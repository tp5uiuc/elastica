#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>

// Add to include path
#include "Vector3.h"
#include "Matrix3.h"
#include "SpeedFunctions.h"

template <typename T> 
using V = std::vector<T>;

struct Holder{
	V<Vector3> w;
	V<Vector3> tempVV3_n;
	V<Matrix3> Q;

	Holder(std::size_t n_values) : w(n_values, Vector3{}),
							tempVV3_n(n_values, Vector3{}),
							Q(n_values, Matrix3{}){}

};

void kernel(Holder * const __restrict__ rod){
	v_a_times_b_equal_c(rod->Q, rod->w, rod->tempVV3_n);                  // in-place	
}

std::string getCmdOption(int argc, const char* argv[], const std::string& option)
{
    std::string cmd;
     for( int i = 0; i < argc; ++i)
     {
          std::string arg = argv[i];
          if(0 == arg.find(option))
          {
               std::size_t found = arg.find_first_of("=");
               cmd =arg.substr(found + 1);
               return cmd;
          }
     }
     return cmd;
}

int main(const int argc, const char **argv) {
	REAL coeffDt = 1.0; 

	std::string n_values_as_str = getCmdOption(argc, argv, "-n_values=");
	std::string n_samples_as_str = getCmdOption(argc, argv, "-n_samples=");

	const std::size_t n_values = std::stoul(n_values_as_str,nullptr,0);
	const std::size_t n_samples = std::stoul(n_samples_as_str,nullptr,0);

	// const std::size_t n_values{1 << 16};  // 4096
	// const std::size_t n_samples{1 << 14};
	std::cout << "(" << n_values << ", " << n_samples << ")" << std::endl;

	Holder holder(n_values);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<REAL> d{1.0, 0.2};

	std::generate(holder.w.begin(), holder.w.end(), [&](){
		return Vector3{d(gen), d(gen), d(gen)};
	});
	std::generate(holder.Q.begin(), holder.Q.end(), [&](){
		const auto angle = d(gen) * M_PI / 3.0;
		const auto c = std::cos(angle);
		const auto s = std::sin(angle);
		const Vector3 r1{c, -s, s + c};
		const Vector3 r2{s, c, c - s};
		const Vector3 r3{s * c, s - c, 2.0 * s};
		return Matrix3{r1, r2, r3};
	});


	for (std::size_t i = 0; i < n_samples; ++i)
		kernel(&holder);
 
#ifdef CHECK_RESULTS
	// Ensure no-nans for sanity
	const bool no_nans = std::none_of(holder.w.cbegin(), holder.w.cend(), [](Vector3 const& v){
		return (std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z));  
	}) && 
	std::none_of(holder.tempVV3_n.cbegin(), holder.tempVV3_n.cend(), [](Vector3 const& v){
		return (std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z));  
	}) && 
	std::none_of(holder.Q.cbegin(), holder.Q.cend(), [](Matrix3 const& m){
		return (std::isnan(m.r1c1)|| std::isnan(m.r1c2)|| std::isnan(m.r1c3)
			|| std::isnan(m.r2c1)|| std::isnan(m.r2c2)|| std::isnan(m.r2c3)|| 
			std::isnan(m.r3c1)|| std::isnan(m.r3c2)|| std::isnan(m.r3c3));
	});

	if (!no_nans) {
		throw std::runtime_error("Nans encountered");
		return EXIT_FAILURE;
	}
#endif

	return EXIT_SUCCESS;
}
