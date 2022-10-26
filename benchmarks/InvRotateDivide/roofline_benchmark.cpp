#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

// Add to include path
#include "Vector3.h"
#include "Matrix3.h"
#include "SpeedFunctions.h"

template <typename T> 
using V = std::vector<T>;

struct Holder{
  	V<Vector3> k0;
  	V<REAL> d0;
	V<Vector3> tempVV3_nminus;
	V<Matrix3> tempVM3_nminus;
  	V<Matrix3> Q;

  	Holder(std::size_t n_values) : k0(n_values - 1, Vector3{}), 
  								   d0(n_values - 1, 0.0),
								   tempVV3_nminus(n_values - 1UL, Vector3{}),
								   tempVM3_nminus(n_values - 1UL, Matrix3{}),
								   Q(n_values, Matrix3{}){}

};

void kernel(Holder * const __restrict__ rod){
	vRotDiff(rod->Q, rod->tempVM3_nminus);
	vLog(rod->tempVM3_nminus, rod->tempVV3_nminus);
	v_a_divide_b_equal_c(rod->tempVV3_nminus, rod->d0, rod->k0);
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
	std::string n_values_as_str = getCmdOption(argc, argv, "-n_values=");
	std::string n_samples_as_str = getCmdOption(argc, argv, "-n_samples=");
	std::string angle_as_str = getCmdOption(argc, argv, "-inclination=");

	const std::size_t n_values = std::stoul(n_values_as_str,nullptr,0);
	const std::size_t n_samples = std::stoul(n_samples_as_str,nullptr,0);
	const std::size_t increment = std::stod(angle_as_str, nullptr);

	// const std::size_t n_values{1 << 16};  // 4096
	// const std::size_t n_samples{1 << 14};

	std::cout << "(" << n_values << ", " << n_samples << ")" << std::endl;

	Holder holder(n_values);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<REAL> d{1.0, 0.1};

	std::generate(holder.d0.begin(), holder.d0.end(), [&](){
		return d(gen);
	});
	std::size_t n = 0UL;
	std::generate(holder.Q.begin(), holder.Q.end(), [&]() mutable {
		const auto angle = REAL(n) * increment;
		const Vector3 r1{std::cos(angle), -std::sin(angle), 0.0};
		const Vector3 r2{std::sin(angle), std::cos(angle), 0.0};
		const Vector3 r3{0.0, 0.0, 1.0};
		++n;
		return Matrix3{r1, r2, r3};
	});


	for (std::size_t i = 0; i < n_samples; ++i)
		kernel(&holder);
 
// #define CHECK_RESULTS
#ifdef CHECK_RESULTS
	// Ensure no-nans for sanity
	const bool no_nans = std::none_of(holder.k0.cbegin(), holder.k0.cend(), [](Vector3 const& v){
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
