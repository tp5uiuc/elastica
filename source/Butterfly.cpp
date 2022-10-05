
#include "Butterfly.h"
#include <string>
#include <algorithm>

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

Rod* butterflyRod(
    const int n, const REAL totalMass, const REAL r0, const Matrix3 _J0,
    const Matrix3 _B0, const Matrix3 _S0, const REAL L0, const REAL inclination,
    const Vector3 origin, const Vector3 direction, const Vector3 normal,
    const REAL nu, const REAL relaxationNu, const bool useSelfContact) {
  // Bunch of sanity checks
  assert(n > 1);
  assert(totalMass > Tolerance::tol());
  assert(r0 > Tolerance::tol());
  assert(L0 > Tolerance::tol());
  assert(nu >= 0.0);
  assert(relaxationNu >= 0.0);
  assert(direction.length() > Tolerance::tol());
  assert(normal.length() > Tolerance::tol());

  assert(_B0.r2c1 == 0.0);
  assert(_B0.r3c1 == 0.0);
  assert(_B0.r1c2 == 0.0);
  assert(_B0.r3c2 == 0.0);
  assert(_B0.r1c3 == 0.0);
  assert(_B0.r2c3 == 0.0);
  assert(_B0.r1c1 >= Tolerance::tol());
  assert(_B0.r2c2 >= Tolerance::tol());
  assert(_B0.r3c3 >= Tolerance::tol());

  assert(_S0.r2c1 == 0.0);
  assert(_S0.r3c1 == 0.0);
  assert(_S0.r1c2 == 0.0);
  assert(_S0.r3c2 == 0.0);
  assert(_S0.r1c3 == 0.0);
  assert(_S0.r2c3 == 0.0);
  assert(_S0.r1c1 >= Tolerance::tol());
  assert(_S0.r2c2 >= Tolerance::tol());
  assert(_S0.r3c3 >= Tolerance::tol());

  // Density
  const REAL density = totalMass / (M_PI * r0 * r0 * L0);


  // u0 is a vector orthogonal to the tangent, .i.e. the direction.
  // If not provided I just compute one, using a random non-zero vector.
  // Note that to be able to compare codes, I fix the random vector and
  // I do not generate one, because that is BAD!
// #ifndef NDEBUG
//   const Vector3 t0 = (x0[1] - x0[0]).unitize();
//   assert((t0 * normal).length() > Tolerance::tol());
// #endif

  const auto n_elem(n);
  const REAL base_length{L0};
  const auto horizontal_direction(direction);
  const auto vertical_direction(normal);
  const REAL bend_angle(inclination);

  // E++ logic
  const REAL dl = base_length / n;
  const auto half_n_elem = n_elem / 2UL;

  auto position_profile = [=](std::size_t index) -> Vector3 {
    const auto indices = std::minmax({index, half_n_elem});
    return origin +
           (dl * static_cast<REAL>(indices.first) *
            (std::cos(bend_angle) * horizontal_direction +
             std::sin(bend_angle) * vertical_direction)) +
           (dl * static_cast<REAL>(indices.second - half_n_elem) *
            (std::cos(bend_angle) * horizontal_direction -
             std::sin(bend_angle) * vertical_direction));
  };


  auto director_profile = [=](std::size_t index) -> Matrix3 {
    // clang-format off
        if (index < half_n_elem) {
          return Matrix3(
              {Vector3{ 0.0            , 1.0, 0.0                 },
               Vector3{ +std::cos(bend_angle), 0.0, -std::sin(bend_angle)},
               Vector3{ +std::sin(bend_angle), 0.0, +std::cos(bend_angle)}}
          );
        } else {
          return Matrix3(
              {Vector3{ 0.0            , 1.0, 0.0                 },
               Vector3{ +std::cos(bend_angle), 0.0, +std::sin(bend_angle)},
               Vector3{ -std::sin(bend_angle), 0.0, +std::cos(bend_angle)}}
          );
        }
    // clang-format on
  };

  // Initialize discretization points
  // TO BE CHANGED
  vector<Vector3> x0(n + 1);
  std::generate(std::begin(x0), std::end(x0), [=, idx = 0]() mutable {
    return position_profile(idx++);
  });

  // Set velocities to zero
  vector<Vector3> v0 = vector<Vector3>(n + 1);

  // Now I can align frames using the orthogonal vector provided/generated
  // TO BE CHANGED
  // vector<Matrix3> Q0 = alignFrames(x0, normal.unitize());
  vector<Matrix3> Q0(n);
  std::generate(std::begin(Q0), std::end(Q0), [=, idx = 0] () mutable {
    return director_profile(idx++);
  });

  // Apply twist about the orthonormal vector d3
  // applyTwists(Q0, vRange((REAL)0.0, totTwist, n));

  // Set angular velocity to zero
  vector<Vector3> w0 = vector<Vector3>(n);

  // Set rest edge lengths
  vector<REAL> l0 = vLength(vDiff(x0));
  const REAL dl0 = l0[0];

  // Set volume discretization elements
  const vector<REAL> V0 = vector<REAL>(n, M_PI * r0 * r0 * dl0);

  // Set shear vector to zero
  vector<Vector3> intrinsicShearStrain0 = vector<Vector3>(n, Vector3(0.0, 0.0, 0.));

  // Mass of vertex point wise element.
  // VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
  // This is why m is obtained dividing the total mass by n, and then to
  // conserve the total mass the masses of the first and last verteces are
  // divideb by 2
  const REAL m = totalMass / (double)(n);
  vector<REAL> masses = vector<REAL>(n + 1, m);
  masses.front() /= 2.0;
  masses.back() /= 2.0;

  // Intrinsic curvature in rest configuration
  const vector<Vector3> intrinsic_k0 = vector<Vector3>(n - 1, Vector3(0.0, 0.0, 0.));

  // Mass second moment of inertia matrix in rest configuration
  vector<Matrix3> J0 = vector<Matrix3>(n, _J0);

  // Bending matrix in reference configuration
  vector<Matrix3> B0 = vector<Matrix3>(n - 1, _B0);

  // Shear matrix in reference configuration
  vector<Matrix3> S0 = vector<Matrix3>(n, _S0);

  return new Rod(n, x0, v0, Q0, w0, l0, intrinsic_k0, intrinsicShearStrain0,
                 masses, V0, density, J0, B0, S0, nu, relaxationNu,
                 useSelfContact);
}

void printEnergies(vector<Rod *>& rodptrs, const int step, const REAL time) {

  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    auto* rod = rodptrs[i];
    rod->computeEnergies();

    std::printf("%1.10e %1.10e %1.10e %1.10e %1.10e %1.10e\n", time,
            rod->totalInternalEnergy, rod->translationalEnergy, rod->rotationalEnergy,
            rod->bendingEnergy, rod->shearEnergy);
  }
}

Butterfly::Butterfly(const int argc, const char **argv) {
  // do parsing here
  std::string n_rods_as_str = getCmdOption(argc, argv, "-n_rods=");
  std::string n_elements_as_str = getCmdOption(argc, argv, "-n_elements=");
  std::string final_time_as_str = getCmdOption(argc, argv, "-final_time=");
  std::string angle_as_str = getCmdOption(argc, argv, "-inclination=");

  n_rod = std::stoul(n_rods_as_str,nullptr,0);
  n_elems_per_rod = std::stoul(n_elements_as_str,nullptr,0);
  final_time = std::stod(final_time_as_str, nullptr);
  inclination = std::stod(angle_as_str, nullptr);

  std::cout << "n_rod: " << n_rod << '\n';
  std::cout << "n_elems_per_rod: " << n_elems_per_rod << '\n';
  std::cout << "final_time: " << final_time << '\n';
  std::cout << "inclination: " << inclination << '\n';

  inclination *= (M_PI / 180.0);

  std::cout << std::endl;
}

Butterfly::~Butterfly() {}

bool Butterfly::_test(const int nEdges, const REAL _dt,
                                      const REAL _L, const REAL _r,
                                      const REAL _P, const REAL _timeSimulation,
                                      const REAL _E, const REAL _G,
                                      const REAL _rho, const REAL _nu,
                                      const REAL _relaxationNu,
                                      const string outfileName) {
  // Input parameters
  const int n = nEdges;  // number of discretization edges (i.e. n+1 points)
                         // along the entire rod
  const REAL timeSimulation = _timeSimulation;  // total simulation time
  const REAL dt = _dt;                          // time step
  const REAL P = _P;
  const REAL L0 = _L;         // total length of rod [m]
  const REAL r0 = _r;         // radius [m]
  const REAL density = _rho;  // [kg/m^3]
  const REAL E = _E;          // GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
  const REAL G = _G;          // GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
  const REAL nu = _nu;        // Numerical damping viscosity [m^2/s]
  const REAL relaxationNu =
      _relaxationNu;  // relaxation time for exponential decay of nu

  // Dumping frequencies (number of frames/dumps per unit time)
  const unsigned int diagPerUnitTime = 5;
  const unsigned int povrayPerUnitTime = 0;

  // Physical parameters
  const REAL dL0 = L0 / (double)n;  // length of cross-section element
  const REAL A0 = M_PI * r0 * r0;
  const REAL Vol = A0 * L0;
  const REAL totalMass = Vol * density;
  // const REAL initialTotalTwist = 0.0;
  const Vector3 originRod = Vector3(0.0, 0.0, 0.0);
  const Vector3 directionRod = Vector3(0.0, 0.0, 1.0);
  const Vector3 normalRod = Vector3(1.0, 0.0, 0.0);

  // Second moment of area for disk cross section
  const REAL I0_1 = A0 * A0 / (4.0 * M_PI);
  const REAL I0_2 = I0_1;
  const REAL I0_3 = 2.0 * I0_1;
  const Matrix3 I0 = Matrix3(I0_1, 0.0, 0.0, 0.0, I0_2, 0.0, 0.0, 0.0, I0_3);

  // Mass inertia matrix for disk cross section
  const Matrix3 J0 = density * dL0 * I0;

  // Bending matrix (TOD: change this is wrong!!)
  Matrix3 B0 =
      Matrix3(E * I0_1, 0.0, 0.0, 0.0, E * I0_2, 0.0, 0.0, 0.0, G * I0_3);

  // Shear matrix
  REAL const shape_factor = REAL(13.5) / REAL(14.0);
  Matrix3 S0 = Matrix3(shape_factor * G * A0, 0.0, 0.0, 0.0,
                       shape_factor * G * A0, 0.0, 0.0, 0.0, E * A0);

  // Initialize straight rod and pack it into a vector of pointers to rod -->
  // Use linear load-strain (hence the true flag at the end)!!!
  const bool useSelfContact = false;
 

  Rod *rod = butterflyRod(
      n, totalMass, r0, J0, B0, S0, L0, inclination, originRod,
      directionRod, normalRod, nu, relaxationNu, useSelfContact);
 

  vector<Rod *> rodPtrs;
  rodPtrs.push_back(rod);
  rod->update(0.0);
  rod->computeEnergies();
  printEnergies(rodPtrs, 0, 0.0);

  // Pack boundary conditions
  FreeBC freeBC = FreeBC();
  vector<RodBC *> boundaryConditionsPtrs;
  boundaryConditionsPtrs.push_back(&freeBC);

  // Pack all forces together
  vector<ExternalForces*> externalForcesPtrs;
  NoForces endpointsForce = NoForces();
  MultipleForces multipleForces;
  multipleForces.add(&endpointsForce);
  MultipleForces* multipleForcesPtr = multipleForces.get();
  externalForcesPtrs.push_back(multipleForcesPtr);

  // Empty interaction forces (no substrate in this case)
  vector<Interaction *> substrateInteractionsPtrs;

  // No external contact function needed
  vector<pair<int, int>> attachpoint;
  vector<ExternalContact *> externalcontactPtrs;
  ExternalContact externalcontact =
      ExternalContact(rodPtrs, 0.0, 0.0, attachpoint);
  externalcontactPtrs.push_back(&externalcontact);

  // No Simple Connection needed
  vector<SimpleConnection *> simpleconnectionPtrs;
  SimpleConnection simpleconnection = SimpleConnection(rodPtrs);
  simpleconnectionPtrs.push_back(&simpleconnection);

  //-----------------------------------------------------------------------------------------------------------------
  // Set up time integrator
  PolymerIntegrator *integrator = new PositionVerlet2nd(
      rodPtrs, externalForcesPtrs, boundaryConditionsPtrs,
      substrateInteractionsPtrs, externalcontactPtrs, simpleconnectionPtrs);

  // Simulate
  Polymer poly = Polymer(integrator);
  // bool goodRun = true;
  const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime,
                                     povrayPerUnitTime, outfileName);

  // Throw exception if something went wrong
  if (!goodRun)
    throw "not good run in localized helical buckling, what is going on?";

  cout << "total internal energy = " << poly.getTotalEnergy() << endl;
  cout << "total translational energy = " << poly.getTotalTranslationalEnergy()
       << endl;
  cout << "total rotational energy = " << poly.getTotalRotationalEnergy()
       << endl;
  cout << "total bending energy = " << poly.getTotalBendingEnergy()
       << endl;
  cout << "total shear energy = " << poly.getTotalShearEnergy()
       << endl;

  return false;
}

void Butterfly::run() {
  // Variable
  const int nEdges = n_elems_per_rod; // Number of degrees of freedom.
  const REAL L = 3.0;

  const REAL rho = 5000;
  const REAL E = 1e4;
  const REAL G = 1e4 / (1.5);
  const REAL nu = 0.0;
  const REAL dL = L / nEdges;
  const REAL r = 0.25;
  const REAL dt = 0.01 * dL;
  const REAL relaxationNu = 0.0;
  const REAL P = 0.0; // unused

  const REAL simTime = final_time;

  _test(nEdges, dt, L, r, P, simTime, E, G, rho, nu, relaxationNu, "timoshenko_final");

  exit(0);
}
