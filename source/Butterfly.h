/*
 * BendingWavesCouple.h
 *
 *  Created on: Jun 21, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_Butterfly_H_
#define SOURCE_Butterfly_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodBoundaryConditions.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "Tolerance.h"
#include "UsualHeaders.h"

using namespace std;

class Butterfly : public Test {
 protected:
  bool _test(const int nEdges, const REAL _dt, const REAL _L,
                    const REAL _r, const REAL _P, const REAL _timeSimulation,
                    const REAL _E, const REAL _G, const REAL _rho,
                    const REAL _nu, const REAL _relaxationNu,
                    const string outfileName);

 public:
  Butterfly(const int argc, const char **argv);
  virtual ~Butterfly();

  void run();
  void paint(){};
  unsigned int n_rod;
  unsigned int n_elems_per_rod;
  REAL final_time;
  REAL inclination{0.0};
};

#endif /* SOURCE_Butterfly_H_ */
