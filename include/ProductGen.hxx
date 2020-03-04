#ifndef ProductGen_H
#define ProductGen_H

#include "Particle.hxx"
#include "CustomRand.hxx"

#include "TF1.h"
#include "TVector3.h"

#include <vector>


using namespace std;

class ProductGen
{

  const double StepSize=0.001;
  const double MTol=0.01; //As fraction of M_p
  const double Start=0.01; //As fraction of photon P
  const double End=0.99; //"

  double W_in_val;

  Particle * Target;
  Particle * Interaction;
  Particle * Initial;

  Particle * Proton;
  Particle * Pion;
  Particle * Final;

  CustomRand * AngleGen;
  TRandom3 * CoinToss;

  TF1 * F;
  TVector3 * UnitVect;

  double Qsq_in();
  double W_in();

  double W_out();

  bool SolnCheck();

public:
  ProductGen(Particle* inInteraction, Particle* inTarget);

  void SetInteraction(Particle * inInteraction);
  void SetTarget(Particle * inTarget);

  int Solve();
  int Solve(double theta, double phi);

  int TestSolve();

  int SolveAnalytic();

  Particle * ProdPion();
  Particle * ProdProton();

  void PrintPars();
};

#endif
