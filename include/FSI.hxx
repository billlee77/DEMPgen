#ifndef FSI_H
#define FSI_H

#include "Particle.hxx"
#include "CustomRand.hxx"
#include "TargetGen.hxx"

#include "TVector3.h"

class FSI
{

private:
  CustomRand* SpherePicker;

  TargetGen* ProtGen;

  double theta_pion, p_pion, phi_pion;
  double a,b,c,x;
  double Z0,Z1,Z2;
  double beta, gamma, beta_pion, g;

  TVector3* CoP;
public:

  Particle* VertInPion;
  Particle* VertTargProt;
  Particle* CMInPion;
  Particle* CMTargProt;

  Particle* VertOutPion;
  Particle* VertOutProt;
  Particle* CMOutPion;
  Particle* CMOutProt;

  double* PhaseShiftWeight;

  double* WilliamsWeight;
  double* DedrickWeight;
  double* CatchenWeight;


  FSI();

  int Generate();
  int CalculateWeights();
  int GenerateNoRandom();

};

#endif
