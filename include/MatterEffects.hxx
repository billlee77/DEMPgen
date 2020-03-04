/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Rory Evans, Feb. 2018.
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef MatterEffects_H
#define MatterEffects_H

#include "Particle.hxx"

class MatterEffects{
private:
  double eta(double aZ);

  double Lrad,Lrad_prime,f_Z, a, aE0, mass;
  double lK, lbetasq, lxi, lhbarwsq, j, Delta_p, lw;
  double result;
  Particle * Pout;
  double beta, mstheta, msphi, mstheta0;

public:
  double b(double aZ);
  double X0(double aZ, double aA);
  void MultiScatter(Particle* P, double radlen);
  void IonLoss(Particle* P, double a, double z, double rho, double t);
  void BremLoss(Particle* P, double abt);
};

#endif
