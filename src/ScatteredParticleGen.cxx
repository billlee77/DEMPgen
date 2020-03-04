#include "Particle.hxx"
#include "ScatteredParticleGen.hxx"
#include "CustomRand.hxx"

#include "TMath.h"

using namespace TMath;

ScatteredParticleGen::ScatteredParticleGen(double mass_in,
                                          double ERange[2],
                                          double ThetaRange[2],
                                          double PhiRange[2])
{
  mass = mass_in;
  char fname[100] = "ScatteredParticle";
  Rand = new CustomRand(fname, ERange, ThetaRange, PhiRange);
  ScatteredParticle = new Particle();
}

Particle * ScatteredParticleGen::GetParticle()
{
  double theta, phi, E; //Randomly generated vars
  double px, py, pz;    //Input vars for Particle cons.
  double p;             //Mom. Mag.

  E = Rand->E();
  theta = Rand->Theta();
  phi = Rand->Phi();

  p = Sqrt(E*E - mass*mass);

  px = p*Sin(theta)*Cos(phi);
  py = p*Sin(theta)*Sin(phi);
  pz = p*Cos(theta);

  ScatteredParticle->SetPxPyPzE(px,py,pz,E);
  return ScatteredParticle;
}
