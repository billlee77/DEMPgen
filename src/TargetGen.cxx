#include "Particle.hxx"
#include "CustomRand.hxx"

#include "TRandom2.h"
#include "TMath.h"

#include "TargetGen.hxx"

using namespace TMath;

TargetGen::TargetGen(double mass_in, bool fermi_in)
{
  mass = mass_in;
  char fname[100] = "TargetParticle";
  double ERange[2] = {0,300};
  double ThetaRange[2] = {0, Pi()};
  double PhiRange[2] = {0, 2*Pi()};
  Rand = new CustomRand(fname, ERange, ThetaRange, PhiRange);

  fermi = fermi_in;

  TargetParticle = new Particle();
  TargetParticle -> SetMass(mass);

  fRandom = new TRandom2(0);

  for(int i=0; i<300; i++) fProb[i]/=2;

}

double TargetGen::FermiMomentum()
{
  fMom;
  kFermi = true;
  while ( kFermi ) {
    fProton_Rand_Mom_Col      = fRandom->Uniform( 0, 300.0);
    fProton_Rand_Mom_Col_Prob = fRandom->Uniform( fProb[299], fProb[0] );
    fProton_Mom_Int           = std::ceil( fProton_Rand_Mom_Col );
    f3He_Value                = fProb[ fProton_Mom_Int - 1 ];

    if ( fProton_Rand_Mom_Col_Prob <= f3He_Value ) {
      fMom = fProton_Rand_Mom_Col;
      kFermi = false;
    }
  }

  return fMom;
}

Particle* TargetGen::GetParticle()
{
  TargetParticle -> SetPxPyPzE(0,0,0,mass);

  if(fermi){
    TargetParticle -> SetThetaPhiP(Rand->Theta(),
                                   Rand->Phi(),
                                   FermiMomentum()
                                   );
  }

  return TargetParticle;
}
