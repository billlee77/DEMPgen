/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Class to generate primary scattered particle
   (electron) within set paramters.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Particle.hxx"
#include "CustomRand.hxx"

class ScatteredParticleGen{
  CustomRand * Rand;
  double mass;
  Particle * ScatteredParticle;

public:
  Particle* GetParticle();

  ScatteredParticleGen(double mass,
                       double ERange[2],
                       double ThetaRange[2],
                       double PhiRange[2]);
};
