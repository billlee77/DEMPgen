#ifndef DEMPEvent_H
#define DEMPEvent_H

#include "TLorentzVector.h"

#include "Particle.hxx"
#include "TVector3.h"

class DEMPEvent
{
public:

  DEMPEvent();
  DEMPEvent(const char* prefix);

  void Boost(TVector3 boostvect);
  void Rotate(double rottheta, double rotphi);
  TVector3 CoM();

  Particle * BeamElec;
  Particle * TargNeut;
  Particle * ScatElec;
  Particle * ProdProt;
  Particle * ProdPion;

  Particle * VirtPhot;

  TLorentzVector* TargPol;

  double* qsq_GeV;
  double* t_GeV;
  double* w_GeV;

  double* t_para_GeV;
  double* t_prime_GeV;
  double* negt;
  double* x_B;

  double* Phi;
  double* Phi_s;
  double* Theta;

  double* Phi_deg;
  double* Phi_s_deg;
  double* Theta_deg;

  double* P_T;

  double* Vertex_x;
  double* Vertex_y;
  double* Vertex_z;

  double STphi; //Angle of the target transverse polarization

  void Update();

  DEMPEvent operator = (const DEMPEvent& q);

};

#endif
