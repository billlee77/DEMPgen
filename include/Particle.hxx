/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Class to contain all information on each particle in the event.
  Includes: TLorentz vector with momenta information.
  Charge, PID, etc.
  Access functions for E, theta, etc.
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef Particle_H
#define Particle_H

#pragma clang diagnostic ignored "-Woverloaded-virtual"

#include "TLorentzVector.h"

class Particle:public TLorentzVector
{
  double vx, vy, vz;

  char identifier[100];


public:
  int pid;
  int charge;
  double proper_mass;

  int GetPid();
  int GetCharge();
  double GetMass();

  void SetPid(int x);
  void SetCharge(int x);
  void SetMass(double x);

  void SetVx(double x);
  void SetVy(double x);
  void SetVz(double x);

  double GetVx();
  double GetVy();
  double GetVz();

  double Pmag();

  int Complete(Particle a, Particle b);
  // Set particle momentum to sum zero with
  // the supplied particle
  // (i.e. p_this + p_a + P_b = 0


  Particle(double m, const char* name, int pid_in);

  Particle(double m, double px, double py, double pz);

  Particle(double m, TVector3& v);

  Particle(double m, Particle a, Particle b);

  Particle():TLorentzVector(){}

  Particle operator + (const Particle& q);
  Particle operator - (const Particle& q);
  Particle operator = (const Particle& q);
  Particle operator -();

  char * GetName();
  void SetName(char * name);

  void SetThetaPhiE(double theta, double phi, double E);
  void SetThetaPhiP(double theta, double phi, double P);

};
#endif
