#include "DEMPEvent.hxx"
#include "Particle.hxx"
#include "Constants.hxx"

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "json/json.h"

using namespace constants;

DEMPEvent::DEMPEvent()
{
  qsq_GeV = new double(0);
  t_GeV = new double(0);
  w_GeV = new double(0);
  t_para_GeV = new double(0);
  t_prime_GeV = new double(0);
  negt = new double(0);
  x_B = new double(0);
  Phi = new double(0);
  Phi_s = new double(0);
  Theta = new double(0);
  Phi_deg = new double(0);
  Phi_s_deg = new double(0);
  Theta_deg = new double(0);
  P_T = new double(0);
  Vertex_x = new double(0);
  Vertex_y = new double(0);
  Vertex_z = new double(0);

  STphi = TMath::Pi();
  //std::cout<<STphi<<std::endl;

  extern Json::Value obj;

  TargPol = new TLorentzVector(obj["targ_pol_x"].asDouble(),
                               obj["targ_pol_y"].asDouble(),
                               obj["targ_pol_z"].asDouble(),
                               0);

}

DEMPEvent::DEMPEvent(const char* prefix)
{
  qsq_GeV = new double(0);
  t_GeV = new double(0);
  w_GeV = new double(0);
  t_para_GeV = new double(0);
  t_prime_GeV = new double(0);
  negt = new double(0);
  x_B = new double(0);
  Phi = new double(0);
  Phi_s = new double(0);
  Theta = new double(0);
  Phi_deg = new double(0);
  Phi_s_deg = new double(0);
  Theta_deg = new double(0);
  P_T = new double(0);
  Vertex_x = new double(0);
  Vertex_y = new double(0);
  Vertex_z = new double(0);

  STphi = TMath::Pi();
  //std::cout<<STphi<<std::endl;

  char be_name[100];
  strcpy(be_name, prefix);
  strcat(be_name, "BeamElec");
  BeamElec = new Particle(electron_mass_mev, be_name, pid_elec);
  BeamElec->SetCharge(-1);

  char tn_name[100];
  strcpy(tn_name, prefix);
  strcat(tn_name, "TargNeut");
  TargNeut = new Particle(neutron_mass_mev, tn_name, pid_neut);
  TargNeut->SetCharge(0);

  char se_name[100];
  strcpy(se_name, prefix);
  strcat(se_name, "ScatElec");
  ScatElec = new Particle(electron_mass_mev, se_name, pid_elec);
  ScatElec->SetCharge(-1);

  char ppi_name[100];
  strcpy(ppi_name, prefix);
  strcat(ppi_name, "ProdPion");
  ProdPion = new Particle(pion_mass_mev, ppi_name, pid_pion);
  ProdPion->SetCharge(-1);

  char ppr_name[100];
  strcpy(ppr_name, prefix);
  strcat(ppr_name, "ProdProt");
  ProdProt = new Particle(proton_mass_mev, ppr_name, pid_prot);
  ProdProt->SetCharge(1);

  /*
  char vp_name[100];
  strcpy(vp_name, prefix);
  strcat(vp_name, "VirtPhot");
  ProdProt = new Particle(0, vp_name, pid_prot);
  ProdProt->SetCharge(0);
  */
  VirtPhot = new Particle();

  extern Json::Value obj;

  TargPol = new TLorentzVector(obj["targ_pol_x"].asDouble(),
                               obj["targ_pol_y"].asDouble(),
                               obj["targ_pol_z"].asDouble(),
                               0);

}

void DEMPEvent::Update()
{
  *VirtPhot = *BeamElec - *ScatElec;

  *qsq_GeV = -(VirtPhot->Mag2())/1000000;

  *w_GeV = (*VirtPhot+*TargNeut).Mag()/1000;

  *t_GeV = (*VirtPhot-*ProdPion).Mag2()/1000000;


  double E1 = VirtPhot->E();
  double P1 = VirtPhot->P();

  double E2 = ProdPion->E();
  double P2 = ProdPion->P();

  *t_para_GeV = (((E1-E2)*(E1-E2))-((P1-P2)*(P1-P2)))/1000000;

  *t_prime_GeV = (*t_GeV) - (*t_para_GeV);

  *negt = (-1.0)*(*t_GeV);

  *x_B = (*qsq_GeV)*1000000/(2*proton_mass_mev *
                             this->VirtPhot->E());

  *Phi_s = TargPol->Phi() - ScatElec->Phi();
  *Phi_s_deg = *Phi_s * DEG;

  *Phi = ProdPion->Phi();
  *Phi_deg = *Phi * DEG;

  *Theta = VirtPhot->Theta();
  *Theta_deg = *Theta * DEG;

  *P_T = TMath::Abs(TargPol->Perp(VirtPhot->Vect()));

}

// Calculate and return the center of momentum 3-vector
TVector3 DEMPEvent::CoM()
{
  return (BeamElec->Vect() + TargNeut->Vect())*(1/(BeamElec->E()+TargNeut->E()));
}

// Copy each particle object from another event.
DEMPEvent DEMPEvent::operator = (const DEMPEvent& q)
{
  *(this->BeamElec) = *(q.BeamElec);
  *(this->TargNeut) = *(q.TargNeut);
  *(this->ScatElec) = *(q.ScatElec);
  *(this->ProdPion) = *(q.ProdPion);
  *(this->ProdProt) = *(q.ProdProt);
  *(this->VirtPhot) = *(q.VirtPhot);
  *(this->TargPol) = *(q.TargPol);

  this->STphi = q.STphi;
}

void DEMPEvent::Boost(TVector3 boostvect)
{
  //std::cout << boostvect.Z() << std::endl;
  //std::cout << BeamElec->T() << std::endl;
  BeamElec->Boost(boostvect);
  //std::cout << BeamElec->T() << std::endl;
  TargNeut->Boost(boostvect);
  ScatElec->Boost(boostvect);
  ProdProt->Boost(boostvect);
  ProdPion->Boost(boostvect);
  VirtPhot->Boost(boostvect);
  TargPol->Boost(boostvect);
}

void DEMPEvent::Rotate(double rottheta, double rotphi)
{

  //std::cout << ScatElec->Phi() << std::endl;

  BeamElec->RotateZ(rotphi);
  TargNeut->RotateZ(rotphi);
  ScatElec->RotateZ(rotphi);
  ProdProt->RotateZ(rotphi);
  ProdPion->RotateZ(rotphi);
  VirtPhot->RotateZ(rotphi);

  TargPol->RotateZ(rotphi);

  //std::cout << ScatElec->Phi() << std::endl;

  //std::cout<< VirtPhot->Theta()*DEG << std::endl;

  BeamElec->RotateY(rottheta);
  TargNeut->RotateY(rottheta);
  ScatElec->RotateY(rottheta);
  ProdProt->RotateY(rottheta);
  ProdPion->RotateY(rottheta);
  VirtPhot->RotateY(rottheta);

  TargPol->RotateY(rottheta);

  //std::cout<< VirtPhot->Theta()*DEG << std::endl;

  //std::cout<<STphi*DEG<<std::endl;

  STphi += rotphi;
  //std::cout<<STphi*DEG<<std::endl;

}
