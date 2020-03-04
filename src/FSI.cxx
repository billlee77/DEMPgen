#include "FSI.hxx"

#include "TMath.h"

#include "CustomRand.hxx"
#include "Particle.hxx"
#include "Constants.hxx"
#include "PhaseShift.h"

#include <stdio.h>
#include "json/json.h"

using namespace std;
using namespace constants;
using namespace TMath;

FSI::FSI()
{
  double ERange[2] = {0,1};
  double ThetaRange[2] = {0,Pi()};
  double PhiRange[2] = {0,2*Pi()};
  char fname[100] = "FSI";
  SpherePicker = new CustomRand(fname, ERange, ThetaRange, PhiRange);

  extern Json::Value obj;
  ProtGen = new TargetGen(proton_mass_mev, obj["fermi_momentum"].asBool());

  CoP = new TVector3();

  VertInPion = new Particle(pion_mass_mev, "", pid_pion);
  VertTargProt = new Particle(proton_mass_mev, "", pid_prot);

  CMInPion = new Particle(pion_mass_mev, "", pid_pion);
  CMTargProt = new Particle(proton_mass_mev, "", pid_prot);

  VertOutPion = new Particle(pion_mass_mev, "", pid_pion);
  VertOutProt = new Particle(proton_mass_mev, "", pid_prot);

  CMOutPion = new Particle(pion_mass_mev, "", pid_pion);
  CMOutProt = new Particle(proton_mass_mev, "", pid_prot);

  PhaseShiftWeight = new double(0);

  WilliamsWeight = new double(0);
  DedrickWeight = new double(0);
  CatchenWeight = new double(0);
}

int FSI::Generate()
{

  *VertTargProt = *ProtGen->GetParticle();

  *CMInPion = *VertInPion;
  *CMTargProt = *VertTargProt;

  //cout << "Pion Momentum Before:\t"<<CMInPion->P() << endl;
  //cout << "Proton Momentum Before:\t"<<CMTargProt->P()<<endl;

  *CoP = ((VertInPion->Vect()+VertTargProt->Vect())*
          (1.0/(VertInPion->E()+VertTargProt->E())));


  CMInPion->Boost(-(*CoP));
  //cout << "Pion Momentum After:\t"<<CMInPion->P() << endl;
  CMTargProt->Boost(-(*CoP));
  //cout << "Proton Momentum After:\t"<<CMTargProt->P()<<endl;

  theta_pion = SpherePicker->Theta();
  phi_pion = SpherePicker->Phi();

  // Some factors to simplify the formula
  a = Sqrt(Power(pion_mass_mev, 2)+Power(CMInPion->P(), 2))
    + Sqrt(Power(proton_mass_mev, 2)+Power(CMTargProt->P(), 2));
  //cout<<a<<endl;
  b = Power(pion_mass_mev, 2);
  //cout<<b<<endl;
  c = Power(proton_mass_mev, 2);
  //cout<<c<<endl;

  // Solution to E_i = E_f, taking advantage of the fact that the momenta of
  // outgoing particles are equal ond opposite in the CoM frame, leaving
  // pion momentum the only unknown.
  x = (a*a*a*a+b*b+c*c-2*a*a*b-2*a*a*c-2*b*c)/(4*a*a);

  //cout<<x<<endl;

  p_pion = Sqrt(x);

  //cout<<p_pion<<endl;

  CMOutPion->SetThetaPhiP(theta_pion, phi_pion, p_pion);
  //cout << "Pion final P:\t" << CMOutPion->P() << endl;
  CMOutProt->SetVectM(-(CMOutPion->Vect()), proton_mass_mev);
  //cout << "Prot final P:\t" << CMOutProt->P() << endl;

  *VertOutPion = *CMOutPion;
  *VertOutProt = *CMOutProt;

  VertOutProt->Boost((*CoP));
  VertOutPion->Boost((*CoP));

 //cout << "Boost x:\t" << CoP->X() <<endl;
 //cout << "Prot x: \t" << CMOutProt->X()<<endl;
 //cout << "Boost y:\t" << CoP->Y() <<endl;
 //cout << "Prot y: \t" << CMOutProt->Y()<<endl;
 //cout << "Boost z:\t" << CoP->Z() <<endl;
 //cout << "Prot z: \t" << CMOutProt->Z()<<endl;
 //cout << "CM E:\t" << CMOutProt->E()<<endl;
 //cout << "Vert E:\t" << VertOutProt->E()<<endl;

  //Check cons laws:

  if (((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).Px() > 1.0){
    cout << "Px: " << ((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).Px() << endl;
    return 1;
  }
  if (((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).Py() > 1.0){
    cout << "Py: " << ((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).Py() << endl;
    return 1;
  }
  if (((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).Pz() > 1.0){
    cout << "Pz: " << ((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).Pz() << endl;
    return 1;
  }
  if (((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).E() > 1.0){
    cout << "E: " << ((*CMInPion+*CMTargProt)-(*CMOutProt+*CMOutPion)).E() << endl;
      return 1;
    }
    else
      return 0;

}

int FSI::CalculateWeights()
{
  phaseshifts(2, CMOutPion->P()/1000, (*CMInPion+*CMTargProt).Mag2()/1000000);

  //cout << CMOutPion->P() << endl;
  //cout << CMOutPion->Px() << endl;

  //cout << "Check 1" << endl;

  Z0 = getZ0();
  //cout << Z0 << endl;
  Z1 = getZ1();
  //cout << Z1 << endl;
  Z2 = getZ2();
  //cout << Z2 << endl;

  //cout << "Check 2" << endl;

  //cout << (Z0 +
  //        (Z1 * CMOutPion->Px()/CMOutPion->P()) +
  //        (Z2 * Power(CMOutPion->Px()/CMOutPion->P(), 2))
  //         ) << endl;

  *PhaseShiftWeight = (Z0 +
                      (Z1 * CMOutPion->Px()/CMOutPion->P()) +
                      (Z2 * Power(CMOutPion->Px()/CMOutPion->P(), 2))
                      );
  //cout << *PhaseShiftWeight << endl;
  *PhaseShiftWeight *= 0.012; //.012 nucleons per half of He_3 nucleus area in milli barns

  //cout << "Check 3" << endl;

  beta = (VertInPion->P()+VertTargProt->P())/(VertInPion->E()+VertTargProt->E());
  gamma = (VertInPion->E()+VertTargProt->E()) / (*VertInPion+*VertTargProt).Mag();

  //cout << "Check 4" << endl;

  // Jacobian by W. S. C. Williams
  *WilliamsWeight = *PhaseShiftWeight * (VertInPion->Vect().Mag2()/
                                        (gamma*CMOutPion->P()*
                                         (VertInPion->P()-(beta*VertInPion->E()*
                                                           VertInPion->CosTheta())
                                          )
                                         )
                                        );

  //cout << "Check 5" << endl;


  beta_pion = CMOutPion->P() / CMOutPion->E();
  g = beta / beta_pion;

  //cout << "Check 6" << endl;


  // Jacobian by K. G. Dedrick. Rev. Mod. Phys. 34, 429 (1962)
  *DedrickWeight = *PhaseShiftWeight * ((Power
                                        (Power(g+CMOutPion->CosTheta(), 2)+
                                         (1 - beta * beta)*
                                         (1 - Power(CMOutPion->CosTheta(), 2)),
                                         1.5))/
                                       ((1-beta*beta)*Abs(1+g*CMOutPion->CosTheta())
                                        )
                                       );

  //cout << "Check 8" << endl;


  // Jacobian by Gary L. Catchen J. Chem. Phys. 69(4), 15 Aug 1978
  *CatchenWeight = *PhaseShiftWeight * (Power(VertInPion->P(),2)*CMOutPion->E()/
                                       (Power(CMOutPion->P(),2)*VertInPion->E())
                                       );

  //cout << "Check 9" << endl;


  return 0;

}

int FSI::GenerateNoRandom(){
  //Debug only
  //For this function, VertInPion and VertOutPion must be specified
  *VertTargProt = *ProtGen->GetParticle();

  *CMInPion = *VertInPion;
  *CMTargProt = *VertTargProt;

  *CoP = ((VertInPion->Vect()+VertTargProt->Vect())*
          (1.0/(VertInPion->E()+VertTargProt->E())));

  CMInPion->Boost(-(*CoP));
  CMTargProt->Boost(-(*CoP));

  *CMOutPion = *VertOutPion;
  CMOutPion->Boost(-(*CoP));

  CMOutProt->SetVectM(-(CMOutPion->Vect()), proton_mass_mev);

  *VertOutProt = *CMOutProt;
  VertOutProt->Boost(*CoP);

  return 0;

}
