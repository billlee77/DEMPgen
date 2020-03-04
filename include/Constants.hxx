/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Declaration of constants for use throughout the program
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef Constants_H
#define Constants_H

#include "TMath.h"

namespace constants{

  const double electron_mass_mev = 0.511;
  const double proton_mass_mev = 938.272;
  const double pion_mass_mev = 139.57;
  const double neutron_mass_mev =  939.565;

  const double electron_mass_gev = 0.511/1000;
  const double proton_mass_gev = 938.272/1000;
  const double pion_mass_gev = 139.57/1000;
  const double neutron_mass_gev =  939.565/1000;

  const double DEG = 180.0/TMath::Pi();

  const int pid_pion = -211;
  const int pid_prot = 2212;
  const int pid_elec = 11;
  const int pid_neut = 0;
  const int pid_phot = 22;

  //const double EBeam = 11000;
  //const double PSF = ((EBeam*0.9-EBeam*0.1)/1000)
  //  * (28.0-5.0)/DEG * 2*TMath::Pi()
  //  *((EBeam*0.9-EBeam*0.1)/1000)
  //  * (60.0-0.0)/DEG * 2*TMath::Pi();
  const double Lumi = 1.0e36;
  const double nBcm2 = 1.0e-33;
  const double uBcm2 = 1.0e-30;

  const double alpha = 1.0/137.0;


  //Helium

  const double Helium_Z = 2.0;
  const double Helium_A = 4.002602;
  const double Helium_Density = 1.345e-3;

  // Taget Window

  const double Window_Density  = 2.76;    // g/cm^3
  const double Window_Thickness  = 19.4246; // g/cm2 this is equal to my X0
  const double Window_Thickness_RadLen = 7.038;    // cm
  const double Window_Z = 17;
  const double Window_A = 35;


  // Air
  const double Airtotal = (0.000124 + 0.755267
                           + 0.231781 + 0.012827);  //total mass of component of air
  const double Air_Z = int( 6 * 0.000124 / Airtotal + 7 * 0.755267
                            / Airtotal + 8 * 0.231781 / Airtotal + 18* 0.012827 / Airtotal);
  const double Air_A = Air_Z / 0.499;
  const double Air_Density = 1.205e-03;


}

#endif
