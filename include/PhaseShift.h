/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  This file is the C++ adaptation of a Fortran routine for SIMC,
  parameterizing the cross-sections for pion-nucleon elastic scattering.
  The adaptation was done by Zafar Ahmed for use in the original DEMP event
  generator and is included here unmodified.
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <complex> 

using namespace std;
using std::setw;
using std::setprecision;
using std::cout;
using std::cin;
using std::endl;

double fZ0 = 0;
double fZ1 = 0;
double fZ2 = 0;


/* ************************************************** */
/* ! written by A. Shinozaki, Oct/00. */

/* ! this function subroutine returns equation 2.25 (p. 812) of Feshbach, */
/* ! "Theoretical Nuclear Physics".  This is the parameterization of the */
/* ! pi-N phase shifts for pion energies less than 400 MeV (LAB). */

/* ! s = pi-N CM energy^2 */
/* ! k = pi-N CM momentum */
/* ! l = pi-N CM angular momentum */
/* ! ij= parameter selecting which set of phaseshift parameters to use */
double formula( double s , double k, int l, int ij ) {

  ij = ij - 1;
    
  /*   ********************************************************************** */
  /* * Phase shift data are from G. Rowe, M. Solomon and R.H. Landau */
  /* ! PRC 18 (1978) 584. */
  /* ! L_2t2j=          S11      S31      P11       P13      P31      P33 */
  double x[6]   = {    0.44,    0.31,    0.61,      0.23,    0.22,    0.99 };
  double ss0[6] = {   1.550,   1.655,   1.435,     1.815,   1.850,   1.233 };
  double k0[6]  = {   0.477,   0.550,   0.393,     0.656,   0.678,   0.228 };
  double G0[6]  = {   0.105,   0.170,   0.230,     0.255,   0.200,   0.116 };
  double b[6]   = {   0.168,  -0.112,  -0.0571,  -0.0131, -0.0291,   0.114 };
  double c[6]   = { -0.0354, -0.0307,   0.0254,  0.00122, 0.00345, -0.0154 };
  double d[6]   = {  27.E-4,  21.E-4,  -29.E-4,  -0.4E-4, -1.5E-4,  7.2E-4 };

  double xmpi = 0.13957;
  double kom = k/xmpi;
  double ss = s;
  double fFormula = 0;

  if ( k > 0.353 ){
    kom = 0.353 / xmpi;
    ss  = pow( 1.383 , 2);
  }

  double kom2 = pow( kom , 2);
  int l21  = 2 * l + 1;
  double s0 = pow( ss0[ ij ] , 2 );
  
  if ( std::abs( s0 - ss ) < 1.e-10 ) {
    fFormula = 1.0;
  }
  else {
    double tangent = ( ( pow( kom , l21 ) ) *
		       ( b[ ij ] + ( c[ ij ] + d[ ij ] * kom2 ) * kom2 ) +
		       x[ ij ] * ( pow ( ( k / k0[ ij ] ) , l21 ) ) *
		       G0[ ij ] * ss0[ ij ] /( s0 - ss ) );
    
    fFormula = tangent / sqrt( 1. + pow( tangent , 2 ) );
    
  }

  /* !     Note on quadrant corrections: */
  /* !     The phaseshifts (delta_Ltj) are from 0-180 degrees. */
  /* !     The tangent function is positive from 0-90 degrees and negative from 90-180 deg. */
  /* !     Sine is positive from 0-180 deg while Cosine is negative from 90-180 deg. */
  /* !     This means that "formula" will give the wrong sign (negative) for sin(delta) */
  /* !     from 90-180 deg.  However, this works out okay, since the cosine calculated  */
  /* !     in the calling routine is always positive, and the real part of the partial  */
  /* !     wave amplitude (cos*sin) will have the correct sign for 0-180 deg. */
  /* !     The imaginary part (sin**2) is always positive.   */
  /* !     Thus, everything works out correctly. */

  return fFormula;

}


/* ! this routine returns the CM angular distribution parameters for pi-N */
/* ! elastic scattering. */

/* ! based on work by A. Shinozaki */
/* ! gh 02.06.20 */

/* ! ich = Isospin channel (1: pure T=3/2,  2: mixed T=1/2,3/2) */
/* ! q  = pi-N CM momentum (GeV/c) */
/* ! s  = pi-N CM energy^2 (GeV) */
/* ! z0 = cos^0 parameter */
/* ! z1 = cos^1 parameter (in pi-N CM frame) */
/* ! z2 = cos^2 parameter */

void phaseshifts( int chan, double q, double s ) {

  double hbarc2 = 0.389379; //  mb(GeV/c)^2
  double reck2 = hbarc2 / pow( q , 2 );

  double ss31 = formula( s , q , 0 , 2 ); // S31
  double cs31 = sqrt( 1.0 - pow( ss31 , 2 ) );
  std::complex<double> fs31 ( ss31 * cs31 , pow( ss31 , 2 ) );
  
  double sp33 = formula( s , q , 1 , 6 ); // P33
  double cp33 = sqrt( 1.0 - pow( sp33 , 2 ) );
  std::complex<double> fp33 ( sp33 * cp33 , pow( sp33 , 2 ) );
  
  double sp31 = formula( s , q , 1 , 5 ); // P31
  double cp31 = sqrt( 1.0 - pow( sp31 , 2 ) );
  std::complex<double> fp31 ( sp31 * cp31 , pow( sp31 , 2 ) );


  if ( chan == 1 ) {
    std::complex<double> fs11 ( 0., 0. );
    std::complex<double> fp11 ( 0., 0. );
    std::complex<double> fp13 ( 0., 0. );

    // ! equations for A,B,C from Appendix F of Aki's thesis.
    // cos^0
    fZ0 = ( reck2 *
	    std::real( ( fs31 + 2.0 * fs11 ) * std::conj( fs31 + 2.0 * fs11 ) +
		       ( ( fp33 + 2.0 * fp13 ) - ( fp31 + 2.0 * fp11 ) ) *
		       std::conj( ( fp33 + 2.0 * fp13 ) -
				  ( fp31 + 2.0 * fp11 ) ) ) );

    // cos^1
    fZ1 = ( reck2 *
	    2.0 * std::real( ( fs31 + 2.0 * fs11 ) *
			     std::conj( 2.0 * ( fp33 + 2.0 * fp13 ) +
					( fp31 + 2.0 * fp11 ) ) ) );


    // cos^2
    fZ2 = ( reck2 *
	    3.0 * std::real( ( fp33 + 2.0 * fp13 ) *
			     std::conj( fp33 + 2.0 * fp13 ) +
			     2.0 * ( fp33 + 2.0 * fp13 ) *
			     std::conj( fp31 + 2.0 * fp11 ) ) );

  }

  
  // For Pion- . mixed isospin 1/2,3/2 channel (e.g. pi-p)
  if ( chan == 2 ) {
    double ss11 = formula( s , q , 0 , 1 ); // S11
    double cs11 = sqrt( 1.0 - pow( ss11 , 2 ) );  
    std::complex<double> fs11 ( ss11 * cs11 , pow( ss11 , 2 ) );

    double sp13 = formula( s , q , 1 , 4 ); // P13
    double cp13 = sqrt( 1.0 - pow( sp13 , 2 ) ); 
    std::complex<double> fp13 ( sp13 * cp13 , pow( sp13  , 2 ) );
  
    double sp11 = formula( s , q , 1 , 3 ); // P11
    double cp11 = sqrt( 1.0 - pow( sp11 , 2 ) );
    std::complex<double> fp11 ( sp11 * cp11 , pow( sp11 , 2 ) );    
  
    
    // ! equations for A,B,C from Appendix F of Aki's thesis.
    // cos^0
    fZ0 = ( reck2 *
	    std::real( ( fs31 + 2.0 * fs11 ) * std::conj( fs31 + 2.0 * fs11 ) +
		       ( ( fp33 + 2.0 * fp13 ) - ( fp31 + 2.0 * fp11 ) ) *
		       std::conj( ( fp33 + 2.0 * fp13 ) -
				  ( fp31 + 2.0 * fp11 ) ) ) );

    // cos^1
    fZ1 = ( reck2 *
	    2.0 * std::real( ( fs31 + 2.0 * fs11 ) *
			     std::conj( 2.0 * ( fp33 + 2.0 * fp13 ) +
					( fp31 + 2.0 * fp11 ) ) ) );


    // cos^2
    fZ2 = ( reck2 *
	    3.0 * std::real( ( fp33 + 2.0 * fp13 ) *
			     std::conj( fp33 + 2.0 * fp13 ) +
			     2.0 * ( fp33 + 2.0 * fp13 ) *
			     std::conj( fp31 + 2.0 * fp11 ) ) );

    fZ0 = fZ0 / 9;
    fZ1 = fZ1 / 9;
    fZ2 = fZ2 / 9;
  }
	   
}

double getZ0(){ return fZ0; }
double getZ1(){ return fZ1; }
double getZ2(){ return fZ2; }
