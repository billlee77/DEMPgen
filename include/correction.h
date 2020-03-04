#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>

using std::setw;
using std::setprecision;
using std::cout;
using std::cin;
using std::endl;
using namespace std;


double correctionToSigTT( double Sig_TT, double Sig_T , double qsq ) { 

  double fCheck = Sig_TT;

  if ( ( qsq < 0.5 ) &&                  ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.977833;  
  if ( ( qsq > 0.5 ) && ( qsq < 1   ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.823989;
  if ( ( qsq > 1   ) && ( qsq < 1.5 ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.669703;
  if ( ( qsq > 1.5 ) && ( qsq < 2   ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.558281;
  if ( ( qsq > 2   ) && ( qsq < 2.5 ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.475644;
  if ( ( qsq > 2.5 ) && ( qsq < 3   ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.411432;
  if ( ( qsq > 3   ) && ( qsq < 3.5 ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.359239;
  if ( ( qsq > 3.5 ) && ( qsq < 4   ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.315289;
  if ( ( qsq > 4   ) && ( qsq < 4.5 ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.27723;
  if ( ( qsq > 4.5 ) && ( qsq < 5   ) && ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.243903;
  if ( ( qsq > 5   ) &&                  ( std::abs( Sig_TT ) > std::abs( Sig_T ) ) ) Sig_TT = Sig_T * 0.219363;
  
  double s =  Sig_TT;

  if ( fCheck < 0 && s > 0 )
    s = -1.0*s;

  return s;
}


double correctionToSigLT( double Sig_LT, double Sig_T , double qsq ) { 

  double fCheck = Sig_LT;

  if ( ( qsq < 0.5 )                  && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.200121;
  if ( ( qsq > 0.5 ) && ( qsq < 1   ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.65716;
  if ( ( qsq > 1   ) && ( qsq < 1.5 ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.744639;
  if ( ( qsq > 1.5 ) && ( qsq < 2   ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.716045;
  if ( ( qsq > 2   ) && ( qsq < 2.5 ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.665329;
  if ( ( qsq > 2.5 ) && ( qsq < 3   ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.615555;
  if ( ( qsq > 3   ) && ( qsq < 3.5 ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.571544;
  if ( ( qsq > 3.5 ) && ( qsq < 4   ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.53366;
  if ( ( qsq > 4   ) && ( qsq < 4.5 ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.513469;
  if ( ( qsq > 4.5 ) && ( qsq < 5   ) && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.498549;
  if ( ( qsq > 5   )                  && ( std::abs( Sig_LT ) > std::abs( Sig_T ) ) ) Sig_LT = Sig_T * 0.47498;

  double s =  Sig_LT;

  if ( fCheck < 0 && s > 0 )
    s = -1.0*s;

  return s;
}
