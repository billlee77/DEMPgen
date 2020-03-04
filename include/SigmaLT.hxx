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

/* double fSigmaLT; */

double MySigmaLT( double qsqgev, double tgev, double wgev ) { 

  double s      =  0.0;
  /* double wgev   = par[0]; */
  /* double qsqgev = par[1]; */
  /* double tgev   = x[0]; */



  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/

  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -1.942 + pow( qsqgev , 1 ) * 2.653 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.551 + pow( qsqgev , 1 ) * 1.036 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 16.84 + pow( qsqgev , 1 ) * -16.19 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -10.18 + pow( qsqgev , 1 ) * 10.7 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -1.155 + pow( qsqgev , 1 ) * 1.895 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.919 + pow( qsqgev , 1 ) * 1.48 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 14.55 + pow( qsqgev , 1 ) * -13.39 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -8.63 + pow( qsqgev , 1 ) * 8.762 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.6923 + pow( qsqgev , 1 ) * 1.457 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.824 + pow( qsqgev , 1 ) * 1.289 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 12.38 + pow( qsqgev , 1 ) * -10.61 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -7.243 + pow( qsqgev , 1 ) * 6.979 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.4203 + pow( qsqgev , 1 ) * 1.193 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.498 + pow( qsqgev , 1 ) * 0.847 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 10.42 + pow( qsqgev , 1 ) * -8.111 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -6.034 + pow( qsqgev , 1 ) * 5.431 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2296 + pow( qsqgev , 1 ) * 0.9805 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.273 + pow( qsqgev , 1 ) * 0.6162 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 9.012 + pow( qsqgev , 1 ) * -6.417 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -5.157 + pow( qsqgev , 1 ) * 4.347 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1191 + pow( qsqgev , 1 ) * 0.8388 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.944 + pow( qsqgev , 1 ) * 0.31 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 7.71 + pow( qsqgev , 1 ) * -4.914 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.375 + pow( qsqgev , 1 ) * 3.412 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.05686 + pow( qsqgev , 1 ) * 0.7382 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.584 + pow( qsqgev , 1 ) * 0.005725 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.546 + pow( qsqgev , 1 ) * -3.639 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -3.693 + pow( qsqgev , 1 ) * 2.627 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.002834 + pow( qsqgev , 1 ) * 0.635 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.357 + pow( qsqgev , 1 ) * -0.09286 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 5.738 + pow( qsqgev , 1 ) * -2.868 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -3.209 + pow( qsqgev , 1 ) * 2.119 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.02147 + pow( qsqgev , 1 ) * 0.5634 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.077 + pow( qsqgev , 1 ) * -0.2457 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.942 + pow( qsqgev , 1 ) * -2.129 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.75 + pow( qsqgev , 1 ) * 1.651 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04124 + pow( qsqgev , 1 ) * 0.4949 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.865 + pow( qsqgev , 1 ) * -0.3078 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.329 + pow( qsqgev , 1 ) * -1.634 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.394 + pow( qsqgev , 1 ) * 1.321 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05085 + pow( qsqgev , 1 ) * 0.439 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.659 + pow( qsqgev , 1 ) * -0.3618 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.783 + pow( qsqgev , 1 ) * -1.225 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.082 + pow( qsqgev , 1 ) * 1.047 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05635 + pow( qsqgev , 1 ) * 0.3898 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.482 + pow( qsqgev , 1 ) * -0.3887 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.323 + pow( qsqgev , 1 ) * -0.9168 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.821 + pow( qsqgev , 1 ) * 0.8354 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05989 + pow( qsqgev , 1 ) * 0.3454 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.332 + pow( qsqgev , 1 ) * -0.39 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.943 + pow( qsqgev , 1 ) * -0.6963 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.606 + pow( qsqgev , 1 ) * 0.676 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05982 + pow( qsqgev , 1 ) * 0.3082 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.192 + pow( qsqgev , 1 ) * -0.3917 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.601 + pow( qsqgev , 1 ) * -0.5115 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.414 + pow( qsqgev , 1 ) * 0.5422 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 0.1  to  1  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.06028 + pow( qsqgev , 1 ) * 0.2734 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.08 + pow( qsqgev , 1 ) * -0.3725 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.326 + pow( qsqgev , 1 ) * -0.3905 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.26 + pow( qsqgev , 1 ) * 0.446 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/



  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.214 + pow( qsqgev , 1 ) * -1.046 + pow( qsqgev , 2 ) * 0.09154 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -21.47 + pow( qsqgev , 1 ) * 25.33 + pow( qsqgev , 2 ) * -6.004 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 72 + pow( qsqgev , 1 ) * -92.44 + pow( qsqgev , 2 ) * 24.21 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -109.6 + pow( qsqgev , 1 ) * 145.4 + pow( qsqgev , 2 ) * -39.74 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 78.76 + pow( qsqgev , 1 ) * -106.3 + pow( qsqgev , 2 ) * 29.79 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -21.52 + pow( qsqgev , 1 ) * 29.37 + pow( qsqgev , 2 ) * -8.372 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.47 + pow( qsqgev , 1 ) * -1.395 + pow( qsqgev , 2 ) * 0.2358 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -18.97 + pow( qsqgev , 1 ) * 23.51 + pow( qsqgev , 2 ) * -5.878 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 60.75 + pow( qsqgev , 1 ) * -83.26 + pow( qsqgev , 2 ) * 22.7 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -89.2 + pow( qsqgev , 1 ) * 128.6 + pow( qsqgev , 2 ) * -36.58 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 62.16 + pow( qsqgev , 1 ) * -92.6 + pow( qsqgev , 2 ) * 27.04 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -16.55 + pow( qsqgev , 1 ) * 25.26 + pow( qsqgev , 2 ) * -7.517 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.616 + pow( qsqgev , 1 ) * -1.541 + pow( qsqgev , 2 ) * 0.3027 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -17.75 + pow( qsqgev , 1 ) * 22.18 + pow( qsqgev , 2 ) * -5.695 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 55.78 + pow( qsqgev , 1 ) * -78.06 + pow( qsqgev , 2 ) * 21.71 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -80.74 + pow( qsqgev , 1 ) * 120.5 + pow( qsqgev , 2 ) * -34.96 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 55.46 + pow( qsqgev , 1 ) * -86.62 + pow( qsqgev , 2 ) * 25.83 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -14.57 + pow( qsqgev , 1 ) * 23.58 + pow( qsqgev , 2 ) * -7.175 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.667 + pow( qsqgev , 1 ) * -1.558 + pow( qsqgev , 2 ) * 0.3226 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -16.98 + pow( qsqgev , 1 ) * 20.91 + pow( qsqgev , 2 ) * -5.42 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 53.48 + pow( qsqgev , 1 ) * -74.22 + pow( qsqgev , 2 ) * 20.75 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -77.73 + pow( qsqgev , 1 ) * 115.7 + pow( qsqgev , 2 ) * -33.72 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 53.5 + pow( qsqgev , 1 ) * -83.77 + pow( qsqgev , 2 ) * 25.08 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -14.07 + pow( qsqgev , 1 ) * 22.93 + pow( qsqgev , 2 ) * -7.008 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.649 + pow( qsqgev , 1 ) * -1.495 + pow( qsqgev , 2 ) * 0.3149 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -16.27 + pow( qsqgev , 1 ) * 19.54 + pow( qsqgev , 2 ) * -5.062 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 51.98 + pow( qsqgev , 1 ) * -70.44 + pow( qsqgev , 2 ) * 19.62 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -76.59 + pow( qsqgev , 1 ) * 111.4 + pow( qsqgev , 2 ) * -32.31 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 53.3 + pow( qsqgev , 1 ) * -81.58 + pow( qsqgev , 2 ) * 24.28 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -14.15 + pow( qsqgev , 1 ) * 22.53 + pow( qsqgev , 2 ) * -6.839 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.581 + pow( qsqgev , 1 ) * -1.385 + pow( qsqgev , 2 ) * 0.292 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -15.48 + pow( qsqgev , 1 ) * 18.05 + pow( qsqgev , 2 ) * -4.642 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 50.36 + pow( qsqgev , 1 ) * -66.22 + pow( qsqgev , 2 ) * 18.28 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -75.46 + pow( qsqgev , 1 ) * 106.4 + pow( qsqgev , 2 ) * -30.53 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 53.24 + pow( qsqgev , 1 ) * -78.83 + pow( qsqgev , 2 ) * 23.2 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -14.3 + pow( qsqgev , 1 ) * 21.98 + pow( qsqgev , 2 ) * -6.589 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.483 + pow( qsqgev , 1 ) * -1.253 + pow( qsqgev , 2 ) * 0.2616 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -14.57 + pow( qsqgev , 1 ) * 16.46 + pow( qsqgev , 2 ) * -4.188 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 48.29 + pow( qsqgev , 1 ) * -61.5 + pow( qsqgev , 2 ) * 16.78 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -73.55 + pow( qsqgev , 1 ) * 100.3 + pow( qsqgev , 2 ) * -28.4 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 52.59 + pow( qsqgev , 1 ) * -75.14 + pow( qsqgev , 2 ) * 21.8 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -14.28 + pow( qsqgev , 1 ) * 21.14 + pow( qsqgev , 2 ) * -6.24 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.369 + pow( qsqgev , 1 ) * -1.114 + pow( qsqgev , 2 ) * 0.2289 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -13.56 + pow( qsqgev , 1 ) * 14.84 + pow( qsqgev , 2 ) * -3.726 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 45.74 + pow( qsqgev , 1 ) * -56.42 + pow( qsqgev , 2 ) * 15.18 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -70.7 + pow( qsqgev , 1 ) * 93.24 + pow( qsqgev , 2 ) * -26.02 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 51.15 + pow( qsqgev , 1 ) * -70.56 + pow( qsqgev , 2 ) * 20.16 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -14.03 + pow( qsqgev , 1 ) * 20 + pow( qsqgev , 2 ) * -5.81 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.248 + pow( qsqgev , 1 ) * -0.9779 + pow( qsqgev , 2 ) * 0.1967 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -12.5 + pow( qsqgev , 1 ) * 13.25 + pow( qsqgev , 2 ) * -3.277 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 42.82 + pow( qsqgev , 1 ) * -51.18 + pow( qsqgev , 2 ) * 13.57 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -67.02 + pow( qsqgev , 1 ) * 85.59 + pow( qsqgev , 2 ) * -23.52 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 48.99 + pow( qsqgev , 1 ) * -65.34 + pow( qsqgev , 2 ) * 18.37 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -13.55 + pow( qsqgev , 1 ) * 18.65 + pow( qsqgev , 2 ) * -5.328 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.129 + pow( qsqgev , 1 ) * -0.8498 + pow( qsqgev , 2 ) * 0.1666 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -11.43 + pow( qsqgev , 1 ) * 11.73 + pow( qsqgev , 2 ) * -2.855 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 39.66 + pow( qsqgev , 1 ) * -45.98 + pow( qsqgev , 2 ) * 12 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -62.75 + pow( qsqgev , 1 ) * 77.7 + pow( qsqgev , 2 ) * -21.02 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 46.25 + pow( qsqgev , 1 ) * -59.76 + pow( qsqgev , 2 ) * 16.53 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -12.87 + pow( qsqgev , 1 ) * 17.15 + pow( qsqgev , 2 ) * -4.82 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.014 + pow( qsqgev , 1 ) * -0.7332 + pow( qsqgev , 2 ) * 0.1397 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -10.37 + pow( qsqgev , 1 ) * 10.32 + pow( qsqgev , 2 ) * -2.469 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 36.42 + pow( qsqgev , 1 ) * -40.99 + pow( qsqgev , 2 ) * 10.53 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -58.13 + pow( qsqgev , 1 ) * 69.89 + pow( qsqgev , 2 ) * -18.61 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 43.14 + pow( qsqgev , 1 ) * -54.1 + pow( qsqgev , 2 ) * 14.73 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -12.08 + pow( qsqgev , 1 ) * 15.61 + pow( qsqgev , 2 ) * -4.314 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.907 + pow( qsqgev , 1 ) * -0.6293 + pow( qsqgev , 2 ) * 0.1161 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -9.368 + pow( qsqgev , 1 ) * 9.033 + pow( qsqgev , 2 ) * -2.122 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 33.21 + pow( qsqgev , 1 ) * -36.31 + pow( qsqgev , 2 ) * 9.175 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -53.41 + pow( qsqgev , 1 ) * 62.4 + pow( qsqgev , 2 ) * -16.35 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 39.86 + pow( qsqgev , 1 ) * -48.58 + pow( qsqgev , 2 ) * 13.01 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -11.21 + pow( qsqgev , 1 ) * 14.07 + pow( qsqgev , 2 ) * -3.827 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.809 + pow( qsqgev , 1 ) * -0.5383 + pow( qsqgev , 2 ) * 0.09588 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -8.428 + pow( qsqgev , 1 ) * 7.879 + pow( qsqgev , 2 ) * -1.817 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 30.13 + pow( qsqgev , 1 ) * -32.02 + pow( qsqgev , 2 ) * 7.96 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -48.75 + pow( qsqgev , 1 ) * 55.42 + pow( qsqgev , 2 ) * -14.29 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 36.56 + pow( qsqgev , 1 ) * -43.36 + pow( qsqgev , 2 ) * 11.43 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -10.32 + pow( qsqgev , 1 ) * 12.61 + pow( qsqgev , 2 ) * -3.375 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7205 + pow( qsqgev , 1 ) * -0.4594 + pow( qsqgev , 2 ) * 0.07875 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -7.563 + pow( qsqgev , 1 ) * 6.857 + pow( qsqgev , 2 ) * -1.552 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 27.22 + pow( qsqgev , 1 ) * -28.14 + pow( qsqgev , 2 ) * 6.882 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -44.28 + pow( qsqgev , 1 ) * 49.02 + pow( qsqgev , 2 ) * -12.45 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 33.33 + pow( qsqgev , 1 ) * -38.51 + pow( qsqgev , 2 ) * 10 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -9.436 + pow( qsqgev , 1 ) * 11.23 + pow( qsqgev , 2 ) * -2.961 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 1  to  2  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6413 + pow( qsqgev , 1 ) * -0.3917 + pow( qsqgev , 2 ) * 0.06443 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -6.776 + pow( qsqgev , 1 ) * 5.96 + pow( qsqgev , 2 ) * -1.323 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 24.53 + pow( qsqgev , 1 ) * -24.68 + pow( qsqgev , 2 ) * 5.94 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -40.08 + pow( qsqgev , 1 ) * 43.24 + pow( qsqgev , 2 ) * -10.81 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 30.27 + pow( qsqgev , 1 ) * -34.1 + pow( qsqgev , 2 ) * 8.724 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -8.59 + pow( qsqgev , 1 ) * 9.973 + pow( qsqgev , 2 ) * -2.59 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/

  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1988 + pow( qsqgev , 1 ) * -0.06115 + pow( qsqgev , 2 ) * -0.149 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 1.085 + pow( qsqgev , 1 ) * 1.19 + pow( qsqgev , 2 ) * 0.4906 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -12.13 + pow( qsqgev , 1 ) * -0.2031 + pow( qsqgev , 2 ) * -1.124 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 23.74 + pow( qsqgev , 1 ) * -2.589 + pow( qsqgev , 2 ) * 1.328 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -19.22 + pow( qsqgev , 1 ) * 3.297 + pow( qsqgev , 2 ) * -0.8253 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 5.7 + pow( qsqgev , 1 ) * -1.233 + pow( qsqgev , 2 ) * 0.2167 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.53 + pow( qsqgev , 1 ) * -0.4321 + pow( qsqgev , 2 ) * -0.01254 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.7113 + pow( qsqgev , 1 ) * 2.137 + pow( qsqgev , 2 ) * -0.05752 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -14.12 + pow( qsqgev , 1 ) * -0.3568 + pow( qsqgev , 2 ) * -0.2638 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 30.49 + pow( qsqgev , 1 ) * -5.395 + pow( qsqgev , 2 ) * 0.896 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -25.91 + pow( qsqgev , 1 ) * 6.761 + pow( qsqgev , 2 ) * -0.9182 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 7.875 + pow( qsqgev , 1 ) * -2.445 + pow( qsqgev , 2 ) * 0.315 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6674 + pow( qsqgev , 1 ) * -0.5521 + pow( qsqgev , 2 ) * 0.04307 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.4416 + pow( qsqgev , 1 ) * 2.349 + pow( qsqgev , 2 ) * -0.276 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -14.44 + pow( qsqgev , 1 ) * -0.163 + pow( qsqgev , 2 ) * 0.09797 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 32.9 + pow( qsqgev , 1 ) * -6.897 + pow( qsqgev , 2 ) * 0.711 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -28.74 + pow( qsqgev , 1 ) * 8.506 + pow( qsqgev , 2 ) * -0.9674 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 8.876 + pow( qsqgev , 1 ) * -3.057 + pow( qsqgev , 2 ) * 0.3632 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.738 + pow( qsqgev , 1 ) * -0.5851 + pow( qsqgev , 2 ) * 0.06653 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.0189 + pow( qsqgev , 1 ) * 2.444 + pow( qsqgev , 2 ) * -0.3794 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -13.14 + pow( qsqgev , 1 ) * -0.4828 + pow( qsqgev , 2 ) * 0.3366 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 31.81 + pow( qsqgev , 1 ) * -6.729 + pow( qsqgev , 2 ) * 0.475 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -28.57 + pow( qsqgev , 1 ) * 8.64 + pow( qsqgev , 2 ) * -0.8701 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 8.973 + pow( qsqgev , 1 ) * -3.154 + pow( qsqgev , 2 ) * 0.3515 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7758 + pow( qsqgev , 1 ) * -0.5833 + pow( qsqgev , 2 ) * 0.07566 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.609 + pow( qsqgev , 1 ) * 2.543 + pow( qsqgev , 2 ) * -0.4363 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -10.79 + pow( qsqgev , 1 ) * -1.271 + pow( qsqgev , 2 ) * 0.5464 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 28.4 + pow( qsqgev , 1 ) * -5.33 + pow( qsqgev , 2 ) * 0.1555 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -26.37 + pow( qsqgev , 1 ) * 7.619 + pow( qsqgev , 2 ) * -0.649 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 8.438 + pow( qsqgev , 1 ) * -2.881 + pow( qsqgev , 2 ) * 0.2936 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7928 + pow( qsqgev , 1 ) * -0.5663 + pow( qsqgev , 2 ) * 0.07798 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.24 + pow( qsqgev , 1 ) * 2.656 + pow( qsqgev , 2 ) * -0.4715 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -7.912 + pow( qsqgev , 1 ) * -2.338 + pow( qsqgev , 2 ) * 0.746 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 23.66 + pow( qsqgev , 1 ) * -3.189 + pow( qsqgev , 2 ) * -0.2195 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -22.91 + pow( qsqgev , 1 ) * 5.881 + pow( qsqgev , 2 ) * -0.3493 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 7.497 + pow( qsqgev , 1 ) * -2.372 + pow( qsqgev , 2 ) * 0.206 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7926 + pow( qsqgev , 1 ) * -0.5406 + pow( qsqgev , 2 ) * 0.07659 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.824 + pow( qsqgev , 1 ) * 2.752 + pow( qsqgev , 2 ) * -0.491 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -4.957 + pow( qsqgev , 1 ) * -3.45 + pow( qsqgev , 2 ) * 0.921 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 18.44 + pow( qsqgev , 1 ) * -0.8004 + pow( qsqgev , 2 ) * -0.5938 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -18.9 + pow( qsqgev , 1 ) * 3.851 + pow( qsqgev , 2 ) * -0.02836 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 6.356 + pow( qsqgev , 1 ) * -1.757 + pow( qsqgev , 2 ) * 0.1079 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7791 + pow( qsqgev , 1 ) * -0.5101 + pow( qsqgev , 2 ) * 0.07321 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.324 + pow( qsqgev , 1 ) * 2.816 + pow( qsqgev , 2 ) * -0.4988 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -2.167 + pow( qsqgev , 1 ) * -4.481 + pow( qsqgev , 2 ) * 1.065 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 13.25 + pow( qsqgev , 1 ) * 1.552 + pow( qsqgev , 2 ) * -0.9373 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -14.75 + pow( qsqgev , 1 ) * 1.779 + pow( qsqgev , 2 ) * 0.282 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 5.143 + pow( qsqgev , 1 ) * -1.115 + pow( qsqgev , 2 ) * 0.009857 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7554 + pow( qsqgev , 1 ) * -0.4772 + pow( qsqgev , 2 ) * 0.06878 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.722 + pow( qsqgev , 1 ) * 2.841 + pow( qsqgev , 2 ) * -0.497 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3193 + pow( qsqgev , 1 ) * -5.357 + pow( qsqgev , 2 ) * 1.175 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 8.401 + pow( qsqgev , 1 ) * 3.688 + pow( qsqgev , 2 ) * -1.231 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -10.76 + pow( qsqgev , 1 ) * -0.1652 + pow( qsqgev , 2 ) * 0.5603 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 3.951 + pow( qsqgev , 1 ) * -0.4985 + pow( qsqgev , 2 ) * -0.08053 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.724 + pow( qsqgev , 1 ) * -0.443 + pow( qsqgev , 2 ) * 0.06379 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.011 + pow( qsqgev , 1 ) * 2.822 + pow( qsqgev , 2 ) * -0.4865 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.406 + pow( qsqgev , 1 ) * -6.022 + pow( qsqgev , 2 ) * 1.246 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 4.141 + pow( qsqgev , 1 ) * 5.466 + pow( qsqgev , 2 ) * -1.457 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -7.159 + pow( qsqgev , 1 ) * -1.843 + pow( qsqgev , 2 ) * 0.7874 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 2.856 + pow( qsqgev , 1 ) * 0.04486 + pow( qsqgev , 2 ) * -0.1566 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6874 + pow( qsqgev , 1 ) * -0.4086 + pow( qsqgev , 2 ) * 0.0586 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.199 + pow( qsqgev , 1 ) * 2.765 + pow( qsqgev , 2 ) * -0.4695 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.085 + pow( qsqgev , 1 ) * -6.478 + pow( qsqgev , 2 ) * 1.281 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.5367 + pow( qsqgev , 1 ) * 6.857 + pow( qsqgev , 2 ) * -1.617 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -4.034 + pow( qsqgev , 1 ) * -3.218 + pow( qsqgev , 2 ) * 0.9612 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 1.89 + pow( qsqgev , 1 ) * 0.5008 + pow( qsqgev , 2 ) * -0.2169 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6476 + pow( qsqgev , 1 ) * -0.3749 + pow( qsqgev , 2 ) * 0.05343 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.299 + pow( qsqgev , 1 ) * 2.674 + pow( qsqgev , 2 ) * -0.4474 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 5.375 + pow( qsqgev , 1 ) * -6.737 + pow( qsqgev , 2 ) * 1.285 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.405 + pow( qsqgev , 1 ) * 7.871 + pow( qsqgev , 2 ) * -1.715 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.41 + pow( qsqgev , 1 ) * -4.283 + pow( qsqgev , 2 ) * 1.083 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 1.065 + pow( qsqgev , 1 ) * 0.8649 + pow( qsqgev , 2 ) * -0.2616 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6061 + pow( qsqgev , 1 ) * -0.3423 + pow( qsqgev , 2 ) * 0.04841 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.325 + pow( qsqgev , 1 ) * 2.557 + pow( qsqgev , 2 ) * -0.4215 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.309 + pow( qsqgev , 1 ) * -6.823 + pow( qsqgev , 2 ) * 1.264 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.712 + pow( qsqgev , 1 ) * 8.532 + pow( qsqgev , 2 ) * -1.757 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.7154 + pow( qsqgev , 1 ) * -5.053 + pow( qsqgev , 2 ) * 1.157 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0.3841 + pow( qsqgev , 1 ) * 1.139 + pow( qsqgev , 2 ) * -0.2916 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.5643 + pow( qsqgev , 1 ) * -0.3114 + pow( qsqgev , 2 ) * 0.04366 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.291 + pow( qsqgev , 1 ) * 2.422 + pow( qsqgev , 2 ) * -0.3936 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.942 + pow( qsqgev , 1 ) * -6.769 + pow( qsqgev , 2 ) * 1.222 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -6.458 + pow( qsqgev , 1 ) * 8.899 + pow( qsqgev , 2 ) * -1.756 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 2.391 + pow( qsqgev , 1 ) * -5.568 + pow( qsqgev , 2 ) * 1.193 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -0.1643 + pow( qsqgev , 1 ) * 1.336 + pow( qsqgev , 2 ) * -0.3094 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 2  to  3  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.5231 + pow( qsqgev , 1 ) * -0.2824 + pow( qsqgev , 2 ) * 0.03921 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.211 + pow( qsqgev , 1 ) * 2.274 + pow( qsqgev , 2 ) * -0.3644 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 7.311 + pow( qsqgev , 1 ) * -6.597 + pow( qsqgev , 2 ) * 1.165 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -7.697 + pow( qsqgev , 1 ) * 9.003 + pow( qsqgev , 2 ) * -1.717 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 3.646 + pow( qsqgev , 1 ) * -5.85 + pow( qsqgev , 2 ) * 1.193 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * -0.5861 + pow( qsqgev , 1 ) * 1.458 + pow( qsqgev , 2 ) * -0.3159 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/


  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -7.057 + pow( qsqgev , 1 ) * 4.6 + pow( qsqgev , 2 ) * -0.8992 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 39.74 + pow( qsqgev , 1 ) * -23.87 + pow( qsqgev , 2 ) * 4.562 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -89.29 + pow( qsqgev , 1 ) * 50.05 + pow( qsqgev , 2 ) * -9.332 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 99.97 + pow( qsqgev , 1 ) * -52.5 + pow( qsqgev , 2 ) * 9.526 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -56.53 + pow( qsqgev , 1 ) * 27.87 + pow( qsqgev , 2 ) * -4.885 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 12.89 + pow( qsqgev , 1 ) * -5.999 + pow( qsqgev , 2 ) * 1.01 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -1.978 + pow( qsqgev , 1 ) * 1.201 + pow( qsqgev , 2 ) * -0.2792 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 15.18 + pow( qsqgev , 1 ) * -7.422 + pow( qsqgev , 2 ) * 1.527 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -42.7 + pow( qsqgev , 1 ) * 18.68 + pow( qsqgev , 2 ) * -3.447 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 57.09 + pow( qsqgev , 1 ) * -23.27 + pow( qsqgev , 2 ) * 3.91 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -37.38 + pow( qsqgev , 1 ) * 14.55 + pow( qsqgev , 2 ) * -2.248 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 9.595 + pow( qsqgev , 1 ) * -3.638 + pow( qsqgev , 2 ) * 0.5223 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.5244 + pow( qsqgev , 1 ) * 0.2375 + pow( qsqgev , 2 ) * -0.08818 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 7.926 + pow( qsqgev , 1 ) * -2.706 + pow( qsqgev , 2 ) * 0.5804 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -28.84 + pow( qsqgev , 1 ) * 9.707 + pow( qsqgev , 2 ) * -1.599 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 44.68 + pow( qsqgev , 1 ) * -15.15 + pow( qsqgev , 2 ) * 2.159 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -32.2 + pow( qsqgev , 1 ) * 11.07 + pow( qsqgev , 2 ) * -1.442 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 8.825 + pow( qsqgev , 1 ) * -3.086 + pow( qsqgev , 2 ) * 0.379 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.031 + pow( qsqgev , 1 ) * -0.06721 + pow( qsqgev , 2 ) * -0.02095 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 5.189 + pow( qsqgev , 1 ) * -1.137 + pow( qsqgev , 2 ) * 0.2379 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -23.26 + pow( qsqgev , 1 ) * 6.628 + pow( qsqgev , 2 ) * -0.9143 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 39.52 + pow( qsqgev , 1 ) * -12.38 + pow( qsqgev , 2 ) * 1.506 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -30.05 + pow( qsqgev , 1 ) * 9.947 + pow( qsqgev , 2 ) * -1.144 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 8.522 + pow( qsqgev , 1 ) * -2.931 + pow( qsqgev , 2 ) * 0.3277 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1664 + pow( qsqgev , 1 ) * -0.1687 + pow( qsqgev , 2 ) * 0.004959 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 3.797 + pow( qsqgev , 1 ) * -0.5151 + pow( qsqgev , 2 ) * 0.09557 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -19.93 + pow( qsqgev , 1 ) * 5.225 + pow( qsqgev , 2 ) * -0.608 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 35.99 + pow( qsqgev , 1 ) * -10.96 + pow( qsqgev , 2 ) * 1.194 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -28.31 + pow( qsqgev , 1 ) * 9.288 + pow( qsqgev , 2 ) * -0.9925 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 8.205 + pow( qsqgev , 1 ) * -2.82 + pow( qsqgev , 2 ) * 0.2995 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2598 + pow( qsqgev , 1 ) * -0.2019 + pow( qsqgev , 2 ) * 0.01553 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 2.857 + pow( qsqgev , 1 ) * -0.1983 + pow( qsqgev , 2 ) * 0.02672 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -17.23 + pow( qsqgev , 1 ) * 4.295 + pow( qsqgev , 2 ) * -0.4347 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 32.63 + pow( qsqgev , 1 ) * -9.796 + pow( qsqgev , 2 ) * 0.9914 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -26.36 + pow( qsqgev , 1 ) * 8.607 + pow( qsqgev , 2 ) * -0.8779 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 7.773 + pow( qsqgev , 1 ) * -2.667 + pow( qsqgev , 2 ) * 0.2743 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3093 + pow( qsqgev , 1 ) * -0.2095 + pow( qsqgev , 2 ) * 0.01976 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 2.115 + pow( qsqgev , 1 ) * 0.005766 + pow( qsqgev , 2 ) * -0.0116 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -14.76 + pow( qsqgev , 1 ) * 3.509 + pow( qsqgev , 2 ) * -0.3146 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 29.2 + pow( qsqgev , 1 ) * -8.624 + pow( qsqgev , 2 ) * 0.8251 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -24.15 + pow( qsqgev , 1 ) * 7.813 + pow( qsqgev , 2 ) * -0.7694 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 7.23 + pow( qsqgev , 1 ) * -2.464 + pow( qsqgev , 2 ) * 0.247 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3367 + pow( qsqgev , 1 ) * -0.2071 + pow( qsqgev , 2 ) * 0.02118 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 1.489 + pow( qsqgev , 1 ) * 0.1608 + pow( qsqgev , 2 ) * -0.03574 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -12.43 + pow( qsqgev , 1 ) * 2.78 + pow( qsqgev , 2 ) * -0.2192 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 25.71 + pow( qsqgev , 1 ) * -7.414 + pow( qsqgev , 2 ) * 0.6738 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -21.75 + pow( qsqgev , 1 ) * 6.928 + pow( qsqgev , 2 ) * -0.6605 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 6.608 + pow( qsqgev , 1 ) * -2.223 + pow( qsqgev , 2 ) * 0.2175 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3518 + pow( qsqgev , 1 ) * -0.2008 + pow( qsqgev , 2 ) * 0.02135 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.9436 + pow( qsqgev , 1 ) * 0.292 + pow( qsqgev , 2 ) * -0.053 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -10.23 + pow( qsqgev , 1 ) * 2.077 + pow( qsqgev , 2 ) * -0.1357 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 22.2 + pow( qsqgev , 1 ) * -6.166 + pow( qsqgev , 2 ) * 0.5276 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -19.23 + pow( qsqgev , 1 ) * 5.969 + pow( qsqgev , 2 ) * -0.548 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 5.927 + pow( qsqgev , 1 ) * -1.952 + pow( qsqgev , 2 ) * 0.1856 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3576 + pow( qsqgev , 1 ) * -0.1923 + pow( qsqgev , 2 ) * 0.02081 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.4817 + pow( qsqgev , 1 ) * 0.3984 + pow( qsqgev , 2 ) * -0.06512 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -8.22 + pow( qsqgev , 1 ) * 1.437 + pow( qsqgev , 2 ) * -0.06484 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 18.85 + pow( qsqgev , 1 ) * -4.969 + pow( qsqgev , 2 ) * 0.3938 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -16.74 + pow( qsqgev , 1 ) * 5.019 + pow( qsqgev , 2 ) * -0.4406 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 5.235 + pow( qsqgev , 1 ) * -1.676 + pow( qsqgev , 2 ) * 0.1541 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3572 + pow( qsqgev , 1 ) * -0.183 + pow( qsqgev , 2 ) * 0.01995 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.0896 + pow( qsqgev , 1 ) * 0.4863 + pow( qsqgev , 2 ) * -0.07403 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -6.389 + pow( qsqgev , 1 ) * 0.8515 + pow( qsqgev , 2 ) * -0.002893 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 15.68 + pow( qsqgev , 1 ) * -3.828 + pow( qsqgev , 2 ) * 0.2698 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -14.31 + pow( qsqgev , 1 ) * 4.086 + pow( qsqgev , 2 ) * -0.3373 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 4.545 + pow( qsqgev , 1 ) * -1.4 + pow( qsqgev , 2 ) * 0.123 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.352 + pow( qsqgev , 1 ) * -0.1731 + pow( qsqgev , 2 ) * 0.01892 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2359 + pow( qsqgev , 1 ) * 0.5557 + pow( qsqgev , 2 ) * -0.08027 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -4.764 + pow( qsqgev , 1 ) * 0.3367 + pow( qsqgev , 2 ) * 0.04949 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 12.76 + pow( qsqgev , 1 ) * -2.786 + pow( qsqgev , 2 ) * 0.1594 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -12.02 + pow( qsqgev , 1 ) * 3.215 + pow( qsqgev , 2 ) * -0.2428 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 3.884 + pow( qsqgev , 1 ) * -1.138 + pow( qsqgev , 2 ) * 0.09414 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3431 + pow( qsqgev , 1 ) * -0.1631 + pow( qsqgev , 2 ) * 0.01779 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.4983 + pow( qsqgev , 1 ) * 0.6063 + pow( qsqgev , 2 ) * -0.0841 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -3.353 + pow( qsqgev , 1 ) * -0.09973 + pow( qsqgev , 2 ) * 0.09199 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 10.15 + pow( qsqgev , 1 ) * -1.867 + pow( qsqgev , 2 ) * 0.06472 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -9.926 + pow( qsqgev , 1 ) * 2.43 + pow( qsqgev , 2 ) * -0.1596 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 3.27 + pow( qsqgev , 1 ) * -0.898 + pow( qsqgev , 2 ) * 0.06828 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3316 + pow( qsqgev , 1 ) * -0.153 + pow( qsqgev , 2 ) * 0.01663 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7069 + pow( qsqgev , 1 ) * 0.6413 + pow( qsqgev , 2 ) * -0.08611 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -2.139 + pow( qsqgev , 1 ) * -0.4658 + pow( qsqgev , 2 ) * 0.1262 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 7.818 + pow( qsqgev , 1 ) * -1.062 + pow( qsqgev , 2 ) * -0.01609 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -8.023 + pow( qsqgev , 1 ) * 1.726 + pow( qsqgev , 2 ) * -0.08652 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 2.704 + pow( qsqgev , 1 ) * -0.6798 + pow( qsqgev , 2 ) * 0.04512 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 3  to  4  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3181 + pow( qsqgev , 1 ) * -0.143 + pow( qsqgev , 2 ) * 0.01547 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.8667 + pow( qsqgev , 1 ) * 0.6618 + pow( qsqgev , 2 ) * -0.08649 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.117 + pow( qsqgev , 1 ) * -0.7603 + pow( qsqgev , 2 ) * 0.1522 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 5.796 + pow( qsqgev , 1 ) * -0.3831 + pow( qsqgev , 2 ) * -0.0819 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -6.34 + pow( qsqgev , 1 ) * 1.12 + pow( qsqgev , 2 ) * -0.02526 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 2.197 + pow( qsqgev , 1 ) * -0.4892 + pow( qsqgev , 2 ) * 0.02541 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/



  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -11.24 + pow( qsqgev , 1 ) * 5.538 + pow( qsqgev , 2 ) * -0.7483 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 44.02 + pow( qsqgev , 1 ) * -20.98 + pow( qsqgev , 2 ) * 2.796 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -64.31 + pow( qsqgev , 1 ) * 29.84 + pow( qsqgev , 2 ) * -3.927 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 41.71 + pow( qsqgev , 1 ) * -18.85 + pow( qsqgev , 2 ) * 2.45 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -10.21 + pow( qsqgev , 1 ) * 4.485 + pow( qsqgev , 2 ) * -0.5748 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -2.627 + pow( qsqgev , 1 ) * 1.241 + pow( qsqgev , 2 ) * -0.1858 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 12.57 + pow( qsqgev , 1 ) * -5.335 + pow( qsqgev , 2 ) * 0.7442 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -21.55 + pow( qsqgev , 1 ) * 8.565 + pow( qsqgev , 2 ) * -1.129 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 16.05 + pow( qsqgev , 1 ) * -6.067 + pow( qsqgev , 2 ) * 0.7618 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -4.455 + pow( qsqgev , 1 ) * 1.619 + pow( qsqgev , 2 ) * -0.1942 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.7055 + pow( qsqgev , 1 ) * 0.2843 + pow( qsqgev , 2 ) * -0.05244 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 5.298 + pow( qsqgev , 1 ) * -1.765 + pow( qsqgev , 2 ) * 0.2487 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -11.43 + pow( qsqgev , 1 ) * 3.631 + pow( qsqgev , 2 ) * -0.4431 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 9.886 + pow( qsqgev , 1 ) * -3.079 + pow( qsqgev , 2 ) * 0.3437 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -3.066 + pow( qsqgev , 1 ) * 0.9465 + pow( qsqgev , 2 ) * -0.09926 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.155 + pow( qsqgev , 1 ) * 0.02308 + pow( qsqgev , 2 ) * -0.01311 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 3.007 + pow( qsqgev , 1 ) * -0.7293 + pow( qsqgev , 2 ) * 0.09719 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -8.037 + pow( qsqgev , 1 ) * 2.136 + pow( qsqgev , 2 ) * -0.2269 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 7.739 + pow( qsqgev , 1 ) * -2.148 + pow( qsqgev , 2 ) * 0.2091 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -2.566 + pow( qsqgev , 1 ) * 0.7328 + pow( qsqgev , 2 ) * -0.06823 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03852 + pow( qsqgev , 1 ) * -0.05691 + pow( qsqgev , 2 ) * 0.0002055 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 2.025 + pow( qsqgev , 1 ) * -0.3578 + pow( qsqgev , 2 ) * 0.04159 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -6.404 + pow( qsqgev , 1 ) * 1.538 + pow( qsqgev , 2 ) * -0.142 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 6.614 + pow( qsqgev , 1 ) * -1.744 + pow( qsqgev , 2 ) * 0.1533 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -2.284 + pow( qsqgev , 1 ) * 0.6335 + pow( qsqgev , 2 ) * -0.0547 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.121 + pow( qsqgev , 1 ) * -0.08264 + pow( qsqgev , 2 ) * 0.005119 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 1.458 + pow( qsqgev , 1 ) * -0.1871 + pow( qsqgev , 2 ) * 0.01718 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -5.313 + pow( qsqgev , 1 ) * 1.205 + pow( qsqgev , 2 ) * -0.09942 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 5.782 + pow( qsqgev , 1 ) * -1.488 + pow( qsqgev , 2 ) * 0.1223 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -2.057 + pow( qsqgev , 1 ) * 0.563 + pow( qsqgev , 2 ) * -0.04645 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1625 + pow( qsqgev , 1 ) * -0.09031 + pow( qsqgev , 2 ) * 0.006984 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 1.053 + pow( qsqgev , 1 ) * -0.08564 + pow( qsqgev , 2 ) * 0.004311 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -4.429 + pow( qsqgev , 1 ) * 0.9609 + pow( qsqgev , 2 ) * -0.07229 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 5.047 + pow( qsqgev , 1 ) * -1.275 + pow( qsqgev , 2 ) * 0.09979 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.842 + pow( qsqgev , 1 ) * 0.4981 + pow( qsqgev , 2 ) * -0.03978 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1853 + pow( qsqgev , 1 ) * -0.09124 + pow( qsqgev , 2 ) * 0.007608 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.7345 + pow( qsqgev , 1 ) * -0.01416 + pow( qsqgev , 2 ) * -0.003504 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -3.662 + pow( qsqgev , 1 ) * 0.7571 + pow( qsqgev , 2 ) * -0.05226 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 4.368 + pow( qsqgev , 1 ) * -1.08 + pow( qsqgev , 2 ) * 0.08123 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.633 + pow( qsqgev , 1 ) * 0.4354 + pow( qsqgev , 2 ) * -0.03385 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1988 + pow( qsqgev , 1 ) * -0.0898 + pow( qsqgev , 2 ) * 0.007731 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.4676 + pow( qsqgev , 1 ) * 0.04379 + pow( qsqgev , 2 ) * -0.009068 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -2.967 + pow( qsqgev , 1 ) * 0.5719 + pow( qsqgev , 2 ) * -0.03538 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 3.72 + pow( qsqgev , 1 ) * -0.8925 + pow( qsqgev , 2 ) * 0.0642 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.425 + pow( qsqgev , 1 ) * 0.372 + pow( qsqgev , 2 ) * -0.02809 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.206 + pow( qsqgev , 1 ) * -0.08724 + pow( qsqgev , 2 ) * 0.007612 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.2428 + pow( qsqgev , 1 ) * 0.09173 + pow( qsqgev , 2 ) * -0.01323 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -2.341 + pow( qsqgev , 1 ) * 0.4042 + pow( qsqgev , 2 ) * -0.02084 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 3.111 + pow( qsqgev , 1 ) * -0.7144 + pow( qsqgev , 2 ) * 0.04856 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.224 + pow( qsqgev , 1 ) * 0.31 + pow( qsqgev , 2 ) * -0.02258 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.209 + pow( qsqgev , 1 ) * -0.08412 + pow( qsqgev , 2 ) * 0.007373 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.05441 + pow( qsqgev , 1 ) * 0.1312 + pow( qsqgev , 2 ) * -0.0164 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.783 + pow( qsqgev , 1 ) * 0.2549 + pow( qsqgev , 2 ) * -0.008343 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 2.549 + pow( qsqgev , 1 ) * -0.55 + pow( qsqgev , 2 ) * 0.03448 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.034 + pow( qsqgev , 1 ) * 0.2514 + pow( qsqgev , 2 ) * -0.01748 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.209 + pow( qsqgev , 1 ) * -0.08076 + pow( qsqgev , 2 ) * 0.007078 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1037 + pow( qsqgev , 1 ) * 0.164 + pow( qsqgev , 2 ) * -0.01887 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.286 + pow( qsqgev , 1 ) * 0.1215 + pow( qsqgev , 2 ) * 0.002565 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 2.031 + pow( qsqgev , 1 ) * -0.398 + pow( qsqgev , 2 ) * 0.02164 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.8535 + pow( qsqgev , 1 ) * 0.1958 + pow( qsqgev , 2 ) * -0.0127 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2067 + pow( qsqgev , 1 ) * -0.07728 + pow( qsqgev , 2 ) * 0.006759 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2351 + pow( qsqgev , 1 ) * 0.1906 + pow( qsqgev , 2 ) * -0.02078 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.8475 + pow( qsqgev , 1 ) * 0.00428 + pow( qsqgev , 2 ) * 0.01197 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 1.561 + pow( qsqgev , 1 ) * -0.2602 + pow( qsqgev , 2 ) * 0.01016 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.6867 + pow( qsqgev , 1 ) * 0.1445 + pow( qsqgev , 2 ) * -0.00833 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2023 + pow( qsqgev , 1 ) * -0.07361 + pow( qsqgev , 2 ) * 0.006413 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3403 + pow( qsqgev , 1 ) * 0.2106 + pow( qsqgev , 2 ) * -0.02207 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.4719 + pow( qsqgev , 1 ) * -0.09402 + pow( qsqgev , 2 ) * 0.01961 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 1.145 + pow( qsqgev , 1 ) * -0.1405 + pow( qsqgev , 2 ) * 0.0004159 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.5366 + pow( qsqgev , 1 ) * 0.09909 + pow( qsqgev , 2 ) * -0.004538 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 4  to  5  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1968 + pow( qsqgev , 1 ) * -0.06998 + pow( qsqgev , 2 ) * 0.00607 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.4252 + pow( qsqgev , 1 ) * 0.2258 + pow( qsqgev , 2 ) * -0.02297 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.1472 + pow( qsqgev , 1 ) * -0.178 + pow( qsqgev , 2 ) * 0.026 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.7752 + pow( qsqgev , 1 ) * -0.03492 + pow( qsqgev , 2 ) * -0.008074 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.4005 + pow( qsqgev , 1 ) * 0.05825 + pow( qsqgev , 2 ) * -0.001164 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/


  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -13.77 + pow( qsqgev , 1 ) * 5.43 + pow( qsqgev , 2 ) * -0.5587 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 36.79 + pow( qsqgev , 1 ) * -14.31 + pow( qsqgev , 2 ) * 1.461 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -32.63 + pow( qsqgev , 1 ) * 12.57 + pow( qsqgev , 2 ) * -1.276 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 9.629 + pow( qsqgev , 1 ) * -3.675 + pow( qsqgev , 2 ) * 0.3711 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -2.217 + pow( qsqgev , 1 ) * 0.859 + pow( qsqgev , 2 ) * -0.09595 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 6.941 + pow( qsqgev , 1 ) * -2.516 + pow( qsqgev , 2 ) * 0.2675 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -6.977 + pow( qsqgev , 1 ) * 2.43 + pow( qsqgev , 2 ) * -0.2494 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 2.293 + pow( qsqgev , 1 ) * -0.7747 + pow( qsqgev , 2 ) * 0.07724 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.4722 + pow( qsqgev , 1 ) * 0.1637 + pow( qsqgev , 2 ) * -0.02166 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 2.291 + pow( qsqgev , 1 ) * -0.6792 + pow( qsqgev , 2 ) * 0.07262 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -2.9 + pow( qsqgev , 1 ) * 0.8258 + pow( qsqgev , 2 ) * -0.07947 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 1.111 + pow( qsqgev , 1 ) * -0.3103 + pow( qsqgev , 2 ) * 0.02805 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.07054 + pow( qsqgev , 1 ) * 0.009612 + pow( qsqgev , 2 ) * -0.004102 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 1.115 + pow( qsqgev , 1 ) * -0.2438 + pow( qsqgev , 2 ) * 0.02459 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.809 + pow( qsqgev , 1 ) * 0.4273 + pow( qsqgev , 2 ) * -0.03618 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.7831 + pow( qsqgev , 1 ) * -0.1913 + pow( qsqgev , 2 ) * 0.01521 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05291 + pow( qsqgev , 1 ) * -0.0321 + pow( qsqgev , 2 ) * 0.001038 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.6679 + pow( qsqgev , 1 ) * -0.1024 + pow( qsqgev , 2 ) * 0.008969 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.343 + pow( qsqgev , 1 ) * 0.2818 + pow( qsqgev , 2 ) * -0.0209 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.6312 + pow( qsqgev , 1 ) * -0.1442 + pow( qsqgev , 2 ) * 0.01039 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1011 + pow( qsqgev , 1 ) * -0.04455 + pow( qsqgev , 2 ) * 0.002737 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.4256 + pow( qsqgev , 1 ) * -0.0396 + pow( qsqgev , 2 ) * 0.002471 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.051 + pow( qsqgev , 1 ) * 0.2039 + pow( qsqgev , 2 ) * -0.0135 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.5266 + pow( qsqgev , 1 ) * -0.1157 + pow( qsqgev , 2 ) * 0.007792 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1225 + pow( qsqgev , 1 ) * -0.04766 + pow( qsqgev , 2 ) * 0.003277 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.2638 + pow( qsqgev , 1 ) * -0.004344 + pow( qsqgev , 2 ) * -0.0007704 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.83 + pow( qsqgev , 1 ) * 0.1505 + pow( qsqgev , 2 ) * -0.00901 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.4411 + pow( qsqgev , 1 ) * -0.09403 + pow( qsqgev , 2 ) * 0.006019 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1332 + pow( qsqgev , 1 ) * -0.04776 + pow( qsqgev , 2 ) * 0.00341 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1402 + pow( qsqgev , 1 ) * 0.02019 + pow( qsqgev , 2 ) * -0.002752 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.644 + pow( qsqgev , 1 ) * 0.1072 + pow( qsqgev , 2 ) * -0.005688 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.3651 + pow( qsqgev , 1 ) * -0.0751 + pow( qsqgev , 2 ) * 0.004578 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1385 + pow( qsqgev , 1 ) * -0.04672 + pow( qsqgev , 2 ) * 0.003381 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.04085 + pow( qsqgev , 1 ) * 0.03904 + pow( qsqgev , 2 ) * -0.004122 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.4822 + pow( qsqgev , 1 ) * 0.06994 + pow( qsqgev , 2 ) * -0.002991 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.2957 + pow( qsqgev , 1 ) * -0.0579 + pow( qsqgev , 2 ) * 0.003321 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1404 + pow( qsqgev , 1 ) * -0.04523 + pow( qsqgev , 2 ) * 0.003286 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.04042 + pow( qsqgev , 1 ) * 0.05414 + pow( qsqgev , 2 ) * -0.005137 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.3403 + pow( qsqgev , 1 ) * 0.0373 + pow( qsqgev , 2 ) * -0.0007094 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.2324 + pow( qsqgev , 1 ) * -0.04216 + pow( qsqgev , 2 ) * 0.002198 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1402 + pow( qsqgev , 1 ) * -0.04355 + pow( qsqgev , 2 ) * 0.003161 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.107 + pow( qsqgev , 1 ) * 0.06628 + pow( qsqgev , 2 ) * -0.005906 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.2163 + pow( qsqgev , 1 ) * 0.008886 + pow( qsqgev , 2 ) * 0.001229 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.1752 + pow( qsqgev , 1 ) * -0.02797 + pow( qsqgev , 2 ) * 0.001202 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1384 + pow( qsqgev , 1 ) * -0.04178 + pow( qsqgev , 2 ) * 0.003024 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1612 + pow( qsqgev , 1 ) * 0.07596 + pow( qsqgev , 2 ) * -0.006488 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.1085 + pow( qsqgev , 1 ) * -0.01573 + pow( qsqgev , 2 ) * 0.002875 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.1237 + pow( qsqgev , 1 ) * -0.01525 + pow( qsqgev , 2 ) * 0.0003214 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1356 + pow( qsqgev , 1 ) * -0.03997 + pow( qsqgev , 2 ) * 0.002883 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2049 + pow( qsqgev , 1 ) * 0.08353 + pow( qsqgev , 2 ) * -0.006919 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.01545 + pow( qsqgev , 1 ) * -0.03683 + pow( qsqgev , 2 ) * 0.004263 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.07791 + pow( qsqgev , 1 ) * -0.00399 + pow( qsqgev , 2 ) * -0.000449 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1319 + pow( qsqgev , 1 ) * -0.03815 + pow( qsqgev , 2 ) * 0.002741 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2393 + pow( qsqgev , 1 ) * 0.08912 + pow( qsqgev , 2 ) * -0.007211 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.06331 + pow( qsqgev , 1 ) * -0.05436 + pow( qsqgev , 2 ) * 0.005388 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.03798 + pow( qsqgev , 1 ) * 0.005684 + pow( qsqgev , 2 ) * -0.001099 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 5  to  6  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1278 + pow( qsqgev , 1 ) * -0.03635 + pow( qsqgev , 2 ) * 0.002601 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2661 + pow( qsqgev , 1 ) * 0.09319 + pow( qsqgev , 2 ) * -0.007403 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1302 + pow( qsqgev , 1 ) * -0.06902 + pow( qsqgev , 2 ) * 0.006311 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.003034 + pow( qsqgev , 1 ) * 0.01406 + pow( qsqgev , 2 ) * -0.001655 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/


  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -704.3 + pow( qsqgev , 1 ) * 225.9 + pow( qsqgev , 2 ) * -18.15 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 1788 + pow( qsqgev , 1 ) * -573.4 + pow( qsqgev , 2 ) * 46.07 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1513 + pow( qsqgev , 1 ) * 485.2 + pow( qsqgev , 2 ) * -38.98 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 426.7 + pow( qsqgev , 1 ) * -136.8 + pow( qsqgev , 2 ) * 11 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -24.95 + pow( qsqgev , 1 ) * 8.227 + pow( qsqgev , 2 ) * -0.6931 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 64.92 + pow( qsqgev , 1 ) * -21.31 + pow( qsqgev , 2 ) * 1.791 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -56.2 + pow( qsqgev , 1 ) * 18.38 + pow( qsqgev , 2 ) * -1.542 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 16.22 + pow( qsqgev , 1 ) * -5.287 + pow( qsqgev , 2 ) * 0.4429 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -3.358 + pow( qsqgev , 1 ) * 1.107 + pow( qsqgev , 2 ) * -0.09881 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 9.529 + pow( qsqgev , 1 ) * -3.044 + pow( qsqgev , 2 ) * 0.2659 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -8.856 + pow( qsqgev , 1 ) * 2.77 + pow( qsqgev , 2 ) * -0.2382 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 2.732 + pow( qsqgev , 1 ) * -0.8385 + pow( qsqgev , 2 ) * 0.0711 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.7278 + pow( qsqgev , 1 ) * 0.226 + pow( qsqgev , 2 ) * -0.02192 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 2.653 + pow( qsqgev , 1 ) * -0.7495 + pow( qsqgev , 2 ) * 0.06618 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -2.907 + pow( qsqgev , 1 ) * 0.7865 + pow( qsqgev , 2 ) * -0.06559 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 1.02 + pow( qsqgev , 1 ) * -0.2684 + pow( qsqgev , 2 ) * 0.02146 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1832 + pow( qsqgev , 1 ) * 0.0462 + pow( qsqgev , 2 ) * -0.005459 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 1.151 + pow( qsqgev , 1 ) * -0.2623 + pow( qsqgev , 2 ) * 0.02222 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.565 + pow( qsqgev , 1 ) * 0.3545 + pow( qsqgev , 2 ) * -0.02684 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.6256 + pow( qsqgev , 1 ) * -0.1418 + pow( qsqgev , 2 ) * 0.01014 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.02693 + pow( qsqgev , 1 ) * -0.001804 + pow( qsqgev , 2 ) * -0.0008344 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.6589 + pow( qsqgev , 1 ) * -0.1176 + pow( qsqgev , 2 ) * 0.009 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -1.091 + pow( qsqgev , 1 ) * 0.2171 + pow( qsqgev , 2 ) * -0.01461 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.4783 + pow( qsqgev , 1 ) * -0.0995 + pow( qsqgev , 2 ) * 0.006429 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0306 + pow( qsqgev , 1 ) * -0.01685 + pow( qsqgev , 2 ) * 0.0006922 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.4303 + pow( qsqgev , 1 ) * -0.06037 + pow( qsqgev , 2 ) * 0.003948 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.8437 + pow( qsqgev , 1 ) * 0.1554 + pow( qsqgev , 2 ) * -0.009445 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.395 + pow( qsqgev , 1 ) * -0.07868 + pow( qsqgev , 2 ) * 0.00474 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05606 + pow( qsqgev , 1 ) * -0.02184 + pow( qsqgev , 2 ) * 0.001232 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.2926 + pow( qsqgev , 1 ) * -0.03124 + pow( qsqgev , 2 ) * 0.001589 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.6757 + pow( qsqgev , 1 ) * 0.1182 + pow( qsqgev , 2 ) * -0.00665 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.334 + pow( qsqgev , 1 ) * -0.06484 + pow( qsqgev , 2 ) * 0.00373 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0695 + pow( qsqgev , 1 ) * -0.02354 + pow( qsqgev , 2 ) * 0.001431 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1926 + pow( qsqgev , 1 ) * -0.0122 + pow( qsqgev , 2 ) * 0.000209 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.5405 + pow( qsqgev , 1 ) * 0.08998 + pow( qsqgev , 2 ) * -0.004714 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.2817 + pow( qsqgev , 1 ) * -0.05337 + pow( qsqgev , 2 ) * 0.002957 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.07673 + pow( qsqgev , 1 ) * -0.02382 + pow( qsqgev , 2 ) * 0.001484 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1158 + pow( qsqgev , 1 ) * 0.001377 + pow( qsqgev , 2 ) * -0.0006792 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.4277 + pow( qsqgev , 1 ) * 0.06713 + pow( qsqgev , 2 ) * -0.003254 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.2358 + pow( qsqgev , 1 ) * -0.0435 + pow( qsqgev , 2 ) * 0.002327 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0808 + pow( qsqgev , 1 ) * -0.02358 + pow( qsqgev , 2 ) * 0.00148 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.05319 + pow( qsqgev , 1 ) * 0.01219 + pow( qsqgev , 2 ) * -0.001337 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.3289 + pow( qsqgev , 1 ) * 0.04718 + pow( qsqgev , 2 ) * -0.002029 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.1939 + pow( qsqgev , 1 ) * -0.03446 + pow( qsqgev , 2 ) * 0.001765 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.08263 + pow( qsqgev , 1 ) * -0.02301 + pow( qsqgev , 2 ) * 0.001445 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.002594 + pow( qsqgev , 1 ) * 0.02068 + pow( qsqgev , 2 ) * -0.00182 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.2436 + pow( qsqgev , 1 ) * 0.03011 + pow( qsqgev , 2 ) * -0.001015 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.1563 + pow( qsqgev , 1 ) * -0.02641 + pow( qsqgev , 2 ) * 0.001277 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.08311 + pow( qsqgev , 1 ) * -0.02232 + pow( qsqgev , 2 ) * 0.001397 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.03919 + pow( qsqgev , 1 ) * 0.0276 + pow( qsqgev , 2 ) * -0.002199 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.1689 + pow( qsqgev , 1 ) * 0.01514 + pow( qsqgev , 2 ) * -0.0001414 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.1223 + pow( qsqgev , 1 ) * -0.01912 + pow( qsqgev , 2 ) * 0.0008393 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.08261 + pow( qsqgev , 1 ) * -0.02155 + pow( qsqgev , 2 ) * 0.001345 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.07361 + pow( qsqgev , 1 ) * 0.03324 + pow( qsqgev , 2 ) * -0.002496 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.1035 + pow( qsqgev , 1 ) * 0.002081 + pow( qsqgev , 2 ) * 0.000611 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.09168 + pow( qsqgev , 1 ) * -0.01255 + pow( qsqgev , 2 ) * 0.000449 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 6  to  7  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0814 + pow( qsqgev , 1 ) * -0.02075 + pow( qsqgev , 2 ) * 0.001289 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1017 + pow( qsqgev , 1 ) * 0.03773 + pow( qsqgev , 2 ) * -0.002724 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.04688 + pow( qsqgev , 1 ) * -0.009171 + pow( qsqgev , 2 ) * 0.00125 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.06433 + pow( qsqgev , 1 ) * -0.006716 + pow( qsqgev , 2 ) * 0.000106 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/


  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( wgev < 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 239.8 + pow( qsqgev , 1 ) * -34.1 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -668.5 + pow( qsqgev , 1 ) * 95.19 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 698.7 + pow( qsqgev , 1 ) * -99.64 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -324.4 + pow( qsqgev , 1 ) * 46.33 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 56.43 + pow( qsqgev , 1 ) * -8.074 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 18.17 + pow( qsqgev , 1 ) * -2.754 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -51.43 + pow( qsqgev , 1 ) * 7.849 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 54.6 + pow( qsqgev , 1 ) * -8.391 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -25.69 + pow( qsqgev , 1 ) * 3.98 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 4.519 + pow( qsqgev , 1 ) * -0.7063 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 3.246 + pow( qsqgev , 1 ) * -0.5446 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -8.954 + pow( qsqgev , 1 ) * 1.545 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 9.278 + pow( qsqgev , 1 ) * -1.647 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.23 + pow( qsqgev , 1 ) * 0.7768 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.7119 + pow( qsqgev , 1 ) * -0.1364 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7931 + pow( qsqgev , 1 ) * -0.1581 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.897 + pow( qsqgev , 1 ) * 0.4258 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.659 + pow( qsqgev , 1 ) * -0.4313 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5799 + pow( qsqgev , 1 ) * 0.1905 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.05779 + pow( qsqgev , 1 ) * -0.03052 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2048 + pow( qsqgev , 1 ) * -0.05651 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2157 + pow( qsqgev , 1 ) * 0.1308 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.1555 + pow( qsqgev , 1 ) * -0.1084 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.2888 + pow( qsqgev , 1 ) * 0.03355 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.09774 + pow( qsqgev , 1 ) * -0.001953 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03703 + pow( qsqgev , 1 ) * -0.02302 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.2411 + pow( qsqgev , 1 ) * 0.03517 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.6341 + pow( qsqgev , 1 ) * -0.004229 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.5116 + pow( qsqgev , 1 ) * -0.017 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.1365 + pow( qsqgev , 1 ) * 0.00723 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.0119 + pow( qsqgev , 1 ) * -0.01037 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.3512 + pow( qsqgev , 1 ) * 0.000849 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.732 + pow( qsqgev , 1 ) * 0.032 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.5492 + pow( qsqgev , 1 ) * -0.03412 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.1416 + pow( qsqgev , 1 ) * 0.01027 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.02295 + pow( qsqgev , 1 ) * -0.005194 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.3519 + pow( qsqgev , 1 ) * -0.01149 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.7094 + pow( qsqgev , 1 ) * 0.04367 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.528 + pow( qsqgev , 1 ) * -0.03906 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.1359 + pow( qsqgev , 1 ) * 0.01104 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.02146 + pow( qsqgev , 1 ) * -0.002997 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.318 + pow( qsqgev , 1 ) * -0.01516 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.6489 + pow( qsqgev , 1 ) * 0.04562 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.4882 + pow( qsqgev , 1 ) * -0.03912 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.1268 + pow( qsqgev , 1 ) * 0.0109 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.01621 + pow( qsqgev , 1 ) * -0.00207 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.2753 + pow( qsqgev , 1 ) * -0.01527 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.5793 + pow( qsqgev , 1 ) * 0.04358 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.4439 + pow( qsqgev , 1 ) * -0.03716 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.1169 + pow( qsqgev , 1 ) * 0.01036 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.01033 + pow( qsqgev , 1 ) * -0.001708 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.2331 + pow( qsqgev , 1 ) * -0.0139 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.5107 + pow( qsqgev , 1 ) * 0.03991 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.4002 + pow( qsqgev , 1 ) * -0.03436 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.107 + pow( qsqgev , 1 ) * 0.00966 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.004953 + pow( qsqgev , 1 ) * -0.001596 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1946 + pow( qsqgev , 1 ) * -0.01202 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.4471 + pow( qsqgev , 1 ) * 0.03571 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.3589 + pow( qsqgev , 1 ) * -0.03128 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.09746 + pow( qsqgev , 1 ) * 0.008896 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.0003834 + pow( qsqgev , 1 ) * -0.001592 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1605 + pow( qsqgev , 1 ) * -0.01001 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.3892 + pow( qsqgev , 1 ) * 0.03145 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.3205 + pow( qsqgev , 1 ) * -0.02815 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.08845 + pow( qsqgev , 1 ) * 0.008117 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.003337 + pow( qsqgev , 1 ) * -0.001626 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1309 + pow( qsqgev , 1 ) * -0.00808 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.3372 + pow( qsqgev , 1 ) * 0.02735 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.2851 + pow( qsqgev , 1 ) * -0.02512 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.07998 + pow( qsqgev , 1 ) * 0.007352 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_LT for Qsq: 7  to  8  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.006297 + pow( qsqgev , 1 ) * -0.00167 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1053 + pow( qsqgev , 1 ) * -0.006303 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.2906 + pow( qsqgev , 1 ) * 0.02352 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.2526 + pow( qsqgev , 1 ) * -0.02223 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.07203 + pow( qsqgev , 1 ) * 0.006611 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/
  /*****************************************************************************************************************************************/


  return s;


}
