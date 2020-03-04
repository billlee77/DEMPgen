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

/* double fSigmaTT; */

double MySigmaTT( double qsqgev, double tgev, double wgev ) {


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
  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -6.039 + pow( qsqgev , 1 ) * 6.791 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 16.34 + pow( qsqgev , 1 ) * -27.83 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -14.97 + pow( qsqgev , 1 ) * 38.41 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 4.309 + pow( qsqgev , 1 ) * -23.07 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1517 + pow( qsqgev , 1 ) * 5.051 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -4.901 + pow( qsqgev , 1 ) * 5.281 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 14.87 + pow( qsqgev , 1 ) * -23.48 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -16.02 + pow( qsqgev , 1 ) * 34.69 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 6.776 + pow( qsqgev , 1 ) * -22.14 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.8159 + pow( qsqgev , 1 ) * 5.149 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -3.998 + pow( qsqgev , 1 ) * 4.139 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 13.28 + pow( qsqgev , 1 ) * -19.69 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -15.9 + pow( qsqgev , 1 ) * 30.66 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 7.96 + pow( qsqgev , 1 ) * -20.42 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.366 + pow( qsqgev , 1 ) * 4.934 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -3.28 + pow( qsqgev , 1 ) * 3.271 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 11.73 + pow( qsqgev , 1 ) * -16.47 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -15.16 + pow( qsqgev , 1 ) * 26.76 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 8.361 + pow( qsqgev , 1 ) * -18.4 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.651 + pow( qsqgev , 1 ) * 4.566 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -2.708 + pow( qsqgev , 1 ) * 2.606 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 10.29 + pow( qsqgev , 1 ) * -13.78 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -14.1 + pow( qsqgev , 1 ) * 23.19 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 8.296 + pow( qsqgev , 1 ) * -16.36 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.768 + pow( qsqgev , 1 ) * 4.14 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -2.249 + pow( qsqgev , 1 ) * 2.094 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 9.005 + pow( qsqgev , 1 ) * -11.55 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -12.92 + pow( qsqgev , 1 ) * 20.02 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 7.964 + pow( qsqgev , 1 ) * -14.41 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.783 + pow( qsqgev , 1 ) * 3.706 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -1.879 + pow( qsqgev , 1 ) * 1.697 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 7.868 + pow( qsqgev , 1 ) * -9.705 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -11.73 + pow( qsqgev , 1 ) * 17.26 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 7.489 + pow( qsqgev , 1 ) * -12.64 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.736 + pow( qsqgev , 1 ) * 3.292 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -1.579 + pow( qsqgev , 1 ) * 1.386 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 6.874 + pow( qsqgev , 1 ) * -8.182 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -10.58 + pow( qsqgev , 1 ) * 14.87 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 6.948 + pow( qsqgev , 1 ) * -11.06 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.653 + pow( qsqgev , 1 ) * 2.912 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -1.335 + pow( qsqgev , 1 ) * 1.14 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 6.011 + pow( qsqgev , 1 ) * -6.924 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -9.508 + pow( qsqgev , 1 ) * 12.83 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 6.39 + pow( qsqgev , 1 ) * -9.667 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.551 + pow( qsqgev , 1 ) * 2.569 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -1.135 + pow( qsqgev , 1 ) * 0.9452 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 5.264 + pow( qsqgev , 1 ) * -5.882 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -8.524 + pow( qsqgev , 1 ) * 11.09 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 5.841 + pow( qsqgev , 1 ) * -8.451 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.442 + pow( qsqgev , 1 ) * 2.264 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.9694 + pow( qsqgev , 1 ) * 0.789 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 4.617 + pow( qsqgev , 1 ) * -5.017 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -7.633 + pow( qsqgev , 1 ) * 9.606 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 5.318 + pow( qsqgev , 1 ) * -7.394 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.331 + pow( qsqgev , 1 ) * 1.995 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.8323 + pow( qsqgev , 1 ) * 0.6629 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 4.058 + pow( qsqgev , 1 ) * -4.296 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -6.832 + pow( qsqgev , 1 ) * 8.339 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 4.829 + pow( qsqgev , 1 ) * -6.478 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.223 + pow( qsqgev , 1 ) * 1.759 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.7179 + pow( qsqgev , 1 ) * 0.5605 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 3.574 + pow( qsqgev , 1 ) * -3.693 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -6.115 + pow( qsqgev , 1 ) * 7.256 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 4.377 + pow( qsqgev , 1 ) * -5.684 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.119 + pow( qsqgev , 1 ) * 1.552 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.6221 + pow( qsqgev , 1 ) * 0.4767 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 3.155 + pow( qsqgev , 1 ) * -3.186 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -5.477 + pow( qsqgev , 1 ) * 6.331 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 3.964 + pow( qsqgev , 1 ) * -4.997 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -1.023 + pow( qsqgev , 1 ) * 1.372 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 0.1  to  1  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 1 )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.5414 + pow( qsqgev , 1 ) * 0.4078 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 2.792 + pow( qsqgev , 1 ) * -2.76 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -4.91 + pow( qsqgev , 1 ) * 5.54 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 3.59 + pow( qsqgev , 1 ) * -4.403 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * -0.9335 + pow( qsqgev , 1 ) * 1.215 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3944 + pow( qsqgev , 1 ) * 0.3564 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -7.876 + pow( qsqgev , 1 ) * 3.071 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 18.85 + pow( qsqgev , 1 ) * -8.218 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -16.26 + pow( qsqgev , 1 ) * 7.355 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 4.712 + pow( qsqgev , 1 ) * -2.176 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.5575 + pow( qsqgev , 1 ) * 0.3911 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -5.523 + pow( qsqgev , 1 ) * 2.027 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 14.68 + pow( qsqgev , 1 ) * -6.212 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -13.15 + pow( qsqgev , 1 ) * 5.802 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 3.89 + pow( qsqgev , 1 ) * -1.754 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.6346 + pow( qsqgev , 1 ) * 0.3945 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.746 + pow( qsqgev , 1 ) * 1.26 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 11.25 + pow( qsqgev , 1 ) * -4.592 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -10.48 + pow( qsqgev , 1 ) * 4.493 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 3.162 + pow( qsqgev , 1 ) * -1.39 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.6576 + pow( qsqgev , 1 ) * 0.3797 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.421 + pow( qsqgev , 1 ) * 0.7051 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 8.473 + pow( qsqgev , 1 ) * -3.312 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -8.232 + pow( qsqgev , 1 ) * 3.417 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 2.536 + pow( qsqgev , 1 ) * -1.083 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.6473 + pow( qsqgev , 1 ) * 0.3549 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.446 + pow( qsqgev , 1 ) * 0.3119 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.261 + pow( qsqgev , 1 ) * -2.319 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -6.381 + pow( qsqgev , 1 ) * 2.552 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 2.01 + pow( qsqgev , 1 ) * -0.8308 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.6176 + pow( qsqgev , 1 ) * 0.3256 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7381 + pow( qsqgev , 1 ) * 0.03944 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.523 + pow( qsqgev , 1 ) * -1.56 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.879 + pow( qsqgev , 1 ) * 1.868 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.575 + pow( qsqgev , 1 ) * -0.6278 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.5773 + pow( qsqgev , 1 ) * 0.2949 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2327 + pow( qsqgev , 1 ) * -0.1439 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.174 + pow( qsqgev , 1 ) * -0.9903 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -3.677 + pow( qsqgev , 1 ) * 1.335 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.221 + pow( qsqgev , 1 ) * -0.4667 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.5324 + pow( qsqgev , 1 ) * 0.2648 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1207 + pow( qsqgev , 1 ) * -0.2624 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.138 + pow( qsqgev , 1 ) * -0.5684 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.724 + pow( qsqgev , 1 ) * 0.9251 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.935 + pow( qsqgev , 1 ) * -0.3405 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.4863 + pow( qsqgev , 1 ) * 0.2364 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.3614 + pow( qsqgev , 1 ) * -0.3344 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.351 + pow( qsqgev , 1 ) * -0.2612 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.976 + pow( qsqgev , 1 ) * 0.6142 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.7071 + pow( qsqgev , 1 ) * -0.2428 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.4413 + pow( qsqgev , 1 ) * 0.2102 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.5191 + pow( qsqgev , 1 ) * -0.3734 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.76 + pow( qsqgev , 1 ) * -0.04156 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.393 + pow( qsqgev , 1 ) * 0.3806 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.5264 + pow( qsqgev , 1 ) * -0.1679 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3987 + pow( qsqgev , 1 ) * 0.1865 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.6163 + pow( qsqgev , 1 ) * -0.3896 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3214 + pow( qsqgev , 1 ) * 0.1119 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.9429 + pow( qsqgev , 1 ) * 0.2075 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.3843 + pow( qsqgev , 1 ) * -0.1109 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3591 + pow( qsqgev , 1 ) * 0.1652 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.6701 + pow( qsqgev , 1 ) * -0.3903 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.0002326 + pow( qsqgev , 1 ) * 0.216 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5977 + pow( qsqgev , 1 ) * 0.08079 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.273 + pow( qsqgev , 1 ) * -0.06809 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3229 + pow( qsqgev , 1 ) * 0.1463 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.6929 + pow( qsqgev , 1 ) * -0.3807 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.2306 + pow( qsqgev , 1 ) * 0.2836 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3352 + pow( qsqgev , 1 ) * -0.01041 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1865 + pow( qsqgev , 1 ) * -0.03619 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.29 + pow( qsqgev , 1 ) * 0.1295 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.6941 + pow( qsqgev , 1 ) * -0.3648 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.3928 + pow( qsqgev , 1 ) * 0.3246 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1373 + pow( qsqgev , 1 ) * -0.07486 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1196 + pow( qsqgev , 1 ) * -0.01266 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 1  to  2  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 1 ) && ( qsqgev <= 2 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2604 + pow( qsqgev , 1 ) * 0.1148 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.6806 + pow( qsqgev , 1 ) * -0.345 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.5031 + pow( qsqgev , 1 ) * 0.3463 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.01007 + pow( qsqgev , 1 ) * -0.1191 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.06822 + pow( qsqgev , 1 ) * 0.004422 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3235 + pow( qsqgev , 1 ) * -0.009602 + pow( qsqgev , 2 ) * -0.01691 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -8.138 + pow( qsqgev , 1 ) * 4.406 + pow( qsqgev , 2 ) * -0.6401 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 16.02 + pow( qsqgev , 1 ) * -9.309 + pow( qsqgev , 2 ) * 1.427 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -12.66 + pow( qsqgev , 1 ) * 7.568 + pow( qsqgev , 2 ) * -1.185 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 3.472 + pow( qsqgev , 1 ) * -2.113 + pow( qsqgev , 2 ) * 0.3351 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04716 + pow( qsqgev , 1 ) * 0.1169 + pow( qsqgev , 2 ) * -0.03257 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -6.214 + pow( qsqgev , 1 ) * 3.276 + pow( qsqgev , 2 ) * -0.4659 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 13.37 + pow( qsqgev , 1 ) * -7.64 + pow( qsqgev , 2 ) * 1.156 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -10.99 + pow( qsqgev , 1 ) * 6.47 + pow( qsqgev , 2 ) * -1.001 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 3.088 + pow( qsqgev , 1 ) * -1.851 + pow( qsqgev , 2 ) * 0.2902 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1314 + pow( qsqgev , 1 ) * 0.1922 + pow( qsqgev , 2 ) * -0.04104 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -4.663 + pow( qsqgev , 1 ) * 2.379 + pow( qsqgev , 2 ) * -0.3291 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 10.99 + pow( qsqgev , 1 ) * -6.156 + pow( qsqgev , 2 ) * 0.918 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -9.373 + pow( qsqgev , 1 ) * 5.425 + pow( qsqgev , 2 ) * -0.8286 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 2.693 + pow( qsqgev , 1 ) * -1.588 + pow( qsqgev , 2 ) * 0.2458 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2425 + pow( qsqgev , 1 ) * 0.233 + pow( qsqgev , 2 ) * -0.0448 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.424 + pow( qsqgev , 1 ) * 1.675 + pow( qsqgev , 2 ) * -0.2231 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 8.887 + pow( qsqgev , 1 ) * -4.869 + pow( qsqgev , 2 ) * 0.7134 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -7.865 + pow( qsqgev , 1 ) * 4.466 + pow( qsqgev , 2 ) * -0.6722 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 2.308 + pow( qsqgev , 1 ) * -1.336 + pow( qsqgev , 2 ) * 0.204 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3069 + pow( qsqgev , 1 ) * 0.2509 + pow( qsqgev , 2 ) * -0.0455 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.446 + pow( qsqgev , 1 ) * 1.129 + pow( qsqgev , 2 ) * -0.1423 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 7.074 + pow( qsqgev , 1 ) * -3.776 + pow( qsqgev , 2 ) * 0.5419 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -6.5 + pow( qsqgev , 1 ) * 3.612 + pow( qsqgev , 2 ) * -0.5348 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.946 + pow( qsqgev , 1 ) * -1.105 + pow( qsqgev , 2 ) * 0.1661 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3396 + pow( qsqgev , 1 ) * 0.2539 + pow( qsqgev , 2 ) * -0.04432 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.681 + pow( qsqgev , 1 ) * 0.7127 + pow( qsqgev , 2 ) * -0.0816 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 5.539 + pow( qsqgev , 1 ) * -2.867 + pow( qsqgev , 2 ) * 0.4011 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -5.294 + pow( qsqgev , 1 ) * 2.872 + pow( qsqgev , 2 ) * -0.4172 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.618 + pow( qsqgev , 1 ) * -0.8981 + pow( qsqgev , 2 ) * 0.1327 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3509 + pow( qsqgev , 1 ) * 0.2475 + pow( qsqgev , 2 ) * -0.04201 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.092 + pow( qsqgev , 1 ) * 0.3994 + pow( qsqgev , 2 ) * -0.03691 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.259 + pow( qsqgev , 1 ) * -2.124 + pow( qsqgev , 2 ) * 0.2876 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.251 + pow( qsqgev , 1 ) * 2.244 + pow( qsqgev , 2 ) * -0.3188 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.327 + pow( qsqgev , 1 ) * -0.7188 + pow( qsqgev , 2 ) * 0.1042 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3482 + pow( qsqgev , 1 ) * 0.2353 + pow( qsqgev , 2 ) * -0.03909 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.6416 + pow( qsqgev , 1 ) * 0.1676 + pow( qsqgev , 2 ) * -0.004703 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.207 + pow( qsqgev , 1 ) * -1.526 + pow( qsqgev , 2 ) * 0.1979 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -3.365 + pow( qsqgev , 1 ) * 1.721 + pow( qsqgev , 2 ) * -0.2381 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.074 + pow( qsqgev , 1 ) * -0.5663 + pow( qsqgev , 2 ) * 0.08032 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3366 + pow( qsqgev , 1 ) * 0.22 + pow( qsqgev , 2 ) * -0.0359 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3032 + pow( qsqgev , 1 ) * -0.0002344 + pow( qsqgev , 2 ) * 0.01786 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.355 + pow( qsqgev , 1 ) * -1.053 + pow( qsqgev , 2 ) * 0.1281 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.623 + pow( qsqgev , 1 ) * 1.292 + pow( qsqgev , 2 ) * -0.1731 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.8587 + pow( qsqgev , 1 ) * -0.439 + pow( qsqgev , 2 ) * 0.06069 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.3195 + pow( qsqgev , 1 ) * 0.2033 + pow( qsqgev , 2 ) * -0.03267 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.05287 + pow( qsqgev , 1 ) * -0.1187 + pow( qsqgev , 2 ) * 0.0331 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.672 + pow( qsqgev , 1 ) * -0.6835 + pow( qsqgev , 2 ) * 0.07484 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.011 + pow( qsqgev , 1 ) * 0.9467 + pow( qsqgev , 2 ) * -0.1216 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.6777 + pow( qsqgev , 1 ) * -0.3343 + pow( qsqgev , 2 ) * 0.04483 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2995 + pow( qsqgev , 1 ) * 0.1862 + pow( qsqgev , 2 ) * -0.02954 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1286 + pow( qsqgev , 1 ) * -0.1995 + pow( qsqgev , 2 ) * 0.04283 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.132 + pow( qsqgev , 1 ) * -0.4001 + pow( qsqgev , 2 ) * 0.03496 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.511 + pow( qsqgev , 1 ) * 0.6717 + pow( qsqgev , 2 ) * -0.08141 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.5276 + pow( qsqgev , 1 ) * -0.2496 + pow( qsqgev , 2 ) * 0.03224 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2781 + pow( qsqgev , 1 ) * 0.1695 + pow( qsqgev , 2 ) * -0.02657 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.2568 + pow( qsqgev , 1 ) * -0.2518 + pow( qsqgev , 2 ) * 0.0485 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7095 + pow( qsqgev , 1 ) * -0.1859 + pow( qsqgev , 2 ) * 0.005676 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.108 + pow( qsqgev , 1 ) * 0.4557 + pow( qsqgev , 2 ) * -0.05054 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4044 + pow( qsqgev , 1 ) * -0.1819 + pow( qsqgev , 2 ) * 0.02235 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2565 + pow( qsqgev , 1 ) * 0.1536 + pow( qsqgev , 2 ) * -0.02382 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.344 + pow( qsqgev , 1 ) * -0.2827 + pow( qsqgev , 2 ) * 0.0512 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3837 + pow( qsqgev , 1 ) * -0.02725 + pow( qsqgev , 2 ) * -0.01521 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7858 + pow( qsqgev , 1 ) * 0.2885 + pow( qsqgev , 2 ) * -0.02727 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.3043 + pow( qsqgev , 1 ) * -0.1283 + pow( qsqgev , 2 ) * 0.01472 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2355 + pow( qsqgev , 1 ) * 0.1387 + pow( qsqgev , 2 ) * -0.0213 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.4003 + pow( qsqgev , 1 ) * -0.2983 + pow( qsqgev , 2 ) * 0.05182 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1353 + pow( qsqgev , 1 ) * 0.08804 + pow( qsqgev , 2 ) * -0.02971 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5307 + pow( qsqgev , 1 ) * 0.1605 + pow( qsqgev , 2 ) * -0.009957 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.2237 + pow( qsqgev , 1 ) * -0.08644 + pow( qsqgev , 2 ) * 0.008893 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 2  to  3  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.2153 + pow( qsqgev , 1 ) * 0.125 + pow( qsqgev , 2 ) * -0.01902 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.4332 + pow( qsqgev , 1 ) * -0.3026 + pow( qsqgev , 2 ) * 0.05099 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.05089 + pow( qsqgev , 1 ) * 0.1693 + pow( qsqgev , 2 ) * -0.03928 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3309 + pow( qsqgev , 1 ) * 0.06411 + pow( qsqgev , 2 ) * 0.002616 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1592 + pow( qsqgev , 1 ) * -0.05411 + pow( qsqgev , 2 ) * 0.004522 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.605 + pow( qsqgev , 1 ) * -0.2134 + pow( qsqgev , 2 ) * 0.0199 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -4.491 + pow( qsqgev , 1 ) * 1.915 + pow( qsqgev , 2 ) * -0.2134 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 7.076 + pow( qsqgev , 1 ) * -3.144 + pow( qsqgev , 2 ) * 0.3622 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.985 + pow( qsqgev , 1 ) * 2.259 + pow( qsqgev , 2 ) * -0.2643 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.259 + pow( qsqgev , 1 ) * -0.5786 + pow( qsqgev , 2 ) * 0.06845 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3895 + pow( qsqgev , 1 ) * -0.1267 + pow( qsqgev , 2 ) * 0.01076 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.664 + pow( qsqgev , 1 ) * 1.542 + pow( qsqgev , 2 ) * -0.1702 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.281 + pow( qsqgev , 1 ) * -2.764 + pow( qsqgev , 2 ) * 0.3162 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.621 + pow( qsqgev , 1 ) * 2.075 + pow( qsqgev , 2 ) * -0.2412 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.204 + pow( qsqgev , 1 ) * -0.5483 + pow( qsqgev , 2 ) * 0.06441 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2338 + pow( qsqgev , 1 ) * -0.06543 + pow( qsqgev , 2 ) * 0.004426 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.959 + pow( qsqgev , 1 ) * 1.228 + pow( qsqgev , 2 ) * -0.1341 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 5.503 + pow( qsqgev , 1 ) * -2.397 + pow( qsqgev , 2 ) * 0.2721 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -4.215 + pow( qsqgev , 1 ) * 1.874 + pow( qsqgev , 2 ) * -0.2163 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.129 + pow( qsqgev , 1 ) * -0.5093 + pow( qsqgev , 2 ) * 0.05939 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1217 + pow( qsqgev , 1 ) * -0.02243 + pow( qsqgev , 2 ) * 6.633e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.362 + pow( qsqgev , 1 ) * 0.9642 + pow( qsqgev , 2 ) * -0.104 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.756 + pow( qsqgev , 1 ) * -2.048 + pow( qsqgev , 2 ) * 0.2306 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -3.783 + pow( qsqgev , 1 ) * 1.665 + pow( qsqgev , 2 ) * -0.1907 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 1.039 + pow( qsqgev , 1 ) * -0.4639 + pow( qsqgev , 2 ) * 0.05367 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04197 + pow( qsqgev , 1 ) * 0.007147 + pow( qsqgev , 2 ) * -0.002843 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.861 + pow( qsqgev , 1 ) * 0.7456 + pow( qsqgev , 2 ) * -0.07922 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.059 + pow( qsqgev , 1 ) * -1.727 + pow( qsqgev , 2 ) * 0.1926 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -3.347 + pow( qsqgev , 1 ) * 1.457 + pow( qsqgev , 2 ) * -0.1654 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.9406 + pow( qsqgev , 1 ) * -0.4153 + pow( qsqgev , 2 ) * 0.04765 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.01399 + pow( qsqgev , 1 ) * 0.02705 + pow( qsqgev , 2 ) * -0.004723 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.444 + pow( qsqgev , 1 ) * 0.5657 + pow( qsqgev , 2 ) * -0.05901 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.42 + pow( qsqgev , 1 ) * -1.435 + pow( qsqgev , 2 ) * 0.1584 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.921 + pow( qsqgev , 1 ) * 1.256 + pow( qsqgev , 2 ) * -0.1413 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.8385 + pow( qsqgev , 1 ) * -0.3658 + pow( qsqgev , 2 ) * 0.0416 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.0524 + pow( qsqgev , 1 ) * 0.03994 + pow( qsqgev , 2 ) * -0.005866 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.099 + pow( qsqgev , 1 ) * 0.419 + pow( qsqgev , 2 ) * -0.04267 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.845 + pow( qsqgev , 1 ) * -1.176 + pow( qsqgev , 2 ) * 0.1282 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.516 + pow( qsqgev , 1 ) * 1.068 + pow( qsqgev , 2 ) * -0.1189 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.737 + pow( qsqgev , 1 ) * -0.3175 + pow( qsqgev , 2 ) * 0.03575 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.0779 + pow( qsqgev , 1 ) * 0.04777 + pow( qsqgev , 2 ) * -0.006485 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.8171 + pow( qsqgev , 1 ) * 0.3005 + pow( qsqgev , 2 ) * -0.02961 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.335 + pow( qsqgev , 1 ) * -0.9487 + pow( qsqgev , 2 ) * 0.1021 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.14 + pow( qsqgev , 1 ) * 0.8952 + pow( qsqgev , 2 ) * -0.09859 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.6393 + pow( qsqgev , 1 ) * -0.2717 + pow( qsqgev , 2 ) * 0.03028 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.094 + pow( qsqgev , 1 ) * 0.05201 + pow( qsqgev , 2 ) * -0.006743 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5876 + pow( qsqgev , 1 ) * 0.2054 + pow( qsqgev , 2 ) * -0.01926 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.888 + pow( qsqgev , 1 ) * -0.7522 + pow( qsqgev , 2 ) * 0.07963 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.797 + pow( qsqgev , 1 ) * 0.7399 + pow( qsqgev , 2 ) * -0.08047 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.5476 + pow( qsqgev , 1 ) * -0.2293 + pow( qsqgev , 2 ) * 0.02526 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1032 + pow( qsqgev , 1 ) * 0.05369 + pow( qsqgev , 2 ) * -0.006747 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.4031 + pow( qsqgev , 1 ) * 0.1303 + pow( qsqgev , 2 ) * -0.01118 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.504 + pow( qsqgev , 1 ) * -0.5852 + pow( qsqgev , 2 ) * 0.06077 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.491 + pow( qsqgev , 1 ) * 0.6031 + pow( qsqgev , 2 ) * -0.06469 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4636 + pow( qsqgev , 1 ) * -0.1911 + pow( qsqgev , 2 ) * 0.02079 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1075 + pow( qsqgev , 1 ) * 0.05362 + pow( qsqgev , 2 ) * -0.006586 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2555 + pow( qsqgev , 1 ) * 0.0713 + pow( qsqgev , 2 ) * -0.004928 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.175 + pow( qsqgev , 1 ) * -0.4444 + pow( qsqgev , 2 ) * 0.04502 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.22 + pow( qsqgev , 1 ) * 0.4838 + pow( qsqgev , 2 ) * -0.05105 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.3878 + pow( qsqgev , 1 ) * -0.157 + pow( qsqgev , 2 ) * 0.01684 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1083 + pow( qsqgev , 1 ) * 0.05233 + pow( qsqgev , 2 ) * -0.006313 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1393 + pow( qsqgev , 1 ) * 0.02585 + pow( qsqgev , 2 ) * -0.0002015 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.8981 + pow( qsqgev , 1 ) * -0.3278 + pow( qsqgev , 2 ) * 0.03214 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.9853 + pow( qsqgev , 1 ) * 0.3818 + pow( qsqgev , 2 ) * -0.03953 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.3207 + pow( qsqgev , 1 ) * -0.1273 + pow( qsqgev , 2 ) * 0.01345 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1067 + pow( qsqgev , 1 ) * 0.05027 + pow( qsqgev , 2 ) * -0.005976 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.04856 + pow( qsqgev , 1 ) * -0.008754 + pow( qsqgev , 2 ) * 0.003319 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.667 + pow( qsqgev , 1 ) * -0.2319 + pow( qsqgev , 2 ) * 0.02169 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7833 + pow( qsqgev , 1 ) * 0.2954 + pow( qsqgev , 2 ) * -0.02988 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.2619 + pow( qsqgev , 1 ) * -0.1017 + pow( qsqgev , 2 ) * 0.01056 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.1033 + pow( qsqgev , 1 ) * 0.04772 + pow( qsqgev , 2 ) * -0.005602 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.02124 + pow( qsqgev , 1 ) * -0.03455 + pow( qsqgev , 2 ) * 0.005871 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.4761 + pow( qsqgev , 1 ) * -0.1542 + pow( qsqgev , 2 ) * 0.01333 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.6115 + pow( qsqgev , 1 ) * 0.2231 + pow( qsqgev , 2 ) * -0.02191 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.2111 + pow( qsqgev , 1 ) * -0.07995 + pow( qsqgev , 2 ) * 0.008122 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 3  to  4  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.09894 + pow( qsqgev , 1 ) * 0.04489 + pow( qsqgev , 2 ) * -0.005213 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.07393 + pow( qsqgev , 1 ) * -0.05325 + pow( qsqgev , 2 ) * 0.007651 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3201 + pow( qsqgev , 1 ) * -0.09193 + pow( qsqgev , 2 ) * 0.006753 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.467 + pow( qsqgev , 1 ) * 0.1633 + pow( qsqgev , 2 ) * -0.01542 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1677 + pow( qsqgev , 1 ) * -0.06164 + pow( qsqgev , 2 ) * 0.006103 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.4969 + pow( qsqgev , 1 ) * -0.1601 + pow( qsqgev , 2 ) * 0.01333 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.363 + pow( qsqgev , 1 ) * 0.8296 + pow( qsqgev , 2 ) * -0.07485 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.161 + pow( qsqgev , 1 ) * -1.139 + pow( qsqgev , 2 ) * 0.1053 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.023 + pow( qsqgev , 1 ) * 0.739 + pow( qsqgev , 2 ) * -0.06911 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4723 + pow( qsqgev , 1 ) * -0.1746 + pow( qsqgev , 2 ) * 0.01648 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3582 + pow( qsqgev , 1 ) * -0.1124 + pow( qsqgev , 2 ) * 0.009139 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.999 + pow( qsqgev , 1 ) * 0.6946 + pow( qsqgev , 2 ) * -0.06221 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.908 + pow( qsqgev , 1 ) * -1.039 + pow( qsqgev , 2 ) * 0.09546 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.952 + pow( qsqgev , 1 ) * 0.7075 + pow( qsqgev , 2 ) * -0.06577 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4735 + pow( qsqgev , 1 ) * -0.1736 + pow( qsqgev , 2 ) * 0.0163 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2546 + pow( qsqgev , 1 ) * -0.07735 + pow( qsqgev , 2 ) * 0.006111 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.685 + pow( qsqgev , 1 ) * 0.5797 + pow( qsqgev , 2 ) * -0.05158 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.656 + pow( qsqgev , 1 ) * -0.9422 + pow( qsqgev , 2 ) * 0.08609 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.863 + pow( qsqgev , 1 ) * 0.6706 + pow( qsqgev , 2 ) * -0.06204 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4677 + pow( qsqgev , 1 ) * -0.1703 + pow( qsqgev , 2 ) * 0.0159 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1758 + pow( qsqgev , 1 ) * -0.05108 + pow( qsqgev , 2 ) * 0.00386 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.406 + pow( qsqgev , 1 ) * 0.4789 + pow( qsqgev , 2 ) * -0.0423 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.393 + pow( qsqgev , 1 ) * -0.8425 + pow( qsqgev , 2 ) * 0.07655 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.749 + pow( qsqgev , 1 ) * 0.6249 + pow( qsqgev , 2 ) * -0.0575 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4524 + pow( qsqgev , 1 ) * -0.1635 + pow( qsqgev , 2 ) * 0.01518 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1158 + pow( qsqgev , 1 ) * -0.0314 + pow( qsqgev , 2 ) * 0.002189 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.161 + pow( qsqgev , 1 ) * 0.391 + pow( qsqgev , 2 ) * -0.03424 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.129 + pow( qsqgev , 1 ) * -0.7433 + pow( qsqgev , 2 ) * 0.06711 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.616 + pow( qsqgev , 1 ) * 0.5731 + pow( qsqgev , 2 ) * -0.05243 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4297 + pow( qsqgev , 1 ) * -0.1541 + pow( qsqgev , 2 ) * 0.01422 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.07066 + pow( qsqgev , 1 ) * -0.01685 + pow( qsqgev , 2 ) * 0.0009727 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.9487 + pow( qsqgev , 1 ) * 0.3158 + pow( qsqgev , 2 ) * -0.02739 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.873 + pow( qsqgev , 1 ) * -0.6483 + pow( qsqgev , 2 ) * 0.05815 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.475 + pow( qsqgev , 1 ) * 0.5188 + pow( qsqgev , 2 ) * -0.04717 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.4019 + pow( qsqgev , 1 ) * -0.143 + pow( qsqgev , 2 ) * 0.01311 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03696 + pow( qsqgev , 1 ) * -0.006213 + pow( qsqgev , 2 ) * 9.968e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7669 + pow( qsqgev , 1 ) * 0.2519 + pow( qsqgev , 2 ) * -0.02161 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.631 + pow( qsqgev , 1 ) * -0.5593 + pow( qsqgev , 2 ) * 0.04982 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.33 + pow( qsqgev , 1 ) * 0.4638 + pow( qsqgev , 2 ) * -0.0419 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.3707 + pow( qsqgev , 1 ) * -0.1307 + pow( qsqgev , 2 ) * 0.01191 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0119 + pow( qsqgev , 1 ) * 0.001511 + pow( qsqgev , 2 ) * -0.0005217 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.6114 + pow( qsqgev , 1 ) * 0.1977 + pow( qsqgev , 2 ) * -0.01675 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.403 + pow( qsqgev , 1 ) * -0.4764 + pow( qsqgev , 2 ) * 0.04211 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.184 + pow( qsqgev , 1 ) * 0.4093 + pow( qsqgev , 2 ) * -0.03672 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.3371 + pow( qsqgev , 1 ) * -0.1179 + pow( qsqgev , 2 ) * 0.01067 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.006484 + pow( qsqgev , 1 ) * 0.007016 + pow( qsqgev , 2 ) * -0.0009525 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.4793 + pow( qsqgev , 1 ) * 0.1521 + pow( qsqgev , 2 ) * -0.01268 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.193 + pow( qsqgev , 1 ) * -0.4007 + pow( qsqgev , 2 ) * 0.03511 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.042 + pow( qsqgev , 1 ) * 0.3566 + pow( qsqgev , 2 ) * -0.03175 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.3027 + pow( qsqgev , 1 ) * -0.1048 + pow( qsqgev , 2 ) * 0.009416 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.01973 + pow( qsqgev , 1 ) * 0.01083 + pow( qsqgev , 2 ) * -0.00124 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.368 + pow( qsqgev , 1 ) * 0.1141 + pow( qsqgev , 2 ) * -0.009323 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.003 + pow( qsqgev , 1 ) * -0.3326 + pow( qsqgev , 2 ) * 0.02886 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.9062 + pow( qsqgev , 1 ) * 0.3069 + pow( qsqgev , 2 ) * -0.02711 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.2685 + pow( qsqgev , 1 ) * -0.09204 + pow( qsqgev , 2 ) * 0.008204 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.02904 + pow( qsqgev , 1 ) * 0.01337 + pow( qsqgev , 2 ) * -0.00142 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2751 + pow( qsqgev , 1 ) * 0.08272 + pow( qsqgev , 2 ) * -0.006567 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.8318 + pow( qsqgev , 1 ) * -0.2721 + pow( qsqgev , 2 ) * 0.02335 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7794 + pow( qsqgev , 1 ) * 0.261 + pow( qsqgev , 2 ) * -0.02284 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.2355 + pow( qsqgev , 1 ) * -0.07986 + pow( qsqgev , 2 ) * 0.007059 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.03535 + pow( qsqgev , 1 ) * 0.01496 + pow( qsqgev , 2 ) * -0.001521 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1978 + pow( qsqgev , 1 ) * 0.05692 + pow( qsqgev , 2 ) * -0.004324 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.68 + pow( qsqgev , 1 ) * -0.2189 + pow( qsqgev , 2 ) * 0.01854 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.6625 + pow( qsqgev , 1 ) * 0.219 + pow( qsqgev , 2 ) * -0.01898 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.2043 + pow( qsqgev , 1 ) * -0.06844 + pow( qsqgev , 2 ) * 0.005995 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.03938 + pow( qsqgev , 1 ) * 0.01583 + pow( qsqgev , 2 ) * -0.001564 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1342 + pow( qsqgev , 1 ) * 0.03597 + pow( qsqgev , 2 ) * -0.002521 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.547 + pow( qsqgev , 1 ) * -0.1728 + pow( qsqgev , 2 ) * 0.01441 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5565 + pow( qsqgev , 1 ) * 0.1814 + pow( qsqgev , 2 ) * -0.01555 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1752 + pow( qsqgev , 1 ) * -0.05797 + pow( qsqgev , 2 ) * 0.005028 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.04169 + pow( qsqgev , 1 ) * 0.01619 + pow( qsqgev , 2 ) * -0.001566 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.08223 + pow( qsqgev , 1 ) * 0.01908 + pow( qsqgev , 2 ) * -0.001083 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.4312 + pow( qsqgev , 1 ) * -0.1332 + pow( qsqgev , 2 ) * 0.01088 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.4613 + pow( qsqgev , 1 ) * 0.148 + pow( qsqgev , 2 ) * -0.01252 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1486 + pow( qsqgev , 1 ) * -0.04848 + pow( qsqgev , 2 ) * 0.004159 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 4  to  5  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.04272 + pow( qsqgev , 1 ) * 0.01615 + pow( qsqgev , 2 ) * -0.001537 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.04017 + pow( qsqgev , 1 ) * 0.005626 + pow( qsqgev , 2 ) * 4.61e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3317 + pow( qsqgev , 1 ) * -0.09944 + pow( qsqgev , 2 ) * 0.007911 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3769 + pow( qsqgev , 1 ) * 0.1188 + pow( qsqgev , 2 ) * -0.009895 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1245 + pow( qsqgev , 1 ) * -0.04001 + pow( qsqgev , 2 ) * 0.00339 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3448 + pow( qsqgev , 1 ) * -0.09885 + pow( qsqgev , 2 ) * 0.007154 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.23 + pow( qsqgev , 1 ) * 0.37 + pow( qsqgev , 2 ) * -0.02817 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.443 + pow( qsqgev , 1 ) * -0.4403 + pow( qsqgev , 2 ) * 0.03409 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.8503 + pow( qsqgev , 1 ) * 0.2612 + pow( qsqgev , 2 ) * -0.0204 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1841 + pow( qsqgev , 1 ) * -0.05703 + pow( qsqgev , 2 ) * 0.004489 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2613 + pow( qsqgev , 1 ) * -0.07348 + pow( qsqgev , 2 ) * 0.005238 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.072 + pow( qsqgev , 1 ) * 0.3185 + pow( qsqgev , 2 ) * -0.02404 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.369 + pow( qsqgev , 1 ) * -0.4136 + pow( qsqgev , 2 ) * 0.0318 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.8506 + pow( qsqgev , 1 ) * 0.2591 + pow( qsqgev , 2 ) * -0.0201 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.193 + pow( qsqgev , 1 ) * -0.05932 + pow( qsqgev , 2 ) * 0.004639 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1989 + pow( qsqgev , 1 ) * -0.05508 + pow( qsqgev , 2 ) * 0.003884 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.9352 + pow( qsqgev , 1 ) * 0.2758 + pow( qsqgev , 2 ) * -0.02074 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.295 + pow( qsqgev , 1 ) * -0.3889 + pow( qsqgev , 2 ) * 0.0298 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.8429 + pow( qsqgev , 1 ) * 0.2554 + pow( qsqgev , 2 ) * -0.01975 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.199 + pow( qsqgev , 1 ) * -0.06085 + pow( qsqgev , 2 ) * 0.004742 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1483 + pow( qsqgev , 1 ) * -0.04023 + pow( qsqgev , 2 ) * 0.002787 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.8033 + pow( qsqgev , 1 ) * 0.2348 + pow( qsqgev , 2 ) * -0.01754 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.2 + pow( qsqgev , 1 ) * -0.3579 + pow( qsqgev , 2 ) * 0.02728 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.8163 + pow( qsqgev , 1 ) * 0.2457 + pow( qsqgev , 2 ) * -0.01891 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1998 + pow( qsqgev , 1 ) * -0.06066 + pow( qsqgev , 2 ) * 0.004704 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1089 + pow( qsqgev , 1 ) * -0.02885 + pow( qsqgev , 2 ) * 0.001956 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.6856 + pow( qsqgev , 1 ) * 0.1987 + pow( qsqgev , 2 ) * -0.01475 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.102 + pow( qsqgev , 1 ) * -0.3265 + pow( qsqgev , 2 ) * 0.02477 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7808 + pow( qsqgev , 1 ) * 0.2336 + pow( qsqgev , 2 ) * -0.0179 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1972 + pow( qsqgev , 1 ) * -0.05951 + pow( qsqgev , 2 ) * 0.004593 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.07799 + pow( qsqgev , 1 ) * -0.02003 + pow( qsqgev , 2 ) * 0.001317 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5798 + pow( qsqgev , 1 ) * 0.1666 + pow( qsqgev , 2 ) * -0.01229 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.001 + pow( qsqgev , 1 ) * -0.2945 + pow( qsqgev , 2 ) * 0.02223 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7369 + pow( qsqgev , 1 ) * 0.2191 + pow( qsqgev , 2 ) * -0.0167 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1915 + pow( qsqgev , 1 ) * -0.05742 + pow( qsqgev , 2 ) * 0.00441 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05384 + pow( qsqgev , 1 ) * -0.01324 + pow( qsqgev , 2 ) * 0.00083 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.486 + pow( qsqgev , 1 ) * 0.1383 + pow( qsqgev , 2 ) * -0.01014 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.8998 + pow( qsqgev , 1 ) * -0.2629 + pow( qsqgev , 2 ) * 0.01975 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.687 + pow( qsqgev , 1 ) * 0.2029 + pow( qsqgev , 2 ) * -0.0154 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1832 + pow( qsqgev , 1 ) * -0.05457 + pow( qsqgev , 2 ) * 0.00417 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.035 + pow( qsqgev , 1 ) * -0.008 + pow( qsqgev , 2 ) * 0.0004581 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.403 + pow( qsqgev , 1 ) * 0.1136 + pow( qsqgev , 2 ) * -0.008257 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.8001 + pow( qsqgev , 1 ) * -0.2321 + pow( qsqgev , 2 ) * 0.01733 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.6326 + pow( qsqgev , 1 ) * 0.1856 + pow( qsqgev , 2 ) * -0.01401 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1727 + pow( qsqgev , 1 ) * -0.05111 + pow( qsqgev , 2 ) * 0.003886 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.02041 + pow( qsqgev , 1 ) * -0.00401 + pow( qsqgev , 2 ) * 0.0001784 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3304 + pow( qsqgev , 1 ) * 0.09213 + pow( qsqgev , 2 ) * -0.006639 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7042 + pow( qsqgev , 1 ) * -0.2027 + pow( qsqgev , 2 ) * 0.01505 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.576 + pow( qsqgev , 1 ) * 0.1678 + pow( qsqgev , 2 ) * -0.0126 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1608 + pow( qsqgev , 1 ) * -0.04724 + pow( qsqgev , 2 ) * 0.003572 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.009165 + pow( qsqgev , 1 ) * -0.0009856 + pow( qsqgev , 2 ) * -3.072e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2671 + pow( qsqgev , 1 ) * 0.07357 + pow( qsqgev , 2 ) * -0.005248 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.613 + pow( qsqgev , 1 ) * -0.175 + pow( qsqgev , 2 ) * 0.01291 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5186 + pow( qsqgev , 1 ) * 0.15 + pow( qsqgev , 2 ) * -0.01119 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1478 + pow( qsqgev , 1 ) * -0.0431 + pow( qsqgev , 2 ) * 0.003241 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0005928 + pow( qsqgev , 1 ) * 0.001274 + pow( qsqgev , 2 ) * -0.0001842 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2124 + pow( qsqgev , 1 ) * 0.05769 + pow( qsqgev , 2 ) * -0.004064 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.5279 + pow( qsqgev , 1 ) * -0.1494 + pow( qsqgev , 2 ) * 0.01094 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.4619 + pow( qsqgev , 1 ) * 0.1325 + pow( qsqgev , 2 ) * -0.00983 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1343 + pow( qsqgev , 1 ) * -0.03887 + pow( qsqgev , 2 ) * 0.002904 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.005862 + pow( qsqgev , 1 ) * 0.002936 + pow( qsqgev , 2 ) * -0.0002945 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1655 + pow( qsqgev , 1 ) * 0.04417 + pow( qsqgev , 2 ) * -0.003063 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.4493 + pow( qsqgev , 1 ) * -0.1259 + pow( qsqgev , 2 ) * 0.00915 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.407 + pow( qsqgev , 1 ) * 0.1158 + pow( qsqgev , 2 ) * -0.008531 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1207 + pow( qsqgev , 1 ) * -0.03465 + pow( qsqgev , 2 ) * 0.002572 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.01065 + pow( qsqgev , 1 ) * 0.004131 + pow( qsqgev , 2 ) * -0.0003715 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1255 + pow( qsqgev , 1 ) * 0.03272 + pow( qsqgev , 2 ) * -0.002221 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3775 + pow( qsqgev , 1 ) * -0.1046 + pow( qsqgev , 2 ) * 0.007536 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3546 + pow( qsqgev , 1 ) * 0.09993 + pow( qsqgev , 2 ) * -0.00731 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.1073 + pow( qsqgev , 1 ) * -0.03052 + pow( qsqgev , 2 ) * 0.00225 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.0141 + pow( qsqgev , 1 ) * 0.004958 + pow( qsqgev , 2 ) * -0.0004223 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.09167 + pow( qsqgev , 1 ) * 0.02316 + pow( qsqgev , 2 ) * -0.001523 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.313 + pow( qsqgev , 1 ) * -0.08568 + pow( qsqgev , 2 ) * 0.00611 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3058 + pow( qsqgev , 1 ) * 0.08532 + pow( qsqgev , 2 ) * -0.006193 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.09441 + pow( qsqgev , 1 ) * -0.0266 + pow( qsqgev , 2 ) * 0.001947 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 5  to  6  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.01653 + pow( qsqgev , 1 ) * 0.005506 + pow( qsqgev , 2 ) * -0.0004537 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.06316 + pow( qsqgev , 1 ) * 0.01517 + pow( qsqgev , 2 ) * -0.000945 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2552 + pow( qsqgev , 1 ) * -0.06886 + pow( qsqgev , 2 ) * 0.00485 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2605 + pow( qsqgev , 1 ) * 0.07188 + pow( qsqgev , 2 ) * -0.005171 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.08217 + pow( qsqgev , 1 ) * -0.02292 + pow( qsqgev , 2 ) * 0.001665 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2485 + pow( qsqgev , 1 ) * -0.066 + pow( qsqgev , 2 ) * 0.004348 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7201 + pow( qsqgev , 1 ) * 0.1961 + pow( qsqgev , 2 ) * -0.01335 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7707 + pow( qsqgev , 1 ) * -0.2111 + pow( qsqgev , 2 ) * 0.01454 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.423 + pow( qsqgev , 1 ) * 0.1156 + pow( qsqgev , 2 ) * -0.007983 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.08539 + pow( qsqgev , 1 ) * -0.02341 + pow( qsqgev , 2 ) * 0.001624 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1776 + pow( qsqgev , 1 ) * -0.04542 + pow( qsqgev , 2 ) * 0.002885 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5806 + pow( qsqgev , 1 ) * 0.153 + pow( qsqgev , 2 ) * -0.0101 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.6693 + pow( qsqgev , 1 ) * -0.1774 + pow( qsqgev , 2 ) * 0.01185 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3905 + pow( qsqgev , 1 ) * 0.1036 + pow( qsqgev , 2 ) * -0.006956 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.08342 + pow( qsqgev , 1 ) * -0.02225 + pow( qsqgev , 2 ) * 0.001503 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1378 + pow( qsqgev , 1 ) * -0.03462 + pow( qsqgev , 2 ) * 0.002169 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5096 + pow( qsqgev , 1 ) * 0.1325 + pow( qsqgev , 2 ) * -0.00866 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.6344 + pow( qsqgev , 1 ) * -0.1661 + pow( qsqgev , 2 ) * 0.011 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3884 + pow( qsqgev , 1 ) * 0.1019 + pow( qsqgev , 2 ) * -0.006784 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.08672 + pow( qsqgev , 1 ) * -0.02289 + pow( qsqgev , 2 ) * 0.001533 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1079 + pow( qsqgev , 1 ) * -0.02675 + pow( qsqgev , 2 ) * 0.00166 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.4517 + pow( qsqgev , 1 ) * 0.1164 + pow( qsqgev , 2 ) * -0.007569 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.6071 + pow( qsqgev , 1 ) * -0.1578 + pow( qsqgev , 2 ) * 0.0104 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3891 + pow( qsqgev , 1 ) * 0.1015 + pow( qsqgev , 2 ) * -0.006727 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.09044 + pow( qsqgev , 1 ) * -0.02373 + pow( qsqgev , 2 ) * 0.001582 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.08361 + pow( qsqgev , 1 ) * -0.02042 + pow( qsqgev , 2 ) * 0.001252 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3962 + pow( qsqgev , 1 ) * 0.1013 + pow( qsqgev , 2 ) * -0.006549 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.5728 + pow( qsqgev , 1 ) * -0.1479 + pow( qsqgev , 2 ) * 0.009709 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3831 + pow( qsqgev , 1 ) * 0.09933 + pow( qsqgev , 2 ) * -0.006558 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.09228 + pow( qsqgev , 1 ) * -0.02407 + pow( qsqgev , 2 ) * 0.001599 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.06371 + pow( qsqgev , 1 ) * -0.0153 + pow( qsqgev , 2 ) * 0.0009242 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3436 + pow( qsqgev , 1 ) * 0.08708 + pow( qsqgev , 2 ) * -0.005598 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.5329 + pow( qsqgev , 1 ) * -0.1367 + pow( qsqgev , 2 ) * 0.008927 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3709 + pow( qsqgev , 1 ) * 0.0956 + pow( qsqgev , 2 ) * -0.006283 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.09228 + pow( qsqgev , 1 ) * -0.02392 + pow( qsqgev , 2 ) * 0.001581 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04783 + pow( qsqgev , 1 ) * -0.01127 + pow( qsqgev , 2 ) * 0.0006692 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.296 + pow( qsqgev , 1 ) * 0.07445 + pow( qsqgev , 2 ) * -0.00476 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.4916 + pow( qsqgev , 1 ) * -0.1254 + pow( qsqgev , 2 ) * 0.00815 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3553 + pow( qsqgev , 1 ) * 0.09108 + pow( qsqgev , 2 ) * -0.005962 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.09098 + pow( qsqgev , 1 ) * -0.02346 + pow( qsqgev , 2 ) * 0.001544 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03503 + pow( qsqgev , 1 ) * -0.008062 + pow( qsqgev , 2 ) * 0.0004675 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2527 + pow( qsqgev , 1 ) * 0.06308 + pow( qsqgev , 2 ) * -0.00401 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.449 + pow( qsqgev , 1 ) * -0.1138 + pow( qsqgev , 2 ) * 0.007366 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3365 + pow( qsqgev , 1 ) * 0.08578 + pow( qsqgev , 2 ) * -0.005593 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.08848 + pow( qsqgev , 1 ) * -0.02269 + pow( qsqgev , 2 ) * 0.001487 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.02472 + pow( qsqgev , 1 ) * -0.005506 + pow( qsqgev , 2 ) * 0.000308 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2134 + pow( qsqgev , 1 ) * 0.05285 + pow( qsqgev , 2 ) * -0.003338 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.4054 + pow( qsqgev , 1 ) * -0.1021 + pow( qsqgev , 2 ) * 0.006578 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3147 + pow( qsqgev , 1 ) * 0.07975 + pow( qsqgev , 2 ) * -0.005176 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.08477 + pow( qsqgev , 1 ) * -0.02161 + pow( qsqgev , 2 ) * 0.00141 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.01651 + pow( qsqgev , 1 ) * -0.003496 + pow( qsqgev , 2 ) * 0.0001838 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1784 + pow( qsqgev , 1 ) * 0.04381 + pow( qsqgev , 2 ) * -0.002749 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3627 + pow( qsqgev , 1 ) * -0.09077 + pow( qsqgev , 2 ) * 0.005818 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2912 + pow( qsqgev , 1 ) * 0.07336 + pow( qsqgev , 2 ) * -0.00474 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.08024 + pow( qsqgev , 1 ) * -0.02033 + pow( qsqgev , 2 ) * 0.001321 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.009958 + pow( qsqgev , 1 ) * -0.00191 + pow( qsqgev , 2 ) * 8.672e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1471 + pow( qsqgev , 1 ) * 0.03581 + pow( qsqgev , 2 ) * -0.002229 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3209 + pow( qsqgev , 1 ) * -0.07977 + pow( qsqgev , 2 ) * 0.005085 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2664 + pow( qsqgev , 1 ) * 0.0667 + pow( qsqgev , 2 ) * -0.004288 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.075 + pow( qsqgev , 1 ) * -0.01888 + pow( qsqgev , 2 ) * 0.001221 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.004813 + pow( qsqgev , 1 ) * -0.0006849 + pow( qsqgev , 2 ) * 1.27e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1198 + pow( qsqgev , 1 ) * 0.02886 + pow( qsqgev , 2 ) * -0.001782 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2814 + pow( qsqgev , 1 ) * -0.06946 + pow( qsqgev , 2 ) * 0.004404 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2415 + pow( qsqgev , 1 ) * 0.06008 + pow( qsqgev , 2 ) * -0.003844 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.06941 + pow( qsqgev , 1 ) * -0.01737 + pow( qsqgev , 2 ) * 0.001117 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0007662 + pow( qsqgev , 1 ) * 0.0002657 + pow( qsqgev , 2 ) * -4.405e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.09573 + pow( qsqgev , 1 ) * 0.0228 + pow( qsqgev , 2 ) * -0.001392 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2439 + pow( qsqgev , 1 ) * -0.05974 + pow( qsqgev , 2 ) * 0.003765 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2166 + pow( qsqgev , 1 ) * 0.05351 + pow( qsqgev , 2 ) * -0.003405 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.0635 + pow( qsqgev , 1 ) * -0.01578 + pow( qsqgev , 2 ) * 0.00101 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.002376 + pow( qsqgev , 1 ) * 0.0009909 + pow( qsqgev , 2 ) * -8.665e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.07476 + pow( qsqgev , 1 ) * 0.01755 + pow( qsqgev , 2 ) * -0.001058 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2089 + pow( qsqgev , 1 ) * -0.05074 + pow( qsqgev , 2 ) * 0.003176 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1922 + pow( qsqgev , 1 ) * 0.04713 + pow( qsqgev , 2 ) * -0.002981 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.05747 + pow( qsqgev , 1 ) * -0.01418 + pow( qsqgev , 2 ) * 0.0009023 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 6  to  7  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * -0.004772 + pow( qsqgev , 1 ) * 0.001531 + pow( qsqgev , 2 ) * -0.0001177 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.05669 + pow( qsqgev , 1 ) * 0.01307 + pow( qsqgev , 2 ) * -0.0007738 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1767 + pow( qsqgev , 1 ) * -0.04254 + pow( qsqgev , 2 ) * 0.002643 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1689 + pow( qsqgev , 1 ) * 0.04109 + pow( qsqgev , 2 ) * -0.002583 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.05153 + pow( qsqgev , 1 ) * -0.01262 + pow( qsqgev , 2 ) * 0.0007985 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( wgev < 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.09313 + pow( qsqgev , 1 ) * -0.02089 + pow( qsqgev , 2 ) * 0.00107 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2065 + pow( qsqgev , 1 ) * 0.04657 + pow( qsqgev , 2 ) * -0.002443 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1962 + pow( qsqgev , 1 ) * -0.04316 + pow( qsqgev , 2 ) * 0.002252 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1135 + pow( qsqgev , 1 ) * 0.02491 + pow( qsqgev , 2 ) * -0.001335 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.02311 + pow( qsqgev , 1 ) * -0.005131 + pow( qsqgev , 2 ) * 0.0002816 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1171 + pow( qsqgev , 1 ) * -0.02781 + pow( qsqgev , 2 ) * 0.001602 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3198 + pow( qsqgev , 1 ) * 0.07705 + pow( qsqgev , 2 ) * -0.004565 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.341 + pow( qsqgev , 1 ) * -0.0817 + pow( qsqgev , 2 ) * 0.004871 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1895 + pow( qsqgev , 1 ) * 0.04498 + pow( qsqgev , 2 ) * -0.002681 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.03832 + pow( qsqgev , 1 ) * -0.0091 + pow( qsqgev , 2 ) * 0.0005442 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.08646 + pow( qsqgev , 1 ) * -0.01995 + pow( qsqgev , 2 ) * 0.00112 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2621 + pow( qsqgev , 1 ) * 0.06145 + pow( qsqgev , 2 ) * -0.003559 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.298 + pow( qsqgev , 1 ) * -0.06943 + pow( qsqgev , 2 ) * 0.004044 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1747 + pow( qsqgev , 1 ) * 0.04041 + pow( qsqgev , 2 ) * -0.002355 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.03717 + pow( qsqgev , 1 ) * -0.008612 + pow( qsqgev , 2 ) * 0.0005039 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0737 + pow( qsqgev , 1 ) * -0.0169 + pow( qsqgev , 2 ) * 0.000951 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2525 + pow( qsqgev , 1 ) * 0.05894 + pow( qsqgev , 2 ) * -0.003423 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3116 + pow( qsqgev , 1 ) * -0.07252 + pow( qsqgev , 2 ) * 0.004242 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1907 + pow( qsqgev , 1 ) * 0.04415 + pow( qsqgev , 2 ) * -0.002584 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.04234 + pow( qsqgev , 1 ) * -0.009824 + pow( qsqgev , 2 ) * 0.0005771 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05929 + pow( qsqgev , 1 ) * -0.01342 + pow( qsqgev , 2 ) * 0.0007489 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2277 + pow( qsqgev , 1 ) * 0.05263 + pow( qsqgev , 2 ) * -0.003039 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3024 + pow( qsqgev , 1 ) * -0.06983 + pow( qsqgev , 2 ) * 0.004065 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1934 + pow( qsqgev , 1 ) * 0.0445 + pow( qsqgev , 2 ) * -0.002595 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.0447 + pow( qsqgev , 1 ) * -0.01031 + pow( qsqgev , 2 ) * 0.0006034 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04668 + pow( qsqgev , 1 ) * -0.01041 + pow( qsqgev , 2 ) * 0.0005738 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2006 + pow( qsqgev , 1 ) * 0.04588 + pow( qsqgev , 2 ) * -0.002629 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2854 + pow( qsqgev , 1 ) * -0.0653 + pow( qsqgev , 2 ) * 0.003775 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1901 + pow( qsqgev , 1 ) * 0.0434 + pow( qsqgev , 2 ) * -0.002515 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.04551 + pow( qsqgev , 1 ) * -0.01042 + pow( qsqgev , 2 ) * 0.000606 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03643 + pow( qsqgev , 1 ) * -0.008009 + pow( qsqgev , 2 ) * 0.0004359 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1757 + pow( qsqgev , 1 ) * 0.03982 + pow( qsqgev , 2 ) * -0.002267 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2673 + pow( qsqgev , 1 ) * -0.0607 + pow( qsqgev , 2 ) * 0.003489 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.185 + pow( qsqgev , 1 ) * 0.04195 + pow( qsqgev , 2 ) * -0.002418 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.04571 + pow( qsqgev , 1 ) * -0.01039 + pow( qsqgev , 2 ) * 0.0006012 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.02825 + pow( qsqgev , 1 ) * -0.006127 + pow( qsqgev , 2 ) * 0.0003295 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1538 + pow( qsqgev , 1 ) * 0.03461 + pow( qsqgev , 2 ) * -0.00196 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.25 + pow( qsqgev , 1 ) * -0.05645 + pow( qsqgev , 2 ) * 0.003231 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1796 + pow( qsqgev , 1 ) * 0.04052 + pow( qsqgev , 2 ) * -0.002327 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.04568 + pow( qsqgev , 1 ) * -0.01034 + pow( qsqgev , 2 ) * 0.0005958 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.02145 + pow( qsqgev , 1 ) * -0.004579 + pow( qsqgev , 2 ) * 0.0002424 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1331 + pow( qsqgev , 1 ) * 0.02973 + pow( qsqgev , 2 ) * -0.001675 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2309 + pow( qsqgev , 1 ) * -0.05182 + pow( qsqgev , 2 ) * 0.002954 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1718 + pow( qsqgev , 1 ) * 0.03856 + pow( qsqgev , 2 ) * -0.002206 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.04487 + pow( qsqgev , 1 ) * -0.0101 + pow( qsqgev , 2 ) * 0.0005798 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.01587 + pow( qsqgev , 1 ) * -0.003326 + pow( qsqgev , 2 ) * 0.0001726 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1141 + pow( qsqgev , 1 ) * 0.02532 + pow( qsqgev , 2 ) * -0.001419 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2113 + pow( qsqgev , 1 ) * -0.04716 + pow( qsqgev , 2 ) * 0.002677 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1627 + pow( qsqgev , 1 ) * 0.03634 + pow( qsqgev , 2 ) * -0.002071 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.04355 + pow( qsqgev , 1 ) * -0.009754 + pow( qsqgev , 2 ) * 0.0005578 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0113 + pow( qsqgev , 1 ) * -0.00231 + pow( qsqgev , 2 ) * 0.0001164 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.09681 + pow( qsqgev , 1 ) * 0.02134 + pow( qsqgev , 2 ) * -0.001189 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1915 + pow( qsqgev , 1 ) * -0.04251 + pow( qsqgev , 2 ) * 0.002402 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1525 + pow( qsqgev , 1 ) * 0.03389 + pow( qsqgev , 2 ) * -0.001924 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.04177 + pow( qsqgev , 1 ) * -0.009308 + pow( qsqgev , 2 ) * 0.0005301 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.007565 + pow( qsqgev , 1 ) * -0.001488 + pow( qsqgev , 2 ) * 7.13e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.08109 + pow( qsqgev , 1 ) * 0.01775 + pow( qsqgev , 2 ) * -0.0009828 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1717 + pow( qsqgev , 1 ) * -0.03789 + pow( qsqgev , 2 ) * 0.002131 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1414 + pow( qsqgev , 1 ) * 0.03124 + pow( qsqgev , 2 ) * -0.001766 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.03957 + pow( qsqgev , 1 ) * -0.008769 + pow( qsqgev , 2 ) * 0.0004974 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.004541 + pow( qsqgev , 1 ) * -0.0008305 + pow( qsqgev , 2 ) * 3.556e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.06702 + pow( qsqgev , 1 ) * 0.01456 + pow( qsqgev , 2 ) * -0.0008009 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1525 + pow( qsqgev , 1 ) * -0.03345 + pow( qsqgev , 2 ) * 0.001873 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1298 + pow( qsqgev , 1 ) * 0.02852 + pow( qsqgev , 2 ) * -0.001605 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.03709 + pow( qsqgev , 1 ) * -0.008174 + pow( qsqgev , 2 ) * 0.0004616 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.002105 + pow( qsqgev , 1 ) * -0.000307 + pow( qsqgev , 2 ) * 7.391e-06 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.05448 + pow( qsqgev , 1 ) * 0.01174 + pow( qsqgev , 2 ) * -0.0006408 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.134 + pow( qsqgev , 1 ) * -0.02921 + pow( qsqgev , 2 ) * 0.001627 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.118 + pow( qsqgev , 1 ) * 0.02578 + pow( qsqgev , 2 ) * -0.001444 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.03441 + pow( qsqgev , 1 ) * -0.00754 + pow( qsqgev , 2 ) * 0.0004239 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_TT for Qsq: 7  to  8  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0001485 + pow( qsqgev , 1 ) * 0.0001082 + pow( qsqgev , 2 ) * -1.473e-05 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.04332 + pow( qsqgev , 1 ) * 0.009244 + pow( qsqgev , 2 ) * -0.0005 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1164 + pow( qsqgev , 1 ) * -0.02519 + pow( qsqgev , 2 ) * 0.001395 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1061 + pow( qsqgev , 1 ) * 0.02304 + pow( qsqgev , 2 ) * -0.001284 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0.03158 + pow( qsqgev , 1 ) * -0.006877 + pow( qsqgev , 2 ) * 0.0003848 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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
