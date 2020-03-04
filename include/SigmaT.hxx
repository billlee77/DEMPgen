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

/* double fSigmaT; */

double MySigmaT( double qsqgev, double tgev, double wgev ) { 


  double s =  0.0;
  /* double wgev   = par[0]; */
  /* double qsqgev = par[1]; */
  /* double tgev   = x[0]; */

   //---------------------------------------------------------------------------------------------//
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //---------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( wgev <= 2.1 )  { 

      s = ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 2.546 + pow( qsqgev , 1 ) * -1.074 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -45.04  + pow( qsqgev , 1 ) * 60.66 + pow( qsqgev , 2 ) * -22.69 ) ) + 
	  (  pow( qsqgev , 0 ) * 6.504 + pow( qsqgev , 1 ) * -7.199 + pow( qsqgev , 2 ) * 2.295 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -4.284  + pow( qsqgev , 1 ) * 5.635 + pow( qsqgev , 2 ) * -1.911 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 2.332 + pow( qsqgev , 1 ) * -1.054 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -40.39  + pow( qsqgev , 1 ) * 51.49 + pow( qsqgev , 2 ) * -18.97 ) ) + 
	  (  pow( qsqgev , 0 ) * 4.724 + pow( qsqgev , 1 ) * -5.13 + pow( qsqgev , 2 ) * 1.628 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -3.525  + pow( qsqgev , 1 ) * 4.489 + pow( qsqgev , 2 ) * -1.495 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 2.12 + pow( qsqgev , 1 ) * -1.023 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -34.61  + pow( qsqgev , 1 ) * 40.38 + pow( qsqgev , 2 ) * -14.55 ) ) + 
	  (  pow( qsqgev , 0 ) * 3.461 + pow( qsqgev , 1 ) * -3.685 + pow( qsqgev , 2 ) * 1.166 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -2.799  + pow( qsqgev , 1 ) * 3.458 + pow( qsqgev , 2 ) * -1.135 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 1.91 + pow( qsqgev , 1 ) * -0.9857 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -27.26  + pow( qsqgev , 1 ) * 26.45 + pow( qsqgev , 2 ) * -9.059 ) ) + 
	  (  pow( qsqgev , 0 ) * 2.538 + pow( qsqgev , 1 ) * -2.63 + pow( qsqgev , 2 ) * 0.8264 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -2.161  + pow( qsqgev , 1 ) * 2.574 + pow( qsqgev , 2 ) * -0.831 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 1.694 + pow( qsqgev , 1 ) * -0.9393 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -18.64  + pow( qsqgev , 1 ) * 10.18 + pow( qsqgev , 2 ) * -2.648 ) ) + 
	  (  pow( qsqgev , 0 ) * 1.771 + pow( qsqgev , 1 ) * -1.691 + pow( qsqgev , 2 ) * 0.5114 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -1.533  + pow( qsqgev , 1 ) * 1.671 + pow( qsqgev , 2 ) * -0.5125 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 1.521 + pow( qsqgev , 1 ) * -0.9198 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -17.9  + pow( qsqgev , 1 ) * 8.338 + pow( qsqgev , 2 ) * -1.914 ) ) + 
	  (  pow( qsqgev , 0 ) * 1.326 + pow( qsqgev , 1 ) * -1.224 + pow( qsqgev , 2 ) * 0.3666 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -1.184  + pow( qsqgev , 1 ) * 1.236 + pow( qsqgev , 2 ) * -0.3718 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 1.358 + pow( qsqgev , 1 ) * -0.9031 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -17.73  + pow( qsqgev , 1 ) * 7.539 + pow( qsqgev , 2 ) * -1.597 ) ) + 
	  (  pow( qsqgev , 0 ) * 1.016 + pow( qsqgev , 1 ) * -0.9168 + pow( qsqgev , 2 ) * 0.2741 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.9309  + pow( qsqgev , 1 ) * 0.9422 + pow( qsqgev , 2 ) * -0.28 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 1.202 + pow( qsqgev , 1 ) * -0.8882 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -17.69  + pow( qsqgev , 1 ) * 6.986 + pow( qsqgev , 2 ) * -1.374 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.7906 + pow( qsqgev , 1 ) * -0.7023 + pow( qsqgev , 2 ) * 0.2104 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.7403  + pow( qsqgev , 1 ) * 0.7301 + pow( qsqgev , 2 ) * -0.2154 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 1.052 + pow( qsqgev , 1 ) * -0.8745 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -17.73  + pow( qsqgev , 1 ) * 6.568 + pow( qsqgev , 2 ) * -1.2 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.6239 + pow( qsqgev , 1 ) * -0.5481 + pow( qsqgev , 2 ) * 0.165 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.5943  + pow( qsqgev , 1 ) * 0.5737 + pow( qsqgev , 2 ) * -0.1686 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 0.9084 + pow( qsqgev , 1 ) * -0.8628 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -17.94  + pow( qsqgev , 1 ) * 6.43 + pow( qsqgev , 2 ) * -1.13 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.4999 + pow( qsqgev , 1 ) * -0.4371 + pow( qsqgev , 2 ) * 0.1326 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.483  + pow( qsqgev , 1 ) * 0.4591 + pow( qsqgev , 2 ) * -0.1349 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 0.7705 + pow( qsqgev , 1 ) * -0.852 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -18.28  + pow( qsqgev , 1 ) * 6.522 + pow( qsqgev , 2 ) * -1.153 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.4062 + pow( qsqgev , 1 ) * -0.3561 + pow( qsqgev , 2 ) * 0.1092 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.3972  + pow( qsqgev , 1 ) * 0.3742 + pow( qsqgev , 2 ) * -0.1106 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 0.6384 + pow( qsqgev , 1 ) * -0.8428 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -18.69  + pow( qsqgev , 1 ) * 6.698 + pow( qsqgev , 2 ) * -1.192 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.3342 + pow( qsqgev , 1 ) * -0.2948 + pow( qsqgev , 2 ) * 0.09148 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.33  + pow( qsqgev , 1 ) * 0.3095 + pow( qsqgev , 2 ) * -0.09206 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 0.5101 + pow( qsqgev , 1 ) * -0.8344 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -19.11  + pow( qsqgev , 1 ) * 6.918 + pow( qsqgev , 2 ) * -1.25 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.2778 + pow( qsqgev , 1 ) * -0.2475 + pow( qsqgev , 2 ) * 0.07769 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.2765  + pow( qsqgev , 1 ) * 0.259 + pow( qsqgev , 2 ) * -0.07767 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 0.3852 + pow( qsqgev , 1 ) * -0.8252 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -19.58  + pow( qsqgev , 1 ) * 7.256 + pow( qsqgev , 2 ) * -1.365 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.2334 + pow( qsqgev , 1 ) * -0.2109 + pow( qsqgev , 2 ) * 0.06704 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.2338  + pow( qsqgev , 1 ) * 0.2199 + pow( qsqgev , 2 ) * -0.06668 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 0.1  to  2  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev <= 2  )  { 
    if ( wgev >= 3.4 )  { 

      s = ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
	    pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) * 
	( exp( (  pow( qsqgev , 0 ) * 0.2654 + pow( qsqgev , 1 ) * -0.8181 + pow( qsqgev , 2 ) * 0 ) + 
	       tgev * ( pow( qsqgev , 0 ) * -20.09  + pow( qsqgev , 1 ) * 7.614 + pow( qsqgev , 2 ) * -1.474 ) ) + 
	  (  pow( qsqgev , 0 ) * 0.1979 + pow( qsqgev , 1 ) * -0.1816 + pow( qsqgev , 2 ) * 0.05842 ) + 
	  ( tgev * ( pow( qsqgev , 0 ) * -0.1993  + pow( qsqgev , 1 ) * 0.1887 + pow( qsqgev , 2 ) * -0.05786 ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


   //---------------------------------------------------------------------------------------------//
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //---------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 3.089 + pow( qsqgev , 1 ) * -0.8796 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -6.194 + pow( qsqgev , 1 ) * 1.932 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.944 + pow( qsqgev , 1 ) * -2.2 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.724 + pow( qsqgev , 1 ) * 0.866 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 2.612 + pow( qsqgev , 1 ) * -0.7368 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -5.789 + pow( qsqgev , 1 ) * 1.774 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.629 + pow( qsqgev , 1 ) * -2.064 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.635 + pow( qsqgev , 1 ) * 0.8255 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 2.219 + pow( qsqgev , 1 ) * -0.6202 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -5.349 + pow( qsqgev , 1 ) * 1.616 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 6.25 + pow( qsqgev , 1 ) * -1.918 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.516 + pow( qsqgev , 1 ) * 0.778 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.893 + pow( qsqgev , 1 ) * -0.5244 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -4.903 + pow( qsqgev , 1 ) * 1.463 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 5.836 + pow( qsqgev , 1 ) * -1.768 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.377 + pow( qsqgev , 1 ) * 0.7259 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.621 + pow( qsqgev , 1 ) * -0.4451 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -4.468 + pow( qsqgev , 1 ) * 1.317 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 5.409 + pow( qsqgev , 1 ) * -1.62 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.224 + pow( qsqgev , 1 ) * 0.6717 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.394 + pow( qsqgev , 1 ) * -0.3794 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -4.055 + pow( qsqgev , 1 ) * 1.183 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.984 + pow( qsqgev , 1 ) * -1.476 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -2.067 + pow( qsqgev , 1 ) * 0.6174 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.202 + pow( qsqgev , 1 ) * -0.3246 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.67 + pow( qsqgev , 1 ) * 1.059 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.573 + pow( qsqgev , 1 ) * -1.341 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.91 + pow( qsqgev , 1 ) * 0.5647 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.041 + pow( qsqgev , 1 ) * -0.2787 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -3.315 + pow( qsqgev , 1 ) * 0.9472 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 4.182 + pow( qsqgev , 1 ) * -1.214 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.758 + pow( qsqgev , 1 ) * 0.5144 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.9044 + pow( qsqgev , 1 ) * -0.2402 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.991 + pow( qsqgev , 1 ) * 0.8465 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.817 + pow( qsqgev , 1 ) * -1.097 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.613 + pow( qsqgev , 1 ) * 0.4674 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7885 + pow( qsqgev , 1 ) * -0.2078 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.697 + pow( qsqgev , 1 ) * 0.7564 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.478 + pow( qsqgev , 1 ) * -0.9906 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.477 + pow( qsqgev , 1 ) * 0.4239 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6898 + pow( qsqgev , 1 ) * -0.1805 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.433 + pow( qsqgev , 1 ) * 0.6762 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 3.167 + pow( qsqgev , 1 ) * -0.8939 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.351 + pow( qsqgev , 1 ) * 0.3842 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6054 + pow( qsqgev , 1 ) * -0.1573 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.195 + pow( qsqgev , 1 ) * 0.6049 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.883 + pow( qsqgev , 1 ) * -0.8068 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.235 + pow( qsqgev , 1 ) * 0.348 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.533 + pow( qsqgev , 1 ) * -0.1375 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.982 + pow( qsqgev , 1 ) * 0.5417 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.625 + pow( qsqgev , 1 ) * -0.7284 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.128 + pow( qsqgev , 1 ) * 0.3153 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.4708 + pow( qsqgev , 1 ) * -0.1207 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.792 + pow( qsqgev , 1 ) * 0.4857 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.391 + pow( qsqgev , 1 ) * -0.6582 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -1.031 + pow( qsqgev , 1 ) * 0.2858 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 2  to  3  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 2 ) && ( qsqgev <= 3 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.4171 + pow( qsqgev , 1 ) * -0.1063 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.621 + pow( qsqgev , 1 ) * 0.4362 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.179 + pow( qsqgev , 1 ) * -0.5953 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.9425 + pow( qsqgev , 1 ) * 0.2592 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------

   //---------------------------------------------------------------------------------------------//
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //---------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.45 + pow( qsqgev , 1 ) * -0.3202 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.038 + pow( qsqgev , 1 ) * 0.5048 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.034 + pow( qsqgev , 1 ) * -0.5094 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7555 + pow( qsqgev , 1 ) * 0.1876 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.259 + pow( qsqgev , 1 ) * -0.2752 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.05 + pow( qsqgev , 1 ) * 0.4916 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.118 + pow( qsqgev , 1 ) * -0.5127 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7953 + pow( qsqgev , 1 ) * 0.1921 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 1.097 + pow( qsqgev , 1 ) * -0.2378 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -2.01 + pow( qsqgev , 1 ) * 0.472 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.142 + pow( qsqgev , 1 ) * -0.5072 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.8164 + pow( qsqgev , 1 ) * 0.1936 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.9593 + pow( qsqgev , 1 ) * -0.2063 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.939 + pow( qsqgev , 1 ) * 0.4481 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.122 + pow( qsqgev , 1 ) * -0.4946 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.8215 + pow( qsqgev , 1 ) * 0.1921 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.8411 + pow( qsqgev , 1 ) * -0.1796 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.849 + pow( qsqgev , 1 ) * 0.4218 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.072 + pow( qsqgev , 1 ) * -0.4768 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.8139 + pow( qsqgev , 1 ) * 0.1881 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.7392 + pow( qsqgev , 1 ) * -0.1568 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.748 + pow( qsqgev , 1 ) * 0.3944 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 2.001 + pow( qsqgev , 1 ) * -0.4555 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7962 + pow( qsqgev , 1 ) * 0.1822 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6512 + pow( qsqgev , 1 ) * -0.1372 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.641 + pow( qsqgev , 1 ) * 0.3667 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.915 + pow( qsqgev , 1 ) * -0.4318 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7707 + pow( qsqgev , 1 ) * 0.1747 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.5749 + pow( qsqgev , 1 ) * -0.1203 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.533 + pow( qsqgev , 1 ) * 0.3395 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.82 + pow( qsqgev , 1 ) * -0.4068 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7398 + pow( qsqgev , 1 ) * 0.1663 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.5085 + pow( qsqgev , 1 ) * -0.1058 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.426 + pow( qsqgev , 1 ) * 0.3132 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.719 + pow( qsqgev , 1 ) * -0.3813 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.7052 + pow( qsqgev , 1 ) * 0.1572 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.4508 + pow( qsqgev , 1 ) * -0.09318 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.323 + pow( qsqgev , 1 ) * 0.2881 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.617 + pow( qsqgev , 1 ) * -0.3558 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.6684 + pow( qsqgev , 1 ) * 0.1478 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.4004 + pow( qsqgev , 1 ) * -0.08226 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.224 + pow( qsqgev , 1 ) * 0.2645 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.516 + pow( qsqgev , 1 ) * -0.3309 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.6306 + pow( qsqgev , 1 ) * 0.1384 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3564 + pow( qsqgev , 1 ) * -0.07279 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.13 + pow( qsqgev , 1 ) * 0.2424 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.416 + pow( qsqgev , 1 ) * -0.3069 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5927 + pow( qsqgev , 1 ) * 0.1291 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3178 + pow( qsqgev , 1 ) * -0.06454 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -1.043 + pow( qsqgev , 1 ) * 0.222 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.32 + pow( qsqgev , 1 ) * -0.2841 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5554 + pow( qsqgev , 1 ) * 0.1201 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.284 + pow( qsqgev , 1 ) * -0.05735 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.9605 + pow( qsqgev , 1 ) * 0.2031 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.229 + pow( qsqgev , 1 ) * -0.2625 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.5193 + pow( qsqgev , 1 ) * 0.1115 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 3  to  4  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 3 ) && ( qsqgev <= 4 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2543 + pow( qsqgev , 1 ) * -0.05108 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.8844 + pow( qsqgev , 1 ) * 0.1857 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 1.142 + pow( qsqgev , 1 ) * -0.2424 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.4847 + pow( qsqgev , 1 ) * 0.1034 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


   //---------------------------------------------------------------------------------------------//
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //---------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6873 + pow( qsqgev , 1 ) * -0.1265 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5864 + pow( qsqgev , 1 ) * 0.1344 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.504 + pow( qsqgev , 1 ) * -0.1185 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.185 + pow( qsqgev , 1 ) * 0.04164 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.6116 + pow( qsqgev , 1 ) * -0.1108 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.6799 + pow( qsqgev , 1 ) * 0.1423 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.632 + pow( qsqgev , 1 ) * -0.1333 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2289 + pow( qsqgev , 1 ) * 0.04735 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.545 + pow( qsqgev , 1 ) * -0.09761 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7289 + pow( qsqgev , 1 ) * 0.1454 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7123 + pow( qsqgev , 1 ) * -0.1425 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2602 + pow( qsqgev , 1 ) * 0.05152 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.4866 + pow( qsqgev , 1 ) * -0.08635 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.75 + pow( qsqgev , 1 ) * 0.1453 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7611 + pow( qsqgev , 1 ) * -0.1475 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2822 + pow( qsqgev , 1 ) * 0.05443 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.4352 + pow( qsqgev , 1 ) * -0.07663 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7528 + pow( qsqgev , 1 ) * 0.1428 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7879 + pow( qsqgev , 1 ) * -0.1496 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.297 + pow( qsqgev , 1 ) * 0.05628 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3898 + pow( qsqgev , 1 ) * -0.06816 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.743 + pow( qsqgev , 1 ) * 0.1388 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7985 + pow( qsqgev , 1 ) * -0.1494 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3058 + pow( qsqgev , 1 ) * 0.05719 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3497 + pow( qsqgev , 1 ) * -0.06075 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7247 + pow( qsqgev , 1 ) * 0.1337 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7972 + pow( qsqgev , 1 ) * -0.1474 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3098 + pow( qsqgev , 1 ) * 0.05733 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3141 + pow( qsqgev , 1 ) * -0.05423 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.7004 + pow( qsqgev , 1 ) * 0.1279 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7867 + pow( qsqgev , 1 ) * -0.1441 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3097 + pow( qsqgev , 1 ) * 0.0568 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2824 + pow( qsqgev , 1 ) * -0.04848 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.6721 + pow( qsqgev , 1 ) * 0.1216 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7692 + pow( qsqgev , 1 ) * -0.1397 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3063 + pow( qsqgev , 1 ) * 0.05573 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2542 + pow( qsqgev , 1 ) * -0.04339 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.6411 + pow( qsqgev , 1 ) * 0.115 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7464 + pow( qsqgev , 1 ) * -0.1345 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.3003 + pow( qsqgev , 1 ) * 0.05424 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2291 + pow( qsqgev , 1 ) * -0.03889 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.6087 + pow( qsqgev , 1 ) * 0.1083 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.7199 + pow( qsqgev , 1 ) * -0.1288 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2922 + pow( qsqgev , 1 ) * 0.05241 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2067 + pow( qsqgev , 1 ) * -0.0349 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5756 + pow( qsqgev , 1 ) * 0.1017 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.6907 + pow( qsqgev , 1 ) * -0.1227 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2826 + pow( qsqgev , 1 ) * 0.05035 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1867 + pow( qsqgev , 1 ) * -0.03136 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5426 + pow( qsqgev , 1 ) * 0.09518 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.6597 + pow( qsqgev , 1 ) * -0.1164 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2719 + pow( qsqgev , 1 ) * 0.04813 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1688 + pow( qsqgev , 1 ) * -0.02821 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.5101 + pow( qsqgev , 1 ) * 0.08887 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.6278 + pow( qsqgev , 1 ) * -0.1101 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2604 + pow( qsqgev , 1 ) * 0.0458 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 4  to  5  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 4 ) && ( qsqgev <= 5 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1528 + pow( qsqgev , 1 ) * -0.02541 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.4784 + pow( qsqgev , 1 ) * 0.08282 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.5956 + pow( qsqgev , 1 ) * -0.1038 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.2484 + pow( qsqgev , 1 ) * 0.04344 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


   //---------------------------------------------------------------------------------------------//
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //---------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.3192 + pow( qsqgev , 1 ) * -0.05197 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.05751 + pow( qsqgev , 1 ) * 0.02695 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.004748 + pow( qsqgev , 1 ) * -0.017 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.01052 + pow( qsqgev , 1 ) * 0.006129 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2946 + pow( qsqgev , 1 ) * -0.04664 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1614 + pow( qsqgev , 1 ) * 0.03711 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1228 + pow( qsqgev , 1 ) * -0.02987 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.04678 + pow( qsqgev , 1 ) * 0.01031 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2697 + pow( qsqgev , 1 ) * -0.04193 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2249 + pow( qsqgev , 1 ) * 0.04321 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1995 + pow( qsqgev , 1 ) * -0.03831 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.0723 + pow( qsqgev , 1 ) * 0.01333 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2464 + pow( qsqgev , 1 ) * -0.03778 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2647 + pow( qsqgev , 1 ) * 0.04686 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.251 + pow( qsqgev , 1 ) * -0.044 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.09101 + pow( qsqgev , 1 ) * 0.01559 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2248 + pow( qsqgev , 1 ) * -0.0341 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2892 + pow( qsqgev , 1 ) * 0.04888 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.2861 + pow( qsqgev , 1 ) * -0.04782 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.105 + pow( qsqgev , 1 ) * 0.01729 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.2052 + pow( qsqgev , 1 ) * -0.03084 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3037 + pow( qsqgev , 1 ) * 0.0498 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3102 + pow( qsqgev , 1 ) * -0.05037 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1155 + pow( qsqgev , 1 ) * 0.01857 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1873 + pow( qsqgev , 1 ) * -0.02794 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3109 + pow( qsqgev , 1 ) * 0.04992 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3261 + pow( qsqgev , 1 ) * -0.05193 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1232 + pow( qsqgev , 1 ) * 0.01949 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1711 + pow( qsqgev , 1 ) * -0.02534 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.313 + pow( qsqgev , 1 ) * 0.04945 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3359 + pow( qsqgev , 1 ) * -0.05273 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1287 + pow( qsqgev , 1 ) * 0.02012 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1563 + pow( qsqgev , 1 ) * -0.023 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3111 + pow( qsqgev , 1 ) * 0.04851 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3408 + pow( qsqgev , 1 ) * -0.05291 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1323 + pow( qsqgev , 1 ) * 0.02048 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1428 + pow( qsqgev , 1 ) * -0.02089 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.3062 + pow( qsqgev , 1 ) * 0.04723 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3418 + pow( qsqgev , 1 ) * -0.05257 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1343 + pow( qsqgev , 1 ) * 0.02061 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1306 + pow( qsqgev , 1 ) * -0.01899 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.299 + pow( qsqgev , 1 ) * 0.04568 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3397 + pow( qsqgev , 1 ) * -0.05182 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1349 + pow( qsqgev , 1 ) * 0.02055 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1194 + pow( qsqgev , 1 ) * -0.01727 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2901 + pow( qsqgev , 1 ) * 0.04394 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3349 + pow( qsqgev , 1 ) * -0.05072 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1342 + pow( qsqgev , 1 ) * 0.02032 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1092 + pow( qsqgev , 1 ) * -0.01572 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2799 + pow( qsqgev , 1 ) * 0.04207 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.328 + pow( qsqgev , 1 ) * -0.04935 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1326 + pow( qsqgev , 1 ) * 0.01995 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.09998 + pow( qsqgev , 1 ) * -0.01431 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2689 + pow( qsqgev , 1 ) * 0.04012 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3196 + pow( qsqgev , 1 ) * -0.04777 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1302 + pow( qsqgev , 1 ) * 0.01947 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 5  to  6  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 5 ) && ( qsqgev <= 6 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.09153 + pow( qsqgev , 1 ) * -0.01304 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.2573 + pow( qsqgev , 1 ) * 0.03813 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.3098 + pow( qsqgev , 1 ) * -0.04603 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.1271 + pow( qsqgev , 1 ) * 0.01889 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


   //---------------------------------------------------------------------------------------------//
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //---------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( wgev <= 2.1 )  { 

      s = ( ( ( pow( ( pow( 2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1287 + pow( qsqgev , 1 ) * -0.01998 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.1505 + pow( qsqgev , 1 ) * -0.008036 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.1683 + pow( qsqgev , 1 ) * 0.01216 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.0459 + pow( qsqgev , 1 ) * -0.003386 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.1 ) && ( wgev <= 2.2 ) )  { 

      s = ( ( ( pow( ( pow( 2.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1341 + pow( qsqgev , 1 ) * -0.01966 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * 0.0373 + pow( qsqgev , 1 ) * 0.003599 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * -0.05297 + pow( qsqgev , 1 ) * -0.0001854 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * 0.01242 + pow( qsqgev , 1 ) * 0.0003032 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.2 ) && ( wgev <= 2.3 ) )  { 

      s = ( ( ( pow( ( pow( 2.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.129 + pow( qsqgev , 1 ) * -0.01827 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.02419 + pow( qsqgev , 1 ) * 0.009363 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.01269 + pow( qsqgev , 1 ) * -0.006786 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.007635 + pow( qsqgev , 1 ) * 0.002411 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.3 ) && ( wgev <= 2.4 ) )  { 

      s = ( ( ( pow( ( pow( 2.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1219 + pow( qsqgev , 1 ) * -0.01685 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.06392 + pow( qsqgev , 1 ) * 0.01302 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.05636 + pow( qsqgev , 1 ) * -0.01116 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.02181 + pow( qsqgev , 1 ) * 0.003912 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.4 ) && ( wgev <= 2.5 ) )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1143 + pow( qsqgev , 1 ) * -0.01551 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.09088 + pow( qsqgev , 1 ) * 0.01546 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.08695 + pow( qsqgev , 1 ) * -0.01424 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.0324 + pow( qsqgev , 1 ) * 0.005049 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.1066 + pow( qsqgev , 1 ) * -0.01428 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1095 + pow( qsqgev , 1 ) * 0.0171 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.109 + pow( qsqgev , 1 ) * -0.01646 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.04052 + pow( qsqgev , 1 ) * 0.005931 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.09933 + pow( qsqgev , 1 ) * -0.01315 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1226 + pow( qsqgev , 1 ) * 0.01821 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1253 + pow( qsqgev , 1 ) * -0.0181 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.04697 + pow( qsqgev , 1 ) * 0.006637 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.09237 + pow( qsqgev , 1 ) * -0.01211 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1316 + pow( qsqgev , 1 ) * 0.01891 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1374 + pow( qsqgev , 1 ) * -0.0193 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.05202 + pow( qsqgev , 1 ) * 0.00719 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.08583 + pow( qsqgev , 1 ) * -0.01116 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1375 + pow( qsqgev , 1 ) * 0.01931 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1464 + pow( qsqgev , 1 ) * -0.02018 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.05605 + pow( qsqgev , 1 ) * 0.007628 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.07969 + pow( qsqgev , 1 ) * -0.01028 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1411 + pow( qsqgev , 1 ) * 0.01945 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1528 + pow( qsqgev , 1 ) * -0.02076 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.05914 + pow( qsqgev , 1 ) * 0.007956 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.07395 + pow( qsqgev , 1 ) * -0.009479 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1427 + pow( qsqgev , 1 ) * 0.01939 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1572 + pow( qsqgev , 1 ) * -0.02112 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.06146 + pow( qsqgev , 1 ) * 0.008194 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.06859 + pow( qsqgev , 1 ) * -0.008737 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1428 + pow( qsqgev , 1 ) * 0.01917 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1598 + pow( qsqgev , 1 ) * -0.02127 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.06307 + pow( qsqgev , 1 ) * 0.008344 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.06359 + pow( qsqgev , 1 ) * -0.008053 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1417 + pow( qsqgev , 1 ) * 0.01883 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.161 + pow( qsqgev , 1 ) * -0.02125 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.06408 + pow( qsqgev , 1 ) * 0.008421 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05894 + pow( qsqgev , 1 ) * -0.007422 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1395 + pow( qsqgev , 1 ) * 0.01837 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1608 + pow( qsqgev , 1 ) * -0.02108 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.06455 + pow( qsqgev , 1 ) * 0.00843 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 6  to  7  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( ( qsqgev >= 6 ) && ( qsqgev <= 7 ) )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05461 + pow( qsqgev , 1 ) * -0.006841 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.1365 + pow( qsqgev , 1 ) * 0.01783 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.1596 + pow( qsqgev , 1 ) * -0.02078 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.06454 + pow( qsqgev , 1 ) * 0.008382 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


   //---------------------------------------------------------------------------------------------//
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //                                                                                             //
   //---------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------
  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 2.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( wgev < 2.5 )  { 

      s = ( ( ( pow( ( pow( 2.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05511 + pow( qsqgev , 1 ) * -0.007004 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.006249 + pow( qsqgev , 1 ) * 0.00326 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.007846 + pow( qsqgev , 1 ) * -0.002827 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.004861 + pow( qsqgev , 1 ) * 0.001075 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 2.5
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.5 ) && ( wgev <= 2.6 ) )  { 

      s = ( ( ( pow( ( pow( 2.5 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05332 + pow( qsqgev , 1 ) * -0.006608 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.0239 + pow( qsqgev , 1 ) * 0.004764 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.02592 + pow( qsqgev , 1 ) * -0.004478 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.01087 + pow( qsqgev , 1 ) * 0.001653 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 2.6
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.6 ) && ( wgev <= 2.7 ) )  { 

      s = ( ( ( pow( ( pow( 2.6 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.05109 + pow( qsqgev , 1 ) * -0.006211 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.03685 + pow( qsqgev , 1 ) * 0.005857 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.03938 + pow( qsqgev , 1 ) * -0.005712 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.01557 + pow( qsqgev , 1 ) * 0.002109 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 2.7
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.7 ) && ( wgev <= 2.8 ) )  { 

      s = ( ( ( pow( ( pow( 2.7 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04867 + pow( qsqgev , 1 ) * -0.005826 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.04658 + pow( qsqgev , 1 ) * 0.006668 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.04974 + pow( qsqgev , 1 ) * -0.006666 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.01935 + pow( qsqgev , 1 ) * 0.00248 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 2.8
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.8 ) && ( wgev <= 2.9 ) )  { 

      s = ( ( ( pow( ( pow( 2.8 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04618 + pow( qsqgev , 1 ) * -0.005458 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.05396 + pow( qsqgev , 1 ) * 0.00727 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.05787 + pow( qsqgev , 1 ) * -0.007417 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.02245 + pow( qsqgev , 1 ) * 0.002787 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 2.9
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 2.9 ) && ( wgev <= 3 ) )  { 

      s = ( ( ( pow( ( pow( 2.9 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.0437 + pow( qsqgev , 1 ) * -0.005109 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.05954 + pow( qsqgev , 1 ) * 0.007709 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.06433 + pow( qsqgev , 1 ) * -0.008012 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.02503 + pow( qsqgev , 1 ) * 0.003042 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3 ) && ( wgev <= 3.1 ) )  { 

      s = ( ( ( pow( ( pow( 3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.04126 + pow( qsqgev , 1 ) * -0.004778 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.06369 + pow( qsqgev , 1 ) * 0.008018 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.06947 + pow( qsqgev , 1 ) * -0.008481 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.02716 + pow( qsqgev , 1 ) * 0.003253 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 3.1
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.1 ) && ( wgev <= 3.2 ) )  { 

      s = ( ( ( pow( ( pow( 3.1 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03888 + pow( qsqgev , 1 ) * -0.004466 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.06666 + pow( qsqgev , 1 ) * 0.008217 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.07349 + pow( qsqgev , 1 ) * -0.008841 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.02891 + pow( qsqgev , 1 ) * 0.003424 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 3.2
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.2 ) && ( wgev <= 3.3 ) )  { 

      s = ( ( ( pow( ( pow( 3.2 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03659 + pow( qsqgev , 1 ) * -0.004172 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.06867 + pow( qsqgev , 1 ) * 0.008324 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.07659 + pow( qsqgev , 1 ) * -0.009106 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.03032 + pow( qsqgev , 1 ) * 0.00356 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 3.3
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( ( wgev >= 3.3 ) && ( wgev <= 3.4 ) )  { 

      s = ( ( ( pow( ( pow( 3.3 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03439 + pow( qsqgev , 1 ) * -0.003895 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.06986 + pow( qsqgev , 1 ) * 0.008353 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.07887 + pow( qsqgev , 1 ) * -0.009287 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.03142 + pow( qsqgev , 1 ) * 0.003663 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 4 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 5 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 6 ) * (  pow( qsqgev , 0 ) * 0 + pow( qsqgev , 1 ) * 0 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) ) );

    }
  }
  //---------------------------------------------------------------------------------------------


  //Paramaterization fo Sigma_T for Qsq: 7  to  8  at w = 3.4
  //---------------------------------------------------------------------------------------------
  if ( qsqgev >= 7 )  { 
    if ( wgev >= 3.4 )  { 

      s = ( ( ( pow( ( pow( 3.4 , 2 ) - pow( 0.938 , 2 ) ) , 2 ) / 
		pow( ( pow( wgev , 2 ) - pow( 0.938 , 2 ) ) , 2 ) ) ) * 
	    ( (pow( tgev  , 0 ) * (  pow( qsqgev , 0 ) * 0.03229 + pow( qsqgev , 1 ) * -0.003634 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 1 ) * (  pow( qsqgev , 0 ) * -0.07035 + pow( qsqgev , 1 ) * 0.008315 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 2 ) * (  pow( qsqgev , 0 ) * 0.08041 + pow( qsqgev , 1 ) * -0.00939 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
	      (  pow( tgev  , 3 ) * (  pow( qsqgev , 0 ) * -0.03225 + pow( qsqgev , 1 ) * 0.003735 + pow( qsqgev , 2 ) * 0 + pow( qsqgev , 3 ) * 0 + pow( qsqgev , 4 ) * 0 ) ) + 
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


  return s; 
}


