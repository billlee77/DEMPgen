#include "reaction_routine.h"
#include "eic.h"
#include "particleType.h"

using namespace std;

Pi0_Production::Pi0_Production() { 

	cout << "Program Start" << endl;


}

/*--------------------------------------------------*/
/// Pi0_Production 

Pi0_Production::Pi0_Production(TString particle_str) { 

	rParticle = particle_str;	

//	cout << rParticle << endl;
//	exit(0);

 	Init();
	Pi0_Decay_Pythia6_Out_Init();

}

/*--------------------------------------------------*/

Pi0_Production::~Pi0_Production() {

//	delete rRand;	

//	cout << "File closed!" << endl;

	ppiOut.close();
	ppiDetails.close();

}

/*--------------------------------------------------*/

void Pi0_Production::process_reaction() {

	if_pi0_decay = gPi0_decay;

	polar_out.open("test.txt", ofstream::out);

   	for( long long int i = 0; i < rNEvents; i++ ) {
 
 		rNEvent_itt = i;
 		fNGenerated ++;
 
 		Progress_Report();  // This is happens at each 10% of the total event is processed
 		Processing_Event();
		
 	}

	polar_out.close();
 
 	Detail_Output();

}

/*--------------------------------------------------*/

void Pi0_Production::Processing_Event() {


	// ----------------------------------------------------
    // Considering Fermi momentum for the proton
    // ----------------------------------------------------

	if( kCalcFermi ) {
		Consider_Proton_Fermi_Momentum();
   	}

	// ----------------------------------------------------
    // Boost vector from collider (lab) frame to protons rest frame (Fix target)
    // ----------------------------------------------------
 
    beta_col_rf = r_lproton.BoostVector();        
    fGamma_Col_RF = 1.0/sqrt( 1 - pow( beta_col_rf.Mag() , 2 ) );


    // ---------------------------------------------------------------------
    // Specify the energy and solid angle of scatterd electron in Collider (lab) frame
    // ---------------------------------------------------------------------
//    fScatElec_Theta_Col  = acos( rRand->Uniform( cos( fScatElec_Theta_I ) , cos( fScatElec_Theta_F ) ) );
//    fScatElec_Phi_Col    = rRand->Uniform( 0 , 2.0 * fPi);

//    fScatElec_Phi_Col    = fRandom->Uniform( 0 , 2.0 * fPi);
//    fScatElec_Energy_Col = rRand->Uniform( fScatElec_E_Lo * fElectron_Energy_Col , fScatElec_E_Hi * fElectron_Energy_Col );
//      fScatElec_Energy_Col = fRandom->Uniform( fScatElec_E_Lo * fElectron_Energy_Col , fScatElec_E_Hi * fElectron_Energy_Col );
 
//    fScatElec_Theta_Col  = acos( rRand->Uniform( cos( 120 *fPi /180) , cos( 180 *fPi /180 ) ) );

//    fScatElec_Theta_Col  = rRand->Uniform(fPi/2, fPi);
//    fScatElec_Energy_Col = rRand->Uniform(5000, 6500);


	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/// Random gen
//      fScatElec_Theta_Col  = fRandom->Uniform(fPi/2, fPi);
//      fScatElec_Energy_Col = fRandom->Uniform(5000, 6200);
//	    fScatElec_Phi_Col    = fRandom->Uniform(fProton_Phi_Col - fPi/6.0, fProton_Phi_Col + fPi/6.0);

//	  cout << "Seed check: " << fRandom->GetSeed() << endl;

//	exit(0);

//	fScatElec_Theta_Col = 134 * fPi /180;
//	fScatElec_Energy_Col = 5911;

/*--------------------------------------------------*/
/// 6 GeV
//	fScatElec_Theta_Col = 152 * fPi /180;
//	fScatElec_Energy_Col = 5303;

////

//	fScatElec_Theta_Col = 152 * fPi /180;
//	fScatElec_Energy_Col = 5336;


//	fScatElec_Theta_Col = 152 * fPi /180;
//	fScatElec_Energy_Col = 5330;


/*--------------------------------------------------*/
/// 7 GeV
//	fScatElec_Theta_Col = 150 * fPi /180;
//	fScatElec_Energy_Col = 5350;

/*--------------------------------------------------*/
/// 8 GeV
// 	fScatElec_Theta_Col = 148 * fPi /180;
// 	fScatElec_Energy_Col = 5403;

 	fScatElec_Theta_Col = 148 * fPi /180;
 	fScatElec_Energy_Col = 5360;




/*--------------------------------------------------*/
/// 9 GeV
//	fScatElec_Theta_Col = 146 * fPi /180;
//	fScatElec_Energy_Col = 5458;

/*--------------------------------------------------*/
/// 10.5 GeV
//	fScatElec_Theta_Col = 144 * fPi /180;
//	fScatElec_Energy_Col = 5518;

/*--------------------------------------------------*/
// Q2=10.5 GeV,  S=25
//	fScatElec_Theta_Col = 144 * fPi /180;
//	fScatElec_Energy_Col = 5477;


///*--------------------------------------------------*/
// For special testing only 
//   fScatElec_Phi_Col    = fPi;
//   fScatElec_Phi_Col    = 0.0;
   fScatElec_Phi_Col    = fProton_Phi_Col;

// 	cout << "Initial and final angle: " << fScatElec_Theta_I*180/fPi << "   " << fScatElec_Theta_F*180/fPi << "   " << fScatElec_Theta_Col*180/fPi << "   " << fScatElec_Energy_Col  << "    " << fElectron_Energy_Col << "   " << fScatElec_E_Lo * fElectron_Energy_Col <<  "     " << fScatElec_E_Hi * fElectron_Energy_Col <<  endl;
// 
// 
// 	cout << " electron:::  " << fScatElec_Phi_Col*180/fPi << endl;


    // ----------------------------------------------------
    // Produced Particle X in Collider frame
    // ----------------------------------------------------  

	/// The generic produced particle in the exclusive reaction is labelled as X 
	// fX_Theta_Col      = acos( rRand->Uniform( cos(fX_Theta_I), cos(fX_Theta_F ) ) ); 

	//fX_Theta_Col      = acos( rRand->Uniform( cos(fX_Theta_I), cos(fX_Theta_F ) ) ); 
    // fX_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi);


    //fX_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi);
    //fX_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi);

    /// ----------------------------------------------------  
    /// ----------------------------------------------------  
	/// Random number generation

    // fX_Phi_Col        = fRandom->Uniform( fProton_Phi_Col - fPi/6, fProton_Phi_Col + fPi/6);
	// fX_Theta_Col      = fRandom->Uniform( 0.01, 0.04); 

//    fX_Phi_Col        = fPi;
//    fX_Phi_Col        = 0;
	
	///*--------------------------------------------------*/
	// For special testing only 
	
	//fX_Theta_Col      = 0.01; // 0.5729 degree
	//fX_Theta_Col      = 0.02; // 0.5729 degree
	//fX_Theta_Col      = 0.0; // 0.5729 degree
	//fX_Theta_Col      = 0.025; // Proton incidence angle

	fX_Theta_Col      = 0.0; // Proton incidence angle
//	fX_Theta_Col      = 0.025; // Proton incidence angle
//	fX_Theta_Col      = 0.05; // Proton incidence angle

//  fX_Phi_Col        = 0.0;
//  fX_Phi_Col        = fPi;
    fX_Phi_Col        = fProton_Phi_Col;

	/*--------------------------------------------------*/


// 	cout << "Angle check: " << fX_Theta_Col*180/fPi << "    " << fX_Theta_I*180/fPi << "   "<<  fX_Theta_F*180/fPi << endl;
// 
// 
// 	
// 	cout << "proton phi " << r_lproton.Vect().Theta() << "     " << r_lproton.Vect().Phi() << endl;
// 
// 	cout << "10m radians: " << 0.01*180/fPi << endl;

//	exit(0); 


	// ---------------------------------------------------------------------
    // Specify the energy and solid angle of scatterd electron in Collider (lab) frame
    // ---------------------------------------------------------------------

    fScatElec_Mom_Col  = sqrt( pow( fScatElec_Energy_Col,2) - pow( fElectron_Mass , 2) );
    fScatElec_MomZ_Col = ( fScatElec_Mom_Col * cos(fScatElec_Theta_Col) );  
    fScatElec_MomX_Col = ( fScatElec_Mom_Col * sin(fScatElec_Theta_Col) * cos(fScatElec_Phi_Col) );
    fScatElec_MomY_Col = ( fScatElec_Mom_Col * sin(fScatElec_Theta_Col) * sin(fScatElec_Phi_Col) );

    r_lscatelec.SetPxPyPzE( fScatElec_MomX_Col, fScatElec_MomY_Col, fScatElec_MomZ_Col, fScatElec_Energy_Col);
 
   	r_lscatelecg = r_lscatelec * fm;

    // ----------------------------------------------------
    // Photon in collider (lab) frame and Qsq
    // ----------------------------------------------------

    r_lphoton  = r_lelectron - r_lscatelec;
    r_lphotong = r_lelectrong - r_lscatelecg;

    fQsq_GeV = -1.* r_lphotong.Mag2();

// 	cout << "e: " << r_lelectron.Px() << "  " << r_lelectron.Py() << "  " <<  r_lelectron.Pz() << "  " << r_lelectron.E() << endl;
// 	cout << "Scattered e: " << r_lscatelec.Px() << "  " << r_lscatelec.Py() << "  " <<  r_lscatelec.Pz() << "  " << r_lscatelec.E() << endl;
// 
// 	cout << "Q2 check: " << fQsq_GeV << endl;
// 
// 	exit(0);


//     if ( fQsq_GeV < 5.0 ) {
//        qsq_ev++;
//        return;
//     }

    // ----------------------------------------------------
    // W square, Invariant Mass (P_g + P_p)^2
    // ----------------------------------------------------
 
    TLorentzVector lwg;
    lwg = r_lprotong + r_lphotong;
    fW_GeV    = lwg.Mag();
    fWSq_GeV  = lwg.Mag2();

//  	cout << endl << endl;
//  	cout << "W: "  << fW_GeV << "   s: " << fW_GeV*fW_GeV << "  "  << "Q2: " << fQsq_GeV << endl;
//  	cout << endl << endl;
//		exit(0);



    if ( fWSq_GeV < 0 ) {
		w_neg_ev++;
		return;
    }    

	
 	if (fW_GeV < 0 || fW_GeV > 40) {
		cout << "W out of range of interests" << endl;
      	return;
	}







    // ---------------------------------------------------------
    // Pion momentum in collider frame, analytic solution starts
    // ---------------------------------------------------------

    double fupx = sin( fX_Theta_Col ) * cos( fX_Phi_Col );
    double fupy = sin( fX_Theta_Col ) * sin( fX_Phi_Col );
    double fupz = cos( fX_Theta_Col );
 
    double fuqx = sin( r_lphoton.Theta() ) * cos( r_lphoton.Phi() );
	double fuqy = sin( r_lphoton.Theta() ) * sin( r_lphoton.Phi() );
    double fuqz = cos( r_lphoton.Theta() );
 
    double fa = -(r_lphoton.Vect()).Mag() * ( fupx * fuqx +  fupy * fuqy +  fupz * fuqz );
    double fb = pow ( (r_lphoton.Vect()).Mag() , 2 );
    double fc = r_lphoton.E() + fProton_Mass;

    fa = ( fa - std::abs( (r_lproton.Vect()).Mag() ) * ( ( ( r_lproton.X() / (r_lproton.Vect()).Mag() ) * fupx ) + 
 						       ( ( r_lproton.Y() / (r_lproton.Vect()).Mag() ) * fupy ) + 
 						       ( ( r_lproton.Z() / (r_lproton.Vect()).Mag() ) * fupz ) ) );
     
    double factor = ( pow( (r_lproton.Vect()).Mag() , 2 ) + 2.0 * (r_lphoton.Vect()).Mag() * (r_lproton.Vect()).Mag() *  
 		      ( ( ( r_lproton.X() / (r_lproton.Vect()).Mag() ) * fuqx ) + 
 			( ( r_lproton.Y() / (r_lproton.Vect()).Mag() ) * fuqy ) + 
 			( ( r_lproton.Z() / (r_lproton.Vect()).Mag() ) * fuqz ) ) );
     
    fb =  fb + factor;  
    fc = r_lphoton.E() + r_lproton.E();
     
//     double ft = fc * fc - fb + fPion_Mass * fPion_Mass - fProton_Mass * fProton_Mass;
//         t_min = -qsq + mass**2 -2.*(e_pi0CM*e_photCM -sqrt( (e_pi0CM**2-mass**2)*(e_photCM**2+qsq) ))
//         t_max = -qsq + mass**2 -2.*(e_pi0CM*e_photCM +sqrt((e_pi0CM**2-mass**2)*(e_photCM**2+qsq)))
//         u_min = -qsq + m_psq -2.*(e_pCM*e_photCM -sqrt( (e_pCM**2-m_psq)*(e_photCM**2+qsq) ))


	double e_X_rf      = lX_rf.E();
	double e_photon_rf = lphoton_rf.E();

//	double e_photon_rf = lproton_rf.E();

//	fu_min = -fQsq + pow(f_Scat_Nucleon_Mass, 2) -2.*(e_pi0CM*e_photCM -sqrt( (e_pi0CM**2-mass**2)*(e_photCM**2+qsq) ))

    double ft = fc * fc - fb + fX_Mass * fX_Mass - fProton_Mass * fProton_Mass;

    double fu = lproton_rf.Dot(lproton_rf) + lX_rf.Dot(lX_rf) - lproton_rf.Dot(lX_rf);



//    double fttt = lproton_rfg * lproton_rfg - 2* lproton_rfg * l_scat_nucleon_rf_g + l_scat_nucleon_rf_g * l_scat_nucleon_rf_g;
//    double fuuu = lproton_rfg * lproton_rfg - 2* lproton_rfg * lX_rfg + lX_rfg * lX_rfg;
	




	ft_min = -fQsq_GeV + pow(fX_Mass, 2) -2.*(e_X_rf*e_photon_rf -sqrt( (pow(e_X_rf, 2) - pow(fX_Mass, 2)) * (pow(e_photon_rf,2)+fQsq) ));

	fu_min = -fQsq_GeV + pow(f_Scat_Nucleon_Mass, 2) -2.*(e_X_rf*e_photon_rf -sqrt( (pow(e_X_rf, 2) - pow(f_Scat_Nucleon_Mass , 2))*(pow(e_photon_rf,2)+fQsq) ));


// 	cout << "asd: " << ft_min  << "    " <<  fu_min  <<  f_Scat_Nucleon_Mass  << "   " << fX_Mass << "    " << fProton_Mass << "     "<< fQsq << "    " << e_photon_rf << "    photon photon:  " <<  lphoton_rf.E() << "    " << fm  << "    " << r_lphoton.E() << "    " << beta_col_rf.Mag() << endl;
// 
//     
// 	cout << r_lprotong.Px() << "  " << r_lprotong.Py() << "  " << r_lprotong.Pz() << "  " << r_lprotong.E() << endl;
// 	cout << lproton_rfg.Px() << "  " << lproton_rfg.Py() << "  " << lproton_rfg.Pz() << "  " << lproton_rfg.E() << endl;
// 	cout << beta_col_rf.Px() << "  " << beta_col_rf.Py() << "  " << beta_col_rf.Pz() << endl;

//	exit(0);


    double fQA = 4.0 * ( fa * fa - fc * fc );
    double fQB = 4.0 * fc * ft;

//     double fQC = -4.0 * fa * fa * fPion_Mass * fPion_Mass - ft * ft;    
    double fQC = -4.0 * fa * fa * fX_Mass * fX_Mass - ft * ft;    
 
    fradical = fQB * fQB - 4.0 * fQA * fQC;
 
    fepi1 = ( -fQB - sqrt( fradical ) ) / ( 2.0 * fQA );
    fepi2 = ( -fQB + sqrt( fradical ) ) / ( 2.0 * fQA );
 
    // ---------------------------------------------------------
    // Particle X momentum in collider frame, analytic solution ends
    // ---------------------------------------------------------
         
    r_lX.SetPxPyPzE( (sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * sin(fX_Theta_Col) * cos(fX_Phi_Col),
 			  ( sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * sin(fX_Theta_Col) * sin(fX_Phi_Col),
 			  ( sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * cos(fX_Theta_Col),
 			  fepi1 );

    r_lX_g = r_lX * fm;

	// ----------------------------------------------------
    // Scattered nucleon collider (lab) frame

    r_l_scat_nucleon.SetPxPyPzE( ( r_lproton + r_lelectron - r_lscatelec - r_lX).X(),
 			     ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Y(),
 			     ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Z(),
 			     sqrt( pow( ( ( ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Vect() ).Mag()),2) +
 				   pow( f_Scat_Nucleon_Mass ,2 ) ) );

    r_l_scat_nucleon_g = r_l_scat_nucleon * fm;

//	cout << "Proton Momentum: " << r_l_scat_nucleon_g.Vect().Mag() << "   " << r_l_scat_nucleon_g.Vect().Theta() << "   Angle: "<< r_l_scat_nucleon_g.Vect().Theta() *180/fPi << endl;
//	cout << "X Momentum: " << r_lX_g.Vect().Mag() << "   Angle: " << r_lX_g.Vect().Theta() << "   " << r_lX_g.Vect().Theta() *180/fPi<< endl ;

//	cout << "Proton-X  Angle: " << r_lX_g.Vect().Angle(r_l_scat_nucleon_g.Vect())*180/fPi << endl ; 

    // ----------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------
    // Calculate w = (proton + photon)^2
    // ----------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------
     
//    // cout << fW_GeV << endl;
//    if ( fW_GeV < 3.0 || fW_GeV > 10.6 ) {
//      w_ev++;
//      return;
//    }

     r_lw = r_lproton + r_lphoton;
     fW = r_lw.Mag();

	 // ----------------------------------------------------------------------------------------------
     // ----------------------------------------------------------------------------------------------
     // Calculate w prime w' = (proton + photon - pion)^2                                             
     // ----------------------------------------------------------------------------------------------
     // ----------------------------------------------------------------------------------------------

	 lwp = r_lprotong + r_lphotong - r_lX_g;
     fW_Prime_GeV = lwp.Mag();    

   	 fsini = r_lelectron + r_lproton;
     fsfin = r_lscatelec + r_lX + r_l_scat_nucleon;
     
     fsinig = fsini * fm;
     fsfing = fsfin * fm; 
 
     fMandSConserve = std::abs( fsinig.Mag() - fsfing.Mag() );

     kSConserve = false;
     if( std::abs( fsinig.Mag() - fsfing.Mag() ) < fDiff ) {
       kSConserve = true;
     }
        
     if ( pd->CheckLaws( r_lelectron, r_lproton, r_lscatelec, r_lX, r_l_scat_nucleon) != 1 )
       return;


    ////////////////////////////////////////////////////////////////////////////////////////////
     //                                          Start                                         //
     // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss //
     ////////////////////////////////////////////////////////////////////////////////////////////
 
     lproton_rf = r_lproton;
     lproton_rf.Boost(-beta_col_rf);
     lproton_rfg = lproton_rf * fm;
 
     lelectron_rf = r_lelectron;
     lelectron_rf.Boost(-beta_col_rf);
     lelectron_rfg = lelectron_rf * fm;
 
     lscatelec_rf = r_lscatelec;
     lscatelec_rf.Boost(-beta_col_rf);
     lscatelec_rfg = lscatelec_rf * fm;
     
     lphoton_rf = r_lphoton;
     lphoton_rf.Boost(-beta_col_rf);
     lphoton_rfg = lphoton_rf * fm;
     
     lX_rf = r_lX;
     lX_rf.Boost(-beta_col_rf);
     lX_rfg = lX_rf * fm;
         
     l_scat_nucleon_rf = r_l_scat_nucleon;
     l_scat_nucleon_rf.Boost(-beta_col_rf);
     l_scat_nucleon_rf_g = l_scat_nucleon_rf * fm;



	 ///*--------------------------------------------------*/
	 /// Doing pi0 decay 

	 if (if_pi0_decay) {
		 Pi0_decay(r_lX);
	 } 


     ////////////////////////////////////////////////////////////////////////////////////////////
     //                                          End                                           //
     // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss //
     ////////////////////////////////////////////////////////////////////////////////////////////

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Polar output

	if (fQsq_GeV > 5 && fQsq_GeV < 12 && fWSq_GeV > 8 && fWSq_GeV < 12) {

// 		cout << fQsq_GeV << "  " <<  fWSq_GeV << "  ";             
// 
// 		cout << r_lscatelecg.Vect().Theta()       << "  " <<  r_lscatelecg.Vect().Mag()    << "  " 
// 			 << r_l_scat_nucleon_g.Vect().Theta() << "  " << r_l_scat_nucleon_g.Vect().Mag() << "  " 
//              << r_lX_g.Vect().Theta()             << "  " << r_lX_g.Vect().Mag()                  
// 			 << endl; 
		
		polar_out << fQsq_GeV << "  " <<  fWSq_GeV << "  ";             

		polar_out << r_lscatelecg.Vect().Theta()  << "  " <<  r_lscatelecg.Vect().Mag()      << "  " 
			 << r_l_scat_nucleon_g.Vect().Theta() << "  " << r_l_scat_nucleon_g.Vect().Mag() << "  " 
             << r_lX_g.Vect().Theta()             << "  " << r_lX_g.Vect().Mag()             << "  "                
             << l_photon_1.Vect().Theta()         << "  " << l_photon_1.Vect().Mag()         << "  "                
             << l_photon_2.Vect().Theta()         << "  " << l_photon_2.Vect().Mag()                  
			 << endl; 

	}

//	exit(0);
//	return;


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
























	e_X_rf      = lX_rf.E();

	double e_p_rf      = lproton_rf.E();

	e_photon_rf = lphoton_rf.E();

     // e_photCM = (Wsq - qsq - m_psq)/invm/2.
     // e_omCM   = (Wsq + mass**2 - m_psq)/invm/2.
     // e_pCM    = (Wsq + m_psq - mass**2)/invm/2.


	

//   double e_photCM = (fWSq_GeV - fQsq_GeV - pow(fProton_Mass/1000, 2))/fW_GeV/2.;

     double e_photCM = (fWSq_GeV - fQsq_GeV - pow(fProton_Mass/1000, 2))/fW_GeV/2.;
     double e_pCM = (fWSq_GeV  + pow(fProton_Mass/1000, 2) - pow(f_Scat_Nucleon_Mass/1000, 2))/fW_GeV/2.;

//	 cout << "*********:  " <<  e_photCM   << "     " << e_pCM << endl;
//	 cout << "*********:  " << e_photon_rf << "      "<< e_p_rf << endl;



	ft_min = -fQsq_GeV + pow(fX_Mass/1000, 2) -2.*(e_X_rf/1000*e_photon_rf/1000 -sqrt( (pow(e_X_rf/1000, 2) - pow(fX_Mass/1000, 2)) * (pow(e_photon_rf/1000,2)+fQsq_GeV) ));


//         u_min = -qsq + m_psq -2.*(e_pCM*e_photCM -sqrt( (e_pCM**2-m_psq)*(e_photCM**2+qsq) ))

//	fu_min = -fQsq_GeV + pow(fProton_Mass/1000, 2) -2.*(e_p_rf/1000*e_photon_rf/1000 -sqrt( fabs(pow(e_p_rf/1000, 2) - pow(fProton_Mass/1000 , 2))*(pow(e_photon_rf/1000,2)+fQsq_GeV) ));

	fu_min = -fQsq_GeV + pow(fProton_Mass/1000, 2) -2.*(e_pCM*e_photCM -sqrt( (pow(e_pCM/1000, 2) - pow(fProton_Mass/1000 , 2))*(pow(e_photCM/1000,2)+fQsq_GeV) ));

	
// 	cout << fQsq_GeV << "   " << "AAAA " <<  (pow(e_p_rf/1000, 2) - pow(f_Scat_Nucleon_Mass/1000 , 2))*(pow(e_photon_rf/1000,2)+fQsq_GeV) << "   ::::  " <<  (pow(e_p_rf/1000, 2) - pow(fProton_Mass/1000 , 2)) << "     " << (pow(e_photon_rf/1000,2)+fQsq_GeV) << "      " <<  pow(e_p_rf/1000, 2) << "     " <<  pow(fProton_Mass/1000 , 2) << "     " << e_p_rf/1000*e_photon_rf/1000 << endl;
// 
//     cout << endl;
//     cout << endl;
//     cout << endl;
// 
// 
// 	cout << "asd: " << ft_min  << "    " <<  fu_min << "   " <<  f_Scat_Nucleon_Mass  << "   " << fX_Mass << "    " << fProton_Mass << "     "<< fQsq << "    " << e_photon_rf << "    photon photon:  " <<  lphoton_rf.E() << "    " << fm  << "    " << r_lphoton.E() << "    " << beta_col_rf.Mag() << endl;
// 
//     
// 	cout << "Photon energy: " << e_X_rf << "    " << e_photon_rf << endl;
// 
// 	cout << r_lprotong.Px() << "  " << r_lprotong.Py() << "  " << r_lprotong.Pz() << "  " << r_lprotong.E() << endl;
// 	cout << lproton_rfg.Px() << "  " << lproton_rfg.Py() << "  " << lproton_rfg.Pz() << "  " << lproton_rfg.E() << endl;
// 	cout << "aaa  " << beta_col_rf.Px() << "  " << beta_col_rf.Py() << "  " << beta_col_rf.Pz() << endl;
//
//	exit(0);

     // -----------------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------------
     // Calculate -t
     // -----------------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------------

     fBeta_CM_RF           = (lphoton_rf.Vect()).Mag() / ( lphoton_rf.E() + fProton_Mass );
     fGamma_CM_RF          = ( lphoton_rf.E() + fProton_Mass ) / fW;
     fX_Energy_CM       = ( pow( fW , 2) + pow(fX_Mass , 2) - pow(f_Scat_Nucleon_Mass , 2) ) / ( 2.0 * fW);    
     fX_Mom_CM          = sqrt( pow(fX_Energy_CM , 2) - pow(fX_Mass , 2));    
     fX_Energy_CM_GeV   = fX_Energy_CM / 1000.0;
     fX_Mom_CM_GeV      = fX_Mom_CM / 1000.0;

     // this equation is valid for parallel kinematics only!
     fT_Para = ( pow(((r_lphoton.Vect()).Mag() - (r_lX.Vect()).Mag()),2) - pow((r_lphoton.E() - r_lX.E()),2));
     fT_Para_GeV = fT_Para/1000000.0;

     lt = r_lphoton - r_lX;
     ltg = lt * fm;
 
     fT = -1.*lt.Mag2();
     fT_GeV = -1.*ltg.Mag2();


//	 cout << "ttt:   " << fT_GeV << endl;



     //if ( gKinematics_type == 1 && fT_GeV > 0.5 ) {
     //  t_ev++;
     //  return;
     //}
     //
     //if ( gKinematics_type == 2 && fT_GeV > 1.3 ) {
     //  t_ev++;
     //  return;
     //}

     fx = fQsq_GeV / ( 2.0 * r_lprotong.Dot( r_lphotong ) );
     fy = r_lprotong.Dot( r_lphotong ) / r_lprotong.Dot( r_lelectrong );
     fz = r_lX.E()/r_lphoton.E();    

     // -------------------------------------------------------------------------------------------------------
     // Calculation of Phi  ( azimuthal angle of pion momentum w.r.t lepton plane in target's rest frame)
     // Calculation of PhiS ( azimuthal angle of target polarization w.r.t lepton plane in target's rest frame)
     // -------------------------------------------------------------------------------------------------------  

     v3Photon.SetX( lphoton_rfg.X() );     
	 v3Photon.SetY( lphoton_rfg.Y() );     
	 v3Photon.SetZ( lphoton_rfg.Z() );    

     v3Electron.SetX( lelectron_rfg.X() ); 
     v3Electron.SetY( lelectron_rfg.Y() ); 
     v3Electron.SetZ( lelectron_rfg.Z() );

     v3X.SetX( lX_rfg.X() ) ;        
	 v3X.SetY( lX_rfg.Y() ) ;        
     v3X.SetZ( lX_rfg.Z() );

     v3S.SetX( -1 );                       
	 v3S.SetY( 0 );                        
     v3S.SetZ( 0 );        

     v3PhotonUnit = v3Photon.Unit();    
     v3QxL        = v3Photon.Cross(v3Electron);
     v3QxP        = v3Photon.Cross(v3X);
     v3QxS        = v3Photon.Cross(v3S);
     v3LxP        = v3Electron.Cross(v3X);
     v3LxS        = v3Electron.Cross(v3S);
     v3PxL        = v3X.Cross(v3Electron);
     v3QUnitxL    = v3PhotonUnit.Cross(v3Electron);
     v3QUnitxP    = v3PhotonUnit.Cross(v3X);
     v3QUnitxS    = v3PhotonUnit.Cross(v3S);

	 /*--------------------------------------------------*/
	 // Get the Phi scattering angle with respect to the electron scattering plane
     fPhi   = Get_Phi_X_LeptonPlane_RF ();

	 /*--------------------------------------------------*/
	 // Get the Phi scattering angle between the target polarization plane and the electron scattering plane
     fPhiS  = Get_Phi_TargPol_LeptonPlane_RF();

     fTheta_X_Photon_RF       = fRAD2DEG * acos( ( v3Photon.Dot( v3X     ) ) / ( v3Photon.Mag()  * v3X.Mag()    ) );
     if ( fTheta_X_Photon_RF < 0 ) { fTheta_X_Photon_RF = 180.0 + fTheta_X_Photon_RF; }

     // -----------------------------------------------------------------------------------
     // If we have fermi momentum then epsilon should be in rest frame 
     // The theta angle of scattered angle used in expression of epsilon is the angle 
     // with respect to direction of incoming electron in the rest frame of target nucleon
     // epsilon=1./(1.+ 2.*(pgam_restg**2)/q2g * *(tand(thscat_rest/2.))**2)
     // -----------------------------------------------------------------------------------
 
     double fTheta_EEp = (lelectron_rf.Vect()).Angle(lscatelec_rf.Vect());

     fEpsilon = 1.0 / ( 1.0 + 2.0 * ( pow( (lphoton_rfg.Vect()).Mag(),2)/fQsq_GeV ) * pow( tan( fTheta_EEp / 2 ) , 2 ) );

	 /// meson produced angle with respect to the Q-vector

     theta_X_rf = (lX_rf.Vect()).Angle(lphoton_rf.Vect());
	

     // ----------------------------------------------------
     // Virtual Photon flux factor in units of 1/(GeV*Sr)
     // ----------------------------------------------------
     fFlux_Factor_Col = (fAlpha/(2.0*pow(fPi,2))) * (r_lscatelecg.E() / r_lelectrong.E()) * 
       ( pow(fW_GeV,2) - pow(fProton_Mass_GeV,2) ) / (2.0*fProton_Mass_GeV*fQsq_GeV*(1.0 - fEpsilon));
         
     fFlux_Factor_RF = ( fAlpha / ( 2.0 * pow( fPi , 2 ) ) ) * ( lscatelec_rfg.E() / lelectron_rfg.E() ) *
       ( pow( fW_GeV , 2 ) - pow( fProton_Mass_GeV , 2 ) ) /
       ( 2.0 * fProton_Mass_GeV * fQsq_GeV * ( 1.0 - fEpsilon ) );
     

     // ----------------------------------------------------
     //  Jacobian  dt/dcos(theta*)dphi in units of GeV2/sr
     // ----------------------------------------------------
     fJacobian_CM = ( (lphoton_rfg.Vect()).Mag() - fBeta_CM_RF * lphoton_rfg.E() ) / ( fGamma_CM_RF * ( 1.0 - pow(fBeta_CM_RF,2) ) );
 
     fA = fJacobian_CM * fX_Mom_CM_GeV / fPi;
 
//	cout << "LXXXX check:   " << ( r_lproton + r_lelectron - r_lscatelec - r_lX).X() << "   " << ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Y() << "    " << ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Z() << "     " << sqrt( pow( ( ( ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Vect() ).Mag()),2) + pow( f_Scat_Nucleon_Mass ,2 ) ) << endl;


    double fttt = r_lprotong * r_lprotong - 2 * r_lprotong * r_l_scat_nucleon_g  + r_l_scat_nucleon_g * r_l_scat_nucleon_g;
    double fuuu = r_lprotong * r_lprotong - 2 * r_lprotong * r_lX_g + r_lX_g * r_lX_g;
 
// 	cout << "check this: " << r_lprotong * r_lprotong << "     " <<  r_l_scat_nucleon_g * r_l_scat_nucleon_g << "    " <<  r_lX_g * r_lX_g << endl;
// 
// 	cout << "t: " << ft << "   " << fttt << "    u: " << fu << "    "  << fuuu << endl;  
// 
// 	cout << "Angle aaa : " << r_lproton.Vect().Angle(r_lX.Vect()) << endl;  
// 
// 
// 	exit(0);

     // ----------------------------------------------------
     // Jacobian dOmega* / dOmega dimensionless
     // ----------------------------------------------------
     fJacobian_CM_RF  = ( pow((lX_rf.Vect()).Mag(),2)*fW) / 
       ( fX_Mom_CM * std::abs( ( fProton_Mass + lphoton_rf.E()) * (lX_rf.Vect()).Mag() - 
 				 ( lX_rf.E() * (lphoton_rf.Vect()).Mag() * cos( lX_rf.Theta() ) ) ) );
 
     fJacobian_CM_Col = ( ( pow((r_lX.Vect()).Mag(),2) * fW ) /
     			 ( fX_Mom_CM * std::abs( ( fProton_Mass + r_lphoton.E() ) * (r_lX.Vect()).Mag() -
 						    ( r_lX.E() * (r_lphoton.Vect()).Mag() * cos( r_lX.Theta() ) ) ) ) );

     // -----------------------------------------------------------------------------------------------------------
     // CKY sigma L and T starts
     // -----------------------------------------------------------------------------------------------------------
     // -------------------------------------------------------------------------------------------
 
     r_fSig = Get_CrossSection();
 
     // -----------------------------------------------------------------------------------------------------------
     // CKY sigma L and T ends
     // -----------------------------------------------------------------------------------------------------------
 
     fSigma_Col = r_fSig * fFlux_Factor_Col * fA * fJacobian_CM_Col;


     if ( ( fSigma_Col <= 0 ) || std::isnan( fSigma_Col ) ) { 
       fNSigmaNeg ++;
       return;
     }
     
     // -----------------------------------------------------------------------------------------------------------
     // -----------------------------------------------------------------------------------------------------------
     //             Lab cross section     Phase Space   Conversion     Luminosity                Total events tried
     // Hz        = ub / ( sr^2 * GeV ) * GeV * sr^2 * ( cm^2 / ub ) * ( # / ( cm^2 * sec ) ) / ( # )
 
     fEventWeight = fSigma_Col * fPSF * fuBcm2 * fLumi / fNEvents;   // in Hz

     
     fNRecorded ++;
     fLundRecorded++;
     fRatio = fNRecorded / fNGenerated;
	
//	 cout << if_pi0_decay << endl;

// 	 if (if_pi0_decay) {
// 	 	Pi0_Decay_Lund_Output();
// 	 } else {
// 	 	Pi0_Lund_Output();
// 	 }
//	 cout << "AAAAAAAAAAAAAA" << endl;

	 Pi0_Decay_Pythia6_Output();

//	exit(0);
	
}

/*--------------------------------------------------*/

void Pi0_Production::Detail_Output() {

   ppiDetails << "Total events tried                           " << setw(50) << fNGenerated   << endl;
   ppiDetails << "Total events recorded                        " << setw(50) << fNRecorded    << endl;

   ppiDetails << "Seed used for the Random Number Generator    " << setw(50) << fSeed         << endl;
   ppiDetails << "Number of lund events                        " << setw(50) << fLundRecorded << endl;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Pi0_Production::Pi0_Lund_Output() {


      ppiOut << "5"
  	   << " \t " << fPhi           // var 1
  	   << " \t " << fPhiS          // var 2
  	   << " \t " << fx             // var 3
  	   << " \t " << "1"	       
  	   << " \t " << fQsq_GeV       // var 4
  	   << " \t " << fT_GeV         // var 5
  	   << " \t " << fW_GeV 	       // var 6
  	   << " \t " << fEpsilon       // var 7
  	   << " \t " << fEventWeight   // var 8	   
  	   << endl;
 
 	///*--------------------------------------------------*/
  	// Initial State
 
      ppiOut << setw(10) << "1" 
  	   << setw(10) << "-1" 
  	   << setw(10) << "0" 
  	   << setw(10) << "11"
  	   << setw(10) << "0" 
  	   << setw(10) << "0" 
  	   << setw(16) << r_lelectrong.X()
  	   << setw(16) << r_lelectrong.Y()   
  	   << setw(16) << r_lelectrong.Z()  
  	   << setw(16) << r_lelectrong.E()
  	   << setw(16) << fElectron_Mass_GeV
  	   << setw(16) << fVertex_X
  	   << setw(16) << fVertex_Y
  	   << setw(16) << fVertex_Z
  	   << endl;

      ppiOut << setw(10) << "2" 
  	   << setw(10) << "1" 
  	   << setw(10) << "0" 
  	   << setw(10) << "2212"
  	   << setw(10) << "0" 
  	   << setw(10) << "0" 
  	   << setw(16) << r_lprotong.X()
  	   << setw(16) << r_lprotong.Y()   
  	   << setw(16) << r_lprotong.Z()  
  	   << setw(16) << r_lprotong.E()
  	   << setw(16) << fProton_Mass_GeV
  	   << setw(16) << fVertex_X
  	   << setw(16) << fVertex_Y
  	   << setw(16) << fVertex_Z
  	   << endl;



 	///*--------------------------------------------------*/
  	// Final State
 	      
      // Produced Particle X
      ppiOut << setw(10) << "3" 
  	   << setw(10) << "1" 
  	   << setw(10) << "1" 
  	   << setw(10) << PDGtype(produced_X)
  	   << setw(10) << "0" 
  	   << setw(10) << "0" 
  	   << setw(16) << r_lX_g.X()
  	   << setw(16) << r_lX_g.Y()   
  	   << setw(16) << r_lX_g.Z()  
  	   << setw(16) << r_lX_g.E()
  	   << setw(16) << fX_Mass_GeV
  	   << setw(16) << fVertex_X
  	   << setw(16) << fVertex_Y
  	   << setw(16) << fVertex_Z
  	   << endl;
      
      // Scattered electron
      ppiOut << setw(10) << "4" 
  	   << setw(10) << "-1" 
  	   << setw(10) << "1" 
  	   << setw(10) << "11" 
  	   << setw(10) << "0" 
  	   << setw(10) << "0" 
  	   << setw(16) << r_lscatelecg.X() 
  	   << setw(16) << r_lscatelecg.Y() 
  	   << setw(16) << r_lscatelecg.Z() 
  	   << setw(16) << r_lscatelecg.E()
  	   << setw(16) << fElectron_Mass_GeV
  	   << setw(16) << fVertex_X
  	   << setw(16) << fVertex_Y
  	   << setw(16) << fVertex_Z
  	   << endl;
  	  
      // Recoiled neutron
      ppiOut << setw(10) << "5" 
  	   << setw(10) << "1" 
  	   << setw(10) << "1" 
  	   << setw(10) << PDGtype(recoil_nucleon)
  	   << setw(10) << "0" 
  	   << setw(10) << "0" 
  	   << setw(16) << r_l_scat_nucleon_g.X() 
  	   << setw(16) << r_l_scat_nucleon_g.Y()
  	   << setw(16) << r_l_scat_nucleon_g.Z()
  	   << setw(16) << r_l_scat_nucleon_g.E()
  	   << setw(16) << f_Scat_Nucleon_Mass_GeV
  	   << setw(16) << fVertex_X
  	   << setw(16) << fVertex_Y
  	   << setw(16) << fVertex_Z
  	   << endl;
 
 	
// 	   cout << "Particle check: " << PDGtype(produced_X) << "       " << PDGtype(recoil_nucleon) << endl;
 




// 	///*--------------------------------------------------*/
// 
//      ppiOut << "3"
//  	   << " \t " << fPhi           // var 1
//  	   << " \t " << fPhiS          // var 2
//  	   << " \t " << fx             // var 3
//  	   << " \t " << "1"	       
//  	   << " \t " << fQsq_GeV       // var 4
//  	   << " \t " << fT_GeV         // var 5
//  	   << " \t " << fW_GeV 	       // var 6
//  	   << " \t " << fEpsilon       // var 7
//  	   << " \t " << fEventWeight   // var 8	   
//  	   << endl;
// 
// 	///*--------------------------------------------------*/
//  	// Final State
// 	      
//      // Produced Particle X
//      ppiOut << setw(10) << "1" 
//  	   << setw(10) << "1" 
//  	   << setw(10) << "1" 
//  	   << setw(10) << PDGtype(produced_X)
//  	   << setw(10) << "0" 
//  	   << setw(10) << "0" 
//  	   << setw(16) << r_lX_g.X()
//  	   << setw(16) << r_lX_g.Y()   
//  	   << setw(16) << r_lX_g.Z()  
//  	   << setw(16) << r_lX_g.E()
//  	   << setw(16) << fX_Mass_GeV
//  	   << setw(16) << fVertex_X
//  	   << setw(16) << fVertex_Y
//  	   << setw(16) << fVertex_Z
//  	   << endl;
//      
//      // Scattered electron
//      ppiOut << setw(10) << "2" 
//  	   << setw(10) << "-1" 
//  	   << setw(10) << "1" 
//  	   << setw(10) << "11" 
//  	   << setw(10) << "0" 
//  	   << setw(10) << "0" 
//  	   << setw(16) << r_lscatelecg.X() 
//  	   << setw(16) << r_lscatelecg.Y() 
//  	   << setw(16) << r_lscatelecg.Z() 
//  	   << setw(16) << r_lscatelecg.E()
//  	   << setw(16) << fElectron_Mass_GeV
//  	   << setw(16) << fVertex_X
//  	   << setw(16) << fVertex_Y
//  	   << setw(16) << fVertex_Z
//  	   << endl;
//  	  
//      // Recoiled neutron
//      ppiOut << setw(10) << "3" 
//  	   << setw(10) << "1" 
//  	   << setw(10) << "1" 
//  	   << setw(10) << PDGtype(recoil_nucleon)
//  	   << setw(10) << "0" 
//  	   << setw(10) << "0" 
//  	   << setw(16) << r_l_scat_nucleon_g.X() 
//  	   << setw(16) << r_l_scat_nucleon_g.Y()
//  	   << setw(16) << r_l_scat_nucleon_g.Z()
//  	   << setw(16) << r_l_scat_nucleon_g.E()
//  	   << setw(16) << f_Scat_Nucleon_Mass_GeV
//  	   << setw(16) << fVertex_X
//  	   << setw(16) << fVertex_Y
//  	   << setw(16) << fVertex_Z
//  	   << endl;
// 
// 	
// 	   cout << "Particle check: " << PDGtype(produced_X) << "       " << PDGtype(recoil_nucleon) << endl;
// 


}

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Pi0_Production::Pi0_Decay_Lund_Output() {

//	 cout << rNEvents << endl;
//	 cout << "-----------+++++++++++++-------------------" << endl;
//	 exit(0);




     ppiOut << "4"
 	   << " \t " << fPhi           // var 1
 	   << " \t " << fPhiS          // var 2
 	   << " \t " << fx             // var 3
 	   << " \t " << "1"	       
 	   << " \t " << fQsq_GeV       // var 4
 	   << " \t " << fT_GeV         // var 5
 	   << " \t " << fW_GeV 	       // var 6
 	   << " \t " << fEpsilon       // var 7
 	   << " \t " << fEventWeight   // var 8	   
 	   << endl;
       
//     // Produced Particle X
//     ppiOut << setw(10) << "1" 
// 	   << setw(10) << "1" 
// 	   << setw(10) << "1" 
// 	   << setw(10) << PDGtype(produced_X)
// 	   << setw(10) << "0" 
// 	   << setw(10) << "0" 
// 	   << setw(16) << r_lX_g.X()
// 	   << setw(16) << r_lX_g.Y()   
// 	   << setw(16) << r_lX_g.Z()  
// 	   << setw(16) << r_lX_g.E()
// 	   << setw(16) << fX_Mass_GeV
// 	   << setw(16) << fVertex_X
// 	   << setw(16) << fVertex_Y
// 	   << setw(16) << fVertex_Z
// 	   << endl;
     
     // Scattered electron
     ppiOut << setw(10) << "2" 
 	   << setw(10) << "-1" 
 	   << setw(10) << "1" 
 	   << setw(10) << "11" 
 	   << setw(10) << "0" 
 	   << setw(10) << "0" 
 	   << setw(16) << r_lscatelecg.X() 
 	   << setw(16) << r_lscatelecg.Y() 
 	   << setw(16) << r_lscatelecg.Z() 
 	   << setw(16) << r_lscatelecg.E()
 	   << setw(16) << fElectron_Mass_GeV
 	   << setw(16) << fVertex_X
 	   << setw(16) << fVertex_Y
 	   << setw(16) << fVertex_Z
 	   << endl;
 	  
      // Recoiled neutron
      ppiOut << setw(10) << "3" 
  	   << setw(10) << "1" 
  	   << setw(10) << "1" 
  	   << setw(10) << PDGtype(recoil_nucleon)
  	   << setw(10) << "0" 
  	   << setw(10) << "0" 
  	   << setw(16) << r_l_scat_nucleon_g.X() 
  	   << setw(16) << r_l_scat_nucleon_g.Y()
  	   << setw(16) << r_l_scat_nucleon_g.Z()
  	   << setw(16) << r_l_scat_nucleon_g.E()
  	   << setw(16) << f_Scat_Nucleon_Mass_GeV
  	   << setw(16) << fVertex_X
  	   << setw(16) << fVertex_Y
  	   << setw(16) << fVertex_Z
  	   << endl;

     // Photon 1
     ppiOut << setw(10) << "3" 
 	   << setw(10) << "1" 
 	   << setw(10) << "1" 
 	   << setw(10) << "22"                       
 	   << setw(10) << "0" 
 	   << setw(10) << "0" 
 	   << setw(16) << l_photon_1.X() 
 	   << setw(16) << l_photon_1.Y()
 	   << setw(16) << l_photon_1.Z()
 	   << setw(16) << l_photon_1.E()
 	   << setw(16) << "0"
 	   << setw(16) << fVertex_X
 	   << setw(16) << fVertex_Y
 	   << setw(16) << fVertex_Z
 	   << endl;

     // Photon 2
     ppiOut << setw(10) << "4" 
 	   << setw(10) << "1" 
 	   << setw(10) << "1" 
 	   << setw(10) << "22"                       
 	   << setw(10) << "0" 
 	   << setw(10) << "0" 
 	   << setw(16) << l_photon_2.X() 
 	   << setw(16) << l_photon_2.Y()
 	   << setw(16) << l_photon_2.Z()
 	   << setw(16) << l_photon_2.E()
 	   << setw(16) << "0"
 	   << setw(16) << fVertex_X
 	   << setw(16) << fVertex_Y
 	   << setw(16) << fVertex_Z
 	   << endl;

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Pi0_Production::Pi0_Decay_Pythia6_Out_Init() {


	print_itt = 0;

//	ppiOut << "PYTHIA EVENT FILE" << endl;
	ppiOut << "SIMPLE Event FILE" << endl;
	ppiOut << "============================================" << endl;
	ppiOut << "I, ievent, nParticles" << endl;
	ppiOut << "============================================" << endl;
	ppiOut << "I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)" << endl;
	ppiOut << "============================================" << endl;

}

/*--------------------------------------------------*/

void Pi0_Production::Pi0_Decay_Pythia6_Output() {



//     ppiOut << "4"
// 	   << " \t " << fPhi           // var 1
// 	   << " \t " << fPhiS          // var 2
// 	   << " \t " << fx             // var 3
// 	   << " \t " << "1"	       
// 	   << " \t " << fQsq_GeV       // var 4
// 	   << " \t " << fT_GeV         // var 5
// 	   << " \t " << fW_GeV 	       // var 6
// 	   << " \t " << fEpsilon       // var 7
// 	   << " \t " << fEventWeight   // var 8	   
// 	   << endl;


    ppiOut << "0" << " \t\t\t\ "  << print_itt << " \t\t\t " << "1" << endl;           // var 1

	print_itt++;

	ppiOut << "============================================" << endl;

 	///*--------------------------------------------------*/
  	// Initial State
 
      ppiOut  << "1" 
  	   << setw(6) << "21" 
  	   << setw(6) << "11"
  	   << setw(6) << "0" 
  	   << setw(6) << "3" 
  	   << setw(6) << "4" 

  	   << setw(14) << r_lelectrong.X()
  	   << setw(14) << r_lelectrong.Y()   
  	   << setw(14) << r_lelectrong.Z()  
  	   << setw(14) << r_lelectrong.E()
  	   << setw(14) << fElectron_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;

      ppiOut << "2" 
  	   << setw(6) << "21" 
  	   << setw(6) << "2212"
  	   << setw(6) << "0" 
  	   << setw(6) << "5" 
  	   << setw(6) << "6" 

  	   << setw(14) << r_lprotong.X()
  	   << setw(14) << r_lprotong.Y()   
  	   << setw(14) << r_lprotong.Z()  
  	   << setw(14) << r_lprotong.E()
  	   << setw(14) << fProton_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;

      ppiOut << "3" 
  	   << setw(6) << "21" 
  	   << setw(6) << "22"
  	   << setw(6) << "1" 
  	   << setw(6) << "0" 
  	   << setw(6) << "0" 

  	   << setw(14) << r_lphotong.X()
  	   << setw(14) << r_lphotong.Y()   
  	   << setw(14) << r_lphotong.Z()  
  	   << setw(14) << r_lphotong.E()
  	   << setw(14) << r_lphotong.M()
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;


 	///*--------------------------------------------------*/
  	// Final State
      
      // Scattered electron
      ppiOut << "4" 
  	   << setw(6) << "1" 
  	   << setw(6) << "11" 
  	   << setw(6) << "1" 
  	   << setw(6) << "0"
  	   << setw(6) << "0"
 
  	   << setw(14) << r_lscatelecg.X() 
  	   << setw(14) << r_lscatelecg.Y() 
  	   << setw(14) << r_lscatelecg.Z() 
  	   << setw(14) << r_lscatelecg.E()
  	   << setw(14) << fElectron_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;
  	  
      // Recoiled nucleon
      ppiOut << "5" 
  	   << setw(6) << "1" 
  	   << setw(6) << PDGtype(recoil_nucleon)
  	   << setw(6) << "2" 
  	   << setw(6) << "0"
  	   << setw(6) << "0"
 
  	   << setw(14) << r_l_scat_nucleon_g.X() 
  	   << setw(14) << r_l_scat_nucleon_g.Y()
  	   << setw(14) << r_l_scat_nucleon_g.Z()
  	   << setw(14) << r_l_scat_nucleon_g.E()
  	   << setw(14) << f_Scat_Nucleon_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;
 
      // Produced Particle X
      ppiOut << "6" 
  	   << setw(6) << "1" 
  	   << setw(6) << PDGtype(produced_X)
  	   << setw(6) << "2" 
  	   << setw(6) << "0" 
  	   << setw(6) << "0"

  	   << setw(14) << r_lX_g.X()
  	   << setw(14) << r_lX_g.Y()   
  	   << setw(14) << r_lX_g.Z()  
  	   << setw(14) << r_lX_g.E()
  	   << setw(14) << fX_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;

	ppiOut << "=============== Event finished ===============" << endl;

}



/*--------------------------------------------------*/
/// Cross section model: Based on 12 GeV Hall C u-channel pi0 proposal
/// Author: Garth Huber
/// Date: April 1st, 2020

Double_t Pi0_Production::Get_CrossSection(){

	double_t sig_total;

	// use fit parameters for xB=0.36 data+model (nb/GeV2)
	Double_t Q2tab[4] = {1.75,  3.00,  4.00,  5.50};
    Double_t Wtab[4]  = {2.00,  2.46,  2.83,  3.26};

	// sigT is parameterized in the form p1+p2/(-t)
    Double_t p1T[4]   = {1577.,  168.,  67.4,  24.7};
    Double_t p2T[4]   = {-125., -11.1,  -4.6, -1.66};

	// sigL is simply constant vs t
    Double_t p1L[4]   = {0.,  2.77,  1.16,  0.43};

	// sigTT  p1+p2/(-t)
    Double_t p1TT[4]  = {-674., -141., -57.6, -21.3};
    Double_t p2TT[4]  = { 102.,  25.4,  10.4,  3.82};

	// sigLT  p1+p2/(-t)
    Double_t p1LT[4]  = { 156., 12.57,  5.17,  1.97};
    Double_t p2LT[4]  = {-60.0, -2.08, -0.85, -0.32};

	double sigThi, sigTlo, sigLhi, sigLlo;
	double sigTThi, sigTTlo, sigLThi, sigLTlo;
	double tmin, tprime;
	double umin, uprime;
	Int_t Q2count =0;

	double thetacm = theta_X_rf;

	double eps = fEpsilon;
	double phicm = fPhi;

	double Q2hi=0.;
    double Q2lo=0.;
   	double delQ2=1.;

	double Whi, Wlo;

	double Q2tmp = fQsq_GeV;
	double m_p = fProton_Mass;
	double W_gev = fW_GeV;
	
	double pi = fPi;
	
    int ndat=4; 

	tprime = ft-ft_min;
	uprime = fu-fu_min;

    if( Q2tmp < Q2tab[0] ) {

       	Q2count = 0;
       	Q2hi = Q2tab[1];
      	Whi  = Wtab[1];
       	Q2lo = Q2tab[0];
      	Wlo  = Wtab[0];
       	delQ2 = (Q2hi - Q2lo);

	} else {

		for(int Q2c=0; Q2c <= ndat-2; Q2c++){

//			cout << Q2c << endl;
			if( (Q2tmp >= Q2tab[Q2c] && Q2tmp < Q2tab[Q2c+1]) || Q2tmp >= Q2tab[Q2c+1] ) {

                Q2count = Q2c ;
                Q2hi = Q2tab[Q2count+1];
                Whi  = Wtab[Q2count+1];
                Q2lo = Q2tab[Q2count];
                Wlo  = Wtab[Q2count];
                delQ2 = (Q2hi - Q2lo);

//				cout << "check this:  " << Q2count << "  "<< Q2tmp << "  " << Q2tab[Q2c] << "  "<< Q2tab[Q2c+1] << endl;

			}


//          do Q2c=1,(ndat-1)
//             if( (Q2tmp.ge.Q2tab(Q2c)).and. (Q2tmp.lt.Q2tab(Q2c+1) )
//      1           .or. Q2tmp.ge.Q2tab(Q2c+1) ) then
//                Q2count = Q2c
//                Q2hi = Q2tab(Q2count+1)
//                Whi  = Wtab(Q2count+1)
//                Q2lo = Q2tab(Q2count)
//                Wlo  = Wtab(Q2count)
//                delQ2 = (Q2hi - Q2lo)
//             endif               !Q2 check
//          enddo                  !Q2
		}
	}




	///*--------------------------------------------------*/
	// t-channel
	if (  thetacm > pi/2.) { 

        sigThi = p1T[Q2count+1] + p2T[Q2count+1]/(tprime+fabs(tmin));
        sigTlo = p1T[Q2count]   + p2T[Q2count]/(tprime+fabs(tmin));
        sigLhi = p1L[Q2count+1];
        sigLlo = p1L[Q2count];

        sigTThi = p1TT[Q2count+1]+p2TT[Q2count+1]/(tprime+fabs(tmin));
        sigTTlo = p1TT[Q2count]  +p2TT[Q2count]/(tprime+fabs(tmin));
        sigLThi = p1LT[Q2count+1]+p2LT[Q2count+1]/(tprime+fabs(tmin));
        sigLTlo = p1LT[Q2count]  +p2LT[Q2count]/(tprime+fabs(tmin));

	} else {
//
	///*--------------------------------------------------*/
	// u-channel
//   christian weiss recommends the following change for u-channel:
//   switch u-slope for t-slope, then divide by 10, since back angle peak 
//   is ~10% of forward angle peak (at least for omega electroproduction)
       
        sigThi = (p1T[Q2count+1]+p2T[Q2count+1]/(uprime+fabs(tmin)))/10.;
        sigTlo = (p1T[Q2count]  +p2T[Q2count]/(uprime+fabs(tmin)))/10.  ;
        sigLhi = p1L[Q2count+1]/10.;
        sigLlo = p1L[Q2count]/10.;
        sigTThi = (p1TT[Q2count+1]+p2TT[Q2count+1]/(uprime+fabs(tmin)))/10.;
        sigTTlo = (p1TT[Q2count]  +p2TT[Q2count]/(uprime+fabs(tmin)))/10.;  
        sigLThi = (p1LT[Q2count+1]+p2LT[Q2count+1]/(uprime+fabs(tmin)))/10.;
        sigLTlo = (p1LT[Q2count]  +p2LT[Q2count]/(uprime+fabs(tmin)))/10.; 

    }

//    double Wfac_hi= ((Whi**2 - m_p**2)**2) / ((W_gev**2-m_p**2)**2);
//    double Wfac_lo= ((Wlo**2 - m_p**2)**2) / ((W_gev**2-m_p**2)**2);

    double Wfac_hi = pow(pow(Whi, 2) - pow(m_p, 2), 2) / pow( pow(W_gev,2)- pow(m_p, 2), 2);
    double Wfac_lo = pow(pow(Wlo, 2) - pow(m_p, 2), 2) / pow( pow(W_gev,2)- pow(m_p, 2), 2);

    double sigThiW = sigThi*Wfac_hi;

    if (sigThiW  <0) {
//       cout << pizero: sigThiW<0 ',sigThiW,uprime, abs(tmin)
       sigThiW=0.;
    }

    double sigTloW = sigTlo*Wfac_lo;

    if (sigTloW < 0) {
//       write(6,*)' pizero: sigTloW<0 ',sigTloW,uprime,abs(tmin)
       sigTloW=0.;
    }

    double sigLhiW = sigLhi*Wfac_hi;

    if (sigLhiW < 0) {
//       write(6,*)' pizero: sigLhiW<0 ',sigLhiW,uprime,abs(tmin)
       sigLhiW=0.;
    }

    double sigLloW = sigLlo*Wfac_lo;

    if (sigLloW < 0) {
//       write(6,*)' pizero: sigLloW<0 ',sigLloW,uprime,abs(tmin)
       sigLloW=0.;
    }

    double sigTThiW = sigTThi*Wfac_hi;
    double sigTTloW;
    double sigLThiW;
    double sigLTloW;

	double sig, sigL, sigT, sigLT, sigTT; 

	if (abs(sigTThiW) < sigThiW) {
		sigTThiW=Sign(sigThiW, sigTThiW);
	}
    sigTTloW = sigTTlo*Wfac_lo;

    if (abs(sigTTloW) < sigTloW) {
		sigTTloW = Sign(sigTloW,sigTTloW);
	}
	sigLThiW = sigLThi*Wfac_hi;
    
	if (abs(sigLThiW) < sigThiW) {
		sigLThiW = Sign(sigThiW,sigLThiW);
	}
    sigLTloW = sigLTlo*Wfac_lo;
    
	if (abs(sigLTloW) < sigTloW) {
		sigLTloW = Sign(sigTloW,sigLTloW);
	}

	/*--------------------------------------------------*/

    if( Q2count <= (ndat-2) && Q2tmp >= Q2tab[Q2count] && Q2tmp < Q2tab[Q2count+1] ) {
         
    	sigT  = ( sigTloW*(Q2hi-Q2tmp) + sigThiW*(Q2tmp-Q2lo))/delQ2;
        sigL  = ( sigLloW*(Q2hi-Q2tmp) + sigLhiW*(Q2tmp-Q2lo))/delQ2;
        sigTT = (sigTTloW*(Q2hi-Q2tmp) + sigTThiW*(Q2tmp-Q2lo))/delQ2;
        sigLT = (sigLTloW*(Q2hi-Q2tmp) + sigLThiW*(Q2tmp-Q2lo))/delQ2;

	} else if (Q2tmp >= Q2tab[ndat-1]) {
         
        sigT  =  sigThiW + (sigThiW-sigTloW) /delQ2;
        sigL  =  sigLhiW + (sigLhiW-sigLloW) /delQ2;
        sigTT = sigTThiW + (sigTThiW-sigTTloW)/delQ2;
        sigLT = sigLThiW + (sigLThiW-sigLTloW)/delQ2;

    } else if (Q2tmp <= Q2tab[1] ) {
         
        sigT  =  sigTloW - (sigThiW-sigTloW) /delQ2;
        sigL  =  sigLloW - (sigLhiW-sigLloW) /delQ2;
        sigTT = sigTTloW - (sigTThiW-sigTTloW)/delQ2;
        sigLT = sigLTloW - (sigLThiW-sigLTloW)/delQ2;
	
	} else {

//         write(6,*)' Q2tmp error ',Q2tmp,Q2count
	}
//      endif


      
    sig = sigT + eps*sigL +eps*cos(2.*phicm)*sigTT +sqrt(2.*eps*(1.+eps))*cos(phicm)*sigLT;

    double sig_pi0gmh = sig/2./pi*1.e-09;  //dsig/dtdphicm in microbarns/MeV^2/rad

	return sig_total;
}




///*--------------------------------------------------*/
// Pi0 decay subroutine

void Pi0_Production::Pi0_decay(TLorentzVector pi0_vec) {

// 	TH1D *h1 = new TH1D("h1", "h1", 100, -100, 100);
// 	TH1D *h2 = new TH1D("h2", "h2", 100, -100, 100);
// 	TH1D *h3 = new TH1D("h3", "h3", 100, -100, 100);

// 	TH2D *photon_2d = new TH2D("photon_2d", "photon_2d", 600, -30, 30, 600, -30, 30);

	TVector3 beta_col_rf_pi0;

    beta_col_rf_pi0 = pi0_vec.BoostVector();        
	
    Double_t fGamma_Col_pi0_RF = 1.0/sqrt( 1 - pow( beta_col_rf_pi0.Mag() , 2 ) );

	TLorentzVector l_photon_1_rf;
	TLorentzVector l_photon_2_rf;

//	TLorentzVector pi0_vec;
	TLorentzVector pi0_vec_rf;
	TLorentzVector pi0_vec_rf_g;

    pi0_vec_rf = pi0_vec;
 
    pi0_vec_rf.Boost(-beta_col_rf_pi0);
    pi0_vec_rf_g = pi0_vec_rf * fm;

	double photon_1_E = pi0_vec_rf_g.E()/2;
    double photon_2_E = pi0_vec_rf_g.E()/2;

//	double photon_1_theta     = rRand->Uniform( 0, 2.0 * fPi ); 
//  double photon_1_phi       = rRand->Uniform( 0 , 2.0 * fPi );

//	for (int i =0;  i<=100000; i++) {
	
//	double photon_1_theta     = fPi/4;

//	double photon_1_theta     = fPi/1000 *i; 
//  double photon_1_phi       = fPi/2;

//	double photon_1_theta = rRand->Uniform( 0, fPi ); 
// 	double photon_1_phi   = rRand->Uniform( 0, 2.0 * fPi);
	
	double photon_1_theta = fRandom->Uniform( 0, fPi ); 
  	double photon_1_phi   = fRandom->Uniform( 0, 2.0 * fPi);

    double photon_1_Mom_Col  = photon_1_E;
    double photon_1_MomZ_Col = ( photon_1_Mom_Col * cos(photon_1_theta) );  
    double photon_1_MomX_Col = ( photon_1_Mom_Col * sin(photon_1_theta) * cos(photon_1_phi) );
    double photon_1_MomY_Col = ( photon_1_Mom_Col * sin(photon_1_theta) * sin(photon_1_phi) );

    l_photon_1_rf.SetPxPyPzE( photon_1_MomX_Col, photon_1_MomY_Col, photon_1_MomZ_Col, photon_1_E );
 
// 	photon_2_theta     = photon_1_theta - fPi; 
//  photon_2_phi       = photon_1_phi - fPi;

    l_photon_2_rf.SetPxPyPzE( -photon_1_MomX_Col, -photon_1_MomY_Col, -photon_1_MomZ_Col, photon_1_E );

	l_photon_1 = l_photon_1_rf;
	l_photon_2 = l_photon_2_rf;

	// TVector3 beta_col_rf_pi0;
    // beta_col_rf_pi0 = pi0_vec.BoostVector();

	TVector3 beta_col_rf_pi0_inverse = -pi0_vec.BoostVector();

	l_photon_1.Boost(beta_col_rf_pi0);
	l_photon_2.Boost(beta_col_rf_pi0);

//	cout << "/*--------------------------------------------------*/" << endl;
//	cout << "pi0 decay" << endl;
//  cout << "pi0 decay"
//	cout << "gamma " << fGamma_Col_pi0_RF << "   " << pi0_vec.Gamma() << endl; 
//	cout << "beta " << pi0_vec.Beta()<< endl; 
//	cout << "pion energy: " << pi0_vec.E() * fm << "  " <<  pi0_vec.P() * fm << "   " << l_photon_1_rf.E() << "    " <<  pi0_vec_rf_g.E() << endl; 
//	cout << "Photon energy 1: " << l_photon_1_rf.E() <<  "  P: " << l_photon_1_rf.P()   << "  After transformation   " << l_photon_1.E() << "   " << l_photon_1.P() << endl; 
//	cout << "Photon energy 2: " << l_photon_2_rf.E() <<  "  P: " <<  l_photon_2_rf.P()  << "  After transformation   " << l_photon_2.E() << "   " << l_photon_2.P() << endl; 	
//	cout << "angle: " << l_photon_1_rf.Vect().Angle(l_photon_2_rf.Vect()) * 180.0/fPi 
//		 << "    " << l_photon_1.Vect().Angle(l_photon_2.Vect()) * 180.0/fPi << endl;
//	exit(0);

	double photons_1_sep = (tan(l_photon_1.Vect().Angle(pi0_vec.Vect())) * 32)*100; 
	double photons_2_sep = (tan(l_photon_2.Vect().Angle(pi0_vec.Vect())) * 32)*100; 

	double photons_seperation = photons_1_sep + photons_2_sep;
 
//	cout << "pi0_vec: "  << pi0_vec.Vect().Theta()    << "    " << pi0_vec.Vect().Theta()  << endl;
//	cout << "Photon 1: " << l_photon_1.Vect().Theta() << "    " << l_photon_1.Vect().Phi() << endl;
//	cout << "Photon 2: " << l_photon_2.Vect().Theta() << "    " << l_photon_2.Vect().Phi() << endl;
//	cout << "photon prjecton at 32 meters: "          << sin(l_photon_1.Vect().Angle(pi0_vec.Vect())) * 32 
//	                               << "    "          << sin(l_photon_2.Vect().Angle(pi0_vec.Vect())) * 32 << "    "
//								   << photons_seperation << endl;
//	exit(0);

	///*--------------------------------------------------*/
	/// Max angle formula: sin(theta) = m_pi/(2*E_gamma)  
	//
	// Decay at 7.8m
	// 2.6033 ± 0.0005 × 10−8 life time

	/// ZDC is at 30 down stream
	// Pi0 drift is 22.2 meters

	// cout << "Maximum angle: " << pi0_vec_rf_g.E() << "   " << pi0_vec.E()*fm << "    Angle in deg: " << asin(pi0_vec_rf_g.E()/(pi0_vec.E() * fm)) * 180/fPi << endl; 

	// 2.6*3 = 7.8m before decay 

	 TLorentzVector photon_1_diff = l_photon_1 - pi0_vec*fm;

//  	cout << "aa" << endl;
//  	cout << pi0_vec.Px() *fm << "  " << pi0_vec.Py() *fm<< "  " << pi0_vec.Pz() *fm<< "  " << pi0_vec.E() *fm<< endl;
//  	cout << l_photon_1.Px() << "  " << l_photon_1.Py() << "  " << l_photon_1.Pz() << "  " << l_photon_1.E() << endl;
//  	cout << photon_1_diff.Px() << "  " << photon_1_diff.Py() << "  " << photon_1_diff.Pz() << "  " << photon_1_diff.E() << endl;
	

	double angle = atan(l_photon_1.Px()/l_photon_1.Py());

//	h1->Fill(photons_seperation);	
//	h2->Fill(photons_1_sep);	
//	h3->Fill(-photons_2_sep);	

//	cout << photons_1_sep  << "  " << l_photon_1.Vect().Phi() << "   " << photons_2_sep << endl;

//	photon_2d->Fill(photons_1_sep *sin(l_photon_1.Vect().Angle(pi0_vec.Vect())) * cos(l_photon_1.Vect().Phi()), photons_1_sep* sin(l_photon_1.Vect().Angle(pi0_vec.Vect())) *sin(l_photon_1.Vect().Phi()) );
//	photon_2d->Fill(photons_2_sep *cos(l_photon_2.Vect().Phi()), photons_2_sep*sin(l_photon_2.Vect().Phi()));

//photons_2_sep, l_photon_2.Vect().Phi()
	
//	cout << l_photon_1.Vect().Phi() << " " << l_photon_2.Vect().Phi() << endl;
//	cout <<"check:  " << photons_1_sep*sin(angle) << "  "  << photons_1_sep*cos(angle) << "    " << angle << endl;
//	photon_2d->Fill(l_photon_1.Vect().Phi(), l_photon_2.Vect().Phi());
//  cout << l_photon_1.Px() << "    " << l_photon_1.Py() << endl;
//	cout <<  photon_1.Gamma() << endl;
//	cout << l_photon_1.Px()**1.1e-7 << "    " << l_photon_1.Py()**1.1e-7 << endl;
//	photon_2d->Fill(l_photon_1.Px()/0.134/l_photon_1.Gamma()*1.1e-7, l_photon_1.Py()/0.134/l_photon_1.Gamma()*1.1e-7);
//	}

// 	TCanvas* c1 = new TCanvas();
// 
// 	h3->SetFillColor(2);
// 	h3->SetLineColor(2);
// 
// 	h2->SetFillColor(3);
// 	h2->SetLineColor(3);
// 
// // h1->SetFillColor(4);
// 	h1->SetLineColor(4);
// 
// 	h1->Draw("hist");
// 	h2->Draw("same");
// 	h3->Draw("same");
// 
// 	h1->Draw("same");
// 
// 	c1->Print("seperation.png");	
//	photon_2d->Draw("colz");
// 	c1->Print("seperation_aaa.png");	
//	cout << "/*--------------------------------------------------*/" << endl;
//	exit(0);
	
}

