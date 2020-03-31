#include "reaction_rountine.h"
#include "eic.h"
#include "particleType.h"

using namespace std;

Pi0_Production::Pi0_Production() { 

	cout << "Program Start" << endl;

}



/*--------------------------------------------------*/
/// PiPlus_Production 


Pi0_Production::Pi0_Production(TString particle_str) { 

	rParticle = particle_str;	
 	Init();

}

Pi0_Production::~Pi0_Production() {

	ppiOut.close();
	ppiDetails.close();

}

void Pi0_Production::process_reaction() {
 
   	for( long long int i = 0; i < rNEvents; i++ ) {
 
 		rNEvent_itt = i;
 		fNGenerated ++;
 
 		Progress_Report();  // This is happens at each 10% of the total event is processed
 		Processing_Event();
		
 	}
 
 	Detail_Output();
 
}



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
    fScatElec_Theta_Col  = acos( fRandom->Uniform( cos( fScatElec_Theta_I ) , cos( fScatElec_Theta_F ) ) );
    fScatElec_Phi_Col    = fRandom->Uniform( 0 , 2.0 * fPi);
    fScatElec_Energy_Col = fRandom->Uniform( fScatElec_E_Lo * fElectron_Energy_Col , fScatElec_E_Hi * fElectron_Energy_Col );

    // ----------------------------------------------------
    // Produced Particle X in Collider frame
    // ----------------------------------------------------  

	/// The generic produced particle in the exclusive reaction is labelled as X 
	fX_Theta_Col      = acos( fRandom->Uniform( cos(fX_Theta_I), cos(fX_Theta_F ) ) ); 
    fX_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi );

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

    if ( fWSq_GeV < 0 ) { 
      w_neg_ev++;
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
	double e_photon_rf = lproton_rf.E();

//	fu_min = -fQsq + pow(f_Scat_Nucleon_Mass, 2) -2.*(e_pi0CM*e_photCM -sqrt( (e_pi0CM**2-mass**2)*(e_photCM**2+qsq) ))

    double ft = fc * fc - fb + fX_Mass * fX_Mass - fProton_Mass * fProton_Mass;

    double fu = lproton_rf.Dot(lproton_rf) + lX_rf.Dot(lX_rf) - lproton_rf.Dot(lX_rf);
 
	ft_min = -fQsq + pow(fX_Mass, 2) -2.*(e_X_rf*e_photon_rf -sqrt( (pow(e_X_rf, 2) - pow(fX_Mass, 2)) * (pow(e_photon_rf,2)+fQsq) ));

	fu_min = -fQsq + pow(f_Scat_Nucleon_Mass, 2) -2.*(e_X_rf*e_photon_rf -sqrt( (pow(e_X_rf, 2) - pow(f_Scat_Nucleon_Mass , 2))*(pow(e_photon_rf,2)+fQsq) ));



    
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
    // Scattered proton collider (lab) frame

    r_l_scat_nucleon.SetPxPyPzE( ( r_lproton + r_lelectron - r_lscatelec - r_lX).X(),
 			     ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Y(),
 			     ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Z(),
 			     sqrt( pow( ( ( ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Vect() ).Mag()),2) +
 				   pow( f_Scat_Nucleon_Mass ,2 ) ) );

    r_l_scat_nucleon_g = r_l_scat_nucleon * fm;

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

     ////////////////////////////////////////////////////////////////////////////////////////////
     //                                          End                                           //
     // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss //
     ////////////////////////////////////////////////////////////////////////////////////////////


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

	 Lund_Output();
	
}

/*--------------------------------------------------*/

void Pi0_Production::Detail_Output() {

   ppiDetails << "Total events tried                           " << setw(50) << fNGenerated   << endl;
   ppiDetails << "Total events recorded                        " << setw(50) << fNRecorded    << endl;

   ppiDetails << "Seed used for the Random Number Generator    " << setw(50) << fSeed         << endl;
   ppiDetails << "Number of lund events                        " << setw(50) << fLundRecorded << endl;

}

/*--------------------------------------------------*/

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

			cout << Q2c << endl;
			if( (Q2tmp >= Q2tab[Q2c] && Q2tmp < Q2tab[Q2c+1]) || Q2tmp >= Q2tab[Q2c+1] ) {

                Q2count = Q2c ;
                Q2hi = Q2tab[Q2count+1];
                Whi  = Wtab[Q2count+1];
                Q2lo = Q2tab[Q2count];
                Wlo  = Wtab[Q2count];
                delQ2 = (Q2hi - Q2lo);

				cout << "check this:  " << Q2count << "  "<< Q2tmp << "  " << Q2tab[Q2c] << "  "<< Q2tab[Q2c+1] << endl;

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

//	cout << "Q2: " << Q2tmp << "   " << Q2count << endl;
//	exit(0);


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
