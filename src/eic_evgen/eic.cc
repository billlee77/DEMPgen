///*--------------------------------------------------*/
/// eic.cc:
/// Original author: Dr. Ahmed Zafar
/// Date: 2015-2018
///
///*--------------------------------------------------*/
/// Modifier: Wenliang (Bill) Li
/// Date: Feb 24 2020
/// Email: wenliang.billlee@gmail.com
///
/// Comment: Feb 24, 2020: the main function is excuted in main.cc


#include "eic.h"

using std::setw;
using std::setprecision;
using std::cout;
using std::cin;
using std::endl;
using namespace std;

//---------------------------------------------------------
// g++ -o pim pimFermi.C `root-config --cflags --glibs`
//---------------------------------------------------------

//int main() {
// 
//  eic();
//  
//  return 0;
//}

void eic() {

    Int_t target_direction, kinematics_type;
 
   	cout << "Target Direction (1->Up, 2->Down): "; cin >> target_direction; cout << endl;
   	cout << "Kinematics type (1->FF, 2->TSSA): ";  cin >> kinematics_type;  cout << endl;
   	cout << "Enter the number of events: ";        cin >> fNEvents;         cout << endl;
   	cout << "Enter the file number: ";             cin >> fNFile;           cout << endl;
 
//	eic(target_direction, kinematics_type, fNEvents);

}

/*--------------------------------------------------*/

void eic(int event_number, int target_direction, int kinematics_type, TString file_name, int fEIC_seed, TString particle, TString det_location) {

   	TString targetname;  
	TString charge;

   	if( target_direction == 1 ) targetname = "up";
  	if( target_direction == 2 ) targetname = "down";

	gKinematics_type = kinematics_type;
	gfile_name = file_name;

	fNFile = 1;

	fNEvents = event_number;

	fSeed = fEIC_seed;

	pim* myPim = new pim(fSeed);
  	myPim->Initilize();
 
//  	TDatime dsTime;
//  	cout << "Start Time:   " << dsTime.GetHour() << ":" << dsTime.GetMinute() << endl;

	particle = ExtractParticle(particle);
	charge = ExtractCharge(particle);

//	if (particle == "pi") {
//		Exclusive_Pion_Prodoction(*myPim);
//	} else if (particle == "omega") {
//		Exclusive_Omega_Prodoction(*myPim);
//	} else {
//		cerr << "Choice of particle is not recognized." << endl;
//		exit(0);
//	}

//    TDatime dsTime;
//  	cout << "Start Time:   " << dsTime.GetHour() << ":" << dsTime.GetMinute() << endl;
	
//	TStopwatch tTime;
//   	tTime.Start();
//
//	Exclusive_Pion_Prodoction(*myPim);
//	
//	tTime.Stop();
//   	tTime.Print();

	Reaction* r1 = new Reaction(particle);
	r1->process_reaction();
	delete r1;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void eic(Json::Value obj) {

   	TString targetname;  
 	TString charge;

	int target_direction = obj["Targ_dir"].asInt();
 	gKinematics_type     = obj["Kinematics_type"].asInt();

   	if( target_direction == 1 ) targetname = "up";
   	if( target_direction == 2 ) targetname = "down";
 
 	gfile_name = obj["file_name"].asString();
 
 	gPi0_decay = obj["pi0_decay"].asBool();

 	fNFile = 1;
// 	fNEvents = obj["n_events"].asInt();
 	fNEvents = obj["n_events"].asUInt64();

 	fSeed = obj["generator_seed"].asInt();
 
 	pim* myPim = new pim(fSeed);
   	myPim->Initilize();
 
//  	TDatime dsTime;
//  	cout << "Start Time:   " << dsTime.GetHour() << ":" << dsTime.GetMinute() << endl;

	TString particle = obj["particle"].asString();

	particle = ExtractParticle(particle);
	charge = ExtractCharge(particle);

	///*--------------------------------------------------*/
	/// The detector selection is determined here
	/// The incidence proton phi angle is 

	gDet_location = obj["det_location"].asString();

	if (gDet_location == "ip8") {

		fProton_incidence_phi = 0.0;

	} else if (gDet_location == "ip6") {

		fProton_incidence_phi = fPi;

	} else {
	
		fProton_incidence_phi = 0.0;
		cout << "The interaction point not recognized!" << endl;
		cout << "Therefore default opition ip6 is used." << endl;

	}
	






//	if (particle == "pi") {
//		Exclusive_Pion_Prodoction(*myPim);
//	} else if (particle == "omega") {
//		Exclusive_Omega_Prodoction(*myPim);
//	} else {
//		cerr << "Choice of particle is not recognized." << endl;
//		exit(0);
//	}

//    TDatime dsTime;
//  	cout << "Start Time:   " << dsTime.GetHour() << ":" << dsTime.GetMinute() << endl;
	
//	TStopwatch tTime;
//   	tTime.Start();
//
//	Exclusive_Pion_Prodoction(*myPim);
//	
//	tTime.Stop();
//   	tTime.Print();

	Reaction* r1 = new Reaction(particle);
	r1->process_reaction();
	delete r1;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/













void SetEICSeed (int seed) {
	fSeed = seed;
}




///*--------------------------------------------------*/
///*--------------------------------------------------*/
///  Some utility functions


///*--------------------------------------------------*/
/// Extracting the particle type

TString ExtractParticle(TString particle) {

	/// Make the input particle case insansitive
	particle.ToLower();

	if (particle.Contains("on")) {
		particle.ReplaceAll("on", "");
	};
	
	if (particle.Contains("plus")) {
		particle.ReplaceAll("plus", "+");
	}

	if (particle.Contains("minus")) {
		particle.ReplaceAll("minus", "-");
	}

	if (particle.Contains("zero")) {
		particle.ReplaceAll("zero", "0");
	}

	particle[0] = toupper(particle[0]);

	cout << "Particle: " << particle << endl;

	return particle;

}

///*--------------------------------------------------*/
/// Extracting the particle charge

TString ExtractCharge(TString particle) {

	TString charge;

	if (particle.Contains("+") || particle.Contains("plus")) {
		charge = "+";
	} else if (particle.Contains("-") || particle.Contains("minus")) {
		charge = "-";
	} else {
		charge = "0";
	}

	return charge;

}
