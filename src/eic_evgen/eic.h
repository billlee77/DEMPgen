# ifndef EIC_H
# define EIC_H

#include <iostream>
#include "TRandom.h"
#include "TRandom2.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TF1.h>

#include "eic_pim.h"

#include "tssa_sig_Para.h"

#include "reaction_rountine.h"
#include "legacy.h"


void eic();
void eic(int, int, int, TString, int, TString);

extern int fSeed;

void SetEICSeed(int);

TString ExtractParticle(TString);
TString ExtractCharge(TString) ;

#endif

