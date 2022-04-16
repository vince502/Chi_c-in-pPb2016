#ifndef ChiFitterInit_h
#define ChiFitterInit_h

 // Declaration of branches that are saved for chi


//#include "TLorentzVector.h"
//#include "TVector3.h"
//#include "TTree.h"
//#include <TClonesArray.h>
//#include <vector>
//#include <sstream>

//#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "FWCore/Common/interface/TriggerNames.h"
//#include "DataFormats/HeavyIonEvent/interface/Centrality.h"



#include "ChiTreeInit.h" //declaration of variables in the trees



const double mass_cutoffJpsi_l = 2.9; //cutoff for Jpsi mass that is accepted to be chic candidate
const double mass_cutoffJpsi_h = 3.25;

//////////////////////////////////////
/////// fitting variables  //////////
/////////////////////////////////////
 
const int PythCode_chic0 = 10441; //Pythia codes
const int PythCode_chic1 = 20443;
const int PythCode_chic2 = 445;

const double k_mass_c0 = 3.4148; //pdg 1. 2018
const double k_mass_c1 = 3.5107;
const double k_mass_c2 = 3.5562;

//const int nMassBins = 80;
//const double mass_window_l = 3.35;
//const double mass_window_h = 3.75;
//const double mass_windowFit_l = 3.35;
//const double mass_windowFit_h = 3.75;
//const string mass_windowFit = "rvmass>3.35 && rvmass<3.75";

//const int nMassBins = 100;
//const double mass_window_l = 3.35;
//const double mass_window_h = 3.85;
//const double mass_windowFit_l = 3.35;
//const double mass_windowFit_h = 3.85;
//const string mass_windowFit = "rvmass>3.35 && rvmass<3.85";

const int nMassBins = 60;
const double mass_window_l = 3.2;
const double mass_window_h = 3.8;
const double mass_windowFit_l = 3.2;
const double mass_windowFit_h = 3.8;
const string mass_windowFit = "rvmass>3.2 && rvmass<3.8";

const int nMassBinsJpsi = 150;
const double mass_windowJpsi_l = 2.5;
const double mass_windowJpsi_h = 4.0;
const double mass_windowFitJpsi_l = 2.5;
const double mass_windowFitJpsi_h = 4.0;
const string mass_windowFitJpsi = "rvmassJpsi>2.5 && rvmassJpsi<4.0";

double bins_pT[] = { 6.5, 9, 12, 18, 30 };
int  nbins_pT = sizeof(bins_pT) / sizeof(double) - 1;

double bins_y[] = { -2.4, -1.6, -1.0, 0.0 , 1.0, 1.6, 2.4 };
//double bins_y[] = { -2.4, -1.6, -1.0, 1.0, 1.6, 2.4 };
int  nbins_y = sizeof(bins_y) / sizeof(double) - 1;

double bins_nTrk[] = { 0, 50, 100, 150, 250, 400 };
int  nbins_nTrk = sizeof(bins_nTrk) / sizeof(double) - 1;









#endif 

