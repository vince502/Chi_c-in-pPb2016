#ifndef ChiFitterInit_h
#define ChiFitterInit_h



#include "ChiTreeInit.h" //declaration of variables in the trees



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

// constraint code is fitted in narrower range, to avoid the turn on which doesn't determine peak shape
const int nMassBinsConstr = 50;
const double mass_windowConstr_l = 3.3;
const double mass_windowConstr_h = 3.8;
const double mass_windowFitConstr_l = 3.3;
const double mass_windowFitConstr_h = 3.8;
const string mass_windowFitConstr = "rvmass>3.3 && rvmass<3.8";

const int nMassBinsJpsi = 100;
const double mass_windowJpsi_l = 2.5;
const double mass_windowJpsi_h = 3.5;
const double mass_windowFitJpsi_l = 2.5;
const double mass_windowFitJpsi_h = 3.5;
const string mass_windowFitJpsi = "rvmassJpsi>2.5 && rvmassJpsi<3.5";

double bins_pT[] = { 6.5, 9, 12, 18, 30 };
int  nbins_pT = sizeof(bins_pT) / sizeof(double) - 1;

double bins_y[] = { -2.4, -1.6, -1.0, 0.0 , 1.0, 1.6, 2.4 };
//double bins_y[] = { -2.4, -1.6, -1.0, 1.0, 1.6, 2.4 };
int  nbins_y = sizeof(bins_y) / sizeof(double) - 1;

double bins_nTrk[] = { 0, 50, 100, 150, 250};
int  nbins_nTrk = sizeof(bins_nTrk) / sizeof(double) - 1;

////// OLD BINS
//std::vector<std::string> fittingSets = { "pt_all", "y", "nTrack", "pt_mid", "pt_fwd", "pt_fwdOnly", "pt_bkwOnly", "pt_fwdOnlyWide", "pt_bkwOnlyWide", "pt_midCMS", "pt_fwdCMS", "pt_bkwCMS", "nTrack_all" };
////const int nFittingSets = fittingSets.size(); //doesn't work for arrays
//const int nFittingSets = 13;
//
////stuff to ease making the plots
//std::vector<std::string> fittingSetsXLabel = { "p_{T} [GeV/c]", "y (J/#psi)", "nTrack", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "p_{T} [GeV/c]", "nTrack" };
//std::vector<double> fittingSetsXLow = { 5, -2.4, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0 };
//std::vector<double> fittingSetsXHigh = { 30, 2.4, 250, 30, 30, 30, 30, 30, 30, 30, 30, 30, 250 };

// rapidity edges for the CM fits, their values in lab-frame rap
double rapCM_Edge1 = -1.535;
double rapCM_Edge2 = -0.535;
double rapCM_Edge3 = 1.465;
double rapCM_Edge4 = 2.4;


//// DEFINE THE KINEMATIC SPLITTING
std::vector<std::string> fittingSets = { "pt_all", "y",  "pt_midCMS", "pt_fwdCMS", "pt_bkwCMS", "nTrack_all" };
//const int nFittingSets = fittingSets.size(); //doesn't work for arrays
const int nFittingSets = 6;

//stuff to ease making the plots
std::vector<std::string> fittingSetsXLabel = { "p_{T}(J/#psi) [GeV/c]", "y(J/#psi)", "p_{T}(J/#psi) [GeV/c]", "p_{T}(J/#psi) [GeV/c]", "p_{T}(J/#psi) [GeV/c]", "N_{tracks}" };
std::vector<double> fittingSetsXLow = { 5, -2.4, 5, 5, 5, 0 };
std::vector<double> fittingSetsXHigh = { 30, 2.4, 30, 30, 30, 250 };


const int nFitFunctionParams = 20; //Real number can be lower, is handled (but not higher)






#endif 

