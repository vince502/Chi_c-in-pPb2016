#include <iostream>
#include <string.h>

#include "TSystem.h"
#include "TTree.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TPave.h"
#include "TText.h"

#include "TROOT.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TLorentzVector.h"

using namespace std;

const int PythCode_chic0 = 10441; //Pythia codes
const int PythCode_chic1 = 20443;
const int PythCode_chic2 = 445;
double binsChiEffpT[] = {0.0, 0.8, 1.6, 2.4, 3.6, 5.0, 7.0, 10.0, 20.0};
//double binsChiEffpT[] = { 0.0, 0.8, 1.6, 2.4, 3.6, 6.0, 8.0, 10.0, 20.0 };
int  nbinsChiEffpT = sizeof(binsChiEffpT) / sizeof(double) - 1;
double binsConvEffpT[] = { 0.0, 0.1, 0.5, 1.0, 1.5, 2.5, 5.0};
int  nbinsConvEffpT = sizeof(binsConvEffpT) / sizeof(double) - 1;
int weird_decay_counter = 0;

//variables

std::vector <bool>* muonIsGlobal = 0;
std::vector <bool>* muonIsTracker = 0;
std::vector <int>* muonTrackerLayersWithMeasurement = 0;
std::vector <bool>* muonIsSoft = 0;
std::vector <double>* muon_eta = 0;
std::vector <double>* muon_pt = 0;

// MC general

std::vector <int>* gen_pdgId = 0;
std::vector <double>* gen_chic_pt = 0;
std::vector <double>* gen_chic_eta = 0;
TClonesArray* gen_chic_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <double>* gen_Jpsi_pt = 0;
std::vector <double>* gen_Jpsi_eta = 0;
std::vector <int>* gen_Jpsi_matchPosition = 0;
std::vector <int>* gen_Jpsi_nMatches = 0;
std::vector <double>* gen_Jpsi_rDelta = 0; //in principle duplicates information
std::vector <double>* gen_Jpsi_ptDeltaRel = 0;//in principle duplicates information
TClonesArray* gen_Jpsi_p4 = new TClonesArray("TLorentzVector", 100);

std::vector <int>* gen_muon_charge = 0;
std::vector <int>* gen_muon_nMatches = 0;
std::vector <int>* gen_muon_matchPosition = 0;
std::vector <double>* gen_muon_eta = 0;
std::vector <double>* gen_muon_pt = 0;
std::vector <double>* gen_muon_rDelta = 0;
std::vector <double>* gen_muon_ptDeltaRel = 0;

std::vector <double>* gen_phot_pt = 0;
std::vector <double>* gen_phot_eta = 0;
std::vector <int>* gen_conv_matchPosition = 0;
std::vector <int>* gen_conv_nMatches = 0;
std::vector <double>* gen_conv_rDelta = 0;
std::vector <double>* gen_conv_ptDeltaRel = 0;



//// Conversions  /////

std::vector <bool>* convQuality_isHighPurity = 0;
std::vector <bool>* convQuality_isGeneralTracksOnly = 0;
std::vector <double>* conv_vertexPositionRho = 0;
std::vector <double>* conv_sigmaTkVtx1 = 0;
std::vector <double>* conv_sigmaTkVtx2 = 0;
std::vector <bool>* conv_tkVtxCompatibilityOK = 0;

std::vector <int>* conv_compatibleInnerHitsOK = 0; //-1: less than 2 tracks, 0: not compatible, 1: yes
std::vector <double>* conv_vertexChi2Prob = 0;
std::vector <double>* conv_zOfPriVtx = 0;
std::vector <double>* conv_zOfPriVtxFromTracks = 0;
std::vector <double>* conv_dzToClosestPriVtx = 0;
std::vector <double>* conv_dxyPriVtx_Tr1 = 0;
std::vector <double>* conv_dxyPriVtx_Tr2 = 0;
std::vector <double>* conv_dxyPriVtxTimesCharge_Tr1 = 0;
std::vector <double>* conv_dxyPriVtxTimesCharge_Tr2 = 0;
std::vector <double>* conv_dxyError_Tr1 = 0;
std::vector <double>* conv_dxyError_Tr2 = 0;

std::vector <int>* conv_tk1NumOfDOF = 0;
std::vector <int>* conv_tk2NumOfDOF = 0;
std::vector <double>* conv_track1Chi2 = 0;
std::vector <double>* conv_track2Chi2 = 0;
std::vector <double>* conv_minDistanceOfApproach = 0;
std::vector <double>* conv_eta = 0;
std::vector <double>* conv_pt = 0;


// chi
std::vector <double>* chi_eta = 0;
std::vector <double>* chi_pt = 0;
std::vector <int>* chi_daughterJpsi_position = 0; //stores position of daughter Jpsi in dimuon collection
std::vector <int>* chi_daughterConv_position = 0; //stores position of daughter photon (conversion)
std::vector <double>* chi_dzPhotToDimuonVtx = 0; //z distance of photon to dimuon vertex when dxy is minimal
std::vector <double>* chi_dxyPhotToDimuonVtx = 0; //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions



bool MuonAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;
	if (fabs(eta) < 0.8 && pt < 3.3) return false;
	if (fabs(eta) >= 0.8 && fabs(eta) < 1.5 && pt < 5.81 - 3.14*fabs(eta)) return false;
	if (fabs(eta) >= 1.5 && (pt < 0.8 || pt < 1.89 - 0.526*fabs(eta))) return false;
	return true;
}

bool PhotAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.5) return false;
	if (pt < 0.2) return false;
	return true;
}

bool MuonSelectionPassMC(int muonPos)  //uses variables loaded in main function
{
	int matchPosition = gen_muon_matchPosition->at(muonPos);
	if (matchPosition < -0.5) return false; //not matched
	if (muonIsSoft->at(matchPosition) != 1) return false;
	return true;
}

bool PhotSelectionPassMC(int photPos)  //uses variables loaded in main function
{
	int matchPosition = gen_conv_matchPosition->at(photPos);
	if (matchPosition < -0.5) return false; //not matched
	if (convQuality_isHighPurity->at(matchPosition) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(matchPosition) != 1) return false;
	if (conv_vertexPositionRho->at(matchPosition) <= 1.5) return false;
	if (conv_sigmaTkVtx1->at(matchPosition) > 10) return false;
	if (conv_sigmaTkVtx2->at(matchPosition) > 10) return false;
	if (conv_tkVtxCompatibilityOK->at(matchPosition) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(matchPosition) != 1) return false;
	if (conv_vertexChi2Prob->at(matchPosition) <= 0.01) return false;
	//if (fabs(conv_zOfPriVtx->at(matchPosition)) >= 20) return false;
	if (fabs(conv_dzToClosestPriVtx->at(matchPosition))>= 10) return false;
	if (conv_tk1NumOfDOF->at(matchPosition) < 2.5) return false;
	if (conv_tk2NumOfDOF->at(matchPosition) < 2.5) return false;
	if (conv_track1Chi2->at(matchPosition) >= 10) return false;
	if (conv_track2Chi2->at(matchPosition) >= 10) return false;
	if (conv_minDistanceOfApproach->at(matchPosition) <=-0.25) return false;
	if (conv_minDistanceOfApproach->at(matchPosition) >= 1.00) return false;

	return true;
}


bool PhotSelectionPassMCLoose(int photPos)  //uses variables loaded in main function
{
	int matchPosition = gen_conv_matchPosition->at(photPos);
	if (matchPosition < -0.5) return false; //not matched
	if (convQuality_isHighPurity->at(matchPosition) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(matchPosition) != 1) return false;
	//if (conv_vertexPositionRho->at(matchPosition) <= 1.5) return false;
	//if (conv_sigmaTkVtx1->at(matchPosition) > 10) return false;
	//if (conv_sigmaTkVtx2->at(matchPosition) > 10) return false;
	//if (conv_tkVtxCompatibilityOK->at(matchPosition) != 1) return false;
	//if (conv_compatibleInnerHitsOK->at(matchPosition) != 1) return false;
	//if (conv_vertexChi2Prob->at(matchPosition) <= 0.01) return false;
	////if (fabs(conv_zOfPriVtx->at(matchPosition)) >= 20) return false;
	//if (fabs(conv_dzToClosestPriVtx->at(matchPosition)) >= 10) return false;
	//if (conv_tk1NumOfDOF->at(matchPosition) < 2.5) return false;
	//if (conv_tk2NumOfDOF->at(matchPosition) < 2.5) return false;
	//if (conv_track1Chi2->at(matchPosition) >= 10) return false;
	//if (conv_track2Chi2->at(matchPosition) >= 10) return false;
	//if (conv_minDistanceOfApproach->at(matchPosition) <= -0.25) return false;
	//if (conv_minDistanceOfApproach->at(matchPosition) >= 1.00) return false;

	return true;
}


int VariablesDistributions(const char* fileIn = "Chi_c_pPb8TeV_MCv5.root", const char* fileOut = "Chi_c_distributionsMC5.root")
{
	gStyle->SetOptStat(1111);
	gROOT->ProcessLine("#include <vector>");

	//TH2F* hdxy_dz = new TH2F("hdxy_dz", "pt vs abseta", 250, 0, 2.5, 200, 0, 50);
	//TH1D* hSignal = new TH1D("hSignal", "", 600, 3, 6);
	//TH1D* hSignal = new TH1D("hSignal", "", 120, 3.2, 3.8);
	TH1D* h_muon_ptRel1 = new TH1D("h_muon_ptRel1", "deltaPtRel from matching", 100, 0.0, 0.5);
	TH1D* h_muon_ptRel2 = new TH1D("h_muon_ptRel2", "deltaPtRel by hand", 100, 0.0, 0.5);

	TH1I* h_muonMatch_isGlobal = new TH1I("h_muonMatch_isGlobal", "h_muonMatch_isGlobal", 2, 0, 2);
	TH1I* h_muonMatch_isTracker = new TH1I("h_muonMatch_isTracker", "h_muonMatch_isTracker", 2, 0, 2);
	TH1I* h_muonMatch_isSoft = new TH1I("h_muonMatch_isSoft", "h_muonMatch_isSoft", 2, 0, 2);
	TH1I* h_muonMatch_trackerLayers = new TH1I("h_muonMatch_trackerLayers", "Tracker layers with measurement", 20, 0, 20);

	TH1I* h_gen_muon_matchPosition = new TH1I("h_gen_muon_matchPosition", "h_gen_muon_matchPosition", 10, 0, 10);
	TH1I* h_gen_muon_nMatches = new TH1I("h_gen_muon_nMatches", "h_gen_muon_nMatches", 10, 0, 10);
	TH1D* h_gen_muon_rDelta = new TH1D("h_gen_muon_rDelta", "h_gen_muon_rDelta", 100, 0, 0.5);


	//conversions
	TH1I* h_convQuality_isHighPurity = new TH1I("h_convQuality_isHighPurity", "h_convQuality_isHighPurity", 2, 0, 2);
	TH1I* h_convQuality_isGeneralTracksOnly = new TH1I("h_convQuality_isGeneralTracksOnly", "h_convQuality_isGeneralTracksOnly", 2, 0, 2);
	TH1D* h_conv_vertexPositionRho = new TH1D("h_conv_vertexPositionRho", "h_conv_vertexPositionRho", 200, 0, 100);
	TH1D* h_conv_sigmaTkVtx1 = new TH1D("h_conv_sigmaTkVtx1", "h_conv_sigmaTkVtx1", 200, 0, 100);
	TH1D* h_conv_sigmaTkVtx2 = new TH1D("h_conv_sigmaTkVtx2", "h_conv_sigmaTkVtx2", 200, 0, 100);
	TH1I* h_conv_tkVtxCompatibilityOK = new TH1I("h_conv_tkVtxCompatibilityOK", "h_conv_tkVtxCompatibilityOK", 2, 0, 2);

	TH1I* h_conv_compatibleInnerHitsOK = new TH1I("h_conv_compatibleInnerHitsOK", "h_conv_compatibleInnerHitsOK", 3, -1, 2); //-1: less than 2 tracks, 0: not compatible, 1: yes
	TH1D* h_conv_vertexChi2Prob = new TH1D("h_conv_vertexChi2Prob", "h_conv_vertexChi2Prob", 200, 0, 1);
	TH1D* h_conv_zOfPriVtx = new TH1D("h_conv_zOfPriVtx", "h_conv_zOfPriVtx", 240, -60, 60);
	TH1D* h_conv_zOfPriVtxFromTracks = new TH1D("h_conv_zOfPriVtxFromTracks", "h_conv_zOfPriVtxFromTracks", 240, -60, 60);
	TH1D* h_conv_dzToClosestPriVtx = new TH1D("h_conv_dzToClosestPriVtx", "h_conv_dzToClosestPriVtx", 240, -60, 60);
	TH1D* h_conv_dxyPriVtx_Tr1 = new TH1D("h_conv_dxyPriVtx_Tr1", "h_conv_dxyPriVtx_Tr1", 400, -100, 100);
	TH1D* h_conv_dxyPriVtx_Tr2 = new TH1D("h_conv_dxyPriVtx_Tr2", "h_conv_dxyPriVtx_Tr2", 400, -100, 100);
	TH1D* h_conv_dxyPriVtxTimesCharge_Tr1 = new TH1D("h_conv_dxyPriVtxTimesCharge_Tr1", "h_conv_dxyPriVtxTimesCharge_Tr1", 400, -100, 100);
	TH1D* h_conv_dxyPriVtxTimesCharge_Tr2 = new TH1D("h_conv_dxyPriVtxTimesCharge_Tr2", "h_conv_dxyPriVtxTimesCharge_Tr2", 400, -100, 100);
	TH1D* h_conv_dxyError_Tr1 = new TH1D("h_conv_dxyError_Tr1", "h_conv_dxyError_Tr1", 200, 0, 10);
	TH1D* h_conv_dxyError_Tr2 = new TH1D("h_conv_dxyError_Tr2", "h_conv_dxyError_Tr2", 200, 0, 10);

	TH1I* h_conv_tk1NumOfDOF = new TH1I("h_conv_tk1NumOfDOF", "h_conv_tk1NumOfDOF", 70, 0, 70);
	TH1I* h_conv_tk2NumOfDOF = new TH1I("h_conv_tk2NumOfDOF", "h_conv_tk2NumOfDOF", 70, 0, 70);
	TH1D* h_conv_track1Chi2 = new TH1D("h_conv_track1Chi2", "h_conv_track1Chi2", 400, 0, 40);
	TH1D* h_conv_track2Chi2 = new TH1D("h_conv_track2Chi2", "h_conv_track2Chi2", 400, 0, 40);
	TH1D* h_conv_minDistanceOfApproach = new TH1D("h_conv_minDistanceOfApproach", "h_conv_minDistanceOfApproach", 400, -2, 2);
	TH1D* h_conv_eta = new TH1D("h_conv_eta", "h_conv_eta", 200, -5, 5);
	TH1D* h_conv_pt = new TH1D("h_conv_pt", "h_conv_pt", 600, 0, 30);

	TH1D* h_gen_phot_pt = new TH1D("h_gen_phot_pt", "h_gen_phot_pt", 600, 0, 30);
	TH1D* h_gen_phot_eta = new TH1D("h_gen_phot_eta", "h_gen_phot_eta", 200, -5, 5);
	TH1I* h_gen_conv_matchPosition = new TH1I("h_gen_conv_matchPosition", "h_gen_conv_matchPosition", 10, 0, 10);
	TH1I* h_gen_conv_nMatches = new TH1I("h_gen_conv_nMatches", "h_gen_conv_nMatches", 10, 0, 10);
	TH1D* h_gen_conv_rDelta = new TH1D("h_gen_conv_rDelta", "h_gen_conv_rDelta", 100, 0, 0.5);
	TH1D* h_gen_conv_ptDeltaRel = new TH1D("h_gen_conv_ptDeltaRel", "h_gen_conv_ptDeltaRel", 200, -1, 1);



	// ACCEPTANCE
	TH2D* h_muonAcceptance2D_den = new TH2D("h_muonAcceptance2D_den", "h_muonAcceptance2D_den; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_muonAcceptance2D_num = new TH2D("h_muonAcceptance2D_num", "h_muonAcceptance2D_num; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_muonAcceptance2D_rat = new TH2D("h_muonAcceptance2D_rat", "h_muonAcceptance2D_rat; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_photAcceptance2D_den = new TH2D("h_photAcceptance2D_den", "h_photAcceptance2D_den; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_photAcceptance2D_num = new TH2D("h_photAcceptance2D_num", "h_photAcceptance2D_num; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_photAcceptance2D_rat = new TH2D("h_photAcceptance2D_rat", "h_photAcceptance2D_rat; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_JpsiAcceptance2D_den = new TH2D("h_JpsiAcceptance2D_den", "h_JpsiAcceptance2D_den; #eta; p_{T}", 60, -6, 6, 100, 0, 20);
	TH2D* h_JpsiAcceptance2D_num = new TH2D("h_JpsiAcceptance2D_num", "h_JpsiAcceptance2D_num; #eta; p_{T}", 60, -6, 6, 100, 0, 20);
	TH2D* h_JpsiAcceptance2D_rat = new TH2D("h_JpsiAcceptance2D_rat", "h_JpsiAcceptance2D_rat; #eta; p_{T}", 60, -6, 6, 100, 0, 20);
	TH2D* h_JpsiAcceptance2D_y_den = new TH2D("h_JpsiAcceptance2D_y_den", "h_JpsiAcceptance2D_y_den; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_JpsiAcceptance2D_y_num = new TH2D("h_JpsiAcceptance2D_y_num", "h_JpsiAcceptance2D_y_num; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_JpsiAcceptance2D_y_rat = new TH2D("h_JpsiAcceptance2D_y_rat", "h_JpsiAcceptance2D_y_rat; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_chiAcceptance2D_den = new TH2D("h_chiAcceptance2D_den", "h_chiAcceptance2D_den; #eta; p_{T}", 60, -6, 6, 100, 0, 20);
	TH2D* h_chiAcceptance2D_num = new TH2D("h_chiAcceptance2D_num", "h_chiAcceptance2D_num; #eta; p_{T}", 60, -6, 6, 100, 0, 20);
	TH2D* h_chiAcceptance2D_rat = new TH2D("h_chiAcceptance2D_rat", "h_chiAcceptance2D_rat; #eta; p_{T}", 60, -6, 6, 100, 0, 20);
	TH2D* h_chiAcceptance2D_y_den = new TH2D("h_chiAcceptance2D_y_den", "h_chiAcceptance2D_y_den; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_chiAcceptance2D_y_num = new TH2D("h_chiAcceptance2D_y_num", "h_chiAcceptance2D_y_num; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_chiAcceptance2D_y_rat = new TH2D("h_chiAcceptance2D_y_rat", "h_chiAcceptance2D_y_rat; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	// EFFICIENCY
	TH2D* h_muonEfficiency2D_den = new TH2D("h_muonEfficiency2D_den", "h_muonEfficiency2D_den; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_muonEfficiency2D_num = new TH2D("h_muonEfficiency2D_num", "h_muonEfficiency2D_num; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_muonEfficiency2D_rat = new TH2D("h_muonEfficiency2D_rat", "h_muonEfficiency2D_rat; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_photEfficiency2D_den = new TH2D("h_photEfficiency2D_den", "h_photEfficiency2D_den; #eta; p_{T}", 10, -3.0, 3.0, nbinsConvEffpT, binsConvEffpT);
	TH2D* h_photEfficiency2D_num = new TH2D("h_photEfficiency2D_num", "h_photEfficiency2D_num; #eta; p_{T}", 10, -3.0, 3.0, nbinsConvEffpT, binsConvEffpT);
	TH2D* h_photEfficiency2D_rat = new TH2D("h_photEfficiency2D_rat", "h_photEfficiency2D_rat; #eta; p_{T}", 10, -3.0, 3.0, nbinsConvEffpT, binsConvEffpT);
	TH2D* h_JpsiEfficiency2D_y_den = new TH2D("h_JpsiEfficiency2D_y_den", "h_JpsiEfficiency2D_y_den; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_JpsiEfficiency2D_y_num = new TH2D("h_JpsiEfficiency2D_y_num", "h_JpsiEfficiency2D_y_num; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_JpsiEfficiency2D_y_rat = new TH2D("h_JpsiEfficiency2D_y_rat", "h_JpsiEfficiency2D_y_rat; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_chiEfficiency2D_y_den = new TH2D("h_chiEfficiency2D_y_den", "h_chiEfficiency2D_y_den; y_{lab}; p_{T}", 10, -3.0, 3.0, nbinsChiEffpT, binsChiEffpT);
	TH2D* h_chiEfficiency2D_y_num = new TH2D("h_chiEfficiency2D_y_num", "h_chiEfficiency2D_y_num; y_{lab}; p_{T}", 10, -3.0, 3.0, nbinsChiEffpT, binsChiEffpT);
	TH2D* h_chiEfficiency2D_y_rat = new TH2D("h_chiEfficiency2D_y_rat", "h_chiEfficiency2D_y_rat; y_{lab}; p_{T}", 10, -3.0, 3.0, nbinsChiEffpT, binsChiEffpT);

	TH1D* h_photEfficiency1D_den = new TH1D("h_photEfficiency1D_den", "h_photEfficiency1D_den; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_photEfficiency1D_num = new TH1D("h_photEfficiency1D_num", "h_photEfficiency1D_num; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_photEfficiency1D_rat = new TH1D("h_photEfficiency1D_rat", "h_photEfficiency1D_rat; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_photEfficiency1D_relative_den = new TH1D("h_photEfficiency1D_relative_den", "h_photEfficiency1D_relative_den; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_photEfficiency1D_relative_num = new TH1D("h_photEfficiency1D_relative_num", "h_photEfficiency1D_relative_num; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_photEfficiency1D_relative_rat = new TH1D("h_photEfficiency1D_relative_rat", "h_photEfficiency1D_relative_rat; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_JpsiEfficiency1D_den = new TH1D("h_JpsiEfficiency1D_den", "h_JpsiEfficiency1D_den; p_{T}", 100, 0, 20);
	TH1D* h_JpsiEfficiency1D_num = new TH1D("h_JpsiEfficiency1D_num", "h_JpsiEfficiency1D_num; p_{T}", 100, 0, 20);
	TH1D* h_JpsiEfficiency1D_rat = new TH1D("h_JpsiEfficiency1D_rat", "h_JpsiEfficiency1D_rat; p_{T}", 100, 0, 20);
	TH1D* h_chiEfficiency1D_den = new TH1D("h_chiEfficiency1D_den", "h_chiEfficiency1D_den; p_{T}", nbinsChiEffpT, binsChiEffpT);
	TH1D* h_chiEfficiency1D_num = new TH1D("h_chiEfficiency1D_num", "h_chiEfficiency1D_num; p_{T}", nbinsChiEffpT, binsChiEffpT);
	TH1D* h_chiEfficiency1D_rat = new TH1D("h_chiEfficiency1D_rat", "h_chiEfficiency1D_rat; p_{T}", nbinsChiEffpT, binsChiEffpT);


	TFile* f1 = new TFile(fileIn, "READ");

	//double Mdiff, dimuonVertexP, EHF_trans;
	//TLorentzVector* LVchic, *LVdimuon, *LVphoton, *LVmuonP, *LVmuonN;
	//TVector3* V3primaryV, *V3secondaryV;
	//int trig;




	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");
	if (!event_tree) { 
		cout << "Problem with event Tree";
		return 1;
	}


	event_tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal);
	event_tree->SetBranchAddress("muonIsTracker", &muonIsTracker);
	event_tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
	event_tree->SetBranchAddress("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);

	event_tree->SetBranchAddress("muon_eta", &muon_eta);
	event_tree->SetBranchAddress("muon_pt", &muon_pt);
	
	
	
	// MC

	event_tree->SetBranchAddress("gen_pdgId", &gen_pdgId);
	event_tree->SetBranchAddress("gen_chic_pt", &gen_chic_pt);
	event_tree->SetBranchAddress("gen_chic_eta", &gen_chic_eta);
	event_tree->SetBranchAddress("gen_chic_p4", &gen_chic_p4);
	event_tree->SetBranchAddress("gen_Jpsi_pt", &gen_Jpsi_pt);
	event_tree->SetBranchAddress("gen_Jpsi_eta", &gen_Jpsi_eta);
	event_tree->SetBranchAddress("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
	event_tree->SetBranchAddress("gen_Jpsi_nMatches", &gen_Jpsi_nMatches);
	event_tree->SetBranchAddress("gen_Jpsi_rDelta", &gen_Jpsi_rDelta);
	event_tree->SetBranchAddress("gen_Jpsi_ptDeltaRel", &gen_Jpsi_ptDeltaRel);
	event_tree->SetBranchAddress("gen_Jpsi_p4", &gen_Jpsi_p4);

	event_tree->SetBranchAddress("gen_muon_nMatches", &gen_muon_nMatches);
	event_tree->SetBranchAddress("gen_muon_matchPosition", &gen_muon_matchPosition);
	event_tree->SetBranchAddress("gen_muon_eta", &gen_muon_eta);
	event_tree->SetBranchAddress("gen_muon_pt", &gen_muon_pt);
	event_tree->SetBranchAddress("gen_muon_rDelta", &gen_muon_rDelta);
	event_tree->SetBranchAddress("gen_muon_ptDeltaRel", &gen_muon_ptDeltaRel);

	event_tree->SetBranchAddress("gen_phot_pt", &gen_phot_pt);
	event_tree->SetBranchAddress("gen_phot_eta", &gen_phot_eta);
	event_tree->SetBranchAddress("gen_conv_matchPosition", &gen_conv_matchPosition);
	event_tree->SetBranchAddress("gen_conv_nMatches", &gen_conv_nMatches);
	event_tree->SetBranchAddress("gen_conv_rDelta", &gen_conv_rDelta);
	event_tree->SetBranchAddress("gen_conv_ptDeltaRel", &gen_conv_ptDeltaRel);




	//Conversions
	event_tree->SetBranchAddress("convQuality_isHighPurity", &convQuality_isHighPurity);
	event_tree->SetBranchAddress("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly);
	event_tree->SetBranchAddress("conv_vertexPositionRho", &conv_vertexPositionRho);
	event_tree->SetBranchAddress("conv_sigmaTkVtx1", &conv_sigmaTkVtx1);
	event_tree->SetBranchAddress("conv_sigmaTkVtx2", &conv_sigmaTkVtx2);
	event_tree->SetBranchAddress("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOK);

	event_tree->SetBranchAddress("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK);
	event_tree->SetBranchAddress("conv_vertexChi2Prob", &conv_vertexChi2Prob);
	event_tree->SetBranchAddress("conv_zOfPriVtx", &conv_zOfPriVtx);
	event_tree->SetBranchAddress("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks);
	event_tree->SetBranchAddress("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx);
	event_tree->SetBranchAddress("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1);
	event_tree->SetBranchAddress("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2);
	event_tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1);
	event_tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2);
	event_tree->SetBranchAddress("conv_dxyError_Tr1", &conv_dxyError_Tr1);
	event_tree->SetBranchAddress("conv_dxyError_Tr2", &conv_dxyError_Tr2);

	event_tree->SetBranchAddress("conv_tk1NumOfDOF", &conv_tk1NumOfDOF);
	event_tree->SetBranchAddress("conv_tk2NumOfDOF", &conv_tk2NumOfDOF);
	event_tree->SetBranchAddress("conv_track1Chi2", &conv_track1Chi2);
	event_tree->SetBranchAddress("conv_track2Chi2", &conv_track2Chi2);
	event_tree->SetBranchAddress("conv_minDistanceOfApproach", &conv_minDistanceOfApproach);
	event_tree->SetBranchAddress("conv_eta", &conv_eta);
	event_tree->SetBranchAddress("conv_pt", &conv_pt);


	// Chic
	event_tree->SetBranchAddress("chi_eta", &chi_eta);
	event_tree->SetBranchAddress("chi_pt", &chi_pt);
	event_tree->SetBranchAddress("chi_daughterJpsi_position", &chi_daughterJpsi_position);
	event_tree->SetBranchAddress("chi_daughterConv_position", &chi_daughterConv_position);
	event_tree->SetBranchAddress("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx);
	event_tree->SetBranchAddress("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx);


	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: "<<nentries << endl;
	//if (nentries > 5000) { nentries = 5000; }


	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);

		//cout << "here" << endl;
		if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }

		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1)
		{
			weird_decay_counter++;
			continue;
		}

		//GEN muons // 2 per event
		for (Long64_t j = 0; j < 2; j++)
		{	
			int matchPosition = gen_muon_matchPosition->at(j);
			if (matchPosition < -0.5) { continue; } //negative means not matched
			double muonPtRel1 = gen_muon_ptDeltaRel->at(j);
			double muonPtRel2 = fabs(muon_pt->at(matchPosition) - gen_muon_pt->at(j)) / gen_muon_pt->at(j);
			
			h_gen_muon_matchPosition->Fill(matchPosition);
			h_gen_muon_nMatches->Fill(gen_muon_nMatches->at(j));
			h_gen_muon_rDelta->Fill(gen_muon_rDelta->at(j));

			h_muon_ptRel1->Fill(muonPtRel1);
			h_muon_ptRel2->Fill(muonPtRel2);
			
			bool muonTracker = muonIsTracker->at(matchPosition);
			bool muonGlobal = muonIsGlobal->at(matchPosition);
			h_muonMatch_isGlobal->Fill(muonGlobal);
			h_muonMatch_isTracker->Fill(muonTracker);
			h_muonMatch_isSoft->Fill(muonIsSoft->at(matchPosition));
			


			//if (!muonGlobal || !muonTracker) { cout << muonGlobal << "  " << muonTracker << endl; }
			if (muonGlobal || muonTracker) {
				h_muonMatch_trackerLayers->Fill(muonTrackerLayersWithMeasurement->at(matchPosition));
			}
			//cout << "gen eta: " << gen_muon_eta->at(j) << "  and muon eta: " << muon_eta->at(matchPosition) <<"   and rDelta: " << gen_muon_rDelta->at(j)<< endl;
			//cout << "gen pt: " << gen_muon_pt->at(j) << "  and muon pt: " << muon_pt->at(matchPosition) << "   and ptDelta: " << gen_muon_ptDeltaRel->at(j) << endl<<endl;
		}
		// Conversions from chic //one per event
		int conv_matchPosition = gen_conv_matchPosition->at(0);
		if (conv_matchPosition > -0.5) { //non-negative means it was matched
			h_convQuality_isHighPurity->Fill(convQuality_isHighPurity->at(conv_matchPosition));
			h_convQuality_isGeneralTracksOnly->Fill(convQuality_isGeneralTracksOnly->at(conv_matchPosition));
			h_conv_vertexPositionRho->Fill(conv_vertexPositionRho->at(conv_matchPosition));
			h_conv_sigmaTkVtx1->Fill(conv_sigmaTkVtx1->at(conv_matchPosition));
			h_conv_sigmaTkVtx2->Fill(conv_sigmaTkVtx2->at(conv_matchPosition));
			h_conv_tkVtxCompatibilityOK->Fill(conv_tkVtxCompatibilityOK->at(conv_matchPosition));

			h_conv_compatibleInnerHitsOK->Fill(conv_compatibleInnerHitsOK->at(conv_matchPosition)); //-1: less than 2 tracks, 0: not compatible, 1: yes
			h_conv_vertexChi2Prob->Fill(conv_vertexChi2Prob->at(conv_matchPosition));
			h_conv_zOfPriVtx->Fill(conv_zOfPriVtx->at(conv_matchPosition));
			h_conv_zOfPriVtxFromTracks->Fill(conv_zOfPriVtxFromTracks->at(conv_matchPosition));
			h_conv_dzToClosestPriVtx->Fill(conv_dzToClosestPriVtx->at(conv_matchPosition));
			h_conv_dxyPriVtx_Tr1->Fill(conv_dxyPriVtx_Tr1->at(conv_matchPosition));
			h_conv_dxyPriVtx_Tr2->Fill(conv_dxyPriVtx_Tr2->at(conv_matchPosition));
			h_conv_dxyPriVtxTimesCharge_Tr1->Fill(conv_dxyPriVtxTimesCharge_Tr1->at(conv_matchPosition));
			h_conv_dxyPriVtxTimesCharge_Tr2->Fill(conv_dxyPriVtxTimesCharge_Tr2->at(conv_matchPosition));
			h_conv_dxyError_Tr1->Fill(conv_dxyError_Tr1->at(conv_matchPosition));
			h_conv_dxyError_Tr2->Fill(conv_dxyError_Tr2->at(conv_matchPosition));

			h_conv_tk1NumOfDOF->Fill(conv_tk1NumOfDOF->at(conv_matchPosition));
			h_conv_tk2NumOfDOF->Fill(conv_tk2NumOfDOF->at(conv_matchPosition));
			h_conv_track1Chi2->Fill(conv_track1Chi2->at(conv_matchPosition));
			h_conv_track2Chi2->Fill(conv_track2Chi2->at(conv_matchPosition));
			h_conv_minDistanceOfApproach->Fill(conv_minDistanceOfApproach->at(conv_matchPosition));
			h_conv_eta->Fill(conv_eta->at(conv_matchPosition));
			h_conv_pt->Fill(conv_pt->at(conv_matchPosition));

			h_gen_phot_pt->Fill(gen_phot_pt->at(0));
			h_gen_phot_eta->Fill(gen_phot_eta->at(0));
			h_gen_conv_matchPosition->Fill(gen_conv_matchPosition->at(0));
			h_gen_conv_nMatches->Fill(gen_conv_nMatches->at(0));
			h_gen_conv_rDelta->Fill(gen_conv_rDelta->at(0));
			h_gen_conv_ptDeltaRel->Fill(gen_conv_ptDeltaRel->at(0));

		}
		//Chic part
		// ACCEPTANCE
		double eta1 = gen_muon_eta->at(0);
		double pt1 = gen_muon_pt->at(0);
		double eta2 = gen_muon_eta->at(1);
		double pt2 = gen_muon_pt->at(1);
		h_muonAcceptance2D_den->Fill(eta1, pt1);
		h_muonAcceptance2D_den->Fill(eta2, pt2);
		if (MuonAcceptance(eta1, pt1) == true) h_muonAcceptance2D_num->Fill(eta1, pt1);
		if (MuonAcceptance(eta2, pt2) == true) h_muonAcceptance2D_num->Fill(eta2, pt2);
		double eta_phot = gen_phot_eta->at(0);
		double pt_phot = gen_phot_pt->at(0);
		h_photAcceptance2D_den->Fill(eta_phot, pt_phot);
		if (PhotAcceptance(eta_phot, pt_phot) == true) h_photAcceptance2D_num->Fill(eta_phot, pt_phot);
		double eta_Jpsi = gen_Jpsi_eta->at(0);
		double pt_Jpsi = gen_Jpsi_pt->at(0);
		h_JpsiAcceptance2D_den->Fill(eta_Jpsi, pt_Jpsi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true) h_JpsiAcceptance2D_num->Fill(eta_Jpsi, pt_Jpsi);
		TLorentzVector* LVJpsi = (TLorentzVector*)gen_Jpsi_p4->At(0);
		h_JpsiAcceptance2D_y_den->Fill(LVJpsi->Rapidity(), pt_Jpsi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) h_JpsiAcceptance2D_y_num->Fill(LVJpsi->Rapidity(), pt_Jpsi);

		double eta_chi = gen_chic_eta->at(0);
		double pt_chi = gen_chic_pt->at(0);
		h_chiAcceptance2D_den->Fill(eta_chi, pt_chi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot)==true) h_chiAcceptance2D_num->Fill(eta_chi, pt_chi);
		TLorentzVector* LVchic = (TLorentzVector*)gen_chic_p4->At(0);
		h_chiAcceptance2D_y_den->Fill(LVchic->Rapidity(), pt_chi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) h_chiAcceptance2D_y_num->Fill(LVchic->Rapidity(), pt_chi);

		// EFFICIENCY
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) //do only if everything in acceptance
		{
			h_muonEfficiency2D_den->Fill(eta1, pt1);
			h_muonEfficiency2D_den->Fill(eta2, pt2);
			if (MuonSelectionPassMC(0) == true) h_muonEfficiency2D_num->Fill(eta1, pt1);
			if (MuonSelectionPassMC(1) == true) h_muonEfficiency2D_num->Fill(eta2, pt2);

			h_photEfficiency2D_den->Fill(eta_phot, pt_phot);
			h_photEfficiency1D_den->Fill(pt_phot);
			if (PhotSelectionPassMC(0) == true) {
				h_photEfficiency2D_num->Fill(eta_phot, pt_phot);
				h_photEfficiency1D_num->Fill(pt_phot);
				h_photEfficiency1D_relative_num->Fill(pt_phot);

			}
			if (PhotSelectionPassMCLoose(0) == true) {
				h_photEfficiency1D_relative_den->Fill(pt_phot);
			}

			h_JpsiEfficiency2D_y_den->Fill(LVJpsi->Rapidity(), pt_Jpsi);
			h_JpsiEfficiency1D_den->Fill(pt_Jpsi);
			if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true) {
				h_JpsiEfficiency2D_y_num->Fill(LVJpsi->Rapidity(), pt_Jpsi);
				h_JpsiEfficiency1D_num->Fill(pt_Jpsi);
			}
			h_chiEfficiency2D_y_den->Fill(LVchic->Rapidity(), pt_chi);
			h_chiEfficiency1D_den->Fill(pt_chi);
			if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && PhotSelectionPassMC(0) == true) {
				h_chiEfficiency2D_y_num->Fill(LVchic->Rapidity(), pt_chi);
				h_chiEfficiency1D_num->Fill(pt_chi);
			}
		}






	} //end of event loop

	cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;


	h_muonAcceptance2D_rat->Divide(h_muonAcceptance2D_num, h_muonAcceptance2D_den);
	h_photAcceptance2D_rat->Divide(h_photAcceptance2D_num, h_photAcceptance2D_den);
	h_JpsiAcceptance2D_rat->Divide(h_JpsiAcceptance2D_num, h_JpsiAcceptance2D_den);
	h_JpsiAcceptance2D_y_rat->Divide(h_JpsiAcceptance2D_y_num, h_JpsiAcceptance2D_y_den);
	h_chiAcceptance2D_rat->Divide(h_chiAcceptance2D_num, h_chiAcceptance2D_den);
	h_chiAcceptance2D_y_rat->Divide(h_chiAcceptance2D_y_num, h_chiAcceptance2D_y_den);

	h_muonEfficiency2D_rat->Divide(h_muonEfficiency2D_num, h_muonEfficiency2D_den);
	h_photEfficiency2D_rat->Divide(h_photEfficiency2D_num, h_photEfficiency2D_den);
	h_JpsiEfficiency2D_y_rat->Divide(h_JpsiEfficiency2D_y_num, h_JpsiEfficiency2D_y_den);
	h_chiEfficiency2D_y_rat->Divide(h_chiEfficiency2D_y_num, h_chiEfficiency2D_y_den);

	h_photEfficiency1D_num->Sumw2();
	h_photEfficiency1D_den->Sumw2();
	h_photEfficiency1D_relative_num->Sumw2();
	h_photEfficiency1D_relative_den->Sumw2();
	h_JpsiEfficiency1D_num->Sumw2();
	h_JpsiEfficiency1D_den->Sumw2();
	h_chiEfficiency1D_num->Sumw2();
	h_chiEfficiency1D_den->Sumw2();

	h_photEfficiency1D_rat->Divide(h_photEfficiency1D_num, h_photEfficiency1D_den, 1, 1, "B");
	h_photEfficiency1D_relative_rat->Divide(h_photEfficiency1D_relative_num, h_photEfficiency1D_relative_den, 1, 1, "B");
	h_JpsiEfficiency1D_rat->Divide(h_JpsiEfficiency1D_num, h_JpsiEfficiency1D_den, 1, 1, "B");
	h_chiEfficiency1D_rat->Divide(h_chiEfficiency1D_num, h_chiEfficiency1D_den, 1, 1, "B");

	
	//TCanvas* can1 = new TCanvas("can1", "plot", 1200, 800);
	//h_muonAcceptance2D_num->Draw("colz");
	//TF1* f_filter = new TF1("f_filter", "2.5/TMath::CosH(x)", -2.5, 2.5);
	//f_filter->SetLineColor(kBlue);
	//f_filter->Draw("same");
	//TF1* f_acceptance = new TF1("f_acceptance","(fabs(x) > 2.4)*100+((fabs(x) < 0.8)*3.3 + (fabs(x)>=0.8)*(fabs(x)<1.5)*(5.81-3.14*fabs(x))+(fabs(x)>=1.5)*(fabs(x)<2.07)*(1.89-0.526*fabs(x))+(fabs(x)>=2.07)*0.8)",-2.5,2.5);
	//f_acceptance->Draw("same");
	//can1->SaveAs("muonAcceptance2D_num.png");
	
	TFile* fout = new TFile(fileOut, "RECREATE");

	h_muon_ptRel1->Write();
	h_muon_ptRel2->Write();

	h_gen_muon_matchPosition->Write();
	h_gen_muon_nMatches->Write();
	h_gen_muon_rDelta->Write();

	h_muonMatch_isGlobal->Write();
	h_muonMatch_isTracker->Write();
	h_muonMatch_isSoft->Write();
	h_muonMatch_trackerLayers->Write();



	//conv

	h_convQuality_isHighPurity->Write();
	h_convQuality_isGeneralTracksOnly->Write();
	h_conv_vertexPositionRho->Write();
	h_conv_sigmaTkVtx1->Write();
	h_conv_sigmaTkVtx2->Write();
	h_conv_tkVtxCompatibilityOK->Write();

	h_conv_compatibleInnerHitsOK->Write(); //-1: less than 2 tracks, 0: not compatible, 1: yes
	h_conv_vertexChi2Prob->Write();
	h_conv_zOfPriVtx->Write();
	h_conv_zOfPriVtxFromTracks->Write();
	h_conv_dzToClosestPriVtx->Write();
	h_conv_dxyPriVtx_Tr1->Write();
	h_conv_dxyPriVtx_Tr2->Write();
	h_conv_dxyPriVtxTimesCharge_Tr1->Write();
	h_conv_dxyPriVtxTimesCharge_Tr2->Write();
	h_conv_dxyError_Tr1->Write();
	h_conv_dxyError_Tr2->Write();

	h_conv_tk1NumOfDOF->Write();
	h_conv_tk2NumOfDOF->Write();
	h_conv_track1Chi2->Write();
	h_conv_track2Chi2->Write();
	h_conv_minDistanceOfApproach->Write();
	h_conv_eta->Write();
	h_conv_pt->Write();

	h_gen_phot_pt->Write();
	h_gen_phot_eta->Write();
	h_gen_conv_matchPosition->Write();
	h_gen_conv_nMatches->Write();
	h_gen_conv_rDelta->Write();
	h_gen_conv_ptDeltaRel->Write();

	h_muonAcceptance2D_den->Write();
	h_muonAcceptance2D_num->Write();
	h_muonAcceptance2D_rat->Write();
	h_photAcceptance2D_den->Write();
	h_photAcceptance2D_num->Write();
	h_photAcceptance2D_rat->Write();
	h_JpsiAcceptance2D_den->Write();
	h_JpsiAcceptance2D_num->Write();
	h_JpsiAcceptance2D_rat->Write();
	h_JpsiAcceptance2D_y_den->Write();
	h_JpsiAcceptance2D_y_num->Write();
	h_JpsiAcceptance2D_y_rat->Write();
	h_chiAcceptance2D_den->Write();
	h_chiAcceptance2D_num->Write();
	h_chiAcceptance2D_rat->Write();
	h_chiAcceptance2D_y_den->Write();
	h_chiAcceptance2D_y_num->Write();
	h_chiAcceptance2D_y_rat->Write();

	h_muonEfficiency2D_den->Write();
	h_muonEfficiency2D_num->Write();
	h_muonEfficiency2D_rat->Write();
	h_photEfficiency2D_den->Write();
	h_photEfficiency2D_num->Write();
	h_photEfficiency2D_rat->Write();
	h_JpsiEfficiency2D_y_den->Write();
	h_JpsiEfficiency2D_y_num->Write();
	h_JpsiEfficiency2D_y_rat->Write();
	h_chiEfficiency2D_y_den->Write();
	h_chiEfficiency2D_y_num->Write();
	h_chiEfficiency2D_y_rat->Write();
	h_photEfficiency1D_den->Write();
	h_photEfficiency1D_num->Write();
	h_photEfficiency1D_rat->Write();
	h_photEfficiency1D_relative_den->Write();
	h_photEfficiency1D_relative_num->Write();
	h_photEfficiency1D_relative_rat->Write();
	h_JpsiEfficiency1D_den->Write();
	h_JpsiEfficiency1D_num->Write();
	h_JpsiEfficiency1D_rat->Write();
	h_chiEfficiency1D_den->Write();
	h_chiEfficiency1D_num->Write();
	h_chiEfficiency1D_rat->Write();

	fout->Close();
	f1->Close();
	cout << "END" << endl;
	return 0;
}
