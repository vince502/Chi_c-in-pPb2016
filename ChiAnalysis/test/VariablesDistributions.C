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
TLorentzVector* LVchic, *LVJpsi, *LVconv, *LVmuon1, *LVmuon2, *LVJpsi_phot;

int Trig_Event_HLTDoubleMuOpen;

//variables
std::vector <bool>* muonIsHLTDoubleMuOpen = 0;
std::vector <bool>* muonIsGlobal = 0;
std::vector <bool>* muonIsTracker = 0;
std::vector <int>* muonTrackerLayersWithMeasurement = 0;
std::vector <bool>* muonIsSoft = 0;
std::vector <double>* muon_eta = 0;
std::vector <double>* muon_pt = 0;

// vertex
std::vector <double>* pvtx_z = 0;
std::vector <double>* pvtx_zError = 0;
std::vector <double>* pvtx_x = 0;
std::vector <double>* pvtx_y = 0;
std::vector <double>* pvtx_nTracks = 0;
std::vector <bool>* pvtx_isFake = 0;


//dimuon
TClonesArray*  dimuon_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <double>* dimuon_eta = 0;
std::vector <double>* dimuon_pt = 0;
std::vector <double>* dimuon_charge = 0;
std::vector <int>*  dimuon_pvtx_index = 0;
std::vector <double>* dimuon_dz_dimuonvtx_pvtx = 0;
std::vector <double>* dimuon_vtxProb = 0;
std::vector <int>* dimuon_muon1_position = 0; //stores position of first muon in muon collection 
std::vector <int>* dimuon_muon2_position = 0; //stores position of second muon in muon collection 
std::vector <double>* dimuon_ctpv = 0;
std::vector <double>* dimuon_ctpvError = 0;


// MC general
std::vector <int>* gen_pdgId = 0;
std::vector <double>* gen_chic_pt = 0;
std::vector <double>* gen_chic_eta = 0;
std::vector <int>* gen_chic_matchPosition = 0;
std::vector <int>* gen_chic_nMatches = 0;
std::vector <double>* gen_chic_rDelta = 0; //in principle duplicates information
std::vector <double>* gen_chic_ptDeltaRel = 0;//in principle duplicates information
TClonesArray* gen_chic_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <double>* gen_Jpsi_pt = 0;
std::vector <double>* gen_Jpsi_eta = 0;
std::vector <int>* gen_Jpsi_matchPosition = 0;
std::vector <int>* gen_Jpsi_nMatches = 0;
std::vector <double>* gen_Jpsi_rDelta = 0; //in principle duplicates information
std::vector <double>* gen_Jpsi_ptDeltaRel = 0;//in principle duplicates information
TClonesArray* gen_Jpsi_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <int>* gen_Jpsi_photon_n = 0;
std::vector <double>* gen_Jpsi_photon_pt = 0;
TClonesArray* gen_Jpsi_photon_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <int>* gen_muon_charge = 0;
std::vector <int>* gen_muon_nMatches = 0;
std::vector <int>* gen_muon_matchPosition = 0;
std::vector <double>* gen_muon_eta = 0;
std::vector <double>* gen_muon_pt = 0;
std::vector <double>* gen_muon_rDelta = 0;
std::vector <double>* gen_muon_ptDeltaRel = 0;
TClonesArray* gen_muon_p4 = new TClonesArray("TLorentzVector", 100);
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
std::vector <int>* chi_kinematicRefitFlag = 0; // -1 kinematic refit not done, 1 done: +2 needed extra refit for photon +4 something wrong with photon at the end +8 something wrong with the final fit 
std::vector <int>* chi_refit_origChicPosition = 0;
std::vector <double>* chi_refit_vprob = 0;
std::vector <double>* chi_refit_ctauPV = 0;
std::vector <double>* chi_refit_ctauErrPV = 0;
std::vector <double>* chi_refit_ctauPV3D = 0;
std::vector <double>* chi_refit_pvtxFromPVwithMuons_x = 0;
std::vector <double>* chi_refit_pvtxFromPVwithMuons_y = 0;
std::vector <double>* chi_refit_pvtxFromPVwithMuons_z = 0;



///////////////////////
//////   D A T A   ////
//////////////////////////

int Trig_Event_HLTDoubleMuOpenRD;

//variables
std::vector <bool>* muonIsHLTDoubleMuOpenRD = 0;
std::vector <bool>* muonIsGlobalRD = 0;
std::vector <bool>* muonIsTrackerRD = 0;
std::vector <int>* muonTrackerLayersWithMeasurementRD = 0;
std::vector <bool>* muonIsSoftRD = 0;
std::vector <double>* muon_etaRD = 0;
std::vector <double>* muon_ptRD = 0;

//dimuon
TClonesArray*  dimuon_p4RD = new TClonesArray("TLorentzVector", 100);
std::vector <double>* dimuon_etaRD = 0;
std::vector <double>* dimuon_ptRD = 0;
std::vector <double>* dimuon_chargeRD = 0;
std::vector <int>*  dimuon_pvtx_indexRD = 0;
std::vector <double>* dimuon_dz_dimuonvtx_pvtxRD = 0;
std::vector <double>* dimuon_vtxProbRD = 0;
std::vector <int>* dimuon_muon1_positionRD = 0; //stores position of first muon in muon collection 
std::vector <int>* dimuon_muon2_positionRD = 0; //stores position of second muon in muon collection 
std::vector <double>* dimuon_ctpvRD = 0;
std::vector <double>* dimuon_ctpvErrorRD = 0;


// MC general
std::vector <int>* gen_pdgIdRD = 0;
std::vector <double>* gen_chic_ptRD = 0;
std::vector <double>* gen_chic_etaRD = 0;
TClonesArray* gen_chic_p4RD = new TClonesArray("TLorentzVector", 100);
std::vector <double>* gen_Jpsi_ptRD = 0;
std::vector <double>* gen_Jpsi_etaRD = 0;
std::vector <int>* gen_Jpsi_matchPositionRD = 0;
std::vector <int>* gen_Jpsi_nMatchesRD = 0;
std::vector <double>* gen_Jpsi_rDeltaRD = 0; //in principle duplicates information
std::vector <double>* gen_Jpsi_ptDeltaRelRD = 0;//in principle duplicates information
TClonesArray* gen_Jpsi_p4RD = new TClonesArray("TLorentzVector", 100);
std::vector <int>* gen_muon_chargeRD = 0;
std::vector <int>* gen_muon_nMatchesRD = 0;
std::vector <int>* gen_muon_matchPositionRD = 0;
std::vector <double>* gen_muon_etaRD = 0;
std::vector <double>* gen_muon_ptRD = 0;
std::vector <double>* gen_muon_rDeltaRD = 0;
std::vector <double>* gen_muon_ptDeltaRelRD = 0;
std::vector <double>* gen_phot_ptRD = 0;
std::vector <double>* gen_phot_etaRD = 0;
std::vector <int>* gen_conv_matchPositionRD = 0;
std::vector <int>* gen_conv_nMatchesRD = 0;
std::vector <double>* gen_conv_rDeltaRD = 0;
std::vector <double>* gen_conv_ptDeltaRelRD = 0;


//// Conversions  /////
std::vector <bool>* convQuality_isHighPurityRD = 0;
std::vector <bool>* convQuality_isGeneralTracksOnlyRD = 0;
std::vector <double>* conv_vertexPositionRhoRD = 0;
std::vector <double>* conv_sigmaTkVtx1RD = 0;
std::vector <double>* conv_sigmaTkVtx2RD = 0;
std::vector <bool>* conv_tkVtxCompatibilityOKRD = 0;
std::vector <int>* conv_compatibleInnerHitsOKRD = 0; //-1: less than 2 tracks, 0: not compatible, 1: yes
std::vector <double>* conv_vertexChi2ProbRD = 0;
std::vector <double>* conv_zOfPriVtxRD = 0;
std::vector <double>* conv_zOfPriVtxFromTracksRD = 0;
std::vector <double>* conv_dzToClosestPriVtxRD = 0;
std::vector <double>* conv_dxyPriVtx_Tr1RD = 0;
std::vector <double>* conv_dxyPriVtx_Tr2RD = 0;
std::vector <double>* conv_dxyPriVtxTimesCharge_Tr1RD = 0;
std::vector <double>* conv_dxyPriVtxTimesCharge_Tr2RD = 0;
std::vector <double>* conv_dxyError_Tr1RD = 0;
std::vector <double>* conv_dxyError_Tr2RD = 0;
std::vector <int>* conv_tk1NumOfDOFRD = 0;
std::vector <int>* conv_tk2NumOfDOFRD = 0;
std::vector <double>* conv_track1Chi2RD = 0;
std::vector <double>* conv_track2Chi2RD = 0;
std::vector <double>* conv_minDistanceOfApproachRD = 0;
std::vector <double>* conv_etaRD = 0;
std::vector <double>* conv_ptRD = 0;

// chi
std::vector <double>* chi_etaRD = 0;
std::vector <double>* chi_ptRD = 0;
std::vector <int>* chi_daughterJpsi_positionRD = 0; //stores position of daughter Jpsi in dimuon collection
std::vector <int>* chi_daughterConv_positionRD = 0; //stores position of daughter photon (conversion)
std::vector <double>* chi_dzPhotToDimuonVtxRD = 0; //z distance of photon to dimuon vertex when dxy is minimal
std::vector <double>* chi_dxyPhotToDimuonVtxRD = 0; //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions





bool MuonAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 1.1 && pt < 3.3) return false;
	if (fabs(eta) >= 1.1 && fabs(eta) < 2.1 && pt < 5.5 - 2 * fabs(eta)) return false;
	if (fabs(eta) >= 2.1 && pt < 1.3) return false;
	return true;
}

bool PhotAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.5) return false;
	if (pt < 0.2) return false;
	return true;
}

bool DimuonAcceptance(double eta, double pt)
{
	//if (fabs(eta) > 2.4) return false;  //2.4
	if (pt < 6.0) return false;
	return true;
}

bool MuonSelectionPass(int muonPos)  //uses variables loaded in main function
{
	if (muonIsSoft->at(muonPos) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(muonPos) != 1) return false;
	return true;
}

bool MuonSelectionPassMC(int muonPos)  //uses variables loaded in main function
{
	int matchPosition = gen_muon_matchPosition->at(muonPos);
	if (matchPosition < -0.5) return false; //not matched
	if (muonIsSoft->at(matchPosition) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(matchPosition) != 1) return false;
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



bool PhotSelectionNMinusOne(int photPos, int selToSkip)  //uses variables loaded in main function
{
	if (photPos < -0.5) return false; //not matched
	if (selToSkip != 0 && convQuality_isHighPurity->at(photPos) != 1) return false;
	if (selToSkip != 1 && convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	if (selToSkip != 2 && conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	if (selToSkip != 3 && conv_sigmaTkVtx1->at(photPos) > 10) return false;
	if (selToSkip != 4 && conv_sigmaTkVtx2->at(photPos) > 10) return false;
	if (selToSkip != 5 && conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (selToSkip != 6 && conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (selToSkip != 7 && conv_vertexChi2Prob->at(photPos) <= 0.01) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	if (selToSkip != 8 && fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	if (selToSkip != 9 && conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	if (selToSkip != 10 && conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	if (selToSkip != 11 && conv_track1Chi2->at(photPos) >= 10) return false;
	if (selToSkip != 12 && conv_track2Chi2->at(photPos) >= 10) return false;
	if (selToSkip != 13 && conv_minDistanceOfApproach->at(photPos) <= -0.25) return false;
	if (selToSkip != 14 && conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;
	return true;
}

bool PhotSelectionSingle(int photPos, int selToDo)  //uses variables loaded in main function
{
	if (photPos < -0.5) return false; //not matched
	if (selToDo == 1 && convQuality_isHighPurity->at(photPos) != 1) return false;
	if (selToDo == 2 && convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	if (selToDo == 3 && conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	if (selToDo == 4 && conv_sigmaTkVtx1->at(photPos) > 10) return false;
	if (selToDo == 5 && conv_sigmaTkVtx2->at(photPos) > 10) return false;
	if (selToDo == 6 && conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (selToDo == 7 && conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (selToDo == 8 && conv_vertexChi2Prob->at(photPos) <= 0.01) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	if (selToDo == 9 && fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	if (selToDo == 10 && conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	if (selToDo == 11 && conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	if (selToDo == 12 && conv_track1Chi2->at(photPos) >= 10) return false;
	if (selToDo == 13 && conv_track2Chi2->at(photPos) >= 10) return false;
	if (selToDo == 14 && conv_minDistanceOfApproach->at(photPos) <= -0.25) return false;
	if (selToDo == 15 && conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;
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

bool DimuonSelectionPass(int dimuonPos)  //uses variables loaded in main function
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	//if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	return true;
}

bool DimuonSelectionPassMC(int dimuonPos)  //uses variables loaded in main function
{
	int matchPosition = gen_Jpsi_matchPosition->at(dimuonPos);
	if (matchPosition < -0.5) return false; //not matched
	if (dimuon_charge->at(matchPosition) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	//if (dimuon_vtxProb->at(matchPosition) < 0.01) return false;
	return true;
}


int VariablesDistributions(const char* fileInMC = "Chi_c_pPb8TeV-MC6v4.root", const char* fileInRD = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-bothDirRW2.root", const char* fileOut = "Chi_c_distributionsMC6Tests.root")
{
	gStyle->SetOptStat(1111);
	gROOT->ProcessLine("#include <vector>");

	TH1D* h_muon_ptRel1 = new TH1D("h_muon_ptRel1", "deltaPtRel from matching", 100, 0.0, 0.5);
	TH1D* h_muon_ptRel2 = new TH1D("h_muon_ptRel2", "deltaPtRel by hand", 100, 0.0, 0.5);

	TH1I* h_muonMatch_isGlobal = new TH1I("h_muonMatch_isGlobal", "h_muonMatch_isGlobal", 2, 0, 2);
	TH1I* h_muonMatch_isTracker = new TH1I("h_muonMatch_isTracker", "h_muonMatch_isTracker", 2, 0, 2);
	TH1I* h_muonMatch_isSoft = new TH1I("h_muonMatch_isSoft", "h_muonMatch_isSoft", 2, 0, 2);
	TH1I* h_muonMatch_trackerLayers = new TH1I("h_muonMatch_trackerLayers", "Tracker layers with measurement", 20, 0, 20);

	TH1I* h_gen_muon_matchPosition = new TH1I("h_gen_muon_matchPosition", "h_gen_muon_matchPosition", 10, 0, 10);
	TH1I* h_gen_muon_nMatches = new TH1I("h_gen_muon_nMatches", "h_gen_muon_nMatches", 10, 0, 10);
	TH1D* h_gen_muon_rDelta = new TH1D("h_gen_muon_rDelta", "h_gen_muon_rDelta", 100, 0, 0.5);


	// Jpsi
	TH1I* h_gen_Jpsi_photon_n = new TH1I("h_gen_Jpsi_photon_n", "h_gen_Jpsi_photon_n", 10, 0, 10);
	TH1D* h_gen_Jpsi_photon_pt = new TH1D("h_gen_Jpsi_photon_pt", "h_gen_Jpsi_photon_pt", 100, 0, 2);
	TH1D* h_gen_Jpsi_mass_GEN = new TH1D("h_gen_Jpsi_mass_GEN", "h_gen_Jpsi_mass_GEN", 300, 3, 3.15);
	TH1D* h_gen_Jpsi_mass_GENMuonOnly = new TH1D("h_gen_Jpsi_mass_GENMuonOnly", "h_gen_Jpsi_mass_GENMuonOnly", 300, 3, 3.15);
	TH1D* h_gen_Jpsi_mass_GENMuonOnlyWide = new TH1D("h_gen_Jpsi_mass_GENMuonOnlyWide", "h_gen_Jpsi_mass_GENMuonOnlyWide", 300, 1, 3.3);
	TH1D* h_gen_Jpsi_mass_GENMuonPhoton = new TH1D("h_gen_Jpsi_mass_GENMuonPhoton", "h_gen_Jpsi_mass_GENMuonPhoton", 300, 3, 3.15);
	TH1D* h_gen_Jpsi_mass_GENMuonPhotonWide = new TH1D("h_gen_Jpsi_mass_GENMuonPhotonWide", "h_gen_Jpsi_mass_GENMuonPhotonWide", 300, 1, 5);

	TH1D* h_dimuon_vtxProb = new TH1D("h_dimuon_vtxProb", "h_dimuon_vtxProb", 200, 0, 1);

	//conversions
	//const char* hname = "h_convQuality_isHighPurity";
	//TH1D* h_convQuality_isHighPurity = new TH1D(hname, hname, 2, 0, 2);
	TH1D* h_convQuality_isHighPurity = new TH1D("h_convQuality_isHighPurity", "h_convQuality_isHighPurity", 2, 0, 2);
	TH1D* h_convQuality_isGeneralTracksOnly = new TH1D("h_convQuality_isGeneralTracksOnly", "h_convQuality_isGeneralTracksOnly", 2, 0, 2);
	TH1D* h_conv_vertexPositionRho = new TH1D("h_conv_vertexPositionRho", "h_conv_vertexPositionRho", 200, 0, 100);
	TH1D* h_conv_sigmaTkVtx1 = new TH1D("h_conv_sigmaTkVtx1", "h_conv_sigmaTkVtx1", 200, 0, 100);
	TH1D* h_conv_sigmaTkVtx2 = new TH1D("h_conv_sigmaTkVtx2", "h_conv_sigmaTkVtx2", 200, 0, 100);
	TH1D* h_conv_tkVtxCompatibilityOK = new TH1D("h_conv_tkVtxCompatibilityOK", "h_conv_tkVtxCompatibilityOK", 2, 0, 2);

	TH1D* h_conv_compatibleInnerHitsOK = new TH1D("h_conv_compatibleInnerHitsOK", "h_conv_compatibleInnerHitsOK", 3, -1, 2); //-1: less than 2 tracks, 0: not compatible, 1: yes
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

	TH1D* h_conv_tk1NumOfDOF = new TH1D("h_conv_tk1NumOfDOF", "h_conv_tk1NumOfDOF", 70, 0, 70);
	TH1D* h_conv_tk2NumOfDOF = new TH1D("h_conv_tk2NumOfDOF", "h_conv_tk2NumOfDOF", 70, 0, 70);
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


	//chic
	TH1I* h_chic_Refit_Flag = new TH1I("h_chic_Refit_Flag", "h_chic_Refit_Flag", 13, 0, 13);
	TH1D* h_chic_Refit_vtxProb = new TH1D("h_chic_Refit_vtxProb", "h_chic_Refit_vtxProb", 250, -1.5, 1);
	TH1D* h_chic_Refit_vtxProbZoom = new TH1D("h_chic_Refit_vtxProbZoom", "h_chic_Refit_vtxProbZoom", 100, 0, 1);
	TH1D* h_chic_Refit_vtxProbAfterSel = new TH1D("h_chic_Refit_vtxProbAfterSel", "h_chic_Refit_vtxProbAfterSel", 250, -1.5, 1);
	TH1D* h_chic_Refit_vtxProbAfterSelZoom = new TH1D("h_chic_Refit_vtxProbAfterSelZoom", "h_chic_Refit_vtxProbAfterSelZoom", 100, 0, 1);
	TH2D* h_chic_Refit_vtxProb2D = new TH2D("h_chic_Refit_vtxProb2D", "h_chic_Refit_vtxProb2D;Jpsi_vtxProb;Chic_vtxProb", 100, 0, 1, 100, 0, 1);
	TH1D* h_chic_Refit_ctau = new TH1D("h_chic_Refit_ctau", "h_chic_Refit_ctau", 200, -0.2, 0.2);
	TH1D* h_chic_Refit_ctau3D = new TH1D("h_chic_Refit_ctau3D", "h_chic_Refit_ctau3D", 200, -0.2, 0.2);
	TH1D* h_chic_Diff = new TH1D("h_chic_Diff", "h_chic_Diff", 200, -0.001, 0.001);

	// Conversion study
	TH1D* hConversionCutsInAcc[16];
	const char* histNameInAcc[16] = { "h_convQuality_isHighPurity", "h_convQuality_isGeneralTracksOnly", "h_conv_vertexPositionRho", "h_conv_sigmaTkVtx1", "h_conv_sigmaTkVtx2", "h_conv_tkVtxCompatibilityOK", "h_conv_compatibleInnerHitsOK", "h_conv_vertexChi2Prob",   };

	

	TFile* fMC = new TFile(fileInMC, "READ");


	TTree* event_tree = (TTree*)fMC->Get("ChiRootuple/event_tree");
	if (!event_tree) { 
		cout << "Problem with event Tree";
		return 1;
	}

	event_tree->SetBranchAddress("Trig_Event_HLTDoubleMuOpen", &Trig_Event_HLTDoubleMuOpen);
	event_tree->SetBranchAddress("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
	event_tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal);
	event_tree->SetBranchAddress("muonIsTracker", &muonIsTracker);
	event_tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
	event_tree->SetBranchAddress("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);

	event_tree->SetBranchAddress("muon_eta", &muon_eta);
	event_tree->SetBranchAddress("muon_pt", &muon_pt);
	
	//vertex

	event_tree->SetBranchAddress("pvtx_z", &pvtx_z);
	event_tree->SetBranchAddress("pvtx_zError", &pvtx_zError);
	event_tree->SetBranchAddress("pvtx_x", &pvtx_x);
	event_tree->SetBranchAddress("pvtx_y", &pvtx_y);
	event_tree->SetBranchAddress("pvtx_nTracks", &pvtx_nTracks);
	event_tree->SetBranchAddress("pvtx_isFake", &pvtx_isFake);

	event_tree->SetBranchAddress("dimuon_p4", &dimuon_p4);
	event_tree->SetBranchAddress("dimuon_eta", &dimuon_eta);
	event_tree->SetBranchAddress("dimuon_pt", &dimuon_pt);
	event_tree->SetBranchAddress("dimuon_charge", &dimuon_charge);
	event_tree->SetBranchAddress("dimuon_pvtx_index", &dimuon_pvtx_index);
	event_tree->SetBranchAddress("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);
	event_tree->SetBranchAddress("dimuon_vtxProb", &dimuon_vtxProb);
	event_tree->SetBranchAddress("dimuon_muon1_position", &dimuon_muon1_position);
	event_tree->SetBranchAddress("dimuon_muon2_position", &dimuon_muon2_position);
	event_tree->SetBranchAddress("dimuon_ctpv", &dimuon_ctpv);
	event_tree->SetBranchAddress("dimuon_ctpvError", &dimuon_ctpvError);
	
	// MC

	event_tree->SetBranchAddress("gen_pdgId", &gen_pdgId);
	event_tree->SetBranchAddress("gen_chic_pt", &gen_chic_pt);
	event_tree->SetBranchAddress("gen_chic_eta", &gen_chic_eta);
	event_tree->SetBranchAddress("gen_chic_matchPosition", &gen_chic_matchPosition);
	event_tree->SetBranchAddress("gen_chic_nMatches", &gen_chic_nMatches);
	event_tree->SetBranchAddress("gen_chic_rDelta", &gen_chic_rDelta);
	event_tree->SetBranchAddress("gen_chic_ptDeltaRel", &gen_chic_ptDeltaRel);
	event_tree->SetBranchAddress("gen_chic_p4", &gen_chic_p4);

	event_tree->SetBranchAddress("gen_Jpsi_pt", &gen_Jpsi_pt);
	event_tree->SetBranchAddress("gen_Jpsi_eta", &gen_Jpsi_eta);
	event_tree->SetBranchAddress("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
	event_tree->SetBranchAddress("gen_Jpsi_nMatches", &gen_Jpsi_nMatches);
	event_tree->SetBranchAddress("gen_Jpsi_rDelta", &gen_Jpsi_rDelta);
	event_tree->SetBranchAddress("gen_Jpsi_ptDeltaRel", &gen_Jpsi_ptDeltaRel);
	event_tree->SetBranchAddress("gen_Jpsi_p4", &gen_Jpsi_p4);

	event_tree->SetBranchAddress("gen_Jpsi_photon_n", &gen_Jpsi_photon_n);
	event_tree->SetBranchAddress("gen_Jpsi_photon_pt", &gen_Jpsi_photon_pt);
	event_tree->SetBranchAddress("gen_Jpsi_photon_p4", &gen_Jpsi_photon_p4);

	event_tree->SetBranchAddress("gen_muon_nMatches", &gen_muon_nMatches);
	event_tree->SetBranchAddress("gen_muon_matchPosition", &gen_muon_matchPosition);
	event_tree->SetBranchAddress("gen_muon_eta", &gen_muon_eta);
	event_tree->SetBranchAddress("gen_muon_pt", &gen_muon_pt);
	event_tree->SetBranchAddress("gen_muon_rDelta", &gen_muon_rDelta);
	event_tree->SetBranchAddress("gen_muon_ptDeltaRel", &gen_muon_ptDeltaRel);
	event_tree->SetBranchAddress("gen_muon_p4", &gen_muon_p4);

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
	event_tree->SetBranchAddress("chi_kinematicRefitFlag", &chi_kinematicRefitFlag);
	event_tree->SetBranchAddress("chi_refit_origChicPosition", &chi_refit_origChicPosition);
	event_tree->SetBranchAddress("chi_refit_vprob", &chi_refit_vprob);
	event_tree->SetBranchAddress("chi_refit_ctauPV", &chi_refit_ctauPV);
	event_tree->SetBranchAddress("chi_refit_ctauErrPV", &chi_refit_ctauErrPV);
	event_tree->SetBranchAddress("chi_refit_ctauPV3D", &chi_refit_ctauPV3D);
	event_tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_x", &chi_refit_pvtxFromPVwithMuons_x);
	event_tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_y", &chi_refit_pvtxFromPVwithMuons_y);
	event_tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_z", &chi_refit_pvtxFromPVwithMuons_z);



	///////////////////////////////////
	////////  S  T  A  R  T  //////////
	////////////////////////////////////


	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: "<<nentries << endl;
	//if (nentries > 5000) { nentries = 5000; }
	int counter = 0;

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




		// Look at J/psi
		LVJpsi = (TLorentzVector*)gen_Jpsi_p4->At(0);
		LVmuon1 = (TLorentzVector*)gen_muon_p4->At(0);
		LVmuon2 = (TLorentzVector*)gen_muon_p4->At(1);
		h_gen_Jpsi_photon_n->Fill(gen_Jpsi_photon_n->at(0));
		for (int i_jpsi_phot = 0; i_jpsi_phot < gen_Jpsi_photon_n->at(0); i_jpsi_phot++)
		{
			h_gen_Jpsi_photon_pt->Fill(gen_Jpsi_photon_pt->at(0));
		}
		h_gen_Jpsi_mass_GEN->Fill(LVJpsi->M());
		double mass_Jpsi_muOnly = (*LVmuon1 + *LVmuon2).M();
		h_gen_Jpsi_mass_GENMuonOnly->Fill(mass_Jpsi_muOnly);
		h_gen_Jpsi_mass_GENMuonOnlyWide->Fill(mass_Jpsi_muOnly);
		if (mass_Jpsi_muOnly < 2.9) { counter++; };
		TLorentzVector LV_aux = (*LVmuon1 + *LVmuon2);
		for (int i_jpsi_phot2 = 0; i_jpsi_phot2 < gen_Jpsi_photon_n->at(0); i_jpsi_phot2++)
		{
		//if (gen_Jpsi_photon_n->at(0) > 0) {
		//	LVJpsi_phot = (TLorentzVector*)gen_Jpsi_photon_p4->At(0);
		//	cout << LVJpsi_phot->Pt();
		//}
			LV_aux += *((TLorentzVector*)gen_Jpsi_photon_p4->At(0));
		}
		h_gen_Jpsi_mass_GENMuonPhoton->Fill(LV_aux.M());
		h_gen_Jpsi_mass_GENMuonPhotonWide->Fill(LV_aux.M());

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


		//////////////
		// chic 
		/////////////////

		double eta1 = gen_muon_eta->at(0);
		double pt1 = gen_muon_pt->at(0);
		double eta2 = gen_muon_eta->at(1);
		double pt2 = gen_muon_pt->at(1);
		double eta_phot = gen_phot_eta->at(0);
		double pt_phot = gen_phot_pt->at(0);
		double eta_Jpsi = gen_Jpsi_eta->at(0);
		double pt_Jpsi = gen_Jpsi_pt->at(0);
		double eta_chi = gen_chic_eta->at(0);
		double pt_chi = gen_chic_pt->at(0);

		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true && DimuonAcceptance(eta_Jpsi, pt_Jpsi)) //do only if everything in acceptance
		{
			if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && PhotSelectionPassMCLoose(0) == true && DimuonSelectionPassMC(0))
			{
				h_dimuon_vtxProb->Fill(dimuon_vtxProb->at(gen_Jpsi_matchPosition->at(0)));

				int chicMatchPos = gen_chic_matchPosition->at(0);
				if (chicMatchPos > -0.5) {
					h_chic_Refit_Flag->Fill(chi_kinematicRefitFlag->at(chicMatchPos));
					h_chic_Refit_vtxProb->Fill(chi_refit_vprob->at(chicMatchPos));
					if (dimuon_vtxProb->at(gen_Jpsi_matchPosition->at(0)) > 0.01) { h_chic_Refit_vtxProbAfterSel->Fill(chi_refit_vprob->at(chicMatchPos)); }
					if(chi_kinematicRefitFlag->at(chicMatchPos) > 0.5 && chi_kinematicRefitFlag->at(chicMatchPos)<3.5) //refit actually fine
					{
						h_chic_Refit_vtxProbZoom->Fill(chi_refit_vprob->at(chicMatchPos));
						if (dimuon_vtxProb->at(gen_Jpsi_matchPosition->at(0)) > 0.01) { h_chic_Refit_vtxProbAfterSelZoom->Fill(chi_refit_vprob->at(chicMatchPos)); }
						h_chic_Refit_vtxProb2D->Fill(dimuon_vtxProb->at(gen_Jpsi_matchPosition->at(0)), chi_refit_vprob->at(chicMatchPos));
						h_chic_Refit_ctau->Fill(chi_refit_ctauPV->at(chicMatchPos));
						h_chic_Refit_ctau3D->Fill(chi_refit_ctauPV3D->at(chicMatchPos));
					}
					h_chic_Diff->Fill(chi_refit_pvtxFromPVwithMuons_x->at(chicMatchPos) - pvtx_x->at(0));
					//event_tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_x", &chi_refit_pvtxFromPVwithMuons_x);
					//event_tree->SetBranchAddress("chi_refit_ctauPV", &chi_refit_ctauPV);
					//event_tree->SetBranchAddress("chi_refit_ctauErrPV", &chi_refit_ctauErrPV);
					//event_tree->SetBranchAddress("chi_refit_ctauPV3D", &chi_refit_ctauPV3D);
				}
			}


		}







		////Chic part
		//// ACCEPTANCE
		//double eta1 = gen_muon_eta->at(0);
		//double pt1 = gen_muon_pt->at(0);
		//double eta2 = gen_muon_eta->at(1);
		//double pt2 = gen_muon_pt->at(1);
		//h_muonAcceptance2D_den->Fill(eta1, pt1);
		//h_muonAcceptance2D_den->Fill(eta2, pt2);
		//if (MuonAcceptance(eta1, pt1) == true) h_muonAcceptance2D_num->Fill(eta1, pt1);
		//if (MuonAcceptance(eta2, pt2) == true) h_muonAcceptance2D_num->Fill(eta2, pt2);
		//double eta_phot = gen_phot_eta->at(0);
		//double pt_phot = gen_phot_pt->at(0);
		//h_photAcceptance2D_den->Fill(eta_phot, pt_phot);
		//if (PhotAcceptance(eta_phot, pt_phot) == true) h_photAcceptance2D_num->Fill(eta_phot, pt_phot);
		//double eta_Jpsi = gen_Jpsi_eta->at(0);
		//double pt_Jpsi = gen_Jpsi_pt->at(0);
		//h_JpsiAcceptance2D_den->Fill(eta_Jpsi, pt_Jpsi);
		//if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true) h_JpsiAcceptance2D_num->Fill(eta_Jpsi, pt_Jpsi);
		//TLorentzVector* LVJpsi = (TLorentzVector*)gen_Jpsi_p4->At(0);
		//h_JpsiAcceptance2D_y_den->Fill(LVJpsi->Rapidity(), pt_Jpsi);
		//if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) h_JpsiAcceptance2D_y_num->Fill(LVJpsi->Rapidity(), pt_Jpsi);

		//double eta_chi = gen_chic_eta->at(0);
		//double pt_chi = gen_chic_pt->at(0);
		//h_chiAcceptance2D_den->Fill(eta_chi, pt_chi);
		//if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot)==true) h_chiAcceptance2D_num->Fill(eta_chi, pt_chi);
		//TLorentzVector* LVchic = (TLorentzVector*)gen_chic_p4->At(0);
		//h_chiAcceptance2D_y_den->Fill(LVchic->Rapidity(), pt_chi);
		//if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) h_chiAcceptance2D_y_num->Fill(LVchic->Rapidity(), pt_chi);

		//// EFFICIENCY
		//if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) //do only if everything in acceptance
		//{
		//	h_muonEfficiency2D_den->Fill(eta1, pt1);
		//	h_muonEfficiency2D_den->Fill(eta2, pt2);
		//	if (MuonSelectionPassMC(0) == true) h_muonEfficiency2D_num->Fill(eta1, pt1);
		//	if (MuonSelectionPassMC(1) == true) h_muonEfficiency2D_num->Fill(eta2, pt2);
		//	h_photEfficiency2D_den->Fill(eta_phot, pt_phot);
		//	h_photEfficiency1D_den->Fill(pt_phot);
		//	if (PhotSelectionPassMC(0) == true) {
		//		h_photEfficiency2D_num->Fill(eta_phot, pt_phot);
		//		h_photEfficiency1D_num->Fill(pt_phot);
		//		h_photEfficiency1D_relative_num->Fill(pt_phot);

		//	}
		//	if (PhotSelectionPassMCLoose(0) == true) {
		//		h_photEfficiency1D_relative_den->Fill(pt_phot);
		//	}

		//	h_JpsiEfficiency2D_y_den->Fill(LVJpsi->Rapidity(), pt_Jpsi);
		//	h_JpsiEfficiency1D_den->Fill(pt_Jpsi);
		//	if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0)==true) {
		//		h_JpsiEfficiency2D_y_num->Fill(LVJpsi->Rapidity(), pt_Jpsi);
		//		h_JpsiEfficiency1D_num->Fill(pt_Jpsi);
		//	}

		//	h_chiEfficiency2D_y_den->Fill(LVchic->Rapidity(), pt_chi);
		//	h_chiEfficiency1D_den->Fill(pt_chi);
		//	if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0) == true && PhotSelectionPassMC(0) == true) {
		//		h_chiEfficiency2D_y_num->Fill(LVchic->Rapidity(), pt_chi);
		//		h_chiEfficiency1D_num->Fill(pt_chi);
		//	}
		//}





	} //end of event loop

	cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;

	cout << "Counter: " << counter << endl;

	
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

	// Jpsi

	h_gen_Jpsi_photon_n->Write();
	h_gen_Jpsi_photon_pt->Write();
	h_gen_Jpsi_mass_GEN->Write();
	h_gen_Jpsi_mass_GENMuonOnly->Write();
	h_gen_Jpsi_mass_GENMuonOnlyWide->Write();
	h_gen_Jpsi_mass_GENMuonPhoton->Write();
	h_gen_Jpsi_mass_GENMuonPhotonWide->Write();

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

	h_dimuon_vtxProb->Write();
	h_chic_Refit_Flag->Write();
	h_chic_Refit_vtxProb->Write();
	h_chic_Refit_vtxProbAfterSel->Write();
	h_chic_Refit_vtxProbZoom->Write();
	h_chic_Refit_vtxProbAfterSelZoom->Write();
	h_chic_Refit_vtxProb2D->Write();
	h_chic_Refit_ctau->Write();
	h_chic_Refit_ctau3D->Write();
	h_chic_Diff->Write();


	fout->Close();
	fMC->Close();
	cout << "END" << endl;
	return 0;
}
