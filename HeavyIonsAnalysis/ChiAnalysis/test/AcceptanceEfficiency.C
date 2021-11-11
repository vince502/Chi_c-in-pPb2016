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

#include "../../Macros/tdrstyle.C"
#include "../../Macros/CMS_lumi.C"
#include "../../Macros/ChiTreeInit.C"


using namespace std;

const int PythCode_chic0 = 10441; //Pythia codes
const int PythCode_chic1 = 20443;
const int PythCode_chic2 = 445;
double binsChiEffpT[] = {0.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};
//double binsChiEffpT[] = { 0.0, 0.8, 1.6, 2.4, 3.6, 6.0, 8.0, 10.0, 20.0 };
int  nbinsChiEffpT = sizeof(binsChiEffpT) / sizeof(double) - 1;
double binsConvEffpT[] = { 0.0, 0.1, 0.5, 1.0, 1.5, 2.5, 5.0};
int  nbinsConvEffpT = sizeof(binsConvEffpT) / sizeof(double) - 1;
int weird_decay_counter = 0;

double bins_Q_pT[] = { 6, 9, 12, 18, 30 };
int  nbins_Q_pT = sizeof(bins_Q_pT) / sizeof(double) - 1;

double bins_Q_y[] = { -2.4, -1.6, -1.2, 0, 1.2, 1.6, 2.4 };
int  nbins_Q_y = sizeof(bins_Q_y) / sizeof(double) - 1;

//double binsWeightChi_pT[] = { 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 30.0 };
//int  nbinsWeightChi_pT = sizeof(binsWeightChi_pT) / sizeof(double) - 1;
//double binsWeightChi_absy[] = { 0.0, 0.4, 0.8, 1.2, 1.6, 2.1, 2.4 };
//int  nbinsWeightChi_absy = sizeof(binsWeightChi_absy) / sizeof(double) - 1;

double binsWeightChi_pT[] = { 6.0, 9.0, 12.0, 18.0, 30.0 };
int  nbinsWeightChi_pT = sizeof(binsWeightChi_pT) / sizeof(double) - 1;
double binsWeightChi_absy[] = { 0.0, 1.0, 1.8, 2.4 };
int  nbinsWeightChi_absy = sizeof(binsWeightChi_absy) / sizeof(double) - 1;

double binsWeightChi_nTrk[] = { 0, 50, 100, 150, 200, 400 };
int  nbinsWeightChi_nTrk = sizeof(binsWeightChi_nTrk) / sizeof(double) - 1;


//int AcceptanceEfficiency(const char* fileInMC = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-MC8_BothDir.root", const char* fileOut = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
int AcceptanceEfficiency(const char* fileInMC = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC9-bothDir.root", const char* fileOut = "Chi_c_WeightsMC9_bothDir.root")
{
	gStyle->SetOptStat(1111);
	gROOT->ProcessLine("#include <vector>");

	//TH1D* h_muon_ptRel1 = new TH1D("h_muon_ptRel1", "deltaPtRel from matching", 100, 0.0, 0.5);
	//TH1D* h_muon_ptRel2 = new TH1D("h_muon_ptRel2", "deltaPtRel by hand", 100, 0.0, 0.5);

	//TH1I* h_muonMatch_isGlobal = new TH1I("h_muonMatch_isGlobal", "h_muonMatch_isGlobal", 2, 0, 2);
	//TH1I* h_muonMatch_isTracker = new TH1I("h_muonMatch_isTracker", "h_muonMatch_isTracker", 2, 0, 2);
	//TH1I* h_muonMatch_isSoft = new TH1I("h_muonMatch_isSoft", "h_muonMatch_isSoft", 2, 0, 2);
	//TH1I* h_muonMatch_trackerLayers = new TH1I("h_muonMatch_trackerLayers", "Tracker layers with measurement", 20, 0, 20);

	//TH1I* h_gen_muon_matchPosition = new TH1I("h_gen_muon_matchPosition", "h_gen_muon_matchPosition", 10, 0, 10);
	//TH1I* h_gen_muon_nMatches = new TH1I("h_gen_muon_nMatches", "h_gen_muon_nMatches", 10, 0, 10);
	//TH1D* h_gen_muon_rDelta = new TH1D("h_gen_muon_rDelta", "h_gen_muon_rDelta", 100, 0, 0.5);


	////conversions
	//TH1I* h_convQuality_isHighPurity = new TH1I("h_convQuality_isHighPurity", "h_convQuality_isHighPurity", 2, 0, 2);
	//TH1I* h_convQuality_isGeneralTracksOnly = new TH1I("h_convQuality_isGeneralTracksOnly", "h_convQuality_isGeneralTracksOnly", 2, 0, 2);
	//TH1D* h_conv_vertexPositionRho = new TH1D("h_conv_vertexPositionRho", "h_conv_vertexPositionRho", 200, 0, 100);
	//TH1D* h_conv_sigmaTkVtx1 = new TH1D("h_conv_sigmaTkVtx1", "h_conv_sigmaTkVtx1", 200, 0, 100);
	//TH1D* h_conv_sigmaTkVtx2 = new TH1D("h_conv_sigmaTkVtx2", "h_conv_sigmaTkVtx2", 200, 0, 100);
	//TH1I* h_conv_tkVtxCompatibilityOK = new TH1I("h_conv_tkVtxCompatibilityOK", "h_conv_tkVtxCompatibilityOK", 2, 0, 2);

	//TH1I* h_conv_compatibleInnerHitsOK = new TH1I("h_conv_compatibleInnerHitsOK", "h_conv_compatibleInnerHitsOK", 3, -1, 2); //-1: less than 2 tracks, 0: not compatible, 1: yes
	//TH1D* h_conv_vertexChi2Prob = new TH1D("h_conv_vertexChi2Prob", "h_conv_vertexChi2Prob", 200, 0, 1);
	//TH1D* h_conv_zOfPriVtx = new TH1D("h_conv_zOfPriVtx", "h_conv_zOfPriVtx", 240, -60, 60);
	//TH1D* h_conv_zOfPriVtxFromTracks = new TH1D("h_conv_zOfPriVtxFromTracks", "h_conv_zOfPriVtxFromTracks", 240, -60, 60);
	//TH1D* h_conv_dzToClosestPriVtx = new TH1D("h_conv_dzToClosestPriVtx", "h_conv_dzToClosestPriVtx", 240, -60, 60);
	//TH1D* h_conv_dxyPriVtx_Tr1 = new TH1D("h_conv_dxyPriVtx_Tr1", "h_conv_dxyPriVtx_Tr1", 400, -100, 100);
	//TH1D* h_conv_dxyPriVtx_Tr2 = new TH1D("h_conv_dxyPriVtx_Tr2", "h_conv_dxyPriVtx_Tr2", 400, -100, 100);
	//TH1D* h_conv_dxyPriVtxTimesCharge_Tr1 = new TH1D("h_conv_dxyPriVtxTimesCharge_Tr1", "h_conv_dxyPriVtxTimesCharge_Tr1", 400, -100, 100);
	//TH1D* h_conv_dxyPriVtxTimesCharge_Tr2 = new TH1D("h_conv_dxyPriVtxTimesCharge_Tr2", "h_conv_dxyPriVtxTimesCharge_Tr2", 400, -100, 100);
	//TH1D* h_conv_dxyError_Tr1 = new TH1D("h_conv_dxyError_Tr1", "h_conv_dxyError_Tr1", 200, 0, 10);
	//TH1D* h_conv_dxyError_Tr2 = new TH1D("h_conv_dxyError_Tr2", "h_conv_dxyError_Tr2", 200, 0, 10);

	//TH1I* h_conv_tk1NumOfDOF = new TH1I("h_conv_tk1NumOfDOF", "h_conv_tk1NumOfDOF", 70, 0, 70);
	//TH1I* h_conv_tk2NumOfDOF = new TH1I("h_conv_tk2NumOfDOF", "h_conv_tk2NumOfDOF", 70, 0, 70);
	//TH1D* h_conv_track1Chi2 = new TH1D("h_conv_track1Chi2", "h_conv_track1Chi2", 400, 0, 40);
	//TH1D* h_conv_track2Chi2 = new TH1D("h_conv_track2Chi2", "h_conv_track2Chi2", 400, 0, 40);
	//TH1D* h_conv_minDistanceOfApproach = new TH1D("h_conv_minDistanceOfApproach", "h_conv_minDistanceOfApproach", 400, -2, 2);
	//TH1D* h_conv_eta = new TH1D("h_conv_eta", "h_conv_eta", 200, -5, 5);
	//TH1D* h_conv_pt = new TH1D("h_conv_pt", "h_conv_pt", 600, 0, 30);

	//TH1D* h_gen_phot_pt = new TH1D("h_gen_phot_pt", "h_gen_phot_pt", 600, 0, 30);
	//TH1D* h_gen_phot_eta = new TH1D("h_gen_phot_eta", "h_gen_phot_eta", 200, -5, 5);
	//TH1I* h_gen_conv_matchPosition = new TH1I("h_gen_conv_matchPosition", "h_gen_conv_matchPosition", 10, 0, 10);
	//TH1I* h_gen_conv_nMatches = new TH1I("h_gen_conv_nMatches", "h_gen_conv_nMatches", 10, 0, 10);
	//TH1D* h_gen_conv_rDelta = new TH1D("h_gen_conv_rDelta", "h_gen_conv_rDelta", 100, 0, 0.5);
	//TH1D* h_gen_conv_ptDeltaRel = new TH1D("h_gen_conv_ptDeltaRel", "h_gen_conv_ptDeltaRel", 200, -1, 1);



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

	TH1D* h_JpsiEfficiency1D_Q_den = new TH1D("h_JpsiEfficiency1D_Q_den", "h_JpsiEfficiency1D_Q_den; p_{T}", nbins_Q_pT, bins_Q_pT);
	TH1D* h_JpsiEfficiency1D_Q_num = new TH1D("h_JpsiEfficiency1D_Q_num", "h_JpsiEfficiency1D_Q_num; p_{T}", nbins_Q_pT, bins_Q_pT);
	TH1D* h_JpsiEfficiency1D_Q_rat = new TH1D("h_JpsiEfficiency1D_Q_rat", "h_JpsiEfficiency1D_Q_rat; p_{T}", nbins_Q_pT, bins_Q_pT);
	TH1D* h_chiEfficiency1D_Q_den = new TH1D("h_chiEfficiency1D_Q_den", "h_chiEfficiency1D_Q_den; p_{T}", nbins_Q_pT, bins_Q_pT);
	TH1D* h_chiEfficiency1D_Q_num = new TH1D("h_chiEfficiency1D_Q_num", "h_chiEfficiency1D_Q_num; p_{T}", nbins_Q_pT, bins_Q_pT);
	TH1D* h_chiEfficiency1D_Q_rat = new TH1D("h_chiEfficiency1D_Q_rat", "h_chiEfficiency1D_Q_rat; p_{T}", nbins_Q_pT, bins_Q_pT);
	TH1D* h_chiEfficiency1D_Q_ratRel = new TH1D("h_chiEfficiency1D_Q_ratRel", "h_chiEfficiency1D_Q_ratRel; p_{T}", nbins_Q_pT, bins_Q_pT);

	TH1D* h_JpsiEfficiency1D_Q_y_den = new TH1D("h_JpsiEfficiency1D_Q_y_den", "h_JpsiEfficiency1D_Q_y_den; y", nbins_Q_y, bins_Q_y);
	TH1D* h_JpsiEfficiency1D_Q_y_num = new TH1D("h_JpsiEfficiency1D_Q_y_num", "h_JpsiEfficiency1D_Q_y_num; y", nbins_Q_y, bins_Q_y);
	TH1D* h_JpsiEfficiency1D_Q_y_rat = new TH1D("h_JpsiEfficiency1D_Q_y_rat", "h_JpsiEfficiency1D_Q_y_rat; y", nbins_Q_y, bins_Q_y);
	TH1D* h_chiEfficiency1D_Q_y_den = new TH1D("h_chiEfficiency1D_Q_y_den", "h_chiEfficiency1D_Q_y_den; y", nbins_Q_y, bins_Q_y);
	TH1D* h_chiEfficiency1D_Q_y_num = new TH1D("h_chiEfficiency1D_Q_y_num", "h_chiEfficiency1D_Q_y_num; y", nbins_Q_y, bins_Q_y);
	TH1D* h_chiEfficiency1D_Q_y_rat = new TH1D("h_chiEfficiency1D_Q_y_rat", "h_chiEfficiency1D_Q_y_rat; y", nbins_Q_y, bins_Q_y);
	TH1D* h_chiEfficiency1D_Q_y_ratRel = new TH1D("h_chiEfficiency1D_Q_y_ratRel", "h_chiEfficiency1D_Q_y_ratRel; y", nbins_Q_y, bins_Q_y);


	// weights 
	TH2D* h_JpsiAccEff_den = new TH2D("h_JpsiAccEff_den", "h_AccEff_den; y_{lab}; p_{T}", 30, 0.0, 3.0, 100, 0, 30);
	TH2D* h_JpsiAccEff_num = new TH2D("h_JpsiAccEff_num", "h_AccEff_num; y_{lab}; p_{T}", 30, 0.0, 3.0, 100, 0, 30);
	TH2D* h_JpsiAccEff_rat = new TH2D("h_JpsiAccEff_rat", "h_AccEff_rat; y_{lab}; p_{T}", 30, 0.0, 3.0, 100, 0, 30);

	TH2D* h_chiAccEff_den = new TH2D("h_chiAccEff_den", "h_chiAccEff_den; |y_{lab}|; p_{T}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_pT, binsWeightChi_pT);
	TH2D* h_chiAccEff_num = new TH2D("h_chiAccEff_num", "h_chiAccEff_num; |y_{lab}|; p_{T}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_pT, binsWeightChi_pT);
	TH2D* h_chiAccEff_rat = new TH2D("h_chiAccEff_rat", "h_chiAccEff_rat; |y_{lab}|; p_{T}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_pT, binsWeightChi_pT);

	TH2D* h_chiAccEff_den_nTrk = new TH2D("h_chiAccEff_den_nTrk", "h_chiAccEff_den_nTrk; |y_{lab}|; N_{trk}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_nTrk, binsWeightChi_nTrk);
	TH2D* h_chiAccEff_num_nTrk = new TH2D("h_chiAccEff_num_nTrk", "h_chiAccEff_num_nTrk; |y_{lab}|; N_{trk}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_nTrk, binsWeightChi_nTrk);
	TH2D* h_chiAccEff_rat_nTrk = new TH2D("h_chiAccEff_rat_nTrk", "h_chiAccEff_rat_nTrk; |y_{lab}|; N_{trk}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_nTrk, binsWeightChi_nTrk);



	TFile* fMC = new TFile(fileInMC, "READ");


	TTree* event_tree = (TTree*)fMC->Get("ChiRootuple/event_tree");
	if (!event_tree) { 
		cout << "Problem with event Tree";
		return 1;
	}

	LoadChiBranches(event_tree, true);



	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: "<<nentries << endl;
	//if (nentries > 5000) { nentries = 5000; }


	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);

		//cout << "here" << endl;
		if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }

		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			weird_decay_counter++;
			continue;
		}


		////GEN muons // 2 per event
		//for (Long64_t j = 0; j < 2; j++)
		//{	
		//	int matchPosition = gen_muon_matchPosition->at(j);
		//	if (matchPosition < -0.5) { continue; } //negative means not matched
		//	double muonPtRel1 = gen_muon_ptDeltaRel->at(j);
		//	double muonPtRel2 = fabs(muon_pt->at(matchPosition) - gen_muon_pt->at(j)) / gen_muon_pt->at(j);
		//	
		//	h_gen_muon_matchPosition->Fill(matchPosition);
		//	h_gen_muon_nMatches->Fill(gen_muon_nMatches->at(j));
		//	h_gen_muon_rDelta->Fill(gen_muon_rDelta->at(j));

		//	h_muon_ptRel1->Fill(muonPtRel1);
		//	h_muon_ptRel2->Fill(muonPtRel2);
		//	
		//	bool muonTracker = muonIsTracker->at(matchPosition);
		//	bool muonGlobal = muonIsGlobal->at(matchPosition);
		//	h_muonMatch_isGlobal->Fill(muonGlobal);
		//	h_muonMatch_isTracker->Fill(muonTracker);
		//	h_muonMatch_isSoft->Fill(muonIsSoft->at(matchPosition));
		//	


		//	//if (!muonGlobal || !muonTracker) { cout << muonGlobal << "  " << muonTracker << endl; }
		//	if (muonGlobal || muonTracker) {
		//		h_muonMatch_trackerLayers->Fill(muonTrackerLayersWithMeasurement->at(matchPosition));
		//	}
		//	//cout << "gen eta: " << gen_muon_eta->at(j) << "  and muon eta: " << muon_eta->at(matchPosition) <<"   and rDelta: " << gen_muon_rDelta->at(j)<< endl;
		//	//cout << "gen pt: " << gen_muon_pt->at(j) << "  and muon pt: " << muon_pt->at(matchPosition) << "   and ptDelta: " << gen_muon_ptDeltaRel->at(j) << endl<<endl;
		//}
		//// Conversions from chic //one per event
		//int conv_matchPosition = gen_conv_matchPosition->at(0);
		//if (conv_matchPosition > -0.5) { //non-negative means it was matched
		//	h_convQuality_isHighPurity->Fill(convQuality_isHighPurity->at(conv_matchPosition));
		//	h_convQuality_isGeneralTracksOnly->Fill(convQuality_isGeneralTracksOnly->at(conv_matchPosition));
		//	h_conv_vertexPositionRho->Fill(conv_vertexPositionRho->at(conv_matchPosition));
		//	h_conv_sigmaTkVtx1->Fill(conv_sigmaTkVtx1->at(conv_matchPosition));
		//	h_conv_sigmaTkVtx2->Fill(conv_sigmaTkVtx2->at(conv_matchPosition));
		//	h_conv_tkVtxCompatibilityOK->Fill(conv_tkVtxCompatibilityOK->at(conv_matchPosition));

		//	h_conv_compatibleInnerHitsOK->Fill(conv_compatibleInnerHitsOK->at(conv_matchPosition)); //-1: less than 2 tracks, 0: not compatible, 1: yes
		//	h_conv_vertexChi2Prob->Fill(conv_vertexChi2Prob->at(conv_matchPosition));
		//	h_conv_zOfPriVtx->Fill(conv_zOfPriVtx->at(conv_matchPosition));
		//	h_conv_zOfPriVtxFromTracks->Fill(conv_zOfPriVtxFromTracks->at(conv_matchPosition));
		//	h_conv_dzToClosestPriVtx->Fill(conv_dzToClosestPriVtx->at(conv_matchPosition));
		//	h_conv_dxyPriVtx_Tr1->Fill(conv_dxyPriVtx_Tr1->at(conv_matchPosition));
		//	h_conv_dxyPriVtx_Tr2->Fill(conv_dxyPriVtx_Tr2->at(conv_matchPosition));
		//	h_conv_dxyPriVtxTimesCharge_Tr1->Fill(conv_dxyPriVtxTimesCharge_Tr1->at(conv_matchPosition));
		//	h_conv_dxyPriVtxTimesCharge_Tr2->Fill(conv_dxyPriVtxTimesCharge_Tr2->at(conv_matchPosition));
		//	h_conv_dxyError_Tr1->Fill(conv_dxyError_Tr1->at(conv_matchPosition));
		//	h_conv_dxyError_Tr2->Fill(conv_dxyError_Tr2->at(conv_matchPosition));

		//	h_conv_tk1NumOfDOF->Fill(conv_tk1NumOfDOF->at(conv_matchPosition));
		//	h_conv_tk2NumOfDOF->Fill(conv_tk2NumOfDOF->at(conv_matchPosition));
		//	h_conv_track1Chi2->Fill(conv_track1Chi2->at(conv_matchPosition));
		//	h_conv_track2Chi2->Fill(conv_track2Chi2->at(conv_matchPosition));
		//	h_conv_minDistanceOfApproach->Fill(conv_minDistanceOfApproach->at(conv_matchPosition));
		//	h_conv_eta->Fill(conv_eta->at(conv_matchPosition));
		//	h_conv_pt->Fill(conv_pt->at(conv_matchPosition));

		//	h_gen_phot_pt->Fill(gen_phot_pt->at(0));
		//	h_gen_phot_eta->Fill(gen_phot_eta->at(0));
		//	h_gen_conv_matchPosition->Fill(gen_conv_matchPosition->at(0));
		//	h_gen_conv_nMatches->Fill(gen_conv_nMatches->at(0));
		//	h_gen_conv_rDelta->Fill(gen_conv_rDelta->at(0));
		//	h_gen_conv_ptDeltaRel->Fill(gen_conv_ptDeltaRel->at(0));

		//}

	
		// ACCEPTANCE
		double eta1 = gen_muon_eta->at(0);
		double pt1 = gen_muon_pt->at(0);
		double eta2 = gen_muon_eta->at(1);
		double pt2 = gen_muon_pt->at(1);
		double eta_phot = gen_phot_eta->at(0);
		double pt_phot = gen_phot_pt->at(0);
		double eta_Jpsi = gen_Jpsi_eta->at(0);
		double pt_Jpsi = gen_Jpsi_pt->at(0);
		TLorentzVector* LVJpsi = (TLorentzVector*)gen_Jpsi_p4->At(0);
		double rap_Jpsi = LVJpsi->Rapidity();

		h_muonAcceptance2D_den->Fill(eta1, pt1);
		h_muonAcceptance2D_den->Fill(eta2, pt2);
		if (MuonAcceptance(eta1, pt1) == true) h_muonAcceptance2D_num->Fill(eta1, pt1);
		if (MuonAcceptance(eta2, pt2) == true) h_muonAcceptance2D_num->Fill(eta2, pt2);

		h_photAcceptance2D_den->Fill(eta_phot, pt_phot);
		if (PhotAcceptance(eta_phot, pt_phot) == true) h_photAcceptance2D_num->Fill(eta_phot, pt_phot);

		h_JpsiAcceptance2D_den->Fill(eta_Jpsi, pt_Jpsi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true) h_JpsiAcceptance2D_num->Fill(eta_Jpsi, pt_Jpsi);
		h_JpsiAcceptance2D_y_den->Fill(rap_Jpsi, pt_Jpsi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true) h_JpsiAcceptance2D_y_num->Fill(rap_Jpsi, pt_Jpsi);

		double eta_chi = gen_chic_eta->at(0);
		double pt_chi = gen_chic_pt->at(0);
		TLorentzVector* LVchic = (TLorentzVector*)gen_chic_p4->At(0);
		double rap_chi = LVchic->Rapidity();

		h_chiAcceptance2D_den->Fill(eta_chi, pt_chi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot)==true) h_chiAcceptance2D_num->Fill(eta_chi, pt_chi);
		
		h_chiAcceptance2D_y_den->Fill(rap_chi, pt_chi);
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) h_chiAcceptance2D_y_num->Fill(rap_chi, pt_chi);

		// EFFICIENCY
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && PhotAcceptance(eta_phot, pt_phot) == true) //do only if J/psi in acceptance
		{
			h_muonEfficiency2D_den->Fill(eta1, pt1);
			h_muonEfficiency2D_den->Fill(eta2, pt2);
			if (MuonSelectionPassMC(0) == true) h_muonEfficiency2D_num->Fill(eta1, pt1);
			if (MuonSelectionPassMC(1) == true) h_muonEfficiency2D_num->Fill(eta2, pt2);

			h_JpsiEfficiency2D_y_den->Fill(rap_Jpsi, pt_Jpsi);
			h_JpsiEfficiency1D_den->Fill(pt_Jpsi);
			h_JpsiEfficiency1D_Q_den->Fill(pt_Jpsi);
			if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0) == true) {
				h_JpsiEfficiency2D_y_num->Fill(rap_Jpsi, pt_Jpsi);
				h_JpsiEfficiency1D_num->Fill(pt_Jpsi);
				h_JpsiEfficiency1D_Q_num->Fill(pt_Jpsi);
			}
			h_chiEfficiency1D_Q_den->Fill(pt_chi);// is used as a quick ratio of the efficiencies, needs to take into account photon

			if (PhotAcceptance(eta_phot, pt_phot) == true) {
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

				h_chiEfficiency2D_y_den->Fill(rap_chi, pt_chi);
				h_chiEfficiency1D_den->Fill(pt_chi);
				
				if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0) == true && PhotSelectionPassMC(0) == true) {
					h_chiEfficiency2D_y_num->Fill(rap_chi, pt_chi);
					h_chiEfficiency1D_num->Fill(pt_chi);
					h_chiEfficiency1D_Q_num->Fill(pt_chi);
				}
			}
		}


		//weights (<acc*eff>, then to be used as 1/w)
		h_JpsiAccEff_den->Fill(abs(rap_Jpsi), pt_Jpsi);
		h_chiAccEff_den->Fill(abs(rap_chi), pt_chi);
		h_chiAccEff_den_nTrk->Fill(abs(rap_chi), ntracks_inEvent);

		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0) == true)
		{
			h_JpsiAccEff_num->Fill(abs(rap_Jpsi), pt_Jpsi);
			if (PhotAcceptance(eta_phot, pt_phot) == true && PhotSelectionPassMC(0) == true && ChiSelectionPassMC(0))
			{
				h_chiAccEff_num->Fill(abs(rap_chi), pt_chi);
				h_chiAccEff_num_nTrk->Fill(abs(rap_chi), ntracks_inEvent);
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
	h_JpsiEfficiency1D_Q_num->Sumw2();
	h_JpsiEfficiency1D_Q_den->Sumw2();
	h_chiEfficiency1D_Q_num->Sumw2();
	h_chiEfficiency1D_Q_den->Sumw2();

	h_photEfficiency1D_rat->Divide(h_photEfficiency1D_num, h_photEfficiency1D_den, 1, 1, "B");
	h_photEfficiency1D_relative_rat->Divide(h_photEfficiency1D_relative_num, h_photEfficiency1D_relative_den, 1, 1, "B");
	h_JpsiEfficiency1D_rat->Divide(h_JpsiEfficiency1D_num, h_JpsiEfficiency1D_den, 1, 1, "B");
	h_chiEfficiency1D_rat->Divide(h_chiEfficiency1D_num, h_chiEfficiency1D_den, 1, 1, "B");
	h_JpsiEfficiency1D_Q_rat->Divide(h_JpsiEfficiency1D_Q_num, h_JpsiEfficiency1D_Q_den, 1, 1, "B");
	h_chiEfficiency1D_Q_rat->Divide(h_chiEfficiency1D_Q_num, h_chiEfficiency1D_Q_den, 1, 1, "B");
	h_chiEfficiency1D_Q_ratRel->Divide(h_chiEfficiency1D_Q_rat, h_JpsiEfficiency1D_Q_rat, 1, 1, "B");
	

	h_JpsiAccEff_den->Sumw2();
	h_JpsiAccEff_num->Sumw2();
	h_chiAccEff_den->Sumw2();
	h_chiAccEff_num->Sumw2();
	h_chiAccEff_den_nTrk->Sumw2();
	h_chiAccEff_num_nTrk->Sumw2();

	h_JpsiAccEff_rat->Divide(h_JpsiAccEff_num, h_JpsiAccEff_den);
	h_chiAccEff_rat->Divide(h_chiAccEff_num, h_chiAccEff_den);
	h_chiAccEff_rat_nTrk->Divide(h_chiAccEff_num_nTrk, h_chiAccEff_den_nTrk);

	TCanvas* can1 = new TCanvas("can1", "plot", 1200, 800);
	h_muonAcceptance2D_num->Draw("colz");
	TF1* f_filter = new TF1("f_filter", "3.3/TMath::CosH(x)", -2.5, 2.5);
	f_filter->SetLineColor(kBlue);
	f_filter->Draw("same");
	TF1* f_acceptance = new TF1("f_acceptance","(fabs(x) > 2.4)*100+((fabs(x) < 0.3)*3.4 + (fabs(x)>=0.3)*(fabs(x)<1.1)*3.3 + (fabs(x)>=1.1)*(fabs(x)<2.1)*(5.5-2.0*fabs(x))+(fabs(x)>=2.1)*1.3)",-2.5,2.5);
	//TF1* f_acceptance = new TF1("f_acceptance","(fabs(x) > 2.4)*100+((fabs(x) < 0.8)*3.3 + (fabs(x)>=0.8)*(fabs(x)<1.5)*(5.81-3.14*fabs(x))+(fabs(x)>=1.5)*(fabs(x)<2.07)*(1.89-0.526*fabs(x))+(fabs(x)>=2.07)*0.8)",-2.5,2.5);
	f_acceptance->Draw("same");
	//can1->SaveAs("muonAcceptance2D_num.png");
	
	TFile* fout = new TFile(fileOut, "RECREATE");

	//h_muon_ptRel1->Write();
	//h_muon_ptRel2->Write();

	//h_gen_muon_matchPosition->Write();
	//h_gen_muon_nMatches->Write();
	//h_gen_muon_rDelta->Write();

	//h_muonMatch_isGlobal->Write();
	//h_muonMatch_isTracker->Write();
	//h_muonMatch_isSoft->Write();
	//h_muonMatch_trackerLayers->Write();



	////conv

	//h_convQuality_isHighPurity->Write();
	//h_convQuality_isGeneralTracksOnly->Write();
	//h_conv_vertexPositionRho->Write();
	//h_conv_sigmaTkVtx1->Write();
	//h_conv_sigmaTkVtx2->Write();
	//h_conv_tkVtxCompatibilityOK->Write();

	//h_conv_compatibleInnerHitsOK->Write(); //-1: less than 2 tracks, 0: not compatible, 1: yes
	//h_conv_vertexChi2Prob->Write();
	//h_conv_zOfPriVtx->Write();
	//h_conv_zOfPriVtxFromTracks->Write();
	//h_conv_dzToClosestPriVtx->Write();
	//h_conv_dxyPriVtx_Tr1->Write();
	//h_conv_dxyPriVtx_Tr2->Write();
	//h_conv_dxyPriVtxTimesCharge_Tr1->Write();
	//h_conv_dxyPriVtxTimesCharge_Tr2->Write();
	//h_conv_dxyError_Tr1->Write();
	//h_conv_dxyError_Tr2->Write();

	//h_conv_tk1NumOfDOF->Write();
	//h_conv_tk2NumOfDOF->Write();
	//h_conv_track1Chi2->Write();
	//h_conv_track2Chi2->Write();
	//h_conv_minDistanceOfApproach->Write();
	//h_conv_eta->Write();
	//h_conv_pt->Write();

	//h_gen_phot_pt->Write();
	//h_gen_phot_eta->Write();
	//h_gen_conv_matchPosition->Write();
	//h_gen_conv_nMatches->Write();
	//h_gen_conv_rDelta->Write();
	//h_gen_conv_ptDeltaRel->Write();

	can1->Write();
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
	h_chiEfficiency1D_Q_ratRel->Write();

	//weights
	h_JpsiAccEff_den->Write();
	h_JpsiAccEff_num->Write();
	h_chiAccEff_den->Write();
	h_chiAccEff_num->Write();
	h_chiAccEff_den_nTrk->Write();
	h_chiAccEff_num_nTrk->Write();

	h_JpsiAccEff_rat->Write();
	h_chiAccEff_rat->Write();
	h_chiAccEff_rat_nTrk->Write();


	fout->Close();
	fMC->Close();
	cout << "END" << endl;
	return 0;
}
