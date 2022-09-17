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
#include "../../Macros/ChiFitterInit.h"


using namespace std;

//const int PythCode_chic0 = 10441; //Pythia codes
//const int PythCode_chic1 = 20443;
//const int PythCode_chic2 = 445;
//double binsChiEffpT[] = {6.5, 8.0, 10.0, 13.0, 17.0, 23.0, 30.0};
double binsChiEffpT[] = { 0.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 17.0, 20.0, 23.0, 26.0, 30.0 };
int  nbinsChiEffpT = sizeof(binsChiEffpT) / sizeof(double) - 1;
//double binsConvEffpT[] = { 0.0, 0.1, 0.5, 1.0, 1.5, 2.5, 5.0};
double binsConvEffpT[] = { 0.0, 0.3, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 7.0, 10.0 };
int  nbinsConvEffpT = sizeof(binsConvEffpT) / sizeof(double) - 1;
double binsConvEffy[] = { -5.0, -4.0, -3.0, -2.6, -2.5, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 3.0, 4.0, 5.0 };
int  nbinsConvEffy = sizeof(binsConvEffy) / sizeof(double) - 1;

int weird_decay_counter = 0;

double bins_Q_pT[] = { 6.5, 9, 12, 18, 30 };
//double bins_Q_pT[] = { 6.5, 9, 12, 16, 22, 30 };
int  nbins_Q_pT = sizeof(bins_Q_pT) / sizeof(double) - 1;

double bins_Q_y[] = { -2.4, -1.6, -1.0, 0, 1.0, 1.6, 2.4 };
int  nbins_Q_y = sizeof(bins_Q_y) / sizeof(double) - 1;

//double binsWeightChi_pT[] = { 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 30.0 };
//int  nbinsWeightChi_pT = sizeof(binsWeightChi_pT) / sizeof(double) - 1;
//double binsWeightChi_absy[] = { 0.0, 0.4, 0.8, 1.2, 1.6, 2.1, 2.4 };
//int  nbinsWeightChi_absy = sizeof(binsWeightChi_absy) / sizeof(double) - 1;

double binsWeightChi_pT[] = { 6.5, 9, 12, 18, 30 };
int  nbinsWeightChi_pT = sizeof(binsWeightChi_pT) / sizeof(double) - 1;
double binsWeightChi_absy[] = { 0.0, 1.0, 1.6, 2.4 };
int  nbinsWeightChi_absy = sizeof(binsWeightChi_absy) / sizeof(double) - 1;

double binsWeightChi_nTrk[] = { 0, 50, 100, 150, 250, 400 };
int  nbinsWeightChi_nTrk = sizeof(binsWeightChi_nTrk) / sizeof(double) - 1;


int AcceptanceEfficiency(const char* fileInMC = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC-Official_v3-bothDir.root", const char* fileOut = "Chi_c_WeightsMC_Official_v3-bothDirTest.root", const char* fileMCWeight = "MCWeight_v2.root")
//int AcceptanceEfficiency(const char* fileInMC = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC-Official_v2-Pbp.root", const char* fileOut = "Chi_c_WeightsMC_Official-PbpT.root", const char* fileMCWeight = "MCWeight_v1.root")
{
	gStyle->SetOptStat(1111);
	gROOT->ProcessLine("#include <vector>");

	// load the MC weights
	TH1D* h_weightPVtrk;
	TFile* fMCWeight = new TFile(fileMCWeight, "READ");
	h_weightPVtrk= (TH1D*)fMCWeight->Get("h_weightPVtrk");
	h_weightPVtrk->SetDirectory(0);

	fMCWeight->Close();

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

	// photon disctribution, for the j/psi bins
	TH3D* h_photDistribution3D_JpsipT = new TH3D("h_photDistribution3D_JpsipT", "h_photDistribution3D_JpsipT; p_{T}; #eta; p_{T}", nbins_pT, bins_pT, nbinsConvEffy, binsConvEffy, nbinsConvEffpT, binsConvEffpT);
	TH3D* h_photDistribution3D_Jpsiy = new TH3D("h_photDistribution3D_Jpsiy", "h_photDistribution3D_Jpsiy; y; #eta; p_{T}", nbins_y, bins_y, nbinsConvEffy, binsConvEffy, nbinsConvEffpT, binsConvEffpT);


	// EFFICIENCY
	TH2D* h_muonEfficiency2D_den = new TH2D("h_muonEfficiency2D_den", "h_muonEfficiency2D_den; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_muonEfficiency2D_num = new TH2D("h_muonEfficiency2D_num", "h_muonEfficiency2D_num; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_muonEfficiency2D_rat = new TH2D("h_muonEfficiency2D_rat", "h_muonEfficiency2D_rat; #eta; p_{T}", 60, -3.0, 3.0, 100, 0, 10);
	TH2D* h_photEfficiency2D_den = new TH2D("h_photEfficiency2D_den", "h_photEfficiency2D_den; #eta; p_{T}", 60, -3.0, 3.0, nbinsConvEffpT, binsConvEffpT);
	TH2D* h_photEfficiency2D_num = new TH2D("h_photEfficiency2D_num", "h_photEfficiency2D_num; #eta; p_{T}", 60, -3.0, 3.0, nbinsConvEffpT, binsConvEffpT);
	TH2D* h_photEfficiency2D_rat = new TH2D("h_photEfficiency2D_rat", "h_photEfficiency2D_rat; #eta; p_{T}", 60, -3.0, 3.0, nbinsConvEffpT, binsConvEffpT);
	TH2D* h_JpsiEfficiency2D_y_den = new TH2D("h_JpsiEfficiency2D_y_den", "h_JpsiEfficiency2D_y_den; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_JpsiEfficiency2D_y_num = new TH2D("h_JpsiEfficiency2D_y_num", "h_JpsiEfficiency2D_y_num; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_JpsiEfficiency2D_y_rat = new TH2D("h_JpsiEfficiency2D_y_rat", "h_JpsiEfficiency2D_y_rat; y_{lab}; p_{T}", 60, -3.0, 3.0, 100, 0, 20);
	TH2D* h_chiEfficiency2D_y_den = new TH2D("h_chiEfficiency2D_y_den", "h_chiEfficiency2D_y_den; y_{lab}; p_{T}", 60, -3.0, 3.0, nbinsChiEffpT, binsChiEffpT);
	TH2D* h_chiEfficiency2D_y_num = new TH2D("h_chiEfficiency2D_y_num", "h_chiEfficiency2D_y_num; y_{lab}; p_{T}", 60, -3.0, 3.0, nbinsChiEffpT, binsChiEffpT);
	TH2D* h_chiEfficiency2D_y_rat = new TH2D("h_chiEfficiency2D_y_rat", "h_chiEfficiency2D_y_rat; y_{lab}; p_{T}", 60, -3.0, 3.0, nbinsChiEffpT, binsChiEffpT);

	TH1D* h_photEfficiency1D_den = new TH1D("h_photEfficiency1D_den", "h_photEfficiency1D_den; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_photEfficiency1D_num = new TH1D("h_photEfficiency1D_num", "h_photEfficiency1D_num; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_photEfficiency1D_rat = new TH1D("h_photEfficiency1D_rat", "h_photEfficiency1D_rat; p_{T}", nbinsConvEffpT, binsConvEffpT);
	TH1D* h_JpsiEfficiency1D_den = new TH1D("h_JpsiEfficiency1D_den", "h_JpsiEfficiency1D_den; p_{T}", 100, 0, 20);
	TH1D* h_JpsiEfficiency1D_num = new TH1D("h_JpsiEfficiency1D_num", "h_JpsiEfficiency1D_num; p_{T}", 100, 0, 20);
	TH1D* h_JpsiEfficiency1D_rat = new TH1D("h_JpsiEfficiency1D_rat", "h_JpsiEfficiency1D_rat; p_{T}", 100, 0, 20);
	TH1D* h_chiEfficiency1D_den = new TH1D("h_chiEfficiency1D_den", "h_chiEfficiency1D_den; p_{T}", nbinsChiEffpT, binsChiEffpT);
	TH1D* h_chiEfficiency1D_num = new TH1D("h_chiEfficiency1D_num", "h_chiEfficiency1D_num; p_{T}", nbinsChiEffpT, binsChiEffpT);
	TH1D* h_chiEfficiency1D_rat = new TH1D("h_chiEfficiency1D_rat", "h_chiEfficiency1D_rat; p_{T}", nbinsChiEffpT, binsChiEffpT);

	// y

	TH1D* h_photEfficiency1D_y_den = new TH1D("h_photEfficiency1D_y_den", "h_photEfficiency1D_y_den; y", nbins_y, bins_y);
	TH1D* h_photEfficiency1D_y_denJpsiVar = new TH1D("h_photEfficiency1D_y_denJpsiVar", "h_photEfficiency1D_y_denJpsiVar; y", nbins_y, bins_y);
	TH1D* h_photEfficiency1D_y_num = new TH1D("h_photEfficiency1D_y_num", "h_photEfficiency1D_y_num; y", nbins_y, bins_y);
	TH1D* h_photEfficiency1D_y_rat = new TH1D("h_photEfficiency1D_y_rat", "h_photEfficiency1D_y_rat; y", nbins_y, bins_y);
	TH1D* h_JpsiEfficiency1D_y_den = new TH1D("h_JpsiEfficiency1D_y_den", "h_JpsiEfficiency1D_y_den; y", nbins_y, bins_y);
	TH1D* h_JpsiEfficiency1D_y_num = new TH1D("h_JpsiEfficiency1D_y_num", "h_JpsiEfficiency1D_y_num; y", nbins_y, bins_y);
	TH1D* h_JpsiEfficiency1D_y_rat = new TH1D("h_JpsiEfficiency1D_y_rat", "h_JpsiEfficiency1D_y_rat; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_y_den = new TH1D("h_chiEfficiency1D_y_den", "h_chiEfficiency1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_y_num = new TH1D("h_chiEfficiency1D_y_num", "h_chiEfficiency1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_y_rat = new TH1D("h_chiEfficiency1D_y_rat", "h_chiEfficiency1D_y_rat; y", nbins_y, bins_y);


	// ntracks


	TH1D* h_photEfficiency1D_nTrack_den = new TH1D("h_photEfficiency1D_nTrack_den", "h_photEfficiency1D_nTrack_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_photEfficiency1D_nTrack_num = new TH1D("h_photEfficiency1D_nTrack_num", "h_photEfficiency1D_nTrack_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_photEfficiency1D_nTrack_rat = new TH1D("h_photEfficiency1D_nTrack_rat", "h_photEfficiency1D_nTrack_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_JpsiEfficiency1D_nTrack_den = new TH1D("h_JpsiEfficiency1D_nTrack_den", "h_JpsiEfficiency1D_nTrack_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_JpsiEfficiency1D_nTrack_num = new TH1D("h_JpsiEfficiency1D_nTrack_num", "h_JpsiEfficiency1D_nTrack_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_JpsiEfficiency1D_nTrack_rat = new TH1D("h_JpsiEfficiency1D_nTrack_rat", "h_JpsiEfficiency1D_nTrack_rat; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_nTrack_den = new TH1D("h_chiEfficiency1D_nTrack_den", "h_chiEfficiency1D_nTrack_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_nTrack_num = new TH1D("h_chiEfficiency1D_nTrack_num", "h_chiEfficiency1D_nTrack_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_nTrack_rat = new TH1D("h_chiEfficiency1D_nTrack_rat", "h_chiEfficiency1D_nTrack_rat; nTrk", nbins_nTrk, bins_nTrk);





	// ratio of conversion eff //in bins of chic y  - crosscheck only
	TH1D* h_photEfficiency_y_relative_den = new TH1D("h_photEfficiency_y_relative_den", "h_photEfficiency_y_relative_den; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_num_loose = new TH1D("h_photEfficiency_y_relative_num_loose", "h_photEfficiency_y_relative_num_loose; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_rat_loose = new TH1D("h_photEfficiency_y_relative_rat_loose", "h_photEfficiency_y_relative_rat_loose; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_num_medium = new TH1D("h_photEfficiency_y_relative_num_medium", "h_photEfficiency_y_relative_num_medium; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_rat_medium = new TH1D("h_photEfficiency_y_relative_rat_medium", "h_photEfficiency_y_relative_rat_medium; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_num_tight = new TH1D("h_photEfficiency_y_relative_num_tight", "h_photEfficiency_y_relative_num_tight; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_rat_tight = new TH1D("h_photEfficiency_y_relative_rat_tight", "h_photEfficiency_y_relative_rat_tight; y", nbins_y, bins_y);

	TH1D* h_photEfficiency_y_relative_num_noCut = new TH1D("h_photEfficiency_y_relative_num_noCut", "h_photEfficiency_y_relative_num_noCut; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_rat_noCut = new TH1D("h_photEfficiency_y_relative_rat_noCut", "h_photEfficiency_y_relative_rat_noCut; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_num_fullCut = new TH1D("h_photEfficiency_y_relative_num_fullCut", "h_photEfficiency_y_relative_num_fullCut; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_rat_fullCut = new TH1D("h_photEfficiency_y_relative_rat_fullCut", "h_photEfficiency_y_relative_rat_fullCut; y", nbins_y, bins_y);

	TH1D* h_photEfficiency_y_relative_num_noCutJpsiVar = new TH1D("h_photEfficiency_y_relative_num_noCutJpsiVar", "h_photEfficiency_y_relative_num_noCutJpsiVar; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_rat_noCutJpsiVar = new TH1D("h_photEfficiency_y_relative_rat_noCutJpsiVar", "h_photEfficiency_y_relative_rat_noCutJpsiVar; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_num_fullCutJpsiVar = new TH1D("h_photEfficiency_y_relative_num_fullCutJpsiVar", "h_photEfficiency_y_relative_num_fullCutJpsiVar; y", nbins_y, bins_y);
	TH1D* h_photEfficiency_y_relative_rat_fullCutJpsiVar = new TH1D("h_photEfficiency_y_relative_rat_fullCutJpsiVar", "h_photEfficiency_y_relative_rat_fullCutJpsiVar; y", nbins_y, bins_y);




	// weights 
	TH2D* h_JpsiAccEff_den = new TH2D("h_JpsiAccEff_den", "h_AccEff_den; y_{lab}; p_{T}", 30, 0.0, 3.0, 100, 0, 30);
	TH2D* h_JpsiAccEff_num = new TH2D("h_JpsiAccEff_num", "h_AccEff_num; y_{lab}; p_{T}", 30, 0.0, 3.0, 100, 0, 30);
	TH2D* h_JpsiAccEff_rat = new TH2D("h_JpsiAccEff_rat", "h_AccEff_rat; y_{lab}; p_{T}", 30, 0.0, 3.0, 100, 0, 30);

	TH2D* h_JpsiAccEff_den_nTrk = new TH2D("h_JpsiAccEff_den_nTrk", "h_JpsiAccEff_den_nTrk; |y_{lab}|; N_{trk}", 30, 0.0, 3.0, nbinsWeightChi_nTrk, binsWeightChi_nTrk);
	TH2D* h_JpsiAccEff_num_nTrk = new TH2D("h_JpsiAccEff_num_nTrk", "h_JpsiAccEff_num_nTrk; |y_{lab}|; N_{trk}", 30, 0.0, 3.0, nbinsWeightChi_nTrk, binsWeightChi_nTrk);
	TH2D* h_JpsiAccEff_rat_nTrk = new TH2D("h_JpsiAccEff_rat_nTrk", "h_JpsiAccEff_rat_nTrk; |y_{lab}|; N_{trk}", 30, 0.0, 3.0, nbinsWeightChi_nTrk, binsWeightChi_nTrk);

	TH2D* h_chiAccEff_den = new TH2D("h_chiAccEff_den", "h_chiAccEff_den; |y_{lab}|; p_{T}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_pT, binsWeightChi_pT);
	TH2D* h_chiAccEff_num = new TH2D("h_chiAccEff_num", "h_chiAccEff_num; |y_{lab}|; p_{T}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_pT, binsWeightChi_pT);
	TH2D* h_chiAccEff_rat = new TH2D("h_chiAccEff_rat", "h_chiAccEff_rat; |y_{lab}|; p_{T}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_pT, binsWeightChi_pT);

	TH2D* h_chiAccEff_den_nTrk = new TH2D("h_chiAccEff_den_nTrk", "h_chiAccEff_den_nTrk; |y_{lab}|; N_{trk}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_nTrk, binsWeightChi_nTrk);
	TH2D* h_chiAccEff_num_nTrk = new TH2D("h_chiAccEff_num_nTrk", "h_chiAccEff_num_nTrk; |y_{lab}|; N_{trk}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_nTrk, binsWeightChi_nTrk);
	TH2D* h_chiAccEff_rat_nTrk = new TH2D("h_chiAccEff_rat_nTrk", "h_chiAccEff_rat_nTrk; |y_{lab}|; N_{trk}", nbinsWeightChi_absy, binsWeightChi_absy, nbinsWeightChi_nTrk, binsWeightChi_nTrk);


	// 1D weights for correcting the ratio of the chic to Jpsi raw yields

	TH1D* h_chiEfficiency1D_Q_pT_all_den = new TH1D("h_chiEfficiency1D_Q_pT_all_den", "h_chiEfficiency1D_Q_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_all_num = new TH1D("h_chiEfficiency1D_Q_pT_all_num", "h_chiEfficiency1D_Q_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_all_rat = new TH1D("h_chiEfficiency1D_Q_pT_all_rat", "h_chiEfficiency1D_Q_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_y_den = new TH1D("h_chiEfficiency1D_Q_y_den", "h_chiEfficiency1D_Q_y_den; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_Q_y_num = new TH1D("h_chiEfficiency1D_Q_y_num", "h_chiEfficiency1D_Q_y_num; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_Q_y_rat = new TH1D("h_chiEfficiency1D_Q_y_rat", "h_chiEfficiency1D_Q_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiEfficiency1D_QGen_y_den = new TH1D("h_chiEfficiency1D_QGen_y_den", "h_chiEfficiency1D_QGen_y_den; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_QGen_y_num = new TH1D("h_chiEfficiency1D_QGen_y_num", "h_chiEfficiency1D_QGen_y_num; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_QGen_y_rat = new TH1D("h_chiEfficiency1D_QGen_y_rat", "h_chiEfficiency1D_QGen_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiEfficiency1D_Q_nTrk_den = new TH1D("h_chiEfficiency1D_Q_nTrk_den", "h_chiEfficiency1D_Q_nTrk_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_Q_nTrk_num = new TH1D("h_chiEfficiency1D_Q_nTrk_num", "h_chiEfficiency1D_Q_nTrk_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_Q_nTrk_rat = new TH1D("h_chiEfficiency1D_Q_nTrk_rat", "h_chiEfficiency1D_Q_nTrk_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiEfficiency1D_Q_nTrk_all_den = new TH1D("h_chiEfficiency1D_Q_nTrk_all_den", "h_chiEfficiency1D_Q_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_Q_nTrk_all_num = new TH1D("h_chiEfficiency1D_Q_nTrk_all_num", "h_chiEfficiency1D_Q_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_Q_nTrk_all_rat = new TH1D("h_chiEfficiency1D_Q_nTrk_all_rat", "h_chiEfficiency1D_Q_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiEfficiency1D_Q_pT_mid_den = new TH1D("h_chiEfficiency1D_Q_pT_mid_den", "h_chiEfficiency1D_Q_pT_mid_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_mid_num = new TH1D("h_chiEfficiency1D_Q_pT_mid_num", "h_chiEfficiency1D_Q_pT_mid_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_mid_rat = new TH1D("h_chiEfficiency1D_Q_pT_mid_rat", "h_chiEfficiency1D_Q_pT_mid_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_fwd_den = new TH1D("h_chiEfficiency1D_Q_pT_fwd_den", "h_chiEfficiency1D_Q_pT_fwd_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwd_num = new TH1D("h_chiEfficiency1D_Q_pT_fwd_num", "h_chiEfficiency1D_Q_pT_fwd_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwd_rat = new TH1D("h_chiEfficiency1D_Q_pT_fwd_rat", "h_chiEfficiency1D_Q_pT_fwd_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_fwdOnly_den = new TH1D("h_chiEfficiency1D_Q_pT_fwdOnly_den", "h_chiEfficiency1D_Q_pT_fwdOnly_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwdOnly_num = new TH1D("h_chiEfficiency1D_Q_pT_fwdOnly_num", "h_chiEfficiency1D_Q_pT_fwdOnly_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwdOnly_rat = new TH1D("h_chiEfficiency1D_Q_pT_fwdOnly_rat", "h_chiEfficiency1D_Q_pT_fwdOnly_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_bkwOnly_den = new TH1D("h_chiEfficiency1D_Q_pT_bkwOnly_den", "h_chiEfficiency1D_Q_pT_bkwOnly_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_bkwOnly_num = new TH1D("h_chiEfficiency1D_Q_pT_bkwOnly_num", "h_chiEfficiency1D_Q_pT_bkwOnly_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_bkwOnly_rat = new TH1D("h_chiEfficiency1D_Q_pT_bkwOnly_rat", "h_chiEfficiency1D_Q_pT_bkwOnly_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_fwdOnlyWide_den = new TH1D("h_chiEfficiency1D_Q_pT_fwdOnlyWide_den", "h_chiEfficiency1D_Q_pT_fwdOnlyWide_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwdOnlyWide_num = new TH1D("h_chiEfficiency1D_Q_pT_fwdOnlyWide_num", "h_chiEfficiency1D_Q_pT_fwdOnlyWide_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwdOnlyWide_rat = new TH1D("h_chiEfficiency1D_Q_pT_fwdOnlyWide_rat", "h_chiEfficiency1D_Q_pT_fwdOnlyWide_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_bkwOnlyWide_den = new TH1D("h_chiEfficiency1D_Q_pT_bkwOnlyWide_den", "h_chiEfficiency1D_Q_pT_bkwOnlyWide_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_bkwOnlyWide_num = new TH1D("h_chiEfficiency1D_Q_pT_bkwOnlyWide_num", "h_chiEfficiency1D_Q_pT_bkwOnlyWide_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_bkwOnlyWide_rat = new TH1D("h_chiEfficiency1D_Q_pT_bkwOnlyWide_rat", "h_chiEfficiency1D_Q_pT_bkwOnlyWide_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_midCMS_den = new TH1D("h_chiEfficiency1D_Q_pT_midCMS_den", "h_chiEfficiency1D_Q_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_midCMS_num = new TH1D("h_chiEfficiency1D_Q_pT_midCMS_num", "h_chiEfficiency1D_Q_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_midCMS_rat = new TH1D("h_chiEfficiency1D_Q_pT_midCMS_rat", "h_chiEfficiency1D_Q_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_fwdCMS_den = new TH1D("h_chiEfficiency1D_Q_pT_fwdCMS_den", "h_chiEfficiency1D_Q_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwdCMS_num = new TH1D("h_chiEfficiency1D_Q_pT_fwdCMS_num", "h_chiEfficiency1D_Q_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_fwdCMS_rat = new TH1D("h_chiEfficiency1D_Q_pT_fwdCMS_rat", "h_chiEfficiency1D_Q_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_Q_pT_bkwCMS_den = new TH1D("h_chiEfficiency1D_Q_pT_bkwCMS_den", "h_chiEfficiency1D_Q_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_bkwCMS_num = new TH1D("h_chiEfficiency1D_Q_pT_bkwCMS_num", "h_chiEfficiency1D_Q_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_Q_pT_bkwCMS_rat = new TH1D("h_chiEfficiency1D_Q_pT_bkwCMS_rat", "h_chiEfficiency1D_Q_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);
	 
	h_chiEfficiency1D_Q_pT_all_den->Sumw2();
	h_chiEfficiency1D_Q_y_den->Sumw2();
	h_chiEfficiency1D_Q_nTrk_den->Sumw2();
	h_chiEfficiency1D_Q_nTrk_all_den->Sumw2();
	h_chiEfficiency1D_Q_pT_mid_den->Sumw2();
	h_chiEfficiency1D_Q_pT_fwd_den->Sumw2();
	h_chiEfficiency1D_Q_pT_fwdOnly_den->Sumw2();
	h_chiEfficiency1D_Q_pT_bkwOnly_den->Sumw2();
	h_chiEfficiency1D_Q_pT_fwdOnlyWide_den->Sumw2();
	h_chiEfficiency1D_Q_pT_bkwOnlyWide_den->Sumw2();
	h_chiEfficiency1D_Q_pT_midCMS_den->Sumw2();
	h_chiEfficiency1D_Q_pT_fwdCMS_den->Sumw2();
	h_chiEfficiency1D_Q_pT_bkwCMS_den->Sumw2();

	h_chiEfficiency1D_Q_pT_all_num->Sumw2();
	h_chiEfficiency1D_Q_y_num->Sumw2();
	h_chiEfficiency1D_Q_nTrk_num->Sumw2();
	h_chiEfficiency1D_Q_nTrk_all_num->Sumw2();
	h_chiEfficiency1D_Q_pT_mid_num->Sumw2();
	h_chiEfficiency1D_Q_pT_fwd_num->Sumw2();
	h_chiEfficiency1D_Q_pT_fwdOnly_num->Sumw2();
	h_chiEfficiency1D_Q_pT_bkwOnly_num->Sumw2();
	h_chiEfficiency1D_Q_pT_fwdOnlyWide_num->Sumw2();
	h_chiEfficiency1D_Q_pT_bkwOnlyWide_num->Sumw2();
	h_chiEfficiency1D_Q_pT_midCMS_num->Sumw2();
	h_chiEfficiency1D_Q_pT_fwdCMS_num->Sumw2();
	h_chiEfficiency1D_Q_pT_bkwCMS_num->Sumw2();

	h_chiEfficiency1D_QGen_y_den->Sumw2();
	h_chiEfficiency1D_QGen_y_num->Sumw2();


	// Same but no weights (systematics)

	TH1D* h_chiEfficiency1D_QnoW_pT_all_den = new TH1D("h_chiEfficiency1D_QnoW_pT_all_den", "h_chiEfficiency1D_QnoW_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_QnoW_pT_all_num = new TH1D("h_chiEfficiency1D_QnoW_pT_all_num", "h_chiEfficiency1D_QnoW_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_QnoW_pT_all_rat = new TH1D("h_chiEfficiency1D_QnoW_pT_all_rat", "h_chiEfficiency1D_QnoW_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_QnoW_y_den = new TH1D("h_chiEfficiency1D_QnoW_y_den", "h_chiEfficiency1D_QnoW_y_den; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_QnoW_y_num = new TH1D("h_chiEfficiency1D_QnoW_y_num", "h_chiEfficiency1D_QnoW_y_num; y", nbins_y, bins_y);
	TH1D* h_chiEfficiency1D_QnoW_y_rat = new TH1D("h_chiEfficiency1D_QnoW_y_rat", "h_chiEfficiency1D_QnoW_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiEfficiency1D_QnoW_nTrk_den = new TH1D("h_chiEfficiency1D_QnoW_nTrk_den", "h_chiEfficiency1D_QnoW_nTrk_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_QnoW_nTrk_num = new TH1D("h_chiEfficiency1D_QnoW_nTrk_num", "h_chiEfficiency1D_QnoW_nTrk_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiEfficiency1D_QnoW_nTrk_rat = new TH1D("h_chiEfficiency1D_QnoW_nTrk_rat", "h_chiEfficiency1D_QnoW_nTrk_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiEfficiency1D_QnoW_pT_mid_den = new TH1D("h_chiEfficiency1D_QnoW_pT_mid_den", "h_chiEfficiency1D_QnoW_pT_mid_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_QnoW_pT_mid_num = new TH1D("h_chiEfficiency1D_QnoW_pT_mid_num", "h_chiEfficiency1D_QnoW_pT_mid_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_QnoW_pT_mid_rat = new TH1D("h_chiEfficiency1D_QnoW_pT_mid_rat", "h_chiEfficiency1D_QnoW_pT_mid_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiEfficiency1D_QnoW_pT_fwd_den = new TH1D("h_chiEfficiency1D_QnoW_pT_fwd_den", "h_chiEfficiency1D_QnoW_pT_fwd_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_QnoW_pT_fwd_num = new TH1D("h_chiEfficiency1D_QnoW_pT_fwd_num", "h_chiEfficiency1D_QnoW_pT_fwd_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiEfficiency1D_QnoW_pT_fwd_rat = new TH1D("h_chiEfficiency1D_QnoW_pT_fwd_rat", "h_chiEfficiency1D_QnoW_pT_fwd_rat; p_{T}", nbins_pT, bins_pT);



	///////////////
	//   Acceptance corrections
	////////////

	TH1D* h_photAcceptanceCor1D_Q_pT_all_den = new TH1D("h_photAcceptanceCor1D_Q_pT_all_den", "h_photAcceptanceCor1D_Q_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_all_num = new TH1D("h_photAcceptanceCor1D_Q_pT_all_num", "h_photAcceptanceCor1D_Q_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_all_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_all_rat", "h_photAcceptanceCor1D_Q_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_y_den = new TH1D("h_photAcceptanceCor1D_Q_y_den", "h_photAcceptanceCor1D_Q_y_den; y", nbins_y, bins_y);
	TH1D* h_photAcceptanceCor1D_Q_y_num = new TH1D("h_photAcceptanceCor1D_Q_y_num", "h_photAcceptanceCor1D_Q_y_num; y", nbins_y, bins_y);
	TH1D* h_photAcceptanceCor1D_Q_y_rat = new TH1D("h_photAcceptanceCor1D_Q_y_rat", "h_photAcceptanceCor1D_Q_y_rat; y", nbins_y, bins_y);

	TH1D* h_photAcceptanceCor1D_QGen_y_den = new TH1D("h_photAcceptanceCor1D_QGen_y_den", "h_photAcceptanceCor1D_QGen_y_den; y", nbins_y, bins_y);
	TH1D* h_photAcceptanceCor1D_QGen_y_num = new TH1D("h_photAcceptanceCor1D_QGen_y_num", "h_photAcceptanceCor1D_QGen_y_num; y", nbins_y, bins_y);
	TH1D* h_photAcceptanceCor1D_QGen_y_rat = new TH1D("h_photAcceptanceCor1D_QGen_y_rat", "h_photAcceptanceCor1D_QGen_y_rat; y", nbins_y, bins_y);

	TH1D* h_photAcceptanceCor1D_Q_nTrk_den = new TH1D("h_photAcceptanceCor1D_Q_nTrk_den", "h_photAcceptanceCor1D_Q_nTrk_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_photAcceptanceCor1D_Q_nTrk_num = new TH1D("h_photAcceptanceCor1D_Q_nTrk_num", "h_photAcceptanceCor1D_Q_nTrk_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_photAcceptanceCor1D_Q_nTrk_rat = new TH1D("h_photAcceptanceCor1D_Q_nTrk_rat", "h_photAcceptanceCor1D_Q_nTrk_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_photAcceptanceCor1D_Q_nTrk_all_den = new TH1D("h_photAcceptanceCor1D_Q_nTrk_all_den", "h_photAcceptanceCor1D_Q_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_photAcceptanceCor1D_Q_nTrk_all_num = new TH1D("h_photAcceptanceCor1D_Q_nTrk_all_num", "h_photAcceptanceCor1D_Q_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_photAcceptanceCor1D_Q_nTrk_all_rat = new TH1D("h_photAcceptanceCor1D_Q_nTrk_all_rat", "h_photAcceptanceCor1D_Q_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_photAcceptanceCor1D_Q_pT_mid_den = new TH1D("h_photAcceptanceCor1D_Q_pT_mid_den", "h_photAcceptanceCor1D_Q_pT_mid_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_mid_num = new TH1D("h_photAcceptanceCor1D_Q_pT_mid_num", "h_photAcceptanceCor1D_Q_pT_mid_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_mid_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_mid_rat", "h_photAcceptanceCor1D_Q_pT_mid_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_fwd_den = new TH1D("h_photAcceptanceCor1D_Q_pT_fwd_den", "h_photAcceptanceCor1D_Q_pT_fwd_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwd_num = new TH1D("h_photAcceptanceCor1D_Q_pT_fwd_num", "h_photAcceptanceCor1D_Q_pT_fwd_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwd_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_fwd_rat", "h_photAcceptanceCor1D_Q_pT_fwd_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_fwdOnly_den = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdOnly_den", "h_photAcceptanceCor1D_Q_pT_fwdOnly_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwdOnly_num = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdOnly_num", "h_photAcceptanceCor1D_Q_pT_fwdOnly_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwdOnly_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdOnly_rat", "h_photAcceptanceCor1D_Q_pT_fwdOnly_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_bkwOnly_den = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwOnly_den", "h_photAcceptanceCor1D_Q_pT_bkwOnly_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_bkwOnly_num = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwOnly_num", "h_photAcceptanceCor1D_Q_pT_bkwOnly_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_bkwOnly_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwOnly_rat", "h_photAcceptanceCor1D_Q_pT_bkwOnly_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_den = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_den", "h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_num = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_num", "h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_rat", "h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_den = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_den", "h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_num = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_num", "h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_rat", "h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_midCMS_den = new TH1D("h_photAcceptanceCor1D_Q_pT_midCMS_den", "h_photAcceptanceCor1D_Q_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_midCMS_num = new TH1D("h_photAcceptanceCor1D_Q_pT_midCMS_num", "h_photAcceptanceCor1D_Q_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_midCMS_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_midCMS_rat", "h_photAcceptanceCor1D_Q_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_fwdCMS_den = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdCMS_den", "h_photAcceptanceCor1D_Q_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwdCMS_num = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdCMS_num", "h_photAcceptanceCor1D_Q_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_fwdCMS_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_fwdCMS_rat", "h_photAcceptanceCor1D_Q_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_photAcceptanceCor1D_Q_pT_bkwCMS_den = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwCMS_den", "h_photAcceptanceCor1D_Q_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_bkwCMS_num = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwCMS_num", "h_photAcceptanceCor1D_Q_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_photAcceptanceCor1D_Q_pT_bkwCMS_rat = new TH1D("h_photAcceptanceCor1D_Q_pT_bkwCMS_rat", "h_photAcceptanceCor1D_Q_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_photAcceptanceCor1D_Q_pT_all_den->Sumw2();
	h_photAcceptanceCor1D_Q_y_den->Sumw2();
	h_photAcceptanceCor1D_Q_nTrk_den->Sumw2();
	h_photAcceptanceCor1D_Q_nTrk_all_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_mid_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwd_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwdOnly_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_bkwOnly_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_midCMS_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwdCMS_den->Sumw2();
	h_photAcceptanceCor1D_Q_pT_bkwCMS_den->Sumw2();

	h_photAcceptanceCor1D_Q_pT_all_num->Sumw2();
	h_photAcceptanceCor1D_Q_y_num->Sumw2();
	h_photAcceptanceCor1D_Q_nTrk_num->Sumw2();
	h_photAcceptanceCor1D_Q_nTrk_all_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_mid_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwd_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwdOnly_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_bkwOnly_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_midCMS_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_fwdCMS_num->Sumw2();
	h_photAcceptanceCor1D_Q_pT_bkwCMS_num->Sumw2();


	TFile* fMC = new TFile(fileInMC, "READ");


	TTree* event_tree = (TTree*)fMC->Get("ChiRootuple/event_tree");
	if (!event_tree) { 
		cout << "Problem with event Tree";
		return 1;
	}

	LoadChiBranches(event_tree, true);

	long counterNormal = 0;
	long counterNormal2 = 0;
	long counterVeryLoose = 0;
	long counterMedium = 0;

	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: "<<nentries << endl;
	if (nentries > 500000) { nentries = 1000000; }


	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);

		//cout << "ispPb: " << ispPb << endl;
		//cout << "here" << endl;
		if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }

		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			weird_decay_counter++;
			continue;
		}

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
		// if (ispPb==false) rap_Jpsi = -rap_Jpsi;
		int muonPos1 = MuonMCMatched(0);
		int muonPos2 = MuonMCMatched(1);
		int JpsiPos = DimuonMCMatched(0);
		int matchPositionConv = PhotMCMatched(0);

		// quick correction ratio - for nominal results
		// logic: Only fiducial region (acceptance cuts passed). For muons the acceptance in reco, for photon the acceptance in gen


		bool isPassedOverall = false;
		bool isPassedDetailed = false;
		if (ChiPassAllCutsVeryLooseConversionMC(0) > -0.5) { counterVeryLoose++; }
		if (ChiPassAllCutsMediumConversionMC(0) > -0.5) { counterMedium++; }
		if (ChiPassAllCutsMC(0) > -0.5) {
			counterNormal2++; isPassedOverall = true;
		}



		if ((muonPos1 < 0 || muonPos2 < 0 || JpsiPos < 0)==false) {

			if (DimuonPassAllCuts(JpsiPos)) //Jpsi passed
			{

				double eta1Reco = muon_eta->at(muonPos1);
				double pt1Reco = muon_pt->at(muonPos1);
				double eta2Reco = muon_eta->at(muonPos2);
				double pt2Reco = muon_pt->at(muonPos2);
				//cout << "here " << muonPos1 << muonPos2 << JpsiPos<< endl;
				TLorentzVector* LVJpsiReco = (TLorentzVector*)dimuon_p4->At(JpsiPos);
				double rap_JpsiReco = LVJpsiReco->Rapidity();
				double pt_JpsiReco = dimuon_pt->at(JpsiPos);
				double eta_JpsiReco = dimuon_eta->at(JpsiPos);
				double nTrack_inPV = pvtx_nTracks->at(dimuon_pvtx_indexFromOniaMuMu->at(JpsiPos));

				double MCweight = h_weightPVtrk->GetBinContent(h_weightPVtrk->FindBin(nTrack_inPV));
				//cout << "Ntrack " << nTrack_inPV << " and weight: " << MCweight << endl;

				if (ispPb == true) { //direction dependent stuff
					//MCweight = 1.28*MCweight;
				}
				else
				{
					//MCweight = 0.72*MCweight;
					rap_JpsiReco = -rap_JpsiReco;
				}

				//MCweight = 1; //REMOVE, TEST


				if (PhotAcceptance(eta_phot, pt_phot) == true) { //cut on GEN rather than reco


					h_chiEfficiency1D_Q_pT_all_den->Fill(pt_JpsiReco, MCweight);
					h_chiEfficiency1D_Q_y_den->Fill(rap_JpsiReco, MCweight);
					h_chiEfficiency1D_QGen_y_den->Fill(rap_Jpsi, MCweight);
					if (fabs(rap_JpsiReco) < 1.0) { h_chiEfficiency1D_Q_nTrk_den->Fill(nTrack_inPV, MCweight); }
					h_chiEfficiency1D_Q_nTrk_all_den->Fill(nTrack_inPV, MCweight);
					if (fabs(rap_JpsiReco) < 1.0) { h_chiEfficiency1D_Q_pT_mid_den->Fill(pt_JpsiReco, MCweight); }
					else { h_chiEfficiency1D_Q_pT_fwd_den->Fill(pt_JpsiReco, MCweight); }

					if (rap_JpsiReco > 1.6) { h_chiEfficiency1D_Q_pT_fwdOnly_den->Fill(pt_JpsiReco, MCweight); }
					if (rap_JpsiReco > 1.0) { h_chiEfficiency1D_Q_pT_fwdOnlyWide_den->Fill(pt_JpsiReco, MCweight); }
					if (rap_JpsiReco < -1.6) { h_chiEfficiency1D_Q_pT_bkwOnly_den->Fill(pt_JpsiReco, MCweight); }
					if (rap_JpsiReco < -1.0) { h_chiEfficiency1D_Q_pT_bkwOnlyWide_den->Fill(pt_JpsiReco, MCweight); }
					if (rap_JpsiReco > -0.566 && rap_JpsiReco < 1.434) { h_chiEfficiency1D_Q_pT_midCMS_den->Fill(pt_JpsiReco, MCweight); }
					if (rap_JpsiReco < 2.4 && rap_JpsiReco > 1.434) { h_chiEfficiency1D_Q_pT_fwdCMS_den->Fill(pt_JpsiReco, MCweight); }
					if (rap_JpsiReco > -1.566 && rap_JpsiReco < -0.566) { h_chiEfficiency1D_Q_pT_bkwCMS_den->Fill(pt_JpsiReco, MCweight); }

					h_chiEfficiency1D_QnoW_pT_all_den->Fill(pt_JpsiReco);
					h_chiEfficiency1D_QnoW_y_den->Fill(rap_JpsiReco);
					h_chiEfficiency1D_QnoW_nTrk_den->Fill(nTrack_inPV);
					if (fabs(rap_JpsiReco) < 1.0) { h_chiEfficiency1D_QnoW_pT_mid_den->Fill(pt_JpsiReco); }
					else { h_chiEfficiency1D_QnoW_pT_fwd_den->Fill(pt_JpsiReco); }

				}

				if (matchPositionConv > -1) {
					if (PhotAcceptance(conv_eta->at(matchPositionConv), conv_pt->at(matchPositionConv)) == true) {

						int chicPos = ChiMCMatched(0);
						if (chicPos > -0.5) {
							if (ChiPassAllCuts(chicPos) == true) //conversion cuts + dimuon mass + any other cuts on chic 
							{
								counterNormal++;
								isPassedDetailed = true;
								//cout << counterNormal << "   weight: " << MCweight << endl;

								TLorentzVector* LVChicReco = (TLorentzVector*)chi_p4->At(chicPos);
								double rap_chicReco = LVChicReco->Rapidity();
								double pt_chicReco = chi_pt->at(chicPos);
								h_chiEfficiency1D_Q_pT_all_num->Fill(pt_JpsiReco, MCweight);
								h_chiEfficiency1D_Q_y_num->Fill(rap_JpsiReco, MCweight);
								h_chiEfficiency1D_QGen_y_num->Fill(rap_Jpsi, MCweight);
								if (fabs(rap_JpsiReco) < 1.0) { h_chiEfficiency1D_Q_nTrk_num->Fill(nTrack_inPV, MCweight); }
								h_chiEfficiency1D_Q_nTrk_all_num->Fill(nTrack_inPV, MCweight);
								if (fabs(rap_JpsiReco) < 1.0) { h_chiEfficiency1D_Q_pT_mid_num->Fill(pt_JpsiReco, MCweight); }
								else { h_chiEfficiency1D_Q_pT_fwd_num->Fill(pt_JpsiReco, MCweight); }

								if (rap_JpsiReco > 1.6) { h_chiEfficiency1D_Q_pT_fwdOnly_num->Fill(pt_JpsiReco, MCweight); }
								if (rap_JpsiReco < -1.6) { h_chiEfficiency1D_Q_pT_bkwOnly_num->Fill(pt_JpsiReco, MCweight); }
								if (rap_JpsiReco > 1.0) { h_chiEfficiency1D_Q_pT_fwdOnlyWide_num->Fill(pt_JpsiReco, MCweight); }
								if (rap_JpsiReco < -1.0) { h_chiEfficiency1D_Q_pT_bkwOnlyWide_num->Fill(pt_JpsiReco, MCweight); }
								if (rap_JpsiReco > -0.566 && rap_JpsiReco < 1.434) { h_chiEfficiency1D_Q_pT_midCMS_num->Fill(pt_JpsiReco, MCweight); }
								if (rap_JpsiReco < 2.4 && rap_JpsiReco > 1.434) { h_chiEfficiency1D_Q_pT_fwdCMS_num->Fill(pt_JpsiReco, MCweight); }
								if (rap_JpsiReco > -1.566 && rap_JpsiReco < -0.566) { h_chiEfficiency1D_Q_pT_bkwCMS_num->Fill(pt_JpsiReco, MCweight); }

								h_chiEfficiency1D_QnoW_pT_all_num->Fill(pt_JpsiReco);
								h_chiEfficiency1D_QnoW_y_num->Fill(rap_JpsiReco);
								h_chiEfficiency1D_QnoW_nTrk_num->Fill(nTrack_inPV);
								if (fabs(rap_JpsiReco) < 1.0) { h_chiEfficiency1D_QnoW_pT_mid_num->Fill(pt_JpsiReco); }
								else { h_chiEfficiency1D_QnoW_pT_fwd_num->Fill(pt_JpsiReco); }
							}
						}
					}
				}
			}
		}

		if (isPassedDetailed == false && isPassedOverall == true)
		{
			cout << "Dimuon cuts: " << DimuonPassAllCuts(JpsiPos) << "   phot acceptance: " << PhotAcceptance(eta_phot, pt_phot) << " and cut: " << PhotSelectionPassMC(0) << endl;
		}
		
	
		// ACCEPTANCE

		double nTrack_inPV = ntracks_inEvent; //here use total ntrack, since unless there is reconstructed J/psi, we don't know what vertex to chose (total is the greatest PV)

		h_muonAcceptance2D_den->Fill(eta1, pt1);
		h_muonAcceptance2D_den->Fill(eta2, pt2);
		if (MuonAcceptance(eta1, pt1) == true) h_muonAcceptance2D_num->Fill(eta1, pt1);
		if (MuonAcceptance(eta2, pt2) == true) h_muonAcceptance2D_num->Fill(eta2, pt2);

		h_photAcceptance2D_den->Fill(eta_phot, pt_phot);
		if (PhotAcceptance(eta_phot, pt_phot) == true) h_photAcceptance2D_num->Fill(eta_phot, pt_phot);

		double eta_chi = gen_chic_eta->at(0);
		double pt_chi = gen_chic_pt->at(0);
		TLorentzVector* LVchic = (TLorentzVector*)gen_chic_p4->At(0);
		double rap_chi = LVchic->Rapidity();

		h_JpsiAcceptance2D_den->Fill(eta_Jpsi, pt_Jpsi);
		h_JpsiAcceptance2D_y_den->Fill(rap_Jpsi, pt_Jpsi);
		h_chiAcceptance2D_den->Fill(eta_chi, pt_chi);
		h_chiAcceptance2D_y_den->Fill(rap_chi, pt_chi);

		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && DimuonAcceptance(rap_Jpsi, pt_Jpsi) == true)
		{
			h_JpsiAcceptance2D_num->Fill(eta_Jpsi, pt_Jpsi);
			h_JpsiAcceptance2D_y_num->Fill(rap_Jpsi, pt_Jpsi);

			h_photDistribution3D_JpsipT->Fill(pt_Jpsi, eta_phot, pt_phot);
			h_photDistribution3D_Jpsiy->Fill(rap_Jpsi, eta_phot, pt_phot);

			// acceptance corrections
			
			double rap_Jpsi_corrections = rap_Jpsi;
			double MCweight = h_weightPVtrk->GetBinContent(h_weightPVtrk->FindBin(nTrack_inPV));
			//cout << "Ntrack " << nTrack_inPV << " and weight: " << MCweight << endl;

			if (ispPb == true) { //direction dependent stuff
				MCweight = 1.28*MCweight;
			}
			else
			{
				MCweight = 0.72*MCweight;
				rap_Jpsi_corrections = -rap_Jpsi_corrections;
			}


			h_photAcceptanceCor1D_Q_pT_all_den->Fill(pt_Jpsi, MCweight);
			h_photAcceptanceCor1D_Q_y_den->Fill(rap_Jpsi_corrections, MCweight);
			if (fabs(rap_Jpsi_corrections) < 1.0) {h_photAcceptanceCor1D_Q_nTrk_den->Fill(nTrack_inPV, MCweight);	}
			h_photAcceptanceCor1D_Q_nTrk_all_den->Fill(nTrack_inPV, MCweight);
			if (fabs(rap_Jpsi_corrections) < 1.0) {h_photAcceptanceCor1D_Q_pT_mid_den->Fill(pt_Jpsi, MCweight); }
			else {	h_photAcceptanceCor1D_Q_pT_fwd_den->Fill(pt_Jpsi, MCweight); }
			if (rap_Jpsi_corrections > 1.6) { h_photAcceptanceCor1D_Q_pT_fwdOnly_den->Fill(pt_Jpsi, MCweight); }
			if (rap_Jpsi_corrections < -1.6) { h_photAcceptanceCor1D_Q_pT_bkwOnly_den->Fill(pt_Jpsi, MCweight); }
			if (rap_Jpsi_corrections > 1.0) { h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_den->Fill(pt_Jpsi, MCweight); }
			if (rap_Jpsi_corrections < -1.0) { h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_den->Fill(pt_Jpsi, MCweight); }
			if (rap_Jpsi_corrections > -0.566 && rap_Jpsi_corrections < 1.434) { h_photAcceptanceCor1D_Q_pT_midCMS_den->Fill(pt_Jpsi, MCweight); }
			if (rap_Jpsi_corrections < 2.4 && rap_Jpsi_corrections > 1.434) { h_photAcceptanceCor1D_Q_pT_fwdCMS_den->Fill(pt_Jpsi, MCweight); }
			if (rap_Jpsi_corrections > -1.566 && rap_Jpsi_corrections < -0.566) { h_photAcceptanceCor1D_Q_pT_bkwCMS_den->Fill(pt_Jpsi, MCweight); }



			if (PhotAcceptance(eta_phot, pt_phot) == true)
			{
				h_chiAcceptance2D_num->Fill(eta_chi, pt_chi);
				h_chiAcceptance2D_y_num->Fill(rap_chi, pt_chi);

				h_photAcceptanceCor1D_Q_pT_all_num->Fill(pt_Jpsi, MCweight);
				h_photAcceptanceCor1D_Q_y_num->Fill(rap_Jpsi_corrections, MCweight);
				if (fabs(rap_Jpsi_corrections) < 1.0) h_photAcceptanceCor1D_Q_nTrk_num->Fill(nTrack_inPV, MCweight);
				h_photAcceptanceCor1D_Q_nTrk_all_num->Fill(nTrack_inPV, MCweight);
				if (fabs(rap_Jpsi_corrections) < 1.0) { h_photAcceptanceCor1D_Q_pT_mid_num->Fill(pt_Jpsi, MCweight); }
				else { h_photAcceptanceCor1D_Q_pT_fwd_num->Fill(pt_Jpsi, MCweight); }
				if (rap_Jpsi_corrections > 1.6) h_photAcceptanceCor1D_Q_pT_fwdOnly_num->Fill(pt_Jpsi, MCweight);
				if (rap_Jpsi_corrections < -1.6) h_photAcceptanceCor1D_Q_pT_bkwOnly_num->Fill(pt_Jpsi, MCweight);
				if (rap_Jpsi_corrections > 1.0) h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_num->Fill(pt_Jpsi, MCweight);
				if (rap_Jpsi_corrections < -1.0) h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_num->Fill(pt_Jpsi, MCweight);
				if (rap_Jpsi_corrections > -0.566 && rap_Jpsi_corrections < 1.434) h_photAcceptanceCor1D_Q_pT_midCMS_num->Fill(pt_Jpsi, MCweight);
				if (rap_Jpsi_corrections < 2.4 && rap_Jpsi_corrections > 1.434) h_photAcceptanceCor1D_Q_pT_fwdCMS_num->Fill(pt_Jpsi, MCweight);
				if (rap_Jpsi_corrections > -1.566 && rap_Jpsi_corrections < -0.566) h_photAcceptanceCor1D_Q_pT_bkwCMS_num->Fill(pt_Jpsi, MCweight);

			}
		}


		// EFFICIENCY
		if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && DimuonAcceptance(rap_Jpsi, pt_Jpsi) == true) //do only if J/psi in acceptance
		{
			h_muonEfficiency2D_den->Fill(eta1, pt1);
			h_muonEfficiency2D_den->Fill(eta2, pt2);
			if (MuonSelectionPassMC(0) == true) h_muonEfficiency2D_num->Fill(eta1, pt1);
			if (MuonSelectionPassMC(1) == true) h_muonEfficiency2D_num->Fill(eta2, pt2);

			h_JpsiEfficiency2D_y_den->Fill(rap_Jpsi, pt_Jpsi);
			h_JpsiEfficiency1D_den->Fill(pt_Jpsi);
			h_JpsiEfficiency1D_y_den->Fill(rap_Jpsi);
			h_JpsiEfficiency1D_nTrack_den->Fill(nTrack_inPV);
			if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0) == true) {
				h_JpsiEfficiency2D_y_num->Fill(rap_Jpsi, pt_Jpsi);
				h_JpsiEfficiency1D_num->Fill(pt_Jpsi);
				h_JpsiEfficiency1D_y_num->Fill(rap_Jpsi);
				h_JpsiEfficiency1D_nTrack_num->Fill(nTrack_inPV);
			}

			if (PhotAcceptance(eta_phot, pt_phot) == true) {
				h_photEfficiency2D_den->Fill(eta_phot, pt_phot);
				h_photEfficiency1D_den->Fill(pt_phot);
				h_photEfficiency1D_y_den->Fill(eta_phot);
				h_photEfficiency1D_nTrack_den->Fill(nTrack_inPV);
				h_photEfficiency1D_y_denJpsiVar->Fill(rap_Jpsi);

				if (matchPositionConv > -0.5)
				{
					h_photEfficiency_y_relative_num_noCut->Fill(eta_phot);
					h_photEfficiency_y_relative_num_noCutJpsiVar->Fill(rap_Jpsi);

					if (PhotSelectionPass(matchPositionConv) == true)
					{
						h_photEfficiency_y_relative_num_fullCut->Fill(eta_phot);
						h_photEfficiency_y_relative_num_fullCutJpsiVar->Fill(rap_Jpsi);
					}
				}

				if (PhotSelectionPassMC(0) == true) {
					h_photEfficiency2D_num->Fill(eta_phot, pt_phot);
					h_photEfficiency1D_num->Fill(pt_phot);
					h_photEfficiency1D_y_num->Fill(eta_phot);
					h_photEfficiency1D_nTrack_num->Fill(nTrack_inPV);

				}


				
				if (MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0) == true)
				{
					if (matchPositionConv > -0.5) {// y effect crosscheck
						h_photEfficiency_y_relative_den->Fill(rap_chi);
						if (PhotSelectionPassLoose(matchPositionConv) == true) h_photEfficiency_y_relative_num_loose->Fill(rap_chi);
						if (PhotSelectionPassMedium(matchPositionConv) == true) h_photEfficiency_y_relative_num_medium->Fill(rap_chi);
						if (PhotSelectionPassTight(matchPositionConv) == true) h_photEfficiency_y_relative_num_tight->Fill(rap_chi);
					}

					if (PhotSelectionPassMC(0) == true) {
						h_chiEfficiency2D_y_num->Fill(rap_chi, pt_chi);
						h_chiEfficiency1D_num->Fill(pt_chi);
						h_chiEfficiency1D_y_num->Fill(rap_chi);
						h_chiEfficiency1D_nTrack_num->Fill(nTrack_inPV);
					}
				}

				h_chiEfficiency2D_y_den->Fill(rap_chi, pt_chi);
				h_chiEfficiency1D_den->Fill(pt_chi);
				h_chiEfficiency1D_y_den->Fill(rap_chi);
				h_chiEfficiency1D_nTrack_den->Fill(nTrack_inPV);
			}
		}

		//*/
		////weights (<acc*eff>, then to be used as 1/w)
		//h_JpsiAccEff_den->Fill(abs(rap_Jpsi), pt_Jpsi);
		//h_JpsiAccEff_den_nTrk->Fill(abs(rap_Jpsi), ntracks_inEvent);
		//h_chiAccEff_den->Fill(abs(rap_chi), pt_chi);
		//h_chiAccEff_den_nTrk->Fill(abs(rap_chi), ntracks_inEvent);

		//if (MuonAcceptance(eta1, pt1) == true && MuonAcceptance(eta2, pt2) == true && MuonSelectionPassMC(0) == true && MuonSelectionPassMC(1) == true && DimuonSelectionPassMC(0) == true)
		//{
		//	h_JpsiAccEff_num->Fill(abs(rap_Jpsi), pt_Jpsi);
		//	h_JpsiAccEff_num_nTrk->Fill(abs(rap_Jpsi), ntracks_inEvent);

		//	if (PhotAcceptance(eta_phot, pt_phot) == true && PhotSelectionPassMC(0) == true)
		//	{
		//		h_chiAccEff_num->Fill(abs(rap_chi), pt_chi);
		//		h_chiAccEff_num_nTrk->Fill(abs(rap_chi), ntracks_inEvent);
		//	}
		//}


		

	} //end of event loop

	cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;


	cout << "Counters - very loose " << counterVeryLoose << "  normal: " << counterNormal << "  medium: " << counterMedium << endl;
	cout << "Counters - normal 2 : " << counterNormal2  << endl;



	// quick nominal corrections:
	

	h_chiEfficiency1D_Q_pT_all_rat->Divide(h_chiEfficiency1D_Q_pT_all_num, h_chiEfficiency1D_Q_pT_all_den, 1, 1, "B");
	h_chiEfficiency1D_Q_y_rat->Divide(h_chiEfficiency1D_Q_y_num, h_chiEfficiency1D_Q_y_den, 1, 1, "B");
	h_chiEfficiency1D_Q_nTrk_rat->Divide(h_chiEfficiency1D_Q_nTrk_num, h_chiEfficiency1D_Q_nTrk_den, 1, 1, "B");
	h_chiEfficiency1D_Q_nTrk_all_rat->Divide(h_chiEfficiency1D_Q_nTrk_all_num, h_chiEfficiency1D_Q_nTrk_all_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_mid_rat->Divide(h_chiEfficiency1D_Q_pT_mid_num, h_chiEfficiency1D_Q_pT_mid_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_fwd_rat->Divide(h_chiEfficiency1D_Q_pT_fwd_num, h_chiEfficiency1D_Q_pT_fwd_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_fwdOnly_rat->Divide(h_chiEfficiency1D_Q_pT_fwdOnly_num, h_chiEfficiency1D_Q_pT_fwdOnly_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_bkwOnly_rat->Divide(h_chiEfficiency1D_Q_pT_bkwOnly_num, h_chiEfficiency1D_Q_pT_bkwOnly_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_fwdOnlyWide_rat->Divide(h_chiEfficiency1D_Q_pT_fwdOnlyWide_num, h_chiEfficiency1D_Q_pT_fwdOnlyWide_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_bkwOnlyWide_rat->Divide(h_chiEfficiency1D_Q_pT_bkwOnlyWide_num, h_chiEfficiency1D_Q_pT_bkwOnlyWide_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_midCMS_rat->Divide(h_chiEfficiency1D_Q_pT_midCMS_num, h_chiEfficiency1D_Q_pT_midCMS_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_fwdCMS_rat->Divide(h_chiEfficiency1D_Q_pT_fwdCMS_num, h_chiEfficiency1D_Q_pT_fwdCMS_den, 1, 1, "B");
	h_chiEfficiency1D_Q_pT_bkwCMS_rat->Divide(h_chiEfficiency1D_Q_pT_bkwCMS_num, h_chiEfficiency1D_Q_pT_bkwCMS_den, 1, 1, "B");

	h_chiEfficiency1D_QGen_y_rat->Divide(h_chiEfficiency1D_QGen_y_num, h_chiEfficiency1D_QGen_y_den, 1, 1, "B");

	h_chiEfficiency1D_QnoW_pT_all_den->Sumw2();
	h_chiEfficiency1D_QnoW_y_den->Sumw2();
	h_chiEfficiency1D_QnoW_nTrk_den->Sumw2();
	h_chiEfficiency1D_QnoW_pT_mid_den->Sumw2();
	h_chiEfficiency1D_QnoW_pT_fwd_den->Sumw2();

	h_chiEfficiency1D_QnoW_pT_all_num->Sumw2();
	h_chiEfficiency1D_QnoW_y_num->Sumw2();
	h_chiEfficiency1D_QnoW_nTrk_num->Sumw2();
	h_chiEfficiency1D_QnoW_pT_mid_num->Sumw2();
	h_chiEfficiency1D_QnoW_pT_fwd_num->Sumw2();

	h_chiEfficiency1D_QnoW_pT_all_rat->Divide(h_chiEfficiency1D_QnoW_pT_all_num, h_chiEfficiency1D_QnoW_pT_all_den, 1, 1, "B");
	h_chiEfficiency1D_QnoW_y_rat->Divide(h_chiEfficiency1D_QnoW_y_num, h_chiEfficiency1D_QnoW_y_den, 1, 1, "B");
	h_chiEfficiency1D_QnoW_nTrk_rat->Divide(h_chiEfficiency1D_QnoW_nTrk_num, h_chiEfficiency1D_QnoW_nTrk_den, 1, 1, "B");
	h_chiEfficiency1D_QnoW_pT_mid_rat->Divide(h_chiEfficiency1D_QnoW_pT_mid_num, h_chiEfficiency1D_QnoW_pT_mid_den, 1, 1, "B");
	h_chiEfficiency1D_QnoW_pT_fwd_rat->Divide(h_chiEfficiency1D_QnoW_pT_fwd_num, h_chiEfficiency1D_QnoW_pT_fwd_den, 1, 1, "B");


	// quick acceptance corrections:

	h_photAcceptanceCor1D_Q_pT_all_rat->Divide(h_photAcceptanceCor1D_Q_pT_all_num, h_photAcceptanceCor1D_Q_pT_all_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_y_rat->Divide(h_photAcceptanceCor1D_Q_y_num, h_photAcceptanceCor1D_Q_y_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_nTrk_rat->Divide(h_photAcceptanceCor1D_Q_nTrk_num, h_photAcceptanceCor1D_Q_nTrk_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_nTrk_all_rat->Divide(h_photAcceptanceCor1D_Q_nTrk_all_num, h_photAcceptanceCor1D_Q_nTrk_all_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_mid_rat->Divide(h_photAcceptanceCor1D_Q_pT_mid_num, h_photAcceptanceCor1D_Q_pT_mid_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_fwd_rat->Divide(h_photAcceptanceCor1D_Q_pT_fwd_num, h_photAcceptanceCor1D_Q_pT_fwd_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_fwdOnly_rat->Divide(h_photAcceptanceCor1D_Q_pT_fwdOnly_num, h_photAcceptanceCor1D_Q_pT_fwdOnly_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_bkwOnly_rat->Divide(h_photAcceptanceCor1D_Q_pT_bkwOnly_num, h_photAcceptanceCor1D_Q_pT_bkwOnly_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_rat->Divide(h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_num, h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_rat->Divide(h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_num, h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_midCMS_rat->Divide(h_photAcceptanceCor1D_Q_pT_midCMS_num, h_photAcceptanceCor1D_Q_pT_midCMS_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_fwdCMS_rat->Divide(h_photAcceptanceCor1D_Q_pT_fwdCMS_num, h_photAcceptanceCor1D_Q_pT_fwdCMS_den, 1, 1, "B");
	h_photAcceptanceCor1D_Q_pT_bkwCMS_rat->Divide(h_photAcceptanceCor1D_Q_pT_bkwCMS_num, h_photAcceptanceCor1D_Q_pT_bkwCMS_den, 1, 1, "B");







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

	h_JpsiEfficiency1D_num->Sumw2();
	h_JpsiEfficiency1D_den->Sumw2();
	h_chiEfficiency1D_num->Sumw2();
	h_chiEfficiency1D_den->Sumw2();

	h_photEfficiency1D_y_num->Sumw2();
	h_photEfficiency1D_y_den->Sumw2();
	h_photEfficiency1D_y_denJpsiVar->Sumw2();
	h_JpsiEfficiency1D_y_num->Sumw2();
	h_JpsiEfficiency1D_y_den->Sumw2();
	h_chiEfficiency1D_y_num->Sumw2();
	h_chiEfficiency1D_y_den->Sumw2();

	h_photEfficiency1D_rat->Divide(h_photEfficiency1D_num, h_photEfficiency1D_den, 1, 1, "B");
	h_JpsiEfficiency1D_rat->Divide(h_JpsiEfficiency1D_num, h_JpsiEfficiency1D_den, 1, 1, "B");
	h_chiEfficiency1D_rat->Divide(h_chiEfficiency1D_num, h_chiEfficiency1D_den, 1, 1, "B");

	h_photEfficiency1D_y_rat->Divide(h_photEfficiency1D_y_num, h_photEfficiency1D_y_den, 1, 1, "B");
	h_JpsiEfficiency1D_y_rat->Divide(h_JpsiEfficiency1D_y_num, h_JpsiEfficiency1D_y_den, 1, 1, "B");
	h_chiEfficiency1D_y_rat->Divide(h_chiEfficiency1D_y_num, h_chiEfficiency1D_y_den, 1, 1, "B");
	

	h_photEfficiency1D_nTrack_den->Sumw2();
	h_photEfficiency1D_nTrack_num->Sumw2();
	h_JpsiEfficiency1D_nTrack_den->Sumw2();
	h_JpsiEfficiency1D_nTrack_num->Sumw2();
	h_chiEfficiency1D_nTrack_den->Sumw2();
	h_chiEfficiency1D_nTrack_num->Sumw2();

	h_photEfficiency1D_nTrack_rat->Divide(h_photEfficiency1D_nTrack_num, h_photEfficiency1D_nTrack_den, 1, 1, "B");
	h_JpsiEfficiency1D_nTrack_rat->Divide(h_JpsiEfficiency1D_nTrack_num, h_JpsiEfficiency1D_nTrack_den, 1, 1, "B");
	h_chiEfficiency1D_nTrack_rat->Divide(h_chiEfficiency1D_nTrack_num, h_chiEfficiency1D_nTrack_den, 1, 1, "B");


	h_photEfficiency_y_relative_den->Sumw2();
	h_photEfficiency_y_relative_num_loose->Sumw2();
	h_photEfficiency_y_relative_num_medium->Sumw2();
	h_photEfficiency_y_relative_num_tight->Sumw2();

	h_photEfficiency_y_relative_num_noCut->Sumw2();
	h_photEfficiency_y_relative_num_noCutJpsiVar->Sumw2();
	h_photEfficiency_y_relative_num_fullCut->Sumw2();
	h_photEfficiency_y_relative_num_fullCutJpsiVar->Sumw2();

	h_photEfficiency_y_relative_rat_loose->Divide(h_photEfficiency_y_relative_num_loose, h_photEfficiency_y_relative_den, 1, 1, "B");
	h_photEfficiency_y_relative_rat_medium->Divide(h_photEfficiency_y_relative_num_medium, h_photEfficiency_y_relative_den, 1, 1, "B");
	h_photEfficiency_y_relative_rat_tight->Divide(h_photEfficiency_y_relative_num_tight, h_photEfficiency_y_relative_den, 1, 1, "B");

	h_photEfficiency_y_relative_rat_noCut->Divide(h_photEfficiency_y_relative_num_noCut, h_photEfficiency1D_y_den, 1, 1, "B");
	h_photEfficiency_y_relative_rat_fullCut->Divide(h_photEfficiency_y_relative_num_fullCut, h_photEfficiency_y_relative_num_noCut, 1, 1, "B");

	h_photEfficiency_y_relative_rat_noCutJpsiVar->Divide(h_photEfficiency_y_relative_num_noCutJpsiVar, h_photEfficiency1D_y_denJpsiVar, 1, 1, "B");
	h_photEfficiency_y_relative_rat_fullCutJpsiVar->Divide(h_photEfficiency_y_relative_num_fullCutJpsiVar, h_photEfficiency_y_relative_num_noCutJpsiVar, 1, 1, "B");

	h_JpsiAccEff_den->Sumw2();
	h_JpsiAccEff_num->Sumw2();
	h_JpsiAccEff_den_nTrk->Sumw2();
	h_JpsiAccEff_num_nTrk->Sumw2();
	h_chiAccEff_den->Sumw2();
	h_chiAccEff_num->Sumw2();
	h_chiAccEff_den_nTrk->Sumw2();
	h_chiAccEff_num_nTrk->Sumw2();

	h_JpsiAccEff_rat->Divide(h_JpsiAccEff_num, h_JpsiAccEff_den);
	h_JpsiAccEff_rat_nTrk->Divide(h_JpsiAccEff_num_nTrk, h_JpsiAccEff_den_nTrk);
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
	

	TCanvas* can2 = new TCanvas("can2", "plot", 1200, 800);
	//TH1D* hOne = new TH1D("hOne", "hOne", nbins_y, bins_y);
	//hOne->SetMaximum(1.1);
	//hOne->SetMinimum(0);
	//hOne->Draw("");
	h_photEfficiency_y_relative_rat_loose->SetLineColor(kBlue);
	h_photEfficiency_y_relative_rat_medium->SetLineColor(kGreen);
	h_photEfficiency_y_relative_rat_tight->SetLineColor(kRed);
	h_photEfficiency_y_relative_rat_loose->SetMarkerColor(kBlue);
	h_photEfficiency_y_relative_rat_medium->SetMarkerColor(kGreen);
	h_photEfficiency_y_relative_rat_tight->SetMarkerColor(kRed);
	h_photEfficiency_y_relative_rat_loose->Draw();
	h_photEfficiency_y_relative_rat_medium->Draw("same");
	h_photEfficiency_y_relative_rat_tight->Draw("same");








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



	h_chiEfficiency1D_Q_pT_all_den->Write();
	h_chiEfficiency1D_Q_y_den->Write();
	h_chiEfficiency1D_Q_nTrk_den->Write();
	h_chiEfficiency1D_Q_nTrk_all_den->Write();
	h_chiEfficiency1D_Q_pT_mid_den->Write();
	h_chiEfficiency1D_Q_pT_fwd_den->Write();
	h_chiEfficiency1D_Q_pT_fwdOnly_den->Write();
	h_chiEfficiency1D_Q_pT_bkwOnly_den->Write();
	h_chiEfficiency1D_Q_pT_fwdOnlyWide_den->Write();
	h_chiEfficiency1D_Q_pT_bkwOnlyWide_den->Write();
	h_chiEfficiency1D_Q_pT_midCMS_den->Write();
	h_chiEfficiency1D_Q_pT_fwdCMS_den->Write();
	h_chiEfficiency1D_Q_pT_bkwCMS_den->Write();

	h_chiEfficiency1D_Q_pT_all_num->Write();
	h_chiEfficiency1D_Q_y_num->Write();
	h_chiEfficiency1D_Q_nTrk_num->Write();
	h_chiEfficiency1D_Q_nTrk_all_num->Write();
	h_chiEfficiency1D_Q_pT_mid_num->Write();
	h_chiEfficiency1D_Q_pT_fwd_num->Write();
	h_chiEfficiency1D_Q_pT_fwdOnly_num->Write();
	h_chiEfficiency1D_Q_pT_bkwOnly_num->Write();
	h_chiEfficiency1D_Q_pT_fwdOnlyWide_num->Write();
	h_chiEfficiency1D_Q_pT_bkwOnlyWide_num->Write();
	h_chiEfficiency1D_Q_pT_midCMS_num->Write();
	h_chiEfficiency1D_Q_pT_fwdCMS_num->Write();
	h_chiEfficiency1D_Q_pT_bkwCMS_num->Write();

	h_chiEfficiency1D_Q_pT_all_rat->Write();
	h_chiEfficiency1D_Q_y_rat->Write();
	h_chiEfficiency1D_Q_nTrk_rat->Write();
	h_chiEfficiency1D_Q_nTrk_all_rat->Write();
	h_chiEfficiency1D_Q_pT_mid_rat->Write();
	h_chiEfficiency1D_Q_pT_fwd_rat->Write();
	h_chiEfficiency1D_Q_pT_fwdOnly_rat->Write();
	h_chiEfficiency1D_Q_pT_bkwOnly_rat->Write();
	h_chiEfficiency1D_Q_pT_fwdOnlyWide_rat->Write();
	h_chiEfficiency1D_Q_pT_bkwOnlyWide_rat->Write();
	h_chiEfficiency1D_Q_pT_midCMS_rat->Write();
	h_chiEfficiency1D_Q_pT_fwdCMS_rat->Write();
	h_chiEfficiency1D_Q_pT_bkwCMS_rat->Write();

	h_chiEfficiency1D_QGen_y_den->Write();
	h_chiEfficiency1D_QGen_y_num->Write();
	h_chiEfficiency1D_QGen_y_rat->Write();


	h_chiEfficiency1D_QnoW_pT_all_den->Write();
	h_chiEfficiency1D_QnoW_y_den->Write();
	h_chiEfficiency1D_QnoW_nTrk_den->Write();
	h_chiEfficiency1D_QnoW_pT_mid_den->Write();
	h_chiEfficiency1D_QnoW_pT_fwd_den->Write();

	h_chiEfficiency1D_QnoW_pT_all_num->Write();
	h_chiEfficiency1D_QnoW_y_num->Write();
	h_chiEfficiency1D_QnoW_nTrk_num->Write();
	h_chiEfficiency1D_QnoW_pT_mid_num->Write();
	h_chiEfficiency1D_QnoW_pT_fwd_num->Write();

	h_chiEfficiency1D_QnoW_pT_all_rat->Write();
	h_chiEfficiency1D_QnoW_y_rat->Write();
	h_chiEfficiency1D_QnoW_nTrk_rat->Write();
	h_chiEfficiency1D_QnoW_pT_mid_rat->Write();
	h_chiEfficiency1D_QnoW_pT_fwd_rat->Write();

	h_photAcceptanceCor1D_Q_pT_all_den->Write();
	h_photAcceptanceCor1D_Q_y_den->Write();
	h_photAcceptanceCor1D_Q_nTrk_den->Write();
	h_photAcceptanceCor1D_Q_nTrk_all_den->Write();
	h_photAcceptanceCor1D_Q_pT_mid_den->Write();
	h_photAcceptanceCor1D_Q_pT_fwd_den->Write();
	h_photAcceptanceCor1D_Q_pT_fwdOnly_den->Write();
	h_photAcceptanceCor1D_Q_pT_bkwOnly_den->Write();
	h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_den->Write();
	h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_den->Write();
	h_photAcceptanceCor1D_Q_pT_midCMS_den->Write();
	h_photAcceptanceCor1D_Q_pT_fwdCMS_den->Write();
	h_photAcceptanceCor1D_Q_pT_bkwCMS_den->Write();

	h_photAcceptanceCor1D_Q_pT_all_num->Write();
	h_photAcceptanceCor1D_Q_y_num->Write();
	h_photAcceptanceCor1D_Q_nTrk_num->Write();
	h_photAcceptanceCor1D_Q_nTrk_all_num->Write();
	h_photAcceptanceCor1D_Q_pT_mid_num->Write();
	h_photAcceptanceCor1D_Q_pT_fwd_num->Write();
	h_photAcceptanceCor1D_Q_pT_fwdOnly_num->Write();
	h_photAcceptanceCor1D_Q_pT_bkwOnly_num->Write();
	h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_num->Write();
	h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_num->Write();
	h_photAcceptanceCor1D_Q_pT_midCMS_num->Write();
	h_photAcceptanceCor1D_Q_pT_fwdCMS_num->Write();
	h_photAcceptanceCor1D_Q_pT_bkwCMS_num->Write();

	h_photAcceptanceCor1D_Q_pT_all_rat->Write();
	h_photAcceptanceCor1D_Q_y_rat->Write();
	h_photAcceptanceCor1D_Q_nTrk_rat->Write();
	h_photAcceptanceCor1D_Q_nTrk_all_rat->Write();
	h_photAcceptanceCor1D_Q_pT_mid_rat->Write();
	h_photAcceptanceCor1D_Q_pT_fwd_rat->Write();
	h_photAcceptanceCor1D_Q_pT_fwdOnly_rat->Write();
	h_photAcceptanceCor1D_Q_pT_bkwOnly_rat->Write();
	h_photAcceptanceCor1D_Q_pT_fwdOnlyWide_rat->Write();
	h_photAcceptanceCor1D_Q_pT_bkwOnlyWide_rat->Write();
	h_photAcceptanceCor1D_Q_pT_midCMS_rat->Write();
	h_photAcceptanceCor1D_Q_pT_fwdCMS_rat->Write();
	h_photAcceptanceCor1D_Q_pT_bkwCMS_rat->Write();



	can1->Write();
	can2->Write();
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

	h_photDistribution3D_JpsipT->Write();
	h_photDistribution3D_Jpsiy->Write();

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
	h_JpsiEfficiency1D_den->Write();
	h_JpsiEfficiency1D_num->Write();
	h_JpsiEfficiency1D_rat->Write();
	h_chiEfficiency1D_den->Write();
	h_chiEfficiency1D_num->Write();
	h_chiEfficiency1D_rat->Write();
	//h_chiEfficiency1D_Q_ratRel->Write();

	h_photEfficiency1D_y_den->Write();
	h_photEfficiency1D_y_denJpsiVar->Write();
	h_photEfficiency1D_y_num->Write();
	h_photEfficiency1D_y_rat->Write();
	h_JpsiEfficiency1D_y_den->Write();
	h_JpsiEfficiency1D_y_num->Write();
	h_JpsiEfficiency1D_y_rat->Write();
	h_chiEfficiency1D_y_den->Write();
	h_chiEfficiency1D_y_num->Write();
	h_chiEfficiency1D_y_rat->Write();

	h_photEfficiency1D_nTrack_den->Write();
	h_photEfficiency1D_nTrack_num->Write();
	h_photEfficiency1D_nTrack_rat->Write();

	h_JpsiEfficiency1D_nTrack_den->Write();
	h_JpsiEfficiency1D_nTrack_num->Write();
	h_JpsiEfficiency1D_nTrack_rat->Write();
	h_chiEfficiency1D_nTrack_den->Write();
	h_chiEfficiency1D_nTrack_num->Write();
	h_chiEfficiency1D_nTrack_rat->Write();







	h_photEfficiency_y_relative_num_noCut->Write();
	h_photEfficiency_y_relative_num_noCutJpsiVar->Write();
	h_photEfficiency_y_relative_num_fullCut->Write();
	h_photEfficiency_y_relative_num_fullCutJpsiVar->Write();
	h_photEfficiency_y_relative_rat_noCut->Write();
	h_photEfficiency_y_relative_rat_noCutJpsiVar->Write();
	h_photEfficiency_y_relative_rat_fullCut->Write();
	h_photEfficiency_y_relative_rat_fullCutJpsiVar->Write();



	h_photEfficiency_y_relative_den->Write();
	h_photEfficiency_y_relative_num_loose->Write();
	h_photEfficiency_y_relative_num_medium->Write();
	h_photEfficiency_y_relative_num_tight->Write();
	h_photEfficiency_y_relative_rat_loose->Write();
	h_photEfficiency_y_relative_rat_medium->Write();
	h_photEfficiency_y_relative_rat_tight->Write();


	//weights
	h_JpsiAccEff_den->Write();
	h_JpsiAccEff_num->Write();
	h_JpsiAccEff_den_nTrk->Write();
	h_JpsiAccEff_num_nTrk->Write();
	h_chiAccEff_den->Write();
	h_chiAccEff_num->Write();
	h_chiAccEff_den_nTrk->Write();
	h_chiAccEff_num_nTrk->Write();

	h_JpsiAccEff_rat->Write();
	h_JpsiAccEff_rat_nTrk->Write();
	h_chiAccEff_rat->Write();
	h_chiAccEff_rat_nTrk->Write();


	fout->Close();
	fMC->Close();
	cout << "END" << endl;
	return 0;
}
