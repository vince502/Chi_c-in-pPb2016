
// Macro to analyze the chi_c. Starting with the event tree


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
#include "TAxis.h"
#include "TLegend.h"
#include "TF1.h"
#include "TClonesArray.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooTFnBinding.h"
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"


//#include "../ChiTreeInit.C"


using namespace std;
using namespace RooFit;

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
const int nMassBins = 160;
const double mass_window_l = 3.25;
const double mass_window_h = 4.05;
const double mass_windowFit_l = 3.25;
const double mass_windowFit_h = 4.05;
const string mass_windowFit = "rvmass>3.25 && rvmass<4.05";

const int nMassBinsJpsi = 250;
const double mass_windowJpsi_l = 2.0;
const double mass_windowJpsi_h = 4.5;
const double mass_windowFitJpsi_l = 2.0;
const double mass_windowFitJpsi_h = 4.5;
const string mass_windowFitJpsi = "rvmassJpsi>2.0 && rvmassJpsi<4.5";

double bins_pT[] = { 6,8,10,18,30 };
int  nbins_pT = sizeof(bins_pT) / sizeof(double) - 1;






int runNumber, eventNumber;
TLorentzVector* LVchic, *LVdimuon, *LVconv, *LVmuon1, *LVmuon2;
TLorentzVector* LVchic_rot, *LVchic_rotGamma, *LVdimuon_rot, *LVconv_rot, *LVmuon1_rot, *LVmuon2_rot;
TLorentzVector LVaux;

//other
int ntracks_inEvent;
double hfTowerSum_inEvent;
int Trig_Event_HLTDoubleMuOpen;
int nPrimVertices;

//
std::vector <double>* pvtx_z = 0;
std::vector <double>* pvtx_zError = 0;
std::vector <double>* pvtx_x = 0;
std::vector <double>* pvtx_y = 0;
std::vector <double>* pvtx_nTracks = 0;
std::vector <bool>* pvtx_isFake = 0;


//muon info
TClonesArray*  muon_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <bool>* muonIsHLTDoubleMuOpen = 0;
std::vector <bool>* muonIsGlobal = 0;
std::vector <bool>* muonIsTracker = 0;
std::vector <bool>* muonIsPF = 0;
std::vector <bool>* muonIsSoft = 0;
std::vector <bool>* muonIsTight = 0;
std::vector <int>* muonTrackerLayersWithMeasurement = 0;
std::vector <double>* muon_eta = 0;
std::vector <double>* muon_pt = 0;

//dimuon
TClonesArray*  dimuon_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <double>* dimuon_eta = 0;
std::vector <double>* dimuon_pt = 0;
std::vector <double>* dimuon_charge = 0;
std::vector <int>*  dimuon_pvtx_index = 0;
std::vector <double>* dimuon_dz_dimuonvtx_pvtx = 0;
std::vector <double>* dimuon_vtxProb = 0;
std::vector <int>* dimuon_muon1_position = 0; //stores position of first muon in muon collection (no specific order)
std::vector <int>* dimuon_muon2_position = 0; //stores position of second muon in muon collection (no specific order)
std::vector <double>* dimuon_ctpv = 0;
std::vector <double>* dimuon_ctpvError = 0;

// Chi
TClonesArray*  chi_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <double>* chi_eta = 0;
std::vector <double>* chi_pt = 0;
std::vector <int>* chi_daughterJpsi_position = 0; //stores position of daughter Jpsi in dimuon collection
std::vector <int>* chi_daughterConv_position = 0; //stores position of daughter photon (conversion)
std::vector <double>* chi_dzPhotToDimuonVtx = 0; //z distance of photon to dimuon vertex when dxy is minimal
std::vector <double>* chi_dxyPhotToDimuonVtx = 0; //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions


//// Conversions  /////

TClonesArray*  conv_p4 = new TClonesArray("TLorentzVector", 100);
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

//MC
std::vector <int>* gen_Jpsi_matchPosition = 0;
std::vector <int>* gen_conv_matchPosition = 0;
std::vector <int>* convGen_motherCode;


bool MuonAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 1.1 && pt < 3.3) return false;  
	if (fabs(eta) >= 1.1 && fabs(eta) < 2.1 && pt < 5.5 - 2*fabs(eta)) return false;
	if (fabs(eta) >= 2.1 && pt < 1.3) return false;
	return true;
}
bool MuonAcceptanceTight(double eta, double pt)
{
	return true;
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 2.4 && pt < 3.3) return false;  
	return true;
}

bool PhotAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.5) return false; //2.5
	if (pt < 0.2) return false; // 0.2
	return true;
}

bool DimuonSelectionPass(int dimuonPos)  //uses variables loaded in main function
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	if (dimuon_pt->at(dimuonPos) < 6) return false;
	return true;
}

bool DimuonSelectionPassTight(int dimuonPos, double rap)  //uses variables loaded in main function
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	//if (fabs(rap) > 1.0) return false;
	//if (dimuon_pt->at(dimuonPos) > 9.0) return false;
	//if (dimuon_pt->at(dimuonPos) < 6) return false;
	return true;
}

bool MuonSelectionPass(int muonPos)  //uses variables loaded in main function
{
	if (muonIsSoft->at(muonPos) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(muonPos) != 1) return false;
	return true;
}


bool PhotSelectionPass(int photPos)  //uses variables loaded in main function
{
	
	//if (convQuality_isHighPurity->at(photPos) != 1) return false;
	//if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	//if (conv_sigmaTkVtx1->at(photPos) > 10) return false;
	//if (conv_sigmaTkVtx2->at(photPos) > 10) return false;
	//if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	//if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (conv_vertexChi2Prob->at(photPos) <= 0.01) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	//if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	//if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_track1Chi2->at(photPos) >= 10) return false;
	//if (conv_track2Chi2->at(photPos) >= 10) return false;
	//if (conv_minDistanceOfApproach->at(photPos) <= -0.25) return false;
	//if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}

int ImportDatasetFromTree(RooWorkspace& Ws, string label) //to be done later
{
	return 0;
};


bool CreateModelPdf(RooWorkspace& Ws, string pdfName)
{

	//RooRealVar erf_offset("erf_offset", "erf_offset", 3.3, 3.1, 3.5);
	//RooRealVar erf_sigma("erf_sigma", "erf_sigma", 0.1, 0.04, 0.2);
	//RooGenericPdf *ErfPdf = new RooGenericPdf("ErfPdf", "Error Function with offset and mean", "(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5", RooArgList(*(Ws.var("rvmass")), erf_offset, erf_sigma));
	//Ws.import(*ErfPdf);

	//RooRealVar exp_lambda("exp_lambda", "exp_lambda", -0.1, -0.2, 0.0);
	//RooRealVar exp_const("exp_const", "exp_const", 0.8, 0.01, 1.0);
	//Ws.import(exp_lambda);
	//Ws.import(exp_const);
	//Ws.factory("EXPR::expConst('exp(exp_lambda * rvmass) + exp_const', rvmass, exp_lambda, exp_const)");
	//
	//= new RooGenericPdf("ErfPdf", "Error Function with offset and mean", "(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5", RooArgList(*(Ws.var("rvmass")), erf_offset, erf_sigma));
	//Ws.import(*ErfPdf);



	//if (pdfName.find("nominalPdf") != string::npos)
	if (pdfName.compare("nominalPdf") == 0)
	{
		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		Ws.factory("CBShape::chic2(rvmass, mean2[3.5562, 3.54, 3.60], sigma1, alpha, n)");
		Ws.factory("SUM::signal(c2toc1[0.4,0.1,1.]*chic2, chic1)");
		//Ws.factory("Exponential::expbkg(rvmass,e_1[-0.2,-0.5,0])");
		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("RooChebychev::background(rvmass, a1[0,-10,10])");
		Ws.var("alpha")->setConstant(true);
		Ws.var("n")->setConstant(true);
		//Ws.var("mean1")->setConstant(true);
		//Ws.var("mean2")->setConstant(true);
		//Ws.var("erf_offset")->setConstant(true);
		//Ws.var("erf_sigma")->setConstant(true);
		//Ws.factory("")
		Ws.factory("SUM::nominalPdf(nsig[50,0,100000]*signal, nbkg[2000,0,100000]*background)");
	}
	else if (pdfName.compare("nominalPdfJpsi") == 0)
	{
		Ws.factory("CBShape::Jpsi(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08], alphaJpsi[1.85, 0.1, 50], nJpsi[1.7, 0.2, 50])");
		Ws.factory("CBShape::psi2(rvmassJpsi, mean2Jpsi[3.686, 3.50, 3.75], sigma1Jpsi, alphaJpsi, nJpsi)");
		Ws.factory("SUM::signalJpsi(ratioPsi[0.1,0.01,2.0]*psi2, Jpsi)");
		//Ws.factory("Exponential::expbkg(rvmass,e_1[-0.2,-0.5,0])");
		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("RooChebychev::backgroundJpsi(rvmassJpsi, a1Jpsi[0.00,-10,10])");
		Ws.var("alphaJpsi")->setConstant(true);
		Ws.var("nJpsi")->setConstant(true);
		//Ws.var("mean1")->setConstant(true);
		//Ws.var("mean2")->setConstant(true);
		//Ws.factory("")
		Ws.factory("SUM::nominalPdfJpsi(nsigJpsi[500,0,1000000]*signalJpsi, nbkgJpsi[2000,0,1000000]*backgroundJpsi)");
	}


	
	//Ws.pdf(pdfName.c_str())->setNormRange(mass_windowFit_l, mass_windowFit_h);
	//RooRealVar nsig("nsig", "nsig", 500, 0., 100000.);
	//RooRealVar nbkg("nbkg", "nbkg", 500, 0., 100000.);

	//RooAddPdf* modelPdf = new RooAddPdf(pdfName.c_str(), pdfName.c_str(), 
	//	RooArgList(signal),// background),
	//	RooArgList(nsig)//, nbkg)
	//	);
	//Ws.import(*modelPdf);
	return true;
}


//cbGausPlusPol1 = cms.vstring(
//	"CBShape::signal1(mass, mean[3.08,3.00,3.3], sigma1[0.03, 0.01, 0.10], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])",
//	"RooFormulaVar::sigma2('@0*@1',{fracS[1.8,1.2,2.4],sigma1})",
//	"Gaussian::signal2(mass, mean, sigma2)",
//	"SUM::signal(frac[0.8,0.5,1.]*signal1,signal2)",
//	"Chebychev::backgroundPass(mass, {cPass[0.,-1.1,1.1]})",
//	"Chebychev::backgroundFail(mass, {cFail[0.,-1.1,1.1]})",
//	"efficiency[0.9,0,1]",
//	"signalFractionInPassing[0.9]"
//),

///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////

//void Analyze_Chic(bool flagGenerateRds = true, const char* fileIn = "/afs/cern.ch/user/o/okukral/Chic_pPb/CMSSW_8_0_30/src/HeavyIonsAnalysis/ChiAnalysis/test/Chi_c_pPb8TeV_MC5_v2.root", const char* fileOut = "Chi_c_output_MC_test.root", const char* fileRds = "rds_saveMC.root", bool doMC = false)
//void Analyze_Chic(bool flagGenerateRds = true, const char* fileIn = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-PbpRW3.root", const char* fileOut = "Chi_c_output_RW3_test.root", const char* fileRds = "rds_save.root", bool doMC = false)
void Analyze_Chic(bool flagGenerateRds = false, const char* fileIn = "Chi_c_pPb8TeV-PbpRW3.root", const char* fileOut = "Chi_c_output_RW3_test3.root", const char* fileRds = "rds_save.root", bool doMC = false, const char* fileCorrection = "Chi_c_distributionsMC7.root")
//void Analyze_Chic(bool flagGenerateRds = true, const char* fileIn = "Chi_c_pPb8TeV-bothDirRW3.root", const char* fileOut = "Chi_c_output_RW3_testFull.root", const char* fileRds = "rds_save.root", bool doMC = false)
//void Analyze_Chic(const char* fileIn = "Chi_c_pPb8TeV_MCv5.root", const char* fileOut = "Chi_c_output_MCv5_test1.root", const bool doMC = false)
{
	gStyle->SetOptStat(1111);
	//gStyle->SetOptStat(0);


	TH1D* hSignal_pT = new TH1D("hSignal_pT", "", nbins_pT, bins_pT);
	TH1D* hSignalJpsi_pT = new TH1D("hSignalJpsi_pT", "", nbins_pT, bins_pT);
	TH1D* hSignalRatio_pT = new TH1D("hSignalRatio_pT", "", nbins_pT, bins_pT);

	TH1D* hSignal = new TH1D("hSignal", "", 200, 3, 5);
	TH1D* hSignal_SS = new TH1D("hSignal_SS", "", 200, 3, 5);
	TH1D* hSignal_SB = new TH1D("hSignal_SB", "", 200, 3, 5);
	TH1D* hSignal_rot = new TH1D("hSignal_rot", "", 200, 3, 5);
	TH1D* hSignal_rotGamma = new TH1D("hSignal_rotGamma", "", 200, 3, 5);
	TH1D* hSignal2 = new TH1D("hSignal2", "", 120, 3.2, 3.8);
	TH1D* hSignal3 = new TH1D("hSignal3", "", 120, 3.2, 3.8);

	////Quick pT binning
	//TH1D* hSignalpt1 = new TH1D("hSignalpt1", "", 100, 3, 4);
	//TH1D* hSignalpt2 = new TH1D("hSignalpt2", "", 100, 3, 4);
	//TH1D* hSignalpt3 = new TH1D("hSignalpt3", "", 100, 3, 4);
	//TH1D* hSignalpt4 = new TH1D("hSignalpt4", "", 100, 3, 4);
	//TH1D* hSignalpt5 = new TH1D("hSignalpt5", "", 100, 3, 4);
	//TH1D* hSignalpt6 = new TH1D("hSignalpt6", "", 100, 3, 4);
	//TH1D* hSignalpt7 = new TH1D("hSignalpt7", "", 100, 3, 4);

	//TH2D* hSignal_m_pt = new TH2D("hSignal_m_pt", "", 400, 0, 40, 200, 3, 5);
	//TH2D* hSignal_m_rap = new TH2D("hSignal_m_y", "", 60, -3, 3, 200, 3, 5);
	//TH2D* hSignal_m_phi = new TH2D("hSignal_m_phi", "", 70, -3.5, 3.5, 200, 3, 5);
	//TH2D* hSignal_m_nTr = new TH2D("hSignal_m_nTr", "", 300, 0, 600, 200, 3, 5);
	//TH2D* hSignal_m_Et = new TH2D("hSignal_m_Et", "", 400, 0, 800, 200, 3, 5);


	TH1D* hdimuon_M = new TH1D("hdimuon_M", "", 300, 2, 5);
	TH1D* hchic_M = new TH1D("hchic_M", "", 100, 2, 5);
	TH1D* hphoton_M = new TH1D("hphoton_M", "", 1000, 0, 1);

	TH1D* hdimuon_M_SS = new TH1D("hdimuon_M_SS", "", 300, 2, 5);
	TH1D* hchic_M_SS = new TH1D("hchic_M_SS", "", 100, 2, 5);
	TH1D* hdimuon_M_SB = new TH1D("hdimuon_M_SB", "", 300, 2, 5);
	TH1D* hchic_M_SB = new TH1D("hchic_M_SB", "", 100, 2, 5);
	TH1D* hdimuon_M_rot = new TH1D("hdimuon_M_rot", "", 300, 2, 5);
	TH1D* hchic_M_rot = new TH1D("hchic_M_rot", "", 100, 2, 5);
	TH1D* hchic_M_rotGamma = new TH1D("hchic_M_rotGamma", "", 100, 2, 5);

	//	TH1D* hdimuon_pt = new TH1D("hdimuon_pt", "", 1000, 0, 100);
	//	TH1D* hphoton_pt = new TH1D("hphoton_pt", "", 1000, 0, 100);

	//	TH2D* hchicCand_pt_rap = new TH2D("hchicCand_pt_rap", "", 60, -3, 3, 200, 0, 20);
	//	TH2D* hdimuon_pt_rap = new TH2D("hdimuon_pt_rap", "", 60, -3, 3, 200, 0, 20);
	//	TH2D* hphoton_pt_rap = new TH2D("hphoton_pt_rap", "", 60, -3, 3, 100, 0, 10);

		//TH2D* hdimuonAfterCuts_pt_rap = new TH2D("hdimuonAfterCuts_pt_rap", "", 60, -3, 3, 200, 0, 20);
		//TH2D* hmuonPAfterCuts_pt_ps = new TH2D("hmuonPAfterCuts_pt_ps", "", 60, -3, 3, 200, 0, 20);
		//TH2D* hmuonNAfterCuts_pt_ps = new TH2D("hmuonNAfterCuts_pt_ps", "", 60, -3, 3, 200, 0, 20);
		//TH2D* hphotonAfterCuts_pt_rap = new TH2D("hphotonAfterCuts_pt_rap", "", 60, -3, 3, 100, 0, 10);

		//TH2D* hchi_restricted_pt_rap = new TH2D("hchi_restricted_pt_rap", "", 60, -3, 3, 200, 0, 20);
		//TH1D* hchi_restrictedM_pt = new TH1D("hchi_restrictedM_pt", "", 1000, 0, 100);
		//TH1D* hchi_restrictedM_dimuonpt = new TH1D("hchi_restrictedM_dimuonpt", "", 1000, 0, 100);
		//TH1D* hchi_restrictedM_mPpt = new TH1D("hchi_restrictedM_mPpt", "", 1000, 0, 100);
		//TH1D* hchi_restrictedM_mNpt = new TH1D("hchi_restrictedM_mNpt", "", 1000, 0, 100);
		//TH1D* hchi_restrictedM_Photpt = new TH1D("hchi_restrictedM_Photpt", "", 1000, 0, 100);


		// Various
	TH1D* hntracks_inEvent = new TH1D("hntracks_inEvent", "", 400, 0, 400);

	TFile* f1 = new TFile(fileIn, "READ");

	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");
	if (!event_tree) {
		cout << "Problem with event Tree";
		//break;
	}

	//LoadChiBranches(event_tree, false);
	event_tree->SetBranchAddress("nPrimVertices", &nPrimVertices);
	event_tree->SetBranchAddress("ntracks_inEvent", &ntracks_inEvent);

	event_tree->SetBranchAddress("pvtx_z", &pvtx_z);
	event_tree->SetBranchAddress("pvtx_zError", &pvtx_zError);
	event_tree->SetBranchAddress("pvtx_x", &pvtx_x);
	event_tree->SetBranchAddress("pvtx_y", &pvtx_y);
	event_tree->SetBranchAddress("pvtx_nTracks", &pvtx_nTracks);
	event_tree->SetBranchAddress("pvtx_isFake", &pvtx_isFake);

	event_tree->SetBranchAddress("muon_p4", &muon_p4);
	event_tree->SetBranchAddress("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
	event_tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal);
	event_tree->SetBranchAddress("muonIsTracker", &muonIsTracker);
	event_tree->SetBranchAddress("muonIsPF", &muonIsPF);
	event_tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
	event_tree->SetBranchAddress("muonIsTight", &muonIsTight);
	event_tree->SetBranchAddress("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);
	event_tree->SetBranchAddress("muon_eta", &muon_eta);
	event_tree->SetBranchAddress("muon_pt", &muon_pt);

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

	event_tree->SetBranchAddress("chi_p4", &chi_p4);
	event_tree->SetBranchAddress("chi_eta", &chi_eta);
	event_tree->SetBranchAddress("chi_pt", &chi_pt);
	event_tree->SetBranchAddress("chi_daughterJpsi_position", &chi_daughterJpsi_position);
	event_tree->SetBranchAddress("chi_daughterConv_position", &chi_daughterConv_position);
	event_tree->SetBranchAddress("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx);
	event_tree->SetBranchAddress("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx);

	//Conversions
	event_tree->SetBranchAddress("conv_p4", &conv_p4);
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

	// Gen stuff
	if (doMC) {
		event_tree->SetBranchAddress("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
		event_tree->SetBranchAddress("gen_conv_matchPosition", &gen_conv_matchPosition);
		event_tree->SetBranchAddress("convGen_motherCode", &convGen_motherCode);
	}

	TRandom3* rand = new TRandom3();

	////////////////////////////////////////  
///////   R  o  o  f  i  t   //////////
/////////////////////////////////////////

	RooWorkspace myWs("myWs", "main workspace");
	//ImportDatasetFromTree(myWs, "nominal");
	//RooRealVar* rvmass = new RooRealVar("rvmass", "#mu#mu mass", 3.0, 5.0, "GeV/c^{2}");
	RooRealVar* rvmass = new RooRealVar("rvmass", "#mu#mu#gamma - #mu#mu + 3.097", mass_window_l, mass_window_h, "GeV/c^{2}");
	//RooRealVar* ctau = new RooRealVar("ctau", "c_{#tau}", -100.0, 100.0, "mm");
	//RooRealVar* ctauErr = new RooRealVar("ctauErr", "#sigma_{c#tau}", -100.0, 100.0, "mm");
	RooRealVar* rvpt = new RooRealVar("rvpt", "#mu#mu p_{T}", 0.0, 50.0, "GeV/c");
	RooRealVar* rvrap = new RooRealVar("rvrap", "#mu#mu y", -2.4, 2.4, "");
	RooRealVar* rvntrack = new RooRealVar("rvntrack", "ntrack", 0.0, 200.0, "");
	RooRealVar* rvweight = new RooRealVar("rvweight", "weight", 0.0, 10000.0, "");
	RooArgSet*  cols = new RooArgSet(*rvmass, *rvpt, *rvrap, *rvntrack, *rvweight);
	//RooArgSet*  cols = new RooArgSet(*rvmass);
	RooDataSet* rdsNominal = new RooDataSet("rdsNominal", "rdsNominal", *cols, WeightVar(*rvweight), StoreAsymError(*rvmass));
	//RooDataSet* rdsNominal = new RooDataSet("nominal", "nominal", *cols);

	RooRealVar* rvmassJpsi = new RooRealVar("rvmassJpsi", "Inv mass #mu#mu", mass_windowJpsi_l, mass_windowJpsi_h, "GeV/c^{2}");
	RooRealVar* rvptJpsi = new RooRealVar("rvptJpsi", "#mu#mu p_{T}", 0.0, 50.0, "GeV/c");
	RooRealVar* rvrapJpsi = new RooRealVar("rvrapJpsi", "#mu#mu y", -2.4, 2.4, "");
	RooRealVar* rvntrackJpsi = new RooRealVar("rvntrackJpsi", "ntrack", 0.0, 200.0, "");
	RooRealVar* rvweightJpsi = new RooRealVar("rvweightJpsi", "weight", 0.0, 10000.0, "");
	RooArgSet*  colsJpsi = new RooArgSet(*rvmassJpsi, *rvptJpsi, *rvrapJpsi, *rvntrackJpsi, *rvweightJpsi);
	RooDataSet* rdsNominalJpsi = new RooDataSet("rdsNominalJpsi", "rdsNominalJpsi", *colsJpsi, WeightVar(*rvweightJpsi), StoreAsymError(*rvmassJpsi));


	if (flagGenerateRds == true) {

		long nchicCounter = 0, nchicCounterPass = 0, muon1ptCounter = 0, muon2ptCounter = 0;
		bool passDimSel = false;
		bool passDimSelTight = false;
		Long64_t nentries = event_tree->GetEntries();
		cout << nentries << endl;
		if (nentries > 1000) { nentries = 1000; }
		for (Long64_t i = 0; i < nentries; i++)
		{

			event_tree->GetEntry(i);
			if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }




			////////////////////////
			//// V A R I O U S  ////
			///////////////////////
			hntracks_inEvent->Fill(ntracks_inEvent);


			///////////////////////
			/////   C H I C ///////
			/////////////////////

			for (int iChi = 0; iChi < chi_p4->GetEntriesFast(); iChi++)
			{
				++nchicCounter;
				cout << "in" << endl;
				// check Acceptance and Cuts
				int dimuonPos = chi_daughterJpsi_position->at(iChi);
				//if (DimuonSelectionPass(dimuonPos) == false) continue;
				passDimSel = DimuonSelectionPass(dimuonPos);
				passDimSelTight = DimuonSelectionPassTight(dimuonPos, ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity());
				cout << "in" << endl;
				// muon cuts
				int muon1Pos = dimuon_muon1_position->at(dimuonPos);
				int muon2Pos = dimuon_muon2_position->at(dimuonPos);
				cout << "in" << endl;
				if (MuonAcceptance(muon_eta->at(muon1Pos), muon_pt->at(muon1Pos)) == false) continue;
				if (MuonAcceptance(muon_eta->at(muon2Pos), muon_pt->at(muon2Pos)) == false) continue;
				cout << "in" << endl;
				if (MuonSelectionPass(muon1Pos) == false) continue;
				if (MuonSelectionPass(muon2Pos) == false) continue;
				cout << "in" << endl;
				// photon
				int convPos = chi_daughterConv_position->at(iChi);
				if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) continue;
				if (PhotSelectionPass(convPos) == false) continue;
				if (passDimSel == true) { ++nchicCounterPass; } // SelectionsPassed
				// Get Lorentz V
				LVchic = (TLorentzVector*)chi_p4->At(iChi);
				LVdimuon = (TLorentzVector*)dimuon_p4->At(dimuonPos);
				LVconv = (TLorentzVector*)conv_p4->At(convPos);


				//////// Rotational bkg   /////
				LVmuon1_rot = (TLorentzVector*)muon_p4->At(muon1Pos);
				LVmuon2_rot = (TLorentzVector*)muon_p4->At(muon2Pos);
				LVconv_rot = new TLorentzVector(*(TLorentzVector*)conv_p4->At(convPos));
				LVconv_rot->RotateZ(TMath::Pi());
				//if (muon_pt->at(muon1Pos) > muon_pt->at(muon2Pos)) { muon1ptCounter++; } //crosscheck that they have ordering
				//else muon2ptCounter++;
				if (rand->Uniform(0, 1) < 0.5) { LVmuon1_rot->RotateZ(TMath::Pi()); }
				else { LVmuon2_rot->RotateZ(TMath::Pi()); }
				LVdimuon_rot = new TLorentzVector(*LVmuon1_rot + *LVmuon2_rot);
				//cout << LVdimuon_rot->Pt() << endl;
				LVchic_rot = new TLorentzVector(*LVdimuon_rot + *LVconv);
				LVchic_rotGamma = new TLorentzVector(*LVdimuon + *LVconv_rot);


				// Obtain yields

				// fill dimuon stuff

				double dimuonM = LVdimuon->M();
				if (passDimSel == true) {
					hdimuon_M->Fill(dimuonM);
					hdimuon_M_rot->Fill(LVdimuon_rot->M());
				}
				if (dimuon_charge->at(dimuonPos) != 0) {
					hdimuon_M_SS->Fill(dimuonM);
				}


				// Rotational bkg
				if (LVdimuon_rot->M() > 2.95 && LVdimuon_rot->M() < 3.2) {
					if (passDimSel == true) {
						hchic_M_rot->Fill(LVchic_rot->M());
						hSignal_rot->Fill(LVchic_rot->M() - LVdimuon_rot->M() + 3.097);
					}
				}
				// Side-band
				if (((dimuonM > 2.5 && dimuonM < 2.75) || (dimuonM > 3.4 && dimuonM < 3.7)) && passDimSel == true) {
					hchic_M_SB->Fill(LVchic->M());
					hSignal_SB->Fill(LVchic->M() - LVdimuon->M() + 3.097);
				}


				if (dimuonM<2.95 || dimuonM>3.2) continue; //require narrow dimuon mass

				if (passDimSel == true) {
					hchic_M->Fill(LVchic->M()); // just raw M
					hchic_M_rotGamma->Fill(LVchic_rotGamma->M());
				}
				if (dimuon_charge->at(dimuonPos) != 0) {
					hchic_M_SS->Fill(LVchic->M());
				}


				// Assume J/psi mass

				double Mdiff = LVchic->M() - LVdimuon->M() + 3.097;
				if (passDimSel == true) {
					hSignal->Fill(Mdiff);
					hSignal_rotGamma->Fill(LVchic_rotGamma->M() - LVdimuon->M() + 3.097);
					//cout << nchicCounterPass << endl;

					//roofit:
					if (Mdiff > mass_window_l && Mdiff < mass_window_h) {
						rvmass->setVal(Mdiff);
						rvpt->setVal(chi_pt->at(iChi));
						rvrap->setVal(LVchic->Rapidity());
						rvntrack->setVal(ntracks_inEvent); //for now using total
						rvweight->setVal(1); // no weights for now

						rdsNominal->add(*cols);

						
					}

				}
				if (dimuon_charge->at(dimuonPos) != 0) { hSignal_SS->Fill(Mdiff); }

				cout << "i2n" << endl;

			} // end of chic loop

			for (int iJpsi = 0; iJpsi < dimuon_p4->GetEntriesFast(); iJpsi++) // Jpsi loop
			{

				passDimSel = DimuonSelectionPass(iJpsi);

				// muon cuts
				int muon1Pos = dimuon_muon1_position->at(iJpsi);
				int muon2Pos = dimuon_muon2_position->at(iJpsi);
				if (MuonAcceptance(muon_eta->at(muon1Pos), muon_pt->at(muon1Pos)) == false) continue;
				if (MuonAcceptance(muon_eta->at(muon2Pos), muon_pt->at(muon2Pos)) == false) continue;
				if (MuonSelectionPass(muon1Pos) == false) continue;
				if (MuonSelectionPass(muon2Pos) == false) continue;

				// Get Lorentz V
				LVdimuon = (TLorentzVector*)dimuon_p4->At(iJpsi);

				// Obtain yields

				// fill dimuon stuff

				double dimuonM = LVdimuon->M();

				if (passDimSel == true) {

					//roofit:
					if (dimuonM > mass_windowJpsi_l && dimuonM < mass_windowJpsi_h) {
						rvmassJpsi->setVal(dimuonM);
						rvptJpsi->setVal(dimuon_pt->at(iJpsi));
						rvrapJpsi->setVal(LVdimuon->Rapidity());
						rvntrackJpsi->setVal(ntracks_inEvent); //for now using total
						rvweightJpsi->setVal(1); // no weights for now

						rdsNominalJpsi->add(*colsJpsi);


					}

				}

			} //end of Jpsi loop





		} //end of event loop

		cout << "Processed " << nchicCounter << " chic candidates, " << nchicCounterPass << " passed the selections, that is " << ((double)nchicCounterPass / (double)nchicCounter * 100) << " percent" << endl;

		// Save rds
		TFile* DBFile = new TFile(fileRds, "RECREATE");
		DBFile->cd();
		rdsNominal->Write();
		rdsNominalJpsi->Write();
		DBFile->Write(); DBFile->Close(); delete DBFile;
		f1->cd();
		myWs.import(*rdsNominal);

	}
	else //load the rds instead
	{
		cout << "[INFO] Loading RooDataSet from " << fileRds << endl;
		TFile* DBFile = new TFile(fileRds, "OPEN");
		DBFile->cd();
		cout << "content " << endl;
		DBFile->Print();
		DBFile->ls();
		cout << endl;
		//RooDataSet* d = (RooDataSet*)f.FindObject("d")
		rdsNominal = (RooDataSet*)DBFile->Get("rdsNominal");
		//DBFile->GetObject("rdsNomial", rdsNominal2);
		cout << "nEntries: " << rdsNominal->numEntries()<<endl;
		myWs.import(*rdsNominal);
		rdsNominalJpsi = (RooDataSet*)DBFile->Get("rdsNominalJpsi");
		cout << "nEntriesJpsi: " << rdsNominalJpsi->numEntries() << endl;
		myWs.import(*rdsNominalJpsi);

		DBFile->Close(); delete DBFile;
	}

	//////////////////////////////////////////  
	/////////   R  o  o  f  i  t   //////////
	///////////////////////////////////////////


	cout << "nEntries: " << rdsNominal->numEntries()<<endl;
	rdsNominal = (RooDataSet*)rdsNominal->reduce(mass_windowFit.c_str());
	//RooDataSet* dataOS = (RooDataSet*)myWs.data("rdsNominal")->reduce("pt>20");
	//RooDataSet* dataOS2 = (RooDataSet*)myWs.data("rdsNominal")->reduce("pt<5");

	RooPlot *massframe = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins );
	massframe->SetTitle("mass");
	rdsNominal->plotOn(massframe);
	//dataOS->plotOn(massframe);
	//dataOS2->plotOn(massframe);
	
	string myPdfName = "nominalPdf";
	CreateModelPdf(myWs, myPdfName);

	
	RooFitResult* fitResult = myWs.pdf(myPdfName.c_str())->fitTo(*myWs.data("rdsNominal"), Extended(true), SumW2Error(true), Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true));
	//fitResult->Print("v");
	myWs.import(*fitResult, TString::Format("fitResult_%s", myPdfName.c_str()));


	myWs.pdf(myPdfName.c_str())->plotOn(massframe);
	myWs.pdf(myPdfName.c_str())->paramOn(massframe, Layout(0.55));
	myWs.pdf(myPdfName.c_str())->plotOn(massframe, Components("background"), LineStyle(kDashed));
	myWs.pdf(myPdfName.c_str())->plotOn(massframe, Components("chic1"), LineStyle(kDashed), LineColor(kRed));
	myWs.pdf(myPdfName.c_str())->plotOn(massframe, Components("chic2"), LineStyle(kDashed), LineColor(kGreen));
	cout << endl << endl << "HERE" << endl << endl;
	
	fitResult->floatParsFinal().Print("s");
	cout << endl << endl << endl;

	myWs.var("nsig")->Print();
	//fitResult->PrintArgs();

	TCanvas *cTest = new TCanvas("cTest", "cTest", 1000, 480);
	gPad->SetLeftMargin(0.15);
	massframe->GetYaxis()->SetTitleOffset(1.6);
	massframe->Draw();


	cTest->SaveAs("CanvasCTest_RW3.png");


	// pT fitting

	for (int i = 0; i < (sizeof(bins_pT) / sizeof(*bins_pT)-1); i++) {
		RooPlot *massframeBin = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
		massframeBin->SetTitle("mass");
		cout << bins_pT[i] << endl;
		TString TstrCut = TString::Format("rvpt > %f", bins_pT[i]) + " && " + TString::Format("rvpt < %f", bins_pT[i + 1]) + " && rvrap>-1 && rvrap <1";
		cout << TstrCut << endl;
		string strCut = TstrCut.Data();
		RooDataSet* rdsDataBin = (RooDataSet*)myWs.data("rdsNominal")->reduce(strCut.c_str());
		rdsDataBin->plotOn(massframeBin);

		//myWs.var("rvpt")->setMin(6);
		//myWs.var("rvpt")->setMax(8);

		//fitting
		//string myPdfNameBin = TString::Format("nominalPdf_%i",i);
		//string myPdfNameBin = "nominalPdf";
	//cout << endl << endl << endl << myPdfNameBin << endl << endl;
	//CreateModelPdf(myWs, myPdfNameBin);

	RooFitResult* fitResultBin = myWs.pdf(myPdfName.c_str())->fitTo(*rdsDataBin);// , Extended(true), SumW2Error(true), Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true));
	//fitResultBin->Print("v");
	//myWs.import(*fitResultBin, TString::Format("fitResult_%s", myPdfNameBin.c_str()));


	myWs.pdf(myPdfName.c_str())->plotOn(massframeBin);
	myWs.pdf(myPdfName.c_str())->paramOn(massframeBin, Layout(0.55));
	myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("background"), LineStyle(kDashed));
	myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic1"), LineStyle(kDashed), LineColor(kRed));
	myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic2"), LineStyle(kDashed), LineColor(kGreen));


	massframeBin->Draw();
	cTest->SaveAs(Form("CanvasCTest_RW3_pT_%i.png",i));

	//myWs.Delete("nominalPdf");
}


// J/psi fit

rdsNominalJpsi = (RooDataSet*)rdsNominalJpsi->reduce(mass_windowFitJpsi.c_str());
RooPlot *massframeJpsi = rvmassJpsi->frame(mass_windowFitJpsi_l, mass_windowFitJpsi_h, nMassBinsJpsi);
massframeJpsi->SetTitle("massJpsi");
rdsNominalJpsi->plotOn(massframeJpsi);

string myPdfNameJpsi = "nominalPdfJpsi";
CreateModelPdf(myWs, myPdfNameJpsi);

RooFitResult* fitResultJpsi = myWs.pdf(myPdfNameJpsi.c_str())->fitTo(*myWs.data("rdsNominalJpsi"), Extended(true), SumW2Error(true), Range(mass_windowFitJpsi_l, mass_windowFitJpsi_h), NumCPU(1), Save(true));
fitResultJpsi->Print("v");
myWs.import(*fitResultJpsi, TString::Format("fitResultJpsi_%s", myPdfNameJpsi.c_str()));


myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi);
myWs.pdf(myPdfNameJpsi.c_str())->paramOn(massframeJpsi, Layout(0.55));
myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("backgroundJpsi"), LineStyle(kDashed));
myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("Jpsi"), LineStyle(kDashed), LineColor(kRed));
myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("psi2"), LineStyle(kDashed), LineColor(kGreen));
cout << endl << endl << "HEREJpsi" << endl << endl;



TCanvas *cTestJpsi = new TCanvas("cTestJpsi", "cTestJpsi", 1000, 480);
gPad->SetLeftMargin(0.15);
massframeJpsi->GetYaxis()->SetTitleOffset(1.6);
massframeJpsi->Draw();


cTestJpsi->SaveAs("CanvasCTest_RW3_Jpsi.png");


// pT fitting

for (int i = 0; i < (sizeof(bins_pT) / sizeof(*bins_pT) - 1); i++) {
	RooPlot *massframeJpsiBin = rvmassJpsi->frame(mass_windowFitJpsi_l, mass_windowFitJpsi_h, nMassBinsJpsi);
	massframeJpsiBin->SetTitle("massJpsi");
	cout << bins_pT[i] << endl;
	TString TstrCut = TString::Format("rvptJpsi > %f", bins_pT[i]) + " && " + TString::Format("rvptJpsi < %f", bins_pT[i + 1]) + "&& rvrapJpsi>-1 && rvrapJpsi <1";
	cout << TstrCut << endl;
	string strCut = TstrCut.Data();
	RooDataSet* rdsDataJpsiBin = (RooDataSet*)myWs.data("rdsNominalJpsi")->reduce(strCut.c_str());
	rdsDataJpsiBin->plotOn(massframeJpsiBin);

	//myWs.var("rvpt")->setMin(6);
	//myWs.var("rvpt")->setMax(8);

	////fitting
	//string myPdfNameJpsiBin = TString::Format("nominalPdf_%i",i);
	//string myPdfNameJpsiBin = "nominalPdf";
	//cout << endl << endl << endl << myPdfNameJpsiBin << endl << endl;
	//CreateModelPdf(myWs, myPdfNameJpsiBin);

	RooFitResult* fitResultJpsiBin = myWs.pdf(myPdfNameJpsi.c_str())->fitTo(*rdsDataJpsiBin);// , Extended(true), SumW2Error(true), Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true));
	//fitResultJpsiBin->Print("v");
	//myWs.import(*fitResultJpsiBin, TString::Format("fitResultJpsi_%s", myPdfNameJpsiBin.c_str()));


	myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsiBin);
	myWs.pdf(myPdfNameJpsi.c_str())->paramOn(massframeJpsiBin, Layout(0.55));
	myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsiBin, Components("backgroundJpsi"), LineStyle(kDashed));
	myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsiBin, Components("Jpsi"), LineStyle(kDashed), LineColor(kRed));
	myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsiBin, Components("psi2"), LineStyle(kDashed), LineColor(kGreen));


	massframeJpsiBin->Draw();
	cTestJpsi->SaveAs(Form("CanvasCTest_RW3_Jpsi_pT_%i.png", i));

}



	hSignal_pT->SetBinContent(1, 2);
	hSignal_pT->SetBinError(1, 2);
	hSignalJpsi_pT->SetBinContent(1, 2);

	hSignal_pT->SetBinContent(2, 2);
	hSignal_pT->SetBinError(2, 2);
	hSignalJpsi_pT->SetBinContent(2, 2);

	hSignal_pT->SetBinContent(3, 2);
	hSignal_pT->SetBinError(3, 2);
	hSignalJpsi_pT->SetBinContent(3, 2);

	hSignal_pT->SetBinContent(4, 2);
	hSignal_pT->SetBinError(4, 2);
	hSignalJpsi_pT->SetBinContent(4, 2);

	hSignalRatio_pT->Divide(hSignal_pT, hSignalJpsi_pT, 1, 1);

	TFile* fCor = new TFile(fileCorrection, "READ");

	TH1D* hCor = (TH1D*)fCor->Get("h_chiEfficiency1D_Q_ratRel");
	for (int i = 1; i < (hCor->GetNbinsX() + 1); i++)
		{
		cout << hCor->GetBinContent(i) << endl;

		hSignalRatio_pT->SetBinContent(i, hSignalRatio_pT->GetBinContent(i)/hCor->GetBinContent(i));
		hSignalRatio_pT->SetBinError(i, hSignalRatio_pT->GetBinError(i)/hCor->GetBinContent(i));
	}

	TCanvas* can4 = new TCanvas("can4", "plot", 800, 600);
	hSignalRatio_pT->Draw();
	hSignalRatio_pT->GetYaxis()->SetTitle("(Chic_1+Chic_2) / Jpsi");
	hSignalRatio_pT->GetYaxis()->SetTitleOffset(1.3);
	hSignalRatio_pT->GetXaxis()->SetTitle("pT");
	can4->SaveAs("RatioChicJpsi.png");


	///////////////////////////// COMPLEMENT

	//cout << "nEntries: " << rdsNominalCompl->numEntries() << endl;
	//rdsNominalCompl = (RooDataSet*)rdsNominalCompl->reduce(mass_windowFit.c_str());
	////RooDataSet* dataOS = (RooDataSet*)myWs.data("rdsNominalCompl")->reduce("pt>20");
	////RooDataSet* dataOS2 = (RooDataSet*)myWs.data("rdsNominalCompl")->reduce("pt<5");

	//RooPlot *massframe2 = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
	//massframe2->SetTitle("mass");
	//rdsNominalCompl->plotOn(massframe2);
	////dataOS->plotOn(massframe2);
	////dataOS2->plotOn(massframe2);
	///*
	//string myPdfNameCompl = "nominalPdfCompl";

	//CreateModelPdf(myWs, myPdfNameCompl);

	////RooFitResult* fitResult = myWs.pdf(myPdfNameCompl.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(true), SumW2Error(true), Range("MassWindow"), NumCPU(1), Save(true));
	//RooFitResult* fitResult2 = myWs.pdf(myPdfNameCompl.c_str())->fitTo(*myWs.data("rdsNominalCompl"), Extended(true), SumW2Error(true), Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true));
	//fitResult2->Print("v");
	//myWs.import(*fitResult2, Form("fitResult_%s", myPdfNameCompl.c_str()));


	//myWs.pdf(myPdfNameCompl.c_str())->plotOn(massframe2);
	//myWs.pdf(myPdfNameCompl.c_str())->paramOn(massframe2, Layout(0.55, 0.9, 0.99));
	//myWs.pdf(myPdfNameCompl.c_str())->plotOn(massframe2, Components("background"), LineStyle(kDashed));
	//myWs.pdf(myPdfNameCompl.c_str())->plotOn(massframe2, Components("chic1"), LineStyle(kDashed), LineColor(kRed));
	//myWs.pdf(myPdfNameCompl.c_str())->plotOn(massframe2, Components("chic2"), LineStyle(kDashed), LineColor(kGreen));
	//cout << endl << endl << "HERE" << endl << endl;
	//*///


	//TCanvas *cTest2 = new TCanvas("cTest2", "cTest2", 1000, 480);
	//gPad->SetLeftMargin(0.15);
	//massframe2->GetYaxis()->SetTitleOffset(1.6);
	//massframe2->Draw();


	//cTest2->SaveAs("CanvasCTestCompl_MC.png");












	//// Declare variables x,mean,sigma with associated name, title, initial value and allowed range
	//RooRealVar x("x", "x", -10, 10);
	//RooRealVar mean("mean", "mean of gaussian", 1, -10, 10);
	//RooRealVar sigma("sigma", "width of gaussian", 1, 0.1, 10);

	//// Build gaussian pdf in terms of x,mean and sigma
	//RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

	////RooAbsPdf *errorFunc = bindPdf("erf", TMath::Erf, x);



	//// Print erf definition
	////errorFunc->Print();


	//// Construct plot frame in 'x'
	//RooPlot *xframe = x.frame(Title("Gaussian pdf."));

	//// P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
	//// ---------------------------------------------------------------------------

	//// Plot gauss in frame (i.e. in x)
	//gauss.plotOn(xframe);
	////errorFunc->plotOn(xframe);

	//// Change the value of sigma to 3
	//sigma.setVal(3);

	//// Plot gauss in frame (i.e. in x) and draw frame on canvas
	//gauss.plotOn(xframe, LineColor(kRed));

	//// G e n e r a t e   e v e n t s
	//// -----------------------------

	//// Generate a dataset of 1000 events in x from gauss
	//RooDataSet *data = gauss.generate(x, 10000);

	//// Make a second plot frame in x and draw both the
	//// data and the pdf in the frame
	//RooPlot *xframe2 = x.frame(Title("Gaussian pdf with data"));
	//data->plotOn(xframe2);
	//gauss.plotOn(xframe2);

	//// F i t   m o d e l   t o   d a t a
	//// -----------------------------

	//// Fit pdf to data
	//gauss.fitTo(*data);

	//// Print values of mean and sigma (that now reflect fitted values and errors)
	//mean.Print();
	//sigma.Print();

	//// Draw all frames on a canvas
	//TCanvas *c = new TCanvas("rf101_basics", "rf101_basics", 800, 400);
	//c->Divide(2);
	//c->cd(1);
	//gPad->SetLeftMargin(0.15);
	//xframe->GetYaxis()->SetTitleOffset(1.6);
	//xframe->Draw();
	//c->cd(2);
	//gPad->SetLeftMargin(0.15);
	//xframe2->GetYaxis()->SetTitleOffset(1.6);
	//xframe2->Draw();
	//
	//
	//c->SaveAs("CanvasC.png");
	//
	//
	//
	//
	//









	/*
	TCanvas* can1 = new TCanvas("can1", "plot", 1200, 800);

	TH2F* h1 = hdxy_dz->Clone();
	//cout<<h1->Integral(101, 300, 1, 30)<<endl;

	h1->GetXaxis()->SetTitle("abseta");
	h1->GetYaxis()->SetTitle("pt");
	can1->SetLogz();
	h1->Draw("colz");
	*/



	//hSignalpt1->Write();
	//hSignalpt2->Write();
	//hSignalpt3->Write();
	//hSignalpt4->Write();
	//hSignalpt5->Write();
	//hSignalpt6->Write();
	//hSignalpt7->Write();

	//hSignal_m_pt->Write();
	//hSignal_m_rap->Write();
	//hSignal_m_phi->Write();
	//hSignal_m_nTr->Write();
	//hSignal_m_Et->Write();

	//hmuon_M->Write();
	//hmuon_Mzoom->Write();
	//hphoton_M->Write();

	//hdimuon_pt->Write();
	//hphoton_pt->Write();

	//hchicCand_pt_rap->Write();
	//hdimuon_pt_rap->Write();
	//hmuonP_pt_ps->Write();
	//hmuonN_pt_ps->Write();
	//hphoton_pt_rap->Write();

	//hdimuonAfterCuts_pt_rap->Write();
	//hmuonPAfterCuts_pt_ps->Write();
	//hmuonNAfterCuts_pt_ps->Write();
	//hphotonAfterCuts_pt_rap->Write();

	//hchi_restricted_pt_rap->Write();
	//hchi_restrictedM_pt->Write();
	//hchi_restrictedM_dimuonpt->Write();
	//hchi_restrictedM_mPpt->Write();
	//hchi_restrictedM_mNpt->Write();
	//hchi_restrictedM_Photpt->Write();

	//if (flagGenerateRds == true)
	{
		TCanvas* can2 = new TCanvas("can2", "plot", 1200, 800);

		TH1D* h2 = dynamic_cast<TH1D*>(hSignal->Clone());

		h2->GetXaxis()->SetTitle("M_inv [GeV]");
		h2->GetYaxis()->SetTitle("counts");
		//can2->SetLogy();
		h2->Draw("");
		TH1D* h3 = dynamic_cast<TH1D*>(hSignal_SS->Clone());
		TH1D* h4 = dynamic_cast<TH1D*>(hSignal_rot->Clone());
		TH1D* h5 = dynamic_cast<TH1D*>(hSignal_rotGamma->Clone());
		TH1D* h6 = dynamic_cast<TH1D*>(hSignal_SB->Clone());
		h3->SetLineColor(kRed);
		h4->SetLineColor(kGreen + 2);
		h5->SetLineColor(kYellow + 2);
		h6->SetLineColor(kMagenta);
		h3->Draw("same");
		h4->Draw("same");
		h5->Draw("same");
		h6->Draw("same");

		TLegend* leg2 = new TLegend(0.7, 0.2, 0.88, 0.4, "");
		leg2->AddEntry(h3, "Same-sign", "l");
		leg2->AddEntry(h4, "Rotation (muon)", "l");
		leg2->AddEntry(h5, "Rotation (gamma)", "l");
		leg2->AddEntry(h6, "Side-band", "l");
		leg2->Draw();


		TCanvas* can3 = new TCanvas("can3", "plot", 1200, 800);

		TH1D* hs2 = dynamic_cast<TH1D*>(hSignal->Clone());
		hs2->GetXaxis()->SetTitle("M_inv [GeV]");
		hs2->GetYaxis()->SetTitle("counts [Arb. scaling]");
		//can2->SetLogy();
		hs2->Draw("");
		TH1D* hs3 = dynamic_cast<TH1D*>(hSignal_SS->Clone());
		TH1D* hs4 = dynamic_cast<TH1D*>(hSignal_rot->Clone());
		TH1D* hs5 = dynamic_cast<TH1D*>(hSignal_rotGamma->Clone());
		TH1D* hs6 = dynamic_cast<TH1D*>(hSignal_SB->Clone());
		hs3->Scale(hs2->GetBinContent(100) / (hs3->GetBinContent(100) + 0.01));
		hs4->Scale(hs2->GetBinContent(100) / (hs4->GetBinContent(100) + 0.01));
		hs5->Scale(hs2->GetBinContent(100) / (hs5->GetBinContent(100) + 0.01));
		hs6->Scale(hs2->GetBinContent(100) / (hs6->GetBinContent(100) + 0.01));

		hs3->SetLineColor(kRed);
		hs4->SetLineColor(kGreen + 2);
		hs5->SetLineColor(kYellow + 2);
		hs6->SetLineColor(kMagenta);
		hs3->Draw("same");
		hs4->Draw("same");
		hs5->Draw("same");
		hs6->Draw("same");

		TLegend* leg3 = new TLegend(0.7, 0.2, 0.88, 0.4, "");
		leg3->AddEntry(hs3, "Same-sign", "l");
		leg3->AddEntry(hs4, "Rotation (muon)", "l");
		leg3->AddEntry(hs5, "Rotation (gamma)", "l");
		leg3->AddEntry(hs6, "Side-band", "l");
		leg3->Draw();


		//hSignal->Sumw2();
		//hSignal2->Sumw2();
		//hSignal3->Sumw2();

		TFile* fout = new TFile(fileOut, "RECREATE");
		can2->Write();
		can3->Write();
		hSignal->Write();
		hSignal_SS->Write();
		hSignal_rot->Write();
		hSignal_rotGamma->Write();
		hSignal_SB->Write();
		hchic_M->Write();
		hdimuon_M->Write();
		hchic_M_SS->Write();
		hdimuon_M_SS->Write();
		hchic_M_rot->Write();
		hdimuon_M_rot->Write();
		hchic_M_rotGamma->Write();
		hchic_M_SB->Write();
		hntracks_inEvent->Write();
	}

}

