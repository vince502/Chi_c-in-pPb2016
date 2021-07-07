
// Macro to obtain the constraints from the MC.
// Based on macro to analyze the chi_c. Therefore some unnecessary structure. Starting with the event tree


#include <iostream>
#include <fstream>
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
#include "RooAddPdf.h"
#include "RooMsgService.h"

#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../ChiTreeInit.C"


using namespace std;
using namespace RooFit;

ofstream file_log;

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

const int nMassBinsJpsi = 150;
const double mass_windowJpsi_l = 2.5;
const double mass_windowJpsi_h = 4.0;
const double mass_windowFitJpsi_l = 2.5;
const double mass_windowFitJpsi_h = 4.0;
const string mass_windowFitJpsi = "rvmassJpsi>2.5 && rvmassJpsi<4.0";

double bins_pT[] = {6, 9, 12, 18, 30};
int  nbins_pT = sizeof(bins_pT) / sizeof(double) - 1;

double bins_y[] = {-2.4, -1.6, -1.0, 0, 1.0, 1.6, 2.4 };
int  nbins_y = sizeof(bins_y) / sizeof(double) - 1;

double bins_nTrk[] = { 0, 50, 100, 150, 200, 300, 400 };
int  nbins_nTrk = sizeof(bins_nTrk) / sizeof(double) - 1;

// parameters

TGraphAsymmErrors* gAsChic_alpha_pT, *gAsChic_alpha_y, *gAsChic_alpha_nTrk;
TGraphAsymmErrors* gAsChic_n_pT, *gAsChic_n_y, *gAsChic_n_nTrk;
TGraphAsymmErrors* gAsChic_sigmaRat_pT, *gAsChic_sigmaRat_y, *gAsChic_sigmaRat_nTrk;
TGraphAsymmErrors* gAsChic_sigma1_pT, *gAsChic_sigma1_y, *gAsChic_sigma1_nTrk;
TGraphAsymmErrors* gAsChic_mean1_pT, *gAsChic_mean1_y, *gAsChic_mean1_nTrk;


TLorentzVector* LVchic, *LVdimuon, *LVconv, *LVmuon1, *LVmuon2;
TLorentzVector* LVchic_rot, *LVchic_rotGamma, *LVdimuon_rot, *LVconv_rot, *LVmuon1_rot, *LVmuon2_rot;
TLorentzVector LVaux;

int GetRatio (TGraphAsymmErrors* gAsResult, TGraphAsymmErrors* gAsNum, TGraphAsymmErrors* gAsDen) // dependent on root version, this is prior 6.20
{
	if (gAsResult->GetN() != gAsResult->GetN() || gAsResult->GetN() != gAsResult->GetN()) {
		cout << "ERROR: Not the same number of N bins, can't get ratio. Will crash" << endl;
		return -1;
	}
	double den_x, den_y, num_x, num_y;

	for (int i = 0; i < gAsResult->GetN(); i++)
	{
		gAsDen->GetPoint(i, den_x, den_y);
		gAsNum->GetPoint(i, num_x, num_y);

		if (den_y > 0.00001) {
			gAsResult->SetPoint(i, num_x, num_y / den_y);
			gAsResult->SetPointEYhigh(i, sqrt((gAsNum->GetErrorYhigh(i) / num_y)* (gAsNum->GetErrorYhigh(i) / num_y) + (gAsDen->GetErrorYlow(i) / den_y)* (gAsDen->GetErrorYlow(i) / den_y)) * (num_y / den_y)); //treating them as independent
			gAsResult->SetPointEYlow(i, sqrt((gAsNum->GetErrorYlow(i) / num_y)* (gAsNum->GetErrorYlow(i) / num_y) + (gAsDen->GetErrorYhigh(i) / den_y)* (gAsDen->GetErrorYhigh(i) / den_y)) * (num_y / den_y)); //treating them as independent

		}
		else {
			gAsResult->SetPoint(i, num_x, -1);
			gAsResult->SetPointError(i, 0,0,0,0);
		}
	}

	return 0;
}

int SetPointFromFit(TGraphAsymmErrors* gAsResult, string fitVarName, RooWorkspace& Ws, int i, double* bins)
{
	gAsResult->SetPoint(i, ((bins[i] + bins[i + 1]) / 2.0), Ws.var(fitVarName.c_str())->getValV()); //placing the point in the middle of the bin
	gAsResult->SetPointEYhigh(i, Ws.var(fitVarName.c_str())->getErrorHi());
	gAsResult->SetPointEYlow(i, -(Ws.var(fitVarName.c_str())->getErrorLo()));
	return 0;
}


bool CreateModelPdf(RooWorkspace& Ws, string pdfName, int intConstrainedFit = 0)
{


	//if (pdfName.find("nominalPdf") != string::npos)
	if (pdfName.compare("nominalPdf") == 0)
	{

		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		Ws.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		Ws.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		Ws.factory("SUM::signal(c2toc1[0.4,0.1,1.]*chic2, chic1)");
		//Ws.factory("Exponential::expbkg(rvmass,e_1[-0.2,-0.5,0])");
		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("RooChebychev::background(rvmass, a1[0,-10,10])");

		if (intConstrainedFit == 1) {
			Ws.var("alpha")->setConstant(true);
			Ws.var("n")->setConstant(true);
		}
		//Ws.var("mean1")->setConstant(true);
		//Ws.var("mean2")->setConstant(true);
		//Ws.var("erf_offset")->setConstant(true);
		//Ws.var("erf_sigma")->setConstant(true);
		//Ws.factory("")
		Ws.factory("SUM::nominalPdf(nsig[50,0,1000000]*signal, nbkg[2000,0,100000]*background)");
	}
	else if (pdfName.compare("nominalPdfJpsi") == 0)
		//else if (pdfName.find("nominalPdfJpsi") != string::npos)
	{
		//Ws.factory("CBShape::signalJpsi(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08], alphaJpsi[1.85, 0.1, 50], nJpsi[1.7, 0.2, 50])");
		Ws.factory("CBShape::Jpsi1(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08], alphaJpsi[1.85, 0.1, 50], nJpsi[1.7, 0.2, 50])");
		Ws.factory("CBShape::Jpsi2(rvmassJpsi, mean1Jpsi, sigma2Jpsi[0.03, 0.003, 0.08], alphaJpsi, nJpsi)");
		Ws.factory("CBShape::psi2(rvmassJpsi, mean2Jpsi[3.686, 3.50, 3.75], sigma1Jpsi, alphaJpsi, nJpsi)");
		Ws.factory("SUM::Jpsi(ratioJpsi[0.5,0.00,1.0]*Jpsi1, Jpsi2)");
		Ws.factory("SUM::signalJpsi(ratioPsi[0.1,0.00,2.0]*psi2, Jpsi)");

		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("Exponential::backgroundJpsi(rvmassJpsi,e_1Jpsi[-0.2,-2.5,0])");
		//Ws.factory("RooChebychev::backgroundJpsi(rvmassJpsi, a1Jpsi[-0.05,-40.0, 40.0])");
		//Ws.var("alphaJpsi")->setConstant(true);
		//Ws.var("nJpsi")->setConstant(true);
		//Ws.var("mean1")->setConstant(true);
		//Ws.var("mean2")->setConstant(true);
		//Ws.factory("")
		Ws.factory("SUM::nominalPdfJpsi(nsigJpsi[5000,0,2000000]*signalJpsi, nbkgJpsi[2000,0,1000000]*backgroundJpsi)");
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

/*bool WriteConstrainedFit(RooWorkspace& Ws, string pdfName)
{
	if (pdfName.compare("nominalPdf") == 0)
	{
		cout << "signal " << Ws.var("nsig")->getValV() << endl;
		gAsResult->SetPoint(i, ((bins[i] + bins[i + 1]) / 2.0), Ws.var("nsig")->getValV()); //placing the point in the middle of the bin
		gAsResult->SetPointEYhigh(i, Ws.var("nsig")->getErrorHi());
		gAsResult->SetPointEYlow(i, -(Ws.var("nsig")->getErrorLo()));



		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.03], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
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
		Ws.factory("SUM::nominalPdf(nsig[50,0,1000000]*signal, nbkg[2000,0,100000]*background)");
	}
}*/

bool RefreshModel(RooWorkspace& Ws, string pdfName) // attempt to prevent the fits in the differential bins to get stuck in a weird state when they stop converge properly
{
	file_log << Ws.var("mean1Jpsi")->getError() << endl;
	file_log << Ws.var("sigma1Jpsi")->getError() << endl;

	Ws.var("mean1Jpsi")->removeError();
	Ws.var("sigma1Jpsi")->removeError();
	Ws.var("alphaJpsi")->removeError();
	Ws.var("nJpsi")->removeError();
	//Ws.var("a1Jpsi")->removeError();
	Ws.var("e_1Jpsi")->removeError();
	Ws.var("nsigJpsi")->removeError();
	Ws.var("nbkgJpsi")->removeError();

	Ws.var("mean1Jpsi")->setVal(3.097);
	Ws.var("sigma1Jpsi")->setVal(0.03);
	Ws.var("alphaJpsi")->setVal(1.85);
	Ws.var("nJpsi")->setVal(1.7);
	//Ws.var("a1Jpsi")->setVal(-0.05);
	Ws.var("e_1Jpsi")->setVal(-0.2);
	Ws.var("nsigJpsi")->setVal(5000);
	Ws.var("nbkgJpsi")->setVal(2000);

	return true;
}



int FitRooDataSet(TGraphAsymmErrors* gAsResult, double* bins, int nbins, RooRealVar* rvmass, RooWorkspace& Ws, bool isJpsi = false, string binVarName = "", string myPdfName = "", string extraCut = "", string canvasName = "", int intConstrainedFit = 0) { //uses global variables
	//constrained fit: 0 no constraints, 1 yes, 2 no constraints, but write the fit results as constraints (for fitting MC)
	cout << endl << "IN FITTING " << binVarName << "   " << nbins << endl << endl;
	TCanvas *cFit = new TCanvas("cFit", "cFit", 1000, 600);
	gPad->SetLeftMargin(0.15);


	for (int i = 0; i < nbins; i++) {
		cout << "Bins " << bins[i] << "  to  " << bins[i + 1] << endl;

		RooPlot* massframeBin;
		if (isJpsi == false) { massframeBin = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins); }
		else { massframeBin = rvmass->frame(mass_windowFitJpsi_l, mass_windowFitJpsi_h, nMassBinsJpsi); }
		massframeBin->SetTitle("mass");
		TString TstrCut = binVarName + TString::Format(" > %f", bins[i]) + " && " + binVarName + TString::Format(" < %f", bins[i + 1]) + extraCut;
		cout << TstrCut << endl;
		string strCut = TstrCut.Data();
		RooDataSet* rdsDataBin;
		if (isJpsi == false) { rdsDataBin = (RooDataSet*)Ws.data("rdsNominal")->reduce(strCut.c_str()); }
		else {
			rdsDataBin = (RooDataSet*)Ws.data("rdsNominalJpsi")->reduce(strCut.c_str()); 
			RefreshModel(Ws, myPdfName);
		}

		rdsDataBin->plotOn(massframeBin);

		RooFitResult* fitResultBin = Ws.pdf(myPdfName.c_str())->fitTo(*rdsDataBin, Extended(true), SumW2Error(true), NumCPU(1));// Range(mass_windowFit_l, mass_windowFit_h), Save(true));
		//fitResultBin->Print("v");
		//Ws.import(*fitResultBin, TString::Format("fitResult_%s", myPdfNameBin.c_str()));
		if (isJpsi == false) {
			cout << "signal " << Ws.var("nsig")->getValV() << endl;
			SetPointFromFit(gAsResult, "nsig", Ws, i, bins);

			if (intConstrainedFit == 2) {
				if (myPdfName.compare("nominalPdf") == 0) {
					if (binVarName.compare("rvpt") == 0){
						SetPointFromFit(gAsChic_alpha_pT, "alpha", Ws, i, bins);
						SetPointFromFit(gAsChic_n_pT, "n", Ws, i, bins);
						SetPointFromFit(gAsChic_sigmaRat_pT, "sigmaRatio", Ws, i, bins);
						SetPointFromFit(gAsChic_sigma1_pT, "sigma1", Ws, i, bins);
						SetPointFromFit(gAsChic_mean1_pT, "mean1", Ws, i, bins);

					} 
					else if (binVarName.compare("rvrap") == 0) {
						SetPointFromFit(gAsChic_alpha_y, "alpha", Ws, i, bins);
						SetPointFromFit(gAsChic_n_y, "n", Ws, i, bins);
						SetPointFromFit(gAsChic_sigmaRat_y, "sigmaRatio", Ws, i, bins);
						SetPointFromFit(gAsChic_sigma1_y, "sigma1", Ws, i, bins);
						SetPointFromFit(gAsChic_mean1_y, "mean1", Ws, i, bins);

					}
					else if (binVarName.compare("rvntrack") == 0) {
						SetPointFromFit(gAsChic_alpha_nTrk, "alpha", Ws, i, bins);
						SetPointFromFit(gAsChic_n_nTrk, "n", Ws, i, bins);
						SetPointFromFit(gAsChic_sigmaRat_nTrk, "sigmaRatio", Ws, i, bins);
						SetPointFromFit(gAsChic_sigma1_nTrk, "sigma1", Ws, i, bins);
						SetPointFromFit(gAsChic_mean1_nTrk, "mean1", Ws, i, bins);

					}
					else cout << "No match for the variable" << endl;

				}
				else cout << "The pdf not matching" << endl;
			
			}

		}
		else {
			cout << "signal Jpsi " << Ws.var("nsigJpsi")->getValV() << endl;
			file_log << i << " signal Jpsi " << Ws.var("nsigJpsi")->getValV() << endl;
			SetPointFromFit(gAsResult, "nsigJpsi", Ws, i, bins);
		}

		Ws.pdf(myPdfName.c_str())->plotOn(massframeBin);
		Ws.pdf(myPdfName.c_str())->paramOn(massframeBin, Layout(0.55));
		
		if (isJpsi == false) {
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("background"), LineStyle(kDashed));
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic1"), LineStyle(kDashed), LineColor(kRed));
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic2"), LineStyle(kDashed), LineColor(kGreen));
		}
		else {
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("backgroundJpsi"), LineStyle(kDashed));
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("Jpsi"), LineStyle(kDashed), LineColor(kRed));
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("psi2"), LineStyle(kDashed), LineColor(kGreen));

		}
		massframeBin->Draw();
		cFit->SaveAs(((string)"FitterOutput/FitResult_" + (isJpsi?"Jpsi_":"Chic_") + canvasName + "_"  + binVarName + Form("_%.1f_%.1f.png", bins[i], bins[i + 1])).c_str());

	}




	return 0;
}



///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////


void ConstraintProducer(bool flagGenerateRds = true, bool flagRunFits = true,  const char* fileIn = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-MC8_BothDir.root", const char* fileOut = "Chi_c_output_MC8_test.root", const char* fileRds = "rds_MC8_Constraint.root", const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
{

	//gStyle->SetOptStat(1111);
	//gStyle->SetOptStat(0);
	setTDRStyle();
	
	bool isMC = true;
	bool flagConstrainedFit = true;

	int intConstrainedFit; //constrained fit: 0 no constraints, 1 yes, 2 no constraints, but write the fit results as constraints (for fitting MC)
	if (flagConstrainedFit == false) { intConstrainedFit = 0; }
	else if (isMC) { intConstrainedFit = 2; }
	else { intConstrainedFit = 1; }

	file_log.open("TestLog.txt", ios::app);
	file_log << endl << endl << "N E W   R U N " << endl << endl;

	TH1D* hSignal_pT = new TH1D("hSignal_pT", "", nbins_pT, bins_pT);
	TH1D* hSignalJpsi_pT = new TH1D("hSignalJpsi_pT", "", nbins_pT, bins_pT);
	TH1D* hSignalRatio_pT = new TH1D("hSignalRatio_pT", "", nbins_pT, bins_pT);

	TGraphAsymmErrors* gAsChic_pT = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT->SetNameTitle("gAsChic_pT","Chic pT dependence");
	TGraphAsymmErrors* gAsJpsi_pT = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT->SetNameTitle("gAsJpsi_pT", "Jpsi pT dependence");
	TGraphAsymmErrors* gAsRatio_pT = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT->SetNameTitle("gAsRatio_pT", "Ratio pT dependence");

	TGraphAsymmErrors* gAsChic_y = new TGraphAsymmErrors(nbins_y);
	gAsChic_y->SetNameTitle("gAsChic_y", "Chic y dependence");
	TGraphAsymmErrors* gAsJpsi_y = new TGraphAsymmErrors(nbins_y);
	gAsJpsi_y->SetNameTitle("gAsJpsi_y", "Jpsi y dependence");
	TGraphAsymmErrors* gAsRatio_y = new TGraphAsymmErrors(nbins_y);
	gAsRatio_y->SetNameTitle("gAsRatio_y", "Ratio y dependence");

	TGraphAsymmErrors* gAsChic_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_nTrk->SetNameTitle("gAsChic_nTrk", "Chic y dependence");
	TGraphAsymmErrors* gAsJpsi_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsJpsi_nTrk->SetNameTitle("gAsJpsi_nTrk", "Jpsi y dependence");
	TGraphAsymmErrors* gAsRatio_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsRatio_nTrk->SetNameTitle("gAsRatio_nTrk", "Ratio y dependence");


	///

	gAsChic_alpha_pT = new TGraphAsymmErrors(nbins_pT);
	gAsChic_alpha_pT->SetNameTitle("gAsChic_alpha_pT", "Chic pT alpha dependence");
	gAsChic_alpha_y = new TGraphAsymmErrors(nbins_y);
	gAsChic_alpha_y->SetNameTitle("gAsChic_alpha_y", "Chic y alpha dependence");
	gAsChic_alpha_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_alpha_nTrk->SetNameTitle("gAsChic_alpha_nTrk", "Chic alpha ntrk dependence");
	gAsChic_n_pT = new TGraphAsymmErrors(nbins_pT);
	gAsChic_n_pT->SetNameTitle("gAsChic_n_pT", "Chic pT n dependence");
	gAsChic_n_y = new TGraphAsymmErrors(nbins_y);
	gAsChic_n_y->SetNameTitle("gAsChic_n_y", "Chic y n dependence");
	gAsChic_n_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_n_nTrk->SetNameTitle("gAsChic_n_nTrk", "Chic ntrk n dependence");
	gAsChic_sigmaRat_pT = new TGraphAsymmErrors(nbins_pT);
	gAsChic_sigmaRat_pT->SetNameTitle("gAsChic_sigmaRat_pT", "Chic pT sigmaRat dependence");
	gAsChic_sigmaRat_y = new TGraphAsymmErrors(nbins_y);
	gAsChic_sigmaRat_y->SetNameTitle("gAsChic_sigmaRat_y", "Chic y sigmaRat dependence");
	gAsChic_sigmaRat_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_sigmaRat_nTrk->SetNameTitle("gAsChic_sigmaRat_nTrk", "Chic sigmaRat ntrk dependence");
	gAsChic_sigma1_pT = new TGraphAsymmErrors(nbins_pT);
	gAsChic_sigma1_pT->SetNameTitle("gAsChic_sigma1_pT", "Chic pT sigma1 dependence");
	gAsChic_sigma1_y = new TGraphAsymmErrors(nbins_y);
	gAsChic_sigma1_y->SetNameTitle("gAsChic_sigma1_y", "Chic y sigma1 dependence");
	gAsChic_sigma1_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_sigma1_nTrk->SetNameTitle("gAsChic_sigma1_nTrk", "Chic sigma1 ntrk dependence");
	gAsChic_mean1_pT = new TGraphAsymmErrors(nbins_pT);
	gAsChic_mean1_pT->SetNameTitle("gAsChic_mean1_pT", "Chic pT mean1 dependence");
	gAsChic_mean1_y = new TGraphAsymmErrors(nbins_y);
	gAsChic_mean1_y->SetNameTitle("gAsChic_mean1_y", "Chic y mean1 dependence");
	gAsChic_mean1_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_mean1_nTrk->SetNameTitle("gAsChic_mean1_nTrk", "Chic mean1 ntrk dependence");


	TH1D* hSignal = new TH1D("hSignal", "", 200, 3, 5);
	TH1D* hSignal_SS = new TH1D("hSignal_SS", "", 200, 3, 5);
	TH1D* hSignal_SB = new TH1D("hSignal_SB", "", 200, 3, 5);
	TH1D* hSignal_rot = new TH1D("hSignal_rot", "", 200, 3, 5);
	TH1D* hSignal_rotGamma = new TH1D("hSignal_rotGamma", "", 200, 3, 5);
	TH1D* hSignal2 = new TH1D("hSignal2", "", 120, 3.2, 3.8);
	TH1D* hSignal3 = new TH1D("hSignal3", "", 120, 3.2, 3.8);



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



		// Various
	TH1D* hntracks_inEvent = new TH1D("hntracks_inEvent", "", 400, 0, 400);

	TFile* f1 = new TFile(fileIn, "READ");

	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");
	if (!event_tree) {
		cout << "Problem with event Tree";
		//break;
	}

	LoadChiBranches(event_tree, isMC);


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
	RooRealVar* rvntrack = new RooRealVar("rvntrack", "ntrack", 0.0, 500.0, "");
	RooRealVar* rvweight = new RooRealVar("rvweight", "weight", 0.0, 10000.0, "");
	RooArgSet*  cols = new RooArgSet(*rvmass, *rvpt, *rvrap, *rvntrack, *rvweight);
	//RooArgSet*  cols = new RooArgSet(*rvmass);
	RooDataSet* rdsNominal = new RooDataSet("rdsNominal", "rdsNominal", *cols, WeightVar(*rvweight), StoreAsymError(*rvmass));
	//RooDataSet* rdsNominal = new RooDataSet("nominal", "nominal", *cols);

	RooDataSet* rdsNominal_chicOne = new RooDataSet("rdsNominal_chicOne", "rdsNominal_chicOne", *cols, WeightVar(*rvweight), StoreAsymError(*rvmass));
	RooDataSet* rdsNominal_chicTwo = new RooDataSet("rdsNominal_chicTwo", "rdsNominal_chicTwo", *cols, WeightVar(*rvweight), StoreAsymError(*rvmass));


	RooRealVar* rvmassJpsi = new RooRealVar("rvmassJpsi", "Inv mass #mu#mu", mass_windowJpsi_l, mass_windowJpsi_h, "GeV/c^{2}");
	RooRealVar* rvptJpsi = new RooRealVar("rvptJpsi", "#mu#mu p_{T}", 0.0, 50.0, "GeV/c");
	RooRealVar* rvrapJpsi = new RooRealVar("rvrapJpsi", "#mu#mu y", -2.4, 2.4, "");
	RooRealVar* rvntrackJpsi = new RooRealVar("rvntrackJpsi", "ntrack", 0.0, 500.0, "");
	RooRealVar* rvweightJpsi = new RooRealVar("rvweightJpsi", "weight", 0.0, 10000.0, "");
	RooArgSet*  colsJpsi = new RooArgSet(*rvmassJpsi, *rvptJpsi, *rvrapJpsi, *rvntrackJpsi, *rvweightJpsi);
	RooDataSet* rdsNominalJpsi = new RooDataSet("rdsNominalJpsi", "rdsNominalJpsi", *colsJpsi, WeightVar(*rvweightJpsi), StoreAsymError(*rvmassJpsi));
	


	if (flagGenerateRds == true) {


		//load corrections
		TFile* fCor = new TFile(fileCorrection, "READ");

		TH2D* hWeightChic = (TH2D*)fCor->Get("h_chiAccEff_rat");
		TH2D* hWeightJpsi = (TH2D*)fCor->Get("h_JpsiAccEff_rat");


		long nchicCounter = 0, nchicCounterPass = 0, nchicCounterPassMass = 0, muon1ptCounter = 0, muon2ptCounter = 0;
		bool passDimSel = false;
		bool passDimSelTight = false;
		Long64_t nentries = event_tree->GetEntries();
		//if (nentries > 50000) { nentries = 50000; }
		cout << nentries << endl;
		for (Long64_t i = 0; i < nentries; i++)
		{

			event_tree->GetEntry(i);
			if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }
			if (isMC==true && (gen_pdgId->at(0) != PythCode_chic1 && gen_pdgId->at(0) != PythCode_chic2)) { continue; } //remove chic0 from MC, expecting one gen chic per event
			//if (isMC == true && (gen_pdgId->at(0) != PythCode_chic1)) { continue; } // only chic1
			//if (isMC == true && (gen_pdgId->at(0) != PythCode_chic2)) { continue; } // only chic2

			////////////////////////
			//// V A R I O U S  ////
			///////////////////////
			hntracks_inEvent->Fill(ntracks_inEvent);
			//cout << ntracks_inEvent << endl;

			///////////////////////
			/////   C H I C ///////
			/////////////////////

			for (int iChi = 0; iChi < chi_p4->GetEntriesFast(); iChi++)
			{
				++nchicCounter;
				// check Acceptance and Cuts
				int dimuonPos = chi_daughterJpsi_position->at(iChi);
				//if (DimuonSelectionPass(dimuonPos) == false) continue;
				passDimSel = DimuonSelectionPass(dimuonPos);
				passDimSelTight = DimuonSelectionPassTight(dimuonPos, ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity());
				// muon cuts
				int muon1Pos = dimuon_muon1_position->at(dimuonPos);
				int muon2Pos = dimuon_muon2_position->at(dimuonPos);
				if (MuonAcceptance(muon_eta->at(muon1Pos), muon_pt->at(muon1Pos)) == false) continue;
				if (MuonAcceptance(muon_eta->at(muon2Pos), muon_pt->at(muon2Pos)) == false) continue;
				if (MuonSelectionPass(muon1Pos) == false) continue;
				if (MuonSelectionPass(muon2Pos) == false) continue;
				// photon
				int convPos = chi_daughterConv_position->at(iChi);
				if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) continue;
				if (PhotSelectionPass(convPos) == false) continue;
				if (passDimSel == true) { ++nchicCounterPass; } // SelectionsPassed
				// Get Lorentz V
				LVchic = (TLorentzVector*)chi_p4->At(iChi);
				LVdimuon = (TLorentzVector*)dimuon_p4->At(dimuonPos);
				LVconv = (TLorentzVector*)conv_p4->At(convPos);

				double pT_chi = chi_pt->at(iChi);
				double rap_chi = LVchic->Rapidity();
				double m_chi = LVchic->M();
				double dimuonM = LVdimuon->M();
				double Mdiff = m_chi - dimuonM + 3.097;// Assume J/psi mass

				if (passDimSel == true && dimuonM>2.95 && dimuonM<3.2 && Mdiff > mass_window_l && Mdiff < mass_window_h) ++nchicCounterPassMass;


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
					hchic_M_SB->Fill(m_chi);
					hSignal_SB->Fill(m_chi - dimuonM + 3.097);
				}


				if (dimuonM<2.95 || dimuonM>3.2) continue; //require narrow dimuon mass

				if (passDimSel == true) {
					hchic_M->Fill(m_chi); // just raw M
					hchic_M_rotGamma->Fill(LVchic_rotGamma->M());
				}
				if (dimuon_charge->at(dimuonPos) != 0) {
					hchic_M_SS->Fill(m_chi);
				}


				if (passDimSel == true) {
					hSignal->Fill(Mdiff);
					hSignal_rotGamma->Fill(LVchic_rotGamma->M() - LVdimuon->M() + 3.097);
					//cout << nchicCounterPass << endl;

					//roofit:
					if (Mdiff > mass_window_l && Mdiff < mass_window_h) {
						rvmass->setVal(Mdiff);
						rvpt->setVal(pT_chi);
						rvrap->setVal(rap_chi);
						rvntrack->setVal(ntracks_inEvent); //for now using total
						//rvweight->setVal(1); // no weights for now
						double accEff_chi = hWeightChic->GetBinContent(hWeightChic->FindBin(abs(rap_chi), pT_chi));
						//cout << "Test rap: " << rap_chi << "  pt  " << pT_chi << " and value accEff " << accEff_chi << endl;
						if (accEff_chi > 0.00001) {// binning should ensure enough statistics, but if we have empty bin, just set weight to be 1
							rvweight->setVal(1 / accEff_chi);
						}
						else rvweight->setVal(1); 
						if (isMC) {
							rvweight->setVal(1); //no weights if MC
						}

						rdsNominal->add(*cols, rvweight->getVal());
						if (gen_pdgId->at(0) == PythCode_chic1) { rdsNominal_chicOne->add(*cols, rvweight->getVal()); }
						if (gen_pdgId->at(0) == PythCode_chic2) { rdsNominal_chicTwo->add(*cols, rvweight->getVal()); }
					}

				}
				if (dimuon_charge->at(dimuonPos) != 0) { hSignal_SS->Fill(Mdiff); }

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

						double accEff_Jpsi = hWeightJpsi->GetBinContent(hWeightJpsi->FindBin(abs(LVdimuon->Rapidity()), dimuon_pt->at(iJpsi)));
						//cout << "Test rap: " << rap_chi << "  pt  " << pT_chi << " and value accEff " << accEff_chi << endl;
						if (accEff_Jpsi > 0.00001) {// binning should ensure enough statistics, but if we have empty bin, just set weight to be 1
							rvweightJpsi->setVal(1 / accEff_Jpsi);
						}
						else rvweightJpsi->setVal(1);

						//rvweightJpsi->setVal(1); // no weights for now

						if (isMC) {
							rvweightJpsi->setVal(1); //no weights if MC
						}


						rdsNominalJpsi->add(*colsJpsi, rvweightJpsi->getVal());


					}

				}

			} //end of Jpsi loop





		} //end of event loop

		cout << "Processed " << nchicCounter << " chic candidates, " << nchicCounterPass << " passed the selections, that is " << ((double)nchicCounterPass / (double)nchicCounter * 100) << " percent" << endl;
		cout << "Processed " << nchicCounter << " chic candidates, " << nchicCounterPassMass << " passed the selections and mass window requirement, that is " << ((double)nchicCounterPassMass / (double)nchicCounter * 100) << " percent" << endl;

		// Save rds
		TFile* DBFile = new TFile(fileRds, "RECREATE");
		DBFile->cd();
		rdsNominal->Write();
		rdsNominal_chicOne->Write();
		rdsNominal_chicTwo->Write();
		rdsNominalJpsi->Write();
		DBFile->Write(); DBFile->Close(); delete DBFile;
		f1->cd();
		myWs.import(*rdsNominal);
		myWs.import(*rdsNominal_chicOne);
		myWs.import(*rdsNominal_chicTwo);
		myWs.import(*rdsNominalJpsi);
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
		cout << "nEntries: " << rdsNominal->numEntries() << endl;
		myWs.import(*rdsNominal);
		rdsNominalJpsi = (RooDataSet*)DBFile->Get("rdsNominalJpsi");
		cout << "nEntriesJpsi: " << rdsNominalJpsi->numEntries() << endl;
		myWs.import(*rdsNominalJpsi);

		rdsNominal_chicOne = (RooDataSet*)DBFile->Get("rdsNominal_chicOne");
		myWs.import(*rdsNominal_chicOne);
		rdsNominal_chicTwo = (RooDataSet*)DBFile->Get("rdsNominal_chicTwo");
		myWs.import(*rdsNominal_chicTwo);


		DBFile->Close(); delete DBFile;
	}


	if (flagRunFits == true) {


		//////////////////////////////////////////  
		/////////   R  o  o  f  i  t   //////////
		///////////////////////////////////////////

		//RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
		//ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
		
		cout << "nEntries: " << rdsNominal->numEntries() << endl;
		rdsNominal = (RooDataSet*)rdsNominal->reduce(mass_windowFit.c_str());
		rdsNominal_chicOne = (RooDataSet*)rdsNominal_chicOne->reduce(mass_windowFit.c_str());
		rdsNominal_chicTwo = (RooDataSet*)rdsNominal_chicTwo->reduce(mass_windowFit.c_str());

		// construct the overall dataset for simultaneous fits

		RooCategory rCatChi("rCatChi", "rCatChi");
		rCatChi.defineType("ChicOne");
		rCatChi.defineType("ChicTwo");

		cout << "HERE" << endl;

		//RooDataSet* rdsNominal_chicBoth = new RooDataSet("rdsNominal_chicBoth", "rdsNominal_chicBoth", RooArgSet(*(myWs.var("rvmass")), *(myWs.var("rvpt")), *(myWs.var("rvrap")), *(myWs.var("rvntrack"))), Index(rCatChi), Import("ChicOne", *rdsNominal_chicOne), Import("ChicTwo", *rdsNominal_chicTwo));
		RooDataSet* rdsNominal_chicBoth = new RooDataSet("rdsNominal_chicBoth", "rdsNominal_chicBoth", *rvmass, Index(rCatChi), Import("ChicOne", *rdsNominal_chicOne), Import("ChicTwo", *rdsNominal_chicTwo));

		//some issue with weight var, as it is not needed, dropping it
		cout << "******** New Combined Dataset ***********" << endl;
		rdsNominal_chicBoth->Print();
		myWs.import(*rdsNominal_chicBoth);

		RooPlot *massframeCheck = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
		massframeCheck->SetTitle("mass");
		rdsNominal->plotOn(massframeCheck, MarkerStyle(21), MarkerColor(kBlack));
		rdsNominal_chicBoth->plotOn(massframeCheck, MarkerStyle(20),	MarkerColor(kBlue));
		rdsNominal_chicOne->plotOn(massframeCheck, MarkerStyle(20), MarkerColor(kRed));
		rdsNominal_chicTwo->plotOn(massframeCheck, MarkerStyle(20), MarkerColor(kGreen));

		TCanvas *cCrosscheck = new TCanvas("cCrosscheck", "cCrosscheck", 1000, 600);
		gPad->SetLeftMargin(0.15);
		massframeCheck->GetYaxis()->SetTitleOffset(1.3);
		massframeCheck->Draw();
		
		cCrosscheck->SaveAs("FitterOutput/CanvasMassComparison.png");


		// C r e a t e   m o d e l   f o r   p h y s i c s   s a m p l e
	// -------------------------------------------------------------

	// Create observables
		//RooRealVar x("x", "x", -8, 8);

		// Construct signal pdf
		RooRealVar mean("mean", "mean", 3.5, 3, 4);
		RooRealVar sigma("sigma", "sigma", 0.03, 0.01, 0.1);
		RooGaussian gx("gx", "gx", *rvmass, mean, sigma);

		// Construct background pdf
		RooRealVar a0("a0", "a0", -0.1, -1, 1);
		RooRealVar a1("a1", "a1", 0.004, -1, 1);
		RooChebychev px("px", "px", *rvmass, RooArgSet(a0, a1));

		// Construct composite pdf
		RooRealVar f("f", "f", 0.2, 0., 1.);
		RooAddPdf model("model", "model", RooArgList(gx, px), f);



		// C r e a t e   m o d e l   f o r   c o n t r o l   s a m p l e
		// --------------------------------------------------------------

		// Construct signal pdf. 
		// NOTE that sigma is shared with the signal sample model
		RooRealVar mean_ctl("mean_ctl", "mean_ctl", 3.5, 3, 4);
		RooGaussian gx_ctl("gx_ctl", "gx_ctl", *rvmass, mean_ctl, sigma);

		// Construct the background pdf
		RooRealVar a0_ctl("a0_ctl", "a0_ctl", -0.1, -1, 1);
		RooRealVar a1_ctl("a1_ctl", "a1_ctl", 0.5, -0.1, 1);
		RooChebychev px_ctl("px_ctl", "px_ctl", *rvmass, RooArgSet(a0_ctl, a1_ctl));

		// Construct the composite model
		RooRealVar f_ctl("f_ctl", "f_ctl", 0.5, 0., 1.);
		RooAddPdf model_ctl("model_ctl", "model_ctl", RooArgList(gx_ctl, px_ctl), f_ctl);



		// C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s 
		// ---------------------------------------------------------------------------

		//// Define category to distinguish physics and control samples events
		//RooCategory sample("sample", "sample");
		//sample.defineType("physics");
		//sample.defineType("control");

		// Construct combined dataset in (*rvmass,sample)
		RooDataSet combData("combData", "combined data", *rvmass, Index(rCatChi), Import("ChicOne", *rdsNominal_chicOne), Import("ChicTwo", *rdsNominal_chicTwo));


		// C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
		// -----------------------------------------------------------------------------------

		// Construct a simultaneous pdf using category sample as index
		RooSimultaneous simPdf("simPdf", "simultaneous pdf", rCatChi);

		// Associate model with the physics state and model_ctl with the control state
		simPdf.addPdf(model, "ChicOne");
		simPdf.addPdf(model_ctl, "ChicTwo");



		// P e r f o r m   a   s i m u l t a n e o u s   f i t
		// ---------------------------------------------------

		// Perform simultaneous fit of model to data and model_ctl to data_ctl
		simPdf.fitTo(combData);



		// P l o t   m o d e l   s l i c e s   o n   d a t a    s l i c e s 
		// ----------------------------------------------------------------

		// Make a frame for the physics sample
		RooPlot* frame1 = rvmass->frame(Bins(30), Title("Physics sample"));

		// Plot all data tagged as physics sample
		combData.plotOn(frame1, Cut("rCatChi==rCatChi::ChicOne"));

		// Plot "physics" slice of simultaneous pdf. 
		// NBL You _must_ project the sample index category with data using ProjWData 
		// as a RooSimultaneous makes no prediction on the shape in the index category 
		// and can thus not be integrated
		simPdf.plotOn(frame1, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, combData), LineColor(kRed));
		simPdf.plotOn(frame1, Slice(rCatChi, "ChicOne"), Components("px"), ProjWData(rCatChi, combData), LineStyle(kDashed), LineColor(kRed));

		// The same plot for the control sample slice
		RooPlot* frame2 = rvmass->frame(Bins(30), Title("Control sample"));
		combData.plotOn(frame2, Cut("rCatChi==rCatChi::ChicTwo"));
		simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, combData));
		simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), Components("px_ctl"), ProjWData(rCatChi, combData), LineStyle(kDashed));
		simPdf.paramOn(frame2, Layout(0.55));


		TCanvas* c = new TCanvas("rf501_simultaneouspdf", "rf501_simultaneouspdf", 800, 400);
		c->Divide(2);
		c->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
		c->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw();

		c->SaveAs("Test.png");






		




			   
		RooPlot *massframe = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
		massframe->SetTitle("mass");
		//rdsNominal->plotOn(massframe);

		string myPdfName = "nominalPdf";
		//CreateModelPdf(myWs, myPdfName);

		
		//////////////////////////////
		// Construct simultaneous PDF
		////////////////////////////

		//myWs.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		//myWs.factory("RooChebychev::background_chi1(rvmass, a1_chi1[0,-10,10])");
		//myWs.factory("SUM::nominalPdf_chi1(nsig_chi1[50,0,1000000]*chic1, nbkg_chi1[2000,0,100000]*background_chi1)");

		//myWs.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		//myWs.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		//myWs.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		//myWs.factory("RooChebychev::background_chi2(rvmass, a1_chi2[0,-10,10])");
		//myWs.factory("SUM::nominalPdf_chi2(nsig_chi2[50,0,1000000]*chic2, nbkg_chi2[2000,0,100000]*background_chi2)");

		myWs.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		myWs.factory("RooChebychev::background_chi1(rvmass, a1_chi1[0,-10,10])");
		myWs.factory("SUM::nominalPdf_chi1(nsig_chi1[50,0,1000000]*chic1, nbkg_chi1[2000,0,100000]*background_chi1)");


		myWs.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		myWs.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		myWs.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		myWs.factory("RooChebychev::background_chi2(rvmass, a1_chi2[0,-10,10])");
		myWs.factory("SUM::nominalPdf_chi2(nsig_chi2[50,0,1000000]*chic2, nbkg_chi2[2000,0,100000]*background_chi2)");

		// C r e a t e   m o d e l   f o r   chic1
// -------------------------------------------------------------

// Create observables
	//RooRealVar x("x", "x", -8, 8);

	// Construct signal pdf
		RooRealVar meanO("meanO", "meanO", 3.5, 3, 4);
		RooRealVar sigmaO("sigmaO", "sigmaO", 0.03, 0.01, 0.1);
		RooGaussian gxO("gxO", "gxO", *rvmass, meanO, sigmaO);

		// Construct background pdf
		RooRealVar a0O("a0O", "a0O", -0.1, -1, 1);
		RooRealVar a1O("a1O", "a1", 0.004, -1, 1);
		RooChebychev pxO("pxO", "pxO", *rvmass, RooArgSet(a0O, a1O));

		// Construct composite pdf
		RooRealVar fO("fO", "fO", 0.2, 0., 1.);
		RooAddPdf modelO("modelO", "modelO", RooArgList(gxO, pxO), fO);



		// C r e a t e   m o d e l   f o r   chic2
		// --------------------------------------------------------------

		// Construct signal pdf. 
		// NOTE that sigma is shared with the signal sample model
		RooRealVar mean_ctlO("mean_ctlO", "mean_ctlO", 3.5, 3, 4);
		RooGaussian gx_ctlO("gx_ctlO", "gx_ctlO", *rvmass, mean_ctlO, sigmaO);

		// Construct the background pdf
		RooRealVar a0_ctlO("a0_ctlO", "a0_ctlO", -0.1, -1, 1);
		RooRealVar a1_ctlO("a1_ctlO", "a1_ctlO", 0.5, -0.1, 1);
		RooChebychev px_ctlO("px_ctlO", "px_ctlO", *rvmass, RooArgSet(a0_ctlO, a1_ctlO));

		// Construct the composite model
		RooRealVar f_ctlO("f_ctlO", "f_ctlO", 0.5, 0., 1.);
		RooAddPdf model_ctlO("model_ctlO", "model_ctlO", RooArgList(gx_ctlO, px_ctlO), f_ctlO);





		cout << "Model done" << endl;


		RooSimultaneous simPdfO("simPdfO", "simultaneous pdf", rCatChi);

		// Associate model with the physics state and model_ctl with the control state
		simPdfO.addPdf(*myWs.pdf("nominalPdf_chi1"), "ChicOne");
		simPdfO.addPdf(*myWs.pdf("nominalPdf_chi2"), "ChicTwo");

		//myWs.import(simPdfO);


		// P e r f o r m   a   s i m u l t a n e o u s   f i t
		// ---------------------------------------------------

		// Perform simultaneous fit of model to data and model_ctl to data_ctl
		simPdfO.fitTo(*rdsNominal_chicBoth);
		//myWs.pdf("simPdfO")->fitTo(*rdsNominal_chicBoth, Extended(true));

		rdsNominal_chicBoth->plotOn(massframe, Cut("rCatChi==rCatChi::ChicOne"), MarkerStyle(20), MarkerColor(kRed));
		rdsNominal_chicBoth->plotOn(massframe, Cut("rCatChi==rCatChi::ChicTwo"), MarkerStyle(20), MarkerColor(kGreen));

		simPdfO.plotOn(massframe, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kRed));
		simPdfO.plotOn(massframe, Slice(rCatChi, "ChicOne"), Components("chic1"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kRed));


		simPdfO.plotOn(massframe, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kGreen));
		simPdfO.plotOn(massframe, Slice(rCatChi, "ChicTwo"), Components("chic2"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kGreen));

		simPdfO.paramOn(massframe, Layout(0.55));

		TCanvas *cTest = new TCanvas("cTest", "cTest", 1000, 600);
		gPad->SetLeftMargin(0.15);
		massframe->GetYaxis()->SetTitleOffset(1.3);
		massframe->Draw();


		cTest->SaveAs("CanvasCTest_RW3.png");









		//////////////////////////////
		//Construct simultaneous PDF
		////////////////////////////

		//myWs.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		//myWs.factory("RooChebychev::background_chi1(rvmass, a1_chi1[0,-10,10])");
		//myWs.factory("SUM::nominalPdf_chi1(nsig_chi1[50,0,1000000]*chic1, nbkg_chi1[2000,0,100000]*background_chi1)");

		//myWs.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		//myWs.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		//myWs.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		//myWs.factory("RooChebychev::background_chi2(rvmass, a1_chi2[0,-10,10])");
		//myWs.factory("SUM::nominalPdf_chi2(nsig_chi2[50,0,1000000]*chic2, nbkg_chi2[2000,0,100000]*background_chi2)");


		
		RooSimultaneous* simPdfD = new RooSimultaneous("simPdfD", "simPdfD", rCatChi);
		simPdfD->addPdf(*myWs.pdf("nominalPdf_chi1"), "ChicOne");
		simPdfD->addPdf(*myWs.pdf("nominalPdf_chi2"), "ChicTwo");


		myWs.import(*simPdfD);

		cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
		//RooFitResult* fitResSim = myWs.pdf("simPdfD")->fitTo(*rdsNominal_chicBoth , Extended(true), SumW2Error(true), NumCPU(1), Save(true));
		myWs.pdf("simPdfD")->fitTo(*rdsNominal_chicBoth, Extended(true), SumW2Error(true), NumCPU(1), Save(true));
		cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;
		//cout << endl << "Importing fit result..." << endl;
		//myWs.import(*fitResSim);

		RooPlot *massframe2 = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
		massframe2->SetTitle("mass");


		//myWs.pdf("simPdf")->plotOn(massframe);
		//rdsNominal_chicOne->plotOn(massframe, MarkerStyle(20), MarkerColor(kRed));
		//rdsNominal_chicTwo->plotOn(massframe, MarkerStyle(20), MarkerColor(kGreen));
		//rdsNominal_chicBoth->plotOn(massframe, MarkerStyle(20), MarkerColor(kBlue));
		rdsNominal_chicBoth->plotOn(massframe2, Cut("rCatChi==rCatChi::ChicOne"), MarkerStyle(20), MarkerColor(kRed));
		rdsNominal_chicBoth->plotOn(massframe2, Cut("rCatChi==rCatChi::ChicTwo"), MarkerStyle(20), MarkerColor(kGreen));
		//myWs.pdf("simPdf")->paramOn(massframe, Layout(0.55));
		//myWs.pdf("simPdf")->plotOn(massframe);//, Components("background"), LineStyle(kDashed));
	//	myWs.pdf("nominalPdf_chi1")->plotOn(massframe, LineColor(kRed));
		//myWs.pdf("nominalPdf_chi2")->plotOn(massframe, LineColor(kGreen));
		//myWs.pdf("nominalPdf_chi1")->plotOn(massframe, Components("chic1"), Range(mass_windowFit_l, mass_windowFit_h), LineStyle(kDashed), LineColor(kRed));
		//myWs.pdf("nominalPdf_chi2")->plotOn(massframe, Components("chic2"), Range(mass_windowFit_l, mass_windowFit_h), LineStyle(kDashed), LineColor(kGreen));


		simPdfD->plotOn(massframe2, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kRed));
		simPdfD->plotOn(massframe2, Slice(rCatChi, "ChicOne"), Components("chic1"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kRed));


		simPdfD->plotOn(massframe2, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kGreen));
		simPdfD->plotOn(massframe2, Slice(rCatChi, "ChicTwo"), Components("chic2"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kGreen));

		simPdfD->paramOn(massframe2, Layout(0.55));


		//simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, combData));
		//simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), Components("px_ctl"), ProjWData(rCatChi, combData), LineStyle(kDashed));
		//simPdf.paramOn(frame2, Layout(0.55));

		//myWs.pdf("simPdfD")->plotOn(massframe, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kRed));
		//myWs.pdf("simPdfD")->plotOn(massframe, Slice(rCatChi, "ChicOne"), Components("chic2"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kRed));
		//myWs.pdf("simPdf")->plotOn(massframe, Slice(rCatChi, "ChicTwo"), Components("chicdgsdgs"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kRed));

		//myWs.pdf("simPdf")->plotOn(massframe, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kGreen));
		//myWs.pdf("simPdfO")->paramOn(massframe, Layout(0.55));

		TCanvas *cTest2 = new TCanvas("cTest2", "cTest2", 1000, 600);
		gPad->SetLeftMargin(0.15);
		massframe2->GetYaxis()->SetTitleOffset(1.3);
		massframe2->Draw();


		cTest2->SaveAs("CanvasCTest2_RW3.png");
		//*/
		
		
		// pT fitting
		/*
		for (int i = 0; i < nbins_pT; i++) {
			RooPlot *massframeBin = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
			massframeBin->SetTitle("mass");
			cout << bins_pT[i] << endl;
			TString TstrCut = TString::Format("rvpt > %f", bins_pT[i]) + " && " + TString::Format("rvpt < %f", bins_pT[i + 1]) + " && rvrap>-1 && rvrap <1";
			cout << TstrCut << endl;
			string strCut = TstrCut.Data();
			RooDataSet* rdsDataBin = (RooDataSet*)myWs.data("rdsNominal")->reduce(strCut.c_str());
			rdsDataBin->plotOn(massframeBin);


			RooFitResult* fitResultBin = myWs.pdf(myPdfName.c_str())->fitTo(*rdsDataBin, SumW2Error(true));// , Extended(true), SumW2Error(true), Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true));
			//fitResultBin->Print("v");
			//myWs.import(*fitResultBin, TString::Format("fitResult_%s", myPdfNameBin.c_str()));
			cout << "signal " << myWs.var("nsig")->getValV() << endl;
			gAsChic_pT->SetPoint(i, ((bins_pT[i] + bins_pT[i + 1]) / 2.0), myWs.var("nsig")->getValV()); //placing the point in the middle of the bin
			gAsChic_pT->SetPointEYhigh(i, myWs.var("nsig")->getErrorHi());
			gAsChic_pT->SetPointEYlow(i, -(myWs.var("nsig")->getErrorLo()));

			myWs.pdf(myPdfName.c_str())->plotOn(massframeBin);
			myWs.pdf(myPdfName.c_str())->paramOn(massframeBin, Layout(0.55));
			myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("background"), LineStyle(kDashed));
			myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic1"), LineStyle(kDashed), LineColor(kRed));
			myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic2"), LineStyle(kDashed), LineColor(kGreen));


			massframeBin->Draw();
			cTest->SaveAs(Form("CanvasCTest_RW3_pT_%.1f_%.1f.png", bins_pT[i], bins_pT[i + 1]));

		}*/

		//FitRooDataSet(gAsChic_pT, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "midrap", intConstrainedFit);
		//FitRooDataSet(gAsChic_y, bins_y, nbins_y, rvmass, myWs, false, "rvrap", myPdfName.c_str(), " && rvpt>6 && rvpt <30", "all", intConstrainedFit);
		//FitRooDataSet(gAsChic_nTrk, bins_nTrk, nbins_nTrk, rvmass, myWs, false, "rvntrack", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "midrap", intConstrainedFit);


	}

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

	if (flagRunFits == true)
	{
		gAsChic_pT->Write();
		gAsJpsi_pT->Write();
		gAsRatio_pT->Write();
		//can_pT->Write();
		gAsChic_y->Write();
		gAsJpsi_y->Write();
		gAsRatio_y->Write();
		//can_y->Write();
		gAsChic_nTrk->Write();
		gAsJpsi_nTrk->Write();
		gAsRatio_nTrk->Write();
	}

	fout->Close();

	if (intConstrainedFit == 2) {
		TFile* fileConstr = new TFile(fileConstraints, "RECREATE");
		gAsChic_alpha_pT->Write();
		gAsChic_alpha_y->Write();
		gAsChic_alpha_nTrk->Write();
		gAsChic_n_pT->Write();
		gAsChic_n_y->Write();
		gAsChic_n_nTrk->Write();
		gAsChic_sigmaRat_pT->Write();
		gAsChic_sigmaRat_y->Write();
		gAsChic_sigmaRat_nTrk->Write();
		gAsChic_sigma1_pT->Write();
		gAsChic_sigma1_y->Write();
		gAsChic_sigma1_nTrk->Write();
		gAsChic_mean1_pT->Write();
		gAsChic_mean1_y->Write();
		gAsChic_mean1_nTrk->Write();
		fileConstr->Close();
	}


	file_log.close();

}


