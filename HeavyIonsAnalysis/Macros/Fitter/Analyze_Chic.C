
// Macro to analyze the chi_c. Starting with the event tree


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
#include "RooMsgService.h"

#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../ChiTreeInit.C"

//#define UsesTMVA // if defined, the code needs to be run in CMSSW_10_3_X, otherwise CMSSW_8_0_X is good enough (production release)

#ifdef UsesTMVA
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif



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
//double bins_y[] = { -2.4, -1.6, -1.0, 1.0, 1.6, 2.4 };
int  nbins_y = sizeof(bins_y) / sizeof(double) - 1;

double bins_nTrk[] = { 0, 50, 100, 150, 200, 300, 400 };
int  nbins_nTrk = sizeof(bins_nTrk) / sizeof(double) - 1;

// constraints

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


bool CreateModelPdf(RooWorkspace& Ws, string pdfName, bool bConstrainedFit = false)
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
		Ws.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		Ws.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		Ws.factory("SUM::signal(c2toc1[0.4,0.1,1.]*chic2, chic1)");
		//Ws.factory("Exponential::expbkg(rvmass,e_1[-0.2,-0.5,0])");
		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("RooChebychev::background(rvmass, a1[0,-10,10])");
		
		if (bConstrainedFit == true) {

			Ws.var("alpha")->setConstant(true);
			Ws.var("n")->setConstant(true);
			Ws.var("sigmaRatio")->setConstant(true);
			Ws.var("sigma1")->setConstant(true);
			Ws.var("mean1")->setConstant(true);
		}
		//Ws.var("mean1")->setConstant(true);
		//Ws.var("mean2")->setConstant(true);
		//Ws.var("erf_offset")->setConstant(true);
		//Ws.var("erf_sigma")->setConstant(true);
		//Ws.factory("")
		Ws.factory("SUM::nominalPdf(nsig[50,0,10000000]*signal, nbkg[2000,0,10000000]*background)");
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



bool RefreshModel(RooWorkspace& Ws, string pdfName, bool isJpsi) // attempt to prevent the fits in the differential bins to get stuck in a weird state when they stop converge properly
{
	if (isJpsi == true) {
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

	}
	else
	{
		Ws.var("c2toc1")->removeError();
		Ws.var("sigma1")->removeError();
		Ws.var("a1")->removeError();
		Ws.var("nsig")->removeError();
		Ws.var("nbkg")->removeError();

		Ws.var("c2toc1")->setVal(0.4);
		Ws.var("sigma1")->setVal(0.005);
		Ws.var("a1")->setVal(0);
		Ws.var("nsig")->setVal(50);
		Ws.var("nbkg")->setVal(2000);
	}




	return true;
}

int SetConstraints(RooWorkspace& Ws, double* bins, int bin, string binVarName = "", string pdfName = "", const char* fileConstraints = "") {
	TFile* fileConstr = new TFile(fileConstraints, "READ"); //opening the file for every fit is probably not the most efficient, but fitting itself is much more time consuming, so it is okay
	if (fileConstr->IsOpen() == false) {
		cout << "Constraint file not properly opened!" << endl << endl;
		return -1;
	}
	string binNaming;
	if (binVarName.compare("rvpt") == 0) {
		binNaming = "pt";
	}
	else if (binVarName.compare("rvrap") == 0) {
		binNaming = "y";
	}
	else if (binVarName.compare("rvntrack") == 0) {
		binNaming = "nTrack";
		if (bin > 2) bin = 2; //WARNING - for now bin 2 used for higher ntracks bins (those constraints are bad due to stats)
	}
	if (pdfName.compare("nominalPdf") == 0)
	{
		TGraphAsymmErrors* gAsChic_sigma1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_sigma1_" + binNaming).c_str());
		TGraphAsymmErrors* gAsChic_mean1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_mean1_" + binNaming).c_str());
		TGraphAsymmErrors* gAsChic_alpha = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_alpha_" + binNaming).c_str());
		TGraphAsymmErrors* gAsChic_n = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_n_" + binNaming).c_str());
		TGraphAsymmErrors* gAsChic_sigmaRatio = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_sigmaRatio_" + binNaming).c_str());		

		double xVal, yVal;
		gAsChic_sigma1->GetPoint(bin, xVal, yVal);
		if (std::abs(xVal- (bins[bin] + bins[bin + 1]) / 2.0)>0.01) cout <<endl<< "WARNING: x values of binning between fit and constraints don't agree! " << xVal << "   " << (bins[bin] + bins[bin + 1]) / 2.0 << endl<<endl;
		cout << "sigma1 " << yVal << endl;
		Ws.var("sigma1")->setVal(yVal);
		gAsChic_mean1->GetPoint(bin, xVal, yVal);
		Ws.var("mean1")->setVal(yVal);
		gAsChic_alpha->GetPoint(bin, xVal, yVal);
		Ws.var("alpha")->setVal(yVal);
		gAsChic_n->GetPoint(bin, xVal, yVal);
		Ws.var("n")->setVal(yVal);
		gAsChic_sigmaRatio->GetPoint(bin, xVal, yVal);
		Ws.var("sigmaRatio")->setVal(yVal);

		Ws.var("alpha")->setConstant(true);
		Ws.var("n")->setConstant(true);
		Ws.var("sigmaRatio")->setConstant(true);
		Ws.var("sigma1")->setConstant(true);
		Ws.var("mean1")->setConstant(true);



	
	}


	fileConstr->Close();
	return 0;
}



int FitRooDataSet(TGraphAsymmErrors* gAsResult, double* bins, int nbins, RooRealVar* rvmass, RooWorkspace& Ws, bool isJpsi = false, string binVarName = "", string myPdfName = "", string extraCut = "", string canvasName = "", bool bConstrainedFit = false, const char* fileConstraints = "") { //uses global variables
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
		if (isJpsi == false) { rdsDataBin = (RooDataSet*)Ws.data("rdsNominal")->reduce(strCut.c_str()); 
			RefreshModel(Ws, myPdfName, false);
		}else {
			rdsDataBin = (RooDataSet*)Ws.data("rdsNominalJpsi")->reduce(strCut.c_str()); 
			RefreshModel(Ws, myPdfName, true);
		}

		rdsDataBin->plotOn(massframeBin);

		if (bConstrainedFit == true) {
			SetConstraints(Ws, bins, i, binVarName, myPdfName, fileConstraints);
		}


		RooFitResult* fitResultBin = Ws.pdf(myPdfName.c_str())->fitTo(*rdsDataBin, Extended(true), SumW2Error(true), NumCPU(1), PrintLevel(-1), Save(true));
		//fitResultBin->Print("v");
		//Ws.import(*fitResultBin, TString::Format("fitResult_%s", myPdfNameBin.c_str()));
		if (isJpsi == false) {
			cout << "signal " << Ws.var("nsig")->getValV() << endl;
			SetPointFromFit(gAsResult, "nsig", Ws, i, bins);

			//if (intConstrainedFit == 2) {
			//	if (myPdfName.compare("nominalPdf") == 0) {
			//		if (binVarName.compare("rvpt") == 0){
			//			SetPointFromFit(gAsChic_alpha_pT, "alpha", Ws, i, bins);
			//			SetPointFromFit(gAsChic_n_pT, "n", Ws, i, bins);
			//			SetPointFromFit(gAsChic_sigmaRat_pT, "sigmaRatio", Ws, i, bins);
			//			SetPointFromFit(gAsChic_sigma1_pT, "sigma1", Ws, i, bins);
			//			SetPointFromFit(gAsChic_mean1_pT, "mean1", Ws, i, bins);

			//		} 
			//		else if (binVarName.compare("rvrap") == 0) {
			//			SetPointFromFit(gAsChic_alpha_y, "alpha", Ws, i, bins);
			//			SetPointFromFit(gAsChic_n_y, "n", Ws, i, bins);
			//			SetPointFromFit(gAsChic_sigmaRat_y, "sigmaRatio", Ws, i, bins);
			//			SetPointFromFit(gAsChic_sigma1_y, "sigma1", Ws, i, bins);
			//			SetPointFromFit(gAsChic_mean1_y, "mean1", Ws, i, bins);

			//		}
			//		else if (binVarName.compare("rvntrack") == 0) {
			//			SetPointFromFit(gAsChic_alpha_nTrk, "alpha", Ws, i, bins);
			//			SetPointFromFit(gAsChic_n_nTrk, "n", Ws, i, bins);
			//			SetPointFromFit(gAsChic_sigmaRat_nTrk, "sigmaRatio", Ws, i, bins);
			//			SetPointFromFit(gAsChic_sigma1_nTrk, "sigma1", Ws, i, bins);
			//			SetPointFromFit(gAsChic_mean1_nTrk, "mean1", Ws, i, bins);

			//		}
			//		else cout << "No match for the variable" << endl;

			//	}
			//	else cout << "The pdf not matching" << endl;
			//
			//}

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
		cFit->SaveAs(((string)"FitterOutput/FitResult_" + (isJpsi?"Jpsi_":"Chic_") + canvasName + "_" + binVarName + Form("_%ibin_", i) + Form("_%.1f_%.1f.png", bins[i], bins[i + 1])).c_str());

	}

	delete cFit;


	return 0;
}



///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////

//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true, const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV-bothDirRW4.root", const char* fileOut = "Chi_c_output_RW5_ComparisonAlberto2M.root", const char* fileRds = "rds_save_test.root", bool isMC = false, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true, const char* fileIn = "/afs/cern.ch/user/o/okukral/Chic_pPb/CMSSW_8_0_30/src/HeavyIonsAnalysis/ChiAnalysis/test/Chi_c_CompOneRun_RW6b.root", const char* fileOut = "Chi_c_output_RW6_ComparisonAlberto285718.root", const char* fileRds = "rds_save_test.root", bool isMC = false, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true,  const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV-bothDirRW4.root", const char* fileOut = "Chi_c_output_RW4_testCutTightRefitAlberto.root", const char* fileRds = "rds_RW4_Full_noWeights_CutTightRefitAlberto.root", bool isMC = false, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true,  const char* fileIn = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-MC8_BothDir.root", const char* fileOut = "Chi_c_output_MC8_test.root", const char* fileRds = "rds_MC8_test.root", bool isMC = true, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
{
	//gStyle->SetOptStat(1111);
	//gStyle->SetOptStat(0);
	setTDRStyle();
	

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

	TH1D* hphoton_MVA_response = new TH1D("hphoton_MVA_response", "", 1000, -2, 2);

		// Various
	TH1D* hntracks_inEvent = new TH1D("hntracks_inEvent", "", 400, 0, 400);
	TH1D* hconv_deltaEta = new TH1D("conv_deltaEta", "", 1000, 0, 1);
	TH2D* hconv_deltapT_Eta = new TH2D("conv_deltapT_Eta", "", 200, 0, 0.02, 200, 0, 0.2);
	TH1D* hconv_m = new TH1D("conv_m", "", 1000, 0, 0.5);







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

	RooRealVar* rvmassJpsi = new RooRealVar("rvmassJpsi", "Inv mass #mu#mu", mass_windowJpsi_l, mass_windowJpsi_h, "GeV/c^{2}");
	RooRealVar* rvptJpsi = new RooRealVar("rvptJpsi", "#mu#mu p_{T}", 0.0, 50.0, "GeV/c");
	RooRealVar* rvrapJpsi = new RooRealVar("rvrapJpsi", "#mu#mu y", -2.4, 2.4, "");
	RooRealVar* rvntrackJpsi = new RooRealVar("rvntrackJpsi", "ntrack", 0.0, 500.0, "");
	RooRealVar* rvweightJpsi = new RooRealVar("rvweightJpsi", "weight", 0.0, 10000.0, "");
	RooArgSet*  colsJpsi = new RooArgSet(*rvmassJpsi, *rvptJpsi, *rvrapJpsi, *rvntrackJpsi, *rvweightJpsi);
	RooDataSet* rdsNominalJpsi = new RooDataSet("rdsNominalJpsi", "rdsNominalJpsi", *colsJpsi, WeightVar(*rvweightJpsi), StoreAsymError(*rvmassJpsi));
	


	if (flagGenerateRds == true) {


		file_log << "index runNumber  eventNumber  (chiCandPerEvent-ignore)  rap_chi  pT_chi  Mass(mchic - mJpsi +3.097) muon1_eta muon1_pt  muon2_eta muon2_pt (conv position) conv_eta conv_pt " << endl;
		//TMVA 
		#ifdef UsesTMVA


		TMVA::Tools::Instance();
		TMVA::Reader *TMWAreader = new TMVA::Reader("!Color:!Silent");

		float convQuality_isHighPurityValue, convQuality_isGeneralTracksOnlyValue, conv_tkVtxCompatibilityOKValue, conv_compatibleInnerHitsOKValue, conv_tk1NumOfDOFValue, conv_tk2NumOfDOFValue; //all variables need to be floats
		float conv_vertexPositionRhoValue, conv_sigmaTkVtx1Value, conv_sigmaTkVtx2Value, conv_vertexChi2ProbValue ;
		float conv_dzToClosestPriVtxValue, conv_dxyPriVtxTimesCharge_Tr1Value, conv_dxyPriVtxTimesCharge_Tr2Value, conv_minDistanceOfApproachValue, conv_etaValue, conv_ptValue;

		TMWAreader->AddVariable("convQuality_isHighPurity", &convQuality_isHighPurityValue);
		TMWAreader->AddVariable("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnlyValue);
		TMWAreader->AddSpectator("conv_vertexPositionRho", &conv_vertexPositionRhoValue);
		TMWAreader->AddVariable("conv_sigmaTkVtx1", &conv_sigmaTkVtx1Value);
		TMWAreader->AddVariable("conv_sigmaTkVtx2", &conv_sigmaTkVtx2Value);
		TMWAreader->AddVariable("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOKValue);
		TMWAreader->AddVariable("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOKValue);
		TMWAreader->AddVariable("conv_vertexChi2Prob", &conv_vertexChi2ProbValue);
		TMWAreader->AddVariable("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtxValue);
		TMWAreader->AddVariable("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1Value);
		TMWAreader->AddVariable("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2Value);
		TMWAreader->AddVariable("conv_tk1NumOfDOF", &conv_tk1NumOfDOFValue);
		TMWAreader->AddVariable("conv_tk2NumOfDOF", &conv_tk2NumOfDOFValue);
		TMWAreader->AddVariable("conv_minDistanceOfApproach", &conv_minDistanceOfApproachValue);
		TMWAreader->AddSpectator("conv_eta", &conv_etaValue);
		TMWAreader->AddSpectator("conv_pt", &conv_ptValue);

		TMWAreader->BookMVA("BDT", "BDTWeight/TMVAClassification_BDT2.weights.xml");
		#endif


		//load corrections
		TFile* fCor = new TFile(fileCorrection, "READ");

		TH2D* hWeightChic = (TH2D*)fCor->Get("h_chiAccEff_rat");
		TH2D* hWeightJpsi = (TH2D*)fCor->Get("h_JpsiAccEff_rat");

		int weird_decay_counter = 0;
		long nchicCounter = 0, nchicCounterPass = 0, nchicCounterPassMass = 0, muon1ptCounter = 0, muon2ptCounter = 0;
		long nConvCounter = 0, nConvSplitCounter = 0, nConvCounterPeak = 0, nConvSplitCounterPeak = 0;
		bool passDimSel = false;
		bool passDimSelTight = false;
		Long64_t nentries = event_tree->GetEntries();
		if (nentries > 2000000) { nentries = 2000000; }
		cout << nentries << endl;
		for (Long64_t i = 0; i < nentries; i++)
		{

			event_tree->GetEntry(i);
			if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }
			if (isMC == true) {
				if (gen_pdgId->at(0) != PythCode_chic1 && gen_pdgId->at(0) != PythCode_chic2) { continue; } //remove chic0 from MC, expecting one gen chic per event
				if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
				{
					weird_decay_counter++;
					continue;
				}
			}


			if (eventNumber != 180043931 && eventNumber != 62726391 && eventNumber != 2372749 && eventNumber != 92140595 && eventNumber != 144554927 && eventNumber != 80353502 && eventNumber != 128861276 && eventNumber != 103338128)
			{
				continue;
			}

			int nRefitNumber=-1; //needed for RW4 only, fix to non-ideal storing of info

			//if (isMC == true && (gen_pdgId->at(0) != PythCode_chic1)) { continue; } // only chic1
			//if (isMC == true && (gen_pdgId->at(0) != PythCode_chic2)) { continue; } // only chic2

			////////////////////////
			//// V A R I O U S  ////
			///////////////////////
			hntracks_inEvent->Fill(ntracks_inEvent);
			//cout << ntracks_inEvent << endl;

			// check the conversions - if they are many identical ones
			int nConv = conv_eta->size();
			if (nConv > 1) {
				//cout << "nConv: " << nConv << endl;
				for (int iconv = 0; iconv < nConv; iconv++)
				{
					//cout << conv_eta->at(iconv) << endl;
					for (int iconv2 = 0; iconv2 < iconv; iconv2++)
					{
						double deltaEta = abs(conv_eta->at(iconv) - conv_eta->at(iconv2));
						hconv_deltaEta->Fill(deltaEta);
						double deltapT = abs(conv_pt->at(iconv) - conv_pt->at(iconv2));
						hconv_deltapT_Eta->Fill(deltaEta, deltapT);

						//if (conv_vertexPositionRho->at(iconv) < 10 && conv_vertexPositionRho->at(iconv2) < 10) {
						//}

							nConvCounter++;
							TVector3* TVconv1 = (TVector3*)conv_vtx->At(iconv);
							TVector3* TVconv2 = (TVector3*)conv_vtx->At(iconv2);
							double convVtx_deltaX = abs(TVconv1->X() - TVconv2->X());
							double convVtx_deltaY = abs(TVconv1->Y() - TVconv2->Y());
							double convVtx_deltaZ = abs(TVconv1->Z() - TVconv2->Z());

							if (deltaEta < 0.02 && deltapT < 0.2 && convVtx_deltaX<2.0 && convVtx_deltaY < 2.0 && convVtx_deltaZ < 2.0) {
								TLorentzVector* LVconv1 = (TLorentzVector*)conv_p4->At(iconv);
								TLorentzVector* LVconv2 = (TLorentzVector*)conv_p4->At(iconv2);
								TLorentzVector LVconvTot = *LVconv1;
								LVconvTot += *LVconv2;
								hconv_m->Fill(LVconvTot.M());
								nConvSplitCounter++;
							}
						
					}
				}

			}








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
				if (MuonAcceptanceTight(muon_eta->at(muon1Pos), muon_pt->at(muon1Pos)) == false) continue;
				if (MuonAcceptanceTight(muon_eta->at(muon2Pos), muon_pt->at(muon2Pos)) == false) continue;
				if (MuonSelectionPass(muon1Pos) == false) continue;
				if (MuonSelectionPass(muon2Pos) == false) continue;
				// photon
				int convPos = chi_daughterConv_position->at(iChi);
				if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) continue;
				
				//cuts or MVA for photons
				#ifdef UsesTMVA
				convQuality_isHighPurityValue = convQuality_isHighPurity->at(convPos);
				convQuality_isGeneralTracksOnlyValue = convQuality_isGeneralTracksOnly->at(convPos);
				conv_vertexPositionRhoValue = conv_vertexPositionRho->at(convPos);
				conv_sigmaTkVtx1Value = conv_sigmaTkVtx1->at(convPos);
				conv_sigmaTkVtx2Value = conv_sigmaTkVtx2->at(convPos);
				conv_tkVtxCompatibilityOKValue = conv_tkVtxCompatibilityOK->at(convPos);
				conv_compatibleInnerHitsOKValue = conv_compatibleInnerHitsOK->at(convPos);
				conv_vertexChi2ProbValue = conv_vertexChi2Prob->at(convPos);
				conv_dzToClosestPriVtxValue = conv_dzToClosestPriVtx->at(convPos);
				conv_dxyPriVtxTimesCharge_Tr1Value = conv_dxyPriVtxTimesCharge_Tr1->at(convPos);
				conv_dxyPriVtxTimesCharge_Tr2Value = conv_dxyPriVtxTimesCharge_Tr2->at(convPos);
				conv_tk1NumOfDOFValue = conv_tk1NumOfDOF->at(convPos);
				conv_tk2NumOfDOFValue = conv_tk2NumOfDOF->at(convPos);
				conv_minDistanceOfApproachValue = conv_minDistanceOfApproach->at(convPos);
				conv_etaValue = conv_eta->at(convPos);
				conv_ptValue = conv_pt->at(convPos);

				
				double photMVA = TMWAreader->EvaluateMVA("BDT");
				hphoton_MVA_response->Fill(photMVA);
				//cout << "Values MVA: HP: " << convQuality_isHighPurityValue << "   TracksOnly: " << convQuality_isGeneralTracksOnlyValue << "    Vprob: " << conv_vertexChi2ProbValue << "    and the response: "<< photMVA << endl;
				//if (photMVA < -0.2) continue;
				#endif


				if (PhotSelectionPassTight(convPos) == false) continue;
				passDimSel = passDimSelTight;

				if (passDimSel == true) { ++nchicCounterPass; } // SelectionsPassed
				// Get Lorentz V
				LVchic = (TLorentzVector*)chi_p4->At(iChi);
				LVdimuon = (TLorentzVector*)dimuon_p4->At(dimuonPos);
				LVconv = (TLorentzVector*)conv_p4->At(convPos);

				double pT_chi = chi_pt->at(iChi);
				double rap_chi = LVchic->Rapidity();
				double m_chi = LVchic->M();
				double refit_vProb=0;
				/////////
				//refit
				if (chi_kinematicRefitFlag->at(iChi) == 1 || chi_kinematicRefitFlag->at(iChi) == 3) { //good refits
					nRefitNumber++; //needed for RW4 only, fix to non-ideal storing of info
					nRefitNumber = iChi; //overwrite previous line

					//cout << chi_refit_vprob->at(nRefitNumber) << endl;
					//cout << chi_refit_ctauPV->at(nRefitNumber) << endl;
					//cout << "RefitStored: " << chi_refitStored->at(nRefitNumber).mass() << endl;
					m_chi = chi_refitStored->at(nRefitNumber).mass();//use refit mass
				    //refit_vProb = chi_refit_vprob->at(nRefitNumber);
					refit_vProb = chi_refit_vprob->at(iChi);
				}
				else continue;//skip those that don't have it
				if (refit_vProb < 0.01) continue;
				////////
							   

				double dimuonM = LVdimuon->M();
				double Mdiff = m_chi - dimuonM + 3.097;// Assume J/psi mass

				//if (passDimSel == true && dimuonM>2.95 && dimuonM<3.2 && Mdiff > mass_window_l && Mdiff < mass_window_h) ++nchicCounterPassMass;
				if (passDimSel == true && dimuonM > 2.9 && dimuonM<3.3 && Mdiff > mass_window_l && Mdiff < mass_window_h) ++nchicCounterPassMass;


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


				//if (dimuonM<2.95 || dimuonM>3.2) continue; //require narrow dimuon mass
				if (dimuonM<2.90 || dimuonM>3.3) continue; //require narrow dimuon mass

				if (passDimSel == true) {
					hchic_M->Fill(m_chi); // just raw M
					hchic_M_rotGamma->Fill(LVchic_rotGamma->M());
				}
				if (dimuon_charge->at(dimuonPos) != 0) {
					hchic_M_SS->Fill(m_chi);
				}


				// nominal

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
							//rvweight->setVal(1 / accEff_chi);
							rvweight->setVal(1);
						}
						else rvweight->setVal(1); 
						if (isMC) {
							rvweight->setVal(1); //no weights if MC
						}

						rdsNominal->add(*cols, rvweight->getVal());

						///////////////////
						// crosschecks
						/////////////////
						
						//TVector3* TVconv1 = (TVector3*)conv_vtx->At(convPos);
						//file_log << nchicCounterPassMass << " " << runNumber << " " << eventNumber << " " << chiCandPerEvent << " " << rap_chi << " " << pT_chi << " " << Mdiff << " m1 " << muon_eta->at(muon1Pos) << " " << muon_pt->at(muon1Pos)<< " m2 " << muon_eta->at(muon2Pos)<< " " << muon_pt->at(muon2Pos) << " conv " << convPos << " " << conv_eta->at(convPos) << " " << conv_pt->at(convPos) << " " << conv_vertexPositionRho->at(convPos) << " "  << TVconv1->X() << " " << TVconv1->Y() << " " << TVconv1->Z() << endl;

						int nRefitNumber2 = -1;

						// check the conversions - if they are many identical ones
						for (int iChi2 = 0; iChi2 < iChi; iChi2++)
						{
							int dimuonPos2 = chi_daughterJpsi_position->at(iChi2);
							double passDimSel2 = DimuonSelectionPass(dimuonPos2);
							double passDimSelTight2 = DimuonSelectionPassTight(dimuonPos2, ((TLorentzVector*)dimuon_p4->At(dimuonPos2))->Rapidity());
							// muon cuts
							int muon1Pos2 = dimuon_muon1_position->at(dimuonPos2);
							int muon2Pos2 = dimuon_muon2_position->at(dimuonPos2);
							if (MuonAcceptanceTight(muon_eta->at(muon1Pos2), muon_pt->at(muon1Pos2)) == false) continue;
							if (MuonAcceptanceTight(muon_eta->at(muon2Pos2), muon_pt->at(muon2Pos2)) == false) continue;
							if (MuonSelectionPass(muon1Pos2) == false) continue;
							if (MuonSelectionPass(muon2Pos2) == false) continue;
							// photon
							int convPos2 = chi_daughterConv_position->at(iChi2);
							if (PhotAcceptance(conv_eta->at(convPos2), conv_pt->at(convPos2)) == false) continue;
							if (PhotSelectionPassTight(convPos2) == false) continue;
							passDimSel2 = passDimSelTight2;

							TLorentzVector* LVchic2 = (TLorentzVector*)chi_p4->At(iChi2);
							TLorentzVector* LVdimuon2 = (TLorentzVector*)dimuon_p4->At(dimuonPos2);

							double pT_chi2 = chi_pt->at(iChi2);
							double rap_chi2 = LVchic2->Rapidity();
							double m_chi2 = LVchic2->M();
							double refit_vProb2 = 0;
							/////////
							//refit
							if (chi_kinematicRefitFlag->at(iChi2) == 1 || chi_kinematicRefitFlag->at(iChi2) == 3) { //good refits
								nRefitNumber2++; //needed for RW4 only, fix to non-ideal storing of info
								nRefitNumber2 = iChi2; //overwrite previous line

								//cout << chi_refit_vprob->at(nRefitNumber) << endl;
								//cout << chi_refit_ctauPV->at(nRefitNumber) << endl;
								//cout << "RefitStored: " << chi_refitStored->at(nRefitNumber).mass() << endl;
								m_chi2 = chi_refitStored->at(nRefitNumber2).mass();//use refit mass
								//refit_vProb = chi_refit_vprob->at(nRefitNumber);
								refit_vProb2 = chi_refit_vprob->at(iChi2);
							}
							else continue;//skip those that don't have it
							if (refit_vProb2 < 0.01) continue;
							////////


							double dimuonM2 = LVdimuon2->M();
							double Mdiff2 = m_chi2 - dimuonM2 + 3.097;// Assume J/psi mass
							if (passDimSel2 == true && dimuonM2 > 2.9 && dimuonM2<3.3 && Mdiff2 > mass_window_l && Mdiff < mass_window_h)

							{
									double deltaEta = abs(conv_eta->at(iChi) - conv_eta->at(iChi2));
									double deltapT = abs(conv_pt->at(iChi) - conv_pt->at(iChi2));


									nConvCounterPeak++;
									TVector3* TVconv1 = (TVector3*)conv_vtx->At(iChi);
									TVector3* TVconv2 = (TVector3*)conv_vtx->At(iChi2);
									double convVtx_deltaX = abs(TVconv1->X() - TVconv2->X());
									double convVtx_deltaY = abs(TVconv1->Y() - TVconv2->Y());
									double convVtx_deltaZ = abs(TVconv1->Z() - TVconv2->Z());

									if (deltaEta < 0.02 && deltapT < 0.2 && convVtx_deltaX < 2.0 && convVtx_deltaY < 2.0 && convVtx_deltaZ < 2.0) {
										//TLorentzVector* LVconv1 = (TLorentzVector*)conv_p4->At(iChi);
										//TLorentzVector* LVconv2 = (TLorentzVector*)conv_p4->At(iChi2);
										//TLorentzVector LVconvTot = *LVconv1;
										//LVconvTot += *LVconv2;
										//hconv_m->Fill(LVconvTot.M());
										nConvSplitCounterPeak++;
										file_log << nchicCounterPassMass << " " << runNumber << " " << eventNumber << " " << chiCandPerEvent << " " << rap_chi << " " << pT_chi << " " << Mdiff << " m1 " << muon_eta->at(muon1Pos) << " " << muon_pt->at(muon1Pos)<< " m2 " << muon_eta->at(muon2Pos)<< " " << muon_pt->at(muon2Pos) << " conv " << convPos << " " << conv_eta->at(convPos) << " " << conv_pt->at(convPos) << " " << conv_vertexPositionRho->at(convPos) << " "  << TVconv1->X() << " " << TVconv1->Y() << " " << TVconv1->Z() << endl;
										file_log << nchicCounterPassMass << " " << runNumber << " " << eventNumber << " " << chiCandPerEvent << " " << rap_chi2 << " " << pT_chi2 << " " << Mdiff2 << " m1 " << muon_eta->at(muon1Pos2) << " " << muon_pt->at(muon1Pos2) << " m2 " << muon_eta->at(muon2Pos2) << " " << muon_pt->at(muon2Pos2) << " conv " << convPos2 << " " << conv_eta->at(convPos2) << " " << conv_pt->at(convPos2) << " " << conv_vertexPositionRho->at(convPos2) << " " << TVconv2->X() << " " << TVconv2->Y() << " " << TVconv2->Z() << endl << endl;

									}

								
							}

						}









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
		cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;
		cout << endl;
		cout << "Nconversions: " << nConvCounter << "  of which " << nConvSplitCounter << " was split, that is " << ((double)nConvSplitCounter / (double)nConvCounter * 100) << " %." <<endl;
		cout << "Nconversions from the candidates: " << nConvCounterPeak << "  of which " << nConvSplitCounterPeak << " was split, that is " << ((double)nConvSplitCounterPeak / (double)nConvCounterPeak * 100) << " %.  Or " << ((double)nConvSplitCounterPeak / (double)nchicCounterPassMass * 100) << " % of total chic" << endl << endl;
		cout << endl;

		// Save rds
		TFile* DBFile = new TFile(fileRds, "RECREATE");
		DBFile->cd();
		rdsNominal->Write();
		rdsNominalJpsi->Write();
		DBFile->Write(); DBFile->Close(); delete DBFile;
		f1->cd();
		myWs.import(*rdsNominal);
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

		DBFile->Close(); delete DBFile;
	}


	if (flagRunFits == true) {


		//////////////////////////////////////////  
		/////////   R  o  o  f  i  t   //////////
		///////////////////////////////////////////

		RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
		//ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
		
		cout << "nEntries: " << rdsNominal->numEntries() << endl;
		rdsNominal = (RooDataSet*)rdsNominal->reduce(mass_windowFit.c_str());

		RooPlot *massframe = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
		massframe->SetTitle("mass");
		rdsNominal->plotOn(massframe);

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
		cout << "signal " << myWs.var("nsig")->getValV() << " + " << myWs.var("nsig")->getErrorHi() << " - " << myWs.var("nsig")->getErrorLo() << endl;

		TCanvas *cTest = new TCanvas("cTest", "cTest", 1000, 600);
		gPad->SetLeftMargin(0.15);
		massframe->GetYaxis()->SetTitleOffset(1.3);
		massframe->Draw();


		cTest->SaveAs("FitterOutput/CanvasCTest_RW3.png");
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

		FitRooDataSet(gAsChic_pT, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "midrap", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_y, bins_y, nbins_y, rvmass, myWs, false, "rvrap", myPdfName.c_str(), " && rvpt>6 && rvpt <30", "all", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_nTrk, bins_nTrk, nbins_nTrk, rvmass, myWs, false, "rvntrack", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "midrap", flagConstrainedFit, fileConstraints);


		
		/////////////////////
		//// J/psi fit  /////
		/////////////////////

		
		rdsNominalJpsi = (RooDataSet*)rdsNominalJpsi->reduce(mass_windowFitJpsi.c_str());
		RooPlot *massframeJpsi = rvmassJpsi->frame(mass_windowFitJpsi_l, mass_windowFitJpsi_h, nMassBinsJpsi);
		massframeJpsi->SetTitle("massJpsi");
		rdsNominalJpsi->plotOn(massframeJpsi);

		string myPdfNameJpsi = "nominalPdfJpsi";
		CreateModelPdf(myWs, myPdfNameJpsi);


		//RooFitResult* fitResultJpsi = myWs.pdf(myPdfNameJpsi.c_str())->fitTo(*myWs.data("rdsNominalJpsi"), Extended(true), SumW2Error(true), Range(mass_windowFitJpsi_l, mass_windowFitJpsi_h), NumCPU(1), Save(true));
		//myWs.import(*fitResultJpsi, TString::Format("fitResultJpsi_%s", myPdfNameJpsi.c_str()));


		myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi);
		myWs.pdf(myPdfNameJpsi.c_str())->paramOn(massframeJpsi, Layout(0.55));
		myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("backgroundJpsi"), LineStyle(kDashed));
		myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("Jpsi"), LineStyle(kDashed), LineColor(kRed));
		myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("psi2"), LineStyle(kDashed), LineColor(kGreen));
		cout << endl << endl << "HEREJpsi" << endl << endl;



		TCanvas *cTestJpsi = new TCanvas("cTestJpsi", "cTestJpsi", 1000, 600);
		gPad->SetLeftMargin(0.15);
		massframeJpsi->GetYaxis()->SetTitleOffset(1.3);
		massframeJpsi->Draw();


		cTestJpsi->SaveAs("FitterOutput/CanvasCTest_RW3_Jpsi.png");



		//FitRooDataSet(gAsJpsi_pT, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && rvrapJpsi>-1 && rvrapJpsi <1", "midrap", false);
		//FitRooDataSet(gAsJpsi_y, bins_y, nbins_y, rvmassJpsi, myWs, true, "rvrapJpsi", myPdfNameJpsi.c_str(), " && rvptJpsi>6 && rvptJpsi <30", "all", false);
		//FitRooDataSet(gAsJpsi_nTrk, bins_nTrk, nbins_nTrk, rvmassJpsi, myWs, true, "rvntrackJpsi", myPdfNameJpsi.c_str(), " && rvrapJpsi>-1 && rvrapJpsi <1", "midrap", false);

		/*
		for (int i = 0; i < nbins_y; i++) {
			RooPlot *massframeJpsiBin = rvmassJpsi->frame(mass_windowFitJpsi_l, mass_windowFitJpsi_h, nMassBinsJpsi);
			massframeJpsiBin->SetTitle("massJpsi");
			file_log << bins_y[i] << endl;
			TString TstrCut = TString::Format("rvrapJpsi > %f", bins_y[i]) + " && " + TString::Format("rvrapJpsi < %f", bins_y[i + 1]) + " && rvptJpsi>6 && rvptJpsi <30";
			file_log << TstrCut << endl;
			string strCut = TstrCut.Data();
			RooDataSet* rdsDataJpsiBin = (RooDataSet*)myWs.data("rdsNominalJpsi")->reduce(strCut.c_str());
			rdsDataJpsiBin->plotOn(massframeJpsiBin);

			//myWs.var("rvpt")->setMin(6);
			//myWs.var("rvpt")->setMax(8);

			////fitting
			//string myPdfNameJpsiBin = (TString::Format("nominalPdfJpsi_%i",i)).Data();
			//string myPdfNameJpsiBin = "nominalPdf";
			//cout << endl << endl << endl << myPdfNameJpsi << endl << endl;
			//CreateModelPdf(myWs, myPdfNameJpsi);

			RefreshModel(myWs, myPdfNameJpsi);

			RooFitResult* fitResultJpsiBin = myWs.pdf(myPdfNameJpsi.c_str())->fitTo(*rdsDataJpsiBin, Extended(true), SumW2Error(true));// , Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true), PrintLevel(-1));
			////fitResultJpsiBin->Print("v");
			////myWs.import(*fitResultJpsiBin, TString::Format("fitResultJpsi_%s", myPdfNameJpsiBin.c_str()));
			//gAsJpsi_y->SetPoint(i, ((bins_y[i] + bins_y[i + 1]) / 2.0), myWs.var("nsigJpsi")->getValV()); //placing the point in the middle of the bin
			//gAsJpsi_y->SetPointEYhigh(i, myWs.var("nsigJpsi")->getErrorHi());
			//gAsJpsi_y->SetPointEYlow(i, -(myWs.var("nsigJpsi")->getErrorLo()));

			//file_log << "signal Jpsi " << myWs.var("nsigJpsi")->getValV() << " + " << myWs.var("nsigJpsi")->getErrorHi() << " - " << myWs.var("nsigJpsi")->getErrorLo() << endl;

			myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsiBin);
			myWs.pdf(myPdfNameJpsi.c_str())->paramOn(massframeJpsiBin, Layout(0.55));
			myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsiBin, Components("backgroundJpsi"), LineStyle(kDashed));
			myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsiBin, Components("Jpsi"), LineStyle(kDashed), LineColor(kRed));
			//myWs.pdf(myPdfNameJpsiBin.c_str())->plotOn(massframeJpsiBin, Components("psi2"), LineStyle(kDashed), LineColor(kGreen));


			massframeJpsiBin->Draw();
			cTestJpsi->SaveAs(Form("CanvasCTest_RW3_Jpsi_y_%.1f_%.1f.png", bins_y[i], bins_y[i + 1]));

		}
		//*/
		

		GetRatio(gAsRatio_pT, gAsChic_pT, gAsJpsi_pT);



		// for now quick correction

		//TFile* fCor = new TFile(fileCorrection, "READ");

		//TH1D* hCor = (TH1D*)fCor->Get("h_chiEfficiency1D_Q_ratRel");
		//for (int i = 0; i < hCor->GetNbinsX(); i++)
		//{
		//	double corr = hCor->GetBinContent(i + 1); //TH1 numbering offset by 1
		//	cout << corr << endl; 
		//	gAsRatio_pT->GetY()[i] *= (1/corr);
		//	gAsRatio_pT->SetPointEYhigh(i, gAsRatio_pT->GetErrorYhigh(i)/corr);
		//	gAsRatio_pT->SetPointEYlow(i, gAsRatio_pT->GetErrorYlow(i)/corr);
		//}



		TCanvas* can_pT = new TCanvas("can_pT", "Ratio_pT", 800, 600);
		gAsRatio_pT->SetMarkerStyle(21);
		gAsRatio_pT->SetMarkerSize(2.0);
		gAsRatio_pT->SetMarkerColor(kBlue);
		gAsRatio_pT->SetLineWidth(3);
		gAsRatio_pT->Draw("AP");
		gAsRatio_pT->GetYaxis()->SetTitle("(Chic_1+Chic_2) / Jpsi");
		gAsRatio_pT->GetYaxis()->SetTitleOffset(1.3);
		gAsRatio_pT->GetXaxis()->SetTitle("pT");
		can_pT->SaveAs("FitterOutput/RatioChicJpsi_pT.png");

		////////  y

		GetRatio(gAsRatio_y, gAsChic_y, gAsJpsi_y);

		//TH1D* hCor_y = (TH1D*)fCor->Get("h_chiEfficiency1D_Q_ratRel");
		//for (int i = 0; i < hCor_y->GetNbinsX(); i++)
		//{
		//	double corr = hCor_y->GetBinContent(i + 1); //TH1 numbering offset by 1
		//	cout << corr << endl;
		//	gAsRatio_y->GetY()[i] *= (1 / corr);
		//	gAsRatio_y->SetPointEYhigh(i, gAsRatio_y->GetErrorYhigh(i) / corr);
		//	gAsRatio_y->SetPointEYlow(i, gAsRatio_y->GetErrorYlow(i) / corr);
		//}



		TCanvas* can_y = new TCanvas("can_y", "Ratio_y", 800, 600);
		gAsRatio_y->SetMarkerStyle(21);
		gAsRatio_y->SetMarkerSize(2.0);
		gAsRatio_y->SetMarkerColor(kRed);
		gAsRatio_y->SetLineWidth(3);
		gAsRatio_y->Draw("AP");
		gAsRatio_y->GetYaxis()->SetTitle("(Chic_1+Chic_2) / Jpsi");
		gAsRatio_y->GetYaxis()->SetTitleOffset(1.3);
		gAsRatio_y->GetXaxis()->SetTitle("rapidity");
		can_y->SaveAs("FitterOutput/RatioChicJpsi_y.png");



		GetRatio(gAsRatio_nTrk, gAsChic_nTrk, gAsJpsi_nTrk);

		TCanvas* can_nTrk = new TCanvas("can_nTrk", "Ratio_nTrk", 800, 600);
		gAsRatio_nTrk->SetMarkerStyle(21);
		gAsRatio_nTrk->SetMarkerSize(2.0);
		gAsRatio_nTrk->SetMarkerColor(kRed);
		gAsRatio_nTrk->SetLineWidth(3);
		gAsRatio_nTrk->Draw("AP");
		gAsRatio_nTrk->GetYaxis()->SetTitle("(Chic_1+Chic_2) / Jpsi");
		gAsRatio_nTrk->GetYaxis()->SetTitleOffset(1.3);
		gAsRatio_nTrk->GetXaxis()->SetTitle("nTrack");
		can_nTrk->SaveAs("FitterOutput/RatioChicJpsi_nTrk.png");



		/*
		TCanvas* can1 = new TCanvas("can1", "plot", 1200, 800);

		TH2F* h1 = hdxy_dz->Clone();
		//cout<<h1->Integral(101, 300, 1, 30)<<endl;

		h1->GetXaxis()->SetTitle("abseta");
		h1->GetYaxis()->SetTitle("pt");
		can1->SetLogz();
		h1->Draw("colz");
		*/

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
	hconv_deltaEta->Write();
	hconv_deltapT_Eta->Write();
	hconv_m->Write();
	hphoton_MVA_response->Write();

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

	//if (intConstrainedFit == 2) {
	//	TFile* fileConstr = new TFile(fileConstraints, "RECREATE");
	//	gAsChic_alpha_pT->Write();
	//	gAsChic_alpha_y->Write();
	//	gAsChic_alpha_nTrk->Write();
	//	gAsChic_n_pT->Write();
	//	gAsChic_n_y->Write();
	//	gAsChic_n_nTrk->Write();
	//	gAsChic_sigmaRat_pT->Write();
	//	gAsChic_sigmaRat_y->Write();
	//	gAsChic_sigmaRat_nTrk->Write();
	//	gAsChic_sigma1_pT->Write();
	//	gAsChic_sigma1_y->Write();
	//	gAsChic_sigma1_nTrk->Write();
	//	gAsChic_mean1_pT->Write();
	//	gAsChic_mean1_y->Write();
	//	gAsChic_mean1_nTrk->Write();
	//	fileConstr->Close();
	//}


	file_log.close();

}


