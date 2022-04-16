
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
#include "../ChiFitterInit.h"

using namespace std;
using namespace RooFit;

ofstream file_log;

/*const int PythCode_chic0 = 10441; //Pythia codes
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

*/



/*
double bins_pT[] = {6, 9, 12, 18, 30};
int  nbins_pT = sizeof(bins_pT) / sizeof(double) - 1;

double bins_y[] = {-2.4, -1.6, -1.0, 0, 1.0, 1.6, 2.4 };
int  nbins_y = sizeof(bins_y) / sizeof(double) - 1;

double bins_nTrk[] = { 0, 50, 100, 150, 200, 300, 400 };
int  nbins_nTrk = sizeof(bins_nTrk) / sizeof(double) - 1;

*/

int nParamOneFit = 5;//Since this is simultaneous fit, it is better to get ndf by hand


// output
const int nBinSets = 5;
const string nBinSetNames[nBinSets] = { "pt_all", "y", "nTrack", "pt_mid", "pt_fwd" };
const int nFitFunctionParams = 15; //Real number can be lower, is handled (but not higher)
TGraphAsymmErrors* gAsOutputArray [nBinSets][nFitFunctionParams] ; //stores the output of the fits

//below: old, replaced with the array above
TGraphAsymmErrors* gAsChic_alpha_pT, *gAsChic_alpha_y, *gAsChic_alpha_nTrk;
TGraphAsymmErrors* gAsChic_n_pT, *gAsChic_n_y, *gAsChic_n_nTrk;
TGraphAsymmErrors* gAsChic_sigmaRat_pT, *gAsChic_sigmaRat_y, *gAsChic_sigmaRat_nTrk;
TGraphAsymmErrors* gAsChic_sigma1_pT, *gAsChic_sigma1_y, *gAsChic_sigma1_nTrk;
TGraphAsymmErrors* gAsChic_mean1_pT, *gAsChic_mean1_y, *gAsChic_mean1_nTrk;


TLorentzVector* LVchic, *LVdimuon, *LVconv, *LVmuon1, *LVmuon2;
TLorentzVector* LVchic_rot, *LVchic_rotGamma, *LVdimuon_rot, *LVconv_rot, *LVmuon1_rot, *LVmuon2_rot;
TLorentzVector LVaux;



int GetRatio (TGraphAsymmErrors* gAsResult, TGraphAsymmErrors* gAsNum, TGraphAsymmErrors* gAsDen) // dependent on root version, this is prior 6.20 (should work anywhere)
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


bool CreateModelPdf(RooWorkspace& Ws, string pdfName)
{
	//if (pdfName.find("nominalPdf") != string::npos)
	if (pdfName.compare("nominalPdf") == 0)
	{

		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.05], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		Ws.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		Ws.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		Ws.factory("SUM::signal(c2toc1[0.4,0.1,1.]*chic2, chic1)");
		//Ws.factory("Exponential::expbkg(rvmass,e_1[-0.2,-0.5,0])");
		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("RooChebychev::background(rvmass, a1[0,-10,10])");


		Ws.factory("SUM::nominalPdf(nsig[50,0,1000000]*signal, nbkg[2000,0,100000]*background)");
	}
	else if (pdfName.compare("nominalPdfSim") == 0)
	{
		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.46, 3.53], sigma1[0.01, 0.003, 0.035], alpha[1.85, 0.1, 50], n[2.7, 1.1, 50])");
		Ws.factory("RooChebychev::background_chi1(rvmass, a1_chi1[0,-10,10])");
		Ws.factory("SUM::nominalPdf_chi1(nsig_chi1[50,0,100000]*chic1, nbkg_chi1[200,0,10000]*background_chi1)");
		
		Ws.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1.05, 1.00, 1.3])");
		Ws.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		Ws.factory("RooChebychev::background_chi2(rvmass, a1_chi2[0,-10,10])");
		Ws.factory("SUM::nominalPdf_chi2(nsig_chi2[50,0,100000]*chic2, nbkg_chi2[200,0,10000]*background_chi2)");
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

		Ws.factory("Exponential::backgroundJpsi(rvmassJpsi,e_1Jpsi[-0.2,-2.5,0])");

		Ws.factory("SUM::nominalPdfJpsi(nsigJpsi[5000,0,2000000]*signalJpsi, nbkgJpsi[2000,0,1000000]*backgroundJpsi)");
	}


	return true;
}



bool RefreshModel(RooWorkspace& Ws) // attempt to prevent the fits in the differential bins to get stuck in a weird state when they stop converge properly
{

	Ws.var("mean1")->removeError();
	Ws.var("sigma1")->removeError();
	Ws.var("alpha")->removeError();
	Ws.var("n")->removeError();
	Ws.var("sigmaRatio")->removeError();
	Ws.var("a1_chi1")->removeError();
	Ws.var("a1_chi2")->removeError();
	Ws.var("nsig_chi1")->removeError();
	Ws.var("nsig_chi2")->removeError();
	Ws.var("nbkg_chi1")->removeError();
	Ws.var("nbkg_chi2")->removeError();

	Ws.var("mean1")->setVal(3.5107);
	Ws.var("sigma1")->setVal(0.01);
	Ws.var("alpha")->setVal(1.85);
	Ws.var("n")->setVal(2.7);
	Ws.var("sigmaRatio")->setVal(1.05);
	Ws.var("a1_chi1")->setVal(0);
	Ws.var("a1_chi2")->setVal(0);
	Ws.var("nsig_chi1")->setVal(50);
	Ws.var("nsig_chi2")->setVal(50);
	Ws.var("nbkg_chi1")->setVal(200);
	Ws.var("nbkg_chi2")->setVal(200);

	return true;
}



int FitRooDataSetSim(TGraphAsymmErrors* gAsResult, double* bins, int nbins, RooRealVar* rvmass, RooWorkspace& Ws, RooCategory& rCatChi, string binVarName = "", RooSimultaneous* myPdf = NULL, string extraCut = "", string regionName = "") { //uses global variables
	//constrained fit: 0 no constraints, 1 yes, 2 no constraints, but write the fit results as constraints (for fitting MC)
	cout << endl << "IN FITTING " << binVarName << "   " << nbins << endl << endl;
	TCanvas *cFit = new TCanvas("cFit", "cFit", 1000, 600);
	//gPad->SetLeftMargin(0.15);
	cFit->cd();
	//cFit->Divide(1, 2);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.275, 0.98, 1.0);
	pad1->SetTicks(1, 1);
	pad1->Draw();

	//pull pad
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.006, 0.98, 0.360);
	pad2->SetTopMargin(0); // Upper and lower plot are joined
	pad2->SetBottomMargin(0.67);
	pad2->SetTicks(1, 1);


	for (int i = 0; i < nbins; i++) {
		cout << "Bins " << bins[i] << "  to  " << bins[i + 1] << endl;
		pad1->cd();
		RooPlot* massframeBin;
		massframeBin = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
		massframeBin->SetTitle("mass");
		massframeBin->GetYaxis()->SetTitleOffset(1.0);

		TString TstrCut = binVarName + TString::Format(" > %f", bins[i]) + " && " + binVarName + TString::Format(" < %f", bins[i + 1]) + extraCut;
		cout << TstrCut << endl;
		string strCut = TstrCut.Data();
		RooDataSet* rdsDataBin = (RooDataSet*)Ws.data("rdsNominal_chicBoth")->reduce(strCut.c_str()); 

		
		RefreshModel(Ws);
		cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
		RooFitResult* fitResultBin;
		if (i != 8) {
			fitResultBin = myPdf->fitTo(*rdsDataBin, Extended(true), SumW2Error(true), NumCPU(1), PrintLevel(-1), Save(true));// Range(mass_windowFit_l, mass_windowFit_h), Save(true));
		}
		//fitResultBin = myPdf->fitTo(*rdsDataBin, Extended(true), SumW2Error(true), NumCPU(1), PrintLevel(-1), Save(true));// Range(mass_windowFit_l, mass_windowFit_h), Save(true));
		cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;
		//cout << endl << "Importing fit result..." << endl;
		//myWs.import(*fitResSim);

		/////////////////////////
		//   R E A D   O U T  ///
		/////////////////////////

		int nBinSet; int nFitFunctionParam = 0;
		if (binVarName.compare("rvpt") == 0) {
			if (regionName.compare("all") == 0)
			{
				nBinSet = 0;
			}
			if (regionName.compare("midrap") == 0)
			{
				nBinSet = 3;
			}
			if (regionName.compare("fwdrap") == 0)
			{
				nBinSet = 4;
			}
			
		}
		else if (binVarName.compare("rvrap") == 0) {
			nBinSet = 1;
		}
		else if (binVarName.compare("rvntrack") == 0) {
			nBinSet = 2;
		}

		if (i != 2) {
			RooArgList paramList = fitResultBin->floatParsFinal();
			RooRealVar *par;
			for (int j = 0; j < paramList.getSize(); j++)
			{
				par = (RooRealVar*)paramList.at(j);
				//cout << par->getTitle() << endl;
				//cout << par->getVal() << endl;
				//cout << par->getError() << endl;

				// if the first bin, create the graph
				if (i == 0) {
					gAsOutputArray[nBinSet][j] = new TGraphAsymmErrors(nbins);
					cout << "gAsChic_Output_" + par->getTitle() + "_" + nBinSetNames[nBinSet] << endl;
					cout << "Chic " + nBinSetNames[nBinSet] + " " + par->getTitle() + " dependence" << endl;
					gAsOutputArray[nBinSet][j]->SetNameTitle("gAsChic_Output_" + par->getTitle() + "_" + nBinSetNames[nBinSet], "Chic " + nBinSetNames[nBinSet] + " " + par->getTitle() + " dependence");
				}

				// save as tgraphAsymmErrors
				SetPointFromFit(gAsOutputArray[nBinSet][j], (string)par->getTitle(), Ws, i, bins);
			}
		}
		//if (binVarName.compare("rvpt") == 0) {
		//	int nbin

		//gAsOutputArray []
		//		TGraphAsymmErrors* gAsOutputArray[nBinSets][nFitFunctionParams]; //stores the output of the fits

		//if (myPdfName.compare("nominalPdfSim") == 0) {
		//	if (binVarName.compare("rvpt") == 0) {
		//		SetPointFromFit(gAsChic_alpha_pT, "alpha", Ws, i, bins);
		//		SetPointFromFit(gAsChic_n_pT, "n", Ws, i, bins);
		//		SetPointFromFit(gAsChic_sigmaRat_pT, "sigmaRatio", Ws, i, bins);
		//		SetPointFromFit(gAsChic_sigma1_pT, "sigma1", Ws, i, bins);
		//		SetPointFromFit(gAsChic_mean1_pT, "mean1", Ws, i, bins);

		//	}
		//	else if (binVarName.compare("rvrap") == 0) {
		//		SetPointFromFit(gAsChic_alpha_y, "alpha", Ws, i, bins);
		//		SetPointFromFit(gAsChic_n_y, "n", Ws, i, bins);
		//		SetPointFromFit(gAsChic_sigmaRat_y, "sigmaRatio", Ws, i, bins);
		//		SetPointFromFit(gAsChic_sigma1_y, "sigma1", Ws, i, bins);
		//		SetPointFromFit(gAsChic_mean1_y, "mean1", Ws, i, bins);

		//	}
		//	else if (binVarName.compare("rvntrack") == 0) {
		//		SetPointFromFit(gAsChic_alpha_nTrk, "alpha", Ws, i, bins);
		//		SetPointFromFit(gAsChic_n_nTrk, "n", Ws, i, bins);
		//		SetPointFromFit(gAsChic_sigmaRat_nTrk, "sigmaRatio", Ws, i, bins);
		//		SetPointFromFit(gAsChic_sigma1_nTrk, "sigma1", Ws, i, bins);
		//		SetPointFromFit(gAsChic_mean1_nTrk, "mean1", Ws, i, bins);

		//	}
		//	else cout << "No match for the variable" << endl;

		//}
		//else cout << "The pdf not matching" << endl;
		//


		// PLOT

		rdsDataBin->plotOn(massframeBin, Cut("rCatChi==rCatChi::ChicOne"), Name("dataHistChicOne"), MarkerStyle(20), MarkerColor(kRed));
		rdsDataBin->plotOn(massframeBin, Cut("rCatChi==rCatChi::ChicTwo"), Name("dataHistChicTwo"), MarkerStyle(20), MarkerColor(kGreen+2));

		myPdf->plotOn(massframeBin, Slice(rCatChi, "ChicOne"), Name("pdfChicOneFull"), ProjWData(rCatChi, *rdsDataBin), LineColor(kRed));
		myPdf->plotOn(massframeBin, Slice(rCatChi, "ChicOne"), Components("chic1"), ProjWData(rCatChi, *rdsDataBin), LineStyle(kDashed), LineColor(kRed));

		myPdf->plotOn(massframeBin, Slice(rCatChi, "ChicTwo"), Name("pdfChicTwoFull"), ProjWData(rCatChi, *rdsDataBin), LineColor(kGreen+2));
		myPdf->plotOn(massframeBin, Slice(rCatChi, "ChicTwo"), Components("chic2"), ProjWData(rCatChi, *rdsDataBin), LineStyle(kDashed), LineColor(kGreen+2));

		myPdf->paramOn(massframeBin, Layout(0.75));


		//cout << "Print" << endl;
		//massframeBin->Print();
		double chi2_fit1 = massframeBin->chiSquare("pdfChicOneFull", "dataHistChicOne", nParamOneFit);
		double chi2_fit2 = massframeBin->chiSquare("pdfChicTwoFull", "dataHistChicTwo", nParamOneFit);
		//cout<< chi2_fit1 << endl;
		//double chi2_fit1Unreduced = massframeBin->chiSquare("pdfChicOneFull", "dataHistChicOne");
		//cout << "Unreduced: " << chi2_fit1Unreduced << endl;

		massframeBin->Draw();


		TPaveText *textchi1 = new TPaveText(.23, .75, .38, .83, "brNDC");
		textchi1->SetFillColor(kWhite);
		textchi1->SetTextColor(kRed);
		textchi1->SetTextSize(0.05);
		textchi1->SetBorderSize(0);
		string txtchi1 = "Chic 1 #chi^{2}/ndf: " + to_string(chi2_fit1);
		textchi1->AddText(txtchi1.c_str());
		textchi1->Draw();

		TPaveText *textchi2 = new TPaveText(.23, .65, .38, .73, "brNDC");
		textchi2->SetFillColor(kWhite);
		textchi2->SetTextColor(kGreen+2);
		textchi2->SetTextSize(0.05);
		textchi2->SetBorderSize(0);
		string txtchi2 = "Chic 2 #chi^{2}/ndf: " + to_string(chi2_fit2);
		textchi2->AddText(txtchi2.c_str());
		textchi2->Draw();

		//txt1->Draw();

		//PULLS
		pad2->cd();
		//RooHist* hpull = massframeBin->pullHist("rPullHist", myPdfName.c_str());
		RooHist* hpull = massframeBin->pullHist("dataHistChicOne","pdfChicOneFull");
		hpull->SetMarkerSize(0.8);
		hpull->SetMarkerColor(kRed);
		RooHist* hpull2 = massframeBin->pullHist("dataHistChicTwo", "pdfChicTwoFull");
		hpull2->SetMarkerSize(0.8);
		hpull2->SetMarkerColor(kGreen+2);
		RooPlot* pullFrame = rvmass->frame(Title("Pull Distribution"));
		pullFrame->addPlotable(hpull, "P");
		pullFrame->addPlotable(hpull2, "P");
		pullFrame->SetTitleSize(0);
		pullFrame->GetYaxis()->SetTitle("Pull");
		pullFrame->GetYaxis()->SetTitleSize(0.19);
		pullFrame->GetYaxis()->SetLabelSize(0.14);
		//pullFrame->GetYaxis()->SetLabelOffset(1.1);
		pullFrame->GetYaxis()->SetRangeUser(-5.0, 5.0);
		pullFrame->GetYaxis()->SetNdivisions(502, kTRUE);
		pullFrame->GetYaxis()->CenterTitle();
		pullFrame->GetXaxis()->SetLabelSize(0.20);
		pullFrame->GetXaxis()->SetTitle("#mu#mu#gamma - #mu#mu + 3.097 [GeV/c^{2}]");
		pullFrame->GetXaxis()->SetTitleOffset(1.05);
		pullFrame->GetXaxis()->SetTitleSize(0.20);

		//// Add text to frame
		//TText* txt = new TText(2, 100, "Signal");
		//txt->SetTextSize(0.04);
		//txt->SetTextColor(kRed);
		//txt->Draw("same");
		////massframeBin->addObject(txt);


		pullFrame->Draw();

		//Chic2
		//double chi2_fit1 = massframeBin->chiSquare("dataHistChicOne",0,20);
		cout << "    CHI2-1     " << massframeBin->chiSquare() << "   and calculated chi2-1:  " << chi2_fit1 << "   and calculated chi2-2:  " << chi2_fit2<< endl;
		cout << "    Unreduced chi2-1:  " << massframeBin->chiSquare("pdfChicOneFull", "dataHistChicOne", 0) << "   and unreduced chi2-2:  " << massframeBin->chiSquare("pdfChicTwoFull", "dataHistChicTwo", 0) << endl;
		//double chi2_fit2 = massframeBin->chiSquare("dataHistChicTwo", "pdfChicTwoFull", (nMassBins - nParamOneFit));








		TLine *l1 = new TLine(mass_windowFit_l, 0, mass_windowFit_h, 0);
		l1->SetLineStyle(9);
		l1->Draw("same");
		//pad1->Update();

		cFit->cd();
		pad1->Draw();
		pad2->Draw();

		pad1->Update();
		pad2->Update();


		cFit->SaveAs(((string)"FitterOutput/FitResultConstraint_Chic_" + regionName + "_" + binVarName + Form("_%ibin_", i) + Form("_%.1f_%.1f.png", bins[i], bins[i + 1])).c_str());

	}




	return 0;
}



///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////


//void ConstraintProducer(bool flagGenerateRds = true, bool flagRunFits = true,  const char* fileIn = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-MC8_BothDir.root", const char* fileOut = "Chi_c_output_MC8_test.root", const char* fileRds = "rds_MC8_Constraint.root", const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
void ConstraintProducer(bool flagGenerateRds = true, bool flagRunFits = true, const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC9-bothDir.root", const char* fileOut = "Chi_c_outConstr_MC9_WideRange.root", const char* fileRds = "rds_MC9_ConstraintWideRange.root", const char* fileConstraints = "Chi_c_constraints_MC9_WideRange.root", const char* fileCorrection = "Chi_c_WeightsMC9_bothDir.root")
{

	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
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
	//	if (nentries > 50000) { nentries = 5000; }
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
				if (ChiSelectionPassMC(iChi, 0) == false)continue; //check it is matched and passed
				// check Acceptance and Cuts
				int dimuonPos = chi_daughterJpsi_position->at(iChi);
				//if (DimuonSelectionPass(dimuonPos) == false) continue;
				passDimSel = DimuonSelectionPass(dimuonPos);
				passDimSelTight = DimuonSelectionPassTight(dimuonPos);
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



				////////////////////////
				////   refit     // comment out if not done
				///////////////////////


				double refit_vProb = 0, ctauPV = 0, ctauPVError = 0, ctauSig = 0, ctauPV3D = 0;

				//if (chi_kinematicRefitFlag->at(iChi) == 1 || chi_kinematicRefitFlag->at(iChi) == 3) { //good refits

				//	//cout << chi_refit_vprob->at(nRefitNumber) << endl;
				//	//cout << chi_refit_ctauPV->at(nRefitNumber) << endl;
				//	//cout << "RefitStored: " << chi_refitStored->at(nRefitNumber).mass() << endl;
				//	//m_chi = chi_refitStored->at(iChi).mass();//use refit mass
				//	ctauPV = chi_refit_ctauPV->at(iChi);
				//	ctauPVError = chi_refit_ctauErrPV->at(iChi);
				//	ctauSig = ctauPV / ctauPVError;
				//	ctauPV3D = chi_refit_ctauPV3D->at(iChi);
				//	refit_vProb = chi_refit_vprob->at(iChi);
				//}
				//else continue;//skip those that don't have it
				//if (refit_vProb < 0.01) continue;
				//Mdiff = chi_refitStored->at(iChi).mass();//use refit mass



				////////










				if (passDimSel == true && dimuonM > mass_cutoffJpsi_l && dimuonM < mass_cutoffJpsi_h && Mdiff > mass_window_l && Mdiff < mass_window_h) ++nchicCounterPassMass;



				// Obtain yields


				if (dimuonM<mass_cutoffJpsi_l || dimuonM>mass_cutoffJpsi_h) continue; //require narrow dimuon mass


				if (passDimSel == true) {


					//roofit:
					if (Mdiff > mass_window_l && Mdiff < mass_window_h) {
						rvmass->setVal(Mdiff);
						rvpt->setVal(pT_chi);
						rvrap->setVal(rap_chi);
						rvntrack->setVal(ntracks_inEvent); //for now using total
						//rvweight->setVal(1); // no weights for now
						double accEff_chi = hWeightChic->GetBinContent(hWeightChic->FindBin(abs(rap_chi), pT_chi));
 
						if (isMC) {
							rvweight->setVal(1); //no weights if MC
						}
						else if (accEff_chi > 0.00001) {// binning should ensure enough statistics, but if we have empty bin, just set weight to be 1
							rvweight->setVal(1 / accEff_chi);
						}
						else rvweight->setVal(1); 

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

		RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
		ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
		
		cout << "nEntries: " << rdsNominal->numEntries() << endl;
		rdsNominal = (RooDataSet*)rdsNominal->reduce(mass_windowFit.c_str());
		rdsNominal_chicOne = (RooDataSet*)rdsNominal_chicOne->reduce(mass_windowFit.c_str());
		rdsNominal_chicTwo = (RooDataSet*)rdsNominal_chicTwo->reduce(mass_windowFit.c_str());

		// construct the overall dataset for simultaneous fits

		RooCategory rCatChi("rCatChi", "rCatChi");
		rCatChi.defineType("ChicOne");
		rCatChi.defineType("ChicTwo");

		cout << "HERE" << endl;

		RooDataSet* rdsNominal_chicBoth = new RooDataSet("rdsNominal_chicBoth", "rdsNominal_chicBoth", RooArgSet(*rvmass, *rvpt, *rvrap, *rvntrack), Index(rCatChi), Import("ChicOne", *rdsNominal_chicOne), Import("ChicTwo", *rdsNominal_chicTwo));
		//RooDataSet* rdsNominal_chicBoth = new RooDataSet("rdsNominal_chicBoth", "rdsNominal_chicBoth", *rvmass, Index(rCatChi), Import("ChicOne", *rdsNominal_chicOne), Import("ChicTwo", *rdsNominal_chicTwo));

		//some issue with weight var, as it is not needed, dropping it
		cout << "******** New Combined Dataset ***********" << endl;
		rdsNominal_chicBoth->Print();
		myWs.import(*rdsNominal_chicBoth);

	//	RooPlot *massframeCheck = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
	//	massframeCheck->SetTitle("mass");
	//	rdsNominal->plotOn(massframeCheck, MarkerStyle(21), MarkerColor(kBlack));
	//	rdsNominal_chicBoth->plotOn(massframeCheck, MarkerStyle(20),	MarkerColor(kBlue));
	//	rdsNominal_chicOne->plotOn(massframeCheck, MarkerStyle(20), MarkerColor(kRed));
	//	rdsNominal_chicTwo->plotOn(massframeCheck, MarkerStyle(20), MarkerColor(kGreen+2));

	//	TCanvas *cCrosscheck = new TCanvas("cCrosscheck", "cCrosscheck", 1000, 600);
	//	gPad->SetLeftMargin(0.15);
	//	massframeCheck->GetYaxis()->SetTitleOffset(1.3);
	//	massframeCheck->Draw();
	//	
	//	cCrosscheck->SaveAs("FitterOutput/CanvasMassComparison.png");


	//	// C r e a t e   m o d e l   f o r   p h y s i c s   s a m p l e
	//// -------------------------------------------------------------

	//// Create observables
	//	//RooRealVar x("x", "x", -8, 8);

	//	// Construct signal pdf
	//	RooRealVar mean("mean", "mean", 3.5, 3, 4);
	//	RooRealVar sigma("sigma", "sigma", 0.03, 0.01, 0.1);
	//	RooGaussian gx("gx", "gx", *rvmass, mean, sigma);

	//	// Construct background pdf
	//	RooRealVar a0("a0", "a0", -0.1, -1, 1);
	//	RooRealVar a1("a1", "a1", 0.004, -1, 1);
	//	RooChebychev px("px", "px", *rvmass, RooArgSet(a0, a1));

	//	// Construct composite pdf
	//	RooRealVar f("f", "f", 0.2, 0., 1.);
	//	RooAddPdf model("model", "model", RooArgList(gx, px), f);



	//	// C r e a t e   m o d e l   f o r   c o n t r o l   s a m p l e
	//	// --------------------------------------------------------------

	//	// Construct signal pdf. 
	//	// NOTE that sigma is shared with the signal sample model
	//	RooRealVar mean_ctl("mean_ctl", "mean_ctl", 3.5, 3, 4);
	//	RooGaussian gx_ctl("gx_ctl", "gx_ctl", *rvmass, mean_ctl, sigma);

	//	// Construct the background pdf
	//	RooRealVar a0_ctl("a0_ctl", "a0_ctl", -0.1, -1, 1);
	//	RooRealVar a1_ctl("a1_ctl", "a1_ctl", 0.5, -0.1, 1);
	//	RooChebychev px_ctl("px_ctl", "px_ctl", *rvmass, RooArgSet(a0_ctl, a1_ctl));

	//	// Construct the composite model
	//	RooRealVar f_ctl("f_ctl", "f_ctl", 0.5, 0., 1.);
	//	RooAddPdf model_ctl("model_ctl", "model_ctl", RooArgList(gx_ctl, px_ctl), f_ctl);



	//	// C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s 
	//	// ---------------------------------------------------------------------------

	//	//// Define category to distinguish physics and control samples events
	//	//RooCategory sample("sample", "sample");
	//	//sample.defineType("physics");
	//	//sample.defineType("control");

	//	// Construct combined dataset in (*rvmass,sample)
	//	RooDataSet combData("combData", "combined data", *rvmass, Index(rCatChi), Import("ChicOne", *rdsNominal_chicOne), Import("ChicTwo", *rdsNominal_chicTwo));


	//	// C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
	//	// -----------------------------------------------------------------------------------

	//	// Construct a simultaneous pdf using category sample as index
	//	RooSimultaneous simPdf("simPdf", "simultaneous pdf", rCatChi);

	//	// Associate model with the physics state and model_ctl with the control state
	//	simPdf.addPdf(model, "ChicOne");
	//	simPdf.addPdf(model_ctl, "ChicTwo");



	//	// P e r f o r m   a   s i m u l t a n e o u s   f i t
	//	// ---------------------------------------------------

	//	// Perform simultaneous fit of model to data and model_ctl to data_ctl
	//	simPdf.fitTo(combData);



	//	// P l o t   m o d e l   s l i c e s   o n   d a t a    s l i c e s 
	//	// ----------------------------------------------------------------

	//	// Make a frame for the physics sample
	//	RooPlot* frame1 = rvmass->frame(Bins(30), Title("Physics sample"));

	//	// Plot all data tagged as physics sample
	//	combData.plotOn(frame1, Cut("rCatChi==rCatChi::ChicOne"));

	//	// Plot "physics" slice of simultaneous pdf. 
	//	// NBL You _must_ project the sample index category with data using ProjWData 
	//	// as a RooSimultaneous makes no prediction on the shape in the index category 
	//	// and can thus not be integrated
	//	simPdf.plotOn(frame1, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, combData), LineColor(kRed));
	//	simPdf.plotOn(frame1, Slice(rCatChi, "ChicOne"), Components("px"), ProjWData(rCatChi, combData), LineStyle(kDashed), LineColor(kRed));

	//	// The same plot for the control sample slice
	//	RooPlot* frame2 = rvmass->frame(Bins(30), Title("Control sample"));
	//	combData.plotOn(frame2, Cut("rCatChi==rCatChi::ChicTwo"));
	//	simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, combData));
	//	simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), Components("px_ctl"), ProjWData(rCatChi, combData), LineStyle(kDashed));
	//	simPdf.paramOn(frame2, Layout(0.55));


	//	TCanvas* c = new TCanvas("rf501_simultaneouspdf", "rf501_simultaneouspdf", 800, 400);
	//	c->Divide(2);
	//	c->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
	//	c->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw();

	//	c->SaveAs("Test.png");




		//////////////////////////////
		//Construct simultaneous PDF
		////////////////////////////


		string myPdfName = "nominalPdfSim";
		CreateModelPdf(myWs, myPdfName);


		RooSimultaneous* simPdfChi = new RooSimultaneous("simPdfChi", "simPdfChi", rCatChi);
		simPdfChi->addPdf(*myWs.pdf("nominalPdf_chi1"), "ChicOne");
		simPdfChi->addPdf(*myWs.pdf("nominalPdf_chi2"), "ChicTwo");
		myWs.import(*simPdfChi); //import doesn't work well for some reason, issue with plotting the individual components
















		int i = 2;
		TString TstrCut = TString::Format("rvrap > %f", bins_y[i]) + " && " + TString::Format("rvrap < %f", bins_y[i + 1]) + " && rvpt>6.5 && rvpt <30";
		cout << TstrCut << endl;
		string strCut = TstrCut.Data();
		RooDataSet* rdsDataBin = (RooDataSet*)myWs.data("rdsNominal_chicBoth")->reduce(strCut.c_str());


		cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
		RooFitResult* fitResSim = myWs.pdf("simPdfChi")->fitTo(*rdsDataBin, Extended(true), SumW2Error(true), NumCPU(1), Save(true));
		//myWs.pdf("simPdfChi")->fitTo(*rdsDataBin, Extended(true), SumW2Error(true), NumCPU(1), Save(true));
		cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;
		//cout << endl << "Importing fit result..." << endl;
		//myWs.import(*fitResSim);


		RooPlot *massframe = rvmass->frame(mass_windowFit_l, mass_windowFit_h, nMassBins);
		massframe->SetTitle("mass");


		//myWs.pdf("simPdf")->plotOn(massframe);
		//rdsNominal_chicOne->plotOn(massframe, MarkerStyle(20), MarkerColor(kRed));
		//rdsNominal_chicTwo->plotOn(massframe, MarkerStyle(20), MarkerColor(kGreen+2));
		//rdsNominal_chicBoth->plotOn(massframe, MarkerStyle(20), MarkerColor(kBlue));
		rdsDataBin->plotOn(massframe, Cut("rCatChi==rCatChi::ChicOne"), MarkerStyle(20), MarkerColor(kRed));
		rdsDataBin->plotOn(massframe, Cut("rCatChi==rCatChi::ChicTwo"), MarkerStyle(20), MarkerColor(kGreen+2));
		//myWs.pdf("simPdf")->paramOn(massframe, Layout(0.55));
		//myWs.pdf("simPdf")->plotOn(massframe);//, Components("background"), LineStyle(kDashed));
	//	myWs.pdf("nominalPdf_chi1")->plotOn(massframe, LineColor(kRed));
		//myWs.pdf("nominalPdf_chi2")->plotOn(massframe, LineColor(kGreen+2));
		//myWs.pdf("nominalPdf_chi1")->plotOn(massframe, Components("chic1"), Range(mass_windowFit_l, mass_windowFit_h), LineStyle(kDashed), LineColor(kRed));
		//myWs.pdf("nominalPdf_chi2")->plotOn(massframe, Components("chic2"), Range(mass_windowFit_l, mass_windowFit_h), LineStyle(kDashed), LineColor(kGreen+2));


		simPdfChi->plotOn(massframe, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, *rdsDataBin), LineColor(kRed));
		simPdfChi->plotOn(massframe, Slice(rCatChi, "ChicOne"), Components("chic1"), ProjWData(rCatChi, *rdsDataBin), LineStyle(kDashed), LineColor(kRed));


		simPdfChi->plotOn(massframe, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, *rdsDataBin), LineColor(kGreen+2));
		simPdfChi->plotOn(massframe, Slice(rCatChi, "ChicTwo"), Components("chic2"), ProjWData(rCatChi, *rdsDataBin), LineStyle(kDashed), LineColor(kGreen+2));

		simPdfChi->paramOn(massframe, Layout(0.55));


		//myWs.pdf("simPdfChi")->plotOn(massframe, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, *rdsDataBin), LineColor(kRed));
		//myWs.pdf("simPdfChi")->plotOn(massframe, Slice(rCatChi, "ChicOne"), Components("chic1"), ProjWData(rCatChi, *rdsDataBin), LineStyle(kDashed), LineColor(kRed));


		//myWs.pdf("simPdfChi")->plotOn(massframe, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, *rdsDataBin), LineColor(kGreen+2));
		//myWs.pdf("simPdfChi")->plotOn(massframe, Slice(rCatChi, "ChicTwo"), Components("chic2"), ProjWData(rCatChi, *rdsDataBin), LineStyle(kDashed), LineColor(kGreen+2));

		//myWs.pdf("simPdfChi")->paramOn(massframe, Layout(0.55));

		//simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, combData));
		//simPdf.plotOn(frame2, Slice(rCatChi, "ChicTwo"), Components("px_ctl"), ProjWData(rCatChi, combData), LineStyle(kDashed));
		//simPdf.paramOn(frame2, Layout(0.55));

		//myWs.pdf("simPdfChi")->plotOn(massframe, Slice(rCatChi, "ChicOne"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kRed));
		//myWs.pdf("simPdfChi")->plotOn(massframe, Slice(rCatChi, "ChicOne"), Components("chic2"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kRed));
		//myWs.pdf("simPdf")->plotOn(massframe, Slice(rCatChi, "ChicTwo"), Components("chicdgsdgs"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineStyle(kDashed), LineColor(kRed));

		//myWs.pdf("simPdf")->plotOn(massframe, Slice(rCatChi, "ChicTwo"), ProjWData(rCatChi, *rdsNominal_chicBoth), LineColor(kGreen+2));
		//myWs.pdf("simPdfO")->paramOn(massframe, Layout(0.55));

		TCanvas *cTest = new TCanvas("cTest", "cTest", 1000, 600);
		gPad->SetLeftMargin(0.15);
		massframe->GetYaxis()->SetTitleOffset(1.3);
		massframe->Draw();


		cTest->SaveAs("CanvasCTest_RW3.png");
		//*/
		



		FitRooDataSetSim(gAsChic_pT, bins_pT, nbins_pT, rvmass, myWs, rCatChi, "rvpt", simPdfChi, " && rvrap>-2.4 && rvrap <2.4", "all");
		//FitRooDataSetSim(gAsChic_y, bins_y, nbins_y, rvmass, myWs, rCatChi, "rvrap", simPdfChi, " && rvpt>6.5 && rvpt <30", "all");
		
		//FitRooDataSetSim(gAsChic_pT, bins_pT, nbins_pT, rvmass, myWs, rCatChi, "rvpt", simPdfChi, " && rvrap>-1 && rvrap <1", "midrap");
		//FitRooDataSetSim(gAsChic_pT, bins_pT, nbins_pT, rvmass, myWs, rCatChi, "rvpt", simPdfChi, " && (rvrap<-1 || rvrap >1)", "fwdrap");
		//FitRooDataSetSim(gAsChic_nTrk, bins_nTrk, nbins_nTrk, rvmass, myWs, rCatChi, "rvntrack", simPdfChi, " && rvrap>-1 && rvrap <1", "midrap");

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
			myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic2"), LineStyle(kDashed), LineColor(kGreen+2));


			massframeBin->Draw();
			cTest->SaveAs(Form("CanvasCTest_RW3_pT_%.1f_%.1f.png", bins_pT[i], bins_pT[i + 1]));

		}*/

		//FitRooDataSet(gAsChic_pT, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "midrap", intConstrainedFit);
		//FitRooDataSet(gAsChic_y, bins_y, nbins_y, rvmass, myWs, false, "rvrap", myPdfName.c_str(), " && rvpt>6 && rvpt <30", "all", intConstrainedFit);
		//FitRooDataSet(gAsChic_nTrk, bins_nTrk, nbins_nTrk, rvmass, myWs, false, "rvntrack", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "midrap", intConstrainedFit);


	}



	TFile* fout = new TFile(fileOut, "RECREATE");


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

	
	TFile* fileConstr = new TFile(fileConstraints, "RECREATE");

	// save the values of the params
	for (int i = 0; i < nBinSets; i++)
	{
		for (int j = 0; j < nFitFunctionParams; j++)
		{
			if (gAsOutputArray[i][j] != NULL)
			{
				gAsOutputArray[i][j]->GetXaxis()->SetTitle(nBinSetNames[i].c_str());
				gAsOutputArray[i][j]->GetXaxis()->SetLabelSize(0.05);
				gAsOutputArray[i][j]->GetXaxis()->SetLabelSize(0.05);
				gAsOutputArray[i][j]->GetXaxis()->SetTitleSize(0.05);
				gAsOutputArray[i][j]->GetXaxis()->SetTitleOffset(1.05);
				gAsOutputArray[i][j]->SetMarkerStyle(21);
				gAsOutputArray[i][j]->SetMarkerSize(1.1);
				gAsOutputArray[i][j]->SetMarkerColor(kBlue);

				gAsOutputArray[i][j]->Write();
			}
		}
	}

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





	file_log.close();

}


