
// Macro to analyze the chi_c. Starting with the event tree


#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>

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
#include "TPaveText.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TClonesArray.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"
#include "TMath.h"
#include "Math/SpecFunc.h" //for Bechsel
#include "Math/Minimizer.h" // minimizer settings
#include "RooNumIntConfig.h"

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
#include "RooHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooMsgService.h"



#include "../../tdrstyle.C"
#include "../../CMS_lumi.C"
#include "../../ChiTreeInit.C"
#include "../../ChiFitterInit.h"




using namespace std;
using namespace RooFit;

ofstream file_log;

string myPdfName = "JpsiTwoCB_Exp";
string myPdfNameSig = "JpsiHypatiaGaus_Exp";
string myPdfNameBkg = "JpsiTwoCB_Pol";

///*
TGraphAsymmErrors* gAsChic_alpha_pT, *gAsChic_alpha_y, *gAsChic_alpha_nTrk;
TGraphAsymmErrors* gAsChic_n_pT, *gAsChic_n_y, *gAsChic_n_nTrk;
TGraphAsymmErrors* gAsChic_sigmaRat_pT, *gAsChic_sigmaRat_y, *gAsChic_sigmaRat_nTrk;
TGraphAsymmErrors* gAsChic_sigma1_pT, *gAsChic_sigma1_y, *gAsChic_sigma1_nTrk;
TGraphAsymmErrors* gAsChic_mean1_pT, *gAsChic_mean1_y, *gAsChic_mean1_nTrk;
//*/


TGraph* gSystOutputArraySig[nFittingSets]; // stores the systematic uncertainty - total
TGraph* gSystOutputArraySig_mean[nFittingSets]; // stores the systematic uncertainty - mean only
TGraph* gSystOutputArraySig_rms[nFittingSets]; // stores the systematic uncertainty - RMS only
TGraph* gSystOutputArraySig_failed[nFittingSets]; // stores the number of failed toys

TGraph* gSystOutputArrayBkg[nFittingSets]; // stores the systematic uncertainty - total
TGraph* gSystOutputArrayBkg_mean[nFittingSets]; // stores the systematic uncertainty - mean only
TGraph* gSystOutputArrayBkg_rms[nFittingSets]; // stores the systematic uncertainty - RMS only
TGraph* gSystOutputArrayBkg_failed[nFittingSets]; // stores the number of failed toys



class RooDoubleCB : public RooAbsPdf {  //double sided cb, taken from Alberto's code
public:
	RooDoubleCB() {};
	RooDoubleCB(const char *name, const char *title,
		RooAbsReal& _x,
		RooAbsReal& _mu,
		RooAbsReal& _sig,
		RooAbsReal& _a1,
		RooAbsReal& _n1,
		RooAbsReal& _a2,
		RooAbsReal& _n2);
	RooDoubleCB(const RooDoubleCB& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const { return new RooDoubleCB(*this, newname); }
	inline virtual ~RooDoubleCB() { }

protected:

	RooRealProxy x;
	RooRealProxy mu;
	RooRealProxy sig;
	RooRealProxy a1;
	RooRealProxy n1;
	RooRealProxy a2;
	RooRealProxy n2;

	Double_t evaluate() const;

private:

	ClassDef(RooDoubleCB, 1)
};

///// Hypatia  - from roofit team - copied here since it requires root 6.2x

class RooHypatia2 : public RooAbsPdf {
public:
	RooHypatia2() {};
	RooHypatia2(const char *name, const char *title,
		RooAbsReal& x, RooAbsReal& lambda, RooAbsReal& zeta, RooAbsReal& beta,
		RooAbsReal& sigma, RooAbsReal& mu, RooAbsReal& a, RooAbsReal& n, RooAbsReal& a2, RooAbsReal& n2);
	RooHypatia2(const RooHypatia2& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const override { return new RooHypatia2(*this, newname); }
	inline virtual ~RooHypatia2() { }


private:
	RooRealProxy _x;
	RooRealProxy _lambda;
	RooRealProxy _zeta;
	RooRealProxy _beta;
	RooRealProxy _sigma;
	RooRealProxy _mu;
	RooRealProxy _a;
	RooRealProxy _n;
	RooRealProxy _a2;
	RooRealProxy _n2;

	Double_t evaluate() const override;

	/// \cond CLASS_DEF_DOXY
	ClassDefOverride(RooHypatia2, 1);
	/// \endcond
};

int GetRatio(TGraphAsymmErrors* gAsResult, TGraphAsymmErrors* gAsNum, TGraphAsymmErrors* gAsDen) // dependent on root version, this is prior 6.20
{
	if (gAsResult->GetN() != gAsNum->GetN() || gAsResult->GetN() != gAsDen->GetN()) {
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
			gAsResult->SetPointError(i, 0, 0, 0, 0);
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

	if (pdfName.compare("DoubleSidedCB_BPHTreshold") == 0)
	{
		Ws.factory("RooDoubleCB::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50], alphaH[1.85, 0.1, 10], nH[2.7, 1.1, 50])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1.11])");
		Ws.factory("RooDoubleCB::chic2(rvmass, mean2[3.556, 3.53, 3.58], sigma2, alpha, n, alphaH, nH)");
		Ws.factory("SUM::signal(chic1, c2ratio[0.3,0.05,1.0]*chic2)");

		Ws.factory("expr::delta('rvmass-qzero', rvmass, qzero[3.2])");
		Ws.factory("expr::b1('beta_bkg*(rvmass-qzero)', rvmass, qzero, beta_bkg[-2, -3.3, -0.5])");
		Ws.factory("EXPR::background('pow(delta,alpha_bkg)*exp(b1)',delta, b1, alpha_bkg[0.2, 0.1, 0.8])");

		if (bConstrainedFit == true) {

			Ws.var("alpha")->setConstant(true);
			Ws.var("n")->setConstant(true);
			Ws.var("alphaH")->setConstant(true);
			Ws.var("nH")->setConstant(true);
			Ws.var("sigma1")->setConstant(true);
			Ws.var("mean1")->setConstant(true);
			Ws.var("mean2")->setConstant(true);
		}

		Ws.factory("SUM::DoubleSidedCB_BPHTreshold(nsig[5000,2,20000]*signal, nbkg[2000,0,50000]*background)");
		return true;
	}
	else if (pdfName.compare("DoubleSidedCB_D0BG") == 0)
	{
		Ws.factory("RooDoubleCB::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50], alphaH[1.85, 0.1, 10], nH[2.7, 1.1, 50])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1.11])");
		Ws.factory("RooDoubleCB::chic2(rvmass, mean2[3.556, 3.53, 3.58], sigma2, alpha, n, alphaH, nH)");
		Ws.factory("SUM::signal(chic1, c2ratio[0.3,0.05,1.0]*chic2)");
		Ws.factory("RooDstD0BG::background(rvmass, massOffset[3.15,3.1,3.2], A_bkg[0.04,0.015,0.4], B_bkg[0], C_bkg[2,-4,4])");
		Ws.factory("SUM::DoubleSidedCB_D0BG(nsig[5000,2,20000]*signal, nbkg[2000,0,50000]*background)");
		return true;
	}
	else if (pdfName.compare("Hypatia_BPHTreshold") == 0)
	{
		Ws.factory("RooHypatia2::chic1(rvmass, lambda[-1, -25,0], zeta[0], beta[-0.001, -15, -0.0001], sigma1[0.02, 0.007, 0.034], mean1[3.5107, 3.46, 3.53], aSig1[1.85, 0.1, 10], nSig1[2.7, 1.1, 50], aSig2[1.85, 0.1, 10], nSig2[2.7, 1.1, 50])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1.11])");
		Ws.factory("RooHypatia2::chic2(rvmass, lambda, zeta, beta, sigma2, mean2[3.556, 3.53, 3.58], aSig1,nSig1,aSig2,nSig2)");
		Ws.factory("SUM::signal(chic1, c2ratio[0.3,0.05,1.0]*chic2)");

		Ws.factory("expr::delta('rvmass-qzero', rvmass, qzero[3.2])");
		Ws.factory("expr::b1('beta_bkg*(rvmass-qzero)', rvmass, qzero, beta_bkg[-2, -3.3, -0.5])");
		Ws.factory("EXPR::background('pow(delta,alpha_bkg)*exp(b1)',delta, b1, alpha_bkg[0.2, 0.1, 0.8])");

		if (bConstrainedFit == true) {

			Ws.var("lambda")->setConstant(true);
			Ws.var("zeta")->setConstant(true);
			Ws.var("beta")->setConstant(true);
			Ws.var("sigma1")->setConstant(true);
			Ws.var("mean1")->setConstant(true);
			Ws.var("aSig1")->setConstant(true);
			Ws.var("nSig1")->setConstant(true);
			Ws.var("aSig2")->setConstant(true);
			Ws.var("nSig2")->setConstant(true);
			Ws.var("mean2")->setConstant(true);
		}

		Ws.factory("SUM::Hypatia_BPHTreshold(nsig[50,0,20000]*signal, nbkg[2000,0,10000]*background)");
		return true;
	}
	else if (pdfName.compare("SingleCB_BPHTreshold") == 0)
	{
		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		Ws.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		Ws.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		Ws.factory("SUM::signal(chic1, c2ratio[0.3,0.05,1.0]*chic2)");

		Ws.factory("expr::delta('rvmass-qzero', rvmass, qzero[3.2])");
		Ws.factory("expr::b1('beta_bkg*(rvmass-qzero)', rvmass, qzero, beta_bkg[-2, -3.3, -0.5])");
		Ws.factory("EXPR::background('pow(delta,alpha_bkg)*exp(b1)',delta, b1, alpha_bkg[0.2, 0.1, 0.8])");

		if (bConstrainedFit == true) {

			Ws.var("alpha")->setConstant(true);
			Ws.var("n")->setConstant(true);
			Ws.var("sigmaRatio")->setConstant(true);
			Ws.var("sigma1")->setConstant(true);
			Ws.var("mean1")->setConstant(true);
		}

		Ws.factory("SUM::SingleCB_BPHTreshold(nsig[5000,2,20000]*signal, nbkg[2000,0,50000]*background)");
		return true;
	}
	else if (pdfName.compare("SingleCB_D0BG") == 0)
	{
		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		Ws.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		Ws.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		Ws.factory("SUM::signal(chic1, c2ratio[0.3,0.05,1.0]*chic2)");
		Ws.factory("RooDstD0BG::background(rvmass, massOffset[3.15,3.1,3.2], A_bkg[0.04,0.015,0.4], B_bkg[0], C_bkg[2,-4,4])");

		if (bConstrainedFit == true) {

			Ws.var("alpha")->setConstant(true);
			Ws.var("n")->setConstant(true);
			Ws.var("sigmaRatio")->setConstant(true);
			Ws.var("sigma1")->setConstant(true);
			Ws.var("mean1")->setConstant(true);
		}

		Ws.factory("SUM::SingleCB_D0BG(nsig[5000,20,10000]*signal, nbkg[2000,0,10000]*background)");
		return true;
	}
	else if (pdfName.compare("nominalFlatBkg") == 0)
	{
		Ws.factory("CBShape::chic1(rvmass, mean1[3.5107, 3.48, 3.52], sigma1[0.005, 0.003, 0.08], alpha[1.85, 0.1, 50], n[1.7, 0.2, 50])");
		Ws.factory("prod::mean2(mean1, massRatioPDG[1.01296])");
		Ws.factory("prod::sigma2(sigma1, sigmaRatio[1, 0.95,1.2])");
		Ws.factory("CBShape::chic2(rvmass, mean2, sigma2, alpha, n)");
		Ws.factory("SUM::signal(chic1, c2ratio[0.3,0.05,1.0]*chic2)");
		Ws.factory("RooChebychev::background(rvmass, a1[0,-10,10])");

		if (bConstrainedFit == true) {

			Ws.var("alpha")->setConstant(true);
			Ws.var("n")->setConstant(true);
			Ws.var("sigmaRatio")->setConstant(true);
			Ws.var("sigma1")->setConstant(true);
			Ws.var("mean1")->setConstant(true);
		}

		Ws.factory("SUM::nominalFlatBkg(nsig[50,0,10000000]*signal, nbkg[2000,0,10000000]*background)");
		return true;
	}
	else if (pdfName.compare("JpsiTwoCB_Exp") == 0)
		//else if (pdfName.find("nominalPdfJpsi") != string::npos)
	{
		//Ws.factory("CBShape::signalJpsi(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08], alphaJpsi[1.85, 0.1, 50], nJpsi[1.7, 0.2, 50])");
		Ws.factory("CBShape::Jpsi1(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.01, 0.08], alphaJpsi[1.85, 0.1, 5], nJpsi[1.7, 0.5, 5])");
		Ws.factory("prod::sigma2Jpsi(sigma1Jpsi, sigmaRatioJpsi[2.0, 1.2, 3.0])");// allowing sigmas to merge (i.e. 1.0) leads to fit instabilities, and is useless (same as one peak - so the shape can be fitted here too). In any case, is 1.5-1.6 in the fits
		Ws.factory("CBShape::Jpsi2(rvmassJpsi, mean1Jpsi, sigma2Jpsi, alphaJpsi, nJpsi)");
		//Ws.factory("CBShape::psi2(rvmassJpsi, mean2Jpsi[3.686, 3.50, 3.75], sigma3Jpsi[0.03, 0.003, 0.08], alphaJpsi, nJpsi)");
		Ws.factory("SUM::signalJpsi(ratioJpsi[0.5,0.00,1.0]*Jpsi1, Jpsi2)");
		//Ws.factory("SUM::signalJpsi(ratioPsi[0.1,0.00,2.0]*psi2, Jpsi)");

		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("Exponential::backgroundJpsi(rvmassJpsi,e_1Jpsi[-0.2,-0.7,0])");

		Ws.factory("SUM::JpsiTwoCB_Exp(nsigJpsi[5000,0,400000]*signalJpsi, nbkgJpsi[2000,0,200000]*backgroundJpsi)");
		return true;
	}
	else if (pdfName.compare("JpsiHypatiaGaus_Exp") == 0)
	{
		//Ws.factory("RooGaussian::Jpsi1(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08])");
		//Ws.factory("RooHypatia2::signalJpsi(rvmassJpsi, lambdaJpsi[-1, -25,0], zetaJpsi[0], betaJpsi[-0.001, -15, -0.0001], sigma1Jpsi[0.08, 0.12, 0.18], mean1Jpsi[3.097, 3.05, 3.15], aSig1Jpsi[1.85, 0.1, 10], nSig1Jpsi[2.7, 1.1, 50], aSig2Jpsi[1.85, 0.1, 10], nSig2Jpsi[2.7, 1.1, 50])");
		//Ws.factory("SUM::signalJpsi(ratioJpsi[0.5,0.00,1.0]*Jpsi1, Jpsi2)");
		Ws.factory("RooGaussian::Jpsi1(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.04, 0.10])");
		//Ws.factory("RooHypatia2::Jpsi2(rvmassJpsi, lambdaJpsi[-1, -30, 0], zetaJpsi[0], betaJpsi[-0.001, -30, -0.0001], sigma2Jpsi[0.08, 0.12, 0.16], mean1Jpsi, aSig1Jpsi[0.5, 0.1, 4], nSig1Jpsi[1.7, 1.1, 5], aSig2Jpsi[1.85, 0.1, 10], nSig1Jpsi)");
		Ws.factory("RooHypatia2::Jpsi2(rvmassJpsi, lambdaJpsi[-1, -30, 0], zetaJpsi[0], betaJpsi[-2, -20, -0.0001], sigma2Jpsi[0.12, 0.08, 0.16], mean1Jpsi, aSig1Jpsi[0.5, 0.1, 4], nSig1Jpsi[2.0], aSig2Jpsi[1.85, 0.1, 10], nSig1Jpsi)");

		Ws.factory("SUM::signalJpsi(ratioJpsi[0.5,0.00,1.0]*Jpsi1, Jpsi2)");

		Ws.factory("Exponential::backgroundJpsi(rvmassJpsi, e_1Jpsi[-0.2,-0.9,0])");

		//Ws.factory("SUM::JpsiHypatiaGaus_Exp(nsigJpsi[5000,0,400000]*signalJpsi, nbkgJpsi[2000,0,200000]*backgroundJpsi)"); //worked for bins 3,4 (fwdCMS, bkwCMS) 
		Ws.factory("SUM::JpsiHypatiaGaus_Exp(nsigJpsi[100000,0,500000]*signalJpsi, nbkgJpsi[6000,0,200000]*backgroundJpsi)"); //worked for bins 3,4 (fwdCMS, bkwCMS)
		return true;
	}
	else if (pdfName.compare("JpsiTwoCB_Pol") == 0)
	{
		Ws.factory("CBShape::Jpsi1(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08], alphaJpsi[1.85, 0.1, 5], nJpsi[1.7, 0.5, 5])");
		Ws.factory("prod::sigma2Jpsi(sigma1Jpsi, sigmaRatioJpsi[2.0, 1.2, 3.0])");
		Ws.factory("CBShape::Jpsi2(rvmassJpsi, mean1Jpsi, sigma2Jpsi, alphaJpsi, nJpsi)");
		Ws.factory("SUM::signalJpsi(ratioJpsi[0.5,0.00,1.0]*Jpsi1, Jpsi2)");

		Ws.factory("RooChebychev::backgroundJpsi(rvmassJpsi, a1Jpsi[-0.1,-10,10])");
		//Ws.factory("RooChebychev::backgroundJpsi(rvmassJpsi, {a1Jpsi[0,-10,10], a2Jpsi[0,-10,10]})");

		Ws.factory("SUM::JpsiTwoCB_Pol(nsigJpsi[5000,0,400000]*signalJpsi, nbkgJpsi[2000,0,200000]*backgroundJpsi)");
		return true;
	}

	cout << " ERROR: MODEL CREATION FAILED, POSSIBLY UNDEFINED MODEL WAS CHOSEN" << endl;
	return false;
}



bool RefreshModel(RooWorkspace& Ws, string pdfName, bool isJpsi) // attempt to prevent the fits in the differential bins to get stuck in a weird state when they stop converging properly (empty bins making errors too small, and those errors used to determine range in the next bin)
{
	cout << "Refreshing model...";
	if (isJpsi == true) {

		if (pdfName.compare("JpsiTwoCB_Exp") == 0)
		{
			Ws.var("mean1Jpsi")->removeError();
			Ws.var("sigma1Jpsi")->removeError();
			Ws.var("sigmaRatioJpsi")->removeError();
			Ws.var("alphaJpsi")->removeError();
			Ws.var("nJpsi")->removeError();
			Ws.var("e_1Jpsi")->removeError();
			Ws.var("ratioJpsi")->removeError();
			Ws.var("nsigJpsi")->removeError();
			Ws.var("nbkgJpsi")->removeError();

			Ws.var("mean1Jpsi")->setVal(3.097);
			Ws.var("sigma1Jpsi")->setVal(0.03);
			Ws.var("sigmaRatioJpsi")->setVal(2.0);
			Ws.var("alphaJpsi")->setVal(1.85);
			Ws.var("nJpsi")->setVal(1.7);
			Ws.var("e_1Jpsi")->setVal(-0.4);
			Ws.var("ratioJpsi")->setVal(0.5);
			Ws.var("nsigJpsi")->setVal(50000);
			Ws.var("nbkgJpsi")->setVal(2000);
		}
		else if (pdfName.compare("JpsiHypatiaGaus_Exp") == 0)
		{

			Ws.var("mean1Jpsi")->removeError();
			Ws.var("sigma1Jpsi")->removeError();
			Ws.var("sigma2Jpsi")->removeError();
			Ws.var("lambdaJpsi")->removeError();
			Ws.var("betaJpsi")->removeError();
			Ws.var("aSig1Jpsi")->removeError();
			//Ws.var("nSig1Jpsi")->removeError();
			Ws.var("aSig2Jpsi")->removeError();
			//Ws.var("nSig2Jpsi")->removeError();
			Ws.var("e_1Jpsi")->removeError();
			Ws.var("ratioJpsi")->removeError();
			Ws.var("nsigJpsi")->removeError();
			Ws.var("nbkgJpsi")->removeError();

			Ws.var("mean1Jpsi")->setVal(3.097);
			Ws.var("sigma1Jpsi")->setVal(0.04);
			Ws.var("sigma2Jpsi")->setVal(0.12);
			Ws.var("lambdaJpsi")->setVal(-5);
			Ws.var("betaJpsi")->setVal(-2);
			Ws.var("aSig1Jpsi")->setVal(0.5);
			//Ws.var("nSig1Jpsi")->setVal(1.7);
			Ws.var("aSig2Jpsi")->setVal(2);
			//Ws.var("nSig2Jpsi")->setVal(2);
			Ws.var("e_1Jpsi")->setVal(-0.4);
			Ws.var("ratioJpsi")->setVal(0.5);
			Ws.var("nsigJpsi")->setVal(100000); //50000 worked for fwd and bkg (bins 3,4)
			Ws.var("nbkgJpsi")->setVal(6000);//2000 worked for fwd and bkg (bins 3,4)
		}
		else if (pdfName.compare("JpsiTwoCB_Pol") == 0)
		{
			Ws.var("mean1Jpsi")->removeError();
			Ws.var("sigma1Jpsi")->removeError();
			Ws.var("sigmaRatioJpsi")->removeError();
			Ws.var("alphaJpsi")->removeError();
			Ws.var("nJpsi")->removeError();
			Ws.var("a1Jpsi")->removeError();
			//Ws.var("a2Jpsi")->removeError();
			Ws.var("ratioJpsi")->removeError();
			Ws.var("nsigJpsi")->removeError();
			Ws.var("nbkgJpsi")->removeError();

			Ws.var("mean1Jpsi")->setVal(3.097);
			Ws.var("sigma1Jpsi")->setVal(0.03);
			Ws.var("sigmaRatioJpsi")->setVal(2.0);
			Ws.var("alphaJpsi")->setVal(1.85);
			Ws.var("nJpsi")->setVal(1.7);
			Ws.var("a1Jpsi")->setVal(-0.10);
			//Ws.var("a2Jpsi")->setVal(0.0);
			Ws.var("ratioJpsi")->setVal(0.5);
			Ws.var("nsigJpsi")->setVal(50000);
			Ws.var("nbkgJpsi")->setVal(2000);
		}
	}
	
	else
	{
		if (pdfName.compare("DoubleSidedCB_BPHTreshold") == 0)
		{
			Ws.var("c2ratio")->removeError();
			//Ws.var("sigma1")->removeError();
			Ws.var("alpha_bkg")->removeError();
			Ws.var("beta_bkg")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			//Ws.var("sigma1")->setVal(0.005);
			Ws.var("alpha_bkg")->setVal(0.2);
			Ws.var("beta_bkg")->setVal(-2);
			Ws.var("nsig")->setVal(500);
			Ws.var("nbkg")->setVal(2000);
		}
		else if (pdfName.compare("DoubleSidedCB_D0BG") == 0)
		{
			Ws.var("c2ratio")->removeError();
			//Ws.var("sigma1")->removeError();
			Ws.var("A_bkg")->removeError();
			//Ws.var("B_bkg")->removeError();
			Ws.var("C_bkg")->removeError();
			Ws.var("massOffset")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			//Ws.var("sigma1")->setVal(0.005);
			Ws.var("A_bkg")->setVal(0.04);
			//Ws.var("B_bkg")->setVal(-2);
			Ws.var("C_bkg")->setVal(0);
			Ws.var("massOffset")->setVal(3.15);
			Ws.var("nsig")->setVal(500);
			Ws.var("nbkg")->setVal(2000);
		}
		else if (pdfName.compare("Hypatia_BPHTreshold") == 0)
		{
			Ws.var("c2ratio")->removeError();
			//Ws.var("sigma1")->removeError();
			Ws.var("alpha_bkg")->removeError();
			Ws.var("beta_bkg")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			//Ws.var("sigma1")->setVal(0.005);
			Ws.var("alpha_bkg")->setVal(0.2);
			Ws.var("beta_bkg")->setVal(-2);
			Ws.var("nsig")->setVal(500);
			Ws.var("nbkg")->setVal(2000);
		}
		else if (pdfName.compare("SingleCB_BPHTreshold") == 0)
		{
			Ws.var("c2ratio")->removeError();
			Ws.var("sigma1")->removeError();
			Ws.var("alpha_bkg")->removeError();
			Ws.var("beta_bkg")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			Ws.var("sigma1")->setVal(0.01);
			Ws.var("alpha_bkg")->setVal(0.2);
			Ws.var("beta_bkg")->setVal(-2);
			Ws.var("nsig")->setVal(500);
			Ws.var("nbkg")->setVal(2000);
		}
		else if (pdfName.compare("SingleCB_D0BG") == 0)
		{
			Ws.var("c2ratio")->removeError();
			Ws.var("sigma1")->removeError();
			Ws.var("A_bkg")->removeError();
			//Ws.var("B_bkg")->removeError();
			Ws.var("C_bkg")->removeError();
			Ws.var("massOffset")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			Ws.var("sigma1")->setVal(0.005);
			Ws.var("A_bkg")->setVal(0.04);
			//Ws.var("B_bkg")->setVal(-2);
			Ws.var("C_bkg")->setVal(0);
			Ws.var("massOffset")->setVal(3.15);
			Ws.var("nsig")->setVal(500);
			Ws.var("nbkg")->setVal(2000);
		}

		else if (pdfName.compare("nominalPdfFlatBkg") == 0)
		{
			Ws.var("c2ratio")->removeError();
			Ws.var("sigma1")->removeError();
			Ws.var("a1")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			Ws.var("sigma1")->setVal(0.005);
			Ws.var("a1")->setVal(0);
			Ws.var("nsig")->setVal(50);
			Ws.var("nbkg")->setVal(2000);
		}
	}

	cout << " Done" << endl;


	return true;
}


int SetConstraints(RooWorkspace& Ws, double* bins, int bin, string fittingSetName = "", string pdfName = "", const char* fileConstraints = "") {
	TFile* fileConstr = new TFile(fileConstraints, "READ"); //opening the file for every fit is probably not the most efficient, but fitting itself is much more time consuming, so it is okay
	if (fileConstr->IsOpen() == false) {
		cout << "ERROR: Constraint file not properly opened!" << endl << endl;
		return -1;
	}
	bool isInTheHeader = false;
	for (int k = 0; k < nFittingSets; k++) { //check if we have the fitting set in the header
		if (fittingSets[k] == fittingSetName)
		{
			isInTheHeader = true;
			break;
		}
	}
	if (isInTheHeader == false) { cout << "ERROR: This sets of bins wasn't found, check fittingSetName and fittingSets definition in ChiFitterInit.h     STOPPING HERE, THIS SET WILL BE SKIPPED " << endl; return -1; }


	if (pdfName.compare("DoubleSidedCB_BPHTreshold") == 0 || pdfName.compare("DoubleSidedCB_D0BG") == 0)
	{
		TGraphAsymmErrors* gAsChic_sigma1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_sigma1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_mean1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_mean1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_alpha = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_alpha_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_n = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_n_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_alphaH = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_alphaH_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_nH = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_nH_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_mean2 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_mean2_" + fittingSetName).c_str());

		double xVal, yVal;
		gAsChic_sigma1->GetPoint(bin, xVal, yVal);
		if (std::abs(xVal - (bins[bin] + bins[bin + 1]) / 2.0) > 0.01) cout << endl << "WARNING: x values of binning between fit and constraints don't agree! " << xVal << "   " << (bins[bin] + bins[bin + 1]) / 2.0 << endl << endl;
		Ws.var("sigma1")->setVal(yVal);
		gAsChic_mean1->GetPoint(bin, xVal, yVal);
		Ws.var("mean1")->setVal(yVal);
		gAsChic_alpha->GetPoint(bin, xVal, yVal);
		Ws.var("alpha")->setVal(yVal);
		gAsChic_n->GetPoint(bin, xVal, yVal);
		Ws.var("n")->setVal(yVal);
		gAsChic_alphaH->GetPoint(bin, xVal, yVal);
		Ws.var("alphaH")->setVal(yVal);
		gAsChic_nH->GetPoint(bin, xVal, yVal);
		Ws.var("nH")->setVal(yVal);
		gAsChic_mean2->GetPoint(bin, xVal, yVal);
		Ws.var("mean2")->setVal(yVal);

		Ws.var("alpha")->setConstant(true);
		Ws.var("n")->setConstant(true);
		Ws.var("alphaH")->setConstant(true);
		Ws.var("nH")->setConstant(true);
		Ws.var("sigma1")->setConstant(true);
		Ws.var("mean1")->setConstant(true);
		Ws.var("mean2")->setConstant(true);
	}
	if (pdfName.compare("Hypatia_BPHTreshold") == 0)
	{
		TGraphAsymmErrors* gAsChic_sigma1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_sigma1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_mean1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_mean1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_aSig1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_aSig1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_nSig1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_nSig1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_aSig2 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_aSig2_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_nSig2 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_nSig2_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_lambda = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_lambda_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_beta = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_beta_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_mean2 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_mean2_" + fittingSetName).c_str());

		double xVal, yVal;
		gAsChic_sigma1->GetPoint(bin, xVal, yVal);
		if (std::abs(xVal - (bins[bin] + bins[bin + 1]) / 2.0) > 0.01) cout << endl << "WARNING: x values of binning between fit and constraints don't agree! " << xVal << "   " << (bins[bin] + bins[bin + 1]) / 2.0 << endl << endl;
		Ws.var("sigma1")->setVal(yVal);
		gAsChic_mean1->GetPoint(bin, xVal, yVal);
		Ws.var("mean1")->setVal(yVal);
		gAsChic_aSig1->GetPoint(bin, xVal, yVal);
		Ws.var("aSig1")->setVal(yVal);
		gAsChic_nSig1->GetPoint(bin, xVal, yVal);
		Ws.var("nSig1")->setVal(yVal);
		gAsChic_aSig2->GetPoint(bin, xVal, yVal);
		Ws.var("aSig2")->setVal(yVal);
		gAsChic_nSig2->GetPoint(bin, xVal, yVal);
		Ws.var("nSig2")->setVal(yVal);
		gAsChic_lambda->GetPoint(bin, xVal, yVal);
		Ws.var("lambda")->setVal(yVal);
		gAsChic_beta->GetPoint(bin, xVal, yVal);
		Ws.var("beta")->setVal(yVal);
		gAsChic_mean2->GetPoint(bin, xVal, yVal);
		Ws.var("mean2")->setVal(yVal);

		Ws.var("sigma1")->setConstant(true);
		Ws.var("mean1")->setConstant(true);
		Ws.var("aSig1")->setConstant(true);
		Ws.var("nSig1")->setConstant(true);
		Ws.var("aSig2")->setConstant(true);
		Ws.var("nSig2")->setConstant(true);
		Ws.var("lambda")->setConstant(true);
		Ws.var("beta")->setConstant(true);
		Ws.var("mean2")->setConstant(true);
	}

	if (pdfName.compare("SingleCB_BPHTreshold") == 0 || pdfName.compare("SingleCB_D0BG") == 0)
	{
		TGraphAsymmErrors* gAsChic_sigma1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_sigma1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_mean1 = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_mean1_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_alpha = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_alpha_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_n = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_n_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_sigmaRatio = (TGraphAsymmErrors*)fileConstr->Get(("gAsChic_Output_sigmaRatio_" + fittingSetName).c_str());

		double xVal, yVal;
		gAsChic_sigma1->GetPoint(bin, xVal, yVal);
		if (std::abs(xVal - (bins[bin] + bins[bin + 1]) / 2.0) > 0.01) cout << endl << "WARNING: x values of binning between fit and constraints don't agree! " << xVal << "   " << (bins[bin] + bins[bin + 1]) / 2.0 << endl << endl;
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




int SetNominalShape(RooWorkspace& Ws, double* bins, int bin, string fittingSetName = "", string pdfName = "", const char* fileInNominal = "") {

	TFile* fNominal = new TFile(fileInNominal, "READ");
	if (fNominal->IsOpen() == false) {
		cout << "ERROR: Nominal file not properly opened!" << endl << endl;
		return -1;
	}


	if (pdfName.compare("DoubleSidedCB_BPHTreshold") == 0 || pdfName.compare("Hypatia_BPHTreshold") == 0)
	{
		TGraphAsymmErrors* gAsChic_nsig = (TGraphAsymmErrors*)fNominal->Get(("gAsChic_FitOutput_nsig_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_c2ratio = (TGraphAsymmErrors*)fNominal->Get(("gAsChic_FitOutput_c2ratio_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_nbkg = (TGraphAsymmErrors*)fNominal->Get(("gAsChic_FitOutput_nbkg_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_alpha_bkg = (TGraphAsymmErrors*)fNominal->Get(("gAsChic_FitOutput_alpha_bkg_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsChic_beta_bkg = (TGraphAsymmErrors*)fNominal->Get(("gAsChic_FitOutput_beta_bkg_" + fittingSetName).c_str());
		if(gAsChic_nsig == nullptr|| gAsChic_c2ratio == nullptr || gAsChic_nbkg == nullptr || gAsChic_alpha_bkg == nullptr || gAsChic_beta_bkg == nullptr) cout << "ERROR, the graphs in setting nominal shape didn't load properly, will crash (maybe bad names, or not in file?)" << endl << "SetName: "<< fittingSetName <<endl<<endl;
		
		double xVal, yVal;
		gAsChic_nsig->GetPoint(bin, xVal, yVal);
		if (std::abs(xVal - (bins[bin] + bins[bin + 1]) / 2.0) > 0.01) cout << endl << "WARNING: x values of binning between fit and nominal setup don't agree! " << xVal << "   " << (bins[bin] + bins[bin + 1]) / 2.0 << endl << endl;
		Ws.var("nsig")->setVal(yVal);
		gAsChic_c2ratio->GetPoint(bin, xVal, yVal);
		Ws.var("c2ratio")->setVal(yVal);
		gAsChic_nbkg->GetPoint(bin, xVal, yVal);
		Ws.var("nbkg")->setVal(yVal);
		gAsChic_alpha_bkg->GetPoint(bin, xVal, yVal);
		Ws.var("alpha_bkg")->setVal(yVal);
		gAsChic_beta_bkg->GetPoint(bin, xVal, yVal);
		Ws.var("beta_bkg")->setVal(yVal);

	}
	else if (pdfName.compare("JpsiTwoCB_Exp") == 0) {
		TGraphAsymmErrors* gAsJpsi_mean1Jpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_mean1Jpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_sigma1Jpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_sigma1Jpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_alphaJpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_alphaJpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_nJpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_nJpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_sigmaRatioJpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_sigmaRatioJpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_ratioJpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_ratioJpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_e_1Jpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_e_1Jpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_nsigJpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_nsigJpsi_" + fittingSetName).c_str());
		TGraphAsymmErrors* gAsJpsi_nbkgJpsi = (TGraphAsymmErrors*)fNominal->Get(("gAsJpsi_FitOutput_nbkgJpsi_" + fittingSetName).c_str());
		if (gAsJpsi_mean1Jpsi == nullptr || gAsJpsi_sigma1Jpsi == nullptr || gAsJpsi_alphaJpsi == nullptr || gAsJpsi_nJpsi == nullptr || gAsJpsi_sigmaRatioJpsi == nullptr || gAsJpsi_ratioJpsi == nullptr || gAsJpsi_e_1Jpsi == nullptr || gAsJpsi_nsigJpsi == nullptr || gAsJpsi_nbkgJpsi == nullptr) cout << "ERROR, the graphs in setting nominal shape didn't load properly, will crash (maybe bad names, or not in file?)" << endl << "SetName: " << fittingSetName << endl << endl;

		double xVal, yVal;
		gAsJpsi_mean1Jpsi->GetPoint(bin, xVal, yVal);
		if (std::abs(xVal - (bins[bin] + bins[bin + 1]) / 2.0) > 0.01) cout << endl << "WARNING: x values of binning between fit and nominal setup don't agree! " << xVal << "   " << (bins[bin] + bins[bin + 1]) / 2.0 << endl << endl;
		Ws.var("mean1Jpsi")->setVal(yVal);
		gAsJpsi_sigma1Jpsi->GetPoint(bin, xVal, yVal);
		Ws.var("sigma1Jpsi")->setVal(yVal);
		gAsJpsi_alphaJpsi->GetPoint(bin, xVal, yVal);
		Ws.var("alphaJpsi")->setVal(yVal);
		gAsJpsi_nJpsi->GetPoint(bin, xVal, yVal);
		Ws.var("nJpsi")->setVal(yVal);
		gAsJpsi_sigmaRatioJpsi->GetPoint(bin, xVal, yVal);
		Ws.var("sigmaRatioJpsi")->setVal(yVal);
		gAsJpsi_ratioJpsi->GetPoint(bin, xVal, yVal);
		Ws.var("ratioJpsi")->setVal(yVal);
		gAsJpsi_e_1Jpsi->GetPoint(bin, xVal, yVal);
		Ws.var("e_1Jpsi")->setVal(yVal);
		gAsJpsi_nsigJpsi->GetPoint(bin, xVal, yVal);
		Ws.var("nsigJpsi")->setVal(yVal);
		gAsJpsi_nbkgJpsi->GetPoint(bin, xVal, yVal);
		Ws.var("nbkgJpsi")->setVal(yVal);

	}
	else {
		cout << " ERROR: MODEL CREATION FAILED, POSSIBLY UNDEFINED MODEL WAS CHOSEN" << endl; return -1;
	}






	fNominal->Close();
	return 0;
}



int DetermineSystematicsBin(double* bins, int binN, bool isSingleBin=false, string fittingSetName = "", int nBinSet=0, const char* fileInNominal = "") { //uses global variables

	RooAbsReal::defaultIntegratorConfig()->setEpsAbs(5e-9); //shouldn't really be here, trying to prevent slow-down
	RooAbsReal::defaultIntegratorConfig()->setEpsRel(5e-9);


	//create the directory:
	string saveDirName = (string)"Output_systJpsi/"+fittingSetName + Form("/bin_%i", binN);
	cout << saveDirName << endl;
	string cmd = string("mkdir ") + string(saveDirName);
	system(cmd.c_str());

	TH1D* hToyDifference_sig = new TH1D("hToyDifference_sig", "Alt signal to nominal - percent difference", 200, -10, 10);
	TH1D* hToyDifference_sig_ZoomedOut = new TH1D("hToyDifference_sig_ZoomedOut", "Alt signal to nominal - percent difference", 100, -50, 50);
	int nFailedSig = 0;

	TH1D* hToyDifference_bkg = new TH1D("hToyDifference_bkg", "Alt background to nominal - percent difference", 200, -10, 10);
	TH1D* hToyDifference_bkg_ZoomedOut = new TH1D("hToyDifference_bkg_ZoomedOut", "Alt background to nominal - percent difference", 100, -50, 50);
	int nFailedBkg = 0;

	

	for (int iToy = 0; iToy < 100; iToy++)
	{
		cout << "Got to iToy " <<iToy<< endl;
		


		//////////////////////////////////////////////////////////////////////////////////

		// WE DEFINE THE WORKSPACE HERE, AS A WORKAROUND, SINCE THERE WAS AN ISSUE WITH THE CODE SLOW-DOWN DUE TO STUFF BEING SAVED SOMEWHERE AND BLOATING (COULDN'T BE TRACED AND DELETED)

		/////////////////////////////////////////////////////////////////////////////////

		RooWorkspace myWs("myWs", "main workspace");
		RooWorkspace myWsNom("myWsNom", "Nominal workspace"); //Maybe not the most elegant, but this let's us to use same names for chic1 and chic2, mean, etc variables without distinguishing among them explicitly (e.g. chic_AltSig)
		RooWorkspace myWsSig("myWsSig", "Alt signal workspace"); //Maybe not the most elegant, but this let's us to use same names for chic1 and chic2, mean, etc variables without distinguishing among them explicitly (e.g. chic_AltSig)
		RooWorkspace myWsBkg("myWsBkg", "Alt bkg workspace");

		//define common variables we use
		RooRealVar* rvmassJpsi = new RooRealVar("rvmassJpsi", "Inv mass #mu#mu", mass_windowJpsi_l, mass_windowJpsi_h, "GeV/c^{2}");
		RooRealVar* rvpt = new RooRealVar("rvpt", "#mu#mu p_{T}", 0.0, 50.0, "GeV/c");
		RooRealVar* rvrap = new RooRealVar("rvrap", "#mu#mu y", -2.4, 2.4, "");
		RooRealVar* rvntrack = new RooRealVar("rvntrack", "ntrack", 0.0, 500.0, "");
		RooArgSet*  cols = new RooArgSet(*rvmassJpsi, *rvpt, *rvrap, *rvntrack);

		myWs.import(*cols);
		myWsNom.import(*cols);
		myWsSig.import(*cols);
		myWsBkg.import(*cols);


		CreateModelPdf(myWs, myPdfName, false);
		CreateModelPdf(myWsNom, myPdfName, false);
		CreateModelPdf(myWsSig, myPdfNameSig, false);
		CreateModelPdf(myWsBkg, myPdfNameBkg, false);

		int flagNominalShape = SetNominalShape(myWs, bins, binN, fittingSetName, myPdfName, fileInNominal);
		if (flagNominalShape == -1) cout << "Setting nominal params failed " << endl;

		TCanvas *cFit = new TCanvas("cFit", "cFit", 1000, 800);
		//gPad->SetLeftMargin(0.15);
		cFit->cd();
		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.275, 0.98, 1.0);
		pad1->SetTicks(1, 1);


		//pull pad
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.006, 0.98, 0.360);
		pad2->SetTopMargin(0); // Upper and lower plot are joined
		pad2->SetBottomMargin(0.67);
		pad2->SetTicks(1, 1);


		RefreshModel(myWsNom, myPdfName, true);
		RefreshModel(myWsSig, myPdfNameSig, true);
		RefreshModel(myWsBkg, myPdfNameBkg, true);




		pad1->cd();

		RooPlot* massframeBin;
		massframeBin = rvmassJpsi->frame(mass_windowFitJpsi_l, mass_windowFitJpsi_h, nMassBinsJpsi);
		massframeBin->SetTitle("mass");

		RooDataSet* rdsToy = myWs.pdf(myPdfName.c_str())->generate(*rvmassJpsi, Extended(kTRUE));
		rdsToy->SetName(TString::Format("rdsToy%i", iToy));
		//rdsToy->SetName(TString::Format("rdsToy%i", 0));
		rdsToy->plotOn(massframeBin, Name("ToyData"));
		//myWs.import(*rdsToy, TString::Format("rdsToy%i", iToy), kTRUE);
		//myWsNom.import(*rdsToy, TString::Format("rdsToy%i", iToy), kTRUE);
		//myWs.import(*rdsToy);
		//myWsNom.import(*rdsToy);
		//myWsSig.import(*rdsToy);

		cout << "Toy generated" << endl;

		myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Name("FullPdf"), LineStyle(kDashed), LineColor(kGray));
		myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("backgroundJpsi"), LineStyle(kDashed), LineColor(kGray));
		myWs.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("signalJpsi"), LineStyle(kDashed), LineColor(kGray));

		cout << "here" << endl;

		//copy initial tail values for better convergence - not fixed, just initial seeds
		myWsNom.var("nJpsi")->setVal(myWs.var("nJpsi")->getValV());
		myWsNom.var("alphaJpsi")->setVal(myWs.var("alphaJpsi")->getValV());
		myWsNom.var("e_1Jpsi")->setVal(myWs.var("e_1Jpsi")->getValV());

		RooFitResult* fitResultNom = myWsNom.pdf(myPdfName.c_str())->fitTo(*rdsToy, Extended(true), SumW2Error(true), PrintLevel(-1), Range(mass_windowFitJpsi_l, mass_windowFitJpsi_h), NumCPU(1), Save(true));
		myWsNom.pdf(myPdfName.c_str())->plotOn(massframeBin, Name("FullPdfNom"), LineColor(kBlue));
		myWsNom.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("backgroundJpsi"), LineStyle(kDashed), LineColor(kBlue));
		myWsNom.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("signalJpsi"), LineStyle(kDashed), LineColor(kBlue));

		RooArgList paramList = fitResultNom->floatParsFinal();
		RooRealVar *par;
		for (int j = 0; j < paramList.getSize(); j++)
		{
			par = (RooRealVar*)paramList.at(j);
			cout << par->getTitle() << " = " << par->getValV() << " + " << par->getErrorHi() << "  " << par->getErrorLo() << endl;
		}

		//cout << "NomFitStatus: " << fitResultNom->status() << endl;

		double sigNomCount = myWsNom.var("nsigJpsi")->getValV();
		cout << "signal nominal " << sigNomCount << endl;
		
		//cout << endl;
		//RooAddPdf::pdfList();

		//myWs.pdf(myPdfName.c_str())->fitTo(*rdsToy);
		RooFitResult* fitResultSig = myWsSig.pdf(myPdfNameSig.c_str())->fitTo(*rdsToy, Extended(true), SumW2Error(true), PrintLevel(-1), Range(mass_windowFitJpsi_l, mass_windowFitJpsi_h), NumCPU(1), Save(true));
		//RooFitResult* fitResult = myWs.pdf(myPdfName.c_str())->fitTo(*myWs.data("rdsToy1"));//, Extended(true), SumW2Error(true), Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true));
		////fitResult->Print("v");
		  //  myWs.import(*fitResult, "fitResult");
		
		
		myWsSig.pdf(myPdfNameSig.c_str())->plotOn(massframeBin, Name("FullPdfSig"), LineColor(kRed));
		//massframeBin->getAttText()->SetTextSize(0.03);
		//myWsSig.pdf(myPdfNameSig.c_str())->paramOn(massframeBin, Layout(0.70,0.99));
		myWsSig.pdf(myPdfNameSig.c_str())->plotOn(massframeBin, Components("backgroundJpsi"), LineStyle(kDashed), LineColor(kRed));
		myWsSig.pdf(myPdfNameSig.c_str())->plotOn(massframeBin, Components("signalJpsi"), LineStyle(kDashed), LineColor(kRed));

		RooArgList paramListSig = fitResultSig->floatParsFinal();
		RooRealVar* parSig;
		for (int j = 0; j < paramListSig.getSize(); j++)
		{
			parSig = (RooRealVar*)paramListSig.at(j);
			cout << parSig->getTitle() << " = " << parSig->getValV() << " + " << parSig->getErrorHi() << "  " << parSig->getErrorLo() << endl;
		}

		////cout << "SigFitStatus: " << fitResultSig->status() << endl;
		double sigAltCount = myWsSig.var("nsigJpsi")->getValV();
		cout << "signal alternate " << sigAltCount << endl;

		double percentDiff_sig = (sigAltCount - sigNomCount) / sigNomCount * 100;
		cout << "Percent difference signal: " << percentDiff_sig << endl;


		//if (fitResultNom->status() != 0 || fitResultSig->status() != 0)  //moved after we have chi2
		//{
		//	nFailedSig++;
		//	cout << "A FIT HAS FAILED, RESULT SKIPPED" << endl;
		//}
		//else
		//{
		//	hToyDifference_sig->Fill(percentDiff_sig);
		//	hToyDifference_sig_ZoomedOut->Fill(percentDiff_sig);
		//}

		//myWsBkg.pdf(myPdfNameBkg.c_str())->fitTo(*myWs.data("rdsToy1"));

		///////  BACKGROUND  

		//fix the tail parameters to nominal, otherwise they tend to get messed up and we won't probe background variation
		myWsBkg.var("nJpsi")->setConstant(true);
		myWsBkg.var("nJpsi")->setVal(myWsNom.var("nJpsi")->getValV());
		myWsBkg.var("alphaJpsi")->setConstant(true);
		myWsBkg.var("alphaJpsi")->setVal(myWsNom.var("alphaJpsi")->getValV());
		cout << "Background fitting, constrained parameters nJpsi to " << myWsBkg.var("nJpsi")->getValV() << " from nominal fit: " << myWsNom.var("nJpsi")->getValV() << " and alphaJpsi to " << myWsBkg.var("alphaJpsi")->getValV() << " from nominal fit: " << myWsNom.var("alphaJpsi")->getValV() << endl;

		RooFitResult* fitResultBkg = myWsBkg.pdf(myPdfNameBkg.c_str())->fitTo(*rdsToy, Extended(true), SumW2Error(true), PrintLevel(-1), Range(mass_windowFitJpsi_l, mass_windowFitJpsi_h), NumCPU(1), Save(true));
		myWsBkg.pdf(myPdfNameBkg.c_str())->plotOn(massframeBin, Name("FullPdfBkg"), LineColor(kGreen+2));
		//myWsBkg.pdf(myPdfNameBkg.c_str())->paramOn(massframeBin, Layout(0.05, 0.35));
		myWsBkg.pdf(myPdfNameBkg.c_str())->plotOn(massframeBin, Components("backgroundJpsi"), LineStyle(kDashed), LineColor(kGreen+2));
		myWsBkg.pdf(myPdfNameBkg.c_str())->plotOn(massframeBin, Components("signalJpsi"), LineStyle(kDashed), LineColor(kGreen+2));


		double bkgAltCount = myWsBkg.var("nsigJpsi")->getValV();
		cout << "Background alternate " << bkgAltCount << endl;

		double percentDiff_bkg = (bkgAltCount - sigNomCount) / sigNomCount * 100;
		cout << "Percent difference background: " << percentDiff_bkg << endl;
		

		//if (fitResultNom->status() != 0 || fitResultBkg->status() != 0) //moved after we have chi2
		//{
		//	nFailedBkg++;
		//	cout << "A FIT HAS FAILED, RESULT SKIPPED" << endl;
		//}
		//else
		//{
		//	hToyDifference_bkg->Fill(percentDiff_bkg);
		//	hToyDifference_bkg_ZoomedOut->Fill(percentDiff_bkg);
		//}


		massframeBin->Draw();

		///////////
		//counts and percent difference
		/////////////

		TPaveText *textchiCounts = new TPaveText(.78, .72, .95, .92, "brNDC");
		textchiCounts->SetTextAlign(31);
		textchiCounts->SetFillColor(kWhite);
		textchiCounts->SetTextColor(kGreen+2);
		textchiCounts->SetTextSize(0.04);
		textchiCounts->SetBorderSize(0);
		//char* txtchiCounts1a = Form("Counts nominal: %.1f, alt sig: %.1f", sigNomCount, sigAltCount);
		char* txtchiCounts1a = Form("Counts nominal: %.1f", sigNomCount);
		char* txtchiCounts2a = Form("  Counts alt sig: %.1f", sigAltCount);
		char* txtchiCounts2b = Form("Percent difference %.2f", percentDiff_sig);
		char* txtchiCounts3a = Form("  Counts alt bkg: %.1f", bkgAltCount);
		char* txtchiCounts3b = Form("Percent difference %.2f", percentDiff_bkg);
		TText* t = textchiCounts->AddText(txtchiCounts1a);
		t->Draw();
		t->SetTextColor(kBlue);
		TText* t2 = textchiCounts->AddText(txtchiCounts2a);
		t2->SetTextColor(kRed);
		t2 = textchiCounts->AddText(txtchiCounts2b);
		t2->Draw();
		t2->SetTextColor(kRed);
		textchiCounts->AddText(txtchiCounts3a);
		textchiCounts->AddText(txtchiCounts3b);
		textchiCounts->Draw();



//		chi2

		int numFitParNom = fitResultNom->floatParsFinal().getSize();
		//cout << "Npar: " << numFitPar << endl;
		double chi2_fitNom = massframeBin->chiSquare("FullPdfNom", "ToyData", numFitParNom);
		TPaveText *textchiNom = new TPaveText(.23, .85, .38, .92, "brNDC");
		textchiNom->SetFillColor(kWhite);
		textchiNom->SetTextColor(kBlue);
		textchiNom->SetTextSize(0.05);
		textchiNom->SetBorderSize(0);
		//string txtchiNom = "#chi^{2}/ndf: " + to_string(chi2_fitNom);
		string txtchiNom = "#chi^{2}/ndf: " + to_string(chi2_fitNom).substr(0, std::to_string(chi2_fitNom).find(".") + 3);// output to string, to 2 dec places (truncated and not rounded, but that should not matter)
		textchiNom->AddText(txtchiNom.c_str());
		textchiNom->Draw();

		if (chi2_fitNom > 5 || fitResultNom->status() != 0) {
			TPaveText *textErrorNom = new TPaveText(.35, .25, .80, .45, "brNDC");
			textErrorNom->SetFillColor(0);
			textErrorNom->SetTextColorAlpha(2, 0.1);
			textErrorNom->SetTextSize(0.30);
			textErrorNom->SetBorderSize(0);
			textErrorNom->AddText("FAILED");
			textErrorNom->Draw();
		}

		int numFitParSig = fitResultSig->floatParsFinal().getSize();
		//cout << "Npar: " << numFitPar << endl;
		double chi2_fitSig = massframeBin->chiSquare("FullPdfSig", "ToyData", numFitParSig);
		TPaveText *textchiSig = new TPaveText(.23, .77, .38, .84, "brNDC");
		textchiSig->SetFillColor(kWhite);
		textchiSig->SetTextColor(kRed);
		textchiSig->SetTextSize(0.05);
		textchiSig->SetBorderSize(0);
		//string txtchiSig = "#chi^{2}/ndf: " + to_string(chi2_fitSig);
		string txtchiSig = "#chi^{2}/ndf: " + to_string(chi2_fitSig).substr(0, std::to_string(chi2_fitSig).find(".") + 3);// output to string, to 2 dec places (truncated and not rounded, but that should not matter)
		textchiSig->AddText(txtchiSig.c_str());
		textchiSig->Draw();

		if (chi2_fitSig > 5 || fitResultSig->status() != 0) {
			TPaveText *textErrorSig = new TPaveText(.35, .25, .80, .45, "brNDC");
			textErrorSig->SetFillColor(0);
			textErrorSig->SetTextColorAlpha(2, 0.1);
			textErrorSig->SetTextSize(0.30);
			textErrorSig->SetBorderSize(0);
			textErrorSig->AddText("FAILED");
			textErrorSig->Draw();
		}

		int numFitParBkg = fitResultBkg->floatParsFinal().getSize();
		//cout << "Npar: " << numFitPar << endl;
		double chi2_fitBkg = massframeBin->chiSquare("FullPdfBkg", "ToyData", numFitParBkg);
		TPaveText *textchiBkg = new TPaveText(.23, .69, .38, .76, "brNDC");
		textchiBkg->SetFillColor(kWhite);
		textchiBkg->SetTextColor(kGreen+2);
		textchiBkg->SetTextSize(0.05);
		textchiBkg->SetBorderSize(0);
		//string txtchiBkg = "#chi^{2}/ndf: " + to_string(chi2_fitBkg);
		string txtchiBkg = "#chi^{2}/ndf: " + to_string(chi2_fitBkg).substr(0, std::to_string(chi2_fitBkg).find(".") + 3);// output to string, to 2 dec places (truncated and not rounded, but that should not matter)
		textchiBkg->AddText(txtchiBkg.c_str());
		textchiBkg->Draw();

		if (chi2_fitBkg > 5 || fitResultBkg->status() != 0) {
			TPaveText *textErrorBkg = new TPaveText(.35, .25, .80, .45, "brNDC");
			textErrorBkg->SetFillColor(0);
			textErrorBkg->SetTextColorAlpha(2, 0.1);
			textErrorBkg->SetTextSize(0.30);
			textErrorBkg->SetBorderSize(0);
			textErrorBkg->AddText("FAILED");
			textErrorBkg->Draw();
		}

		//PULLS

		pad2->cd();
		//RooHist* hpull = massframeBin->pullHist("rPullHist", myPdfName.c_str());
		RooHist* hpullNom = massframeBin->pullHist("ToyData", "FullPdfNom");
		hpullNom->SetMarkerSize(0.8);
		hpullNom->SetMarkerColor(kBlue);
		RooHist* hpullSig = massframeBin->pullHist("ToyData", "FullPdfSig");
		hpullSig->SetMarkerSize(0.8);
		hpullSig->SetMarkerColor(kRed);
		RooHist* hpullBkg = massframeBin->pullHist("ToyData", "FullPdfBkg");
		hpullBkg->SetMarkerSize(0.8);
		hpullBkg->SetMarkerColor(kGreen+2);
		RooPlot* pullFrame = rvmassJpsi->frame(Title("Pull Distribution"));
		pullFrame->addPlotable(hpullNom, "P");
		pullFrame->addPlotable(hpullSig, "P");
		pullFrame->addPlotable(hpullBkg, "P");
		pullFrame->SetTitleSize(0);
		pullFrame->GetYaxis()->SetTitle("Pull");
		pullFrame->GetYaxis()->SetTitleSize(0.19);
		pullFrame->GetYaxis()->SetLabelSize(0.14);
		//pullFrame->GetYaxis()->SetLabelOffset(1.1);
		pullFrame->GetYaxis()->SetRangeUser(-5.0, 5.0);
		pullFrame->GetYaxis()->SetNdivisions(502, kTRUE);
		pullFrame->GetYaxis()->CenterTitle();
		pullFrame->GetXaxis()->SetLabelSize(0.15);
		pullFrame->GetXaxis()->SetTitle("M (#mu#mu#gamma - #mu#mu + 3.097) [GeV/c^{2}]"); 
		pullFrame->GetXaxis()->SetTitleOffset(1.05);
		pullFrame->GetXaxis()->SetTitleSize(0.15);
		pullFrame->Draw();

		TLine *l1 = new TLine(mass_windowFit_l, 0, mass_windowFit_h, 0);
		l1->SetLineStyle(9);
		l1->Draw("same");


		cout << "Done with pulls"<< endl;
		

		//TObject* testObj = FindObject("DoubleSidedCB_BPHTreshold_Int[rvmassJpsi|fit_nll_DoubleSidedCB_BPHTreshold_rdsToy0]_Norm[rvmassJpsi]");
		
		cFit->cd();
		pad1->Draw();
		pad2->Draw();

		pad1->Update();
		pad2->Update();

		cFit->SaveAs(((string) saveDirName+"/JpsiToys-" + fittingSetName +"_bin_"+ Form("%i_toy_%i.png",binN, iToy)).c_str());
		cFit->SaveAs(((string)saveDirName + "/JpsiToys-" + fittingSetName + "_bin_" + Form("%i_toy_%i.pdf", binN, iToy)).c_str());
		//cFit->SaveAs(((string)"FitterOutput/FitResult_" + (isJpsi ? "Jpsi_" : "Chic_") + fittingSetName + "_" + binVarName + Form("_%ibin_", i) + Form("_%.1f_%.1f.png", bins[i], bins[i + 1])).c_str());
		cout << endl;



		if (fitResultNom->status() != 0 || fitResultSig->status() != 0 || chi2_fitSig > 5.0 || chi2_fitNom > 5.0)
		{
			nFailedSig++;
			cout << "A FIT HAS FAILED, RESULT SKIPPED" << endl;
		}
		else
		{
			hToyDifference_sig->Fill(percentDiff_sig);
			hToyDifference_sig_ZoomedOut->Fill(percentDiff_sig);
		}


		if (fitResultNom->status() != 0 || fitResultBkg->status() != 0 || chi2_fitBkg > 5.0 || chi2_fitNom > 5.0)
		{
			nFailedBkg++;
			cout << "A FIT HAS FAILED, RESULT SKIPPED" << endl;
		}
		else
		{
			hToyDifference_bkg->Fill(percentDiff_bkg);
			hToyDifference_bkg_ZoomedOut->Fill(percentDiff_bkg);
		}








		//clean up
		delete rdsToy;
		delete massframeBin;

		delete cFit;
		delete fitResultNom;
		delete fitResultSig;
		delete fitResultBkg;
		//delete hpullNom;
		//delete hpullSig;
		//delete hpullBkg;
		delete pullFrame;

	}

	cout << "Done with toys in the bin" << endl;

	TCanvas *cDiff = new TCanvas("cDiff", "cDiff", 800, 600);
	//hToyDifference_sig->SetStats(true);
	hToyDifference_sig->SetLineColor(kRed);
	hToyDifference_bkg->SetLineColor(kGreen+2);
	hToyDifference_sig->Draw();
	//TPaveStats *st = (TPaveStats*)cDiff->GetPrimitive("stats");
	//st->SetName("Stat_sig");
	//st->SetX1NDC(0.2); //new x start position
	//st->SetX2NDC(0.4); //new x end position

	hToyDifference_bkg->Draw("same");
	cDiff->SaveAs(((string)saveDirName + "/JpsiToys" + fittingSetName + Form("Result_bin%i", binN) + ".png").c_str());
	cDiff->SaveAs(((string)saveDirName + "/JpsiToys" + fittingSetName + Form("Result_bin%i", binN) + ".root").c_str());
	cDiff->SaveAs(((string)saveDirName + "/JpsiToys" + fittingSetName + Form("Result_bin%i", binN) + ".pdf").c_str());

	
	TCanvas *cDiff2 = new TCanvas("cDiff2", "cDiff2", 800, 600);
	hToyDifference_bkg_ZoomedOut->SetLineColor(kGreen+2);
	hToyDifference_bkg_ZoomedOut->Draw();
	hToyDifference_sig_ZoomedOut->SetLineColor(kRed);
	hToyDifference_sig_ZoomedOut->Draw("same");
	cDiff2->SaveAs(((string)saveDirName + "/JpsiToys" + fittingSetName + Form("ResultZoomedOut_bin%i", binN) + ".png").c_str());
	cDiff2->SaveAs(((string)saveDirName + "/JpsiToys" + fittingSetName + Form("ResultZoomedOut_bin%i", binN) + ".root").c_str());
	cDiff2->SaveAs(((string)saveDirName + "/JpsiToys" + fittingSetName + Form("ResultZoomedOut_bin%i", binN) + ".pdf").c_str());


	double rmsSig = hToyDifference_sig_ZoomedOut->GetRMS();
	double meanSig = hToyDifference_sig_ZoomedOut->GetMean();

	cout << "Sig RMS: " << rmsSig << " zoomed out value: " << hToyDifference_sig_ZoomedOut->GetRMS() << endl;
	cout << "Sig Mean: " << meanSig << endl << endl;

	double rmsBkg = hToyDifference_bkg_ZoomedOut->GetRMS();
	double meanBkg = hToyDifference_bkg_ZoomedOut->GetMean();

	cout << "Bkg RMS: " << rmsBkg << " zoomed out value: " << hToyDifference_bkg_ZoomedOut->GetRMS() << endl;
	cout << "Bkg Mean: " << meanBkg << endl << endl;

	if (isSingleBin == false) {
		gSystOutputArraySig[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), sqrt(rmsSig*rmsSig + meanSig * meanSig));
		gSystOutputArraySig_rms[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), rmsSig);
		gSystOutputArraySig_mean[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), meanSig);
		gSystOutputArraySig_failed[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), nFailedSig);

		gSystOutputArrayBkg[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), sqrt(rmsBkg*rmsBkg + meanBkg * meanBkg));
		gSystOutputArrayBkg_rms[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), rmsBkg);
		gSystOutputArrayBkg_mean[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), meanBkg);
		gSystOutputArrayBkg_failed[nBinSet]->SetPoint(binN, ((bins[binN] + bins[binN + 1]) / 2.0), nFailedBkg);
	}
	else {
		gSystOutputArraySig[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), sqrt(rmsSig*rmsSig + meanSig * meanSig));
		gSystOutputArraySig_rms[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), rmsSig);
		gSystOutputArraySig_mean[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), meanSig);
		gSystOutputArraySig_failed[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), nFailedSig);

		gSystOutputArrayBkg[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), sqrt(rmsBkg*rmsBkg + meanBkg * meanBkg));
		gSystOutputArrayBkg_rms[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), rmsBkg);
		gSystOutputArrayBkg_mean[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), meanBkg);
		gSystOutputArrayBkg_failed[nBinSet]->SetPoint(0, ((bins[binN] + bins[binN + 1]) / 2.0), nFailedBkg);
	}
	//delete cFit;
	delete cDiff;
	delete cDiff2;
	delete hToyDifference_sig;
	delete hToyDifference_sig_ZoomedOut;
	delete hToyDifference_bkg;
	delete hToyDifference_bkg_ZoomedOut;

	return 0;
	
}

	
int DetermineSystematicsSet(double* bins, int nbins, string fittingSetName = "", int BinToFit=-1, bool bConstrainedFit = false, const char* fileInNominal = "")
{
	string saveDirName = (string)"Output_systJpsi/" + fittingSetName;
	cout << saveDirName << endl;
	string cmd = string("mkdir ") + string(saveDirName);
	system(cmd.c_str());

	int nBinSet = -1;
	for (int k = 0; k < nFittingSets; k++) { //find the corresponding set index for readout
		if (fittingSets[k] == fittingSetName)
		{
			nBinSet = k;
			break;
		}
	}
	if (nBinSet == -1) { cout << "ERROR: This sets of bins wasn't found, check fittingSetName and fittingSets definition in ChiFitterInit.h     STOPPING HERE, THIS SET WILL BE SKIPPED " << endl; return -1; }

	if (BinToFit == -1) {
		gSystOutputArraySig[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - total
		gSystOutputArraySig_mean[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - mean only
		gSystOutputArraySig_rms[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - RMS only
		gSystOutputArraySig_failed[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - RMS only
	}
	else //we are doing only one point
	{
		gSystOutputArraySig[nBinSet] = new TGraph(1); // stores the systematic uncertainty - total
		gSystOutputArraySig_mean[nBinSet] = new TGraph(1); // stores the systematic uncertainty - mean only
		gSystOutputArraySig_rms[nBinSet] = new TGraph(1); // stores the systematic uncertainty - RMS only
		gSystOutputArraySig_failed[nBinSet] = new TGraph(1); // stores the systematic uncertainty - RMS only
	}

	cout << "gSystOutputArraySig_" + fittingSets.at(nBinSet) << endl;
	gSystOutputArraySig[nBinSet]->SetNameTitle(("gSystOutputArraySig_" + fittingSets.at(nBinSet)).c_str(), ("Signal: Total syst uncertainty for " + fittingSets.at(nBinSet)).c_str());
	gSystOutputArraySig_mean[nBinSet]->SetNameTitle(("gSystOutputArraySig_mean_" + fittingSets.at(nBinSet)).c_str(), ("Signal: Mean for " + fittingSets.at(nBinSet)).c_str());
	gSystOutputArraySig_rms[nBinSet]->SetNameTitle(("gSystOutputArraySig_rms_" + fittingSets.at(nBinSet)).c_str(), ("Signal: RMS for " + fittingSets.at(nBinSet)).c_str());
	gSystOutputArraySig_failed[nBinSet]->SetNameTitle(("gSystOutputArraySig_failed_" + fittingSets.at(nBinSet)).c_str(), ("Signal: number of failed fits for " + fittingSets.at(nBinSet)).c_str());

	if (BinToFit == -1) {
		gSystOutputArrayBkg[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - total
		gSystOutputArrayBkg_mean[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - mean only
		gSystOutputArrayBkg_rms[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - RMS only
		gSystOutputArrayBkg_failed[nBinSet] = new TGraph(nbins); // stores the systematic uncertainty - RMS only
	}
	else //we are doing only one point
	{
		gSystOutputArrayBkg[nBinSet] = new TGraph(1); // stores the systematic uncertainty - total
		gSystOutputArrayBkg_mean[nBinSet] = new TGraph(1); // stores the systematic uncertainty - mean only
		gSystOutputArrayBkg_rms[nBinSet] = new TGraph(1); // stores the systematic uncertainty - RMS only
		gSystOutputArrayBkg_failed[nBinSet] = new TGraph(1); // stores the systematic uncertainty - RMS only
	}

	cout << "gSystOutputArrayBkg_" + fittingSets.at(nBinSet) << endl;
	gSystOutputArrayBkg[nBinSet]->SetNameTitle(("gSystOutputArrayBkg_" + fittingSets.at(nBinSet)).c_str(), ("Bkg: Total syst uncertainty for " + fittingSets.at(nBinSet)).c_str());
	gSystOutputArrayBkg_mean[nBinSet]->SetNameTitle(("gSystOutputArrayBkg_mean_" + fittingSets.at(nBinSet)).c_str(), ("Bkg: Mean for " + fittingSets.at(nBinSet)).c_str());
	gSystOutputArrayBkg_rms[nBinSet]->SetNameTitle(("gSystOutputArrayBkg_rms_" + fittingSets.at(nBinSet)).c_str(), ("Bkg: RMS for " + fittingSets.at(nBinSet)).c_str());
	gSystOutputArrayBkg_failed[nBinSet]->SetNameTitle(("gSystOutputArrayBkg_failed_" + fittingSets.at(nBinSet)).c_str(), ("Bkg: number of failed fits for " + fittingSets.at(nBinSet)).c_str());


	if (BinToFit == -1) { // do all
		for (int binN = 0; binN < nbins; binN++)
		{
			DetermineSystematicsBin(bins, binN, false, fittingSetName, nBinSet, fileInNominal);
		}
	}
	else
	{
		if (BinToFit < nbins) { DetermineSystematicsBin(bins, BinToFit, true, fittingSetName, nBinSet, fileInNominal); }
		else { cout << "THIS BIN ISN'T PRESENT, SKIPPED" << endl; } //should be checked before
	}

	return 0;
}
	
	


///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////

//void SystFits_Jpsi(int SetToFit = 0, int BinToFit = -1, const char* fileInNominal = "Chi_c_output_Nominal_v2-bothDir_DCB.root", const char* fileOut = "Chi_c_Jpsi_Syst.root", bool flagConstrainedFit = false)
void SystFits_Jpsi(int SetToFit = 0, int BinToFit = -1, const char* fileInNominal = "Chi_c_output_Nominal_vDissertation-bothDir_DCB_NewBins.root", const char* fileOut = "Chi_c_Jpsi_Syst_vDissertation.root", bool flagConstrainedFit = false)
{
	
	
	//gStyle->SetOptStat(0);
	setTDRStyle();
	
	gStyle->SetOptStat(1111);


	file_log.open("TestLog.txt", ios::app);
	//file_log << endl << endl << "N E W   R U N " << endl << endl;



	////////////////////////////// 
	///   SETUP FOR CONDOR   ////
	/////////////////////////////
	   
	string myFittingSetName = fittingSets[SetToFit];
	
	cout << endl << "GOING TO FIT " << myFittingSetName << endl << endl;

	if (BinToFit == -1) {
		cout << "Fitting all bins for the set" << endl << endl;
	}
	else {
		cout << "Bin: " << BinToFit << endl << endl;
		// skip those jobs that would be off-limit
		if (SetToFit == 1) {
			if (BinToFit >= nbins_y) { cout << "Job would be out of limits, done" << endl; return; }
		}
		//else if (SetToFit == 2 || SetToFit == 12) {
		else if (SetToFit == 5) {
			if (BinToFit >= nbins_nTrk) { cout << "Job would be out of limits, done" << endl; return; }
		}
		else {
			if (BinToFit >= nbins_pT) { cout << "Job would be out of limits, done" << endl; return; }
		}
	}


	string fileOutCon = fileOut;
	fileOutCon.erase(fileOutCon.end() - 5, fileOutCon.end()); //erase .root
	fileOutCon = fileOutCon + Form("_%i_%i.root", SetToFit, BinToFit);
	cout << "File to save is: " << fileOutCon << endl;

	// create the output structure, in case it doesn't exist yet
	string cmd = string("mkdir Output_systJpsi/");
	system(cmd.c_str());

	////////////////////////////////////////  
	///////   R  o  o  f  i  t   //////////
	/////////////////////////////////////////


	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
	ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);


	if (SetToFit==1){ DetermineSystematicsSet(bins_y, nbins_y, myFittingSetName, BinToFit, flagConstrainedFit, fileInNominal); }
	//else if (SetToFit==2 || SetToFit==12){ DetermineSystematicsSet(bins_nTrk, nbins_nTrk, myFittingSetName, BinToFit, flagConstrainedFit, fileInNominal); }
	else if (SetToFit == 5) { DetermineSystematicsSet(bins_nTrk, nbins_nTrk, myFittingSetName, BinToFit, flagConstrainedFit, fileInNominal); }
	else { DetermineSystematicsSet(bins_pT, nbins_pT, myFittingSetName, BinToFit, flagConstrainedFit, fileInNominal); }

	gDirectory->ls();




	TFile* fout = new TFile(fileOutCon.c_str(), "RECREATE");


	for (int i = 0; i < nFittingSets; i++)
	{

		if (gSystOutputArraySig[i] != NULL)
		{

			gSystOutputArraySig[i]->Write();
			gSystOutputArraySig_mean[i]->Write();
			gSystOutputArraySig_rms[i]->Write();
			gSystOutputArraySig_failed[i]->Write();

			gSystOutputArrayBkg[i]->Write();
			gSystOutputArrayBkg_mean[i]->Write();
			gSystOutputArrayBkg_rms[i]->Write();
			gSystOutputArrayBkg_failed[i]->Write();

		}

	}

	fout->Close();

}



RooDoubleCB::RooDoubleCB(const char *name, const char *title,
	RooAbsReal& _x,
	RooAbsReal& _mu,
	RooAbsReal& _sig,
	RooAbsReal& _a1,
	RooAbsReal& _n1,
	RooAbsReal& _a2,
	RooAbsReal& _n2) :
	RooAbsPdf(name, title),
	x("x", "x", this, _x),
	mu("mu", "mu", this, _mu),
	sig("sig", "sig", this, _sig),
	a1("a1", "a1", this, _a1),
	n1("n1", "n1", this, _n1),
	a2("a2", "a2", this, _a2),
	n2("n2", "n2", this, _n2)
{
}


RooDoubleCB::RooDoubleCB(const RooDoubleCB& other, const char* name) :
	RooAbsPdf(other, name),
	x("x", this, other.x),
	mu("mu", this, other.mu),
	sig("sig", this, other.sig),
	a1("a1", this, other.a1),
	n1("n1", this, other.n1),
	a2("a2", this, other.a2),
	n2("n2", this, other.n2)
{
}



Double_t RooDoubleCB::evaluate() const
{
	double u = (x - mu) / sig;
	double A1 = TMath::Power(n1 / TMath::Abs(a1), n1)*TMath::Exp(-a1 * a1 / 2);
	double A2 = TMath::Power(n2 / TMath::Abs(a2), n2)*TMath::Exp(-a2 * a2 / 2);
	double B1 = n1 / TMath::Abs(a1) - TMath::Abs(a1);
	double B2 = n2 / TMath::Abs(a2) - TMath::Abs(a2);

	double result(1);
	if (u < -a1) result *= A1 * TMath::Power(B1 - u, -n1);
	else if (u < a2)  result *= TMath::Exp(-u * u / 2);
	else            result *= A2 * TMath::Power(B2 + u, -n2);
	return result;
}


/// \param[in] n2 Shape parameter of right tail (\f$ n2 \ge 0 \f$). With \f$ n2 = 0 \f$, the function is constant.
RooHypatia2::RooHypatia2(const char *name, const char *title, RooAbsReal& x, RooAbsReal& lambda,
	RooAbsReal& zeta, RooAbsReal& beta, RooAbsReal& sigm, RooAbsReal& mu, RooAbsReal& a,
	RooAbsReal& n, RooAbsReal& a2, RooAbsReal& n2) :
	RooAbsPdf(name, title),
	_x("x", "x", this, x),
	_lambda("lambda", "Lambda", this, lambda),
	_zeta("zeta", "zeta", this, zeta),
	_beta("beta", "Asymmetry parameter beta", this, beta),
	_sigma("sigma", "Width parameter sigma", this, sigm),
	_mu("mu", "Location parameter mu", this, mu),
	_a("a", "Left tail location a", this, a),
	_n("n", "Left tail parameter n", this, n),
	_a2("a2", "Right tail location a2", this, a2),
	_n2("n2", "Right tail parameter n2", this, n2)
{
	//	RooHelpers::checkRangeOfParameters(this, { &sigm }, 0.);
	//	RooHelpers::checkRangeOfParameters(this, { &zeta, &n, &n2, &a, &a2 }, 0., std::numeric_limits<double>::max(), true);
	//	if (zeta.getVal() == 0. && zeta.isConstant()) {
	//		RooHelpers::checkRangeOfParameters(this, { &lambda }, -std::numeric_limits<double>::max(), 0., false,
	//			std::string("Lambda needs to be negative when ") + _zeta.GetName() + " is zero.");
	//	}
	//
	//#ifndef R__HAS_MATHMORE
	//	throw std::logic_error("RooHypatia2 needs ROOT with mathmore enabled to access gsl functions.");
	//#endif
}


///////////////////////////////////////////////////////////////////////////////////////////
/// Copy a new Hypatia2 PDF.
/// \param[in] other Original to copy from.
/// \param[in] name Optional new name.
RooHypatia2::RooHypatia2(const RooHypatia2& other, const char* name) :
	RooAbsPdf(other, name),
	_x("x", this, other._x),
	_lambda("lambda", this, other._lambda),
	_zeta("zeta", this, other._zeta),
	_beta("beta", this, other._beta),
	_sigma("sigma", this, other._sigma),
	_mu("mu", this, other._mu),
	_a("a", this, other._a),
	_n("n", this, other._n),
	_a2("a2", this, other._a2),
	_n2("n2", this, other._n2)
{
#ifndef R__HAS_MATHMORE
	throw std::logic_error("RooHypatia2 needs ROOT with mathmore enabled to access gsl functions.");
#endif
}

namespace {
	const double sq2pi_inv = 1. / std::sqrt(TMath::TwoPi());
	const double logsq2pi = std::log(std::sqrt(TMath::TwoPi()));
	const double ln2 = std::log(2.);

	double low_x_BK(double nu, double x) {
		return TMath::Gamma(nu)*std::pow(2., nu - 1.)*std::pow(x, -nu);
	}

	double low_x_LnBK(double nu, double x) {
		return std::log(TMath::Gamma(nu)) + ln2 * (nu - 1.) - std::log(x) * nu;
	}

	double besselK(double ni, double x) {
		const double nu = std::fabs(ni);
		if ((x < 1.e-06 && nu > 0.) ||
			(x < 1.e-04 && nu > 0. && nu < 55.) ||
			(x < 0.1 && nu >= 55.))
			return low_x_BK(nu, x);

#ifdef R__HAS_MATHMORE
		return ROOT::Math::cyl_bessel_k(nu, x);
#else
		return std::numeric_limits<double>::signaling_NaN();
#endif

	}

	double LnBesselK(double ni, double x) {
		const double nu = std::fabs(ni);
		if ((x < 1.e-06 && nu > 0.) ||
			(x < 1.e-04 && nu > 0. && nu < 55.) ||
			(x < 0.1 && nu >= 55.))
			return low_x_LnBK(nu, x);

#ifdef R__HAS_MATHMORE
		return std::log(ROOT::Math::cyl_bessel_k(nu, x));
#else
		return std::numeric_limits<double>::signaling_NaN();
#endif
	}


	double LogEval(double d, double l, double alpha, double beta, double delta) {
		const double gamma = alpha;//std::sqrt(alpha*alpha-beta*beta);
		const double dg = delta * gamma;
		const double thing = delta * delta + d * d;
		const double logno = l * std::log(gamma / delta) - logsq2pi - LnBesselK(l, dg);

		return std::exp(logno + beta * d
			+ (0.5 - l)*(std::log(alpha) - 0.5*std::log(thing))
			+ LnBesselK(l - 0.5, alpha*std::sqrt(thing)));// + std::log(std::fabs(beta)+0.0001) );

	}


	double diff_eval(double d, double l, double alpha, double beta, double delta) {
		const double gamma = alpha;
		const double dg = delta * gamma;

		const double thing = delta * delta + d * d;
		const double sqrthing = std::sqrt(thing);
		const double alphasq = alpha * sqrthing;
		const double no = std::pow(gamma / delta, l) / besselK(l, dg)*sq2pi_inv;
		const double ns1 = 0.5 - l;

		return no * std::pow(alpha, ns1) * std::pow(thing, l / 2. - 1.25)
			* (-d * alphasq * (besselK(l - 1.5, alphasq)
				+ besselK(l + 0.5, alphasq))
				+ (2.*(beta*thing + d * l) - d) * besselK(ns1, alphasq))
			* std::exp(beta*d) * 0.5;
	}

	/*
	double Gauss2F1(double a, double b, double c, double x){
	  if (fabs(x) <= 1.) {
		return ROOT::Math::hyperg(a, b, c, x);
	  } else {
		return ROOT::Math::hyperg(c-a, b, c, 1-1/(1-x))/std::pow(1-x, b);
	  }
	}

	double stIntegral(double d1, double delta, double l){
	  return d1 * Gauss2F1(0.5, 0.5-l, 3./2, -d1*d1/(delta*delta));
	}
	*/
}

double RooHypatia2::evaluate() const
{
	const double d = _x - _mu;
	const double cons0 = std::sqrt(_zeta);
	const double asigma = _a * _sigma;
	const double a2sigma = _a2 * _sigma;
	const double beta = _beta;
	double out = 0.;

	if (_zeta > 0.) {
		// careful if zeta -> 0. You can implement a function for the ratio,
		// but careful again that |nu + 1 | != |nu| + 1 so you have to deal with the signs
		const double phi = besselK(_lambda + 1., _zeta) / besselK(_lambda, _zeta);
		const double cons1 = _sigma / std::sqrt(phi);
		const double alpha = cons0 / cons1;
		const double delta = cons0 * cons1;

		if (d < -asigma) {
			const double k1 = LogEval(-asigma, _lambda, alpha, beta, delta);
			const double k2 = diff_eval(-asigma, _lambda, alpha, beta, delta);
			const double B = -asigma + _n * k1 / k2;
			const double A = k1 * std::pow(B + asigma, _n);

			out = A * std::pow(B - d, -_n);
		}
		else if (d > a2sigma) {
			const double k1 = LogEval(a2sigma, _lambda, alpha, beta, delta);
			const double k2 = diff_eval(a2sigma, _lambda, alpha, beta, delta);
			const double B = -a2sigma - _n2 * k1 / k2;
			const double A = k1 * std::pow(B + a2sigma, _n2);

			out = A * std::pow(B + d, -_n2);
		}
		else {
			out = LogEval(d, _lambda, alpha, beta, delta);
		}
	}
	else if (_zeta < 0.) {
		coutE(Eval) << "The parameter " << _zeta.GetName() << " of the RooHypatia2 " << GetName() << " cannot be < 0." << std::endl;
	}
	else if (_lambda < 0.) {
		const double delta = _sigma;

		if (d < -asigma) {
			const double cons1 = std::exp(-beta * asigma);
			const double phi = 1. + _a * _a;
			const double k1 = cons1 * std::pow(phi, _lambda - 0.5);
			const double k2 = beta * k1 - cons1 * (_lambda - 0.5) * std::pow(phi, _lambda - 1.5) * 2.*_a / delta;
			const double B = -asigma + _n * k1 / k2;
			const double A = k1 * std::pow(B + asigma, _n);

			out = A * std::pow(B - d, -_n);
		}
		else if (d > a2sigma) {
			const double cons1 = std::exp(beta*a2sigma);
			const double phi = 1. + _a2 * _a2;
			const double k1 = cons1 * std::pow(phi, _lambda - 0.5);
			const double k2 = beta * k1 + cons1 * (_lambda - 0.5) * std::pow(phi, _lambda - 1.5) * 2.*_a2 / delta;
			const double B = -a2sigma - _n2 * k1 / k2;
			const double A = k1 * std::pow(B + a2sigma, _n2);

			out = A * std::pow(B + d, -_n2);
		}
		else {
			out = std::exp(beta*d) * std::pow(1. + d * d / (delta*delta), _lambda - 0.5);
		}
	}
	else {
		coutE(Eval) << "zeta = 0 only supported for lambda < 0. lambda = " << double(_lambda) << std::endl;
	}

	return out;
}




