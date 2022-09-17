
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



#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../ChiTreeInit.C"
#include "../ChiFitterInit.h"

//#define UsesTMVA // if defined, the code needs to be run in CMSSW_10_3_X, otherwise CMSSW_8_0_X is good enough (production release)

#ifdef UsesTMVA
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif


using namespace std;
using namespace RooFit;

ofstream file_log;

string myPdfName = "DoubleSidedCB_BPHTreshold";


///*
TGraphAsymmErrors* gAsChic_alpha_pT, *gAsChic_alpha_y, *gAsChic_alpha_nTrk;
TGraphAsymmErrors* gAsChic_n_pT, *gAsChic_n_y, *gAsChic_n_nTrk;
TGraphAsymmErrors* gAsChic_sigmaRat_pT, *gAsChic_sigmaRat_y, *gAsChic_sigmaRat_nTrk;
TGraphAsymmErrors* gAsChic_sigma1_pT, *gAsChic_sigma1_y, *gAsChic_sigma1_nTrk;
TGraphAsymmErrors* gAsChic_mean1_pT, *gAsChic_mean1_y, *gAsChic_mean1_nTrk;
//*/

TGraphAsymmErrors* gAsFitOutputArray[nFittingSets][nFitFunctionParams]; // stores all the parameters from fits in one giant pile - used in systematics
TGraphAsymmErrors* gAsFitOutputArrayJpsi[nFittingSets][nFitFunctionParams]; // stores all the parameters from fits in one giant pile - used in systematics


TLorentzVector* LVchic, *LVdimuon, *LVconv, *LVmuon1, *LVmuon2;
TLorentzVector* LVchic_rot, *LVchic_rotGamma, *LVdimuon_rot, *LVconv_rot, *LVmuon1_rot, *LVmuon2_rot;
TLorentzVector LVaux;


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

int GetRatio (TGraphAsymmErrors* gAsResult, TGraphAsymmErrors* gAsNum, TGraphAsymmErrors* gAsDen) // dependent on root version, this is prior 6.20
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
			gAsResult->SetPointError(i, 0,0,0,0);
		}
	}

	return 0;
}

int GetCorrectedRatio(TGraphAsymmErrors* gAsResult, TGraphAsymmErrors* gAsNum, TGraphAsymmErrors* gAsDen, const char* fileCorrection, const char* corrName)
{
	GetRatio(gAsResult, gAsNum, gAsDen);
	TFile* fCor = new TFile(fileCorrection, "READ");
	cout << "Opening correction: "<< corrName << endl;
	TH1D* hCor = (TH1D*)fCor->Get(corrName);
	cout << "Opened" << endl;
	for (int i = 0; i < hCor->GetNbinsX(); i++)
	{
		double corr = hCor->GetBinContent(i + 1); //TH1 numbering offset by 1
		cout << corr << endl;
		gAsResult->GetY()[i] *= (1 / corr);
		gAsResult->SetPointEYhigh(i, gAsResult->GetErrorYhigh(i) / corr);
		gAsResult->SetPointEYlow(i, gAsResult->GetErrorYlow(i) / corr);
	}
	cout << "DoneCorrecting" << endl;
	fCor->Close();
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
		Ws.factory("EXPR::background('pow(delta,alpha_bkg)*exp(b1)',delta, b1, alpha_bkg[0.2, 0.1, 1.2])");  //upper was 0.8

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
	 
		Ws.factory("SUM::nominalFlatBkg(nsig[50,0,10000000]*signal, nbkg[2000,0,10000000]*background)");
		return true;
	}
	else if (pdfName.compare("nominalPdfJpsi") == 0)
		//else if (pdfName.find("nominalPdfJpsi") != string::npos)
	{
		//Ws.factory("CBShape::signalJpsi(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08], alphaJpsi[1.85, 0.1, 50], nJpsi[1.7, 0.2, 50])");
		Ws.factory("CBShape::Jpsi1(rvmassJpsi, mean1Jpsi[3.097, 3.05, 3.15], sigma1Jpsi[0.03, 0.003, 0.08], alphaJpsi[1.85, 0.1, 50], nJpsi[1.7, 0.2, 50])");
		Ws.factory("prod::sigma2Jpsi(sigma1Jpsi, sigmaRatioJpsi[2.0, 1.2, 3.0])"); // allowing sigmas to merge (i.e. 1.0) leads to fit instabilities, and is useless (same as one peak - so the shape can be fitted here too). In any case, is 1.5-1.6 in the fits
		Ws.factory("CBShape::Jpsi2(rvmassJpsi, mean1Jpsi, sigma2Jpsi, alphaJpsi, nJpsi)");
		//Ws.factory("CBShape::psi2(rvmassJpsi, mean2Jpsi[3.686, 3.50, 3.75], sigma3Jpsi[0.03, 0.003, 0.08], alphaJpsi, nJpsi)");
		Ws.factory("SUM::signalJpsi(ratioJpsi[0.5,0.00,1.0]*Jpsi1, Jpsi2)");
		//Ws.factory("SUM::signalJpsi(ratioPsi[0.1,0.00,2.0]*psi2, Jpsi)");

		//Ws.factory("Uniform::one(rvmass)");
		//Ws.factory("SUM::expConst(const[0.8,0.01,1]*one, expbkg)");
		//Ws.factory("EXPR::background('expConst*ErfPdf', expConst, ErfPdf)");
		Ws.factory("Exponential::backgroundJpsi(rvmassJpsi,e_1Jpsi[-0.2,-0.7,0])");
 
		Ws.factory("SUM::nominalPdfJpsi(nsigJpsi[5000,0,2000000]*signalJpsi, nbkgJpsi[2000,0,1000000]*backgroundJpsi)");
		return true;
	}

	cout << " ERROR: MODEL CREATION FAILED, POSSIBLY UNDEFINED MODEL WAS CHOSEN" << endl;
	return false;
}



bool RefreshModel(RooWorkspace& Ws, string pdfName, bool isJpsi) // attempt to prevent the fits in the differential bins to get stuck in a weird state when they stop converging properly (empty bins making errors too small, and those errors used to determine range in the next bin)
{
	cout << "Refreshing model...";
	if (isJpsi == true) {
		file_log << Ws.var("mean1Jpsi")->getError() << endl;
		file_log << Ws.var("sigma1Jpsi")->getError() << endl;

		Ws.var("mean1Jpsi")->removeError();
		Ws.var("sigma1Jpsi")->removeError();
		Ws.var("sigmaRatioJpsi")->removeError();
		Ws.var("alphaJpsi")->removeError();
		Ws.var("nJpsi")->removeError();
		Ws.var("e_1Jpsi")->removeError();
		Ws.var("nsigJpsi")->removeError();
		Ws.var("nbkgJpsi")->removeError();

		Ws.var("mean1Jpsi")->setVal(3.097);
		Ws.var("sigma1Jpsi")->setVal(0.03);
		Ws.var("sigmaRatioJpsi")->setVal(2.0);
		Ws.var("alphaJpsi")->setVal(1.85);
		Ws.var("nJpsi")->setVal(1.7);
		Ws.var("e_1Jpsi")->setVal(-0.4);
		Ws.var("nsigJpsi")->setVal(40000); //50000
		Ws.var("nbkgJpsi")->setVal(3000); //2000

	}
	else
	{
		if (pdfName.compare("DoubleSidedCB_BPHTreshold") == 0)
		{
			Ws.var("c2ratio")->removeError();
			Ws.var("sigma1")->removeError();
			Ws.var("alpha_bkg")->removeError();
			Ws.var("beta_bkg")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			Ws.var("sigma1")->setVal(0.005);
			Ws.var("alpha_bkg")->setVal(0.2);
			Ws.var("beta_bkg")->setVal(-2);
			Ws.var("nsig")->setVal(500);
			Ws.var("nbkg")->setVal(2000);
		}
		else if (pdfName.compare("DoubleSidedCB_D0BG") == 0)
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
		else if (pdfName.compare("Hypatia_BPHTreshold") == 0)
		{
			Ws.var("c2ratio")->removeError();
			Ws.var("sigma1")->removeError();
			Ws.var("alpha_bkg")->removeError();
			Ws.var("beta_bkg")->removeError();
			Ws.var("nsig")->removeError();
			Ws.var("nbkg")->removeError();

			Ws.var("c2ratio")->setVal(0.4);
			Ws.var("sigma1")->setVal(0.005);
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
		if (std::abs(xVal- (bins[bin] + bins[bin + 1]) / 2.0)>0.01) cout <<endl<< "WARNING: x values of binning between fit and constraints don't agree! " << xVal << "   " << (bins[bin] + bins[bin + 1]) / 2.0 << endl<<endl;
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



int FitRooDataSet(TGraphAsymmErrors* gAsResult, double* bins, int nbins, RooRealVar* rvmass, RooWorkspace& Ws, bool isJpsi = false, string binVarName = "", string myPdfName = "", string extraCut = "", string fittingSetName = "", bool bConstrainedFit = false, const char* fileConstraints = "") { //uses global variables
	cout << endl << endl << endl << "NEW SET, IN FITTING " << fittingSetName << "   " << nbins << endl << endl;
	TCanvas *cFit = new TCanvas("cFit", "cFit", 1000, 800);
	//gPad->SetLeftMargin(0.15);
	cFit->cd();
	//cFit->Divide(1, 2);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.275, 0.98, 1.0);
	pad1->SetTicks(1, 1);
	pad1->Draw();

	RooAbsReal::defaultIntegratorConfig()->setEpsAbs(2e-8);
	RooAbsReal::defaultIntegratorConfig()->setEpsRel(2e-8);

	//pull pad
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.006, 0.98, 0.360);
	pad2->SetTopMargin(0); // Upper and lower plot are joined
	pad2->SetBottomMargin(0.67);
	pad2->SetTicks(1, 1);

	double mass_windowFitCommon_l, mass_windowFitCommon_h;
	int nMassBinsCommon;

	if (isJpsi == false) {
		mass_windowFitCommon_l = mass_windowFit_l;
		mass_windowFitCommon_h = mass_windowFit_h;
		nMassBinsCommon = nMassBins;
	}
	else {
		mass_windowFitCommon_l = mass_windowFitJpsi_l;
		mass_windowFitCommon_h = mass_windowFitJpsi_h;
		nMassBinsCommon = nMassBinsJpsi;
	}

	int nBinSet = -1;
	for (int k = 0; k < nFittingSets; k++) { //find the corresponding set index for readout
		if (fittingSets[k] == fittingSetName)
		{
			nBinSet = k;
			break;
		}
	}
	if (nBinSet == -1) { cout << "ERROR: This sets of bins wasn't found, check fittingSetName and fittingSets definition in ChiFitterInit.h     STOPPING HERE, THIS SET WILL BE SKIPPED " << endl; return -1; }



	for (int i = 0; i < nbins; i++) {
		cout << endl << "Bins " << bins[i] << "  to  " << bins[i + 1] << endl;
		pad1->cd();
		RooPlot* massframeBin;

		//if (isJpsi == false) {
		massframeBin = rvmass->frame(mass_windowFitCommon_l, mass_windowFitCommon_h, nMassBinsCommon); 
		//}
		//else { massframeBin = rvmass->frame(mass_windowFitJpsi_l, mass_windowFitJpsi_h, nMassBinsJpsi); }
		massframeBin->SetTitle("mass");
		TString TstrCut = binVarName + TString::Format(" > %f", bins[i]) + " && " + binVarName + TString::Format(" < %f", bins[i + 1]) + extraCut;
		cout << TstrCut << endl;
		string strCut = TstrCut.Data();
		RooDataSet* rdsDataBin;
		if (isJpsi == false) { rdsDataBin = (RooDataSet*)Ws.data("rdsNominal")->reduce(strCut.c_str()); 
			//rdsDataBin = (RooDataSet*)rdsDataBin->reduce(mass_windowFit.c_str());
			RefreshModel(Ws, myPdfName, false);
		}else {
			rdsDataBin = (RooDataSet*)Ws.data("rdsNominalJpsi")->reduce(strCut.c_str()); 
			//rdsDataBin = (RooDataSet*)rdsDataBin->reduce(mass_windowFitJpsi.c_str());
			RefreshModel(Ws, myPdfName, true);
		}

		rdsDataBin->plotOn(massframeBin, Name("dataHist"));

		if (bConstrainedFit == true) {
			int flagConstrSet = SetConstraints(Ws, bins, i, fittingSetName, myPdfName, fileConstraints); //sets constraints, returns -1 if a bad set is passed (can still fail in other ways)
			if (flagConstrSet == -1) break; //if it doesn't work, skip all the set
		}


		RooFitResult* fitResultBin = Ws.pdf(myPdfName.c_str())->fitTo(*rdsDataBin, Extended(true), /*Range(mass_windowFitCommon_l, mass_windowFitCommon_h),*/ SumW2Error(true), NumCPU(1), PrintLevel(-1), Save(true));
		//fitResultBin->Print("v");
		//Ws.import(*fitResultBin, TString::Format("fitResult_%s", myPdfNameBin.c_str()));
		if (isJpsi == false) {
			cout << "signal " << Ws.var("nsig")->getValV() << endl;
			SetPointFromFit(gAsResult, "nsig", Ws, i, bins);

			// And also store all the other parameters in the output file
			RooArgList paramList = fitResultBin->floatParsFinal();
			RooRealVar *par;
			for (int j = 0; j < paramList.getSize(); j++)
			{
				par = (RooRealVar*)paramList.at(j);

				// if the first bin, create the graph
				if (i == 0) {
					gAsFitOutputArray[nBinSet][j] = new TGraphAsymmErrors(nbins);
					cout << "gAsChic_FitOutput_" + par->getTitle() + "_" + fittingSets.at(nBinSet) << endl;
					cout << "Chic " + fittingSets.at(nBinSet) + " " + par->getTitle() + " dependence" << endl;
					gAsFitOutputArray[nBinSet][j]->SetNameTitle("gAsChic_FitOutput_" + par->getTitle() + "_" + fittingSets.at(nBinSet), "Chic " + fittingSets.at(nBinSet) + " " + par->getTitle() + " dependence");
				}

				// save as tgraphAsymmErrors
				SetPointFromFit(gAsFitOutputArray[nBinSet][j], (string)par->getTitle(), Ws, i, bins);
			}

		}
		else {
			cout << "signal Jpsi " << Ws.var("nsigJpsi")->getValV() << endl;
			file_log << i << " signal Jpsi " << Ws.var("nsigJpsi")->getValV() << endl;
			SetPointFromFit(gAsResult, "nsigJpsi", Ws, i, bins);

			// And also store all the other parameters in the output file
			RooArgList paramList = fitResultBin->floatParsFinal();
			RooRealVar *par;
			for (int j = 0; j < paramList.getSize(); j++)
			{
				par = (RooRealVar*)paramList.at(j);

				// if the first bin, create the graph
				if (i == 0) {
					gAsFitOutputArrayJpsi[nBinSet][j] = new TGraphAsymmErrors(nbins);
					cout << "gAsJpsi_FitOutput_" + par->getTitle() + "_" + fittingSets.at(nBinSet) << endl;
					cout << "Jpsi " + fittingSets.at(nBinSet) + " " + par->getTitle() + " dependence" << endl;
					gAsFitOutputArrayJpsi[nBinSet][j]->SetNameTitle("gAsJpsi_FitOutput_" + par->getTitle() + "_" + fittingSets.at(nBinSet), "Chic " + fittingSets.at(nBinSet) + " " + par->getTitle() + " dependence");
				}

				// save as tgraphAsymmErrors
				SetPointFromFit(gAsFitOutputArrayJpsi[nBinSet][j], (string)par->getTitle(), Ws, i, bins);
			}


		}

				
		if (isJpsi == false) {
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("background"), LineStyle(kDashed));
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic1"), LineStyle(kDashed), LineColor(kRed));
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("chic2"), LineStyle(kDashed), LineColor(kGreen));
		}
		else {
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("backgroundJpsi"), LineStyle(kDashed));
			Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("signalJpsi"), LineStyle(kDashed), LineColor(kRed));
			//Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Components("psi2"), LineStyle(kDashed), LineColor(kGreen));

		}

		Ws.pdf(myPdfName.c_str())->plotOn(massframeBin, Name("FullPdf"));
		if (isJpsi == false) { Ws.pdf(myPdfName.c_str())->paramOn(massframeBin, Layout(0.65, 0.95)); }
		else{ Ws.pdf(myPdfName.c_str())->paramOn(massframeBin, Layout(0.20, 0.55)); }
		massframeBin->GetXaxis()->SetTitleSize(0);
		massframeBin->GetXaxis()->SetLabelSize(0);
		massframeBin->Draw();

		//chi2

		int numFitPar = fitResultBin->floatParsFinal().getSize();
		//cout << "Npar: " << numFitPar << endl;
		double chi2_fit = massframeBin->chiSquare("FullPdf", "dataHist", numFitPar);
		TPaveText *textchi = new TPaveText(.23, .75, .38, .83, "brNDC");
		if (isJpsi == true) { 
			textchi->SetX1(.68); //move the text to the right
			textchi->SetX2(.88);	
		}
		textchi->SetFillColor(kWhite);
		textchi->SetTextColor(kBlue);
		textchi->SetTextSize(0.05);
		textchi->SetBorderSize(0);
		string txtchi = "#chi^{2}/ndf: " + to_string(chi2_fit);
		textchi->AddText(txtchi.c_str());
		textchi->Draw();


		//PULLS

		pad2->cd();
		//RooHist* hpull = massframeBin->pullHist("rPullHist", myPdfName.c_str());
		RooHist* hpull = massframeBin->pullHist("dataHist", "FullPdf");
		hpull->SetMarkerSize(0.8);
		RooPlot* pullFrame = rvmass->frame(Title("Pull Distribution"));
		pullFrame->addPlotable(hpull, "P");
		pullFrame->SetTitleSize(0);
		pullFrame->GetYaxis()->SetTitle("Pull");
		pullFrame->GetYaxis()->SetTitleSize(0.19);
		pullFrame->GetYaxis()->SetLabelSize(0.14);
		//pullFrame->GetYaxis()->SetLabelOffset(1.1);
		pullFrame->GetYaxis()->SetRangeUser(-5.0, 5.0);
		pullFrame->GetYaxis()->SetNdivisions(502, kTRUE);
		pullFrame->GetYaxis()->CenterTitle();
		pullFrame->GetXaxis()->SetLabelSize(0.20);
		if (isJpsi == false) { pullFrame->GetXaxis()->SetTitle("M (#mu#mu#gamma - #mu#mu + 3.097) [GeV/c^{2}]"); }
		else { pullFrame->GetXaxis()->SetTitle("M(#mu#mu)[GeV/c^{2}]"); }
		pullFrame->GetXaxis()->SetTitleOffset(1.05);
		pullFrame->GetXaxis()->SetTitleSize(0.20);
		pullFrame->Draw();

		TLine *l1 = new TLine(mass_windowFitCommon_l, 0, mass_windowFitCommon_h, 0);
		l1->SetLineStyle(9);
		l1->Draw("same");
		//pad1->Update();

		cFit->cd();
		pad1->Draw();
		pad2->Draw();

		pad1->Update();
		pad2->Update();

		cFit->SaveAs(((string)"FitterOutput/FitResult_" + (isJpsi?"Jpsi_":"Chic_") + fittingSetName + "_" + binVarName + Form("_%ibin_", i) + Form("_%.1f_%.1f.png", bins[i], bins[i + 1])).c_str());
		if (isJpsi == true) cFit->SaveAs(((string)"FitterOutput/FitResult_" + (isJpsi ? "Jpsi_" : "Chic_") + fittingSetName + "_" + binVarName + Form("_%ibin_", i) + Form("_%.1f_%.1f.root", bins[i], bins[i + 1])).c_str());
	}

	delete cFit;


	return 0;
}



///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////

void Analyze_Chic(bool flagGenerateRds = false, bool flagRunFits = true, const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV-Nominal_v2-bothDir.root", const char* fileOut = "Chi_c_output_Nominal_v2-bothDir_DCB_test.root", const char* fileRds = "rds_Nominal_v2_test.root", bool isMC = false, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints_Officialv3_NarrowRangeDCBW.root", const char* fileCorrection = "Chi_c_WeightsMC_Official_v3-bothDir.root")
//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true, const char* fileIn = "/afs/cern.ch/user/o/okukral/Chic_pPb/CMSSW_8_0_30/src/HeavyIonsAnalysis/ChiAnalysis/test/Chi_c_pPb8TeV-Comp285993.root", const char* fileOut = "Chi_c_output_RW6_ComparisonAlberto285993.root", const char* fileRds = "rds_save_test.root", bool isMC = false, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = false,  const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_AOD_Pbp_RW6Comp5/PADoubleMuon/crab_Chi_c_pPb8TeV_AOD_Pbp_RW6Comp5/211202_000528/0000/Chi_c_pPb8TeV-PbpCompTest.root", const char* fileOut = "Chi_c_output_RW4_testCutTightRefitAlberto2.root", const char* fileRds = "rds_RW4_Full_noWeights_CutTightRefitAlberto.root", bool isMC = false, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true,  const char* fileIn = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-MC8_BothDir.root", const char* fileOut = "Chi_c_output_MC8_test.root", const char* fileRds = "rds_MC8_test.root", bool isMC = true, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraints.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true,  const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC9-pPb.root", const char* fileOut = "Chi_c_output_MC9_test.root", const char* fileRds = "rds_MC9_test.root", bool isMC = true, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraintsMC8.root", const char* fileCorrection = "Chi_c_WeightsMC8_pPb_comparisonBothDir.root")
//void Analyze_Chic(bool flagGenerateRds = true, bool flagRunFits = true, const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC9-pPb.root", const char* fileOut = "Chi_c_output_UsingMC9_test.root", const char* fileRds = "rds_UsingMC9_test.root", bool isMC = false, bool flagConstrainedFit = true, const char* fileConstraints = "Chi_c_constraintsMC8.root", const char* fileCorrection = "Chi_c_WeightsMC9_bothDir.root")
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
	TGraphAsymmErrors* gAsChic_pT_mid = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_mid->SetNameTitle("gAsChic_pT_mid", "Chic pT dependence - midrapidity");
	TGraphAsymmErrors* gAsChic_pT_fwd = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_fwd->SetNameTitle("gAsChic_pT_fwd", "Chic pT dependence - forward + backward rapidity");
	TGraphAsymmErrors* gAsChic_pT_fwdOnly = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_fwdOnly->SetNameTitle("gAsChic_pT_fwdOnly", "Chic pT dependence - forward rapidity");
	TGraphAsymmErrors* gAsChic_pT_bkwOnly = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_bkwOnly->SetNameTitle("gAsChic_pT_bkwOnly", "Chic pT dependence - backward rapidity");
	TGraphAsymmErrors* gAsChic_pT_fwdOnlyWide = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_fwdOnlyWide->SetNameTitle("gAsChic_pT_fwdOnlyWide", "Chic pT dependence - forward rapidity");
	TGraphAsymmErrors* gAsChic_pT_bkwOnlyWide = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_bkwOnlyWide->SetNameTitle("gAsChic_pT_bkwOnlyWide", "Chic pT dependence - backward rapidity");

	TGraphAsymmErrors* gAsChic_pT_midCMS = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_midCMS->SetNameTitle("gAsChic_pT_midCMS", "Chic pT dependence - midrapidity CMS");
	TGraphAsymmErrors* gAsChic_pT_fwdCMS = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_fwdCMS->SetNameTitle("gAsChic_pT_fwdCMS", "Chic pT dependence - forward rapidity CMS");
	TGraphAsymmErrors* gAsChic_pT_bkwCMS = new TGraphAsymmErrors(nbins_pT);
	gAsChic_pT_bkwCMS->SetNameTitle("gAsChic_pT_bkwCMS", "Chic pT dependence - backward rapidity CMS");


	TGraphAsymmErrors* gAsJpsi_pT = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT->SetNameTitle("gAsJpsi_pT", "Jpsi pT dependence");
	TGraphAsymmErrors* gAsJpsi_pT_mid = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_mid->SetNameTitle("gAsJpsi_pT_mid", "Jpsi pT dependence - midrapidity");
	TGraphAsymmErrors* gAsJpsi_pT_fwd = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_fwd->SetNameTitle("gAsJpsi_pT_fwd", "Jpsi pT dependence - forward + backward rapidity");
	TGraphAsymmErrors* gAsJpsi_pT_fwdOnly = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_fwdOnly->SetNameTitle("gAsJpsi_pT_fwdOnly", "Jpsi pT dependence - forward rapidity");
	TGraphAsymmErrors* gAsJpsi_pT_bkwOnly = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_bkwOnly->SetNameTitle("gAsJpsi_pT_bkwOnly", "Jpsi pT dependence - backward rapidity");
	TGraphAsymmErrors* gAsJpsi_pT_fwdOnlyWide = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_fwdOnlyWide->SetNameTitle("gAsJpsi_pT_fwdOnlyWide", "Jpsi pT dependence - forward rapidity");
	TGraphAsymmErrors* gAsJpsi_pT_bkwOnlyWide = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_bkwOnlyWide->SetNameTitle("gAsJpsi_pT_bkwOnlyWide", "Jpsi pT dependence - backward rapidity");

	TGraphAsymmErrors* gAsJpsi_pT_midCMS = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_midCMS->SetNameTitle("gAsJpsi_pT_midCMS", "Jpsi pT dependence - midrapidity CMS");
	TGraphAsymmErrors* gAsJpsi_pT_fwdCMS = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_fwdCMS->SetNameTitle("gAsJpsi_pT_fwdCMS", "Jpsi pT dependence - forward rapidity CMS");
	TGraphAsymmErrors* gAsJpsi_pT_bkwCMS = new TGraphAsymmErrors(nbins_pT);
	gAsJpsi_pT_bkwCMS->SetNameTitle("gAsJpsi_pT_bkwCMS", "Jpsi pT dependence - backward rapidity CMS");


	TGraphAsymmErrors* gAsRatio_pT = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT->SetNameTitle("gAsRatio_pT", "Ratio pT dependence");
	TGraphAsymmErrors* gAsRatio_pT_mid = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_mid->SetNameTitle("gAsRatio_pT_mid", "Ratio pT dependence - midrapidity");
	TGraphAsymmErrors* gAsRatio_pT_fwd = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_fwd->SetNameTitle("gAsRatio_pT_fwd", "Ratio pT dependence - forward + backward rapidity");
	TGraphAsymmErrors* gAsRatio_pT_fwdOnly = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_fwdOnly->SetNameTitle("gAsRatio_pT_fwdOnly", "Ratio pT dependence - forward rapidity");
	TGraphAsymmErrors* gAsRatio_pT_bkwOnly = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_bkwOnly->SetNameTitle("gAsRatio_pT_bkwOnly", "Ratio pT dependence - backward rapidity");
	TGraphAsymmErrors* gAsRatio_pT_fwdOnlyWide = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_fwdOnlyWide->SetNameTitle("gAsRatio_pT_fwdOnlyWide", "Ratio pT dependence - forward rapidity");
	TGraphAsymmErrors* gAsRatio_pT_bkwOnlyWide = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_bkwOnlyWide->SetNameTitle("gAsRatio_pT_bkwOnlyWide", "Ratio pT dependence - backward rapidity");

	TGraphAsymmErrors* gAsRatio_pT_midCMS = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_midCMS->SetNameTitle("gAsRatio_pT_midCMS", "Ratio pT dependence - midrapidity CMS");
	TGraphAsymmErrors* gAsRatio_pT_fwdCMS = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_fwdCMS->SetNameTitle("gAsRatio_pT_fwdCMS", "Ratio pT dependence - forward rapidity CMS");
	TGraphAsymmErrors* gAsRatio_pT_bkwCMS = new TGraphAsymmErrors(nbins_pT);
	gAsRatio_pT_bkwCMS->SetNameTitle("gAsRatio_pT_bkwCMS", "Ratio pT dependence - backward rapidity CMS");

	TGraphAsymmErrors* gAsChic_y = new TGraphAsymmErrors(nbins_y);
	gAsChic_y->SetNameTitle("gAsChic_y", "Chic y dependence");
	TGraphAsymmErrors* gAsJpsi_y = new TGraphAsymmErrors(nbins_y);
	gAsJpsi_y->SetNameTitle("gAsJpsi_y", "Jpsi y dependence");
	TGraphAsymmErrors* gAsRatio_y = new TGraphAsymmErrors(nbins_y);
	gAsRatio_y->SetNameTitle("gAsRatio_y", "Ratio y dependence");

	TGraphAsymmErrors* gAsChic_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_nTrk->SetNameTitle("gAsChic_nTrk", "Chic nTrk midrapidity dependence");
	TGraphAsymmErrors* gAsJpsi_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsJpsi_nTrk->SetNameTitle("gAsJpsi_nTrk", "Jpsi nTrk midrapidity dependence");
	TGraphAsymmErrors* gAsRatio_nTrk = new TGraphAsymmErrors(nbins_nTrk);
	gAsRatio_nTrk->SetNameTitle("gAsRatio_nTrk", "Ratio nTrk midrapidity dependence");

	TGraphAsymmErrors* gAsChic_nTrk_all = new TGraphAsymmErrors(nbins_nTrk);
	gAsChic_nTrk_all->SetNameTitle("gAsChic_nTrk_all", "Chic nTrk dependence");
	TGraphAsymmErrors* gAsJpsi_nTrk_all = new TGraphAsymmErrors(nbins_nTrk);
	gAsJpsi_nTrk_all->SetNameTitle("gAsJpsi_nTrk_all", "Jpsi nTrk dependence");
	TGraphAsymmErrors* gAsRatio_nTrk_all = new TGraphAsymmErrors(nbins_nTrk);
	gAsRatio_nTrk_all->SetNameTitle("gAsRatio_nTrk_all", "Ratio nTrk dependence");

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


		//file_log << "index runNumber  eventNumber  (chiCandPerEvent-ignore)  rap_chi  pT_chi  Mass (refit) muon1_eta muon1_pt  muon2_eta muon2_pt muon1_pass muon2_pass (conv position) conv_eta conv_pt vProb Conversion cuts: convHP  convGT  conv_rho  convSig1  convSig2  conv_IHits  conv_vProb  convDOF1  convDOF2  convChi1  convChi2  convDOA ctau/ctauErr ctau3D ctauErr significance(Jpsi) ctauJpsi ctauErr vertexNumber Conv_duplicity_status" << endl;
		
		//TMVA 
		#ifdef UsesTMVA


		TMVA::Tools::Instance();
		TMVA::Reader *TMWAreader = new TMVA::Reader("!Color:!Silent");

		float convQuality_isHighPurityValue, convQuality_isGeneralTracksOnlyValue, conv_tkVtxCompatibilityOKValue, conv_compatibleInnerHitsOKValue, conv_tk1NumOfDOFValue, conv_tk2NumOfDOFValue; //all variables need to be floats
		float conv_vertexPositionRhoValue, conv_sigmaTkVtx1Value, conv_sigmaTkVtx2Value, conv_vertexChi2ProbValue, conv_track1Chi2Value, conv_track2Chi2Value;
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
		TMWAreader->AddVariable("conv_track1Chi2", &conv_track1Chi2Value);
		TMWAreader->AddVariable("conv_track2Chi2", &conv_track2Chi2Value);
		TMWAreader->AddVariable("conv_minDistanceOfApproach", &conv_minDistanceOfApproachValue);
		TMWAreader->AddSpectator("conv_eta", &conv_etaValue);
		TMWAreader->AddSpectator("conv_pt", &conv_ptValue);

		TMWAreader->BookMVA("BDT", "BDTWeight/TMVAClassification_BDT2new.weights.xml");
		#endif


		//load corrections
		TFile* fCor = new TFile(fileCorrection, "READ");

		TH2D* hWeightChic = (TH2D*)fCor->Get("h_chiAccEff_rat");
		TH2D* hWeightJpsi = (TH2D*)fCor->Get("h_JpsiAccEff_rat");

		int weird_decay_counter = 0;
		long nchicCounter = 0, nchicCounterPass = 0, nchicCounterPassMass = 0, muon1ptCounter = 0, muon2ptCounter = 0;
		long nConvCounter = 0, nConvSplitCounter = 0, nConvCounterPeak = 0, nConvSplitCounterPeak = 0;
		int testCounter = 0;
		bool passDimSel = false;
		bool passDimSelTight = false;
		Long64_t nentries = event_tree->GetEntries();
		//if (nentries > 5000) { nentries = 5000; }
		cout << nentries << endl;
		for (Long64_t i = 0; i < nentries; i++)
		{

			event_tree->GetEntry(i);
			if (i % 50000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }
			if (isMC == true) {
				if (gen_pdgId->at(0) != PythCode_chic1 && gen_pdgId->at(0) != PythCode_chic2) { continue; } //remove chic0 from MC, expecting one gen chic per event
				if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
				{
					weird_decay_counter++;
					continue;
				}
			}

			//if (isMC == true && (gen_pdgId->at(0) != PythCode_chic1)) { continue; } // only chic1
			//if (isMC == true && (gen_pdgId->at(0) != PythCode_chic2)) { continue; } // only chic2

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
				/*
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
				conv_track1Chi2Value = conv_track1Chi2->at(convPos);
				conv_track2Chi2Value = conv_track1Chi2->at(convPos);
				conv_minDistanceOfApproachValue = conv_minDistanceOfApproach->at(convPos);
				conv_etaValue = conv_eta->at(convPos);
				conv_ptValue = conv_pt->at(convPos);

				
				double photMVA = TMWAreader->EvaluateMVA("BDT");
				hphoton_MVA_response->Fill(photMVA);
				//cout << "Values MVA: HP: " << convQuality_isHighPurityValue << "   TracksOnly: " << convQuality_isGeneralTracksOnlyValue << "    Vprob: " << conv_vertexChi2ProbValue << "    and the response: "<< photMVA << endl;
				if (photMVA < -0.14) continue;
				#endif


				if (PhotSelectionPass(convPos) == false) continue;

				if (passDimSel == true) { ++nchicCounterPass; } // SelectionsPassed
				
				//*/												
				

				if (ChiPassAllCuts(iChi) == false) { continue; }
				++nchicCounterPass;
				int dimuonPos = chi_daughterJpsi_position->at(iChi);
				int convPos = chi_daughterConv_position->at(iChi);

				// Get Lorentz V
				LVchic = (TLorentzVector*)chi_p4->At(iChi);
				LVdimuon = (TLorentzVector*)dimuon_p4->At(dimuonPos);
				LVconv = (TLorentzVector*)conv_p4->At(convPos);
				double pT_Jpsi = dimuon_pt->at(dimuonPos);
				double rap_Jpsi = LVdimuon->Rapidity();
				if (runNumber < 285833) rap_Jpsi = -rap_Jpsi; //flip Pbp direction

				double pT_chi = chi_pt->at(iChi);
				double rap_chi = LVchic->Rapidity();
				if (runNumber < 285833) rap_chi = -rap_chi; //flip Pbp direction
				double m_chi = LVchic->M();

				double dimuonM = LVdimuon->M();
				double Mdiff = m_chi - dimuonM + 3.097;// Assume J/psi mass 


				if (rap_Jpsi > -0.566 && rap_Jpsi < 1.434)
				{
					if (pT_Jpsi > 18 && pT_Jpsi < 30)
					{
						testCounter++;
						cout << "counter is " << testCounter << " run number and event: " <<  runNumber << "   "  << eventNumber << " is pPb and rapidity: " << ispPb << " " << rap_Jpsi <<  endl;
					}
				}




				////////////////////////
				////   refit     // comment out if not done
				///////////////////////


				double refit_vProb=0, ctauPV=0, ctauPVError=0, ctauSig=0, ctauPV3D=0;

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

				//////////////
				/// Nominal
				//////////////

				if (dimuonM > mass_cutoffJpsi_l && dimuonM<mass_cutoffJpsi_h && Mdiff > mass_window_l && Mdiff < mass_window_h)
				{
					++nchicCounterPassMass;

				//roofit:
					rvmass->setVal(Mdiff);
					rvpt->setVal(pT_Jpsi);
					rvrap->setVal(rap_Jpsi);
					rvntrack->setVal(pvtx_nTracks->at(dimuon_pvtx_indexFromOniaMuMu->at(dimuonPos))); //for now using total
					//rvweight->setVal(1); // no weights for now
					double accEff_chi = 1;// hWeightChic->GetBinContent(hWeightChic->FindBin(abs(rap_chi), pT_chi));
					//cout << "Test rap: " << rap_chi << "  pt  " << pT_chi << " and value accEff " << accEff_chi << endl;
					if (accEff_chi > 0.00001) {// binning should ensure enough statistics, but if we have empty bin, just set weight to be 1
						//rvweight->setVal(1 / accEff_chi);
						rvweight->setVal(1); //switch off weights
					}
					else rvweight->setVal(1);
					if (isMC) {
						rvweight->setVal(1); //no weights if MC
					}

					rdsNominal->add(*cols, rvweight->getVal());

				}
				////////
							   



//				/*
//
//				//if (passDimSel == true && dimuonM>2.95 && dimuonM<3.2 && Mdiff > mass_window_l && Mdiff < mass_window_h) ++nchicCounterPassMass;
//				if (passDimSel == true && dimuonM > mass_cutoffJpsi_l && dimuonM<mass_cutoffJpsi_h && Mdiff > mass_window_l && Mdiff < mass_window_h) ++nchicCounterPassMass;
//
//
//				//////// Rotational bkg   /////
//				LVmuon1_rot = (TLorentzVector*)muon_p4->At(muon1Pos);
//				LVmuon2_rot = (TLorentzVector*)muon_p4->At(muon2Pos);
//				LVconv_rot = new TLorentzVector(*(TLorentzVector*)conv_p4->At(convPos));
//				LVconv_rot->RotateZ(TMath::Pi());
//				//if (muon_pt->at(muon1Pos) > muon_pt->at(muon2Pos)) { muon1ptCounter++; } //crosscheck that they have ordering
//				//else muon2ptCounter++;
//				if (rand->Uniform(0, 1) < 0.5) { LVmuon1_rot->RotateZ(TMath::Pi()); }
//				else { LVmuon2_rot->RotateZ(TMath::Pi()); }
//				LVdimuon_rot = new TLorentzVector(*LVmuon1_rot + *LVmuon2_rot);
//				//cout << LVdimuon_rot->Pt() << endl;
//				LVchic_rot = new TLorentzVector(*LVdimuon_rot + *LVconv);
//				LVchic_rotGamma = new TLorentzVector(*LVdimuon + *LVconv_rot);
//
//
//
//
//				// Obtain yields
//
//				// fill dimuon stuff
//
//				
//				if (passDimSel == true) {
//					hdimuon_M->Fill(dimuonM);
//					hdimuon_M_rot->Fill(LVdimuon_rot->M());
//				}
//				if (dimuon_charge->at(dimuonPos) != 0) {
//					hdimuon_M_SS->Fill(dimuonM);
//				}
//
//
//				// Rotational bkg
//				if (LVdimuon_rot->M() > 2.95 && LVdimuon_rot->M() < 3.2) {
//					if (passDimSel == true) {
//						hchic_M_rot->Fill(LVchic_rot->M());
//						hSignal_rot->Fill(LVchic_rot->M() - LVdimuon_rot->M() + 3.097);
//					}
//				}
//				// Side-band
//				if (((dimuonM > 2.5 && dimuonM < 2.75) || (dimuonM > 3.4 && dimuonM < 3.7)) && passDimSel == true) {
//					hchic_M_SB->Fill(m_chi);
//					hSignal_SB->Fill(m_chi - dimuonM + 3.097);
//				}
//
//
//				if (dimuonM<mass_cutoffJpsi_l || dimuonM>mass_cutoffJpsi_h) continue; //require narrow dimuon mass
//
//				if (passDimSel == true) {
//					hchic_M->Fill(m_chi); // just raw M
//					hchic_M_rotGamma->Fill(LVchic_rotGamma->M());
//				}
//				if (dimuon_charge->at(dimuonPos) != 0) {
//					hchic_M_SS->Fill(m_chi);
//				}
//				*/
//
//				// nominal
//
//				if (passDimSel == true) {
//					hSignal->Fill(Mdiff);
//					//hSignal_rotGamma->Fill(LVchic_rotGamma->M() - LVdimuon->M() + 3.097);
//					//cout << nchicCounterPass << endl;
//
//					//roofit:
//					if (Mdiff > mass_window_l && Mdiff < mass_window_h) {
//						rvmass->setVal(Mdiff);
//						rvpt->setVal(pT_chi);
//						rvrap->setVal(rap_chi);
//						rvntrack->setVal(ntracks_inEvent); //for now using total
//						//rvweight->setVal(1); // no weights for now
//						double accEff_chi = hWeightChic->GetBinContent(hWeightChic->FindBin(abs(rap_chi), pT_chi));
//						//cout << "Test rap: " << rap_chi << "  pt  " << pT_chi << " and value accEff " << accEff_chi << endl;
//						if (accEff_chi > 0.00001) {// binning should ensure enough statistics, but if we have empty bin, just set weight to be 1
//							//rvweight->setVal(1 / accEff_chi);
//							rvweight->setVal(1); //switch off weights
//						}
//						else rvweight->setVal(1); 
//						if (isMC) {
//							rvweight->setVal(1); //no weights if MC
//						}
//
//						rdsNominal->add(*cols, rvweight->getVal());
//
//						///////////////////
//						// crosschecks
//						/////////////////
//						
//						if (muonIsHLTDoubleMuOpen->at(muon1Pos) != 1 /*|| muonIsHLTDoubleMuOpen->at(muon2Pos) != 1*/) {
//							//cout << "This is the event I want: " << runNumber << "  " << eventNumber << endl; 
//							//file_log << "This is the event I want: " << runNumber << "  " << eventNumber << endl;
//						}
//
//						TVector3* TVconv1 = (TVector3*)conv_vtx->At(convPos);
//						//file_log << nchicCounterPassMass << " " << runNumber << " " << eventNumber << " " << chiCandPerEvent << " " << rap_chi << " " << pT_chi << " " << Mdiff << " m1 " << muon_eta->at(muon1Pos) << " " << muon_pt->at(muon1Pos)<< " m2 " << muon_eta->at(muon2Pos)<< " " << muon_pt->at(muon2Pos) 
//							//<< " Muon Pass:" << MuonSelectionPass(muon1Pos) << MuonSelectionPass(muon2Pos) << " Trig " << muonIsHLTDoubleMuOpen->at(muon1Pos) << muonIsHLTDoubleMuOpen->at(muon2Pos) << " conv " << convPos << " " << conv_eta->at(convPos) << " " << conv_pt->at(convPos) << " Refit_v prob: " << refit_vProb
//							//<< " Conversion cut values: "<< convHP <<  convGT <<  conv_rho <<  convSig1 <<  convSig2 <<  conv_IHits <<  conv_vProb <<  convDOF1 <<  convDOF2 <<  convChi1 <<  convChi2 <<  convDOA
//							//<< " ctau/ctauErr: " << ctauSig<< " "<< ctauPV3D << " " << ctauPVError << " " << (dimuon_ctpv->at(dimuonPos)/ dimuon_ctpvError->at(dimuonPos)) << " " << dimuon_ctpv->at(dimuonPos) << " " << dimuon_ctpvError->at(dimuonPos) << " PVindex: " << dimuon_pvtx_index->at(dimuonPos) << " " << dimuon_pvtx_indexFromOniaMuMu->at(dimuonPos)<< " ConvDuplicity: " << conv_duplicityStatus->at(convPos) << endl;
///*
//						int nRefitNumber2 = -1;
//
//						// check the conversions - if they are many identical ones
//						for (int iChi2 = 0; iChi2 < iChi; iChi2++)
//						{
//							int dimuonPos2 = chi_daughterJpsi_position->at(iChi2);
//							double passDimSel2 = DimuonSelectionPass(dimuonPos2);
//							double passDimSelTight2 = DimuonSelectionPassTight(dimuonPos2, ((TLorentzVector*)dimuon_p4->At(dimuonPos2))->Rapidity());
//							// muon cuts
//							int muon1Pos2 = dimuon_muon1_position->at(dimuonPos2);
//							int muon2Pos2 = dimuon_muon2_position->at(dimuonPos2);
//							if (MuonAcceptanceTight(muon_eta->at(muon1Pos2), muon_pt->at(muon1Pos2)) == false) continue;
//							if (MuonAcceptanceTight(muon_eta->at(muon2Pos2), muon_pt->at(muon2Pos2)) == false) continue;
//							if (MuonSelectionPass(muon1Pos2) == false) continue;
//							if (MuonSelectionPass(muon2Pos2) == false) continue;
//							// photon
//							int convPos2 = chi_daughterConv_position->at(iChi2);
//							if (PhotAcceptance(conv_eta->at(convPos2), conv_pt->at(convPos2)) == false) continue;
//							if (PhotSelectionPassTight(convPos2) == false) continue;
//							passDimSel2 = passDimSelTight2;
//
//							TLorentzVector* LVchic2 = (TLorentzVector*)chi_p4->At(iChi2);
//							TLorentzVector* LVdimuon2 = (TLorentzVector*)dimuon_p4->At(dimuonPos2);
//
//							double pT_chi2 = chi_pt->at(iChi2);
//							double rap_chi2 = LVchic2->Rapidity();
//							double m_chi2 = LVchic2->M();
//							double refit_vProb2 = 0;
//							/////////
//							//refit
//							if (chi_kinematicRefitFlag->at(iChi2) == 1 || chi_kinematicRefitFlag->at(iChi2) == 3) { //good refits
//								nRefitNumber2++; //needed for RW4 only, fix to non-ideal storing of info
//								nRefitNumber2 = iChi2; //overwrite previous line
//
//								//cout << chi_refit_vprob->at(nRefitNumber) << endl;
//								//cout << chi_refit_ctauPV->at(nRefitNumber) << endl;
//								//cout << "RefitStored: " << chi_refitStored->at(nRefitNumber).mass() << endl;
//								m_chi2 = chi_refitStored->at(nRefitNumber2).mass();//use refit mass
//								//refit_vProb = chi_refit_vprob->at(nRefitNumber);
//								refit_vProb2 = chi_refit_vprob->at(iChi2);
//							}
//							else continue;//skip those that don't have it
//							if (refit_vProb2 < 0.01) continue;
//							////////
//
//
//							double dimuonM2 = LVdimuon2->M();
//							double Mdiff2 = m_chi2;// -dimuonM2 + 3.097;// Assume J/psi mass
//							if (passDimSel2 == true && dimuonM2 > 2.9 && dimuonM2<3.3 && Mdiff2 > mass_window_l && Mdiff < mass_window_h)
//
//							{
//									double deltaEta = abs(conv_eta->at(iChi) - conv_eta->at(iChi2));
//									double deltapT = abs(conv_pt->at(iChi) - conv_pt->at(iChi2));
//
//
//									nConvCounterPeak++;
//									TVector3* TVconv1 = (TVector3*)conv_vtx->At(iChi);
//									TVector3* TVconv2 = (TVector3*)conv_vtx->At(iChi2);
//									double convVtx_deltaX = abs(TVconv1->X() - TVconv2->X());
//									double convVtx_deltaY = abs(TVconv1->Y() - TVconv2->Y());
//									double convVtx_deltaZ = abs(TVconv1->Z() - TVconv2->Z());
//
//									if (deltaEta < 0.02 && deltapT < 0.2 && convVtx_deltaX < 2.0 && convVtx_deltaY < 2.0 && convVtx_deltaZ < 2.0) {
//										//TLorentzVector* LVconv1 = (TLorentzVector*)conv_p4->At(iChi);
//										//TLorentzVector* LVconv2 = (TLorentzVector*)conv_p4->At(iChi2);
//										//TLorentzVector LVconvTot = *LVconv1;
//										//LVconvTot += *LVconv2;
//										//hconv_m->Fill(LVconvTot.M());
//										nConvSplitCounterPeak++;
//										//file_log << nchicCounterPassMass << " " << runNumber << " " << eventNumber << " " << chiCandPerEvent << " " << rap_chi << " " << pT_chi << " " << Mdiff << " m1 " << muon_eta->at(muon1Pos) << " " << muon_pt->at(muon1Pos)<< " m2 " << muon_eta->at(muon2Pos)<< " " << muon_pt->at(muon2Pos) << " conv " << convPos << " " << conv_eta->at(convPos) << " " << conv_pt->at(convPos) << " " << conv_vertexPositionRho->at(convPos) << " "  << TVconv1->X() << " " << TVconv1->Y() << " " << TVconv1->Z() << endl;
//										//file_log << nchicCounterPassMass << " " << runNumber << " " << eventNumber << " " << chiCandPerEvent << " " << rap_chi2 << " " << pT_chi2 << " " << Mdiff2 << " m1 " << muon_eta->at(muon1Pos2) << " " << muon_pt->at(muon1Pos2) << " m2 " << muon_eta->at(muon2Pos2) << " " << muon_pt->at(muon2Pos2) << " conv " << convPos2 << " " << conv_eta->at(convPos2) << " " << conv_pt->at(convPos2) << " " << conv_vertexPositionRho->at(convPos2) << " " << TVconv2->X() << " " << TVconv2->Y() << " " << TVconv2->Z() << endl << endl;
//
//									}
//
//								
//							}
//
//						}
//
//
//						*/
//
//
//
//
//
//
//					}
//
//				}
//				if (dimuon_charge->at(dimuonPos) != 0) { hSignal_SS->Fill(Mdiff); }
//
//				

			} // end of chic loop




			////////////////////////////////////
			////////      JPSI         /////////
			///////////////////////////////////////

			//for (int iJpsi = 0; iJpsi < dimuon_p4->GetEntriesFast(); iJpsi++) // Jpsi loop
			for (int iJpsi = 0; iJpsi < 0; iJpsi++) // Jpsi loop //debug
			{

				passDimSel = DimuonPassAllCuts(iJpsi);

				if (passDimSel == false)continue;


				// fill dimuon stuff
				LVdimuon = (TLorentzVector*)dimuon_p4->At(iJpsi);
				double dimuonM = LVdimuon->M();
				double rap_Jpsi = LVdimuon->Rapidity();
				if (runNumber < 285833) rap_Jpsi = -rap_Jpsi; //flip Pbp direction

				//roofit:
				if (dimuonM > mass_windowJpsi_l && dimuonM < mass_windowJpsi_h) {
					rvmassJpsi->setVal(dimuonM);
					rvptJpsi->setVal(dimuon_pt->at(iJpsi));
					rvrapJpsi->setVal(rap_Jpsi);
					rvntrackJpsi->setVal(pvtx_nTracks->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi))); //for now using total

					double accEff_Jpsi = 1;// hWeightJpsi->GetBinContent(hWeightJpsi->FindBin(abs(LVdimuon->Rapidity()), dimuon_pt->at(iJpsi)));
					//cout << "Test rap: " << rap_chi << "  pt  " << pT_chi << " and value accEff " << accEff_chi << endl;
					if (accEff_Jpsi > 0.00001) {// binning should ensure enough statistics, but if we have empty bin, just set weight to be 1
						rvweightJpsi->setVal(1 / accEff_Jpsi);
					}
					else rvweightJpsi->setVal(1);

					rvweightJpsi->setVal(1); // no weights for now

					if (isMC) {
						rvweightJpsi->setVal(1); //no weights if MC
					}


					rdsNominalJpsi->add(*colsJpsi, rvweightJpsi->getVal());




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


		CreateModelPdf(myWs, myPdfName);

		cout << "Gothere" << endl;

		RooFitResult* fitResult = myWs.pdf(myPdfName.c_str())->fitTo(*myWs.data("rdsNominal"), Extended(true), SumW2Error(true), Range(mass_windowFit_l, mass_windowFit_h), NumCPU(1), Save(true));
		//fitResult->Print("v");
		myWs.import(*fitResult, TString::Format("fitResult_%s", myPdfName.c_str()));
		cout << "Gothere2" << endl;

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




		FitRooDataSet(gAsChic_pT, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && rvrap>-2.4 && rvrap <2.4", "pt_all", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_pT_mid, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "pt_mid", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_pT_fwd, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap<-1 || rvrap >1)", "pt_fwd", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_y, bins_y, nbins_y, rvmass, myWs, false, "rvrap", myPdfName.c_str(), " && rvpt>6.5 && rvpt <30", "y", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_nTrk, bins_nTrk, nbins_nTrk, rvmass, myWs, false, "rvntrack", myPdfName.c_str(), " && rvrap>-1 && rvrap <1", "nTrack", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_nTrk_all, bins_nTrk, nbins_nTrk, rvmass, myWs, false, "rvntrack", myPdfName.c_str(), " && rvrap>-2.4 && rvrap <2.4", "nTrack_all", flagConstrainedFit, fileConstraints);

		FitRooDataSet(gAsChic_pT_fwdOnly, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap >1.6 && rvrap <2.4)", "pt_fwdOnly", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_pT_bkwOnly, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap <-1.6 && rvrap >-2.4)", "pt_bkwOnly", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_pT_fwdOnlyWide, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap >1.0 && rvrap <2.4)", "pt_fwdOnlyWide", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_pT_bkwOnlyWide, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap <-1.0 && rvrap >-2.4)", "pt_bkwOnlyWide", flagConstrainedFit, fileConstraints);

		FitRooDataSet(gAsChic_pT_midCMS, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap >-0.566 && rvrap <1.434)", "pt_midCMS", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_pT_fwdCMS, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap >1.434 && rvrap <2.4)", "pt_fwdCMS", flagConstrainedFit, fileConstraints);
		FitRooDataSet(gAsChic_pT_bkwCMS, bins_pT, nbins_pT, rvmass, myWs, false, "rvpt", myPdfName.c_str(), " && (rvrap <-0.566 && rvrap >-1.566)", "pt_bkwCMS", flagConstrainedFit, fileConstraints);

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
		myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("signalJpsi"), LineStyle(kDashed), LineColor(kRed));
		//myWs.pdf(myPdfNameJpsi.c_str())->plotOn(massframeJpsi, Components("psi2"), LineStyle(kDashed), LineColor(kGreen));
		cout << endl << endl << "HEREJpsi" << endl << endl;



		TCanvas *cTestJpsi = new TCanvas("cTestJpsi", "cTestJpsi", 1000, 600);
		gPad->SetLeftMargin(0.15);
		massframeJpsi->GetYaxis()->SetTitleOffset(1.3);
		massframeJpsi->Draw();


		cTestJpsi->SaveAs("FitterOutput/CanvasCTest_RW3_Jpsi.png");



		
		FitRooDataSet(gAsJpsi_pT_mid, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && rvrapJpsi>-1.0 && rvrapJpsi <1.0", "pt_mid", false);
		FitRooDataSet(gAsJpsi_pT_fwd, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi<-1.0 || rvrapJpsi >1.0)", "pt_fwd", false);
		FitRooDataSet(gAsJpsi_pT, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && rvrapJpsi>-2.4 && rvrapJpsi <2.4", "pt_all", false);
		FitRooDataSet(gAsJpsi_y, bins_y, nbins_y, rvmassJpsi, myWs, true, "rvrapJpsi", myPdfNameJpsi.c_str(), " && rvptJpsi>6.5 && rvptJpsi < 30", "y", false);
		FitRooDataSet(gAsJpsi_nTrk, bins_nTrk, nbins_nTrk, rvmassJpsi, myWs, true, "rvntrackJpsi", myPdfNameJpsi.c_str(), " && rvrapJpsi>-1 && rvrapJpsi <1", "nTrack", false);
		FitRooDataSet(gAsJpsi_nTrk_all, bins_nTrk, nbins_nTrk, rvmassJpsi, myWs, true, "rvntrackJpsi", myPdfNameJpsi.c_str(), " && rvrapJpsi>-2.4 && rvrapJpsi <2.4", "nTrack_all", false);

		FitRooDataSet(gAsJpsi_pT_fwdOnly, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi >1.6 && rvrapJpsi <2.4)", "pt_fwdOnly", false);
		FitRooDataSet(gAsJpsi_pT_bkwOnly, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi <-1.6 && rvrapJpsi >-2.4)", "pt_bkwOnly", false);
		FitRooDataSet(gAsJpsi_pT_fwdOnlyWide, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi >1.0 && rvrapJpsi <2.4)", "pt_fwdOnlyWide", false);
		FitRooDataSet(gAsJpsi_pT_bkwOnlyWide, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi <-1.0 && rvrapJpsi >-2.4)", "pt_bkwOnlyWide", false);

		FitRooDataSet(gAsJpsi_pT_midCMS, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi >-0.566 && rvrapJpsi <1.434)", "pt_midCMS", false);
		FitRooDataSet(gAsJpsi_pT_fwdCMS, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi >1.434 && rvrapJpsi <2.4)", "pt_fwdCMS", false);
		FitRooDataSet(gAsJpsi_pT_bkwCMS, bins_pT, nbins_pT, rvmassJpsi, myWs, true, "rvptJpsi", myPdfNameJpsi.c_str(), " && (rvrapJpsi <-0.566 && rvrapJpsi >-1.566)", "pt_bkwCMS", false);

	}
	
		////////////// DONE, GET RATIOS

		//pt
	///*
		GetCorrectedRatio(gAsRatio_pT, gAsChic_pT, gAsJpsi_pT, fileCorrection, "h_chiEfficiency1D_Q_pT_all_rat");
		GetCorrectedRatio(gAsRatio_pT_mid, gAsChic_pT_mid, gAsJpsi_pT_mid, fileCorrection, "h_chiEfficiency1D_Q_pT_mid_rat");
		GetCorrectedRatio(gAsRatio_pT_fwd, gAsChic_pT_fwd, gAsJpsi_pT_fwd, fileCorrection, "h_chiEfficiency1D_Q_pT_fwd_rat");

		GetCorrectedRatio(gAsRatio_pT_fwdOnly, gAsChic_pT_fwdOnly, gAsJpsi_pT_fwdOnly, fileCorrection, "h_chiEfficiency1D_Q_pT_fwdOnly_rat");
		GetCorrectedRatio(gAsRatio_pT_bkwOnly, gAsChic_pT_bkwOnly, gAsJpsi_pT_bkwOnly, fileCorrection, "h_chiEfficiency1D_Q_pT_bkwOnly_rat");
		GetCorrectedRatio(gAsRatio_pT_fwdOnlyWide, gAsChic_pT_fwdOnlyWide, gAsJpsi_pT_fwdOnlyWide, fileCorrection, "h_chiEfficiency1D_Q_pT_fwdOnlyWide_rat");
		GetCorrectedRatio(gAsRatio_pT_bkwOnlyWide, gAsChic_pT_bkwOnlyWide, gAsJpsi_pT_bkwOnlyWide, fileCorrection, "h_chiEfficiency1D_Q_pT_bkwOnlyWide_rat");

		GetCorrectedRatio(gAsRatio_pT_midCMS, gAsChic_pT_midCMS, gAsJpsi_pT_midCMS, fileCorrection, "h_chiEfficiency1D_Q_pT_midCMS_rat");
		GetCorrectedRatio(gAsRatio_pT_fwdCMS, gAsChic_pT_fwdCMS, gAsJpsi_pT_fwdCMS, fileCorrection, "h_chiEfficiency1D_Q_pT_fwdCMS_rat");
		GetCorrectedRatio(gAsRatio_pT_bkwCMS, gAsChic_pT_bkwCMS, gAsJpsi_pT_bkwCMS, fileCorrection, "h_chiEfficiency1D_Q_pT_bkwCMS_rat");

		TCanvas* can_pT = new TCanvas("can_pT", "Ratio_pT", 600, 400);
		gAsRatio_pT->SetMarkerStyle(21);
		gAsRatio_pT->SetMarkerSize(1.7);
		gAsRatio_pT->SetMarkerColor(kBlue);
		gAsRatio_pT->SetLineColor(kBlue);
		gAsRatio_pT->SetLineWidth(2);
		gAsRatio_pT->Draw("AP");
		gAsRatio_pT->GetYaxis()->SetRangeUser(0, 0.5);
		gAsRatio_pT_mid->SetMarkerStyle(20);
		gAsRatio_pT_mid->SetMarkerSize(2.0);
		gAsRatio_pT_mid->SetMarkerColor(kRed);
		gAsRatio_pT_mid->SetLineColor(kRed);
		gAsRatio_pT_mid->SetLineWidth(2);
		gAsRatio_pT_mid->Draw("P");
		gAsRatio_pT_fwd->SetMarkerStyle(29);
		gAsRatio_pT_fwd->SetMarkerSize(2.0);
		gAsRatio_pT_fwd->SetMarkerColor(kGreen);
		gAsRatio_pT_fwd->SetLineColor(kGreen);
		gAsRatio_pT_fwd->SetLineWidth(2);
		gAsRatio_pT_fwd->Draw("P");
		gAsRatio_pT->GetYaxis()->SetTitle("(#chic_{1}+#chic_{2}) / J/#psi");
		gAsRatio_pT->GetYaxis()->SetTitleOffset(1.1);
		gAsRatio_pT->GetXaxis()->SetTitle("pT");

		TLegend* leg_pT = new TLegend(0.7, 0.2, 0.88, 0.4, "");
		leg_pT->AddEntry(gAsRatio_pT, "Integrated", "p");
		leg_pT->AddEntry(gAsRatio_pT_mid, "Midrapidity |y|<1", "p");
		leg_pT->AddEntry(gAsRatio_pT_fwd, "Forward |y|>1", "p");
		leg_pT->Draw();

		can_pT->SaveAs("FitterOutput/RatioChicJpsi_pT.png");


		///*


		////////  y

		GetCorrectedRatio(gAsRatio_y, gAsChic_y, gAsJpsi_y, fileCorrection, "h_chiEfficiency1D_Q_y_rat");



		TCanvas* can_y = new TCanvas("can_y", "Ratio_y", 600, 400);
		gAsRatio_y->SetMarkerStyle(21);
		gAsRatio_y->SetMarkerSize(1.7);
		gAsRatio_y->SetMarkerColor(kRed);
		gAsRatio_y->SetLineColor(kRed);
		gAsRatio_y->SetLineWidth(2);
		gAsRatio_y->Draw("AP");
		gAsRatio_y->GetYaxis()->SetRangeUser(0, 0.5);
		gAsRatio_y->GetYaxis()->SetTitle("(#chic_{1}+#chic_{2}) / J/#psi");
		gAsRatio_y->GetYaxis()->SetTitleOffset(1.1);
		gAsRatio_y->GetXaxis()->SetTitle("rapidity");
		can_y->SaveAs("FitterOutput/RatioChicJpsi_y.png");



		GetCorrectedRatio(gAsRatio_nTrk, gAsChic_nTrk, gAsJpsi_nTrk, fileCorrection, "h_chiEfficiency1D_Q_nTrk_rat");
		gAsRatio_nTrk->RemovePoint(4); // not enough statistics there

		GetCorrectedRatio(gAsRatio_nTrk_all, gAsChic_nTrk_all, gAsJpsi_nTrk_all, fileCorrection, "h_chiEfficiency1D_Q_nTrk_all_rat");
		gAsRatio_nTrk_all->RemovePoint(4); // not enough statistics there

		TCanvas* can_nTrk = new TCanvas("can_nTrk", "Ratio_nTrk", 600, 400);
		gAsRatio_nTrk->SetMarkerStyle(21);
		gAsRatio_nTrk->SetMarkerSize(1.7);
		gAsRatio_nTrk->SetMarkerColor(kRed);
		gAsRatio_nTrk->SetLineColor(kRed);
		gAsRatio_nTrk->SetLineWidth(2);
		gAsRatio_nTrk->Draw("AP");
		gAsRatio_nTrk->GetYaxis()->SetRangeUser(0, 0.5);
		gAsRatio_nTrk->GetYaxis()->SetTitle("(#chic_{1}+#chic_{2}) / J/#psi");
		gAsRatio_nTrk->GetYaxis()->SetTitleOffset(1.1);
		gAsRatio_nTrk->GetXaxis()->SetTitle("nTrack");
		gAsRatio_nTrk->GetXaxis()->SetRangeUser(0,250);
		can_nTrk->SaveAs("FitterOutput/RatioChicJpsi_nTrk.png");

		//*/


		/*
		TCanvas* can1 = new TCanvas("can1", "plot", 1200, 800);

		TH2F* h1 = hdxy_dz->Clone();
		//cout<<h1->Integral(101, 300, 1, 30)<<endl;

		h1->GetXaxis()->SetTitle("abseta");
		h1->GetYaxis()->SetTitle("pt");
		can1->SetLogz();
		h1->Draw("colz");
		*/

	

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
		gAsChic_pT_mid->Write();
		gAsChic_pT_fwd->Write();
		gAsChic_pT_fwdOnly->Write();
		gAsChic_pT_bkwOnly->Write();
		gAsChic_pT_fwdOnlyWide->Write();
		gAsChic_pT_bkwOnlyWide->Write();
		gAsChic_pT_midCMS->Write();
		gAsChic_pT_fwdCMS->Write();
		gAsChic_pT_bkwCMS->Write();


		gAsJpsi_pT->Write();
		gAsJpsi_pT_mid->Write();
		gAsJpsi_pT_fwd->Write();
		gAsJpsi_pT_fwdOnly->Write();
		gAsJpsi_pT_bkwOnly->Write();
		gAsJpsi_pT_fwdOnlyWide->Write();
		gAsJpsi_pT_bkwOnlyWide->Write();
		gAsJpsi_pT_midCMS->Write();
		gAsJpsi_pT_fwdCMS->Write();
		gAsJpsi_pT_bkwCMS->Write();

		gAsRatio_pT->Write();
		gAsRatio_pT_mid->Write();
		gAsRatio_pT_fwd->Write();
		gAsRatio_pT_fwdOnly->Write();
		gAsRatio_pT_bkwOnly->Write();
		gAsRatio_pT_fwdOnlyWide->Write();
		gAsRatio_pT_bkwOnlyWide->Write();
		gAsRatio_pT_midCMS->Write();
		gAsRatio_pT_fwdCMS->Write();
		gAsRatio_pT_bkwCMS->Write();

		//can_pT->Write();
		gAsChic_y->Write();
		gAsJpsi_y->Write();
		gAsRatio_y->Write();
		//can_y->Write();
		gAsChic_nTrk->Write();
		gAsJpsi_nTrk->Write();
		gAsRatio_nTrk->Write();

		gAsChic_nTrk_all->Write();
		gAsJpsi_nTrk_all->Write();
		gAsRatio_nTrk_all->Write();
		//can_nTrk->Write();

		// dump all the fit parameters

		for (int i = 0; i < nFittingSets; i++)
		{
			for (int j = 0; j < nFitFunctionParams; j++)
			{
				if (gAsFitOutputArray[i][j] != NULL)
				{
					gAsFitOutputArray[i][j]->GetXaxis()->SetTitle(fittingSets.at(i).c_str());
					gAsFitOutputArray[i][j]->GetXaxis()->SetLabelSize(0.05);
					gAsFitOutputArray[i][j]->GetXaxis()->SetLabelSize(0.05);
					gAsFitOutputArray[i][j]->GetXaxis()->SetTitleSize(0.05);
					gAsFitOutputArray[i][j]->GetXaxis()->SetTitleOffset(1.05);
					gAsFitOutputArray[i][j]->SetMarkerStyle(21);
					gAsFitOutputArray[i][j]->SetMarkerSize(1.1);
					gAsFitOutputArray[i][j]->SetMarkerColor(kBlue);

					gAsFitOutputArray[i][j]->Write();
				}
			}
		}

		for (int i = 0; i < nFittingSets; i++)
		{
			for (int j = 0; j < nFitFunctionParams; j++)
			{
				if (gAsFitOutputArrayJpsi[i][j] != NULL)
				{
					gAsFitOutputArrayJpsi[i][j]->GetXaxis()->SetTitle(fittingSets.at(i).c_str());
					gAsFitOutputArrayJpsi[i][j]->GetXaxis()->SetLabelSize(0.05);
					gAsFitOutputArrayJpsi[i][j]->GetXaxis()->SetLabelSize(0.05);
					gAsFitOutputArrayJpsi[i][j]->GetXaxis()->SetTitleSize(0.05);
					gAsFitOutputArrayJpsi[i][j]->GetXaxis()->SetTitleOffset(1.05);
					gAsFitOutputArrayJpsi[i][j]->SetMarkerStyle(21);
					gAsFitOutputArrayJpsi[i][j]->SetMarkerSize(1.1);
					gAsFitOutputArrayJpsi[i][j]->SetMarkerColor(kBlue);

					gAsFitOutputArrayJpsi[i][j]->Write();
				}
			}
		}
		
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



