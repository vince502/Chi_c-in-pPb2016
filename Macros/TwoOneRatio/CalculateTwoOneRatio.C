
// Macro to combine the fits and correction to calculate the chic2/chic1


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

TGraphAsymmErrors* gAsFitOutputArray_c2ratio[nFittingSets]; //stores the c2 ratio
TGraphAsymmErrors* gAsFitOutputArray_TwoOneRatio[nFittingSets]; // stores the chic2/chic1 ratio


int CorrectRatio_InverseCorrection(TGraphAsymmErrors* gAsResult, const char* fileCorrection, const char* corrName)  // for inverse correction (c1/c2), so it gets multiplied
{
	TFile* fCor = new TFile(fileCorrection, "READ");
	cout << "Opening correction: "<< corrName << endl;
	TH1D* hCor = (TH1D*)fCor->Get(corrName);
	cout << "Opened" << endl;
	for (int i = 0; i < hCor->GetNbinsX(); i++)
	{
		double corr = hCor->GetBinContent(i + 1); //TH1 numbering offset by 1
		cout << corr << endl;
		//gAsResult->GetY()[i] *= (1 / corr);
		//gAsResult->SetPointEYhigh(i, gAsResult->GetErrorYhigh(i) / corr);
		//gAsResult->SetPointEYlow(i, gAsResult->GetErrorYlow(i) / corr);
		gAsResult->GetY()[i] *= corr;
		gAsResult->SetPointEYhigh(i, gAsResult->GetErrorYhigh(i) * corr);
		gAsResult->SetPointEYlow(i, gAsResult->GetErrorYlow(i) * corr);
	}
	cout << "DoneCorrecting" << endl;
	fCor->Close();
	return 0;
}



///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////


void CalculateTwoOneRatio(const char* fileIn = "Chi_c_output_Nominal_vDissertation-bothDir_DCB_NewBins.root", const char* fileOut = "Chi_c_output_TwoOne_Nominal_v4-bothDir.root", const char* fileCorrection = "Chi_c_WeightsMC_Official_v4-bothDir.root")
{
	//gStyle->SetOptStat(1111);
	//gStyle->SetOptStat(0);
	setTDRStyle();





	TFile* f1 = new TFile(fileIn, "READ");

	///////////// LOAD THE RATIO

	for (int nBinSet = 0; nBinSet < nFittingSets; nBinSet++)
	{

		gAsFitOutputArray_c2ratio[nBinSet] = (TGraphAsymmErrors*)f1->Get(("gAsChic_FitOutput_c2ratio_" + fittingSets.at(nBinSet)).c_str()); // load the ratio
		gAsFitOutputArray_TwoOneRatio[nBinSet] = CalculateChicRatioFromC2Ratio(gAsFitOutputArray_c2ratio[nBinSet]); // recalculate the ratio
		cout << ("h_chiTotCorrChic1toChic2_1D_" + fittingSetsCorrection.at(nBinSet) + "_rat").c_str()<<endl;
		CorrectRatio_InverseCorrection(gAsFitOutputArray_TwoOneRatio[nBinSet], fileCorrection, ("h_chiTotCorrChic1toChic2_1D_" + fittingSetsCorrection.at(nBinSet) + "_rat").c_str());
	}

	gAsFitOutputArray_TwoOneRatio[1]->Draw();
	gAsFitOutputArray_c2ratio[1]->Draw("same");
	
	// instead of this, save it in a root file "fileOut"

}