
/////////// Obtain overall systematic uncertainty (combain individual) and make detailed syst plots. Produce output file to be used for final-results plots

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



#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../ChiTreeInit.C"
#include "../ChiFitterInit.h"

const double conversionSelectionSyst = 12;  // the systematic assigned to extra-loose/loose/medium conversion selection (flat)


using namespace std;
using namespace RooFit;


int addUpSystematics(TGraph* gResult, TGraph* gSyst1, TGraph* gSyst2) {
	if (gResult->GetN() != gSyst1->GetN() || gResult->GetN() != gSyst2->GetN()) {
		cout << "Error: Not the same number of N bins" << endl;
		return -1;
	}
	double syst1_x, syst1_y, syst2_x, syst2_y;

	for (int i = 0; i < gResult->GetN(); i++)
	{
		gSyst1->GetPoint(i, syst1_x, syst1_y);
		gSyst2->GetPoint(i, syst2_x, syst2_y);

		
		gResult->SetPoint(i, syst1_x, sqrt(syst1_y*syst1_y + syst2_y * syst2_y));

	}
	return 0;
}

	


///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////

void ObtainTotalSystematic()
{
	gStyle->SetOptStat(0);
	setTDRStyle();	
	gStyle->SetPadRightMargin(0.03);

	const char* fileOut = "TotalSyst_vDissertation.root";
	

	const char* fileInChic = "Chi_c_Syst_Total_vDissertation.root";
	const char* fileInJpsi = "Chi_c_Jpsi_Syst_Total_vDissertation.root";
	const char* fileInMCSyst = "Chi_c_MCSyst_vDissertation.root";

	TFile* fileChic = new TFile(fileInChic, "OPEN");

	TGraph* gSystChic_sig[nFittingSets];
	TGraph* gSystChic_bkg[nFittingSets];

	TFile* fileJpsi = new TFile(fileInJpsi, "OPEN");

	TGraph* gSystJpsi_sig[nFittingSets];
	TGraph* gSystJpsi_bkg[nFittingSets];

	TGraph* gSyst_Total[nFittingSets];
	for (int i = 0; i < nFittingSets; i++) { 
		gSyst_Total[i] = nullptr; //in principle should already be nullptr, but it turns out that it is instead some random stuff. Needs to be set to nullptr in case that we skip sets (for readout)
	}

	TFile* fileMCSyst = new TFile(fileInMCSyst, "OPEN");

	TGraph* gSyst_MC[nFittingSets];

	TGraph* gSyst_Selection[nFittingSets]; // create tgraph with the points for the values


	// MAIN LOOP OVER ALL THE SETS
	
	for (int iSet = 0; iSet < nFittingSets; iSet++) {

		//DON'T DO USELESS BINS
		//if (iSet > 2.5 && iSet < 8.5) continue; no longer in the fitting sets


		gSystChic_sig[iSet] = (TGraph*)fileChic->Get(("gSystOutputArraySig_" + fittingSets[iSet]).c_str());
		if (gSystChic_sig[iSet] == nullptr) cout << "failed to open" << endl;

		gSystChic_bkg[iSet] = (TGraph*)fileChic->Get(("gSystOutputArrayBkg_" + fittingSets[iSet]).c_str());
		if (gSystChic_bkg[iSet] == nullptr) cout << "failed to open" << endl;


		gSystJpsi_sig[iSet] = (TGraph*)fileJpsi->Get(("gSystOutputArraySig_" + fittingSets[iSet]).c_str());
		if (gSystJpsi_sig[iSet] == nullptr) cout << "failed to open" << endl;

		gSystJpsi_bkg[iSet] = (TGraph*)fileJpsi->Get(("gSystOutputArrayBkg_" + fittingSets[iSet]).c_str());
		if (gSystJpsi_bkg[iSet] == nullptr) cout << "failed to open" << endl;



		gSyst_MC[iSet] = (TGraph*)fileMCSyst->Get(("gSystOutputArrayMC_" + fittingSets[iSet]).c_str());
		if (gSyst_MC[iSet] == nullptr) cout << "failed to open" << endl;


		// loading done

		// create total graphs
		int nBinsInSet = gSystChic_sig[iSet]->GetN();


		TGraph* gSystChic_FitTotal[nFittingSets];
		gSystChic_FitTotal[iSet] = new TGraph(nBinsInSet);

		addUpSystematics(gSystChic_FitTotal[iSet], gSystChic_sig[iSet], gSystChic_bkg[iSet]);

		TGraph* gSystJpsi_FitTotal[nFittingSets];
		gSystJpsi_FitTotal[iSet] = new TGraph(nBinsInSet);

		addUpSystematics(gSystJpsi_FitTotal[iSet], gSystJpsi_sig[iSet], gSystJpsi_bkg[iSet]);


		gSyst_Total[iSet] = new TGraph(nBinsInSet);
		cout << iSet<< " gSystOutputArrayTotal_" + fittingSets.at(iSet) << endl;
		gSyst_Total[iSet]->SetNameTitle(("gSystOutputArrayTotal_" + fittingSets.at(iSet)).c_str(), ("Total syst uncertainty for " + fittingSets.at(iSet)).c_str());


		addUpSystematics(gSyst_Total[iSet], gSystChic_FitTotal[iSet], gSystJpsi_FitTotal[iSet]);

		addUpSystematics(gSyst_Total[iSet], gSyst_Total[iSet], gSyst_MC[iSet]); // add MC syst

		// conversion selection
		gSyst_Selection[iSet] = new TGraph(nBinsInSet);
		for (int j = 0; j < nBinsInSet; j++) {
			gSyst_Selection[iSet]->SetPoint(j, gSyst_Total[iSet]->GetPointX(j), conversionSelectionSyst);
		}

		addUpSystematics(gSyst_Total[iSet], gSyst_Total[iSet], gSyst_Selection[iSet]);

		//if (iSet == 5) //ntrack
		//{
		//	gSystChic_sig[iSet]->RemovePoint(4); //remove the test point at 250-400 (no stats)
		//	gSystChic_bkg[iSet]->RemovePoint(4);
		//	gSystChic_FitTotal[iSet]->RemovePoint(4);
		//	gSystJpsi_sig[iSet]->RemovePoint(4);
		//	gSystJpsi_bkg[iSet]->RemovePoint(4);
		//	gSystJpsi_FitTotal[iSet]->RemovePoint(4);
		//	gSyst_Selection[iSet]->RemovePoint(4);
		//	gSyst_MC[iSet]->RemovePoint(4);
		//	gSyst_Total[iSet]->RemovePoint(4);
		//
		//
		//}

		// PLOT

		cout << "Set number, xname, low and high edge: " << iSet << " " << fittingSetsXLabel[iSet] << " " << fittingSetsXLow[iSet] << " " << fittingSetsXHigh[iSet] << endl;

		TCanvas* CanvasResults;

		CanvasResults = new TCanvas("CanvasResult", "Systematic uncertainty results", 800, 640);

		TH1D* hFrame = new TH1D("hFrame", "Systematic uncertainty", 1, fittingSetsXLow[iSet], fittingSetsXHigh[iSet]);
		//	TH1D* hFrame = new TH1D("hFrame", "Systematic uncertainty", 1, 0, 30);

		hFrame->Draw();
		hFrame->GetYaxis()->SetRangeUser(0, 40);
		hFrame->GetYaxis()->SetTitle("Uncertainty [%]");
		//	hFrame->GetXaxis()->SetTitle("p_{T}");
		hFrame->GetXaxis()->SetTitleSize(0.05);
		hFrame->GetYaxis()->SetTitleSize(0.05);
		hFrame->GetXaxis()->SetTitleOffset(1.05);
		hFrame->GetXaxis()->SetTitle(fittingSetsXLabel[iSet].c_str());

		gSystChic_sig[iSet]->SetMarkerStyle(25);
		gSystChic_sig[iSet]->SetMarkerSize(1.5);
		gSystChic_sig[iSet]->SetMarkerColor(kRed);
		gSystChic_sig[iSet]->SetLineColor(kRed);
		gSystChic_sig[iSet]->SetLineWidth(2);
		gSystChic_sig[iSet]->SetLineStyle(2);
		gSystChic_sig[iSet]->Draw("PL");

		gSystChic_bkg[iSet]->SetMarkerStyle(25);
		gSystChic_bkg[iSet]->SetMarkerSize(1.5);
		gSystChic_bkg[iSet]->SetMarkerColor(kRed);
		gSystChic_bkg[iSet]->SetLineColor(kRed);
		gSystChic_bkg[iSet]->SetLineWidth(2);
		gSystChic_bkg[iSet]->SetLineStyle(3);
		gSystChic_bkg[iSet]->Draw("PL");

		gSystChic_FitTotal[iSet]->SetMarkerStyle(21);
		gSystChic_FitTotal[iSet]->SetMarkerSize(1.5);
		gSystChic_FitTotal[iSet]->SetMarkerColor(kRed);
		gSystChic_FitTotal[iSet]->SetLineColor(kRed);
		gSystChic_FitTotal[iSet]->SetLineWidth(3);
		gSystChic_FitTotal[iSet]->Draw("PL");


		gSystJpsi_sig[iSet]->SetMarkerStyle(24);
		gSystJpsi_sig[iSet]->SetMarkerSize(1.5);
		gSystJpsi_sig[iSet]->SetMarkerColor(kBlue);
		gSystJpsi_sig[iSet]->SetLineColor(kBlue);
		gSystJpsi_sig[iSet]->SetLineWidth(2);
		gSystJpsi_sig[iSet]->SetLineStyle(2);
		gSystJpsi_sig[iSet]->Draw("PL");

		gSystJpsi_bkg[iSet]->SetMarkerStyle(24);
		gSystJpsi_bkg[iSet]->SetMarkerSize(1.5);
		gSystJpsi_bkg[iSet]->SetMarkerColor(kBlue);
		gSystJpsi_bkg[iSet]->SetLineColor(kBlue);
		gSystJpsi_bkg[iSet]->SetLineWidth(2);
		gSystJpsi_bkg[iSet]->SetLineStyle(3);
		gSystJpsi_bkg[iSet]->Draw("PL");

		gSystJpsi_FitTotal[iSet]->SetMarkerStyle(20);
		gSystJpsi_FitTotal[iSet]->SetMarkerSize(1.5);
		gSystJpsi_FitTotal[iSet]->SetMarkerColor(kBlue);
		gSystJpsi_FitTotal[iSet]->SetLineColor(kBlue);
		gSystJpsi_FitTotal[iSet]->SetLineWidth(3);
		gSystJpsi_FitTotal[iSet]->Draw("PL");

		gSyst_Selection[iSet]->SetMarkerStyle(33);
		gSyst_Selection[iSet]->SetMarkerSize(1.8);
		gSyst_Selection[iSet]->SetMarkerColor(kOrange);
		gSyst_Selection[iSet]->SetLineColor(kOrange);
		gSyst_Selection[iSet]->SetLineWidth(3);
		gSyst_Selection[iSet]->Draw("PL");

		gSyst_MC[iSet]->SetMarkerStyle(34);
		gSyst_MC[iSet]->SetMarkerSize(1.8);
		gSyst_MC[iSet]->SetMarkerColor(kGreen + 2);
		gSyst_MC[iSet]->SetLineColor(kGreen+2);
		gSyst_MC[iSet]->SetLineWidth(3);
		gSyst_MC[iSet]->Draw("PL");

		gSyst_Total[iSet]->SetMarkerStyle(29);
		gSyst_Total[iSet]->SetMarkerSize(2.7);
		gSyst_Total[iSet]->SetMarkerColor(kBlack);
		gSyst_Total[iSet]->SetLineColor(kBlack);
		gSyst_Total[iSet]->SetLineWidth(3);
		gSyst_Total[iSet]->Draw("PL");

		TLatex latex_text;
		latex_text.SetNDC();
		latex_text.SetTextColor(kBlack);
		latex_text.SetTextFont(42);
		latex_text.SetTextAlign(11);
		latex_text.SetTextSize(0.045);
		if (fittingSets.at(iSet) == "pt_all") { // plotting details
			latex_text.DrawLatex(0.2, 0.76, "-2.4<y_{lab}(J/#psi)<2.4");
		}
		if (fittingSets.at(iSet) == "pt_midCMS") { // plotting details
			latex_text.DrawLatex(0.2, 0.76, "-1.0<y_{CM}(J/#psi)<1.0");
		}
		if (fittingSets.at(iSet) == "pt_bkwCMS") { // plotting details
			latex_text.DrawLatex(0.2, 0.76, "-2.0<y_{CM}(J/#psi)<-1.0");
		}
		if (fittingSets.at(iSet) == "pt_fwdCMS") { // plotting details
			latex_text.DrawLatex(0.2, 0.76, "1.0<y_{CM}(J/#psi)<1.9");
		}

		TLegend* leg1 = new TLegend(0.52, 0.64, 0.86, 0.92, "");
		leg1->SetFillStyle(0);
		leg1->SetFillColor(0);
		leg1->SetBorderSize(0);
		leg1->SetTextSize(0.032);
		leg1->AddEntry(gSystChic_sig[iSet], "Chi fitting - signal", "pl");
		leg1->AddEntry(gSystChic_bkg[iSet], "Chi fitting - background", "pl");
		leg1->AddEntry(gSystChic_FitTotal[iSet], "Chi fitting - overall", "pl");
		leg1->AddEntry(gSystJpsi_sig[iSet], "J/#psi fitting - signal", "pl");
		leg1->AddEntry(gSystJpsi_bkg[iSet], "J/#psi fitting - background", "pl");
		leg1->AddEntry(gSystJpsi_FitTotal[iSet], "J/#psi fitting - overall", "pl");
		leg1->AddEntry(gSyst_Selection[iSet], "Conversion selection", "pl");
		leg1->AddEntry(gSyst_MC[iSet], "Pythia settings", "pl");
		leg1->AddEntry(gSyst_Total[iSet], "Total systematic uncertainty", "pl");
		leg1->Draw();








		CanvasResults->SaveAs(("Plots/OverallSystematicUncertainty_" + fittingSets[iSet] + ".png").c_str());
		CanvasResults->SaveAs(("Plots/OverallSystematicUncertainty_" + fittingSets[iSet] + ".pdf").c_str());
		CanvasResults->SaveAs(("Plots/OverallSystematicUncertainty_" + fittingSets[iSet] + ".root").c_str());

		delete hFrame;
		delete leg1;
		delete CanvasResults;

	}


	fileChic->Close();
	fileJpsi->Close();
	fileMCSyst->Close();


	// READOUT TO FILE
	TFile* fout = new TFile(fileOut, "RECREATE");




	for (int i = 0; i < nFittingSets; i++)

	{
		//cout << i << endl;
		if (gSyst_Total[i] == NULL)
		{
			cout << i<<" Set not selected, not saved" << endl;

		}
		else {
			gSyst_Total[i]->Write();
		}

	}

	fout->Close();

}


