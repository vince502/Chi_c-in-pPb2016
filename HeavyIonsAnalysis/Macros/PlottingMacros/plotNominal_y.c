#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <iostream>
#include "Riostream.h"
#include <math.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <TROOT.h>
#include <cmath>
#include <TMatrixTSym.h>
#include <TGraphAsymmErrors.h>

#include "tdrstyle.C"
#include "CMS_lumi.C"


const TString sEffName= "Nominal";
//const float cLowX = 5, cHighX = 30;
//const TString sNameTag = "_pT_CMS";
const float cLowX = -2.4, cHighX = 2.4;
const TString sNameTag = "_y";
//const float cLowX = 0.0, cHighX = 270;
//const TString sNameTag = "_ntrk";
double systErrorWidth = 0.1;
double pointXShift = 0;

const float _markerSize = 2.4;
void RemoveXError(TGraphAsymmErrors* gAS);
void PrepareSystPlotting(TGraphAsymmErrors* gAS_Result_Syst, TGraphAsymmErrors* gAS_Result, TGraph* g_Systematics, double errorWidth); // adds the systematic uncertainty stored in percents as TGraph to our results (gAS_Result), and stores it in the separate gAS for plotting
void ShiftXPosition(TGraphAsymmErrors* gAS, double shiftSize); // move the points to avoid overlap
void PrintOutTable(TGraphAsymmErrors* gAS, TGraphAsymmErrors* gAS_Syst);

int plotNominal_y()
{

	setTDRStyle();
	gStyle->SetEndErrorSize(1);
	gStyle->SetPadRightMargin(0.03);

	//gStyle->SetOptStat(000000000);
	//gStyle->SetOptFit(0);
	//gStyle->SetEndErrorSize(5);
	//gStyle->SetLineWidth(2);
	////gStyle->SetErrorX(0); //remove bars along X

	//gStyle->SetLabelFont(62, "xyz");
	//gStyle->SetTitleFont(62, "xyzt");
	//gStyle->SetLabelSize(0.042, "xyz");
	//gStyle->SetCanvasBorderMode(0);
	//gStyle->SetCanvasColor(kWhite);
	//gStyle->SetFrameBorderMode(0);
	//gStyle->SetFrameFillColor(kWhite);
	//gStyle->SetPalette(1, 0);
	//gStyle->SetTitleSize(0.07, "t");
	//gStyle->SetLabelSize(0.042, "t");
	
	//TCanvas* cankres1 = new TCanvas("cankres1", "Canvas with results1", 640, 480);
	TCanvas* cankres1 = new TCanvas("cankres1", "Canvas with results1", 800, 600);
	cankres1->SetBottomMargin(0.145);
	TH1F* hframe = new TH1F("hframe", "", 1, cLowX, cHighX);
	hframe->Draw();
	hframe->GetYaxis()->SetRangeUser(0.0, 0.45);
	//hframe->GetXaxis()->SetTitle("p_{T} (J/#psi) [GeV/c]");
	hframe->GetXaxis()->SetTitle("y_{lab,p} (J/#psi)");
	//hframe->GetXaxis()->SetTitle("N_{tracks}");
	//hframe->GetXaxis()->SetTitleSize(0.05);
	hframe->GetXaxis()->SetTitleOffset(1.03);
	hframe->GetYaxis()->SetTitle("(#chi_{c1}+#chi_{c2}) / J/#psi");
	hframe->GetYaxis()->SetTitleSize(0.06);
	//hframe->GetYaxis()->SetTitleOffset(1.12);

	





	//TF1* one = new TF1("one", "1", cLowX, cHighX);
	//one->SetLineStyle(2);
	//one->SetLineWidth(2);
	//one->SetLineColor(kBlack);

	// Load data
	
	TFile* myFile1 = new TFile("Chi_c_output_Nominal_vDissertation-bothDir_DCB_NewBins.root", "READ");
	//TCanvas* c1 = (TCanvas*)myFile1->Get("tpTreeTrk/Trk_ntracksdep/fit_eff_plots/tag_hiBin_PLOT_TrackCuts_true"); //Get canvas with the plot
	//TGraphAsymmErrors* gAS_PbPb= (TGraphAsymmErrors*)c1->GetListOfPrimitives()->FindObject("hxy_fit_eff");

	//TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT");
	//TGraphAsymmErrors* gAS_Result2 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_mid");
	//TGraphAsymmErrors* gAS_Result3 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_fwd");

	//TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_midCMS");
	//TGraphAsymmErrors* gAS_Result2 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_fwdCMS");
	//TGraphAsymmErrors* gAS_Result3 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_bkwCMS");

	TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_y");
	//TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_nTrk");

		/////////// SYSTEMATICS

	TFile* myFile2 = new TFile("TotalSyst_vDissertation.root", "READ");
	TGraph* g_Systematics = (TGraph*)myFile2->Get("gSystOutputArrayTotal_y");
	//TGraph* g_Systematics2 = (TGraph*)myFile2->Get("gSystOutputArrayTotal_nTrack");

	TGraphAsymmErrors* gAS_Result_Syst = (TGraphAsymmErrors*)gAS_Result->Clone();
	//TGraphAsymmErrors* gAS_Result_Syst2 = (TGraphAsymmErrors*)gAS_Result2->Clone();

	PrepareSystPlotting(gAS_Result_Syst, gAS_Result, g_Systematics, systErrorWidth);
	//PrepareSystPlotting(gAS_Result_Syst2, gAS_Result2, g_Systematics2, systErrorWidth);


	cout << "Loading done" << endl;

	//Setup

	gAS_Result_Syst->SetLineWidth(2);
	gAS_Result_Syst->SetLineColor(kRed);
	gAS_Result_Syst->SetMarkerSize(_markerSize);
	gAS_Result_Syst->SetFillColor(kRed);
	//gAS_Result_Syst->SetFillStyle(3002);//3002 looks good as png, but poor in pdf
	gAS_Result_Syst->SetFillStyle(3013);

	gAS_Result->SetMarkerSize(0.9*_markerSize);
	gAS_Result->SetMarkerColor(kRed);
	gAS_Result->SetMarkerStyle(20);
	gAS_Result->SetLineColor(kRed);
	gAS_Result->SetLineWidth(2);

	//gAS_Result2->SetMarkerSize(0.9*_markerSize);
	//gAS_Result2->SetMarkerColor(kBlue);
	//gAS_Result2->SetMarkerStyle(25);
	//gAS_Result2->SetLineColor(kBlue);
	//gAS_Result2->SetLineWidth(2);

	//gAS_Result3->SetMarkerSize(1.1*_markerSize);
	//gAS_Result3->SetMarkerColor(kGreen+3);
	//gAS_Result3->SetMarkerStyle(30);
	//gAS_Result3->SetLineColor(kGreen+3);
	//gAS_Result3->SetLineWidth(2);

	//gAS_Result2->SetPointEYlow(3, 0.08);
	//gAS_Result2->SetPointEYhigh(3, 0.08);

	////kill x errors
	//RemoveXError(gAS_pp);
	//RemoveXError(gAS_pPb);
	//RemoveXError(gAS_PbPb);

	cankres1->cd();
	//cankres1->SetLogx();
	//one->Draw("Same");

	gAS_Result_Syst->Draw("2");

	gAS_Result->Draw("PE");
	//gAS_Result2->Draw("PE");
	//gAS_Result3->Draw("PE");

	// add CMS text
	//CMS_lumi(cankres1, 0, 22); //mid top
	CMS_lumi(cankres1, 0, 10); //left top

	//TLegend*leg = new TLegend(0.50, 0.20, 0.88, 0.45, sEffName+" efficiency");
	//TLegend*leg = new TLegend(0.58, 0.18, 0.88, 0.4, "");
	TLegend*leg = new TLegend(0.20, 0.18, 0.58, 0.32, "");
	leg->SetFillColor(kWhite);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);

	leg->AddEntry(gAS_Result, "6.5<p_{T}<30 GeV/c", "p");
	//leg->AddEntry(gAS_Result, "|y|<1.0, 6.5<p_{T}<30 GeV/c", "p");

	//leg->AddEntry(gAS_Result, "Midrapidity: |y|<1.0", "p");
	//leg->AddEntry(gAS_Result2, "Forward: 1.6<y<2.4", "p");
	//leg->AddEntry(gAS_Result3, "Backward: -2.4<y<-1.6", "p");

	//leg->AddEntry(gAS_Result2, "Forward:      1.0<y_{CMS}<2.0", "p");
	//leg->AddEntry(gAS_Result,  "Midrapidity: -1.0<y_{CMS}<1.0", "p");
	//leg->AddEntry(gAS_Result3, "Backward:   -2.0<y_{CMS}<-1.0", "p");

	//leg->AddEntry(gAS_Result, "|y|<2.4", "p");
	//leg->AddEntry(gAS_Result2, "|y|<1.0", "p");
	//leg->AddEntry(gAS_Result3, "1.0<|y|<2.4", "p");
	leg->Draw("same");

	//TBox* box_mid = new TBox(30, 50, 50, 80);
	//box_mid->SetFillStyle(3002);
	//box_mid->Draw("same");

	TLine* yCMS = new TLine(0.465, 0.3, 0.465, 0.45);
	yCMS->SetLineStyle(7);
	yCMS->SetLineWidth(4);
	yCMS->SetLineColor(kOrange+4);
	yCMS->Draw();

	TLine* yCMS2 = new TLine(0.465, 0.0, 0.465, 0.15);
	yCMS2->SetLineStyle(7);
	yCMS2->SetLineWidth(4);
	yCMS2->SetLineColor(kOrange + 4);
	yCMS2->Draw();

	TPaveText* pText1 = new TPaveText(0.66, 0.7, 0.88, 0.78, "NDC NB");
	pText1->AddText("p-going");
	pText1->SetTextSize(0.05);
	pText1->SetFillColor(0);
	pText1->Draw("");

	TPaveText* pText2 = new TPaveText(0.41, 0.7, 0.63, 0.78, "NDC NB");
	pText2->AddText("Pb-going");
	pText2->SetTextSize(0.05);
	pText2->SetFillColor(0);
	pText2->Draw("");

	// Add the paves displayng the ranges for CM fits
	TPave* pBackward = new TPave(-1.535, 0.12, -0.535, 0.28, 0, "NB");
	pBackward->SetFillColor(kGreen+3);
	pBackward->SetFillStyle(3006);
	pBackward->Draw();

	TPave* pMid = new TPave(-0.535, 0.12, 1.465, 0.28, 0, "NB");
	pMid->SetFillColor(kOrange + 1);
	pMid->SetFillStyle(3007);
	pMid->Draw();

	TPave* pForward = new TPave(1.465, 0.12, 2.4, 0.28, 0, "NB");
	pForward->SetFillColor(kBlue);
	pForward->SetFillStyle(3006);
	pForward->Draw();


	cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + "_boxes.root");
	cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + "_boxes.pdf");
	cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + "_boxes.png");

	PrintOutTable(gAS_Result, gAS_Result_Syst);

	//cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + ".root");
	//cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + ".pdf");
	//cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + ".png");

	return 0;
}


void RemoveXError(TGraphAsymmErrors* gAS)
{
	for (int i = 0; i < gAS->GetN(); i++)
	{
		gAS->SetPointEXlow(i, 0);
		gAS->SetPointEXhigh(i, 0);
	}
}
void PrepareSystPlotting(TGraphAsymmErrors* gAS_Result_Syst, TGraphAsymmErrors* gAS_Result, TGraph* g_Systematics, double errorWidth) {
	if (gAS_Result_Syst->GetN() <= g_Systematics->GetN()) { //smaller, because of nTrk dependence which has an extra bin
		for (int i = 0; i < gAS_Result_Syst->GetN(); i++) {
			if (gAS_Result_Syst->GetPointX(i) != g_Systematics->GetPointX(i)) {
				cout << "SYSTEMATICS AND NOMINAL HAVE DIFFERENT BINNING: " << gAS_Result_Syst->GetPointX(i) << " " << g_Systematics->GetPointX(i) << endl;
			}

			gAS_Result_Syst->SetPointEYhigh(i, g_Systematics->GetPointY(i)*0.01*gAS_Result->GetPointY(i)); //change from percent, and multiply by value
			gAS_Result_Syst->SetPointEYlow(i, g_Systematics->GetPointY(i)*0.01*gAS_Result->GetPointY(i)); //change from percent, and multiply by value
			gAS_Result_Syst->SetPointEXhigh(i, errorWidth);
			gAS_Result_Syst->SetPointEXlow(i, errorWidth);
		}

	}
	else { cout << "DIFFERENT NUMBER OF BINS BETWEEN SYST AND NOMINAL " << gAS_Result_Syst->GetN() << " " << g_Systematics->GetN() << endl; }
}

void ShiftXPosition(TGraphAsymmErrors* gAS, double shiftSize)
{
	for (int i = 0; i < gAS->GetN(); i++) {
		gAS->SetPointX(i, gAS->GetPointX(i) + shiftSize);
	}

}

void PrintOutTable(TGraphAsymmErrors* gAS, TGraphAsymmErrors* gAS_Syst) 
{
	cout << endl << "******************************" << endl;
	cout<< " Printing out table: " << endl;
	cout << "Graph name: " << gAS->GetName() << " and title: " << gAS->GetTitle() << endl;
	cout <<  "******************************" << endl << endl;
	//cout.precision(precision);
	//cout<<std::fixed;

	cout << "Bin  " << " Values" << endl;
	for (int i = 0; i < gAS->GetN(); i++)
	{
		cout << endl<< gAS->GetPointX(i); printf(" & $ %.3f \\pm %.3f \\pm %.3f $ ", gAS->GetPointY(i), gAS->GetErrorY(i), gAS_Syst->GetErrorY(i));
	}
	

	cout << endl << endl;
}