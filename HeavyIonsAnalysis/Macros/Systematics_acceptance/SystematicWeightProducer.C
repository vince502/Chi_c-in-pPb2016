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

#include "../../Macros/ChiTreeInit.C"
#include "../../Macros/ChiFitterInit.h"
#include "../../Macros/tdrstyle.C"
#include "../../Macros/CMS_lumi.C"


using namespace std;

double bins_ntrack[227] = {};
int  nBins_ntrack = sizeof(bins_ntrack) / sizeof(double) - 1;

const float _markerSize = 2.4;
const int testSize = 1000000; // for testing, maximum size of the data set to be used

int ObtainMCSpectrum(const char* fileIn, TH1D* hOutput)
{
	int totalCountMC = 0;
	TFile* fMC = new TFile(fileIn, "READ");

	TTree* event_tree = (TTree*)fMC->Get("ChiRootuple/event_tree");
	if (!event_tree) {
		cout << "Problem with event tree - MC";
		return 1;
	}

	LoadChiBranches(event_tree, true, false);


	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: " << nentries << endl;
	if (nentries > testSize) { nentries = testSize; }
	int weird_decay_counter = 0;

	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);

		//cout << "here" << endl;
		if (i % 10000 == 0) { cout << "MC event: " << i << " done: " << 100 * i / nentries << "%" << endl; }


		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			//cout << gen_Jpsi_pt->size() << gen_muon_pt->size() << gen_phot_pt->size() << gen_isGoodChicDecay->at(0) << endl;
			weird_decay_counter++;
			continue;
		}

		int chicPos = ChiPassAllCutsMC(0);

		if (chicPos > -1) // had a passing chic

		//int JpsiPos = DimuonPassAllCutsMC(0);
		//if (JpsiPos > -1) // had a passing chic
		{
			int JpsiPos = chi_daughterJpsi_position->at(chicPos);
			TLorentzVector* LVJpsiReco = (TLorentzVector*)dimuon_p4->At(JpsiPos);
			double rap_JpsiReco = LVJpsiReco->Rapidity();
			double pt_JpsiReco = dimuon_pt->at(JpsiPos);
			double eta_JpsiReco = dimuon_eta->at(JpsiPos);

			hOutput->Fill(pt_JpsiReco);
			totalCountMC++;
		}


	} //end of event loop

	cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;

	fMC->Close();

	return totalCountMC;
}

void ObtainRDMCRatio(TH1D* hRatioMC, TGraphAsymmErrors* gAS_Result, TH1D* hSpectrumMC)
{
	if (hRatioMC->GetNbinsX() != hSpectrumMC->GetNbinsX() || hRatioMC->GetNbinsX() != gAS_Result->GetN())
	{
		cout << "Different number of points, will crash" << endl;
		return;
	}

	for (int i = 0; i < hRatioMC->GetNbinsX(); i++)
	{
		hRatioMC->SetBinContent(i + 1, gAS_Result->GetY()[i] / hSpectrumMC->GetBinContent(i + 1));
		double relativegASError = (gAS_Result->GetEYhigh()[i] + gAS_Result->GetEYlow()[i]) / (2 * gAS_Result->GetY()[i]); //average it
		//cout << "relative error gAS: " << relativegASError << endl;
		double relativeSpectrumError = hSpectrumMC->GetBinError(i + 1) / hSpectrumMC->GetBinContent(i + 1);
		hRatioMC->SetBinError(i + 1, sqrt(relativegASError*relativegASError + relativeSpectrumError * relativeSpectrumError) * hRatioMC->GetBinContent(i + 1));
	}

}

double ObtainPhotonGenSpectrum(const char* fileIn, TH1D* hOutput, TH1D* h_weight) 
{
	double totalCountMC = 0; //uses weights
	TFile* fMC = new TFile(fileIn, "READ");

	TTree* event_tree = (TTree*)fMC->Get("ChiRootuple/event_tree");
	if (!event_tree) {
		cout << "Problem with event tree - MC";
		return 1;
	}

	LoadChiBranches(event_tree, true, false);


	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: " << nentries << endl;
	if (nentries > testSize) { nentries = testSize; }
	int weird_decay_counter = 0;

	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);

		//cout << "here" << endl;
		if (i % 10000 == 0) { cout << "MC event: " << i << " done: " << 100 * i / nentries << "%" << endl; }


		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			//cout << gen_Jpsi_pt->size() << gen_muon_pt->size() << gen_phot_pt->size() << gen_isGoodChicDecay->at(0) << endl;
			weird_decay_counter++;
			continue;
		}
		double etaM1 = gen_muon_eta->at(0);
		double ptM1 = gen_muon_pt->at(0);
		double etaM2 = gen_muon_eta->at(1);
		double ptM2 = gen_muon_pt->at(1);
		double ptGen_Jpsi = gen_Jpsi_pt->at(0);
		TLorentzVector* LVJpsi = (TLorentzVector*)gen_Jpsi_p4->At(0);
		double rapGen_Jpsi = LVJpsi->Rapidity();

		//if (abs(rapGen_Jpsi) < 1.6) continue;
		//if (ptGen_Jpsi < 12) continue;

		if (MuonAcceptance(etaM1, ptM1)==true && MuonAcceptance(etaM2, ptM2)==true && DimuonAcceptance(rapGen_Jpsi, ptGen_Jpsi)==true) // good J/psi
		{
			
			double eta_phot = gen_phot_eta->at(0);
			double pt_phot = gen_phot_pt->at(0);

			double pTweight = h_weight->GetBinContent(h_weight->FindBin(ptGen_Jpsi));
			//double pTweight = WeightForMC_pTpart(ptGen_Jpsi); Not used, since we want to weight all the MC's
			hOutput->Fill(pt_phot, pTweight);
			totalCountMC+= pTweight;
		}


	} //end of event loop

	cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;

	fMC->Close();
	return totalCountMC;
}


int SystematicWeightProducer()
{
	const char* fileDataResults = "Chi_c_output_Nominal_v2-bothDir_DCB.root";
	const char* fileCorrection = "Chi_c_WeightsMC_Official_v3-bothDir.root";

	const char* fileInMC = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC-Official_v3-bothDir.root";
	const char* fileInMCSyst1 = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Systematics/Chi_c_pPb8TeV_MC-Systematics_3.root";
	const char* fileInMCSyst2 = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Systematics/Chi_c_pPb8TeV_MC-Systematics_6.root";
	const char* fileInMCSyst3 = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Systematics/Chi_c_pPb8TeV_MC-Systematics_cMass.root";
	const char* fileInMCSyst4 = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Systematics/Chi_c_pPb8TeV_MC-Systematics_reNorm.root";

	TH1D* hSpectrumMC_pt = new TH1D("hSpectrumMC_pt", "pT distribution of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hSpectrumMC_ptSyst1 = new TH1D("hSpectrumMC_ptSyst1", "pT distribution of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hSpectrumMC_ptSyst2 = new TH1D("hSpectrumMC_ptSyst2", "pT distribution of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hSpectrumMC_ptSyst3 = new TH1D("hSpectrumMC_ptSyst3", "pT distribution of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hSpectrumMC_ptSyst4 = new TH1D("hSpectrumMC_ptSyst4", "pT distribution of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);

	TH1D* hSpectrumMC_finePt = new TH1D("hSpectrumMC_finePt", "pT distribution of chic; p_{T} (J/#psi)", 60, 0, 30);
	TH1D* hSpectrumMC_finePtSyst1 = new TH1D("hSpectrumMC_finePtSyst1", "pT distribution of chic; p_{T} (J/#psi)", 60, 0, 30);
	TH1D* hSpectrumMC_finePtSyst2 = new TH1D("hSpectrumMC_finePtSyst2", "pT distribution of chic; p_{T} (J/#psi)", 60, 0, 30);
	TH1D* hSpectrumMC_finePtSyst3 = new TH1D("hSpectrumMC_finePtSyst3", "pT distribution of chic; p_{T} (J/#psi)", 60, 0, 30);
	TH1D* hSpectrumMC_finePtSyst4 = new TH1D("hSpectrumMC_finePtSyst4", "pT distribution of chic; p_{T} (J/#psi)", 60, 0, 30);

	TH1D* hRatioMC_pt = new TH1D("hRatioMC_pt", "pT RD/MC ratio of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hRatioMC_ptSyst1 = new TH1D("hRatioMC_ptSyst1", "pT RD/MC ratio of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hRatioMC_ptSyst2 = new TH1D("hRatioMC_ptSyst2", "pT RD/MC ratio of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hRatioMC_ptSyst3 = new TH1D("hRatioMC_ptSyst3", "pT RD/MC ratio of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);
	TH1D* hRatioMC_ptSyst4 = new TH1D("hRatioMC_ptSyst4", "pT RD/MC ratio of chic; p_{T} (J/#psi)", nbins_pT, bins_pT);

	// photon

	TH1D* hSpectrumPhotonMC_pt = new TH1D("hSpectrumPhotonMC_pt", "pT distribution of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hSpectrumPhotonMC_ptSyst1 = new TH1D("hSpectrumPhotonMC_ptSyst1", "pT distribution of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hSpectrumPhotonMC_ptSyst2 = new TH1D("hSpectrumPhotonMC_ptSyst2", "pT distribution of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hSpectrumPhotonMC_ptSyst3 = new TH1D("hSpectrumPhotonMC_ptSyst3", "pT distribution of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hSpectrumPhotonMC_ptSyst4 = new TH1D("hSpectrumPhotonMC_ptSyst4", "pT distribution of photon; p_{T} (#gamma)", 25, 0, 5);

	//TH1D* hRatioPhotonMC_pt = new TH1D("hRatioPhotonMC_pt", "pT RD/MC ratio of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hRatioPhotonMC_ptSyst1 = new TH1D("hRatioPhotonMC_ptSyst1", "pT RD/MC ratio of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hRatioPhotonMC_ptSyst2 = new TH1D("hRatioPhotonMC_ptSyst2", "pT RD/MC ratio of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hRatioPhotonMC_ptSyst3 = new TH1D("hRatioPhotonMC_ptSyst3", "pT RD/MC ratio of photon; p_{T} (#gamma)", 25, 0, 5);
	TH1D* hRatioPhotonMC_ptSyst4 = new TH1D("hRatioPhotonMC_ptSyst4", "pT RD/MC ratio of photon; p_{T} (#gamma)", 25, 0, 5);



	Long64_t totalCountMC = 0;
	Long64_t totalCountRD = 0;

	gStyle->SetOptStat(0);
	setTDRStyle();

	///////////
	// check the spectrum in the fiducial region 
	//////////


	////////  RD   /////////////

	TFile* fRD = new TFile(fileDataResults, "READ");
	TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)fRD->Get("gAsChic_pT");
	//TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)fRD->Get("gAsJpsi_pT");
	for (int i = 0; i < gAS_Result->GetN(); i++)
	{
		totalCountRD += gAS_Result->GetY()[i];
	}

	//TFile* fCor = new TFile(fileCorrection, "READ");
	//const char* corrName = "h_chiEfficiency1D_AnalysisBinning_rat"; 
	//cout << "Opening correction: " << corrName << endl;
	//TH1D* hCor = (TH1D*)fCor->Get(corrName);
	//cout << "Opened" << endl;
	//for (int i = 0; i < hCor->GetNbinsX(); i++)
	//{
	//	double corr = hCor->GetBinContent(i + 1); //TH1 numbering offset by 1
	//	cout << corr << endl;
	//	gAS_Result->GetY()[i] *= (1 / corr);
	//	gAS_Result->SetPointEYhigh(i, gAS_Result->GetErrorYhigh(i) / corr);
	//	gAS_Result->SetPointEYlow(i, gAS_Result->GetErrorYlow(i) / corr);
	//	totalCountRD += gAS_Result->GetY()[i];
	//}
	//cout << "DoneCorrecting, total number: " << totalCountRD << endl;
	//fCor->Close();

	///////////////  MC  ////////////

	totalCountMC = ObtainMCSpectrum(fileInMC, hSpectrumMC_pt);
	hSpectrumMC_pt->Sumw2();
	//hSpectrumMC_pt->Scale((double)gAS_Result->GetY()[1] / (double)hSpectrumMC_pt->GetBinContent(2));
	hSpectrumMC_pt->Scale((double)totalCountRD / (double)totalCountMC);

	totalCountMC = ObtainMCSpectrum(fileInMCSyst1, hSpectrumMC_ptSyst1);
	hSpectrumMC_ptSyst1->Sumw2();
	//hSpectrumMC_ptSyst1->Scale((double)gAS_Result->GetY()[1] / (double)hSpectrumMC_ptSyst1->GetBinContent(2));
	hSpectrumMC_ptSyst1->Scale((double)totalCountRD / (double)totalCountMC);

	totalCountMC = ObtainMCSpectrum(fileInMCSyst2, hSpectrumMC_ptSyst2);
	hSpectrumMC_ptSyst2->Sumw2();
	//hSpectrumMC_ptSyst2->Scale((double)gAS_Result->GetY()[1] / (double)hSpectrumMC_ptSyst2->GetBinContent(2));
	hSpectrumMC_ptSyst2->Scale((double)totalCountRD / (double)totalCountMC);

	totalCountMC = ObtainMCSpectrum(fileInMCSyst3, hSpectrumMC_ptSyst3);
	hSpectrumMC_ptSyst3->Sumw2();
	//hSpectrumMC_ptSyst3->Scale((double)gAS_Result->GetY()[1] / (double)hSpectrumMC_ptSyst3->GetBinContent(2));
	hSpectrumMC_ptSyst3->Scale((double)totalCountRD / (double)totalCountMC);

	totalCountMC = ObtainMCSpectrum(fileInMCSyst4, hSpectrumMC_ptSyst4);
	hSpectrumMC_ptSyst4->Sumw2();
	//hSpectrumMC_ptSyst4->Scale((double)gAS_Result->GetY()[1] / (double)hSpectrumMC_ptSyst4->GetBinContent(2));
	hSpectrumMC_ptSyst4->Scale((double)totalCountRD / (double)totalCountMC);

	ObtainRDMCRatio(hRatioMC_pt, gAS_Result, hSpectrumMC_pt);
	ObtainRDMCRatio(hRatioMC_ptSyst1, gAS_Result, hSpectrumMC_ptSyst1);
	ObtainRDMCRatio(hRatioMC_ptSyst2, gAS_Result, hSpectrumMC_ptSyst2);
	ObtainRDMCRatio(hRatioMC_ptSyst3, gAS_Result, hSpectrumMC_ptSyst3);
	ObtainRDMCRatio(hRatioMC_ptSyst4, gAS_Result, hSpectrumMC_ptSyst4);


	//totalCountMC = ObtainMCSpectrum(fileInMC, hSpectrumMC_pt);
	//hSpectrumMC_pt->Sumw2();
	//hSpectrumMC_pt->Scale((double)totalCountRD / (double)totalCountMC);

	//totalCountMC = ObtainMCSpectrum(fileInMCSyst1, hSpectrumMC_ptSyst1);
	//hSpectrumMC_ptSyst1->Sumw2();
	//hSpectrumMC_ptSyst1->Scale((double)totalCountRD / (double)totalCountMC);

	//totalCountMC = ObtainMCSpectrum(fileInMCSyst2, hSpectrumMC_ptSyst2);
	//hSpectrumMC_ptSyst2->Sumw2();
	//hSpectrumMC_ptSyst2->Scale((double)totalCountRD / (double)totalCountMC);

	//totalCountMC = ObtainMCSpectrum(fileInMCSyst3, hSpectrumMC_ptSyst3);
	//hSpectrumMC_ptSyst3->Sumw2();
	//hSpectrumMC_ptSyst3->Scale((double)totalCountRD / (double)totalCountMC);

	//totalCountMC = ObtainMCSpectrum(fileInMCSyst4, hSpectrumMC_ptSyst4);
	//hSpectrumMC_ptSyst4->Sumw2();
	//hSpectrumMC_ptSyst4->Scale((double)totalCountRD / (double)totalCountMC);

	// plot corrected yields
	TCanvas* cankres1 = new TCanvas("cankres1", "Canvas with results1", 900, 720);
	const float cLowX = 5, cHighX = 30;

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1, 1);
	cankres1->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.45*0.98);
	pad1->SetBottomMargin(0.04);
	pad2->SetTopMargin(0);
	pad2->SetFillColor(0);
	pad2->SetFillStyle(0);
	pad2->SetBottomMargin(gStyle->GetPadBottomMargin() / 0.55);
	pad1->SetTopMargin(gStyle->GetPadTopMargin() / 0.7);

	pad1->Draw();
	pad1->cd();


	TH1F* hframe = new TH1F("hframe", "", 1, cLowX, cHighX);
	hframe->Draw();
	hframe->GetYaxis()->SetRangeUser(0, 1.8*gAS_Result->GetY()[0]); //scale to the first bin in data
	hframe->GetXaxis()->SetTitle("");
	hframe->GetXaxis()->SetLabelSize(0); // switch it off
	//hframe->GetXaxis()->SetTitle("y (J/#psi)");
	//hframe->GetXaxis()->SetTitle("N_{tracks}");
	//hframe->GetXaxis()->SetTitleSize(0.05);
	//hframe->GetXaxis()->SetTitleOffset(1.12);
	hframe->GetYaxis()->SetTitle("Raw (#chi_{c1}+#chi_{c2}) yield ");
	//hframe->GetYaxis()->SetTitle("Raw J/#psi yield ");
	hframe->GetYaxis()->SetTitleSize(0.08);
	hframe->GetYaxis()->SetLabelSize(0.07);
	hframe->GetYaxis()->SetTitleOffset(0.9);

	//cankres1->cd();

	gAS_Result->SetMarkerSize(0.9*_markerSize);
	gAS_Result->SetMarkerColor(kRed);
	gAS_Result->SetMarkerStyle(21);
	gAS_Result->SetLineColor(kRed);
	gAS_Result->SetLineWidth(2);

	gAS_Result->Draw("ZPE");

	hSpectrumMC_pt->SetMarkerSize(0.9*_markerSize);
	hSpectrumMC_pt->SetMarkerColor(kBlue);
	hSpectrumMC_pt->SetMarkerStyle(20);
	hSpectrumMC_pt->SetLineColor(kBlue);
	hSpectrumMC_pt->SetLineWidth(2);
	hSpectrumMC_pt->Draw("same");

	hSpectrumMC_ptSyst1->SetMarkerSize(0.9*_markerSize);
	hSpectrumMC_ptSyst1->SetMarkerColor(kGreen);
	hSpectrumMC_ptSyst1->SetMarkerStyle(24);
	hSpectrumMC_ptSyst1->SetLineColor(kGreen);
	hSpectrumMC_ptSyst1->SetLineWidth(2);
	hSpectrumMC_ptSyst1->Draw("same");

	hSpectrumMC_ptSyst2->SetMarkerSize(0.9*_markerSize);
	hSpectrumMC_ptSyst2->SetMarkerColor(kGreen+3);
	hSpectrumMC_ptSyst2->SetMarkerStyle(24);
	hSpectrumMC_ptSyst2->SetLineColor(kGreen+3);
	hSpectrumMC_ptSyst2->SetLineWidth(2);
	hSpectrumMC_ptSyst2->Draw("same");

	hSpectrumMC_ptSyst3->SetMarkerSize(0.9*_markerSize);
	hSpectrumMC_ptSyst3->SetMarkerColor(kOrange);
	hSpectrumMC_ptSyst3->SetMarkerStyle(24);
	hSpectrumMC_ptSyst3->SetLineColor(kOrange);
	hSpectrumMC_ptSyst3->SetLineWidth(2);
	hSpectrumMC_ptSyst3->Draw("same");

	hSpectrumMC_ptSyst4->SetMarkerSize(0.9*_markerSize);
	hSpectrumMC_ptSyst4->SetMarkerColor(kCyan);
	hSpectrumMC_ptSyst4->SetMarkerStyle(24);
	hSpectrumMC_ptSyst4->SetLineColor(kCyan);
	hSpectrumMC_ptSyst4->SetLineWidth(2);
	hSpectrumMC_ptSyst4->Draw("same");

	//CMS_lumi(cankres1, 0, 10); //left top // or bottom, depending on version

	TLegend*leg = new TLegend(0.50, 0.54, 0.88, 0.90, "");
	leg->SetFillColor(kWhite);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.085);

	leg->AddEntry(gAS_Result, "Raw yield, data", "p");
	leg->AddEntry(hSpectrumMC_pt, "MC nominal", "p");
	leg->AddEntry(hSpectrumMC_ptSyst1, "MC Syst1 - 3", "p");
	leg->AddEntry(hSpectrumMC_ptSyst2, "MC Syst2 - 6 ", "p");
	leg->AddEntry(hSpectrumMC_ptSyst3, "MC Syst3 - cMass", "p");
	leg->AddEntry(hSpectrumMC_ptSyst4, "MC Syst4 - reNorm", "p");

	leg->Draw("same");

	cankres1->cd();
	pad2->Draw();
	pad2->cd();
	

	TH1F* hframe2 = new TH1F("hframe2", "", 1, cLowX, cHighX);
	hframe2->Draw();
	hframe2->GetYaxis()->SetRangeUser(0.0,3.0); //range of ratios
	hframe2->GetXaxis()->SetTitle("p_{T} (J/#psi) [GeV/c]");
	hframe2->GetXaxis()->SetTitleSize(0.10);
	hframe2->GetXaxis()->SetLabelSize(0.10);
	hframe2->GetYaxis()->SetTitle("Ratio RD/MC ");
	hframe2->GetYaxis()->SetTitleSize(0.09);
	hframe2->GetYaxis()->SetLabelSize(0.09);
	hframe2->GetYaxis()->SetTitleOffset(0.8);

	TLine *l1 = new TLine(cLowX, 1, cHighX, 1);
	l1->SetLineStyle(9);
	l1->Draw("same");

	hRatioMC_pt->SetMarkerSize(0.9*_markerSize);
	hRatioMC_pt->SetMarkerColor(kBlue);
	hRatioMC_pt->SetMarkerStyle(20);
	hRatioMC_pt->SetLineColor(kBlue);
	hRatioMC_pt->SetLineWidth(2);
	hRatioMC_pt->Draw("sameP");

	hRatioMC_ptSyst1->SetMarkerSize(0.9*_markerSize);
	hRatioMC_ptSyst1->SetMarkerColor(kGreen);
	hRatioMC_ptSyst1->SetMarkerStyle(24);
	hRatioMC_ptSyst1->SetLineColor(kGreen);
	hRatioMC_ptSyst1->SetLineWidth(2);
	hRatioMC_ptSyst1->Draw("sameP");

	hRatioMC_ptSyst2->SetMarkerSize(0.9*_markerSize);
	hRatioMC_ptSyst2->SetMarkerColor(kGreen + 3);
	hRatioMC_ptSyst2->SetMarkerStyle(24);
	hRatioMC_ptSyst2->SetLineColor(kGreen + 3);
	hRatioMC_ptSyst2->SetLineWidth(2);
	hRatioMC_ptSyst2->Draw("sameP");

	hRatioMC_ptSyst3->SetMarkerSize(0.9*_markerSize);
	hRatioMC_ptSyst3->SetMarkerColor(kOrange);
	hRatioMC_ptSyst3->SetMarkerStyle(24);
	hRatioMC_ptSyst3->SetLineColor(kOrange);
	hRatioMC_ptSyst3->SetLineWidth(2);
	hRatioMC_ptSyst3->Draw("sameP");

	hRatioMC_ptSyst4->SetMarkerSize(0.9*_markerSize);
	hRatioMC_ptSyst4->SetMarkerColor(kCyan);
	hRatioMC_ptSyst4->SetMarkerStyle(24);
	hRatioMC_ptSyst4->SetLineColor(kCyan);
	hRatioMC_ptSyst4->SetLineWidth(2);
	hRatioMC_ptSyst4->Draw("sameP");

	//gStyle->SetOptFit(1011); // doesn't work for some reason

	///////EXP

	TF1* fit_pT_weight = new TF1("fit_pT_weight", "[0]*(1-exp([1]*x))", 6.5, 30);
	fit_pT_weight->SetLineWidth(2);
	fit_pT_weight->SetParNames("Norm", "Exp_scale");
	//basic parameter settings
	fit_pT_weight->SetParameters(1, -0.05);
	TFitResultPtr FitRes = hRatioMC_pt->Fit(fit_pT_weight, "", "", 6.5, 30);
	double normValue = fit_pT_weight->GetParameter(0);
	double normError = fit_pT_weight->GetParError(0);
	double scaleValue = fit_pT_weight->GetParameter(1);
	double scaleError = fit_pT_weight->GetParError(1);
	TPaveText* pText1 = new TPaveText(.25, .74, .45, .96, "brNDC");
	pText1->SetFillColor(kWhite);
	pText1->SetTextSize(0.06);
	pText1->SetTextColor(kRed);
	pText1->SetBorderSize(0);
	pText1->SetTextAlign(32); //right alignment
	pText1->AddText(Form("Normalization: %.2f +- %.2f", normValue, normError));
	pText1->AddText(Form("Scale %.3f +- %.3f", scaleValue, scaleError));
	pText1->Draw();

	/////// Line might be better
	TF1* fit_pT_weight2 = new TF1("fit_pT_weight2", "[0]+[1]*x", 6.5, 30);
	fit_pT_weight2->SetLineWidth(2);
	fit_pT_weight2->SetLineColor(kBlue);
	fit_pT_weight2->SetParNames("A0", "A1");
	//basic parameter settings
	fit_pT_weight2->SetParameters(1, 0.05);
	TFitResultPtr FitRes2 = hRatioMC_pt->Fit(fit_pT_weight2, "S", "", 6.5, 30);
	double A0Value = fit_pT_weight2->GetParameter(0);
	double A0Error = fit_pT_weight2->GetParError(0);
	double A1Value = fit_pT_weight2->GetParameter(1);
	double A1Error = fit_pT_weight2->GetParError(1);
	TPaveText* pText2 = new TPaveText(.65, .24, .85, .46, "brNDC");
	pText2->SetFillColor(kWhite);
	pText2->SetTextSize(0.06);
	pText2->SetTextColor(kBlue);
	pText2->SetBorderSize(0);
	pText2->SetTextAlign(32); //right alignment
	pText2->AddText(Form("Offset: %.2f +- %.2f", A0Value, A0Error));
	pText2->AddText(Form("Slope %.3f +- %.3f", A1Value, A1Error));
	pText2->Draw();

	cankres1->SaveAs((TString)"pTdistribution_fiducial_comparison" + ".root");
	cankres1->SaveAs((TString)"pTdistribution_fiducial_comparison" + ".pdf");
	cankres1->SaveAs((TString)"pTdistribution_fiducial_comparison" + ".png");


	///////////////////////////////////////
	//////////
	/// THE ACTUAL WEIGHTS
	////////////
	///////////////////////////////////////


	double totalCountMCNomPhot = 0;
	double totalCountMCPhot = 0;

	totalCountMCNomPhot = ObtainPhotonGenSpectrum(fileInMC, hSpectrumPhotonMC_pt, hRatioMC_pt);
	//hSpectrumPhotonMC_pt->Sumw2();

	totalCountMCPhot = ObtainPhotonGenSpectrum(fileInMCSyst1, hSpectrumPhotonMC_ptSyst1, hRatioMC_ptSyst1);
	//hSpectrumPhotonMC_ptSyst1->Sumw2();
	hSpectrumPhotonMC_ptSyst1->Scale(totalCountMCNomPhot / totalCountMCPhot);
	//hSpectrumPhotonMC_ptSyst1->Scale((double)hSpectrumPhotonMC_pt->GetBinContent(6) / (double)hSpectrumPhotonMC_ptSyst1->GetBinContent(6));

	totalCountMCPhot = ObtainPhotonGenSpectrum(fileInMCSyst2, hSpectrumPhotonMC_ptSyst2, hRatioMC_ptSyst2);
	//hSpectrumPhotonMC_ptSyst2->Sumw2();
	//hSpectrumPhotonMC_ptSyst2->Scale((double)hSpectrumPhotonMC_pt->GetBinContent(6) / (double)hSpectrumPhotonMC_ptSyst2->GetBinContent(6));
	hSpectrumPhotonMC_ptSyst2->Scale(totalCountMCNomPhot / totalCountMCPhot);

	totalCountMCPhot = ObtainPhotonGenSpectrum(fileInMCSyst3, hSpectrumPhotonMC_ptSyst3, hRatioMC_ptSyst3);
	//hSpectrumPhotonMC_ptSyst3->Sumw2();
	//hSpectrumPhotonMC_ptSyst3->Scale((double)hSpectrumPhotonMC_pt->GetBinContent(6) / (double)hSpectrumPhotonMC_ptSyst3->GetBinContent(6));
	hSpectrumPhotonMC_ptSyst3->Scale(totalCountMCNomPhot / totalCountMCPhot);

	totalCountMCPhot = ObtainPhotonGenSpectrum(fileInMCSyst4, hSpectrumPhotonMC_ptSyst4, hRatioMC_ptSyst4);
	//hSpectrumPhotonMC_ptSyst4->Sumw2();
	//hSpectrumPhotonMC_ptSyst4->Scale((double)hSpectrumPhotonMC_pt->GetBinContent(6) / (double)hSpectrumPhotonMC_ptSyst4->GetBinContent(6));
	hSpectrumPhotonMC_ptSyst4->Scale(totalCountMCNomPhot / totalCountMCPhot);


	// ratios

	hRatioPhotonMC_ptSyst1->Divide(hSpectrumPhotonMC_ptSyst1, hSpectrumPhotonMC_pt);
	hRatioPhotonMC_ptSyst2->Divide(hSpectrumPhotonMC_ptSyst2, hSpectrumPhotonMC_pt);
	hRatioPhotonMC_ptSyst3->Divide(hSpectrumPhotonMC_ptSyst3, hSpectrumPhotonMC_pt);
	hRatioPhotonMC_ptSyst4->Divide(hSpectrumPhotonMC_ptSyst4, hSpectrumPhotonMC_pt);


	// plot corrected yields
	TCanvas* cankres2 = new TCanvas("cankres2", "Canvas with results1", 900, 720);
	const float cLowX_phot = 0, cHighX_phot = 5;

	TPad *pad3 = new TPad("pad3", "pad3", 0, 0.45, 1, 1);
	cankres2->cd();
	TPad *pad4 = new TPad("pad4", "pad4", 0, 0, 1, 0.45*0.98);
	pad3->SetBottomMargin(0.04);
	pad4->SetTopMargin(0);
	pad4->SetFillColor(0);
	pad4->SetFillStyle(0);
	pad4->SetBottomMargin(gStyle->GetPadBottomMargin() / 0.55);
	pad3->SetTopMargin(gStyle->GetPadTopMargin() / 0.7);

	pad3->Draw();
	pad3->cd();


	TH1F* hframe3 = new TH1F("hframe3", "", 1, cLowX_phot, cHighX_phot);
	hframe3->Draw();
	hframe3->GetYaxis()->SetRangeUser(0, 1.3*hSpectrumPhotonMC_pt->GetBinContent(1)); //scale to the bin in nominal MC
	hframe3->GetXaxis()->SetTitle("");
	hframe3->GetXaxis()->SetLabelSize(0); // switch it off
	//hframe3->GetXaxis()->SetTitle("y (J/#psi)");
	//hframe3->GetXaxis()->SetTitle("N_{tracks}");
	//hframe3->GetXaxis()->SetTitleSize(0.05);
	//hframe3->GetXaxis()->SetTitleOffset(1.12);
	//hframe3->GetYaxis()->SetTitle("Raw (#chi_{c1}+#chi_{c2}) yield ");
	hframe3->GetYaxis()->SetTitle("Gen photon spectrum");
	hframe3->GetYaxis()->SetTitleSize(0.08);
	hframe3->GetYaxis()->SetLabelSize(0.07);
	hframe3->GetYaxis()->SetTitleOffset(0.9);

	//cankres2->cd();


	hSpectrumPhotonMC_pt->SetMarkerSize(0.9*_markerSize);
	hSpectrumPhotonMC_pt->SetMarkerColor(kBlue);
	hSpectrumPhotonMC_pt->SetMarkerStyle(20);
	hSpectrumPhotonMC_pt->SetLineColor(kBlue);
	hSpectrumPhotonMC_pt->SetLineWidth(2);
	hSpectrumPhotonMC_pt->Draw("same");

	hSpectrumPhotonMC_ptSyst1->SetMarkerSize(0.9*_markerSize);
	hSpectrumPhotonMC_ptSyst1->SetMarkerColor(kGreen);
	hSpectrumPhotonMC_ptSyst1->SetMarkerStyle(24);
	hSpectrumPhotonMC_ptSyst1->SetLineColor(kGreen);
	hSpectrumPhotonMC_ptSyst1->SetLineWidth(2);
	hSpectrumPhotonMC_ptSyst1->Draw("same");

	hSpectrumPhotonMC_ptSyst2->SetMarkerSize(0.9*_markerSize);
	hSpectrumPhotonMC_ptSyst2->SetMarkerColor(kGreen + 3);
	hSpectrumPhotonMC_ptSyst2->SetMarkerStyle(24);
	hSpectrumPhotonMC_ptSyst2->SetLineColor(kGreen + 3);
	hSpectrumPhotonMC_ptSyst2->SetLineWidth(2);
	hSpectrumPhotonMC_ptSyst2->Draw("same");

	hSpectrumPhotonMC_ptSyst3->SetMarkerSize(0.9*_markerSize);
	hSpectrumPhotonMC_ptSyst3->SetMarkerColor(kOrange);
	hSpectrumPhotonMC_ptSyst3->SetMarkerStyle(24);
	hSpectrumPhotonMC_ptSyst3->SetLineColor(kOrange);
	hSpectrumPhotonMC_ptSyst3->SetLineWidth(2);
	hSpectrumPhotonMC_ptSyst3->Draw("same");

	hSpectrumPhotonMC_ptSyst4->SetMarkerSize(0.9*_markerSize);
	hSpectrumPhotonMC_ptSyst4->SetMarkerColor(kCyan);
	hSpectrumPhotonMC_ptSyst4->SetMarkerStyle(24);
	hSpectrumPhotonMC_ptSyst4->SetLineColor(kCyan);
	hSpectrumPhotonMC_ptSyst4->SetLineWidth(2);
	hSpectrumPhotonMC_ptSyst4->Draw("same");

	//CMS_lumi(cankres2, 0, 10); //left top // or bottom, depending on version

	TLegend*leg2 = new TLegend(0.50, 0.54, 0.88, 0.90, "");
	leg2->SetFillColor(kWhite);
	leg2->SetBorderSize(0);
	leg2->SetTextFont(42);
	leg2->SetTextSize(0.085);

	//leg2->AddEntry(gAS_Result, "Raw yield, data", "p");
	leg2->AddEntry(hSpectrumPhotonMC_pt, "MC nominal", "p");
	leg2->AddEntry(hSpectrumPhotonMC_ptSyst1, "MC Syst1 - 3", "p");
	leg2->AddEntry(hSpectrumPhotonMC_ptSyst2, "MC Syst2 - 6 ", "p");
	leg2->AddEntry(hSpectrumPhotonMC_ptSyst3, "MC Syst3 - cMass", "p");
	leg2->AddEntry(hSpectrumPhotonMC_ptSyst4, "MC Syst4 - reNorm", "p");

	leg2->Draw("same");





	cankres2->cd();
	pad4->Draw();
	pad4->cd();


	TH1F* hframe4 = new TH1F("hframe4", "", 1, cLowX_phot, cHighX_phot);
	hframe4->Draw();
	hframe4->GetYaxis()->SetRangeUser(0.0, 2.0); //range of ratios
	hframe4->GetXaxis()->SetTitle("p_{T} (#gamma) [GeV/c]");
	hframe4->GetXaxis()->SetTitleSize(0.10);
	hframe4->GetXaxis()->SetLabelSize(0.10);
	hframe4->GetYaxis()->SetTitle("Ratio Syst/Nom ");
	hframe4->GetYaxis()->SetTitleSize(0.09);
	hframe4->GetYaxis()->SetLabelSize(0.09);
	hframe4->GetYaxis()->SetTitleOffset(0.8);

	TLine *l2 = new TLine(cLowX_phot, 1, cHighX_phot, 1);
	l2->SetLineStyle(9);
	l2->Draw("same");


	hRatioPhotonMC_ptSyst1->SetMarkerSize(0.9*_markerSize);
	hRatioPhotonMC_ptSyst1->SetMarkerColor(kGreen);
	hRatioPhotonMC_ptSyst1->SetMarkerStyle(24);
	hRatioPhotonMC_ptSyst1->SetLineColor(kGreen);
	hRatioPhotonMC_ptSyst1->SetLineWidth(2);
	hRatioPhotonMC_ptSyst1->Draw("sameP");

	hRatioPhotonMC_ptSyst2->SetMarkerSize(0.9*_markerSize);
	hRatioPhotonMC_ptSyst2->SetMarkerColor(kGreen + 3);
	hRatioPhotonMC_ptSyst2->SetMarkerStyle(24);
	hRatioPhotonMC_ptSyst2->SetLineColor(kGreen + 3);
	hRatioPhotonMC_ptSyst2->SetLineWidth(2);
	hRatioPhotonMC_ptSyst2->Draw("sameP");

	hRatioPhotonMC_ptSyst3->SetMarkerSize(0.9*_markerSize);
	hRatioPhotonMC_ptSyst3->SetMarkerColor(kOrange);
	hRatioPhotonMC_ptSyst3->SetMarkerStyle(24);
	hRatioPhotonMC_ptSyst3->SetLineColor(kOrange);
	hRatioPhotonMC_ptSyst3->SetLineWidth(2);
	hRatioPhotonMC_ptSyst3->Draw("sameP");

	hRatioPhotonMC_ptSyst4->SetMarkerSize(0.9*_markerSize);
	hRatioPhotonMC_ptSyst4->SetMarkerColor(kCyan);
	hRatioPhotonMC_ptSyst4->SetMarkerStyle(24);
	hRatioPhotonMC_ptSyst4->SetLineColor(kCyan);
	hRatioPhotonMC_ptSyst4->SetLineWidth(2);
	hRatioPhotonMC_ptSyst4->Draw("sameP");


	// fitting
	hRatioPhotonMC_ptSyst1->Fit("pol2");
	TF1 *fitFunc_syst1 = hRatioPhotonMC_ptSyst1->GetFunction("pol2");
	fitFunc_syst1->SetLineColor(kGreen);
	fitFunc_syst1->SetLineWidth(2);
	TPaveText* pText_syst1 = new TPaveText(.25, .90, .58, .96, "brNDC");
	pText_syst1->SetFillColor(kWhite);
	pText_syst1->SetTextSize(0.055);
	pText_syst1->SetTextColor(kGreen);
	pText_syst1->SetBorderSize(0);
	pText_syst1->AddText(Form("Syst1: p2 %.3f +- %.3f; p1 %.2f +- %.2f; p0 %.2f +- %.2f", fitFunc_syst1->GetParameter(2), fitFunc_syst1->GetParError(2), fitFunc_syst1->GetParameter(1), fitFunc_syst1->GetParError(1), fitFunc_syst1->GetParameter(0), fitFunc_syst1->GetParError(0)));
	pText_syst1->Draw();

	hRatioPhotonMC_ptSyst2->Fit("pol2");
	TF1 *fitFunc_syst2 = hRatioPhotonMC_ptSyst2->GetFunction("pol2");
	fitFunc_syst2->SetLineColor(kGreen+3);
	fitFunc_syst2->SetLineWidth(2);
	TPaveText* pText_syst2 = new TPaveText(.25, .82, .58, .88, "brNDC");
	pText_syst2->SetFillColor(kWhite);
	pText_syst2->SetTextSize(0.055);
	pText_syst2->SetTextColor(kGreen+3);
	pText_syst2->SetBorderSize(0);
	pText_syst2->AddText(Form("Syst2: p2 %.3f +- %.3f; p1 %.2f +- %.2f; p0 %.2f +- %.2f", fitFunc_syst2->GetParameter(2), fitFunc_syst2->GetParError(2), fitFunc_syst2->GetParameter(1), fitFunc_syst2->GetParError(1), fitFunc_syst2->GetParameter(0), fitFunc_syst2->GetParError(0)));
	pText_syst2->Draw();

	hRatioPhotonMC_ptSyst3->Fit("pol2");
	TF1 *fitFunc_syst3 = hRatioPhotonMC_ptSyst3->GetFunction("pol2");
	fitFunc_syst3->SetLineColor(kOrange);
	fitFunc_syst3->SetLineWidth(2);
	TPaveText* pText_syst3 = new TPaveText(.25, .34, .58, .40, "brNDC");
	pText_syst3->SetFillColor(kWhite);
	pText_syst3->SetTextSize(0.055);
	pText_syst3->SetTextColor(kOrange);
	pText_syst3->SetBorderSize(0);
	pText_syst3->AddText(Form("Syst3: p2 %.3f +- %.3f; p1 %.2f +- %.2f; p0 %.2f +- %.2f", fitFunc_syst3->GetParameter(2), fitFunc_syst3->GetParError(2), fitFunc_syst3->GetParameter(1), fitFunc_syst3->GetParError(1), fitFunc_syst3->GetParameter(0), fitFunc_syst3->GetParError(0)));
	pText_syst3->Draw();

	hRatioPhotonMC_ptSyst4->Fit("pol2");
	TF1 *fitFunc_syst4 = hRatioPhotonMC_ptSyst4->GetFunction("pol2");
	fitFunc_syst4->SetLineColor(kCyan);
	fitFunc_syst4->SetLineWidth(2);
	TPaveText* pText_syst4 = new TPaveText(.25, .26, .58, .32, "brNDC");
	pText_syst4->SetFillColor(kWhite);
	pText_syst4->SetTextSize(0.055);
	pText_syst4->SetTextColor(kCyan);
	pText_syst4->SetBorderSize(0);
	pText_syst4->AddText(Form("Syst4: p2 %.3f +- %.3f; p1 %.2f +- %.2f; p0 %.2f +- %.2f", fitFunc_syst4->GetParameter(2), fitFunc_syst4->GetParError(2), fitFunc_syst4->GetParameter(1), fitFunc_syst4->GetParError(1), fitFunc_syst4->GetParameter(0), fitFunc_syst4->GetParError(0)));
	pText_syst4->Draw();



	cankres2->SaveAs((TString)"pTdistribution_photon_weighted" + ".root");
	cankres2->SaveAs((TString)"pTdistribution_photon_weighted" + ".pdf");
	cankres2->SaveAs((TString)"pTdistribution_photon_weighted" + ".png");


	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_fwdrap1624" + ".root");
	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_fwdrap1624" + ".pdf");
	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_fwdrap1624" + ".png");

	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_highpt1230" + ".root");
	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_highpt1230" + ".pdf");
	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_highpt1230" + ".png");

	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_lowpt6p512" + ".root");
	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_lowpt6p512" + ".pdf");
	//cankres2->SaveAs((TString)"pTdistribution_photon_weighted_lowpt6p512" + ".png");



	return 0;
}