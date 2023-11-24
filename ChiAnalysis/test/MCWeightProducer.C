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

using namespace std;

double bins_ntrack[227] = {};
//double binsChiEffpT[] = { 6.5, 8.0, 10.0, 15.0, 20.0, 25.0 };
////double binsChiEffpT[] = { 0.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0 };
int  nBins_ntrack = sizeof(bins_ntrack) / sizeof(double) - 1;

int MCWeightProducer(const char* fileInMC = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC-Official_v3-bothDir.root", const char* fileInRD = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV-Nominal_v2-bothDir.root", const char* fileOut = "MCWeight_v2.root")
{
	gStyle->SetOptStat(1111);

	Long64_t totalCountMC = 0;
	Long64_t totalCountRD = 0;

	// create the binning
	for (int i = 0; i < 227; i++)
	{
		if (i < 201) // single track bins until 200
		{
			bins_ntrack[i] = i;
		}
		else if (i < 221) // bins of 2 until 240
		{
			bins_ntrack[i] = 200 + 2 * (i - 200);
		}
		else if (i<227) //bins of 5 until cutoff at 270
		{
			bins_ntrack[i] = 240 + 5 * (i - 220);
		}
	}
	for (int i = 0; i < 227; i++)
	{
		cout << bins_ntrack[i] << " ";
	}
	cout << endl;

	////weighting based on the number of tracks in the whole event
	//TH1D* h_distributionNtrackMC = new TH1D("h_distributionNtrackMC", "Distribution of ntrack in MC", 250, 0.0, 250);
	//TH1D* h_distributionNtrackRD = new TH1D("h_distributionNtrackRD", "Distribution of ntrack in RD", 250, 0.0, 250);
	//TH1D* h_weightNtrack = new TH1D("h_weightNtrack", "MC weight in ntrack binning", 250, 0.0, 250);
	//TH1D* h_distributionNtrackMC_weighted = new TH1D("h_distributionNtrackMC_weighted", "Distribution of ntrack in MC after weighting", 250, 0.0, 250); //cross check
	//TH1D* h_distributionNtrackMC_weightedScaled = new TH1D("h_distributionNtrackMC_weightedScaled", "Distribution of ntrack in MC after weighting, scaled to RD counts", 250, 0.0, 250); //cross check

	//// weighting based on the number of tracks in the associated PV
	//TH1D* h_distributionPVtrkMC = new TH1D("h_distributionPVtrkMC", "Distribution of PVtrk in MC", 250, 0.0, 250);
	//TH1D* h_distributionPVtrkRD = new TH1D("h_distributionPVtrkRD", "Distribution of PVtrk in RD", 250, 0.0, 250);
	//TH1D* h_weightPVtrk = new TH1D("h_weightPVtrk", "MC weight in PVtrk binning", 250, 0.0, 250);
	//TH1D* h_distributionPVtrkMC_weighted = new TH1D("h_distributionPVtrkMC_weighted", "Distribution of PVtrk in MC after weighting", 250, 0.0, 250); //cross check
	//TH1D* h_distributionPVtrkMC_weightedScaled = new TH1D("h_distributionPVtrkMC_weightedScaled", "Distribution of PVtrk in MC after weighting, scaled to RD counts", 250, 0.0, 250); //cross check

	//weighting based on the number of tracks in the whole event
	TH1D* h_distributionNtrackMC = new TH1D("h_distributionNtrackMC", "Distribution of ntrack in MC", nBins_ntrack, bins_ntrack);
	TH1D* h_distributionNtrackRD = new TH1D("h_distributionNtrackRD", "Distribution of ntrack in RD", nBins_ntrack, bins_ntrack);
	TH1D* h_weightNtrack = new TH1D("h_weightNtrack", "MC weight in ntrack binning", nBins_ntrack, bins_ntrack);
	TH1D* h_distributionNtrackMC_weighted = new TH1D("h_distributionNtrackMC_weighted", "Distribution of ntrack in MC after weighting", nBins_ntrack, bins_ntrack); //cross check
	TH1D* h_distributionNtrackMC_weightedScaled = new TH1D("h_distributionNtrackMC_weightedScaled", "Distribution of ntrack in MC after weighting, scaled to RD counts", nBins_ntrack, bins_ntrack); //cross check

	// weighting based on the number of tracks in the associated PV
	TH1D* h_distributionPVtrkMC = new TH1D("h_distributionPVtrkMC", "Distribution of PVtrk in MC", nBins_ntrack, bins_ntrack);
	TH1D* h_distributionPVtrkRD = new TH1D("h_distributionPVtrkRD", "Distribution of PVtrk in RD", nBins_ntrack, bins_ntrack);
	TH1D* h_weightPVtrk = new TH1D("h_weightPVtrk", "MC weight in PVtrk binning", nBins_ntrack, bins_ntrack);
	TH1D* h_distributionPVtrkMC_weighted = new TH1D("h_distributionPVtrkMC_weighted", "Distribution of PVtrk in MC after weighting", nBins_ntrack, bins_ntrack); //cross check
	TH1D* h_distributionPVtrkMC_weightedScaled = new TH1D("h_distributionPVtrkMC_weightedScaled", "Distribution of PVtrk in MC after weighting, scaled to RD counts", nBins_ntrack, bins_ntrack); //cross check

	// weighting based on the number of tracks in the associated PV
	TH1D* h_distributionPVtrkUncutMC = new TH1D("h_distributionPVtrkUncutMC", "Distribution of PVtrkUncut in MC", nBins_ntrack, bins_ntrack);
	TH1D* h_distributionPVtrkUncutRD = new TH1D("h_distributionPVtrkUncutRD", "Distribution of PVtrkUncut in RD", nBins_ntrack, bins_ntrack);
	TH1D* h_weightPVtrkUncut = new TH1D("h_weightPVtrkUncut", "MC weight in PVtrkUncut binning", nBins_ntrack, bins_ntrack);
	TH1D* h_distributionPVtrkUncutMC_weighted = new TH1D("h_distributionPVtrkUncutMC_weighted", "Distribution of PVtrkUncut in MC after weighting", nBins_ntrack, bins_ntrack); //cross check
	TH1D* h_distributionPVtrkUncutMC_weightedScaled = new TH1D("h_distributionPVtrkUncutMC_weightedScaled", "Distribution of PVtrkUncut in MC after weighting, scaled to RD counts", nBins_ntrack, bins_ntrack); //cross check





	int weird_decay_counter = 0;
	TLorentzVector* LVdimuon;

	///////////////////////////////////
	////////  S  T  A  R  T  //////////
	////////////////////////////////////




	////////  MC   /////////////


	TFile* fMC = new TFile(fileInMC, "READ");

	TTree* event_tree = (TTree*)fMC->Get("ChiRootuple/event_tree");
	if (!event_tree) { 
		cout << "Problem with event tree - MC";
		return 1;
	}

	LoadChiBranches(event_tree, true, true);


	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: "<<nentries << endl;
	//if (nentries > 2000000) { nentries = 20000; }

	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);

		//cout << "here" << endl;
		if (i % 10000 == 0) { cout << "MC event: " << i << " done: " << 100 * i / nentries << "%" << endl; }


		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			weird_decay_counter++;
			continue;
		}


		for (int iJpsi = 0; iJpsi < dimuon_p4->GetEntriesFast(); iJpsi++) // Jpsi loop
		{
			if (DimuonPassAllCuts(iJpsi) == false) continue;

			LVdimuon = (TLorentzVector*)dimuon_p4->At(iJpsi);
			double dimuonM = LVdimuon->M();

			//limit to likely signal
			if (dimuonM > mass_windowJpsi_l && dimuonM < mass_windowJpsi_h)
			{
				h_distributionNtrackMC->Fill(ntracks_inEvent);
				++totalCountMC;
				h_distributionPVtrkMC->Fill(pvtx_nTracks->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)));
				h_distributionPVtrkUncutMC->Fill(pvtx_nTracksUncut->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)));
			}

		}

	} //end of event loop

	cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;

	fMC->Close();

	////////  RD   /////////////

	TFile* fRD = new TFile(fileInRD, "READ");

	TTree* event_treeRD = (TTree*)fRD->Get("ChiRootuple/event_tree");
	if (!event_treeRD) {
		cout << "Problem with event tree - RD";
		return 1;
	}

	LoadChiBranches(event_treeRD, false, true);


	nentries = event_treeRD->GetEntries();
	cout << "n entries: " << nentries << endl;
	//if (nentries > 2000000) { nentries = 2000000; }

	for (Long64_t i = 0; i < nentries; i++) {

		event_treeRD->GetEntry(i);

		//cout << "here" << endl;
		if (i % 10000 == 0) { cout << "RD event: " << i << " done: " << 100 * i / nentries << "%" << endl; }

		for (int iJpsi = 0; iJpsi < dimuon_p4->GetEntriesFast(); iJpsi++) // Jpsi loop
		{
			if (DimuonPassAllCuts(iJpsi) == false) continue;

			LVdimuon = (TLorentzVector*)dimuon_p4->At(iJpsi);
			double dimuonM = LVdimuon->M();

			//limit to likely signal
			if (dimuonM > mass_windowJpsi_l && dimuonM < mass_windowJpsi_h)
			{
				h_distributionNtrackRD->Fill(ntracks_inEvent);
				++totalCountRD;
				h_distributionPVtrkRD->Fill(pvtx_nTracks->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)));
				h_distributionPVtrkUncutRD->Fill(pvtx_nTracksUncut->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)));
			}

		}

	}
	fRD->Close();
	
	// ntrack weighting
	h_distributionNtrackMC->Sumw2();
	h_distributionNtrackRD->Sumw2();

	h_weightNtrack->Divide(h_distributionNtrackRD, h_distributionNtrackMC, totalCountMC, totalCountRD); // normalized to the same number as in MC (i.e. average weight~1)

	//crosscheck
	for (int i = 1; i < h_distributionNtrackMC_weighted->GetNbinsX() + 2; i++)
	{
		h_distributionNtrackMC_weighted->SetBinContent(i, h_distributionNtrackMC->GetBinContent(i)*h_weightNtrack->GetBinContent(i));
	}
	h_distributionNtrackMC_weightedScaled=(TH1D*)h_distributionNtrackMC_weighted->Clone("h_distributionNtrackMC_weightedScaled");
	h_distributionNtrackMC_weightedScaled->SetTitle("Distribution of ntrack in MC after weighting, scaled to RD counts");
	h_distributionNtrackMC_weightedScaled->Scale((double)totalCountRD / (double) totalCountMC);

	// PV track weighting
	h_distributionPVtrkMC->Sumw2();
	h_distributionPVtrkRD->Sumw2();

	h_weightPVtrk->Divide(h_distributionPVtrkRD, h_distributionPVtrkMC, totalCountMC, totalCountRD); // normalized to the same number as in MC (i.e. average weight~1)

	//crosscheck
	for (int i = 1; i < h_distributionPVtrkMC_weighted->GetNbinsX() + 2; i++)
	{
		h_distributionPVtrkMC_weighted->SetBinContent(i, h_distributionPVtrkMC->GetBinContent(i)*h_weightPVtrk->GetBinContent(i));
	}
	h_distributionPVtrkMC_weightedScaled = (TH1D*)h_distributionPVtrkMC_weighted->Clone("h_distributionPVtrkMC_weightedScaled");
	h_distributionPVtrkMC_weightedScaled->SetTitle("Distribution of PVtrk in MC after weighting, scaled to RD counts");
	h_distributionPVtrkMC_weightedScaled->Scale((double)totalCountRD / (double)totalCountMC);

	// PV uncut track weighting
	h_distributionPVtrkUncutMC->Sumw2();
	h_distributionPVtrkUncutRD->Sumw2();

	h_weightPVtrkUncut->Divide(h_distributionPVtrkUncutRD, h_distributionPVtrkUncutMC, totalCountMC, totalCountRD); // normalized to the same number as in MC (i.e. average weight~1)

	//crosscheck
	for (int i = 1; i < h_distributionPVtrkUncutMC_weighted->GetNbinsX() + 2; i++)
	{
		h_distributionPVtrkUncutMC_weighted->SetBinContent(i, h_distributionPVtrkUncutMC->GetBinContent(i)*h_weightPVtrkUncut->GetBinContent(i));
	}
	h_distributionPVtrkUncutMC_weightedScaled = (TH1D*)h_distributionPVtrkUncutMC_weighted->Clone("h_distributionPVtrkUncutMC_weightedScaled");
	h_distributionPVtrkUncutMC_weightedScaled->SetTitle("Distribution of PVtrkUncut in MC after weighting, scaled to RD counts");
	h_distributionPVtrkUncutMC_weightedScaled->Scale((double)totalCountRD / (double)totalCountMC);



	TFile* fout = new TFile(fileOut, "RECREATE");

	h_distributionNtrackMC->Write();
	h_distributionNtrackRD->Write();
	h_weightNtrack->Write();
	h_distributionNtrackMC_weighted->Write();
	h_distributionNtrackMC_weightedScaled->Write();
	h_distributionPVtrkMC->Write();
	h_distributionPVtrkRD->Write();
	h_weightPVtrk->Write();
	h_distributionPVtrkMC_weighted->Write();
	h_distributionPVtrkMC_weightedScaled->Write();
	h_distributionPVtrkUncutMC->Write();
	h_distributionPVtrkUncutRD->Write();
	h_weightPVtrkUncut->Write();
	h_distributionPVtrkUncutMC_weighted->Write();
	h_distributionPVtrkUncutMC_weightedScaled->Write();
	fout->Close();
	cout << "END" << endl;
	return 0;
}
