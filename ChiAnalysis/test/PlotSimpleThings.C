
// Macro to plot simple distributions and similar stuff, to be run on eos


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

#include "../../Macros/tdrstyle.C"
#include "../../Macros/CMS_lumi.C"
#include "../../Macros/ChiTreeInit.C"
#include "../../Macros/ChiFitterInit.h"

void PlotSimpleThings()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TH1D* h_ntracksInEvent = new TH1D("h_ntracksInEvent", "h_ntracksInEvent", 100, 0, 300);
	TH1D* h_ntracksInPVwithDimuon = new TH1D("h_ntracksInPVwithDimuon", "h_ntracksInPVwithDimuon", 100, 0, 300);
	TH1D* h_ntracksInPVwithDimuonM = new TH1D("h_ntracksInPVwithDimuonM", "h_ntracksInPVwithDimuonM", 100, 0, 300);

	TH1D* h_PV_z = new TH1D("h_PV_z", "h_PV_z; z[cm]", 120, -30, 30);
	TH2D* h_PV_xy = new TH2D("h_PV_xy", "h_PV_xy; x[cm]; y[cm]", 200, -1, 1, 200, -1, 1);

	//TFile* f1 = new TFile("../../../../Data/Chi_c_pPb8TeV-Nominal_v2-bothDir.root", "READ");
	TFile* f1 = new TFile("/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV-Nominal_v2-bothDir.root", "READ");

	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");
	if (!event_tree) {
		cout << "Problem with event Tree";
		return;
	}
	LoadChiBranches(event_tree, false);

	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: " << nentries << endl;
	if (nentries > 500000) { nentries = 3000000; }

	
	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);
		if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }

		for (int j = 0; j < nPrimVertices; j++) {
			h_PV_z->Fill(pvtx_z->at(j));
			h_PV_xy->Fill(pvtx_x->at(j), pvtx_y->at(j));
		}
		for (int iJpsi = 0; iJpsi < dimuon_p4->GetEntriesFast(); iJpsi++) // Jpsi loop
		{
			bool passDimSel = DimuonPassAllCuts(iJpsi);

			if (passDimSel == false) continue;

			TLorentzVector* LVdimuon = (TLorentzVector*)dimuon_p4->At(iJpsi);
			double dimuonM = LVdimuon->M();

			if (dimuonM > mass_windowJpsi_l && dimuonM < mass_windowJpsi_h)
			{
				h_ntracksInEvent->Fill(ntracks_inEvent);
				h_ntracksInPVwithDimuon->Fill(pvtx_nTracks->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)));
				h_ntracksInPVwithDimuonM->Fill(pvtx_nTracks->at(dimuon_pvtx_index->at(iJpsi)));

				//h_PV_z->Fill(pvtx_z->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)));
				//h_PV_xy->Fill(pvtx_x->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)), pvtx_y->at(dimuon_pvtx_indexFromOniaMuMu->at(iJpsi)));
			}
		}


	}

	//*/

	TCanvas* cankres1 = new TCanvas("cankres1", "Canvas with results1", 800, 640);
	cankres1->cd();

	//event_tree->Draw("ntracks_inEvent");
	h_ntracksInEvent->SetLineWidth(2);
	h_ntracksInEvent->SetLineColor(kRed);
	h_ntracksInEvent->SetMarkerColor(kRed);
	h_ntracksInEvent->Draw();

	h_ntracksInEvent->GetYaxis()->SetTitle("Counts");
	h_ntracksInEvent->GetYaxis()->SetTitleSize(0.05);
	h_ntracksInEvent->GetYaxis()->SetLabelSize(0.045);
	h_ntracksInEvent->GetYaxis()->SetTitleOffset(1.05);

	h_ntracksInEvent->GetXaxis()->SetTitle("N_{tracks} in event/dimuon PV");
	h_ntracksInEvent->GetXaxis()->SetRangeUser(0, 300);
	h_ntracksInEvent->GetXaxis()->SetTitleSize(0.05);
	h_ntracksInEvent->GetXaxis()->SetLabelSize(0.045);
	h_ntracksInEvent->GetXaxis()->SetTitleOffset(1.05);


	h_ntracksInPVwithDimuonM->SetLineWidth(2);
	h_ntracksInPVwithDimuonM->SetLineColor(kGreen);
	h_ntracksInPVwithDimuonM->SetMarkerColor(kGreen);
	h_ntracksInPVwithDimuonM->Draw("same");

	h_ntracksInPVwithDimuon->SetLineWidth(2);
	h_ntracksInPVwithDimuon->SetLineColor(kBlue);
	h_ntracksInPVwithDimuon->SetMarkerColor(kBlue);
	h_ntracksInPVwithDimuon->Draw("same");



	TLegend*leg = new TLegend(0.48, 0.78, 0.9, 0.9, "");
	leg->SetFillColor(kWhite);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.038);

	leg->AddEntry(h_ntracksInEvent, "Number of tracks in event", "l");
	leg->AddEntry(h_ntracksInPVwithDimuonM, "N. of tracks in dimuon PV, v1", "l");
	leg->AddEntry(h_ntracksInPVwithDimuon, "N. of tracks in dimuon PV, v2", "l");

	leg->Draw("same");

	cankres1->SaveAs((TString)"SimplePlots/Ntrack_Comparison_eventPV" + ".root");
	cankres1->SaveAs((TString)"SimplePlots/Ntrack_Comparison_eventPV" + ".pdf");
	cankres1->SaveAs((TString)"SimplePlots/Ntrack_Comparison_eventPV" + ".png");


	TCanvas* cankres2 = new TCanvas("cankres2", "Canvas with results2", 800, 640);
	cankres2->cd();
	cankres2->SetRightMargin(0.14);
	cankres2->SetLogz();
	h_PV_xy->Draw("colz");
	cankres2->SaveAs((TString)"SimplePlots/PV_xy" + ".root");
	cankres2->SaveAs((TString)"SimplePlots/PV_xy" + ".pdf");
	cankres2->SaveAs((TString)"SimplePlots/PV_xy" + ".png");

	TCanvas* cankres3 = new TCanvas("cankres3", "Canvas with results3", 800, 640);
	cankres3->cd();
	cankres3->SetRightMargin(0.05);
	cankres3->SetLogy();
	h_PV_z->SetLineWidth(2);
	h_PV_z->GetYaxis()->SetTitle("Counts");
	h_PV_z->Draw("");
	cankres3->SaveAs((TString)"SimplePlots/PV_z" + ".root");
	cankres3->SaveAs((TString)"SimplePlots/PV_z" + ".pdf");
	cankres3->SaveAs((TString)"SimplePlots/PV_z" + ".png");

	TFile* fout = new TFile("SimplePlots/SimplePlots.root", "RECREATE");
	cankres1->Write();
	cankres2->Write();
	cankres3->Write();
	h_ntracksInEvent->Write();
	h_ntracksInPVwithDimuon->Write();
	h_ntracksInPVwithDimuonM->Write();

	h_PV_z->Write();
	h_PV_xy->Write();

	fout->Close();
}