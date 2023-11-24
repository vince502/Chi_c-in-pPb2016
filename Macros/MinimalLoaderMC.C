
//Example of tree loading


#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <iostream>
#include "ChiTreeInit.C" //containts the tree structure. I included the .C file, that way you can run this uncompiled (e.g. root -l -b MinimalLoader.C), as well as compiled. 
//If you intend to run it compiled, header can be included instead, and the ChiTreeInit compiled and loaded separately


using namespace std;

void MinimalLoaderMC(const char* fileInMC = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC-Official_v3-bothDir.root", bool isMC = true, const char* fileMCWeight = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Workspace/CMSSW_8_0_30/src/HeavyIonsAnalysis/ChiAnalysis/test/MCWeight_v2.root")
// there is additional file that is used for weighting the MC to match data
{
	//gROOT->ProcessLine(".L ../ChiTreeInit.C+"); //if the header was included, this makes sure it is compiled and loaded

	// load the Ntracks weights (stored as the histogram)
	TH1D* h_weightPVtrk;
	TFile* fMCWeight = new TFile(fileMCWeight, "READ");
	h_weightPVtrk = (TH1D*)fMCWeight->Get("h_weightPVtrk");
	h_weightPVtrk->SetDirectory(0);

	fMCWeight->Close();



	TFile* f1 = new TFile(fileInMC, "READ");
	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");

	if (!event_tree) {
		cout << "Problem with event Tree";
	}
	LoadChiBranches(event_tree, isMC); //in ChiTreeInit.C, can access all the variables now
	
	Long64_t nentries = event_tree->GetEntries();
	if (nentries > 10000) { nentries = 10000; }//Not needed, cut-off for testing

	int counterNotReconstructed = 0; // not needed, just something for the readout

	for (Long64_t i = 0; i < nentries; i++) //event loop
	{
		event_tree->GetEntry(i); //loaded event
		if (i % 1000 == 0) { cout << endl << "event: " << i << " done: " << 100 * i / nentries << "%" << endl << endl; }




		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			continue; //just removes the suspect MC - tested, not many of those, we limit to events with one good simulated decay
		}
		


		/////////////////////////////////////////////////////////////////////////////////////////////
		// this is MC READER, so instead of looping on the RECO collection, we will start with the GEN
		/////////////////////////////////////////////////////////////////////////////////////////////////


		int iChi = ChiPassAllCutsMC(0); //we have only one generated chic (thus index 0), (or remove check above and do this in a for loop) - check if we have a matched reconstructed chic, and passing the cuts. iChi is then the index of the good reco chi candidate
		if (iChi < 0) { //either not matched or didn't pass the cuts
			counterNotReconstructed++;
			continue;
		}

		// further are those chic that were reconstructed + did pass all the cuts

		TLorentzVector* LVchic = (TLorentzVector*)chi_p4->At(iChi); //iChi is now the position of matching reco object, we can read it mass, rap, etc as above.

		int dimuonPos = chi_daughterJpsi_position->at(iChi); // now we are on the reco portion
		int convPos = chi_daughterConv_position->at(iChi); //may not be needed, example of accesing RECO information for daughters
		int muonPos1 = dimuon_muon1_position->at(dimuonPos); //may not be needed, example
		int muonPos2 = dimuon_muon2_position->at(dimuonPos); //may not be needed, example


		//////////////////////
		////   MC WEIGHTING
		//////////////////////

		//The following is an example of code to weight the MC with three weights
		// 1: pPb/Pbp direction (and flip the rapidity for Pbp)
		// 2: pT - to match pT(Jpsi) distributions - WeightForMC_pTpart(pt_JpsiReco) - defined in the helper header ChiTreeInit
		// 3: Ntrack weight - MCweight - loaded from the weight root file

		TLorentzVector* LVJpsiReco = (TLorentzVector*)dimuon_p4->At(dimuonPos);
		double rap_JpsiReco = LVJpsiReco->Rapidity();
		double pt_JpsiReco = dimuon_pt->at(dimuonPos);
		double nTrack_inPV = pvtx_nTracks->at(dimuon_pvtx_indexFromOniaMuMu->at(dimuonPos));

		double MCweight = h_weightPVtrk->GetBinContent(h_weightPVtrk->FindBin(nTrack_inPV)); //here we use the info from the provided 

		if (ispPb == true) { //direction dependent stuff
			MCweight = 1.28*MCweight*WeightForMC_pTpart(pt_JpsiReco);
		}
		else
		{
			MCweight = 0.72*MCweight*WeightForMC_pTpart(pt_JpsiReco);
			rap_JpsiReco = -rap_JpsiReco;
		}

		// Now we could use the MCweight to weight any output histograms, etc.




		///////////////
		// we can continue as in MinimalLoader.C example
		/////////////
		double YourVariable_pT_chi = chi_pt->at(iChi);
		double YourVariable_pT_Jpsi = dimuon_pt->at(dimuonPos);
		bool YourVariable_isSoftMuon1 = muonIsSoft->at(muonPos1);
			
		//////////////
		// we also have access to the GEN variables
		//////////////

		TLorentzVector* LVJpsiGen = (TLorentzVector*)gen_Jpsi_p4->At(0);
		TLorentzVector* LVchicGen = (TLorentzVector*)gen_chic_p4->At(0);

		cout << "We had " << counterNotReconstructed << " chic that were generated but for various reasons not reconstructed (most often photon didn't convert)" << endl<< endl;
		counterNotReconstructed = 0; // just reset the counter
		cout << "This gen chic was reconstructed! Gen mass:  "<< LVchicGen->M() << " Reco mass: " << LVchic->M() << ", pT(Jpsi) gen: " << gen_Jpsi_pt->at(0) << ", and reco: " << pt_JpsiReco << endl;
		cout << "We can access gen and reco 4-vectors, and get angles, e.g. muon 1 angle phi gen: " << ((TLorentzVector*)gen_muon_p4->At(0))->Phi() << " and reco: " << ((TLorentzVector*)muon_p4->At(MuonMCMatched(0)))->Phi() << endl;
		// here I'm using the information that the chic passed all the cuts, so the muons had to be matched. We don't know by default if the first gen muon (index 0) is also the first muon in the reco chic - so I use the MuonMCMatched function for it
		// these functions are described in the ChiTreeInit header, and there are various other useful ones.
		cout << "This gen chic should have weight: " << MCweight << endl << endl;
			
		

	}
}