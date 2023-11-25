
//Example of tree loading


#include "TROOT.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <iostream>
#include "ChiTreeInit.C" //containts the tree structure. I included the .C file, that way you can run this uncompiled (e.g. root -l -b MinimalLoader.C), as well as compiled. 
//If you intend to run it compiled, header can be included instead, and the ChiTreeInit compiled and loaded separately


using namespace std;

void MinimalLoader(const char* fileIn = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV-Nominal_v2-bothDir.root", bool isMC = false)
{
	//gROOT->ProcessLine(".L ../ChiTreeInit.C+"); //if the header was included, this makes sure it is compiled and loaded
	TFile* f1 = new TFile(fileIn, "READ");
	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");

	if (!event_tree) {
		cout << "Problem with event Tree";
	}
	LoadChiBranches(event_tree, isMC); //in ChiTreeInit.C, can access all the variables now
	
	Long64_t nentries = event_tree->GetEntries();
	if (nentries > 50000) { nentries = 50000; }//Not needed, cut-off for testing

	for (Long64_t i = 0; i < nentries; i++) //event loop
	{
		event_tree->GetEntry(i); //loaded event
		if (i % 5000 == 0) { cout << endl << "event: " << i << " done: " << 100 * i / nentries << "%" << endl << endl; }

		for (int iChi = 0; iChi < chi_p4->GetEntriesFast(); iChi++) //chic loop
		{
			int dimuonPos = chi_daughterJpsi_position->at(iChi);
			int convPos = chi_daughterConv_position->at(iChi); //may not be needed, example of accesing information for daughters
			int muonPos1 = dimuon_muon1_position->at(dimuonPos); //may not be needed, example
			int muonPos2 = dimuon_muon2_position->at(dimuonPos); //may not be needed, example

			///////////////
			// you can load your variables, the variable list is in ChiTreeInit.h - your type should match
			/////////////
			double YourVariable_pT_chi = chi_pt->at(iChi);
			double YourVariable_pT_Jpsi = dimuon_pt->at(dimuonPos);
			bool YourVariable_isSoftMuon1 = muonIsSoft->at(muonPos1);
			
			//////////////
			// Alternatively, you can use the functions that are also declared in the ChiTreeInit.h
			//////////////

			if (ChiPassAllCuts(iChi) == true) //select only those that pass all the cuts (including muon, conversion selection/acceptance)
			{
				// Get Lorentz V
				TLorentzVector* LVchic = (TLorentzVector*)chi_p4->At(iChi);
				TLorentzVector* LVdimuon = (TLorentzVector*)dimuon_p4->At(dimuonPos);
				double pT_Jpsi = dimuon_pt->at(dimuonPos);
				double rap_Jpsi = LVdimuon->Rapidity();
				if (runNumber < 285833) rap_Jpsi = -rap_Jpsi; //flip Pbp direction //run number is again declared in ChiTreeInit.h - it is per event variable, not saved as vector
				//also can use (ispPb==false)
				cout << "This is chic candidate passing all the cuts - mass (direct, not PDG subtracted): " << LVchic->M() << ", pT(Jpsi): " << pT_Jpsi << ", and rapidity (J/psi) in p-going direction: " << rap_Jpsi << endl;
				// we normally subtract mumu mass and replace it with PDG value - better resolution. Also, a lot of candidates have too high mass, and would be outside of the fitting window (and are most likely combinatorial bkg)
			}
		}


		////////////
		// if you want to use MC gen information, it is easier to start with the gen object: (needs isMC=true and MC trees) -> I made a new loader for that, as an example. But you could also use this to look at the reco information in the MC
		//////////////
		/*
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			continue; //just removes suspect MC - tested, not many of those, we limit to events with one good simulated decay
		}

		int iChi = ChiPassAllCutsMC(0); //we have only one generated chic (thus index 0), (or remove check above and do this in a for loop)
		if (iChi < 0) continue; //either not matched or didn't pass the cuts
		TLorentzVector* LVchic = (TLorentzVector*)chi_p4->At(iChi); //iChi is now the position of matching reco object, we can read it mass, rap, etc as above.
		*/
	}
}