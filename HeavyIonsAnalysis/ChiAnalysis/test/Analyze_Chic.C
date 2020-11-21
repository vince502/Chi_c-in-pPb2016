
// Macro to analyze the chi_c. Starting with the event tree


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
#include "TClonesArray.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TLorentzVector.h"

using namespace std;

const int PythCode_chic0 = 10441; //Pythia codes
const int PythCode_chic1 = 20443;
const int PythCode_chic2 = 445;

int runNumber, eventNumber;
TLorentzVector* LVchic, *LVdimuon, *LVconv, *LVmuon1, *LVmuon2;
TLorentzVector* LVchic_rot, *LVchic_rotGamma, *LVdimuon_rot, *LVconv_rot, *LVmuon1_rot, *LVmuon2_rot;
TLorentzVector LVaux;

//other
int ntracks_inEvent;
double hfTowerSum_inEvent;
int Trig_Event_HLTDoubleMuOpen;
int nPrimVertices;

//
std::vector <double>* pvtx_z = 0;
std::vector <double>* pvtx_zError = 0;
std::vector <double>* pvtx_x = 0;
std::vector <double>* pvtx_y = 0;
std::vector <double>* pvtx_nTracks = 0;
std::vector <bool>* pvtx_isFake = 0;


//muon info
TClonesArray*  muon_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <bool>* muonIsHLTDoubleMuOpen = 0;
std::vector <bool>* muonIsGlobal = 0;
std::vector <bool>* muonIsTracker = 0;
std::vector <bool>* muonIsPF = 0;
std::vector <bool>* muonIsSoft = 0;
std::vector <bool>* muonIsTight = 0;
std::vector <int>* muonTrackerLayersWithMeasurement = 0;
std::vector <double>* muon_eta = 0;
std::vector <double>* muon_pt = 0;

//dimuon
TClonesArray*  dimuon_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <double>* dimuon_eta = 0;
std::vector <double>* dimuon_pt = 0;
std::vector <double>* dimuon_charge = 0;
std::vector <int>*  dimuon_pvtx_index = 0;
std::vector <double>* dimuon_dz_dimuonvtx_pvtx = 0;
std::vector <double>* dimuon_vtxProb = 0;
std::vector <int>* dimuon_muon1_position = 0; //stores position of first muon in muon collection (no specific order)
std::vector <int>* dimuon_muon2_position = 0; //stores position of second muon in muon collection (no specific order)
std::vector <double>* dimuon_ctpv = 0;
std::vector <double>* dimuon_ctpvError = 0;

// Chi
TClonesArray*  chi_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <double>* chi_eta = 0;
std::vector <double>* chi_pt = 0;
std::vector <int>* chi_daughterJpsi_position = 0; //stores position of daughter Jpsi in dimuon collection
std::vector <int>* chi_daughterConv_position = 0; //stores position of daughter photon (conversion)
std::vector <double>* chi_dzPhotToDimuonVtx = 0; //z distance of photon to dimuon vertex when dxy is minimal
std::vector <double>* chi_dxyPhotToDimuonVtx = 0; //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions


//// Conversions  /////

TClonesArray*  conv_p4 = new TClonesArray("TLorentzVector", 100);
std::vector <bool>* convQuality_isHighPurity = 0;
std::vector <bool>* convQuality_isGeneralTracksOnly = 0;
std::vector <double>* conv_vertexPositionRho = 0;
std::vector <double>* conv_sigmaTkVtx1 = 0;
std::vector <double>* conv_sigmaTkVtx2 = 0;
std::vector <bool>* conv_tkVtxCompatibilityOK = 0;

std::vector <int>* conv_compatibleInnerHitsOK = 0; //-1: less than 2 tracks, 0: not compatible, 1: yes
std::vector <double>* conv_vertexChi2Prob = 0;
std::vector <double>* conv_zOfPriVtx = 0;
std::vector <double>* conv_zOfPriVtxFromTracks = 0;
std::vector <double>* conv_dzToClosestPriVtx = 0;
std::vector <double>* conv_dxyPriVtx_Tr1 = 0;
std::vector <double>* conv_dxyPriVtx_Tr2 = 0;
std::vector <double>* conv_dxyPriVtxTimesCharge_Tr1 = 0;
std::vector <double>* conv_dxyPriVtxTimesCharge_Tr2 = 0;
std::vector <double>* conv_dxyError_Tr1 = 0;
std::vector <double>* conv_dxyError_Tr2 = 0;

std::vector <int>* conv_tk1NumOfDOF = 0;
std::vector <int>* conv_tk2NumOfDOF = 0;
std::vector <double>* conv_track1Chi2 = 0;
std::vector <double>* conv_track2Chi2 = 0;
std::vector <double>* conv_minDistanceOfApproach = 0;
std::vector <double>* conv_eta = 0;
std::vector <double>* conv_pt = 0;

//MC
std::vector <int>* gen_Jpsi_matchPosition = 0;
std::vector <int>* gen_conv_matchPosition = 0;
std::vector <int>* convGen_motherCode;


bool MuonAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;
	if (fabs(eta) < 0.8 && pt < 3.3) return false;
	if (fabs(eta) >= 0.8 && fabs(eta) < 1.5 && pt < 5.81 - 3.14*fabs(eta)) return false;
	if (fabs(eta) >= 1.5 && (pt < 0.8 || pt < 1.89 - 0.526*fabs(eta))) return false;
	return true;
}

bool PhotAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.5) return false;
	if (pt < 0.2) return false;
	return true;
}

bool DimuonSelectionPass(int dimuonPos)  //uses variables loaded in main function
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	return true;
}

bool MuonSelectionPass(int muonPos)  //uses variables loaded in main function
{
	if (muonIsSoft->at(muonPos) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(muonPos) != 1) return false;
	return true;
}


bool PhotSelectionPass(int photPos)  //uses variables loaded in main function
{
	if (convQuality_isHighPurity->at(photPos) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	if (conv_sigmaTkVtx1->at(photPos) > 10) return false;
	if (conv_sigmaTkVtx2->at(photPos) > 10) return false;
	if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (conv_vertexChi2Prob->at(photPos) <= 0.01) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_track1Chi2->at(photPos) >= 10) return false;
	if (conv_track2Chi2->at(photPos) >= 10) return false;
	if (conv_minDistanceOfApproach->at(photPos) <= -0.25) return false;
	if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}


void Analyze_Chic(const char* fileIn = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-bothDirRW2.root", const char* fileOut = "Chi_c_output_RW2_test4.root", const bool doMC = false)
//void Analyze_Chic(const char* fileIn = "Chi_c_pPb8TeV_MCv5.root", const char* fileOut = "Chi_c_output_MCv5_test1.root", const bool doMC = false)
{
	gStyle->SetOptStat(1111);
	//gStyle->SetOptStat(0);

	TH1D* hSignal = new TH1D("hSignal", "", 200, 3, 5);
	TH1D* hSignal_SS = new TH1D("hSignal_SS", "", 200, 3, 5);
	TH1D* hSignal_SB = new TH1D("hSignal_SB", "", 200, 3, 5);
	TH1D* hSignal_rot = new TH1D("hSignal_rot", "", 200, 3, 5);
	TH1D* hSignal_rotGamma = new TH1D("hSignal_rotGamma", "", 200, 3, 5);
	TH1D* hSignal2 = new TH1D("hSignal2", "", 120, 3.2, 3.8);
	TH1D* hSignal3 = new TH1D("hSignal3", "", 120, 3.2, 3.8);

	////Quick pT binning
	//TH1D* hSignalpt1 = new TH1D("hSignalpt1", "", 100, 3, 4);
	//TH1D* hSignalpt2 = new TH1D("hSignalpt2", "", 100, 3, 4);
	//TH1D* hSignalpt3 = new TH1D("hSignalpt3", "", 100, 3, 4);
	//TH1D* hSignalpt4 = new TH1D("hSignalpt4", "", 100, 3, 4);
	//TH1D* hSignalpt5 = new TH1D("hSignalpt5", "", 100, 3, 4);
	//TH1D* hSignalpt6 = new TH1D("hSignalpt6", "", 100, 3, 4);
	//TH1D* hSignalpt7 = new TH1D("hSignalpt7", "", 100, 3, 4);

	//TH2D* hSignal_m_pt = new TH2D("hSignal_m_pt", "", 400, 0, 40, 200, 3, 5);
	//TH2D* hSignal_m_rap = new TH2D("hSignal_m_y", "", 60, -3, 3, 200, 3, 5);
	//TH2D* hSignal_m_phi = new TH2D("hSignal_m_phi", "", 70, -3.5, 3.5, 200, 3, 5);
	//TH2D* hSignal_m_nTr = new TH2D("hSignal_m_nTr", "", 300, 0, 600, 200, 3, 5);
	//TH2D* hSignal_m_Et = new TH2D("hSignal_m_Et", "", 400, 0, 800, 200, 3, 5);


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

//	TH1D* hdimuon_pt = new TH1D("hdimuon_pt", "", 1000, 0, 100);
//	TH1D* hphoton_pt = new TH1D("hphoton_pt", "", 1000, 0, 100);

//	TH2D* hchicCand_pt_rap = new TH2D("hchicCand_pt_rap", "", 60, -3, 3, 200, 0, 20);
//	TH2D* hdimuon_pt_rap = new TH2D("hdimuon_pt_rap", "", 60, -3, 3, 200, 0, 20);
//	TH2D* hphoton_pt_rap = new TH2D("hphoton_pt_rap", "", 60, -3, 3, 100, 0, 10);

	//TH2D* hdimuonAfterCuts_pt_rap = new TH2D("hdimuonAfterCuts_pt_rap", "", 60, -3, 3, 200, 0, 20);
	//TH2D* hmuonPAfterCuts_pt_ps = new TH2D("hmuonPAfterCuts_pt_ps", "", 60, -3, 3, 200, 0, 20);
	//TH2D* hmuonNAfterCuts_pt_ps = new TH2D("hmuonNAfterCuts_pt_ps", "", 60, -3, 3, 200, 0, 20);
	//TH2D* hphotonAfterCuts_pt_rap = new TH2D("hphotonAfterCuts_pt_rap", "", 60, -3, 3, 100, 0, 10);

	//TH2D* hchi_restricted_pt_rap = new TH2D("hchi_restricted_pt_rap", "", 60, -3, 3, 200, 0, 20);
	//TH1D* hchi_restrictedM_pt = new TH1D("hchi_restrictedM_pt", "", 1000, 0, 100);
	//TH1D* hchi_restrictedM_dimuonpt = new TH1D("hchi_restrictedM_dimuonpt", "", 1000, 0, 100);
	//TH1D* hchi_restrictedM_mPpt = new TH1D("hchi_restrictedM_mPpt", "", 1000, 0, 100);
	//TH1D* hchi_restrictedM_mNpt = new TH1D("hchi_restrictedM_mNpt", "", 1000, 0, 100);
	//TH1D* hchi_restrictedM_Photpt = new TH1D("hchi_restrictedM_Photpt", "", 1000, 0, 100);


	// Various
	TH1D* hntracks_inEvent = new TH1D("hntracks_inEvent", "", 400, 0, 400);


	TFile* f1 = new TFile(fileIn, "READ");
	


	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");
	if (!event_tree) {
		cout << "Problem with event Tree";
		//break;
	}

	event_tree->SetBranchAddress("nPrimVertices", &nPrimVertices);
	event_tree->SetBranchAddress("ntracks_inEvent", &ntracks_inEvent);

	event_tree->SetBranchAddress("pvtx_z", &pvtx_z);
	event_tree->SetBranchAddress("pvtx_zError", &pvtx_zError);
	event_tree->SetBranchAddress("pvtx_x", &pvtx_x);
	event_tree->SetBranchAddress("pvtx_y", &pvtx_y);
	event_tree->SetBranchAddress("pvtx_nTracks", &pvtx_nTracks);
	event_tree->SetBranchAddress("pvtx_isFake", &pvtx_isFake);

	event_tree->SetBranchAddress("muon_p4", &muon_p4);
	event_tree->SetBranchAddress("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
	event_tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal);
	event_tree->SetBranchAddress("muonIsTracker", &muonIsTracker);
	event_tree->SetBranchAddress("muonIsPF", &muonIsPF);
	event_tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
	event_tree->SetBranchAddress("muonIsTight", &muonIsTight);
	event_tree->SetBranchAddress("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);
	event_tree->SetBranchAddress("muon_eta", &muon_eta);
	event_tree->SetBranchAddress("muon_pt", &muon_pt);

	event_tree->SetBranchAddress("dimuon_p4", &dimuon_p4);
	event_tree->SetBranchAddress("dimuon_eta", &dimuon_eta);
	event_tree->SetBranchAddress("dimuon_pt", &dimuon_pt);
	event_tree->SetBranchAddress("dimuon_charge", &dimuon_charge);
   // event_tree->SetBranchAddress("dimuon_pvtx_index", &dimuon_pvtx_index);
	event_tree->SetBranchAddress("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);
	event_tree->SetBranchAddress("dimuon_vtxProb", &dimuon_vtxProb);
	event_tree->SetBranchAddress("dimuon_muon1_position", &dimuon_muon1_position);
	event_tree->SetBranchAddress("dimuon_muon2_position", &dimuon_muon2_position);
	event_tree->SetBranchAddress("dimuon_ctpv", &dimuon_ctpv);
	event_tree->SetBranchAddress("dimuon_ctpvError", &dimuon_ctpvError);

	event_tree->SetBranchAddress("chi_p4", &chi_p4);
	event_tree->SetBranchAddress("chi_eta", &chi_eta);
	event_tree->SetBranchAddress("chi_pt", &chi_pt);
	event_tree->SetBranchAddress("chi_daughterJpsi_position", &chi_daughterJpsi_position);
	event_tree->SetBranchAddress("chi_daughterConv_position", &chi_daughterConv_position);
	event_tree->SetBranchAddress("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx);
	event_tree->SetBranchAddress("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx);

	//Conversions
	event_tree->SetBranchAddress("conv_p4", &conv_p4);
	event_tree->SetBranchAddress("convQuality_isHighPurity", &convQuality_isHighPurity);
	event_tree->SetBranchAddress("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly);
	event_tree->SetBranchAddress("conv_vertexPositionRho", &conv_vertexPositionRho);
	event_tree->SetBranchAddress("conv_sigmaTkVtx1", &conv_sigmaTkVtx1);
	event_tree->SetBranchAddress("conv_sigmaTkVtx2", &conv_sigmaTkVtx2);
	event_tree->SetBranchAddress("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOK);

	event_tree->SetBranchAddress("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK);
	event_tree->SetBranchAddress("conv_vertexChi2Prob", &conv_vertexChi2Prob);
	event_tree->SetBranchAddress("conv_zOfPriVtx", &conv_zOfPriVtx);
	event_tree->SetBranchAddress("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks);
	event_tree->SetBranchAddress("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx);
	event_tree->SetBranchAddress("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1);
	event_tree->SetBranchAddress("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2);
	event_tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1);
	event_tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2);
	event_tree->SetBranchAddress("conv_dxyError_Tr1", &conv_dxyError_Tr1);
	event_tree->SetBranchAddress("conv_dxyError_Tr2", &conv_dxyError_Tr2);

	event_tree->SetBranchAddress("conv_tk1NumOfDOF", &conv_tk1NumOfDOF);
	event_tree->SetBranchAddress("conv_tk2NumOfDOF", &conv_tk2NumOfDOF);
	event_tree->SetBranchAddress("conv_track1Chi2", &conv_track1Chi2);
	event_tree->SetBranchAddress("conv_track2Chi2", &conv_track2Chi2);
	event_tree->SetBranchAddress("conv_minDistanceOfApproach", &conv_minDistanceOfApproach);
	event_tree->SetBranchAddress("conv_eta", &conv_eta);
	event_tree->SetBranchAddress("conv_pt", &conv_pt);

	// Gen stuff
	if (doMC) {
		event_tree->SetBranchAddress("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
		event_tree->SetBranchAddress("gen_conv_matchPosition", &gen_conv_matchPosition);
		event_tree->SetBranchAddress("convGen_motherCode", &convGen_motherCode);
	}

	TRandom3* rand = new TRandom3();

	long nchicCounter = 0, nchicCounterPass = 0, muon1ptCounter = 0, muon2ptCounter = 0;
	bool passDimSel = false;

	Long64_t nentries = event_tree->GetEntries();
	cout << nentries << endl;
	//if (nentries > 10000) { nentries = 10000; }
	for (Long64_t i = 0; i < nentries; i++) 
	{

		event_tree->GetEntry(i);
		if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }



		
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

			// check Acceptance and Cuts
			int dimuonPos = chi_daughterJpsi_position->at(iChi);
			//if (DimuonSelectionPass(dimuonPos) == false) continue;
			passDimSel = DimuonSelectionPass(dimuonPos);

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
		    if (PhotSelectionPass(convPos) == false) continue;
			if (passDimSel == true) { ++nchicCounterPass; } // SelectionsPassed
			
			// Get Lorentz V
			LVchic = (TLorentzVector*)chi_p4->At(iChi);
			LVdimuon = (TLorentzVector*)dimuon_p4->At(dimuonPos);
			LVconv = (TLorentzVector*)conv_p4->At(convPos);


			//////// Rotational bkg   /////
			LVmuon1_rot = (TLorentzVector*)muon_p4->At(muon1Pos);
			LVmuon2_rot = (TLorentzVector*)muon_p4->At(muon2Pos);
			LVconv_rot = new TLorentzVector(*(TLorentzVector*)conv_p4->At(convPos));
			LVconv_rot->RotateZ(TMath::Pi());
			//if (muon_pt->at(muon1Pos) > muon_pt->at(muon2Pos)) { muon1ptCounter++; } //crosscheck that they have ordering
			//else muon2ptCounter++;


			if (rand->Uniform(0,1)<0.5) { LVmuon1_rot->RotateZ(TMath::Pi()); }
			else { LVmuon2_rot->RotateZ(TMath::Pi()); }
			LVdimuon_rot = new TLorentzVector (*LVmuon1_rot + *LVmuon2_rot);
			//cout << LVdimuon_rot->Pt() << endl;
			LVchic_rot = new TLorentzVector(*LVdimuon_rot + *LVconv);
			LVchic_rotGamma = new TLorentzVector(*LVdimuon + *LVconv_rot);


			// Obtain yields
			
			// fill dimuon stuff
			
			double dimuonM = LVdimuon->M();
			if (passDimSel == true) {
				hdimuon_M->Fill(dimuonM);
			}
			if (dimuon_charge->at(dimuonPos) != 0) {
				hdimuon_M_SS->Fill(dimuonM);
			}
			if (passDimSel == true) {
				hdimuon_M_rot->Fill(LVdimuon_rot->M());
			}

			// Rotational bkg
			if (LVdimuon_rot->M() >2.95 && LVdimuon_rot->M() <3.2) {
				if (passDimSel == true) {
					hchic_M_rot->Fill(LVchic_rot->M());
					hSignal_rot->Fill(LVchic_rot->M() - LVdimuon_rot->M() + 3.097);
				}
			}
			// Side-band
			if (((dimuonM > 2.5 && dimuonM < 2.75) || (dimuonM > 3.4 && dimuonM < 3.7)) && passDimSel == true) {
				hchic_M_SB->Fill(LVchic->M()); 
				hSignal_SB->Fill(LVchic->M() - LVdimuon->M() + 3.097);
			}

			if (dimuonM<2.95 || dimuonM>3.2) continue; //require narrow dimuon mass
			
			if (passDimSel == true) {
				hchic_M->Fill(LVchic->M()); // just raw M
				hchic_M_rotGamma->Fill(LVchic_rotGamma->M());
			}
			if (dimuon_charge->at(dimuonPos) != 0) {
				hchic_M_SS->Fill(LVchic->M());
			}


			// Assume J/psi mass

			double Mdiff = LVchic->M() - LVdimuon->M() + 3.097;
			if (passDimSel == true) {
				hSignal->Fill(Mdiff); 
				hSignal_rotGamma->Fill(LVchic_rotGamma->M() - LVdimuon->M() + 3.097);
			}
			if (dimuon_charge->at(dimuonPos) != 0) { hSignal_SS->Fill(Mdiff); }



		//}
		//for (int iDimu = 0; iDimu < dimuon_p4->GetEntriesFast(); iDimu++)
		//{
		//	LVdimuon = (TLorentzVector*)dimuon_p4->At(iDimu);

		//	double dimuonM = LVdimuon->M();
		//	if (dimuon_charge->at(iDimu) != 0) continue;
		//	hdimuon_M->Fill(dimuonM);
		//	//cout << dimuonM << endl;

		//cout << "nextEvent" << endl;




					//hmuon_M->Fill(dimuonM);
					////hmuon_Mzoom->Fill(dimuonM);
					//hphoton_M->Fill(LVphoton->M());

					//hdimuon_pt->Fill(LVdimuon->Pt());
					//hphoton_pt->Fill(LVphoton->Pt());

					//hchicCand_pt_rap->Fill(LVchic->Rapidity(), LVchic->Pt());
					//hdimuon_pt_rap->Fill(LVdimuon->Rapidity(), LVdimuon->Pt());
					//hmuonP_pt_ps->Fill(LVmuonP->PseudoRapidity(), LVmuonP->Pt());
					//hmuonN_pt_ps->Fill(LVmuonN->PseudoRapidity(), LVmuonN->Pt());
					//hphoton_pt_rap->Fill(LVphoton->Rapidity(), LVphoton->Pt());


					////cut section
					//if (dimuonM > 3.3) { continue; }
					//if (dimuonM < 2.9) { continue; }

					//if (LVmuonN->Pt() < 2) { continue; } //2
					//if (LVmuonP->Pt() < 2) { continue; } //2
					//if (dimuonVertexP < 0.05) { continue; }

					////extra cuts

					////if (TMath::Abs(V3primaryV->Z() - V3secondaryV->Z()) > 0.4) { continue; }
					////if (TMath::Abs(LVphoton->Pt())<0.5) { continue; }

					////if (TMath::Sqrt((LVphoton->Eta() - LVchic->Rapidity())*(LVphoton->Eta() - LVchic->Rapidity()) + ((LVphoton->Phi() - LVchic->Phi())*(LVphoton->Phi() - LVchic->Phi()))) < 0.8) { continue; };

					//hmuon_Mzoom->Fill(dimuonM);
					//Mdiff = LVchic->M() - LVdimuon->M() + 3.097;

					//double dchicPt = LVchic->Pt();
					//if (dchicPt > 6.0 && dchicPt < 40.0) { hSignal->Fill(Mdiff); }

					//if (dchicPt < 6.0) { hSignalpt1->Fill(Mdiff); }
					//else if (dchicPt < 9.0) { hSignalpt2->Fill(Mdiff); }
					//else if (dchicPt < 12.0) { hSignalpt3->Fill(Mdiff); }
					//else if (dchicPt < 15.0) { hSignalpt4->Fill(Mdiff); }
					//else if (dchicPt < 20.0) { hSignalpt5->Fill(Mdiff); }
					//else if (dchicPt < 40.0) { hSignalpt6->Fill(Mdiff); }
					//else { hSignalpt7->Fill(Mdiff); }

					//if (dchicPt > 6.0 && dchicPt < 40.0) {
					//	hSignal_m_pt->Fill(dchicPt, Mdiff);
					//	hSignal_m_rap->Fill(LVchic->Rapidity(), Mdiff);
					//	hSignal_m_phi->Fill(LVchic->Phi(), Mdiff);
					//	hSignal_m_nTr->Fill(nTrack, Mdiff);
					//	hSignal_m_Et->Fill(EHF_trans, Mdiff);
					//}

					//hchic_M->Fill(LVchic->M());
					////if (TMath::Abs(LVphoton->Eta())<1.3){ hSignal2->Fill(Mdiff); }
					////if (TMath::Abs(LVphoton->Eta())<1.2) { hSignal3->Fill(Mdiff); }

					//if (Mdiff > 3.5 &&Mdiff < 3.52)
					//{
					//	hchi_restricted_pt_rap->Fill(LVchic->Rapidity(), dchicPt);
					//	hchi_restrictedM_pt->Fill(LVchic->Pt());
					//	hchi_restrictedM_dimuonpt->Fill(LVdimuon->Pt());
					//	hchi_restrictedM_mPpt->Fill(LVmuonP->Pt());
					//	hchi_restrictedM_mNpt->Fill(LVmuonN->Pt());
					//	hchi_restrictedM_Photpt->Fill(LVphoton->Pt());
					//	//cout << "Event number: " << eventNumber << endl;
					//}

					//hdimuonAfterCuts_pt_rap->Fill(LVdimuon->Rapidity(), LVdimuon->Pt());
					//hmuonPAfterCuts_pt_ps->Fill(LVmuonP->PseudoRapidity(), LVmuonP->Pt());
					//hmuonNAfterCuts_pt_ps->Fill(LVmuonN->PseudoRapidity(), LVmuonN->Pt());
					//hphotonAfterCuts_pt_rap->Fill(LVphoton->Rapidity(), LVphoton->Pt());


		} // end of chic loop

	} //end of event loop

	cout << "Processed " << nchicCounter << " chic candidates, " << nchicCounterPass << " passed the selections, that is " << ((double)nchicCounterPass/(double)nchicCounter*100) << " percent" << endl;
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

	//hSignalpt1->Write();
	//hSignalpt2->Write();
	//hSignalpt3->Write();
	//hSignalpt4->Write();
	//hSignalpt5->Write();
	//hSignalpt6->Write();
	//hSignalpt7->Write();

	//hSignal_m_pt->Write();
	//hSignal_m_rap->Write();
	//hSignal_m_phi->Write();
	//hSignal_m_nTr->Write();
	//hSignal_m_Et->Write();

	//hmuon_M->Write();
	//hmuon_Mzoom->Write();
	//hphoton_M->Write();

	//hdimuon_pt->Write();
	//hphoton_pt->Write();

	//hchicCand_pt_rap->Write();
	//hdimuon_pt_rap->Write();
	//hmuonP_pt_ps->Write();
	//hmuonN_pt_ps->Write();
	//hphoton_pt_rap->Write();

	//hdimuonAfterCuts_pt_rap->Write();
	//hmuonPAfterCuts_pt_ps->Write();
	//hmuonNAfterCuts_pt_ps->Write();
	//hphotonAfterCuts_pt_rap->Write();

	//hchi_restricted_pt_rap->Write();
	//hchi_restrictedM_pt->Write();
	//hchi_restrictedM_dimuonpt->Write();
	//hchi_restrictedM_mPpt->Write();
	//hchi_restrictedM_mNpt->Write();
	//hchi_restrictedM_Photpt->Write();


	hntracks_inEvent->Write();
	
}

