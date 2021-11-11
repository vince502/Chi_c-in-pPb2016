
// Macro to create the boosted decision trees to obtain conversion cuts for the chi_c. Starting with the event tree


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
#include "TF1.h"
#include "TClonesArray.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../ChiTreeInit.C"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"


using namespace std;
using namespace RooFit;

const int PythCode_chic0 = 10441; //Pythia codes
const int PythCode_chic1 = 20443;
const int PythCode_chic2 = 445;

const double k_mass_c0 = 3.4148; //pdg 1. 2018
const double k_mass_c1 = 3.5107;
const double k_mass_c2 = 3.5562;


double bins_pT[] = {6, 9, 12, 18, 30};
int  nbins_pT = sizeof(bins_pT) / sizeof(double) - 1;

double bins_y[] = {-2.4, -1.6, -1.0, 0, 1.0, 1.6, 2.4 };
int  nbins_y = sizeof(bins_y) / sizeof(double) - 1;

double bins_nTrk[] = { 0, 50, 100, 150, 200, 300, 400 };
int  nbins_nTrk = sizeof(bins_nTrk) / sizeof(double) - 1;


TLorentzVector* LVchic, *LVdimuon, *LVconv, *LVmuon1, *LVmuon2;
TLorentzVector* LVchic_rot, *LVchic_rotGamma, *LVdimuon_rot, *LVconv_rot, *LVmuon1_rot, *LVmuon2_rot;
TLorentzVector LVaux;


	//conversion info

bool convQuality_isHighPurityOut;
bool convQuality_isGeneralTracksOnlyOut;

double conv_vertexPositionRhoOut;
double conv_sigmaTkVtx1Out;
double conv_sigmaTkVtx2Out;
bool conv_tkVtxCompatibilityOKOut;
bool conv_tkVtxCompatible_bestVertexOut;

int conv_compatibleInnerHitsOKOut; //-1: less than 2 tracks, 0: not compatible, 1: yes

double conv_vertexChi2ProbOut;
int  conv_pvtx_indexOut;
double conv_zOfPriVtxOut; // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
double conv_zOfPriVtxFromTracksOut;
double conv_dzToClosestPriVtxOut;
double conv_dxyPriVtx_Tr1Out;
double conv_dxyPriVtx_Tr2Out;
double conv_dxyPriVtxTimesCharge_Tr1Out;
double conv_dxyPriVtxTimesCharge_Tr2Out;
double conv_dxyError_Tr1Out;
double conv_dxyError_Tr2Out;

int conv_tk1NumOfDOFOut;
int conv_tk2NumOfDOFOut;
double conv_track1Chi2Out;
double conv_track2Chi2Out;
double conv_Tr1_ptOut;
double conv_Tr2_ptOut;

double conv_minDistanceOfApproachOut;
TClonesArray*  conv_p4Out; //TLorentzVector
double conv_etaOut;
double conv_ptOut;

double chic_mDiff, chic_mJpsi;
ofstream file_cutWorkingPoints;

void GenerateSample(bool isSignal, const char* fileIn, const char* fileOut, bool isMC=false)
{
	TFile* f1 = new TFile(fileIn, "READ");

	TTree* event_tree = (TTree*)f1->Get("ChiRootuple/event_tree");
	if (!event_tree) {
		cout << "Problem with event Tree";
		//break;
	}
	LoadChiBranches(event_tree, isMC);

	//counters
	int weird_decay_counter = 0;
	long nchicCounter = 0, nchicCounterPass = 0, nchicCounterPassMass = 0;
	long signalEff_den = 0, signalEff_numLoose = 0, signalEff_numTight = 0;
	long bkgEff_den = 0, bkgEff_numLoose = 0, bkgEff_numTight = 0;

	file_cutWorkingPoints.open("Log_cutBasedWorkingPoints.txt", ios::app);
	file_cutWorkingPoints << endl << endl << "N E W   R U N " << endl << endl;


	bool passDimSel = false;
	bool passDimSelTight = false;

	//create the output tree

	TFile* f2 = new TFile(fileOut, "RECREATE");
	TTree* event_treeOut = new TTree("event_tree", "Tree with events");
	// conversions
	event_treeOut->Branch("convQuality_isHighPurity", &convQuality_isHighPurityOut);
	event_treeOut->Branch("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnlyOut);
	event_treeOut->Branch("conv_vertexPositionRho", &conv_vertexPositionRhoOut);
	event_treeOut->Branch("conv_sigmaTkVtx1", &conv_sigmaTkVtx1Out);
	event_treeOut->Branch("conv_sigmaTkVtx2", &conv_sigmaTkVtx2Out);
	event_treeOut->Branch("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOKOut);
	event_treeOut->Branch("conv_tkVtxCompatible_bestVertex", &conv_tkVtxCompatible_bestVertexOut);
	event_treeOut->Branch("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOKOut);
	event_treeOut->Branch("conv_vertexChi2Prob", &conv_vertexChi2ProbOut);
	event_treeOut->Branch("conv_pvtx_index", &conv_pvtx_indexOut);
	event_treeOut->Branch("conv_zOfPriVtx", &conv_zOfPriVtxOut);
	event_treeOut->Branch("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracksOut);
	event_treeOut->Branch("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtxOut);
	event_treeOut->Branch("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1Out);
	event_treeOut->Branch("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2Out);
	event_treeOut->Branch("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1Out);
	event_treeOut->Branch("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2Out);
	event_treeOut->Branch("conv_dxyError_Tr1", &conv_dxyError_Tr1Out);
	event_treeOut->Branch("conv_dxyError_Tr2", &conv_dxyError_Tr2Out);

	event_treeOut->Branch("conv_tk1NumOfDOF", &conv_tk1NumOfDOFOut);
	event_treeOut->Branch("conv_tk2NumOfDOF", &conv_tk2NumOfDOFOut);
	event_treeOut->Branch("conv_track1Chi2", &conv_track1Chi2Out);
	event_treeOut->Branch("conv_track2Chi2", &conv_track2Chi2Out);
	event_treeOut->Branch("conv_Tr1_pt", &conv_Tr1_ptOut);
	event_treeOut->Branch("conv_Tr2_pt", &conv_Tr2_ptOut);

	event_treeOut->Branch("conv_minDistanceOfApproach", &conv_minDistanceOfApproachOut);

	//event_treeOut->Branch("conv_p4", "TClonesArray", &conv_p4, 32000, 0Out);
	event_treeOut->Branch("conv_eta", &conv_etaOut);
	event_treeOut->Branch("conv_pt", &conv_ptOut);
	event_treeOut->Branch("chic_mJpsi", &chic_mJpsi);
	event_treeOut->Branch("chic_m", &chic_mDiff);


	Long64_t nentries = event_tree->GetEntries();
	//if (nentries > 200000) { nentries = 200000; }
	cout << nentries << endl;
	for (Long64_t i = 0; i < nentries; i++)
	{
		event_tree->GetEntry(i);
		if (i % 10000 == 0) { cout << "Processing " << (isSignal? "signal ":"background ")<< "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }
		if (isMC == true) {
			if (gen_pdgId->at(0) != PythCode_chic1 && gen_pdgId->at(0) != PythCode_chic2) { continue; } //remove chic0 from MC, expecting one gen chic per event
			if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
			{
				weird_decay_counter++;
				continue;
			}
		}

		for (int iChi = 0; iChi < chi_p4->GetEntriesFast(); iChi++) {
			if (isMC == true && isSignal == true) { // if MC signal, use only chic that is matched to the gen, otherwise use all the candidates
				if (ChiIsMatchedAllDaughters(iChi, 0) == false) continue;
			}

			++nchicCounter;
			// check Acceptance and Cuts
			int dimuonPos = chi_daughterJpsi_position->at(iChi);
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
			//if (PhotSelectionPass(convPos) == false) continue;
			if (passDimSel == true) { ++nchicCounterPass; } // SelectionsPassed
			// Get Lorentz V
			LVchic = (TLorentzVector*)chi_p4->At(iChi);
			LVdimuon = (TLorentzVector*)dimuon_p4->At(dimuonPos);
			LVconv = (TLorentzVector*)conv_p4->At(convPos);

			double pT_chi = chi_pt->at(iChi);
			double rap_chi = LVchic->Rapidity();
			double m_chi = LVchic->M();
			double dimuonM = LVdimuon->M();
			double Mdiff = m_chi - dimuonM + 3.097;// Assume J/psi mass

			if (passDimSel == true && dimuonM > 2.95 && dimuonM < 3.2) ++nchicCounterPassMass;

			// Obtain yields

			// fill dimuon stuff



			//// Rotational bkg
			//if (LVdimuon_rot->M() > 2.95 && LVdimuon_rot->M() < 3.2) {
			//	if (passDimSel == true) {
			//		hchic_M_rot->Fill(LVchic_rot->M());
			//		hSignal_rot->Fill(LVchic_rot->M() - LVdimuon_rot->M() + 3.097);
			//	}
			//}
			//// Side-band
			//if (((dimuonM > 2.5 && dimuonM < 2.75) || (dimuonM > 3.4 && dimuonM < 3.7)) && passDimSel == true) {
			//	hchic_M_SB->Fill(m_chi);
			//	hSignal_SB->Fill(m_chi - dimuonM + 3.097);
			//}


			if (dimuonM<2.95 || dimuonM>3.2) continue; //require narrow dimuon mass


			if (isSignal==true && passDimSel == true) {
				//signal

				convQuality_isHighPurityOut = convQuality_isHighPurity->at(convPos);
				convQuality_isGeneralTracksOnlyOut = convQuality_isGeneralTracksOnly->at(convPos);

				conv_vertexPositionRhoOut = conv_vertexPositionRho->at(convPos);
				conv_sigmaTkVtx1Out = conv_sigmaTkVtx1->at(convPos);
				conv_sigmaTkVtx2Out = conv_sigmaTkVtx2->at(convPos);
				conv_tkVtxCompatibilityOKOut = conv_tkVtxCompatibilityOK->at(convPos);
				conv_tkVtxCompatible_bestVertexOut = conv_tkVtxCompatible_bestVertex->at(convPos);


				conv_compatibleInnerHitsOKOut = conv_compatibleInnerHitsOK->at(convPos); //-1: less than 2 tracks, 0: not compatible, 1: yes

				conv_vertexChi2ProbOut = conv_vertexChi2Prob->at(convPos);
				conv_zOfPriVtxOut = conv_zOfPriVtx->at(convPos); // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
				conv_zOfPriVtxFromTracksOut = conv_zOfPriVtxFromTracks->at(convPos);
				conv_dzToClosestPriVtxOut = conv_dzToClosestPriVtx->at(convPos);
				conv_dxyPriVtx_Tr1Out = conv_dxyPriVtx_Tr1->at(convPos);
				conv_dxyPriVtx_Tr2Out = conv_dxyPriVtx_Tr2->at(convPos);
				conv_dxyPriVtxTimesCharge_Tr1Out = conv_dxyPriVtxTimesCharge_Tr1->at(convPos);
				conv_dxyPriVtxTimesCharge_Tr2Out = conv_dxyPriVtxTimesCharge_Tr2->at(convPos);
				conv_dxyError_Tr1Out = conv_dxyError_Tr1->at(convPos);
				conv_dxyError_Tr2Out = conv_dxyError_Tr2->at(convPos);

				conv_tk1NumOfDOFOut = conv_tk1NumOfDOF->at(convPos);
				conv_tk2NumOfDOFOut = conv_tk2NumOfDOF->at(convPos);
				conv_track1Chi2Out = conv_track1Chi2->at(convPos);
				conv_track2Chi2Out = conv_track2Chi2->at(convPos);
				conv_Tr1_ptOut = conv_Tr1_pt->at(convPos);
				conv_Tr2_ptOut = conv_Tr2_pt->at(convPos);

				conv_minDistanceOfApproachOut = conv_minDistanceOfApproach->at(convPos);
				conv_etaOut = conv_eta->at(convPos);
				conv_ptOut = conv_pt->at(convPos);
				chic_mDiff = Mdiff;
				chic_mJpsi = dimuonM;

				event_treeOut->Fill();
				signalEff_den++;
				if (PhotSelectionPass(convPos) == true) {
					signalEff_numLoose++;
				}
				if (PhotSelectionPassTight(convPos) == true) {
					signalEff_numTight++;
				}
				
			}

			if (isSignal == false && dimuon_charge->at(dimuonPos) != 0 && dimuon_pt->at(dimuonPos)>6) {
				//bkg
				convQuality_isHighPurityOut = convQuality_isHighPurity->at(convPos);
				convQuality_isGeneralTracksOnlyOut = convQuality_isGeneralTracksOnly->at(convPos);

				conv_vertexPositionRhoOut = conv_vertexPositionRho->at(convPos);
				conv_sigmaTkVtx1Out = conv_sigmaTkVtx1->at(convPos);
				conv_sigmaTkVtx2Out = conv_sigmaTkVtx2->at(convPos);
				conv_tkVtxCompatibilityOKOut = conv_tkVtxCompatibilityOK->at(convPos);
				conv_tkVtxCompatible_bestVertexOut = conv_tkVtxCompatible_bestVertex->at(convPos);


				conv_compatibleInnerHitsOKOut = conv_compatibleInnerHitsOK->at(convPos); //-1: less than 2 tracks, 0: not compatible, 1: yes

				conv_vertexChi2ProbOut = conv_vertexChi2Prob->at(convPos);
				conv_zOfPriVtxOut = conv_zOfPriVtx->at(convPos); // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
				conv_zOfPriVtxFromTracksOut = conv_zOfPriVtxFromTracks->at(convPos);
				conv_dzToClosestPriVtxOut = conv_dzToClosestPriVtx->at(convPos);
				conv_dxyPriVtx_Tr1Out = conv_dxyPriVtx_Tr1->at(convPos);
				conv_dxyPriVtx_Tr2Out = conv_dxyPriVtx_Tr2->at(convPos);
				conv_dxyPriVtxTimesCharge_Tr1Out = conv_dxyPriVtxTimesCharge_Tr1->at(convPos);
				conv_dxyPriVtxTimesCharge_Tr2Out = conv_dxyPriVtxTimesCharge_Tr2->at(convPos);
				conv_dxyError_Tr1Out = conv_dxyError_Tr1->at(convPos);
				conv_dxyError_Tr2Out = conv_dxyError_Tr2->at(convPos);

				conv_tk1NumOfDOFOut = conv_tk1NumOfDOF->at(convPos);
				conv_tk2NumOfDOFOut = conv_tk2NumOfDOF->at(convPos);
				conv_track1Chi2Out = conv_track1Chi2->at(convPos);
				conv_track2Chi2Out = conv_track2Chi2->at(convPos);
				conv_Tr1_ptOut = conv_Tr1_pt->at(convPos);
				conv_Tr2_ptOut = conv_Tr2_pt->at(convPos);

				conv_minDistanceOfApproachOut = conv_minDistanceOfApproach->at(convPos);
				conv_etaOut = conv_eta->at(convPos);
				conv_ptOut = conv_pt->at(convPos);
				chic_mDiff = Mdiff;
				chic_mJpsi = dimuonM;

				event_treeOut->Fill();

				bkgEff_den++;
				if (PhotSelectionPass(convPos) == true) {
					bkgEff_numLoose++;
				}
				if (PhotSelectionPassTight(convPos) == true) {
					bkgEff_numTight++;
				}

				
			}

		}


	} //end of event loop

	cout << "Processed " << nchicCounter << " chic candidates, " << nchicCounterPass << " passed the selections, that is " << ((double)nchicCounterPass / (double)nchicCounter * 100) << " percent" << endl;
	cout << "Processed " << nchicCounter << " chic candidates, " << nchicCounterPassMass << " passed the selections and mass window requirement, that is " << ((double)nchicCounterPassMass / (double)nchicCounter * 100) << " percent" << endl;

	cout << endl << "Cut based selections: " << endl;
	if (isSignal) {
		cout << "Loose selection - signal eff: " << (double)signalEff_numLoose / (double)signalEff_den * 100 << "% (" << signalEff_numLoose << ")/(" << signalEff_den << ")" << endl; 
		cout << "Tight selection - signal eff: " << (double)signalEff_numTight / (double)signalEff_den * 100 << "% (" << signalEff_numTight << ")/(" << signalEff_den << ")" << endl;
		file_cutWorkingPoints << "Loose selection - signal eff: " << (double)signalEff_numLoose / (double)signalEff_den * 100 << "% (" << signalEff_numLoose << ")/(" << signalEff_den << ")" << endl;
		file_cutWorkingPoints << "Tight selection - signal eff: " << (double)signalEff_numTight / (double)signalEff_den * 100 << "% (" << signalEff_numTight << ")/(" << signalEff_den << ")" << endl;
	}
	else {
		cout << "Loose selection - bkg rejection: " << (1.0 - (double)bkgEff_numLoose / (double)bkgEff_den) * 100 << "% (1-(" << bkgEff_numLoose << ")/(" << bkgEff_den << "))" << endl;
		cout << "Tight selection - bkg rejection: " << (1.0 - (double)bkgEff_numTight / (double)bkgEff_den) * 100 << "% (1-(" << bkgEff_numTight << ")/(" << bkgEff_den << "))" << endl;
		file_cutWorkingPoints << "Loose selection - bkg rejection: " << (1.0 - (double)bkgEff_numLoose / (double)bkgEff_den) * 100 << "% (1-(" << bkgEff_numLoose << ")/(" << bkgEff_den << "))" << endl;
		file_cutWorkingPoints << "Tight selection - bkg rejection: " << (1.0 - (double)bkgEff_numTight / (double)bkgEff_den) * 100 << "% (1-(" << bkgEff_numTight << ")/(" << bkgEff_den << "))" << endl;
	}


	f2->Write();
	f2->Close();
	file_cutWorkingPoints.close();

}




///////////////////////////////////////

/// P R O G R A M   S T A R T   ///////

////////////////////////////////////


void BDT_Chic(bool flagGenerateSBsamples = false, const char* fileSig = "BDT_Signal_test.root", const char* fileBkg = "BDT_Background_test.root", const char* fileInSigTree = "/afs/cern.ch/work/o/okukral/ChicData/Chi_c_pPb8TeV-MC8_BothDir.root", const char* fileInBkgTree = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV-bothDirRW3.root", const char* fileOut = "BDT_OutputWeight.root")
{
	//gStyle->SetOptStat(1111);
	//gStyle->SetOptStat(0);
	setTDRStyle();


	if (flagGenerateSBsamples)
	{
		GenerateSample(true, fileInSigTree, fileSig, true); //signal, fileIn, fileOut, isMC
		GenerateSample(false, fileInBkgTree, fileBkg, false);
	}

	// Load the trees

	TFile* fSig = new TFile(fileSig, "READ");
	TTree* event_treeSig = (TTree*)fSig->Get("event_tree");
	
	TFile* fBkg = new TFile(fileBkg, "READ");
	TTree* event_treeBkg = (TTree*)fBkg->Get("event_tree");



	

	////////////////////////////////////////  
	///////   M V A   //////////                               // based on TMVAClassification.C tutorial
	/////////////////////////////////////////


   // This loads the library
	TMVA::Tools::Instance();


	// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
	TFile* outputFile = TFile::Open(fileOut, "RECREATE");

	// The first argument is the base of the name of all the
// weightfiles in the directory weight/
//
// The second argument is the output file for the training results
// All TMVA output can be suppressed by removing the "!" (not) in
// front of the "Silent" argument in the option string
	TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,
		"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

	TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

	// Define the input variables that shall be used for the MVA training
// note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
// [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
	
	//dataloader->AddVariable("myvar2 := var1-var2", "Expression 2", "", 'F');
	//dataloader->AddVariable("var3", "Variable 3", "units", 'F');
	//dataloader->AddVariable("var4", "Variable 4", "units", 'F');

	dataloader->AddVariable("convQuality_isHighPurity", 'I');
	dataloader->AddVariable("convQuality_isGeneralTracksOnly", 'I');
	dataloader->AddSpectator("conv_vertexPositionRho", 'F');
	dataloader->AddVariable("conv_sigmaTkVtx1", 'F');
	dataloader->AddVariable("conv_sigmaTkVtx2", 'F');
	dataloader->AddVariable("conv_tkVtxCompatibilityOK", 'I');
	dataloader->AddVariable("conv_compatibleInnerHitsOK", 'I');
	dataloader->AddVariable("conv_vertexChi2Prob", 'F');
	//dataloader->AddSpectator("conv_zOfPriVtx", 'F');
	//dataloader->AddSpectator("conv_zOfPriVtxFromTracks", 'F');
	dataloader->AddVariable("conv_dzToClosestPriVtx", 'F');
	dataloader->AddVariable("conv_dxyPriVtxTimesCharge_Tr1", 'F');
	dataloader->AddVariable("conv_dxyPriVtxTimesCharge_Tr2", 'F');
	//dataloader->AddSpectator("conv_dxyError_Tr1", 'F');
	//dataloader->AddSpectator("conv_dxyError_Tr2", 'F');
	dataloader->AddVariable("conv_tk1NumOfDOF", 'I');
	dataloader->AddVariable("conv_tk2NumOfDOF", 'I');
	dataloader->AddVariable("conv_track1Chi2", 'F');
	dataloader->AddVariable("conv_track2Chi2", 'F');
	//dataloader->AddSpectator("conv_Tr1_pt", 'F');
	//dataloader->AddSpectator("conv_Tr2_pt", 'F');

	dataloader->AddVariable("conv_minDistanceOfApproach", 'F');
	dataloader->AddSpectator("conv_eta", 'F');
	dataloader->AddSpectator("conv_pt", 'F');

	//bool convQuality_isHighPurity;
	//bool convQuality_isGeneralTracksOnly;

	//double conv_vertexPositionRho;
	//double conv_sigmaTkVtx1;
	//double conv_sigmaTkVtx2;
	//bool conv_tkVtxCompatibilityOK;
	//bool conv_tkVtxCompatible_bestVertex;

	//int conv_compatibleInnerHitsOK; //-1: less than 2 tracks, 0: not compatible, 1: yes

	//double conv_vertexChi2Prob;
	//int  conv_pvtx_index;
	//double conv_zOfPriVtx; // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
	//double conv_zOfPriVtxFromTracks;
	//double conv_dzToClosestPriVtx;
	//double conv_dxyPriVtx_Tr1;
	//double conv_dxyPriVtx_Tr2;
	//double conv_dxyPriVtxTimesCharge_Tr1;
	//double conv_dxyPriVtxTimesCharge_Tr2;
	//double conv_dxyError_Tr1;
	//double conv_dxyError_Tr2;

	//int conv_tk1NumOfDOF;
	//int conv_tk2NumOfDOF;
	//double conv_track1Chi2;
	//double conv_track2Chi2;
	//double conv_Tr1_pt;
	//double conv_Tr2_pt;

	//double conv_minDistanceOfApproach;
	//TClonesArray*  conv_p4; //TLorentzVector
	//double conv_eta;
	//double conv_pt;

	//double chic_mDiff, chic_mJpsi;








	// You can add so-called "Spectator variables", which are not used in the MVA training,
	// but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
	// input variables, the response values of all trained MVAs, and the spectator variables

	//dataloader->AddSpectator("spec1 := var1*2", "Spectator 1", "units", 'F');
	//dataloader->AddSpectator("spec2 := var1*3", "Spectator 2", "units", 'F');


	// global event weights per tree (see below for setting event-wise weights)
	double signalWeight = 1.0;
	double backgroundWeight = 1.0;

	dataloader->AddSignalTree(event_treeSig, signalWeight);
	dataloader->AddBackgroundTree(event_treeBkg, backgroundWeight);

	dataloader->PrepareTrainingAndTestTree("", "SplitMode=random:!V");

	factory->BookMethod(dataloader, TMVA::Types::kCuts, "Cuts",
		"!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");

	//// Likelihood ("naive Bayes estimator")
	//factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "Likelihood",
	//	"H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");


	// Fisher discriminant (same as LD)
	factory->BookMethod(dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");


	// TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons

	factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator");


	//// Gradient Boost
	//factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
	//	"!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

	// Adaptive Boost // default
	factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
		"!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

	factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT2",
		"!H:!V:NTrees=100:MinNodeSize=0.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

	factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT3",
		"!H:!V:NTrees=300:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

	// Decorrelation + Adaptive Boost
	factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTD",
		"!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");


	factory->TrainAllMethods();

	// Evaluate all MVAs using the set of test events
	factory->TestAllMethods();

	// Evaluate and compare performance of all configured MVAs
	factory->EvaluateAllMethods();

	// --------------------------------------------------------------

	// Save the output
	outputFile->Close();

	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;

	delete factory;
	delete dataloader;
	// Launch the GUI for the root macros
	if (!gROOT->IsBatch()) TMVA::TMVAGui(fileOut);


}







