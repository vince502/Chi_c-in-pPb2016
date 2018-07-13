#ifndef __ChiRootupler_h_
#define __ChiRootupler_h_

  // Declaration of ChiRootupler
  // Description: Saves the muon, dimuon and chi candidate information
  // Implementation:  Ota Kukral based on work of Andre Stahl and Stefano Argiro, Alessandro Degano  and the Torino team, Alberto Sanchez

//O: includes should be cleaned up and moved to cc
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include <TLorentzVector.h>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/PatCandidates/interface/UserData.h> 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
//
// class declaration
//

class ChiRootupler :public edm::EDAnalyzer {
public:
	explicit ChiRootupler(const edm::ParameterSet &);
	~ChiRootupler();
	//static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
	//bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

private:
	
	virtual void analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup);						// ------------ method called for each event  ------------
		
	virtual void beginJob();																				// ------------ method called once each job just before starting event loop  ------------
	virtual void endJob();																					// ------------ method called once each job just after ending the event loop  ------------
	virtual void beginRun(edm::Run const &, edm::EventSetup const &);										// ------------ method called when starting to processes a run  ------------
	virtual void endRun(edm::Run const &, edm::EventSetup const &);											// ------------ method called when ending the processing of a run  ------------
	virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);				// ------------ method called when starting to processes a luminosity block  ------------
	virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);					// ------------ method called when ending the processing of a luminosity block  ------------

	// photon checks
	static bool lt_comparator(std::pair<double, short> a, std::pair<double, short> b); // comparator for checking conversions
	bool Conv_checkTkVtxCompatibility(const reco::Conversion& conv, const reco::VertexCollection&  priVtxs, double sigmaTkVtxComp_, bool& Flag_Best_Out, bool& Flag_SecondBestA_Out, bool& Flag_SecondBestB_Out, double& sigmaMinValue1Out, double& sigmaMinValue2Out);
	bool Conv_foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB);



	// ----------member data ---------------------------
	std::string file_name;

	edm::EDGetTokenT< edm::View <pat::Muon> > muon_label; //is a muon collection
	//edm::EDGetTokenT<reco::MuonCollection>  muon_label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> photon_label;
	edm::EDGetTokenT<reco::ConversionCollection> conversion_label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_label;
	edm::EDGetTokenT<reco::VertexCollection> primaryVertices_label;
	edm::EDGetTokenT<edm::TriggerResults> triggerResults_label;
	edm::EDGetTokenT<reco::GenParticleCollection> genParticles_label;

	bool flag_doMC;
	bool flag_saveExtraThings = false;  //saves all muons, and other not so necessary data
	
	TTree* generalInfo_tree;
	TTree* muon_tree;
	TTree* dimuon_tree;
	TTree* conv_tree;
	TTree* chi_tree;

	//general
	long runNumber;
	long eventNumber;
	long nPrimVertices;
	int muonPerEvent;
	int convPerTriggeredEvent;

	//muon info
	bool muonIsGlobal;
	bool muonIsTracker;
	bool muonIsPF;
	bool muonIsNotGlobalNorTracker;
	bool muonIDHas_TMOneStationTight;
	double muonInnerTrack_dxy;
	double muonInnerTrack_dz;
	int muonTrackerLayersWithMeasurement;
	int muonPixelLayersWithMeasurement;
	bool muonQuality_isHighPurity;
	int muon_charge;
	double muon_eta;
	double muon_pt;
	TLorentzVector muon_p4;

	pat::Muon patMuonStored;
	
	//muon MC
	bool muon_isMatchedMC;
	double muonGen_eta;
	double muonGen_pt;
	TLorentzVector muonGen_p4;
	double muonGen_rDelta;
	double muonGen_ptDelta;
	double muonGen_ptDeltaRel;

	//dimuon info

	//conversion info

	bool convQuality_isHighPurity; 
	bool convQuality_isGeneralTracksOnly;
	double conv_vertexPositionRho; //tbd
	double conv_sigmaTkVtx1;
	double conv_sigmaTkVtx2;
	bool conv_tkVtxCompatibilityOK50;
	bool conv_tkVtxCompatible_bestVertex50;
	bool conv_tkVtxCompatible_secondBestVertexA50;
	bool conv_tkVtxCompatible_secondBestVertexB50;
	bool conv_tkVtxCompatibilityOK3; //test ones
	bool conv_tkVtxCompatible_bestVertex3; //test ones
	bool conv_tkVtxCompatible_secondBestVertexA3; //test
	bool conv_tkVtxCompatible_secondBestVertexB3; //test

	int conv_compatibleInnerHitsOK; //-1: less than 2 tracks, 0: not compatible, 1: yes
	reco::HitPattern conv_hitPat1;
	reco::HitPattern conv_hitPat2;
	bool conv_isCustomHighPurity;//tbd
	double conv_zOfPriVtxFromTracks;
	double conv_dxyPriVtx_Tr1;
	double conv_dxyPriVtx_Tr2;
	double conv_dxyPriVtxTimesCharge_Tr1;
	double conv_dxyPriVtxTimesCharge_Tr2;
	double conv_dxyError_Tr1;
	double conv_dxyError_Tr2;

	int conv_tk1NumOfDOF;
	int conv_tk2NumOfDOF;
	double conv_track1Chi2;
	double conv_track2Chi2;
	double conv_vertexChi2Prob;

	double conv_minDistanceOfApproach;
	TLorentzVector conv_p4;
	double conv_eta;
	double conv_pt;
	
		//convQuality = [''], #O: Changed['highPurity', 'generalTracksOnly']
		//primaryVertexTag = 'offlinePrimaryVertices',
		//convSelection = 'conversionVertex.position.rho>0.0', #O : Changed 1.5
		//wantTkVtxCompatibility = False,
		//sigmaTkVtxComp = 50, #O : Changed 5
		//wantCompatibleInnerHits = True,
		//pfcandidates = 'particleFlow',
		//pi0OnlineSwitch = False,
		//TkMinNumOfDOF = 0, #O : Changed 3
		//wantHighpurity = False,
		//#test
		//vertexChi2ProbCut = 0.0000,
		//trackchi2Cut = 1000,
		//minDistanceOfApproachMinCut = -100.25,
		//minDistanceOfApproachMaxCut = 100.00,
		//#O: Original
		//#vertexChi2ProbCut = 0.0005,
		//#trackchi2Cut = 10,
		//#minDistanceOfApproachMinCut = -0.25,
		//#minDistanceOfApproachMaxCut = 1.00,

	// run and vertex info

	TVector3 primary_v;
	TVector3 secondary_v;
	TVector3 dimuon_v;
	TVector3 conversion_v;

	// Lorentz vectors
	TLorentzVector chi_p4;
	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;
	TLorentzVector photon_p4;

	//Various
	Double_t ele_lowerPt_pt;
	Double_t ele_higherPt_pt;
	Double_t ctpv;
	Double_t ctpv_error;
	Double_t pi0_abs_mass;
	Double_t psi1S_nsigma;
	Double_t psi2S_nsigma;
	Double_t conv_vertex;
	Double_t dz;
	Int_t trigger;


	//MC
	TLorentzVector gen_chic_p4;
	Int_t chic_pdgId;
	TLorentzVector gen_jpsi_p4;
	TLorentzVector gen_photon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;
	edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
	edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

	// static data member definitions
	//
	const double pi0_mass = 0.1349766;
	const Double_t psi1SMass = 3.09691;
	const Double_t psi2SMass = 3.68610;

	/*
	// 2011 par
	static const double Y_sig_par_A = 0.058;
	static const double Y_sig_par_B = 0.047;
	static const double Y_sig_par_C = 0.22;
	*/

	// 2012 par
	const double Y_sig_par_A = 62.62;
	const double Y_sig_par_B = 56.3;
	const double Y_sig_par_C = -20.77;

	// constants for conversion cuts
	const double conv_TkVtxCompSigmaCut = 50.0;

};

#endif 

