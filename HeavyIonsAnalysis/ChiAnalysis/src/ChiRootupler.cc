// -*- C++ -*-
//
// Package:    ChiRootupler
// Class:      ChiRootupler
// 
/**
 Description: Saves the muon, dimuon and chi candidate information

 Implementation:  Ota Kukral based on work of Andre Stahl and Stefano Argiro, Alessandro Degano  and the Torino team, Alberto Sanchez
*/

// system include files
#include <memory>

// user include files
#include <HeavyIonsAnalysis/ChiAnalysis/interface/ChiRootupler.h>
//#include <HeavyIonsAnalysis/ChiAnalysis/interface/mydict.cxx>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/PatCandidates/interface/UserData.h> 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for kinematic refit
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"



#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
//#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include <boost/foreach.hpp>


//constructor
ChiRootupler::ChiRootupler(const edm::ParameterSet & iConfig) :
	//muon_label(consumes<edm::View <pat::Muon> >(iConfig.getParameter < edm::InputTag >("muon_cand"))),
	muon_label(consumes< pat::MuonCollection >(iConfig.getParameter < edm::InputTag >("muon_cand"))),
	dimuon_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("dimuon_cand"))),
	//photon_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("photon_cand"))),
	conversion_label(consumes<reco::ConversionCollection>(iConfig.getParameter < edm::InputTag >("conversions_ch"))),
	primaryVertices_label(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag >("primaryVertices"))),
	triggerResults_label(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag >("TriggerResults"))),
	centrality_label(consumes<reco::Centrality>(iConfig.getParameter < edm::InputTag >("centralityInfo"))),
	genParticles_label(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag >("genParticlesTag"))),
	flag_doMC(iConfig.getParameter < bool >("isMC"))
{

	// initialize arrays
	muon_p4 = new TClonesArray("TLorentzVector", 1000);
	muonGen_p4 = new TClonesArray("TLorentzVector", 1000);

	dimuon_p4 = new TClonesArray("TLorentzVector", 100);
	dimuon_vtx = new TClonesArray("TVector3", 100);
	
	conv_vtx = new TClonesArray("TVector3", 1000);
	conv_p4 = new TClonesArray("TLorentzVector", 1000);
	convGen_p4 = new TClonesArray("TLorentzVector", 1000);

	gen_chic_p4 = new TClonesArray("TLorentzVector", 100);
	gen_Jpsi_p4 = new TClonesArray("TLorentzVector", 100);
	gen_Jpsi_photon_p4 = new TClonesArray("TLorentzVector", 100);
	gen_muon_p4 = new TClonesArray("TLorentzVector", 100);
	gen_phot_p4 = new TClonesArray("TLorentzVector", 100);

	chi_p4 = new TClonesArray("TLorentzVector", 100);


	edm::Service < TFileService > fs;
	event_tree = fs->make < TTree >("event_tree", "Tree with event");

	//general information
	event_tree->Branch("runNumber", &runNumber, "runNumber/L");
	event_tree->Branch("eventNumber", &eventNumber, "eventNumber/L");
	event_tree->Branch("nPrimVertices", &nPrimVertices, "nPrimVertices/I");
	event_tree->Branch("muonPerEvent_noCuts", &muonPerEvent, "muonPerEvent/I");
	event_tree->Branch("convPerEvent_noCuts", &convPerTriggeredEvent, "convPerEvent/I");
	event_tree->Branch("dimuonPerEvent", &dimuonPerEvent, "dimuonPerEvent/I");
	event_tree->Branch("chiCandPerEvent", &chiCandPerEvent, "chiCandPerEvent/I");
	event_tree->Branch("ntracks_inEvent", &ntracks_inEvent, "ntracks_inEvent/I");
	event_tree->Branch("hfTowerSum_inEvent", &hfTowerSum_inEvent, "hfTowerSum_inEvent/D");
	event_tree->Branch("Trig_Event_HLTDoubleMuOpen", &Trig_Event_HLTDoubleMuOpen, "Trig_Event_HLTDoubleMuOpen/I");

	//	PV
	event_tree->Branch("pvtx_z", &pvtx_z);
	event_tree->Branch("pvtx_zError", &pvtx_zError);
	event_tree->Branch("pvtx_x", &pvtx_x);
	event_tree->Branch("pvtx_y", &pvtx_y);
	event_tree->Branch("pvtx_nTracks", &pvtx_nTracks);
	event_tree->Branch("pvtx_isFake", &pvtx_isFake);


	// muon
	event_tree->Branch("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
	event_tree->Branch("muonIsHLTDoubleMuOpenFilter", &muonIsHLTDoubleMuOpenFilter);
	event_tree->Branch("muonIsGlobal", &muonIsGlobal);
	event_tree->Branch("muonIsTracker", &muonIsTracker);
	event_tree->Branch("muonIsPF", &muonIsPF);
	event_tree->Branch("muonIsSoft", &muonIsSoft);
	event_tree->Branch("muonIsTight", &muonIsTight);
	event_tree->Branch("muonIsNotGlobalNorTracker", &muonIsNotGlobalNorTracker);
	event_tree->Branch("muonIDHas_TMOneStationTight", &muonIDHas_TMOneStationTight);
	event_tree->Branch("muon_pvtx_index", &muon_pvtx_index);
	event_tree->Branch("muonInnerTrack_dxy", &muonInnerTrack_dxy);
	event_tree->Branch("muonInnerTrack_dz", &muonInnerTrack_dz);
	event_tree->Branch("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);
	event_tree->Branch("muonPixelLayersWithMeasurement", &muonPixelLayersWithMeasurement);
	event_tree->Branch("muonQuality_isHighPurity", &muonQuality_isHighPurity);
	event_tree->Branch("muon_charge", &muon_charge);
	event_tree->Branch("muon_eta", &muon_eta);
	event_tree->Branch("muon_pt", &muon_pt);
	event_tree->Branch("muon_p4", "TClonesArray", &muon_p4, 32000, 0);
	if (flag_doMC) {
		event_tree->Branch("muon_isMatchedMC", &muon_isMatchedMC);
		event_tree->Branch("muonGen_eta", &muonGen_eta);
		event_tree->Branch("muonGen_pt", &muonGen_pt);
		event_tree->Branch("muonGen_p4", &muonGen_p4, 32000, 0);
		event_tree->Branch("muonGen_rDelta", &muonGen_rDelta);
		event_tree->Branch("muonGen_ptDelta", &muonGen_ptDelta);
		event_tree->Branch("muonGen_ptDeltaRel", &muonGen_ptDeltaRel);
	}
	if (flag_saveExtraThings) {
		event_tree->Branch("patMuon", "std::vector <pat::Muon>", &patMuonStored);
	}

	// dimuon

	event_tree->Branch("dimuon_p4", "TClonesArray", &dimuon_p4, 32000, 0);
	event_tree->Branch("dimuon_eta", &dimuon_eta);
	event_tree->Branch("dimuon_pt", &dimuon_pt);
	event_tree->Branch("dimuon_charge", &dimuon_charge);
	event_tree->Branch("dimuon_vtx", "TClonesArray", &dimuon_vtx, 32000, 0);
	event_tree->Branch("dimuon_pvtx_indexFromOniaMuMu", &dimuon_pvtx_indexFromOniaMuMu);
	event_tree->Branch("dimuon_pvtx_index", &dimuon_pvtx_index);
	event_tree->Branch("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);
	event_tree->Branch("dimuon_vtxProb", &dimuon_vtxProb);

	if (flag_saveExtraThings) {
		event_tree->Branch("dimuonStored", "std::vector <pat::CompositeCandidate>", &dimuonStored);
	}
	event_tree->Branch("dimuon_muon1_position", &dimuon_muon1_position);
	event_tree->Branch("dimuon_muon2_position", &dimuon_muon2_position);
	event_tree->Branch("dimuon_ctpv", &dimuon_ctpv);
	event_tree->Branch("dimuon_ctpvError", &dimuon_ctpvError);





	// conversions
	event_tree->Branch("convQuality_isHighPurity", &convQuality_isHighPurity);
	event_tree->Branch("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly);
	event_tree->Branch("conv_vtx", "TClonesArray", &conv_vtx, 32000, 0);
	event_tree->Branch("conv_vertexPositionRho", &conv_vertexPositionRho);
	event_tree->Branch("conv_sigmaTkVtx1", &conv_sigmaTkVtx1);
	event_tree->Branch("conv_sigmaTkVtx2", &conv_sigmaTkVtx2);
	event_tree->Branch("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOK);
	event_tree->Branch("conv_tkVtxCompatible_bestVertex", &conv_tkVtxCompatible_bestVertex);
	event_tree->Branch("conv_tkVtxCompatible_secondBestVertexA", &conv_tkVtxCompatible_secondBestVertexA);
	event_tree->Branch("conv_tkVtxCompatible_secondBestVertexB", &conv_tkVtxCompatible_secondBestVertexB);
	event_tree->Branch("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK);
	if (flag_saveExtraThings) {
		event_tree->Branch("conv_hitPat1", "reco::HitPattern", &conv_hitPat1);
		event_tree->Branch("conv_hitPat2", "reco::HitPattern", &conv_hitPat2);
		event_tree->Branch("conv_tkVtxCompatibilityOK_test", &conv_tkVtxCompatibilityOK_test);
		event_tree->Branch("conv_tkVtxCompatible_bestVertex_test", &conv_tkVtxCompatible_bestVertex_test);
		event_tree->Branch("conv_tkVtxCompatible_secondBestVertexA_test", &conv_tkVtxCompatible_secondBestVertexA_test);
		event_tree->Branch("conv_tkVtxCompatible_secondBestVertexB_test", &conv_tkVtxCompatible_secondBestVertexB_test);
	}
	event_tree->Branch("conv_vertexChi2Prob", &conv_vertexChi2Prob);
	event_tree->Branch("conv_pvtx_index", &conv_pvtx_index);
	event_tree->Branch("conv_zOfPriVtx", &conv_zOfPriVtx);
	event_tree->Branch("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks);
	event_tree->Branch("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx);
	event_tree->Branch("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1);
	event_tree->Branch("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2);
	event_tree->Branch("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1);
	event_tree->Branch("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2);
	event_tree->Branch("conv_dxyError_Tr1", &conv_dxyError_Tr1);
	event_tree->Branch("conv_dxyError_Tr2", &conv_dxyError_Tr2);

	event_tree->Branch("conv_tk1NumOfDOF", &conv_tk1NumOfDOF);
	event_tree->Branch("conv_tk2NumOfDOF", &conv_tk2NumOfDOF);
	event_tree->Branch("conv_track1Chi2", &conv_track1Chi2);
	event_tree->Branch("conv_track2Chi2", &conv_track2Chi2);
	event_tree->Branch("conv_Tr1_pt", &conv_Tr1_pt);
	event_tree->Branch("conv_Tr2_pt", &conv_Tr2_pt);

	event_tree->Branch("conv_minDistanceOfApproach", &conv_minDistanceOfApproach);

	event_tree->Branch("conv_p4", "TClonesArray", &conv_p4, 32000, 0);
	event_tree->Branch("conv_eta", &conv_eta);
	event_tree->Branch("conv_pt", &conv_pt);



	if (flag_doMC) {
		event_tree->Branch("conv_isMatchedMC", &conv_isMatchedMC);
		event_tree->Branch("convGen_eta", &convGen_eta);
		event_tree->Branch("convGen_pt", &convGen_pt);
		event_tree->Branch("convGen_p4", "TClonesArray", &convGen_p4, 32000, 0);
		event_tree->Branch("convGen_rDelta", &convGen_rDelta);
		event_tree->Branch("convGen_ptDelta", &convGen_ptDelta);
		event_tree->Branch("convGen_ptDeltaRel", &convGen_ptDeltaRel);
		event_tree->Branch("convGen_motherCode", &convGen_motherCode);
	}

	
	// MC general
	if (flag_doMC) {
		event_tree->Branch("gen_isGoodChicDecay", &gen_isGoodChicDecay);
		event_tree->Branch("gen_pdgId", &gen_pdgId);

		event_tree->Branch("gen_chic_pt", &gen_chic_pt);
		event_tree->Branch("gen_chic_eta", &gen_chic_eta);
		event_tree->Branch("gen_chic_p4", "TClonesArray", &gen_chic_p4, 32000, 0);
		event_tree->Branch("gen_chic_matchPosition", &gen_chic_matchPosition);
		event_tree->Branch("gen_chic_nMatches", &gen_chic_nMatches);
		event_tree->Branch("gen_chic_rDelta", &gen_chic_rDelta);
		event_tree->Branch("gen_chic_ptDeltaRel", &gen_chic_ptDeltaRel);

		event_tree->Branch("gen_Jpsi_pt", &gen_Jpsi_pt);
		event_tree->Branch("gen_Jpsi_eta", &gen_Jpsi_eta);
		event_tree->Branch("gen_Jpsi_p4", "TClonesArray", &gen_Jpsi_p4, 32000, 0);
		event_tree->Branch("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
		event_tree->Branch("gen_Jpsi_nMatches", &gen_Jpsi_nMatches);
		event_tree->Branch("gen_Jpsi_rDelta", &gen_Jpsi_rDelta);
		event_tree->Branch("gen_Jpsi_ptDeltaRel", &gen_Jpsi_ptDeltaRel);

		event_tree->Branch("gen_Jpsi_photon_n", &gen_Jpsi_photon_n);
		event_tree->Branch("gen_Jpsi_photon_pt", &gen_Jpsi_photon_pt);
		event_tree->Branch("gen_Jpsi_photon_p4", "TClonesArray", &gen_Jpsi_photon_p4, 32000, 0);

		event_tree->Branch("gen_muon_charge", &gen_muon_charge);
		event_tree->Branch("gen_muon_pt", &gen_muon_pt);
		event_tree->Branch("gen_muon_eta", &gen_muon_eta);
		event_tree->Branch("gen_muon_p4", "TClonesArray", &gen_muon_p4, 32000, 0);
		event_tree->Branch("gen_muon_matchPosition", &gen_muon_matchPosition);
		event_tree->Branch("gen_muon_nMatches", &gen_muon_nMatches);
		event_tree->Branch("gen_muon_rDelta", &gen_muon_rDelta);
		event_tree->Branch("gen_muon_ptDeltaRel", &gen_muon_ptDeltaRel);


		event_tree->Branch("gen_phot_pt", &gen_phot_pt);
		event_tree->Branch("gen_phot_eta", &gen_phot_eta);
		event_tree->Branch("gen_phot_p4", "TClonesArray", &gen_phot_p4, 32000, 0);
		event_tree->Branch("gen_conv_matchPosition", &gen_conv_matchPosition);
		event_tree->Branch("gen_conv_nMatches", &gen_conv_nMatches);
		event_tree->Branch("gen_conv_rDelta", &gen_conv_rDelta);
		event_tree->Branch("gen_conv_ptDeltaRel", &gen_conv_ptDeltaRel);
	}


	//chi

	event_tree->Branch("chi_p4", "TClonesArray", &chi_p4, 32000, 0);
	event_tree->Branch("chi_eta", &chi_eta);
	event_tree->Branch("chi_pt", &chi_pt);
	event_tree->Branch("chi_daughterJpsi_position", &chi_daughterJpsi_position);
	event_tree->Branch("chi_daughterConv_position", &chi_daughterConv_position);
	event_tree->Branch("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx);
	event_tree->Branch("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx);
	event_tree->Branch("chiStored", "std::vector <pat::CompositeCandidate>", &chiStored);
	event_tree->Branch("chi_kinematicRefitFlag", &chi_kinematicRefitFlag);
	event_tree->Branch("chi_refit_origChicPosition", &chi_refit_origChicPosition);
	event_tree->Branch("chi_refitStored", "std::vector <pat::CompositeCandidate>", &chi_refitStored);
	event_tree->Branch("chi_refit_vprob", &chi_refit_vprob);
	event_tree->Branch("chi_refit_ctauPV", &chi_refit_ctauPV);
	event_tree->Branch("chi_refit_ctauErrPV", &chi_refit_ctauErrPV);
	event_tree->Branch("chi_refit_ctauPV3D", &chi_refit_ctauPV3D);
	event_tree->Branch("chi_refit_pvtxFromPVwithMuons_x", &chi_refit_pvtxFromPVwithMuons_x);
	event_tree->Branch("chi_refit_pvtxFromPVwithMuons_y", &chi_refit_pvtxFromPVwithMuons_y);
	event_tree->Branch("chi_refit_pvtxFromPVwithMuons_z", &chi_refit_pvtxFromPVwithMuons_z);

}

ChiRootupler::~ChiRootupler() {}

//
// member functions
//


// ------------ method called for each event  ------------
void ChiRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
	using namespace edm;
	using namespace std;

	edm::Handle  < std::vector <pat::Muon> > muon_handle;
	iEvent.getByToken(muon_label, muon_handle);

	edm::Handle < std::vector < pat::CompositeCandidate > >dimuon_handle;
	iEvent.getByToken(dimuon_label, dimuon_handle);

	//edm::Handle < std::vector < pat::CompositeCandidate > >photon_handle;
	//iEvent.getByToken(photon_label, photon_handle);

	edm::Handle < std::vector <reco::Conversion>> conversion_handle;
	iEvent.getByToken(conversion_label, conversion_handle);

	edm::Handle  < reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_label, primaryVertices_handle);

	edm::Handle < edm::TriggerResults > triggerResults_handle;
	iEvent.getByToken(triggerResults_label, triggerResults_handle);

	edm::Handle <reco::Centrality> centrality_handle;
	iEvent.getByToken(centrality_label, centrality_handle);

	edm::Handle <reco::GenParticleCollection> genParticles_handle;
	iEvent.getByToken(genParticles_label, genParticles_handle);



	/////////////////////
	//// S T A R T //////
	////////////////////


	//general info
	runNumber = iEvent.id().run();
	eventNumber = iEvent.id().event();

	///////////////////////////
	/////////  TRIGGER   /////
	//////////////////////////

	// Event trigger
	if (triggerResults_handle.isValid()) {
		const auto triggerIndex = hltConfig.triggerIndex(triggerName);
		if (triggerIndex < triggerResults_handle->size())
		{
			Trig_Event_HLTDoubleMuOpen = triggerResults_handle->accept(triggerIndex);
		}
		else {
			cout << "Trigger not found" << endl;
			return;
		}
	}
	else {
		std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
		return;
	}
	// if not trigger, skip the event (should not happen due to preselection in the main cfg) // Leave for MC
	if (Trig_Event_HLTDoubleMuOpen == 0 && flag_doMC==false) return;



	nPrimVertices = 0;
	if (primaryVertices_handle.isValid()) {
		nPrimVertices = primaryVertices_handle->size();
	}
	muonPerEvent = 0;
	if (muon_handle.isValid()) {
		muonPerEvent = muon_handle->size(); //all muons without any cuts
	}
	else cout << "Problem with muon handle" << endl;

	convPerTriggeredEvent = 0;
	if (conversion_handle.isValid())
	{
		convPerTriggeredEvent = conversion_handle->size(); //all conversions without any cuts
	}
	else cout << "Problem with conv handle" << endl;

	//int photonSize = photon_handle->size();
	//if (photonSize != convPerTriggeredEvent) cout << "Different photon and conversion counts " << photonSize << "    " << convPerTriggeredEvent << endl;


	dimuonPerEvent = 0;
	if (dimuon_handle.isValid())
	{
		dimuonPerEvent = dimuon_handle->size(); //all conversions without any cuts
	}
	else cout << "Problem with dimuon handle" << endl;

	///// CENTRALITY

	if (centrality_handle.isValid())
	{
		ntracks_inEvent = centrality_handle->Ntracks();
		hfTowerSum_inEvent = centrality_handle->EtHFtowerSum();
	}
	else cout << "Problem with centrality handle" << endl;


	if (dimuonPerEvent == 0 && flag_doMC == false) //don't save events without J/psi in them, unless MC
	{
		Clear();
		return;
	}

	//////////////////////
	//    C H I       ////
	/////////////////////

	/// create the collection 
	//cout << "At chic" << endl;
	int n_chic = 0;  // counter for n chic in event
	pat::CompositeCandidateCollection* chiCandColl = new pat::CompositeCandidateCollection;
	for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon_handle->begin(); dimuonCand != dimuon_handle->end(); ++dimuonCand)
	{

		// loop on conversion candidates, make chi cand
		for (reco::ConversionCollection::const_iterator photCand = conversion_handle->begin(); photCand != conversion_handle->end(); ++photCand) {

			if (ConvSelection(*photCand) == false) { continue; } //pre-select photons here (only pT>0.2 and |eta|<2.5)

			chi_cand = makeChiCandidate(*dimuonCand, *photCand);

			//chi cuts
			if (chi_cand.mass() < 2.0 || chi_cand.mass() > 6.0) continue;


			TLorentzVector chi_p4_aux;
			chi_p4_aux.SetPtEtaPhiM(chi_cand.pt(), chi_cand.eta(), chi_cand.phi(), chi_cand.mass());
			new ((*chi_p4)[n_chic]) TLorentzVector(chi_p4_aux);
			++n_chic;
			chi_eta.push_back(chi_p4_aux.Eta());
			chi_pt.push_back(chi_p4_aux.Pt());
			chi_daughterJpsi_position.push_back(std::distance(dimuon_handle->begin(), dimuonCand));
			chi_daughterConv_position.push_back(std::distance(conversion_handle->begin(), photCand));
			chi_dzPhotToDimuonVtx.push_back(Getdz(*photCand, dimuonCand->vertex()));
			chi_dxyPhotToDimuonVtx.push_back(Getdxy(*photCand, dimuonCand->vertex()));
			chiStored.push_back(chi_cand);

			chiCandColl->push_back(chi_cand);

			// KinematicRefit - the following is based on github.com/alberto-sanchez/chi-analysis-miniaod/blob/master/src/OniaPhotonKinematicFit.cc from Alberto  
			int refitStatusFlag = 0;
			bool refitGoodFitFlag = false;
			if (flag_doKinematicRefit) {
				if ((dimuonCand->daughter("muon1")->isTrackerMuon()|| dimuonCand->daughter("muon1")->isGlobalMuon()) && (dimuonCand->daughter("muon2")->isTrackerMuon() || dimuonCand->daughter("muon2")->isGlobalMuon())) //if not, it doesn't have inner track, no way to refit anything
				{
					//cout << "in kinematic refit  " << n_chic<<endl;
					refitStatusFlag += 1;
					edm::ESHandle<TransientTrackBuilder> theB;
					iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);

					reco::TrackRef JpsiTk[2] = { (dynamic_cast<const pat::Muon*>(chi_cand.daughter("dimuon")->daughter("muon1")))->innerTrack(), (dynamic_cast<const pat::Muon*>(chi_cand.daughter("dimuon")->daughter("muon2")))->innerTrack() };

					const reco::Track conv_tk1 = *(photCand->tracks().at(0));
					const reco::Track conv_tk2 = *(photCand->tracks().at(1));
					//photCand->
					//conv.tracks()[0]
					//const reco::Track tk0 = *(dynamic_cast<const pat::CompositeCandidate*>(chi_cand.daughter("photon"))->userData<reco::Track>("track0"));
					//const reco::Track tk1 = *(dynamic_cast<const pat::CompositeCandidate*>(chi_cand.daughter("photon"))->userData<reco::Track>("track1"));
					//convTracks.push_back(conv_tk1);
					//convTracks.push_back(conv_tk2);
					const reco::Vertex thePrimaryV = *(dimuonCand->userData<reco::Vertex>("PVwithmuons"));
					const reco::Vertex& pvtx = primaryVertices_handle->at(0);
					if (abs(thePrimaryV.x() - pvtx.x()) > 0.01) { cout << "primaryV x   " << thePrimaryV.x() << "    " << pvtx.x() << endl; } //check if it is always 0th vertex
					chi_refit_pvtxFromPVwithMuons_x.push_back(thePrimaryV.x());
					chi_refit_pvtxFromPVwithMuons_y.push_back(thePrimaryV.y());
					chi_refit_pvtxFromPVwithMuons_z.push_back(thePrimaryV.z());

					//cout << "dimuon  " << dimuonCand->userInt("muonPosition1") << endl;
					//cout << "dimuon  " << dimuonCand->userInt("muonPosition2") << endl;
					//const edm::RefToBase<reco::Track> conv_tk1 = photCand->tracks().at(0);
					//cout<<"mom  " << conv_tk1.outerP()<<endl;
					//cout << chi_cand.daughter("photon")->pt() << endl;
					//cout << chi_cand.daughter("photon")->userData<reco::Track>("track1")->outerPt() << endl;

					std::vector<reco::TransientTrack> MuMuTT;
					MuMuTT.push_back((*theB).build(&JpsiTk[0]));
					MuMuTT.push_back((*theB).build(&JpsiTk[1]));

					std::vector<reco::TransientTrack> EETT;
					EETT.push_back((*theB).build(conv_tk1));
					EETT.push_back((*theB).build(conv_tk2));

					//values copied from Alberto
					const ParticleMass photon_mass(0);
					float photon_sigmaM = 1E-6;
					const ParticleMass ele_mass(0.000511);
					float ele_sigmaM = 1E-6;
					const ParticleMass muon_mass(0.1056584);
					float muon_sigmaM = muon_mass * 1E-6;

					KinematicParticleFactoryFromTransientTrack pFactory;
					std::vector<RefCountedKinematicParticle> PhotonParticles;
					PhotonParticles.push_back(pFactory.particle(EETT[0], ele_mass, float(0), float(0), ele_sigmaM));
					PhotonParticles.push_back(pFactory.particle(EETT[1], ele_mass, float(0), float(0), ele_sigmaM));

					KinematicParticleVertexFitter fitter;
					RefCountedKinematicTree photonVertexFitTree;
					photonVertexFitTree = fitter.fit(PhotonParticles);

					if (!photonVertexFitTree->isValid()) {
						refitStatusFlag += 2;
						edm::ParameterSet pSet;
						pSet.addParameter<double>("maxDistance", 3);
						pSet.addParameter<int>("maxNbrOfIterations", 10000);
						KinematicParticleVertexFitter fitter2(pSet);
						photonVertexFitTree = fitter2.fit(PhotonParticles);
					}

					if (photonVertexFitTree->isValid()) {
						// now apply Photon mass constraint
						KinematicParticleFitter csFitterPhoton;
						KinematicConstraint * pho_c = new MassKinematicConstraint(photon_mass, photon_sigmaM);

						// add mass constraint to the photon fit to do a constrained fit:
						photonVertexFitTree->movePointerToTheTop();
						photonVertexFitTree = csFitterPhoton.fit(pho_c, photonVertexFitTree);

						photonVertexFitTree->movePointerToTheTop();
						RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
						std::vector<RefCountedKinematicParticle> allChiDaughters;
						allChiDaughters.push_back(pFactory.particle(MuMuTT[0], muon_mass, float(0), float(0), muon_sigmaM));
						allChiDaughters.push_back(pFactory.particle(MuMuTT[1], muon_mass, float(0), float(0), muon_sigmaM));
						allChiDaughters.push_back(fitPhoton);

						KinematicConstrainedVertexFitter constVertexFitter;

						MultiTrackKinematicConstraint *chic_mtc = new  TwoTrackMassKinematicConstraint(jpsi_mass);
						RefCountedKinematicTree ChiKinTree = constVertexFitter.fit(allChiDaughters, chic_mtc);

						if (!ChiKinTree->isEmpty()) {

							ChiKinTree->movePointerToTheTop();
							RefCountedKinematicParticle fitChi = ChiKinTree->currentParticle();
							RefCountedKinematicVertex ChiDecayVertex = ChiKinTree->currentDecayVertex();

							if (fitChi->currentState().isValid()) { //Get Chi         

								refitGoodFitFlag = true;
								float ChiM_fit = fitChi->currentState().mass();
								float ChiPx_fit = fitChi->currentState().kinematicParameters().momentum().x();
								float ChiPy_fit = fitChi->currentState().kinematicParameters().momentum().y();
								float ChiPz_fit = fitChi->currentState().kinematicParameters().momentum().z();
								float ChiVtxX_fit = ChiDecayVertex->position().x();
								float ChiVtxY_fit = ChiDecayVertex->position().y();
								float ChiVtxZ_fit = ChiDecayVertex->position().z();
								float ChiVtxP_fit = ChiSquaredProbability((double)(ChiDecayVertex->chiSquared()),
									(double)(ChiDecayVertex->degreesOfFreedom()));

								reco::CompositeCandidate recoChiRefit(0, math::XYZTLorentzVector(ChiPx_fit, ChiPy_fit, ChiPz_fit,
									sqrt(ChiM_fit*ChiM_fit + ChiPx_fit * ChiPx_fit + ChiPy_fit * ChiPy_fit + ChiPz_fit * ChiPz_fit)), 
									math::XYZPoint(ChiVtxX_fit,	ChiVtxY_fit, ChiVtxZ_fit), PythCode_chic1); //code is for pythyia, chosen without meaning

								chi_refit_origChicPosition.push_back(n_chic-1); //stores position of the chic candidate for the refit (there will be gaps if refit fails)
								pat::CompositeCandidate patChiRefit(recoChiRefit);
								chi_refitStored.push_back(patChiRefit);								
								chi_refit_vprob.push_back(ChiVtxP_fit);

								// lifetime using PV
								TVector3 vtx;
								TVector3 pvtx;
								VertexDistanceXY vdistXY;
								reco::Vertex myVertex = *ChiDecayVertex;

								vtx.SetXYZ(ChiVtxX_fit, ChiVtxY_fit, 0);
								pvtx.SetXYZ(thePrimaryV.position().x(), thePrimaryV.position().y(), 0);
								TVector3 pperp(ChiPx_fit, ChiPy_fit, 0);
								AlgebraicVector3 vpperp(pperp.x(), pperp.y(), 0);

								TVector3 vdiff = vtx - pvtx;
								double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
								Measurement1D distXY = vdistXY.distance(myVertex, thePrimaryV);
								double ctauPV = distXY.value() * cosAlpha * ChiM_fit / pperp.Perp();
								GlobalError v1e = myVertex.error();
								GlobalError v2e = thePrimaryV.error();
								AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
								double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp, vXYe)) * ChiM_fit / (pperp.Perp2());

								chi_refit_ctauPV.push_back(ctauPV);
								chi_refit_ctauErrPV.push_back(ctauErrPV);

								// same but 3D
								vtx.SetXYZ(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit);
								pvtx.SetXYZ(thePrimaryV.position().x(), thePrimaryV.position().y(), thePrimaryV.position().z());
								TVector3 pperp3D(ChiPx_fit, ChiPy_fit, ChiPz_fit);
								TVector3 vdiff3D = vtx - pvtx;
								cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
								VertexDistance3D vdistXYZ;
								distXY = vdistXYZ.distance(myVertex, thePrimaryV);
								ctauPV = distXY.value() * cosAlpha * ChiM_fit / pperp.Perp();
								chi_refit_ctauPV3D.push_back(ctauPV);
							}
							else refitStatusFlag += 16;
						}
						else refitStatusFlag += 8;
					}
					else refitStatusFlag += 4;



				}
				else
				{
				refitStatusFlag = -1;
				}
			}

			if (refitGoodFitFlag == false) //fill some default value for direct looking at the branches
			{
				chi_refit_vprob.push_back(-1);
				chi_refit_ctauPV.push_back(-100);
				chi_refit_ctauErrPV.push_back(-100);
				chi_refit_ctauPV3D.push_back(-100);
			}
			chi_kinematicRefitFlag.push_back(refitStatusFlag);

		}

	}
	chiCandPerEvent = n_chic;
	


	//PV
	if (primaryVertices_handle.isValid()) {
		//if (primaryVertices_handle->size() == 1) return; //test, to be deleted
		for (uint i = 0; i < primaryVertices_handle->size(); i++) {
			const reco::Vertex& pvtx = primaryVertices_handle->at(i);
			pvtx_z.push_back(pvtx.z());
			pvtx_zError.push_back(pvtx.zError());
			pvtx_x.push_back(pvtx.x());
			pvtx_y.push_back(pvtx.y());
			pvtx_nTracks.push_back(pvtx.nTracks());
			pvtx_isFake.push_back(pvtx.isFake());

		}

	}
	else cout << "Problem with PV handle" << endl;

	int pvtx_index = 0; // 0: top primary vertex. Set in dimuon loop to vertex that is "closest" in z to dimuon vertex  (with optional cut-off, see ChiRootupler::SelectVertex) //changed to match vertex from HiOnia2MuMuPAT (mostly same version of the closest in dz)

	////////////////////////
	////  D I M U O N   ///
	//////////////////////

	if (dimuon_handle.isValid())
	{
		//cout << "At Dimuon" << endl;
		for (uint i = 0; i < dimuon_handle->size(); i++) {
			const pat::CompositeCandidate& dimuon = dimuon_handle->at(i);
			TLorentzVector dimuon_p4_aux;
			dimuon_p4_aux.SetPtEtaPhiM(dimuon.pt(), dimuon.eta(), dimuon.phi(), dimuon.mass());
			new ((*dimuon_p4)[i]) TLorentzVector(dimuon_p4_aux);

			dimuon_charge.push_back(dimuon.charge());
			dimuon_eta.push_back(dimuon.eta());
			dimuon_pt.push_back(dimuon.pt());

			if (flag_saveExtraThings)
			{
				dimuonStored.push_back(dimuon);
			}

			const reco::Vertex* dimuon_recovtx = dimuon.userData<reco::Vertex>("commonVertex");
			TVector3 dimuon_vtx_aux;
			dimuon_vtx_aux.SetXYZ(dimuon_recovtx->x(), dimuon_recovtx->y(), dimuon_recovtx->z());
			new ((*dimuon_vtx)[i]) TVector3(dimuon_vtx_aux);
			dimuon_pvtx_indexFromOniaMuMu.push_back(dimuon.userInt("PVposition"));//PVposition as got from the HIOnia2MuMuPAT
			pvtx_index = SelectVertex(primaryVertices_handle, dimuon_recovtx->z());
			dimuon_pvtx_index.push_back(pvtx_index);
			pvtx_index = dimuon.userInt("PVposition"); //use the one selected by HiOnia2MuMuPAT
			dimuon_dz_dimuonvtx_pvtx.push_back(dimuon_recovtx->z() - (*primaryVertices_handle.product())[pvtx_index].z());
			dimuon_vtxProb.push_back(dimuon.userFloat("vProb"));

			int muonPos1 = dimuon.userInt("muonPosition1");
			int muonPos2 = dimuon.userInt("muonPosition2");
			if (fabs(muon_handle->at(muonPos2).pt() - dimuon.daughter("muon2")->pt()) > 0.01) { cout << "SOMETHING WRONG WITH THE MATCHING FROM DIMUON TO MUON - possibly different muon collections used" << endl; }//muon matching assumes that muon collection going to dimuon producer and here is the same
			dimuon_muon1_position.push_back(muonPos1);
			dimuon_muon2_position.push_back(muonPos2);
			dimuon_ctpv.push_back(dimuon.userFloat("ppdlPV"));
			dimuon_ctpvError.push_back(dimuon.userFloat("ppdlErrPV"));
			//cout << "Muon positions  " << muonPos1 << "   " << muonPos2 << endl;

		}
	}
	else cout << "Problem with dimuon handle" << endl;


	/////////////////////
	//   M U O N S   ////
	////////////////////
	if (muon_handle.isValid()) {
		//cout << "At muons" << endl;
		for (uint i = 0; i < muon_handle->size(); i++) {
			const pat::Muon& patMuon = muon_handle->at(i);
			// pick the closest vertex of those that were selected by dimuon (and ignore the vertices that had no dimuon in them)
						
			int bestPvtx_index = 0;
			float minDz = 9999;
			if (patMuon.isTrackerMuon() == true || patMuon.isGlobalMuon() == true) //If the muon has an inner track, we can see if it comes from a different vertex than the one with index 0. If it is just standalone muon, it doesn't change anything.
			{
				for (uint iVtx = 0; iVtx < dimuon_pvtx_index.size(); iVtx++)
				{
					float deltaZ = fabs(patMuon.innerTrack()->dz() - primaryVertices_handle->at(dimuon_pvtx_index.at(iVtx)).z());
					if (deltaZ < minDz) {
						minDz = deltaZ;
						bestPvtx_index = dimuon_pvtx_index.at(iVtx);
					}
				}
			}
			//bestPvtx_index = 0; //disable

			muon_pvtx_index.push_back(bestPvtx_index);

			
			const reco::Vertex& pvtx = primaryVertices_handle->at(bestPvtx_index); //use our selected PV
			if (patMuon.triggerObjectMatchByPath(triggerName) != 0) {
				muonIsHLTDoubleMuOpen.push_back(true);
			}
			else { muonIsHLTDoubleMuOpen.push_back(false); }
			if (patMuon.triggerObjectMatchByFilter(triggerFilter) != 0) {
				muonIsHLTDoubleMuOpenFilter.push_back(true);
			}
			else { muonIsHLTDoubleMuOpenFilter.push_back(false); }
			muonIsGlobal.push_back(patMuon.isGlobalMuon());
			muonIsTracker.push_back(patMuon.isTrackerMuon());
			muonIsPF.push_back(patMuon.isPFMuon());
			muonIsSoft.push_back(patMuon.isSoftMuon(pvtx));
			muonIsTight.push_back(patMuon.isTightMuon(pvtx));
			if (!patMuon.isGlobalMuon() && !patMuon.isTrackerMuon()) { muonIsNotGlobalNorTracker.push_back(true); }
			else muonIsNotGlobalNorTracker.push_back(false); // just for convenience

			muonIDHas_TMOneStationTight.push_back(patMuon.muonID("TMOneStationTight"));
			if (patMuon.isGlobalMuon() || patMuon.isTrackerMuon()) {
				muonInnerTrack_dxy.push_back(patMuon.innerTrack()->dxy());
				muonInnerTrack_dz.push_back(patMuon.innerTrack()->dz());
				muonTrackerLayersWithMeasurement.push_back(patMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement());
				muonPixelLayersWithMeasurement.push_back(patMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement());
				reco::TrackBase::TrackQuality tq = reco::TrackBase::qualityByName("highPurity");//high purity=2 //see DataFormats/TrackReco/interface/TrackBase.h
				muonQuality_isHighPurity.push_back(patMuon.innerTrack()->quality(tq));
			}
			else
			{
				muonInnerTrack_dxy.push_back(-100);
				muonInnerTrack_dz.push_back(-100);
				muonTrackerLayersWithMeasurement.push_back(-1);
				muonPixelLayersWithMeasurement.push_back(-1);
				muonQuality_isHighPurity.push_back(0);
			}
			muon_charge.push_back(patMuon.charge());
			muon_eta.push_back(patMuon.eta());
			muon_pt.push_back(patMuon.pt());
			TLorentzVector muon_p4_aux;
			muon_p4_aux.SetPtEtaPhiM(patMuon.pt(), patMuon.eta(), patMuon.phi(), patMuon.mass());
			new ((*muon_p4)[i]) TLorentzVector(muon_p4_aux);
			// MC generated information
			//cout << "At muonMC" << endl;
			if (flag_doMC) {
				bool muon_isMatchedMC_flag = patMuon.genLepton(); //false for null pointer, true if match exists, from pat::Lepton.h //O: Stopped working in MC v5, since not really needed, left out
				//cout << muon_isMatchedMC_flag << endl;
				muon_isMatchedMC.push_back(muon_isMatchedMC_flag);
				if (muon_isMatchedMC_flag) {
					const reco::GenParticle genMuon = *patMuon.genLepton();
					//cout << "Got gen particle " << endl;
					muonGen_eta.push_back(genMuon.eta());
					muonGen_pt.push_back(genMuon.pt());
					TLorentzVector muonGen_p4_aux;
					muonGen_p4_aux.SetPtEtaPhiM(genMuon.pt(), genMuon.eta(), genMuon.phi(), genMuon.mass());
					new ((*muonGen_p4)[i]) TLorentzVector(muonGen_p4_aux);
					muonGen_rDelta.push_back(reco::deltaR(patMuon, genMuon)); //sqrt(phi^2+eta^2)
					muonGen_ptDelta.push_back(patMuon.pt() - genMuon.pt());
					muonGen_ptDeltaRel.push_back((patMuon.pt() - genMuon.pt()) / genMuon.pt()); //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
				}
				else { //default values to store if no match - in principle can be ommitted, just for direct looking at branches
					muonGen_eta.push_back(-5);
					muonGen_pt.push_back(0);
					new ((*muonGen_p4)[i]) TLorentzVector(TLorentzVector());
					muonGen_rDelta.push_back(-5);
					muonGen_ptDelta.push_back(-5);
					muonGen_ptDeltaRel.push_back(-5);
				}
			}
			//cout << "out of muonMC" << endl;

			if (flag_saveExtraThings)
			{
				patMuonStored.push_back(patMuon);
			}

		}
	}
	else cout << "Problem with muon handle" << endl;


	///////////////////////////////////////////
	///////   C O N V E R S I O N S   ////////
	/////////////////////////////////////////
	if (conversion_handle.isValid())
	{
		//cout << "At conversions" << endl;
		for (uint i = 0; i < conversion_handle->size(); i++) {
			const reco::Conversion& candPhoton = conversion_handle->at(i);
			TLorentzVector conv_p4_aux;
			conv_p4_aux.SetXYZT(candPhoton.refittedPair4Momentum().x(), candPhoton.refittedPair4Momentum().y(), candPhoton.refittedPair4Momentum().z(), candPhoton.refittedPair4Momentum().t());
			new ((*conv_p4)[i]) TLorentzVector(conv_p4_aux);
			conv_eta.push_back(conv_p4_aux.Eta());
			conv_pt.push_back(conv_p4_aux.Pt());

			convQuality_isHighPurity.push_back(candPhoton.quality((reco::Conversion::ConversionQuality)(8))); //8 is high purity, see reco::Conversion Class Reference
			convQuality_isGeneralTracksOnly.push_back(candPhoton.quality((reco::Conversion::ConversionQuality)(0))); //0 is general tracks only, see reco::Conversion Class Reference
			const reco::Vertex conv_recovtx = candPhoton.conversionVertex();
			TVector3 conv_vtx_aux;
			conv_vtx_aux.SetXYZ(conv_recovtx.x(), conv_recovtx.y(), conv_recovtx.z());
			new ((*conv_vtx)[i]) TVector3(conv_vtx_aux);
			conv_vertexPositionRho.push_back(candPhoton.conversionVertex().position().rho());
			bool conv_tkVtxCompatible_bestVertex_aux, conv_tkVtxCompatible_secondBestVertexA_aux, conv_tkVtxCompatible_secondBestVertexB_aux;
			double conv_sigmaTkVtx1_aux, conv_sigmaTkVtx2_aux;
			if (flag_saveExtraThings)
			{
				conv_tkVtxCompatibilityOK_test.push_back(Conv_checkTkVtxCompatibility(candPhoton, *primaryVertices_handle.product(), 20, conv_tkVtxCompatible_bestVertex_aux, conv_tkVtxCompatible_secondBestVertexA_aux, conv_tkVtxCompatible_secondBestVertexB_aux, conv_sigmaTkVtx1_aux, conv_sigmaTkVtx2_aux));
				conv_tkVtxCompatible_bestVertex_test.push_back(conv_tkVtxCompatible_bestVertex_aux);
				conv_tkVtxCompatible_secondBestVertexA_test.push_back(conv_tkVtxCompatible_secondBestVertexA_aux);
				conv_tkVtxCompatible_secondBestVertexB_test.push_back(conv_tkVtxCompatible_secondBestVertexB_aux);
			}
			conv_tkVtxCompatibilityOK.push_back(Conv_checkTkVtxCompatibility(candPhoton, *primaryVertices_handle.product(), conv_TkVtxCompSigmaCut, conv_tkVtxCompatible_bestVertex_aux, conv_tkVtxCompatible_secondBestVertexA_aux, conv_tkVtxCompatible_secondBestVertexB_aux, conv_sigmaTkVtx1_aux, conv_sigmaTkVtx2_aux));
			conv_tkVtxCompatible_bestVertex.push_back(conv_tkVtxCompatible_bestVertex_aux);
			conv_tkVtxCompatible_secondBestVertexA.push_back(conv_tkVtxCompatible_secondBestVertexA_aux);
			conv_tkVtxCompatible_secondBestVertexB.push_back(conv_tkVtxCompatible_secondBestVertexB_aux);
			conv_sigmaTkVtx1.push_back(conv_sigmaTkVtx1_aux);
			conv_sigmaTkVtx2.push_back(conv_sigmaTkVtx2_aux);

			if (candPhoton.tracks().size() == 2) {
				const edm::RefToBase<reco::Track> conv_tk1 = candPhoton.tracks().at(0);
				const edm::RefToBase<reco::Track> conv_tk2 = candPhoton.tracks().at(1);

				reco::HitPattern hitPatA = conv_tk1->hitPattern();
				reco::HitPattern hitPatB = conv_tk2->hitPattern();
				conv_hitPat1.push_back(hitPatA);
				conv_hitPat2.push_back(hitPatB);
				conv_compatibleInnerHitsOK.push_back((Conv_foundCompatibleInnerHits(hitPatA, hitPatB) && Conv_foundCompatibleInnerHits(hitPatB, hitPatA)));


				// pick the closest vertex of those that were selected by dimuon (and ignore the vertices that had no dimuon in them)
				int bestPvtx_index = 0;
				float minDz = 9999;
				for (uint iVtx = 0; iVtx < dimuon_pvtx_index.size(); iVtx++)
				{
					float deltaZ = fabs(candPhoton.zOfPrimaryVertexFromTracks(primaryVertices_handle->at(dimuon_pvtx_index.at(iVtx)).position()) - primaryVertices_handle->at(dimuon_pvtx_index.at(iVtx)).z());
					if (deltaZ < minDz) {
						minDz = deltaZ;
						bestPvtx_index = dimuon_pvtx_index.at(iVtx);
					}
				}
				conv_pvtx_index.push_back(bestPvtx_index);

				conv_zOfPriVtx.push_back((*primaryVertices_handle.product())[bestPvtx_index].z());
				conv_zOfPriVtxFromTracks.push_back(candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[bestPvtx_index].position()));
				conv_dzToClosestPriVtx.push_back(candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[bestPvtx_index].position()) - (*primaryVertices_handle.product())[bestPvtx_index].z());
				// Now check impact parameter wrt primary vertex
				conv_dxyPriVtx_Tr1.push_back(conv_tk1->dxy((*primaryVertices_handle.product())[bestPvtx_index].position()));
				conv_dxyPriVtx_Tr2.push_back(conv_tk2->dxy((*primaryVertices_handle.product())[bestPvtx_index].position()));
				conv_dxyPriVtxTimesCharge_Tr1.push_back(conv_tk1->dxy((*primaryVertices_handle.product())[bestPvtx_index].position())*conv_tk1->charge());
				conv_dxyPriVtxTimesCharge_Tr2.push_back(conv_tk2->dxy((*primaryVertices_handle.product())[bestPvtx_index].position())*conv_tk2->charge());
				conv_dxyError_Tr1.push_back(conv_tk1->dxyError());
				conv_dxyError_Tr2.push_back(conv_tk2->dxyError());

				conv_tk1NumOfDOF.push_back(conv_tk1->ndof());
				conv_tk2NumOfDOF.push_back(conv_tk2->ndof());
				conv_track1Chi2.push_back(candPhoton.tracks().at(0)->normalizedChi2());
				conv_track2Chi2.push_back(candPhoton.tracks().at(1)->normalizedChi2());
				conv_Tr1_pt.push_back(conv_tk1->pt());
				conv_Tr2_pt.push_back(conv_tk2->pt());

			}
			else conv_compatibleInnerHitsOK.push_back(-1);

			conv_vertexChi2Prob.push_back(ChiSquaredProbability(candPhoton.conversionVertex().chi2(), candPhoton.conversionVertex().ndof()));
			conv_minDistanceOfApproach.push_back(candPhoton.distOfMinimumApproach());
			//if (candPhoton.distOfMinimumApproach() > -10 && candPhoton.distOfMinimumApproach() < 10) { conv_minDistanceOfApproach = candPhoton.distOfMinimumApproach(); }
			//else conv_minDistanceOfApproach = 0;

			//MC for conversions // for testing, normally useless information
			if (flag_doMC)
			{
				if (genParticles_handle.isValid()) {
					reco::GenParticle genConv_best = reco::GenParticle();
					bool conv_isMatchedMC_aux = false;
					for (uint i = 0; i < genParticles_handle->size(); i++) {
						const reco::GenParticle& genParticle = genParticles_handle->at(i);

						int pdgId = genParticle.pdgId();
						if (pdgId != 22) { continue; } //if not photon, don't bother
						if (genParticle.status() != 1) {
							//cout << "notStable" << endl;
							//cout << genParticle.daughter(0)->pdgId() << endl;
							//cout << "ptdif " << genParticle.daughter(0)->pt() - genParticle.pt() << endl;
							continue;
						} //if not stable, don't bother

						bool genParticleMatched = false; //this particular - is it matched?
						genParticleMatched = Conv_isMatched(candPhoton.refittedPair4Momentum(), genParticle, conv_maxDeltaR, conv_maxDPtRel);
						if (genParticleMatched == true) {
							if (conv_isMatchedMC_aux == false) {//first one found
								conv_isMatchedMC_aux = true;
								genConv_best = genParticle;
							}
							else { //check whether the second match is better by deltaR than the first, save the better one
								if (reco::deltaR(candPhoton.refittedPair4Momentum(), genConv_best) > reco::deltaR(candPhoton.refittedPair4Momentum(), genParticle)) { genConv_best = genParticle; }
							}
						}
					}
					if (conv_isMatchedMC_aux) {
						conv_isMatchedMC.push_back(true);
						convGen_eta.push_back(genConv_best.eta());
						convGen_pt.push_back(genConv_best.pt());
						TLorentzVector convGen_p4_aux;
						convGen_p4_aux.SetPtEtaPhiM(genConv_best.pt(), genConv_best.eta(), genConv_best.phi(), genConv_best.mass());
						new ((*convGen_p4)[i]) TLorentzVector(convGen_p4_aux);

						convGen_rDelta.push_back(reco::deltaR(candPhoton.refittedPair4Momentum(), genConv_best)); //sqrt(phi^2+eta^2)
						convGen_ptDelta.push_back(candPhoton.refittedPair4Momentum().pt() - genConv_best.pt());
						convGen_ptDeltaRel.push_back((candPhoton.refittedPair4Momentum().pt() - genConv_best.pt()) / genConv_best.pt()); //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
						convGen_motherCode.push_back(genConv_best.mother()->pdgId());
					}
					else { //default values to store if no match - in principle can be ommitted, just for direct looking at branches
						conv_isMatchedMC.push_back(false);
						convGen_eta.push_back(-5);
						convGen_pt.push_back(0);
						new ((*convGen_p4)[i]) TLorentzVector(TLorentzVector());
						convGen_rDelta.push_back(-5);
						convGen_ptDelta.push_back(-5);
						convGen_ptDeltaRel.push_back(-5);
						convGen_motherCode.push_back(-1);
					}
				}
				else cout << "Problem with gen handle" << endl;
			}
		}
	}
	else { cout << "Conversions handle problem" << endl; }


	////////////////
	////// GEN ////
	///////////////

	if (flag_doMC) {
		if (genParticles_handle.isValid()) {
			//cout << endl<< "entered Gen part, evt: " << eventNumber<< endl;
			//for (uint i = 0; i < genParticles_handle->size(); i++) {
			//	const reco::GenParticle& genParticle = genParticles_handle->at(i);
			//	int gen_pdgId_aux = genParticle.pdgId();
			//	if (gen_pdgId_aux == 443 || gen_pdgId_aux == -13 || gen_pdgId_aux == 13)
			//	{
			//		cout << "Gen pdg ID: " << gen_pdgId_aux << endl;
			//		if (gen_pdgId_aux == 443)
			//		{
			//			cout << "  Is last copy: " << genParticle.isLastCopy() << " and has daughters " << genParticle.numberOfDaughters() << endl;
			//		}
			//	}
			//}
			//bool foundChicInEvent=false;
			for (uint i = 0; i < genParticles_handle->size(); i++) {
				const reco::GenParticle& genParticle = genParticles_handle->at(i);
				int gen_pdgId_aux = genParticle.pdgId();
				if (gen_pdgId_aux != PythCode_chic0 && gen_pdgId_aux != PythCode_chic1 && gen_pdgId_aux != PythCode_chic2) { continue; } //if not chi, skip
				if (genParticle.isLastCopy() == false) { continue; } // if not last copy of a given particle, skip
				//foundChicInEvent
				// it is the chic
				bool isGoodMCchic = true; //we have some weird decays, skip those
				int nDaughters = genParticle.numberOfDaughters();
				if (nDaughters != 2) { //all of them should be J/psi + gamma
					//cout << "weird gen decay, number of daughters is: " << nDaughters << "  in event number : " << eventNumber << endl;
					isGoodMCchic = false; 
				} 
				for (int j = 0; j < nDaughters; j++) {
					const reco::Candidate& gen_chiDaughter = *genParticle.daughter(j);
					int dauId = gen_chiDaughter.pdgId();
					if (dauId != 443 && dauId != 22) { isGoodMCchic = false; } //not J/psi gamma decay

					if (dauId == 443) {
						int nDaughtJpsi = gen_chiDaughter.numberOfDaughters();
						for (int k = 0; k < nDaughtJpsi; k++) { // loop over jpsi daughters
							const reco::Candidate& gen_chiMuon = *gen_chiDaughter.daughter(k);
							if (abs(gen_chiMuon.pdgId()) != 13 && gen_chiMuon.pdgId() != 22)
							{
								//cout << "Warning, jpsi decay daughter is not a muon nor photon" << endl; //continue; 
								//cout << "ID: " << gen_chiMuon.pdgId() << "  and pt:  " << gen_chiMuon.pt() << "  in event number : " << eventNumber << endl;
								isGoodMCchic = false;
								continue;
							}
						}
					}
				}
				gen_isGoodChicDecay.push_back(isGoodMCchic);
				if (isGoodMCchic == false) { continue; } //only few bad MC decays (as of v6), thus we just skip them


				gen_pdgId.push_back(gen_pdgId_aux);
				gen_chic_pt.push_back(genParticle.pt());
				gen_chic_eta.push_back(genParticle.eta());
				TLorentzVector gen_chic_p4_aux;
				gen_chic_p4_aux.SetXYZT(genParticle.p4().x(), genParticle.p4().y(), genParticle.p4().z(), genParticle.p4().t());
				new ((*gen_chic_p4)[0]) TLorentzVector(gen_chic_p4_aux); //assuming one chic in the event

				int nMatches_chic_aux;
				double chic_rDelta_aux, chic_ptDeltaRel_aux;
				gen_chic_matchPosition.push_back(MatchGen(genParticle, chiCandColl, jpsi_maxDeltaR, jpsi_maxDPtRel, nMatches_chic_aux, chic_rDelta_aux, chic_ptDeltaRel_aux));
				gen_chic_nMatches.push_back(nMatches_chic_aux);
				gen_chic_rDelta.push_back(chic_rDelta_aux); 
				gen_chic_ptDeltaRel.push_back(chic_ptDeltaRel_aux);


				//loop over the daughters
				for (int j = 0; j < nDaughters; j++) {
					const reco::Candidate& gen_chiDaughter = *genParticle.daughter(j);
					//if (gen_chiDaughter.isLastCopy() == false) { cout << "Warning: daughter is not a last copy!" << endl; } //doesn't work
					int dauId = gen_chiDaughter.pdgId();
					if (dauId == 443) {//Jpsi
						gen_Jpsi_pt.push_back(gen_chiDaughter.pt());
						gen_Jpsi_eta.push_back(gen_chiDaughter.eta());
						TLorentzVector gen_Jpsi_p4_aux;
						gen_Jpsi_p4_aux.SetXYZT(gen_chiDaughter.p4().x(), gen_chiDaughter.p4().y(), gen_chiDaughter.p4().z(), gen_chiDaughter.p4().t());
						new ((*gen_Jpsi_p4)[0]) TLorentzVector(gen_Jpsi_p4_aux);
						int nMatches_aux;
						double jpsi_rDelta_aux, jpsi_ptDeltaRel_aux;
						gen_Jpsi_matchPosition.push_back(MatchGen(gen_chiDaughter, dimuon_handle, jpsi_maxDeltaR, jpsi_maxDPtRel, nMatches_aux, jpsi_rDelta_aux, jpsi_ptDeltaRel_aux));
						gen_Jpsi_nMatches.push_back(nMatches_aux);
						gen_Jpsi_rDelta.push_back(jpsi_rDelta_aux);
						gen_Jpsi_ptDeltaRel.push_back(jpsi_ptDeltaRel_aux);

						//check the decay muons

						int gen_Jpsi_nPhoton = 0;
						int nDaughtJpsi = gen_chiDaughter.numberOfDaughters();
						for (int k = 0; k < nDaughtJpsi; k++) { // loop over jpsi daughters
							const reco::Candidate& gen_chiMuon = *gen_chiDaughter.daughter(k);
							//cout << "ID: "<<gen_chiMuon.pdgId() << endl;

							if (abs(gen_chiMuon.pdgId()) != 13 && gen_chiMuon.pdgId() != 22)
							{ 
								//cout << "Warning, jpsi decay daughter is not a muon nor photon" << endl; //continue; 
								//cout << "ID: " << gen_chiMuon.pdgId()<< "  and pt:  " << gen_chiMuon.pt() << "  in event number : "<< eventNumber << endl;
								continue;
							}

							if (gen_chiMuon.pdgId() == 22) // if it is the extra photon
							{
								gen_Jpsi_photon_pt.push_back(gen_chiMuon.pt());
								TLorentzVector gen_phot_p4_aux;
								gen_phot_p4_aux.SetXYZT(gen_chiMuon.p4().x(), gen_chiMuon.p4().y(), gen_chiMuon.p4().z(), gen_chiMuon.p4().t());
								new ((*gen_Jpsi_photon_p4)[gen_Jpsi_nPhoton]) TLorentzVector(gen_phot_p4_aux);
								gen_Jpsi_nPhoton++;
								continue;
							}
							
							gen_muon_charge.push_back(gen_chiMuon.charge()); 
							gen_muon_pt.push_back(gen_chiMuon.pt());
							gen_muon_eta.push_back(gen_chiMuon.eta());
							TLorentzVector gen_muon_p4_aux;
							gen_muon_p4_aux.SetXYZT(gen_chiMuon.p4().x(), gen_chiMuon.p4().y(), gen_chiMuon.p4().z(), gen_chiMuon.p4().t());
							new ((*gen_muon_p4)[k]) TLorentzVector(gen_muon_p4_aux); //Assumes that the muon pair is the first two daughters, otherwise would need a separate counter as is done for photons from Jpsi
							int nMatchesMuon_aux;
							double muon_rDelta_aux, muon_ptDeltaRel_aux;
							gen_muon_matchPosition.push_back(MatchGen(gen_chiMuon, muon_handle, muon_maxDeltaR, muon_maxDPtRel, nMatchesMuon_aux, muon_rDelta_aux, muon_ptDeltaRel_aux));
							gen_muon_nMatches.push_back(nMatchesMuon_aux);
							gen_muon_rDelta.push_back(muon_rDelta_aux);
							gen_muon_ptDeltaRel.push_back(muon_ptDeltaRel_aux);
						}
						gen_Jpsi_photon_n.push_back(gen_Jpsi_nPhoton);
					}
					else if (dauId == 22) {//photon
						gen_phot_pt.push_back(gen_chiDaughter.pt());
						gen_phot_eta.push_back(gen_chiDaughter.eta());
						TLorentzVector gen_phot_p4_aux;
						gen_phot_p4_aux.SetXYZT(gen_chiDaughter.p4().x(), gen_chiDaughter.p4().y(), gen_chiDaughter.p4().z(), gen_chiDaughter.p4().t());
						new ((*gen_phot_p4)[0]) TLorentzVector(gen_phot_p4_aux);
						int nMatches_aux;
						double conv_rDelta_aux, conv_ptDeltaRel_aux;
						gen_conv_matchPosition.push_back(MatchGen(gen_chiDaughter, conversion_handle, conv_maxDeltaR, conv_maxDPtRel, nMatches_aux, conv_rDelta_aux, conv_ptDeltaRel_aux));
						gen_conv_nMatches.push_back(nMatches_aux);
						gen_conv_rDelta.push_back(conv_rDelta_aux);
						gen_conv_ptDeltaRel.push_back(conv_ptDeltaRel_aux);

						//MatchGen(gen_chiDaughter, conv_maxDeltaR, conv_maxDPtRel);

						//int photGen_motherCode = gen_chiDaughter->mother()->pdgId();
						//int photGen_status = gen_chiDaughter->status();
						//cout << "Photon " << photGen_motherCode << " " << photGen_status << endl;
					}
					else {
						cout << "Different decay: " << dauId << endl;
					}
				}

				//{ cout << " event and GEN ID, status and n daughters, first id " << eventNumber << "   " << genParticle.pdgId()<< "   " <<st << "    " <<n << "   " << genParticle.daughter(0)->pdgId()<< "  pt:" << genParticle.pt()<< endl; }
				//if (nChicPerEvent > 1) {
					//cout << endl << endl << endl;
				//}

				break; //We found chic, no need to go through the rest of the collection. However, very few events (<1%) exist with multiple chic in them - this takes only the first one into account
			}
			//cout << "Done with Gen part, evt: " << eventNumber << endl;
		}
		else cout << "Problem with gen handle" << endl;
	}


	event_tree->Fill();
	Clear();

}

/////////////////////////////////////////////
///// G E N E R A L   F U N C T I O N S ///// 
/////////////////////////////////////////////

// ------------ method called once each job just before starting event loop  ------------
void ChiRootupler::beginJob() 
{

}

// ------------ method called once each job just after ending the event loop  ------------
void ChiRootupler::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void ChiRootupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) 
{
	bool changed = true;
	hltConfig.init(iRun, iSetup, "HLT", changed);
}

// ------------ method called when ending the processing of a run  ------------
void ChiRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void ChiRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void ChiRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ clear the vectors after each event
void ChiRootupler::Clear()
{
	//PV

	pvtx_z.clear();
	pvtx_zError.clear();
	pvtx_x.clear();
	pvtx_y.clear();
	pvtx_nTracks.clear();
	pvtx_isFake.clear();


	//muon
	muonIsHLTDoubleMuOpen.clear();
	muonIsHLTDoubleMuOpenFilter.clear();
	muonIsGlobal.clear();
	muonIsTracker.clear();
	muonIsPF.clear();
	muonIsSoft.clear();
	muonIsTight.clear();
	muonIsNotGlobalNorTracker.clear();
	muonIDHas_TMOneStationTight.clear();
	muon_pvtx_index.clear();
	muonInnerTrack_dxy.clear();
	muonInnerTrack_dz.clear();
	muonTrackerLayersWithMeasurement.clear();
	muonPixelLayersWithMeasurement.clear();
	muonQuality_isHighPurity.clear();
	muon_charge.clear();
	muon_eta.clear();
	muon_pt.clear();
	muon_p4->Clear(); 
	patMuonStored.clear();

	//muon MC
	muon_isMatchedMC.clear();
	muonGen_eta.clear();
	muonGen_pt.clear();
	muonGen_p4->Clear();
	muonGen_rDelta.clear();
	muonGen_ptDelta.clear();
	muonGen_ptDeltaRel.clear();
	
	//dimuon 
	dimuon_p4->Clear();
	dimuon_eta.clear();
	dimuon_pt.clear();
	dimuon_charge.clear(); 
	dimuon_vtx->Clear();
	dimuon_pvtx_indexFromOniaMuMu.clear();
	dimuon_pvtx_index.clear();
	dimuon_dz_dimuonvtx_pvtx.clear();
	dimuon_vtxProb.clear();
	dimuonStored.clear();
	dimuon_muon1_position.clear();
	dimuon_muon2_position.clear();
	dimuon_ctpv.clear();
	dimuon_ctpvError.clear();

	//conv
	convQuality_isHighPurity.clear();
	convQuality_isGeneralTracksOnly.clear();
	conv_vtx->Clear();

	conv_vertexPositionRho.clear();
	conv_sigmaTkVtx1.clear();
	conv_sigmaTkVtx2.clear();
	conv_tkVtxCompatibilityOK.clear();
	conv_tkVtxCompatible_bestVertex.clear();
	conv_tkVtxCompatible_secondBestVertexA.clear();
	conv_tkVtxCompatible_secondBestVertexB.clear();
	conv_tkVtxCompatibilityOK_test.clear(); 
	conv_tkVtxCompatible_bestVertex_test.clear(); 
	conv_tkVtxCompatible_secondBestVertexA_test.clear(); 
	conv_tkVtxCompatible_secondBestVertexB_test.clear(); 

	conv_compatibleInnerHitsOK.clear();
	conv_hitPat1.clear();
	conv_hitPat2.clear();
	conv_isCustomHighPurity.clear();
	conv_pvtx_index.clear();
	conv_zOfPriVtx.clear();
	conv_zOfPriVtxFromTracks.clear();
	conv_dzToClosestPriVtx.clear();
	conv_dxyPriVtx_Tr1.clear();
	conv_dxyPriVtx_Tr2.clear();
	conv_dxyPriVtxTimesCharge_Tr1.clear();
	conv_dxyPriVtxTimesCharge_Tr2.clear();
	conv_dxyError_Tr1.clear();
	conv_dxyError_Tr2.clear();

	conv_tk1NumOfDOF.clear();
	conv_tk2NumOfDOF.clear();
	conv_track1Chi2.clear();
	conv_track2Chi2.clear();
	conv_Tr1_pt.clear();
	conv_Tr2_pt.clear();
	conv_vertexChi2Prob.clear();
	conv_minDistanceOfApproach.clear();
	conv_p4->Clear(); 
	conv_eta.clear();
	conv_pt.clear();


	//conv MC
	conv_isMatchedMC.clear();
	convGen_eta.clear();
	convGen_pt.clear();
	convGen_p4->Clear();
	convGen_rDelta.clear();
	convGen_ptDelta.clear();
	convGen_ptDeltaRel.clear();
	convGen_motherCode.clear();

	//MC
	gen_pdgId.clear();
	gen_isGoodChicDecay.clear();
	gen_chic_pt.clear();
	gen_chic_eta.clear();
	gen_chic_p4->Clear();
	gen_chic_matchPosition.clear();
	gen_chic_nMatches.clear();
	gen_chic_rDelta.clear();
	gen_chic_ptDeltaRel.clear();
	gen_Jpsi_pt.clear();
	gen_Jpsi_eta.clear();
	gen_Jpsi_p4->Clear();
	gen_Jpsi_matchPosition.clear();
	gen_Jpsi_nMatches.clear();
	gen_Jpsi_rDelta.clear();
	gen_Jpsi_ptDeltaRel.clear();
	gen_Jpsi_photon_n.clear();
	gen_Jpsi_photon_pt.clear();
	gen_Jpsi_photon_p4->Clear();
	gen_muon_charge.clear();
	gen_muon_pt.clear();
	gen_muon_eta.clear();
	gen_muon_p4->Clear();
	gen_muon_matchPosition.clear();
	gen_muon_nMatches.clear();
	gen_muon_rDelta.clear();
	gen_muon_ptDeltaRel.clear();
	gen_phot_pt.clear();
	gen_phot_eta.clear();
	gen_phot_p4->Clear();
	gen_conv_matchPosition.clear();
	gen_conv_nMatches.clear();
	gen_conv_rDelta.clear();
	gen_conv_ptDeltaRel.clear();

	//chi
	chi_p4->Clear();
	chi_eta.clear();
	chi_pt.clear();
	chi_daughterJpsi_position.clear();
	chi_daughterConv_position.clear();
	chi_dzPhotToDimuonVtx.clear();
	chi_dxyPhotToDimuonVtx.clear();
	chiStored.clear();
	chi_kinematicRefitFlag.clear();
	chi_refit_origChicPosition.clear();
	chi_refitStored.clear();
	chi_refit_vprob.clear();
	chi_refit_ctauPV.clear();
	chi_refit_ctauErrPV.clear();
	chi_refit_ctauPV3D.clear();
	chi_refit_pvtxFromPVwithMuons_x.clear();
	chi_refit_pvtxFromPVwithMuons_y.clear();
	chi_refit_pvtxFromPVwithMuons_z.clear();

}


int ChiRootupler::SelectVertex(edm::Handle<std::vector<reco::Vertex> >& priVtxs, double zPos)
{
	float minDz = 9999;
	//int bestIdx = -1;
	int bestIdx = 0;
	for (uint i = 0; i < priVtxs->size(); i++) {
		const reco::Vertex& pvtx = priVtxs->at(i);
		if (i==0){
			if (fabs(zPos - pvtx.position().z()) < cPVMatching_cutoff) return 0; //if the difference is less than the value, don't do anything
		}
		if (pvtx.isFake() || fabs(pvtx.position().z()) > 50 || pvtx.position().Rho() > 2 || pvtx.tracksSize() < 2)
		{
			std::cout << "Vertex is bad, this shouldn't happen" << std::endl;
			continue;
		}
		float deltaZ = fabs(zPos - pvtx.position().z());
		if (deltaZ < minDz) {
			//std::cout << "deltaZ my code  " << deltaZ << " and the z pos "<< pvtx.position().z() << std::endl;
			minDz = deltaZ;
			bestIdx = i;
		}
	}
	return bestIdx;
}


//////////////////////////////////////////////
///// P H O T O N   C H E C K S    //////////  based on HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer
////////////////////////////////////////////

bool ChiRootupler::lt_comparator(std::pair<double, short> a, std::pair<double, short> b) {
	return a.first < b.first;
}

bool ChiRootupler::Conv_checkTkVtxCompatibility(const reco::Conversion& conv, const reco::VertexCollection& priVtxs, double sigmaTkVtxComp_, bool& Flag_Best_Out, bool& Flag_SecondBestA_Out, bool& Flag_SecondBestB_Out, double& sigmaMinValue1Out, double& sigmaMinValue2Out)
{
	sigmaMinValue1Out = 100; //starting values
	sigmaMinValue2Out = 100;
	Flag_Best_Out = false;
	Flag_SecondBestA_Out = false;
	Flag_SecondBestB_Out = false;

	std::vector< std::pair< double, short> > idx[2];  //dz distance, vertex number
	short ik = -1;
	BOOST_FOREACH(edm::RefToBase<reco::Track> tk, conv.tracks()) {  //ik - going over conversion tracks
		ik++;
		short count = -1;
		BOOST_FOREACH(const reco::Vertex& vtx, priVtxs) {  //count - going over vertices
			count++;
			double dz_ = tk->dz(vtx.position());
			double dzError_ = tk->dzError();
			dzError_ = sqrt(dzError_*dzError_ + vtx.covariance(2, 2));

			if ((ik==0)&&((fabs(dz_) / dzError_) < sigmaMinValue1Out)) {//save lowest value of sigma
				sigmaMinValue1Out = fabs(dz_) / dzError_;
			}
			if ((ik == 1) && ((fabs(dz_) / dzError_) < sigmaMinValue2Out)) {
				sigmaMinValue2Out = fabs(dz_) / dzError_;
			}

			if (fabs(dz_) / dzError_ > sigmaTkVtxComp_) continue;

			idx[ik].push_back(std::pair<double, short>(fabs(dz_), count));
		}
		if (idx[ik].size() == 0) { return false; }

		std::stable_sort(idx[ik].begin(), idx[ik].end(), lt_comparator);
	}

	Flag_Best_Out = (idx[0][0].second == idx[1][0].second);
	Flag_SecondBestA_Out = (idx[0][1].second == idx[1][0].second); //second best vertex is OK
	Flag_SecondBestB_Out = (idx[0][0].second == idx[1][1].second);

	if (idx[0][0].second == idx[1][0].second || idx[0][1].second == idx[1][0].second || idx[0][0].second == idx[1][1].second) return true;  //original return
	return false;
}

bool ChiRootupler::Conv_foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB) { //directly copied from HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer
	size_t count = 0;
	uint32_t oldSubStr = 0;
	for (int i = 0; i<hitPatA.numberOfHits(reco::HitPattern::HitCategory::TRACK_HITS) && count<2; i++) {
		uint32_t hitA = hitPatA.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS, i);
		if (!hitPatA.validHitFilter(hitA) || !hitPatA.trackerHitFilter(hitA)) continue;

		if (hitPatA.getSubStructure(hitA) == oldSubStr && hitPatA.getLayer(hitA) == oldSubStr)
			continue;

		if (hitPatB.getTrackerMonoStereo(reco::HitPattern::HitCategory::TRACK_HITS, hitPatA.getSubStructure(hitA), hitPatA.getLayer(hitA)) != 0)
			return true;

		oldSubStr = hitPatA.getSubStructure(hitA);
		count++;
	}
	return false;
}

//bool OniaPhotonConversionProducer::HighpuritySubset(const reco::Conversion& conv, const reco::VertexCollection& priVtxs) {	// select high purity conversions same way as OniaPhotonConversionProducer:
//	// vertex chi2 cut
//	if (ChiSquaredProbability(conv.conversionVertex().chi2(), conv.conversionVertex().ndof())< _vertexChi2ProbCut) return false;
//
//	// d0 cut
//	// Find closest primary vertex
//	int closest_pv_index = 0;
//	int i = 0;
//	BOOST_FOREACH(const reco::Vertex& vtx, priVtxs) {
//		if (conv.zOfPrimaryVertexFromTracks(vtx.position()) < conv.zOfPrimaryVertexFromTracks(priVtxs[closest_pv_index].position())) closest_pv_index = i;
//		i++;
//	}
//	// Now check impact parameter wtr with the just found closest primary vertex
//	BOOST_FOREACH(const edm::RefToBase<reco::Track> tk, conv.tracks()) if (-tk->dxy(priVtxs[closest_pv_index].position())*tk->charge() / tk->dxyError()<0) return false;
//
//	// chi2 of single tracks
//	BOOST_FOREACH(const edm::RefToBase<reco::Track> tk, conv.tracks()) if (tk->normalizedChi2() > _trackchi2Cut) return false;
//
//	// dof for each track  
//	BOOST_FOREACH(const edm::RefToBase<reco::Track> tk, conv.tracks()) if (tk->ndof()< TkMinNumOfDOF_) return false;
//
//	// distance of approach cut
//	if (conv.distOfMinimumApproach() < _minDistanceOfApproachMinCut || conv.distOfMinimumApproach() > _minDistanceOfApproachMaxCut) return false;
//
//	return true;
//}

bool ChiRootupler::Conv_isMatched(const math::XYZTLorentzVectorF& reco_conv, const reco::GenParticle& gen_phot, double maxDeltaR, double maxDPtRel)
{
	return reco::deltaR(reco_conv, gen_phot) < maxDeltaR && ((reco_conv.pt() - gen_phot.pt()) / (gen_phot.pt()+ 1E-9)) < maxDPtRel;
}



//////////////////////////////////
////  C H I   P A R T S   ////
//////////////////////////////


const pat::CompositeCandidate ChiRootupler::makeChiCandidate(const pat::CompositeCandidate& dimuon,	const pat::CompositeCandidate& photon) {
	pat::CompositeCandidate chiCand;
	chiCand.addDaughter(dimuon, "dimuon");
	chiCand.addDaughter(photon, "photon");
	chiCand.setVertex(dimuon.vertex());
	reco::Candidate::LorentzVector vChic = dimuon.p4() + photon.p4();
	chiCand.setP4(vChic);
	return chiCand;
}

double ChiRootupler::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {
	reco::Candidate::LorentzVector mom = c.p4();
	reco::Candidate::Point vtx = c.vertex();
	double dz = (vtx.Z() - p.Z()) - ((vtx.X() - p.X())*mom.X() + (vtx.Y() - p.Y())*mom.Y()) / mom.Rho() * mom.Z() / mom.Rho();
	return dz;
}

double ChiRootupler::Getdxy(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {
	reco::Candidate::LorentzVector mom = c.p4();
	reco::Candidate::Point vtx = c.vertex();
	double dx = (vtx.X() - p.X()) - ((vtx.Z() - p.Z()) * mom.X() / mom.Z());
	double dy = (vtx.Y() - p.Y()) - ((vtx.Z() - p.Z()) * mom.Y() / mom.Z());
	double dxy = sqrt(dx*dx+dy*dy);
	return dxy;
}

const pat::CompositeCandidate ChiRootupler::makeChiCandidate(const pat::CompositeCandidate& dimuon, const reco::Conversion& photon) {
	pat::CompositeCandidate chiCand;
	chiCand.addDaughter(dimuon, "dimuon");
	pat::CompositeCandidate photonComp;

	//chiCand.addDaughter(photon, "photon");
	chiCand.setVertex(dimuon.vertex());
	reco::Candidate::LorentzVector vChic = dimuon.p4() + photon.refittedPair4Momentum();  //checked, works
	chiCand.setP4(vChic);
	return chiCand;
}

double ChiRootupler::Getdz(const reco::Conversion& c, const reco::Candidate::Point &p) {
	math::XYZTLorentzVectorF mom = c.refittedPair4Momentum();
	reco::Candidate::Point vtx = c.conversionVertex().position();
	double dz = (vtx.Z() - p.Z()) - ((vtx.X() - p.X())*mom.X() + (vtx.Y() - p.Y())*mom.Y()) / mom.Rho() * mom.Z() / mom.Rho();
	return dz;
}

double ChiRootupler::Getdxy(const reco::Conversion& c, const reco::Candidate::Point &p) {
	math::XYZTLorentzVectorF mom = c.refittedPair4Momentum();
	reco::Candidate::Point vtx = c.conversionVertex().position();
	double dx = (vtx.X() - p.X()) - ((vtx.Z() - p.Z()) * mom.X() / mom.Z());
	double dy = (vtx.Y() - p.Y()) - ((vtx.Z() - p.Z()) * mom.Y() / mom.Z());
	double dxy = sqrt(dx*dx + dy * dy);
	return dxy;
}

bool ChiRootupler::ConvSelection(const reco::Conversion& conv)
{
	math::XYZTLorentzVectorF mom = conv.refittedPair4Momentum();
	if (mom.eta() > 2.5) return false;
	if (mom.pt() < 0.2) return false;
	return true;
}


///////////////////////////////////
////  G E N   P A R T S  //////////
//////////////////////////////////////

template <typename particle_in, typename particle_type>
int ChiRootupler::MatchGen(particle_in& myGenParticle, edm::Handle< std::vector <particle_type>>& collToBeMatched, double maxDeltaR, double maxDPtRel, int& nMatches_Out, double& rDelta_Out, double& ptDeltaRel_Out) //returns position of the best match
{
	nMatches_Out = 0;
	rDelta_Out = -5;
	ptDeltaRel_Out = -5;
	bool isMatched = false;
	int bestMatchPosition = -1;

	if (collToBeMatched.isValid())
	{
		for (uint j = 0; j < collToBeMatched->size(); j++) { //loop over the collection, find the best match
			const particle_type& myParticle = collToBeMatched->at(j);
			const reco::Candidate& recoParticle = GetRecoCandidate(myParticle); //get reco candidate of whatever type went in
			const reco::Candidate& genParticle = GetRecoCandidate(myGenParticle);
			double deltaR = reco::deltaR(recoParticle, genParticle);// reco::deltaR(inCand.eta(), inCand.phi(), mCand.eta(), mCand.phi());
			double dPtRel = std::fabs(recoParticle.pt() - genParticle.pt()) / (double)(genParticle.pt() + 1E-9); //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
			//std::cout << "pt reco "<< recoParticle.pt()<<" pt gen "<< genParticle.pt()<<"  and relative "<<dPtRel << "   difference "<< std::fabs(recoParticle.pt() - genParticle.pt())<<"  and denom "<< (genParticle.pt() + 1E-9) <<std::endl;
			
			if ((deltaR < maxDeltaR) && (dPtRel < maxDPtRel) && (recoParticle.charge() == genParticle.charge())) //this particular is matched
			{
				if (isMatched == false)//first match
				{
					isMatched = true;
					nMatches_Out++;
					bestMatchPosition = j;
					rDelta_Out = deltaR;
					ptDeltaRel_Out = dPtRel;
				}
				else //it has been matched already
				{
					nMatches_Out++;
					//compare with the best match by R, save the better one
					const reco::Candidate& prevMatch = GetRecoCandidate(collToBeMatched->at(bestMatchPosition));
					if (deltaR < reco::deltaR(prevMatch, genParticle)) {
						bestMatchPosition = j;
						rDelta_Out = deltaR;
						ptDeltaRel_Out = dPtRel;
					}
				}

			}

		}
		return bestMatchPosition;
	}
	else {
		std::cout << "Passed handle invalid" << std::endl;
		return -2; //-2 is a problem
	}
}

template <typename particle_in, typename particle_type>
int ChiRootupler::MatchGen(particle_in& myGenParticle, std::vector <particle_type>* collToBeMatched, double maxDeltaR, double maxDPtRel, int& nMatches_Out, double& rDelta_Out, double& ptDeltaRel_Out) //returns position of the best match
{
	nMatches_Out = 0;
	rDelta_Out = -5;
	ptDeltaRel_Out = -5;
	bool isMatched = false;
	int bestMatchPosition = -1;


	for (uint j = 0; j < collToBeMatched->size(); j++) { //loop over the collection, find the best match
		const particle_type& myParticle = collToBeMatched->at(j);
		const reco::Candidate& recoParticle = GetRecoCandidate(myParticle); //get reco candidate of whatever type went in
		const reco::Candidate& genParticle = GetRecoCandidate(myGenParticle);
		double deltaR = reco::deltaR(recoParticle, genParticle);// reco::deltaR(inCand.eta(), inCand.phi(), mCand.eta(), mCand.phi());
		double dPtRel = std::fabs(recoParticle.pt() - genParticle.pt()) / (double)(genParticle.pt() + 1E-9); //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
		//std::cout << "pt reco "<< recoParticle.pt()<<" pt gen "<< genParticle.pt()<<"  and relative "<<dPtRel << "   difference "<< std::fabs(recoParticle.pt() - genParticle.pt())<<"  and denom "<< (genParticle.pt() + 1E-9) <<std::endl;

		if ((deltaR < maxDeltaR) && (dPtRel < maxDPtRel) && (recoParticle.charge() == genParticle.charge())) //this particular is matched
		{
			if (isMatched == false)//first match
			{
				isMatched = true;
				nMatches_Out++;
				bestMatchPosition = j;
				rDelta_Out = deltaR;
				ptDeltaRel_Out = dPtRel;
			}
			else //it has been matched already
			{
				nMatches_Out++;
				//compare with the best match by R, save the better one
				const reco::Candidate& prevMatch = GetRecoCandidate(collToBeMatched->at(bestMatchPosition));
				if (deltaR < reco::deltaR(prevMatch, genParticle)) {
					bestMatchPosition = j;
					rDelta_Out = deltaR;
					ptDeltaRel_Out = dPtRel;
				}
			}

		}

	}
	return bestMatchPosition;

}

template<typename T>  reco::LeafCandidate ChiRootupler::GetRecoCandidate(const T& a) { return (reco::LeafCandidate) a; }

reco::LeafCandidate ChiRootupler::GetRecoCandidate(const reco::Conversion& conv) 
{
	reco::Candidate::LorentzVector p4;
	if (conv.conversionVertex().isValid()) {
		p4.SetPxPyPzE(conv.refittedPair4Momentum().Px(), conv.refittedPair4Momentum().Py(), conv.refittedPair4Momentum().Pz(), conv.refittedPair4Momentum().E());
	}
	else { 
		p4.SetPxPyPzE(conv.pairMomentum().X(), conv.pairMomentum().Y(), conv.pairMomentum().Z(), conv.pairMomentum().R()); //assumes no mass
	}
	reco::LeafCandidate out;
	out.setP4(p4);
	out.setCharge(0.0);
	return out;
}

reco::LeafCandidate ChiRootupler::GetRecoCandidate(const pat::CompositeCandidate& comp)
{
	reco::Candidate::LorentzVector p4 = comp.p4();
	reco::LeafCandidate out;
	out.setP4(p4);
	out.setCharge(comp.charge());
	return out;
}


reco::LeafCandidate ChiRootupler::GetRecoCandidate(const pat::Muon& comp)
{

	reco::Candidate::LorentzVector p4 = comp.p4();
	reco::LeafCandidate out;
	out.setP4(p4);
	out.setCharge(comp.charge());
	return out;
}

reco::LeafCandidate ChiRootupler::GetRecoCandidate(const reco::GenParticle& comp)
{
	reco::Candidate::LorentzVector p4 = comp.p4();
	reco::LeafCandidate out;
	out.setP4(p4);
	out.setCharge(comp.charge());
	return out;
}











////////////////////////////////////////////////

//define this as a plug-in
DEFINE_FWK_MODULE(ChiRootupler);
