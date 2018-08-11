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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/PatCandidates/interface/UserData.h> 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <boost/foreach.hpp>



//constructor
ChiRootupler::ChiRootupler(const edm::ParameterSet & iConfig) :
	muon_label(consumes<edm::View <pat::Muon> >(iConfig.getParameter < edm::InputTag >("muon_cand"))),
	dimuon_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("dimuon_cand"))),
	photon_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("photon_cand"))),
	conversion_label(consumes<reco::ConversionCollection>(iConfig.getParameter < edm::InputTag >("conversions_ch"))),
	//chi_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("chi_cand"))),
	primaryVertices_label(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag >("primaryVertices"))),
	triggerResults_label(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag >("TriggerResults"))),
	genParticles_label(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag >("genParticlesTag"))),
	flag_doMC(iConfig.getParameter < bool >("isMC"))
{

	edm::Service < TFileService > fs;

	//general information
	generalInfo_tree = fs->make < TTree >("generalInfo_tree", "General Information");
	generalInfo_tree->Branch("run", &runNumber, "run/L");
	generalInfo_tree->Branch("event", &eventNumber, "event/L");
	generalInfo_tree->Branch("nPrimVertices", &nPrimVertices, "nPrimVertices/L");
	generalInfo_tree->Branch("muonPerEvent_noCuts", &muonPerEvent, "muonPerEvent/I");
	generalInfo_tree->Branch("convPerEvent_noCuts", &convPerTriggeredEvent, "convPerEvent/I");

	// muon
	muon_tree = fs->make < TTree >("muonTree", "Muon tree");

	muon_tree->Branch("muonPerEvent_noCuts", &muonPerEvent, "muonPerEvent/I");
	muon_tree->Branch("muonIsGlobal", &muonIsGlobal, "muonIsGlobal/B");
	muon_tree->Branch("muonIsTracker", &muonIsTracker, "muonIsTracker/B");
	muon_tree->Branch("muonIsPF", &muonIsPF, "muonIsPF/B");
	muon_tree->Branch("muonIsNotGlobalNorTracker", &muonIsNotGlobalNorTracker, "muonIsNotGlobalNorTracker/B");
	muon_tree->Branch("muonIDHas_TMOneStationTight", &muonIDHas_TMOneStationTight, "muonIDHas_TMOneStationTight/B");
	muon_tree->Branch("muonInnerTrack_dxy", &muonInnerTrack_dxy, "muonInnerTrack_dxy/D");
	muon_tree->Branch("muonInnerTrack_dz", &muonInnerTrack_dz, "muonInnerTrack_dz/D");
	muon_tree->Branch("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement, "muonTrackerLayersWithMeasurement/I");
	muon_tree->Branch("muonPixelLayersWithMeasurement", &muonPixelLayersWithMeasurement, "muonpixelLayersWithMeasurement/I");
	muon_tree->Branch("muonQuality_isHighPurity", &muonQuality_isHighPurity, "muonQuality_isHighPurity/B");
	muon_tree->Branch("muon_charge", &muon_charge, "muon_charge/I");
	muon_tree->Branch("muon_eta", &muon_eta, "muon_eta/D");
	muon_tree->Branch("muon_pt", &muon_pt, "muon_pt/D");
	muon_tree->Branch("muon_p4", "TLorentzVector", &muon_p4);
	if (flag_doMC) {
		muon_tree->Branch("muon_isMatchedMC", &muon_isMatchedMC, "muon_isMatchedMC/B");
		muon_tree->Branch("muonGen_eta", &muonGen_eta, "muonGen_eta/D");
		muon_tree->Branch("muonGen_pt", &muonGen_pt, "muonGen_pt/D");
		muon_tree->Branch("muonGen_p4", "TLorentzVector", &muonGen_p4);
		muon_tree->Branch("muonGen_rDelta", &muonGen_rDelta, "muonGen_rDelta/D");
		muon_tree->Branch("muonGen_ptDelta", &muonGen_ptDelta, "muonGen_ptDelta/D");
		muon_tree->Branch("muonGen_ptDeltaRel", &muonGen_ptDeltaRel, "muonGen_ptDeltaRel/D");
	}
	if (flag_saveExtraThings) {
		muon_tree->Branch("patMuon", "pat::Muon", &patMuonStored);
	}

	// dimuon - TBD
	dimuon_tree = fs->make < TTree >("dimuonTree", "Tree of dimuon");
	dimuon_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
	dimuon_tree->Branch("dimuon_pt", &dimuon_pt);

	dimuon_tree->Branch("run", &runNumber, "run/L");
	dimuon_tree->Branch("event", &eventNumber, "event/L");
	dimuon_tree->Branch("primary_v", "TVector3", &primary_v);
	dimuon_tree->Branch("secondary_v", "TVector3", &secondary_v);
	dimuon_tree->Branch("dimuon_v", "TVector3", &dimuon_v);


	dimuon_tree->Branch("dimuonStored", "pat::CompositeCandidate", &dimuonStored);


	// conversions
	conv_tree = fs->make < TTree >("convTree", "Tree of conversions");
	
	conv_tree->Branch("convPerEvent_noCuts2", &convPerTriggeredEvent2);
	conv_tree->Branch("convQuality_isHighPurity", &convQuality_isHighPurity, "convQuality_isHighPurity/B");
	conv_tree->Branch("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly, "convQuality_isGeneralTracksOnly/B");
	conv_tree->Branch("conv_vertexPositionRho", &conv_vertexPositionRho, "conv_vertexPositionRho/D");
	conv_tree->Branch("conv_sigmaTkVtx1", &conv_sigmaTkVtx1, "conv_sigmaTkVtx1/D");
	conv_tree->Branch("conv_sigmaTkVtx2", &conv_sigmaTkVtx2, "conv_sigmaTkVtx2/D");
	conv_tree->Branch("conv_tkVtxCompatibilityOK50", &conv_tkVtxCompatibilityOK50, "conv_tkVtxCompatibilityOK50/B");
	conv_tree->Branch("conv_tkVtxCompatible_bestVertex50", &conv_tkVtxCompatible_bestVertex50, "conv_tkVtxCompatible_bestVertex50/B");
	conv_tree->Branch("conv_tkVtxCompatible_secondBestVertexA50", &conv_tkVtxCompatible_secondBestVertexA50, "conv_tkVtxCompatible_secondBestVertexA50/B");
	conv_tree->Branch("conv_tkVtxCompatible_secondBestVertexB50", &conv_tkVtxCompatible_secondBestVertexB50, "conv_tkVtxCompatible_secondBestVertexB50/B");
	conv_tree->Branch("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK, "conv_compatibleInnerHitsOK/I");
	if (flag_saveExtraThings) {
		conv_tree->Branch("conv_hitPat1", "reco::HitPattern", &conv_hitPat1);
		conv_tree->Branch("conv_hitPat2", "reco::HitPattern", &conv_hitPat2);
		conv_tree->Branch("conv_tkVtxCompatibilityOK3", &conv_tkVtxCompatibilityOK3, "conv_tkVtxCompatibilityOK3/B");
		conv_tree->Branch("conv_tkVtxCompatible_bestVertex3", &conv_tkVtxCompatible_bestVertex3, "conv_tkVtxCompatible_bestVertex3/B");
		conv_tree->Branch("conv_tkVtxCompatible_secondBestVertexA3", &conv_tkVtxCompatible_secondBestVertexA3, "conv_tkVtxCompatible_secondBestVertexA3/B");
		conv_tree->Branch("conv_tkVtxCompatible_secondBestVertexB3", &conv_tkVtxCompatible_secondBestVertexB3, "conv_tkVtxCompatible_secondBestVertexB3/B");
	}
	conv_tree->Branch("conv_vertexChi2Prob", &conv_vertexChi2Prob, "conv_vertexChi2Prob/D");
	conv_tree->Branch("conv_zOfPriVtx", &conv_zOfPriVtx, "conv_zOfPriVtx/D");
	conv_tree->Branch("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks, "conv_zOfPriVtxFromTracks/D");
	conv_tree->Branch("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx, "conv_dzToClosestPriVtx/D");
	conv_tree->Branch("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1, "conv_dxyPriVtx_Tr1/D");
	conv_tree->Branch("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2, "conv_dxyPriVtx_Tr2/D");
	conv_tree->Branch("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1, "conv_dxyPriVtxTimesCharge_Tr1/D");
	conv_tree->Branch("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2, "conv_dxyPriVtxTimesCharge_Tr2/D");
	conv_tree->Branch("conv_dxyError_Tr1", &conv_dxyError_Tr1, "conv_dxyError_Tr1/D");
	conv_tree->Branch("conv_dxyError_Tr2", &conv_dxyError_Tr2, "conv_dxyError_Tr2/D");

	conv_tree->Branch("conv_tk1NumOfDOF", &conv_tk1NumOfDOF, "conv_tk1NumOfDOF/I");
	conv_tree->Branch("conv_tk2NumOfDOF", &conv_tk2NumOfDOF, "conv_tk2NumOfDOF/I");
	conv_tree->Branch("conv_track1Chi2", &conv_track1Chi2, "conv_track1Chi2/D");
	conv_tree->Branch("conv_track2Chi2", &conv_track2Chi2, "conv_track2Chi2/D");
	conv_tree->Branch("conv_minDistanceOfApproach", &conv_minDistanceOfApproach, "conv_minDistanceOfApproach/D");

	conv_tree->Branch("conv_p4", "TLorentzVector", &conv_p4);
	conv_tree->Branch("conv_eta", &conv_eta, "conv_eta/D");
	conv_tree->Branch("conv_pt", &conv_pt, "conv_pt/D");
	conv_tree->Branch("conv_pt2", &conv_pt2);



	if (flag_doMC) {
		conv_tree->Branch("conv_isMatchedMC", &conv_isMatchedMC, "conv_isMatchedMC/B");
		conv_tree->Branch("convGen_eta", &convGen_eta, "convGen_eta/D");
		conv_tree->Branch("convGen_pt", &convGen_pt, "convGen_pt/D");
		conv_tree->Branch("convGen_p4", "TLorentzVector", &convGen_p4);
		conv_tree->Branch("convGen_rDelta", &convGen_rDelta, "convGen_rDelta/D");
		conv_tree->Branch("convGen_ptDelta", &convGen_ptDelta, "convGen_ptDelta/D");
		conv_tree->Branch("convGen_ptDeltaRel", &convGen_ptDeltaRel, "convGen_ptDeltaRel/D");
		conv_tree->Branch("convGen_motherCode", &convGen_motherCode, "convGen_motherCode/I");
	}


	// MC general

	if (flag_doMC) {
		gen_tree = fs->make < TTree >("genTree", "MC general Tree");
		gen_tree->Branch("gen_pdgId", &gen_pdgId, "gen_pdgId/I");

		gen_tree->Branch("gen_chic_pt", &gen_chic_pt, "gen_chic_pt/D");
		gen_tree->Branch("gen_chic_eta", &gen_chic_eta, "gen_chic_eta/D");
		gen_tree->Branch("gen_chic_p4", "TLorentzVector", &gen_chic_p4);

		gen_tree->Branch("gen_Jpsi_pt", &gen_Jpsi_pt, "gen_Jpsi_pt/D");
		gen_tree->Branch("gen_Jpsi_eta", &gen_Jpsi_eta, "gen_Jpsi_eta/D");
		gen_tree->Branch("gen_Jpsi_p4", "TLorentzVector", &gen_Jpsi_p4);

		gen_tree->Branch("gen_phot_pt", &gen_phot_pt, "gen_phot_pt/D");
		gen_tree->Branch("gen_phot_eta", &gen_phot_eta, "gen_phot_eta/D");
		gen_tree->Branch("gen_phot_p4", "TLorentzVector", &gen_phot_p4);

	}

	// chi - TBD
	chi_tree = fs->make < TTree >("chiTree", "Tree of chi");

	chi_tree->Branch("run", &runNumber, "run/L");
	chi_tree->Branch("event", &eventNumber, "event/L");
	chi_tree->Branch("chi_cand", "pat::CompositeCandidate", &chi_cand);
	chi_tree->Branch("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx, "chi_dzPhotToDimuonVtx/D");

	chi_tree->Branch("primary_v", "TVector3", &primary_v);
	chi_tree->Branch("secondary_v", "TVector3", &secondary_v);
	chi_tree->Branch("dimuon_v", "TVector3", &dimuon_v);
	chi_tree->Branch("conversion_v", "TVector3", &conversion_v);

	chi_tree->Branch("chi_p4", "TLorentzVector", &chi_p4);
	chi_tree->Branch("chi_dimuon_p4", "TLorentzVector", &chi_dimuon_p4);
	chi_tree->Branch("muonP_p4", "TLorentzVector", &muonP_p4);
	chi_tree->Branch("muonN_p4", "TLorentzVector", &muonN_p4);
	chi_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);

	chi_tree->Branch("ele_lowerPt_pt", &ele_lowerPt_pt, "ele_lowerPt_pt/D");
	chi_tree->Branch("ele_higherPt_pt", &ele_higherPt_pt, "ele_higherPt_pt/D");
	chi_tree->Branch("ctpv", &ctpv, "ctpv/D");
	chi_tree->Branch("ctpv_error", &ctpv_error, "ctpv_error/D");
	chi_tree->Branch("pi0_abs_mass", &pi0_abs_mass, "pi0_abs_mass/D");
	chi_tree->Branch("psi1S_nsigma", &psi1S_nsigma, "psi1S_nsigma/D");
	chi_tree->Branch("psi2S_nsigma", &psi2S_nsigma, "psi2S_nsigma/D");

	chi_tree->Branch("conv_vertex", &conv_vertex, "conv_vertex/D");
	chi_tree->Branch("trigger", &trigger, "trigger/I");

	if (flag_doMC) {
		chi_tree->Branch("gen_chic_p4", "TLorentzVector", &gen_chic_p4);
		chi_tree->Branch("chic_pdgId", &chic_pdgId, "chic_pdgId/I");
		chi_tree->Branch("gen_jpsi_p4", "TLorentzVector", &gen_jpsi_p4);
		chi_tree->Branch("gen_photon_p4", "TLorentzVector", &gen_photon_p4);
		chi_tree->Branch("gen_muonP_p4", "TLorentzVector", &gen_muonP_p4);
		chi_tree->Branch("gen_muonM_p4", "TLorentzVector", &gen_muonM_p4);
	}
	genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
	packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");

}

ChiRootupler::~ChiRootupler() {}

//
// member functions
//

////Check recursively if any ancestor of particle is the given one
//bool ChiRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
//   if (ancestor == particle ) return true;
//   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
//   return false;
//}

// ------------ method called for each event  ------------
void ChiRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
	using namespace edm;
	using namespace std;

	edm::Handle  < edm::View <pat::Muon> > muon_handle;
	iEvent.getByToken(muon_label, muon_handle);

	edm::Handle < std::vector < pat::CompositeCandidate > >dimuon_handle;
	iEvent.getByToken(dimuon_label, dimuon_handle);

	edm::Handle < std::vector < pat::CompositeCandidate > >photon_handle;
	iEvent.getByToken(photon_label, photon_handle);

	edm::Handle< std::vector <reco::Conversion>> conversion_handle;
	iEvent.getByToken(conversion_label, conversion_handle);

	//edm::Handle < std::vector < pat::CompositeCandidate >>chi_handle;
	//iEvent.getByToken(chi_label, chi_handle);

	edm::Handle  < reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_label, primaryVertices_handle);

	edm::Handle < edm::TriggerResults > triggerResults_handle;
	iEvent.getByToken(triggerResults_label, triggerResults_handle);

	edm::Handle <reco::GenParticleCollection> genParticles_handle;
	iEvent.getByToken(genParticles_label, genParticles_handle);

	

	//edm::Handle<reco::VertexCollection> priVtxs;
	//event.getByToken(thePVsToken_, priVtxs);



	//edm::Handle<reco::PFCandidateCollection> pfcandidates;
	//event.getByToken(pfCandidateCollectionToken_, pfcandidates);

	//const reco::PFCandidateCollection pfphotons = selectPFPhotons(*pfcandidates);






	/////////////////////
	//// S T A R T //////
	////////////////////

	//general info
	runNumber = iEvent.id().run();
	eventNumber = iEvent.id().event();
	nPrimVertices = 0;
	if (primaryVertices_handle.isValid()) {
		nPrimVertices = primaryVertices_handle->size();
	}
	muonPerEvent = 0;
	if (muon_handle.isValid()) {
		muonPerEvent = muon_handle->size(); //all muons without any cuts
	}
	convPerTriggeredEvent = 0;
	if (conversion_handle.isValid())
	{
		convPerTriggeredEvent = conversion_handle->size(); //all conversions without any cuts
	}

	generalInfo_tree->Fill();

	//muons
	if (muon_handle.isValid()) {
		for (uint i = 0; i < muon_handle->size(); i++) {
			const pat::Muon& patMuon = muon_handle->at(i);

			muonIsGlobal = patMuon.isGlobalMuon();
			muonIsTracker = patMuon.isTrackerMuon();
			muonIsPF = patMuon.isPFMuon();
			if (!patMuon.isGlobalMuon() && !patMuon.isTrackerMuon()) { muonIsNotGlobalNorTracker = true; }
			else muonIsNotGlobalNorTracker = false; // just for convenience

			muonIDHas_TMOneStationTight = patMuon.muonID("TMOneStationTight");
			if (patMuon.isGlobalMuon() || patMuon.isTrackerMuon()) {
				muonInnerTrack_dxy = patMuon.innerTrack()->dxy();
				muonInnerTrack_dz = patMuon.innerTrack()->dz();
				muonTrackerLayersWithMeasurement = patMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement();
				muonPixelLayersWithMeasurement = patMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement();
				reco::TrackBase::TrackQuality tq = reco::TrackBase::qualityByName("highPurity");//high purity=2 //see DataFormats/TrackReco/interface/TrackBase.h
				muonQuality_isHighPurity = patMuon.innerTrack()->quality(tq);
			}
			muon_charge = patMuon.charge();
			muon_eta = patMuon.eta();
			muon_pt = patMuon.pt();
			muon_p4.SetPtEtaPhiM(patMuon.pt(), patMuon.eta(), patMuon.phi(), patMuon.mass());

			// MC generated information
			if (flag_doMC) {
				muon_isMatchedMC = patMuon.genLepton(); //false for null pointer, true if match exists, from pat::Lepton.h
				if (muon_isMatchedMC) {
					const reco::GenParticle genMuon = *patMuon.genLepton();
					muonGen_eta = genMuon.eta();
					muonGen_pt = genMuon.pt();
					muonGen_p4.SetPtEtaPhiM(genMuon.pt(), genMuon.eta(), genMuon.phi(), genMuon.mass());
					muonGen_rDelta = reco::deltaR(patMuon, genMuon); //sqrt(phi^2+eta^2)
					muonGen_ptDelta = patMuon.pt() - genMuon.pt();
					muonGen_ptDeltaRel = muonGen_ptDelta / muonGen_pt; //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
				}
				else { //default values to store if no match - in principle can be ommitted, just for direct looking at branches
					muonGen_eta = -5;
					muonGen_pt = 0;
					muonGen_p4 = TLorentzVector();
					muonGen_rDelta = -5;
					muonGen_ptDelta = -5;
					muonGen_ptDeltaRel = -5;
				}
			}


			if (flag_saveExtraThings)
			{
				patMuonStored = patMuon;
			}

			muon_tree->Fill();
		}
	}
	else cout << "Problem with muon handle" << endl;

	//
	//if (genParticles_handle.isValid()) {
	//	//cout << "Gen: "<<genParticles_handle->size() << endl;
	//	/*for (uint i = 0; i < genParticles_handle->size(); i++) {
	//		const reco::GenParticle& genParticle = genParticles_handle->at(i);// genParticle(size_t idx = 0) const;
	//		if (i == 100) { cout << "GEN pT " << genParticle.pt() << endl; }
	//	}*/
	//}
	//else cout << "Problem with gen handle" << endl;


	if (dimuon_handle.isValid()) //O: add if store dimuons
	{
		for (uint i = 0; i < dimuon_handle->size(); i++) {
			const pat::CompositeCandidate& dimuon = dimuon_handle->at(i);
			dimuon_p4.SetPtEtaPhiM(dimuon.pt(), dimuon.eta(), dimuon.phi(), dimuon.mass());
			dimuon_pt = dimuon.pt();
			//dimuon.userData<reco::Vertex>("PVwithmuons");
			//const reco::Vertex* opv = dimuon.userData<reco::Vertex>("PVwithmuons");
			const reco::Vertex* opv = dimuon.userData<reco::Vertex>("commonVertex");
			//const reco::Vertex* opv = (dynamic_cast < pat::CompositeCandidate*> (&dimuon))->userData<reco::Vertex>("PVwithmuons");
			//cout << opv << endl;
			//primary_v.SetXYZ(opv->x(), opv->y(), opv->z());
			cout << eventNumber << "   " << i << " " << dimuon.userFloat("DCA")<< endl;
			//cout << opv->x() << endl;
			dimuonStored = dimuon;
			dimuon_tree->Fill();
		}
	}

	// conversions

	if (conversion_handle.isValid())
	{
		for (uint i = 0; i < conversion_handle->size(); i++) {
			const reco::Conversion& candPhoton = conversion_handle->at(i);
			conv_p4.SetXYZT(candPhoton.refittedPair4Momentum().x(), candPhoton.refittedPair4Momentum().y(), candPhoton.refittedPair4Momentum().z(), candPhoton.refittedPair4Momentum().t());
			conv_eta = conv_p4.Eta();
			conv_pt = conv_p4.Pt();

			convQuality_isHighPurity = candPhoton.quality((reco::Conversion::ConversionQuality)(8)); //8 is high purity, see reco::Conversion Class Reference
			convQuality_isGeneralTracksOnly = candPhoton.quality((reco::Conversion::ConversionQuality)(0)); //0 is general tracks only, see reco::Conversion Class Reference
			conv_vertexPositionRho = candPhoton.conversionVertex().position().rho();
			conv_tkVtxCompatibilityOK3 = Conv_checkTkVtxCompatibility(candPhoton, *primaryVertices_handle.product(), 20, conv_tkVtxCompatible_bestVertex3, conv_tkVtxCompatible_secondBestVertexA3, conv_tkVtxCompatible_secondBestVertexB3, conv_sigmaTkVtx1, conv_sigmaTkVtx2);
			conv_tkVtxCompatibilityOK50 = Conv_checkTkVtxCompatibility(candPhoton, *primaryVertices_handle.product(), conv_TkVtxCompSigmaCut, conv_tkVtxCompatible_bestVertex50, conv_tkVtxCompatible_secondBestVertexA50, conv_tkVtxCompatible_secondBestVertexB50, conv_sigmaTkVtx1, conv_sigmaTkVtx2);
			if (candPhoton.tracks().size() == 2) {
				const edm::RefToBase<reco::Track> conv_tk1 = candPhoton.tracks().at(0);
				const edm::RefToBase<reco::Track> conv_tk2 = candPhoton.tracks().at(1);

				reco::HitPattern hitPatA = conv_tk1->hitPattern();
				reco::HitPattern hitPatB = conv_tk2->hitPattern();
				conv_hitPat1 = hitPatA;
				conv_hitPat2 = hitPatB;
				conv_compatibleInnerHitsOK = (Conv_foundCompatibleInnerHits(hitPatA, hitPatB) && Conv_foundCompatibleInnerHits(hitPatB, hitPatA));

				//find vertex that points closest
				int closest_pv_index = 0;
				int i = 0;
				BOOST_FOREACH(const reco::Vertex& vtx, *primaryVertices_handle.product()) {
					if (fabs(candPhoton.zOfPrimaryVertexFromTracks(vtx.position()) - vtx.z()) < fabs(candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[closest_pv_index].position()) - (*primaryVertices_handle.product())[closest_pv_index].z())) { closest_pv_index = i; }
					i++;
				}
				conv_zOfPriVtx = (*primaryVertices_handle.product())[closest_pv_index].z();
				conv_zOfPriVtxFromTracks = candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[closest_pv_index].position());
				conv_dzToClosestPriVtx = candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[closest_pv_index].position()) - (*primaryVertices_handle.product())[closest_pv_index].z();
				// Now check impact parameter wtr with the just found closest primary vertex
				conv_dxyPriVtx_Tr1 = conv_tk1->dxy((*primaryVertices_handle.product())[closest_pv_index].position());
				conv_dxyPriVtx_Tr2 = conv_tk2->dxy((*primaryVertices_handle.product())[closest_pv_index].position());
				conv_dxyPriVtxTimesCharge_Tr1 = conv_tk1->dxy((*primaryVertices_handle.product())[closest_pv_index].position())*conv_tk1->charge();
				conv_dxyPriVtxTimesCharge_Tr2 = conv_tk2->dxy((*primaryVertices_handle.product())[closest_pv_index].position())*conv_tk2->charge();
				conv_dxyError_Tr1 = conv_tk1->dxyError();
				conv_dxyError_Tr2 = conv_tk2->dxyError();

				conv_tk1NumOfDOF = conv_tk1->ndof();
				conv_tk2NumOfDOF = conv_tk2->ndof();
				conv_track1Chi2 = candPhoton.tracks().at(0)->normalizedChi2();
				conv_track2Chi2 = candPhoton.tracks().at(1)->normalizedChi2();

			}
			else conv_compatibleInnerHitsOK = -1;

			conv_vertexChi2Prob = ChiSquaredProbability(candPhoton.conversionVertex().chi2(), candPhoton.conversionVertex().ndof());
			conv_minDistanceOfApproach = candPhoton.distOfMinimumApproach();
			//if (candPhoton.distOfMinimumApproach() > -10 && candPhoton.distOfMinimumApproach() < 10) { conv_minDistanceOfApproach = candPhoton.distOfMinimumApproach(); }
			//else conv_minDistanceOfApproach = 0;

			//MC for conversions
			if (flag_doMC)
			{
				if (genParticles_handle.isValid()) {
					reco::GenParticle genConv_best = reco::GenParticle();
					conv_isMatchedMC = false;
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
							if (conv_isMatchedMC == false) {//first one found
								conv_isMatchedMC = true;
								genConv_best = genParticle;
							}
							else { //check whether the second match is better by deltaR than the first, save the better one
								if (reco::deltaR(candPhoton.refittedPair4Momentum(), genConv_best) > reco::deltaR(candPhoton.refittedPair4Momentum(), genParticle)) { genConv_best = genParticle; }
							}
						}
					}
					if (conv_isMatchedMC) {
						convGen_eta = genConv_best.eta();
						convGen_pt = genConv_best.pt();
						convGen_p4.SetPtEtaPhiM(genConv_best.pt(), genConv_best.eta(), genConv_best.phi(), genConv_best.mass());
						convGen_rDelta = reco::deltaR(candPhoton.refittedPair4Momentum(), genConv_best); //sqrt(phi^2+eta^2)
						convGen_ptDelta = candPhoton.refittedPair4Momentum().pt() - genConv_best.pt();
						convGen_ptDeltaRel = convGen_ptDelta / convGen_pt; //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
						convGen_motherCode = genConv_best.mother()->pdgId();
					}
					else { //default values to store if no match - in principle can be ommitted, just for direct looking at branches
						convGen_eta = -5;
						convGen_pt = 0;
						convGen_p4 = TLorentzVector();
						convGen_rDelta = -5;
						convGen_ptDelta = -5;
						convGen_ptDeltaRel = -5;
						convGen_motherCode = -1;
					}
				}
				else cout << "Problem with gen handle" << endl;
			}
			conv_tree->Fill();
		}
	}
	else { cout << "Conversions handle problem" << endl; }

//if (conversion_handle.isValid())
//{
//	convPerTriggeredEvent2.push_back(conversion_handle->size());
//	for (uint i = 0; i < conversion_handle->size(); i++) {
//		const reco::Conversion& candPhoton = conversion_handle->at(i);
//		conv_pt2.push_back(candPhoton.refittedPair4Momentum().y());
//
//	}
//	conv_tree->Fill();
//	convPerTriggeredEvent2.clear();
//	conv_pt2.clear();
//}
//else { cout << "Conversions handle problem" << endl; }



	// GEN 
	if (flag_doMC) {
		if (genParticles_handle.isValid()) {
			//int nChicPerEvent=0;
			for (uint i = 0; i < genParticles_handle->size(); i++) {
				gen_Jpsi_pt = -1;
				gen_phot_pt = -1;
				const reco::GenParticle& genParticle = genParticles_handle->at(i);
				gen_pdgId = genParticle.pdgId();
				if (gen_pdgId != PythCode_chic0 && gen_pdgId != PythCode_chic1 && gen_pdgId != PythCode_chic2) { continue; } //if not chi, skip
				if (genParticle.isLastCopy() == false) { continue; } // if not last copy of a given particle, skip
				//nChicPerEvent++;
				gen_chic_pt = genParticle.pt();
				gen_chic_eta = genParticle.eta();
				gen_chic_p4.SetXYZT(genParticle.p4().x(), genParticle.p4().y(), genParticle.p4().z(), genParticle.p4().t());

				int nDaughters = genParticle.numberOfDaughters();
				if (nDaughters != 2) { cout << "weird gen decay" << endl; } //all of them should be J/psi + gamma
				for (int j = 0; j < nDaughters; j++) {
					const reco::Candidate * gen_chiDaughter = genParticle.daughter(j);
					int dauId = gen_chiDaughter->pdgId();
					if (dauId == 443) {//Jpsi
						gen_Jpsi_pt = gen_chiDaughter->pt();
						gen_Jpsi_eta = gen_chiDaughter->eta();
						gen_Jpsi_p4.SetXYZT(gen_chiDaughter->p4().x(), gen_chiDaughter->p4().y(), gen_chiDaughter->p4().z(), gen_chiDaughter->p4().t());
					}
					else if (dauId == 22) {//photon
						gen_phot_pt = gen_chiDaughter->pt();
						gen_phot_eta = gen_chiDaughter->eta();
						gen_phot_p4.SetXYZT(gen_chiDaughter->p4().x(), gen_chiDaughter->p4().y(), gen_chiDaughter->p4().z(), gen_chiDaughter->p4().t());
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
				gen_tree->Fill();
			}
		}
		else cout << "Problem with gen handle" << endl;
	}




	//if (flag_doMC) {
	//    // Pruned particles are the one containing "important" stuff
	//    edm::Handle<std::vector<reco::GenParticle> > pruned;
	//    iEvent.getByToken(genCands_, pruned);
	//    // Packed particles are all the status 1, so usable to remake jets
	//    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
	//    edm::Handle<std::vector<pat::PackedGenParticle> > packed;
	//    iEvent.getByToken(packCands_, packed);

	//    //let's try to find all status1 originating directly from a Chi_c1 meson decay

	//    chic_pdgId = 0;
	//    int foundit = 0;
	//    for (size_t i = 0; i < pruned->size(); i++) {
	//	    int p_id = abs((*pruned)[i].pdgId());
	//	    if (p_id == 20443 || p_id == 445 || p_id == 10443) {
	//	  	  chic_pdgId = p_id;
	//	  	  foundit++;
	//	  	  const reco::Candidate * chic = &(*pruned)[i];
	//	  	  MC_chic_p4.SetPtEtaPhiM(chic->pt(), chic->eta(), chic->phi(), chic->mass());

	//	  	  const reco::Candidate * j1 = nullptr;
	//	  	  const reco::Candidate * j2 = nullptr;

	//	  	  for (size_t j = 0; j < packed->size(); j++) {
	//	  		  //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
	//	  		  const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
	//	  		  if (motherInPrunedCollection != nullptr && isAncestor(chic, motherInPrunedCollection)) {
	//	  			  const reco::Candidate * d = &(*packed)[j];
	//	  			  int dauId = d->pdgId();

	//	  			  if (dauId == 13 && motherInPrunedCollection->pdgId() == 443) {
	//	  				  gen_muonM_p4.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	//	  				  foundit++;
	//	  				  j1 = motherInPrunedCollection;
	//	  			  }
	//	  			  if (dauId == -13 && motherInPrunedCollection->pdgId() == 443) {
	//	  				  gen_muonP_p4.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	//	  				  foundit++;
	//	  				  j2 = motherInPrunedCollection;
	//	  			  }
	//	  			  if (dauId == 22) {
	//	  				  gen_photon_p4.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	//	  				  foundit++;
	//	  			  }
	//	  		  }
	//	  		  if (foundit == 4) break;
	//	  	  }
	//	  	  if (foundit == 4) {
	//	  		  if (j1 != nullptr && j2 != nullptr && j1 == j2) {
	//	  			  gen_jpsi_p4.SetPtEtaPhiM(j1->pt(), j1->eta(), j1->phi(), j1->mass());
	//	  			  break;
	//	  		  }
	//	  		  else {
	//	  			  std::cout << "Mother of muons does not match (" << j1->pdgId() << "," << j2->pdgId() << ")" << std::endl;
	//	  			  foundit = 0;
	//	  			  chic_pdgId = 0;
	//	  		  }
	//	  	  }
	//	  	  else {
	//	  		  std::cout << "Found just " << foundit << " out of 4 particles" << std::endl;
	//	  		  foundit = 0;
	//	  		  chic_pdgId = 0;
	//	  	  }
	//	    }  // if ( p_id
	//    } // for (size

	//    if (!chic_pdgId) { // sanity check
	//	    std::cout << "Rootupler does not found the given decay " <<
	//	  	  iEvent.id().run() << "," << iEvent.id().event() << std::endl;
	//    }
	//}


	//grab Trigger informations
	// save it in variable trigger, trigger is an int between 0 and 15, in binary it is:
	// (pass 11)(pass 8)(pass 7)(pass 5)
	// es. 11 = pass 5, 7 and 11
	// es. 4 = pass only 8

	trigger = 0;
	if (triggerResults_handle.isValid()) {

		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);

		unsigned int NTRIGGERS = 9;
		//  string TriggersToTest[NTRIGGERS] = {"HLT_Dimuon10_Jpsi_Barrel","HLT_Dimuon16_Jpsi","HLT_Dimuon20_Jpsi","HLT_Dimuon8_Upsilon_Barrel","HLT_Dimuon13_Upsilon","HLT_Dimuon8_PsiPrime_Barrel","HLT_Dimuon13_PsiPrime","HLT_Mu16_TkMu0_dEta18_Onia","HLT_Mu25_TkMu0_dEta18_Onia"};
		string TriggersToTest[NTRIGGERS] = { "HLT_PAL1DoubleMuOpen","HLT_Dimuon16_Jpsi","HLT_PAL1DoubleMuOpen_OS","HLT_Dimuon8_Upsilon_Barrel","HLT_Dimuon13_Upsilon","HLT_PAL1DoubleMuOpen_SS","HLT_Dimuon13_PsiPrime","HLT_Mu16_TkMu0_dEta18_Onia","HLT_Mu25_TkMu0_dEta18_Onia" };

		for (unsigned int i = 0; i < NTRIGGERS; i++) {
			for (int version = 1; version < 5; version++) {
				stringstream ss;
				ss << TriggersToTest[i] << "_v" << version;
				unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label().c_str());
				if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					trigger += (1 << i);
					break;
				}
			}
		}

	}
	else {
		std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
	} // if (trigger...

	//////////////////////
	//    C H I       ////
	/////////////////////

	pat::CompositeCandidateCollection* chiCandColl=new pat::CompositeCandidateCollection;
	// Note: since Dimuon cand are sorted by decreasing vertex probability then the first chi cand is the one associated with the "best" dimuon 
	for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon_handle->begin(); dimuonCand != dimuon_handle->end(); ++dimuonCand) 
	{
		// use only trigger-matched Jpsi or Upsilon if so requested 
		//if (triggerMatch_) {
		//if (!dimuonCand->userInt("isTriggerMatched")) continue;
		//}

		// loop on conversion candidates, make chi cand
		for (pat::CompositeCandidateCollection::const_iterator photCand = photon_handle->begin(); photCand != photon_handle->end(); ++photCand) {

			pat::CompositeCandidate chiCand = makeChiCandidate(*dimuonCand, *photCand);

			//if (!cutDeltaMass(chiCand, *dimuonCand)) {
			//	delta_mass_fail++;
			//	continue;
			//}
			chi_dzPhotToDimuonVtx = fabs(Getdz(*photCand, dimuonCand->vertex())); // double dz = fabs( photCand->dz(dimuonCand->vertex()) );
			

			//if (!cutdz(dz)) {
			//	dz_cut_fail++;
			//	continue;
			//}


			chiCandColl->push_back(chiCand);
		}
	}
		
	

	if (chiCandColl) {

		unsigned int csize = chiCandColl->size();
		for (unsigned int i = 0; i < csize; i++) {
			chi_cand = chiCandColl->at(i);
			chi_p4.SetPtEtaPhiM(chi_cand.pt(), chi_cand.eta(), chi_cand.phi(), chi_cand.mass());
			chi_dimuon_p4.SetPtEtaPhiM(chi_cand.daughter("dimuon")->pt(), chi_cand.daughter("dimuon")->eta(),
				chi_cand.daughter("dimuon")->phi(), chi_cand.daughter("dimuon")->mass());
			photon_p4.SetPtEtaPhiM(chi_cand.daughter("photon")->pt(), chi_cand.daughter("photon")->eta(),
				chi_cand.daughter("photon")->phi(), chi_cand.daughter("photon")->mass());
			reco::Candidate::LorentzVector vP = chi_cand.daughter("dimuon")->daughter("muon1")->p4();
			reco::Candidate::LorentzVector vM = chi_cand.daughter("dimuon")->daughter("muon2")->p4();

			if (chi_cand.daughter("dimuon")->daughter("muon1")->charge() < 0) {
				vP = chi_cand.daughter("dimuon")->daughter("muon2")->p4();
				vM = chi_cand.daughter("dimuon")->daughter("muon1")->p4();
			}

			muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
			muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

			Double_t ele1_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->userData<reco::Track>("track0"))->pt();
			Double_t ele2_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->userData<reco::Track>("track1"))->pt();

			if (ele1_pt > ele2_pt) {
				ele_higherPt_pt = ele1_pt;
				ele_lowerPt_pt = ele2_pt;
			}
			else {
				ele_higherPt_pt = ele2_pt;
				ele_lowerPt_pt = ele1_pt;
			}

			//float ctpv = (dynamic_cast <pat::CompositeCandidate *>(chi_cand.daughter("dimuon")))->userFloat("ppdlPV");
			//float ctpv_error = (dynamic_cast <pat::CompositeCandidate *>(chi_cand.daughter("dimuon")))->userFloat("ppdlErrPV");
			pi0_abs_mass = 1.;//pi0_abs_values[0];

			reco::Candidate::Point vtx = chi_cand.daughter("photon")->vertex();
			conversion_v.SetXYZ(vtx.X(), vtx.Y(), vtx.Z());
			conv_vertex = vtx.rho(); //chi_cand.daughter("photon")->vertex().rho();

			reco::Candidate::Point dimuon_vtx = chi_cand.daughter("dimuon")->vertex();
			dimuon_v.SetXYZ(dimuon_vtx.X(), dimuon_vtx.Y(), dimuon_vtx.Z());

			//const reco::Vertex *opv = (dynamic_cast <pat::CompositeCandidate *>(chi_cand.daughter("dimuon")))->userData<reco::Vertex>("PVwithmuons");
			//const reco::Vertex *opv = chi_cand.daughter("dimuon")->userData<reco::Vertex>("PVwithmuons");
			//primary_v.SetXYZ(opv->x(), opv->y(), opv->z());

			// 2012 parameterization
			double sigma = Y_sig_par_A + Y_sig_par_B * pow(fabs(dimuon_p4.Rapidity()), 2) +
				Y_sig_par_C * pow(fabs(dimuon_p4.Rapidity()), 3);

			psi1S_nsigma = fabs(dimuon_p4.M() - psi1SMass) / sigma;
			psi2S_nsigma = fabs(dimuon_p4.M() - psi2SMass) / sigma;

			runNumber = iEvent.id().run();
			eventNumber = iEvent.id().event();
			chi_tree->Fill();
		}

	}
	else { std::cout << "no valid chi handle" << std::endl; }

}

// ------------ method called once each job just before starting event loop  ------------
void ChiRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void ChiRootupler::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void ChiRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void ChiRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void ChiRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void ChiRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//void ChiRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
//	edm::ParameterSetDescription desc;
//	desc.setUnknown();
//	descriptions.addDefault(desc);
//}


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

	std::vector< std::pair< double, short> > idx[2];
	short ik = -1;
	BOOST_FOREACH(edm::RefToBase<reco::Track> tk, conv.tracks()) {
		ik++;
		short count = -1;
		BOOST_FOREACH(const reco::Vertex& vtx, priVtxs) {
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

////////////////////////////////////////////////

//define this as a plug-in
DEFINE_FWK_MODULE(ChiRootupler);
