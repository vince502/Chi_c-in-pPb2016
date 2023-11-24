#ifndef ChiRootupler_h
#define ChiRootupler_h

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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <TClonesArray.h>
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

//
// class declaration
//

class ChiRootupler : public edm::one::EDAnalyzer<> {
public:
	explicit ChiRootupler(const edm::ParameterSet &);
	~ChiRootupler();
	//static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
	//bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

private:
	
	virtual void analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup);						// ------------ method called for each event  ------------
	virtual void beginJob();																				// ------------ method called once each job just before starting event loop  ------------
	virtual void endJob();																					// ------------ method called once each job just after ending the event loop  ------------
	virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);								// ------------ method called when starting to processes a run  ------------
	virtual void endRun(edm::Run const &, edm::EventSetup const &);											// ------------ method called when ending the processing of a run  ------------
	virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);				// ------------ method called when starting to processes a luminosity block  ------------
	virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);					// ------------ method called when ending the processing of a luminosity block  ------------

	// general functions
	void Clear();
	int SelectVertex(edm::Handle<std::vector<reco::Vertex> >& priVtxs, double zPos); //returns the index (in the collection) of vertex matched to the z position

	// photon functions

	void Conv_removeDuplicates(edm::Handle <std::vector <reco::Conversion>> convCollIn, reco::ConversionCollection* convCollOut); // Removes the duplicate conversions, since they mess up matching. We remove them in an early pass, so all the other steps are done with a cleaned-up collection

	static bool lt_comparator(std::pair<double, short> a, std::pair<double, short> b); // comparator for checking conversions
	bool Conv_checkTkVtxCompatibility(const reco::Conversion& conv, const reco::VertexCollection&  priVtxs, double sigmaTkVtxComp_, bool& Flag_Best_Out, bool& Flag_SecondBestA_Out, bool& Flag_SecondBestB_Out, double& sigmaMinValue1Out, double& sigmaMinValue2Out);
	bool Conv_foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB);
	int Conv_checkDuplicity(const reco::Conversion& conv, int convPos, edm::Handle < std::vector <reco::Conversion>>& conversion_handle, double& dROut, double& dpTOut);
	int Conv_checkDuplicityAlbertoVersion(const reco::Conversion& conv, int convPos, edm::Handle < std::vector <reco::Conversion>>& conversion_handle, double& dROut, double& dpTOut);
	//photon MC
	bool Conv_isMatched(const math::XYZTLorentzVectorF& reco_conv, const reco::GenParticle& gen_phot, double maxDeltaR, double maxDPtRel);
	bool ConvSelection(const reco::Conversion& conv);


	// chi functions
	const pat::CompositeCandidate makeChiCandidate(const pat::CompositeCandidate& dimuon, const pat::CompositeCandidate& photon);
	double Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point& p); //gets dz when dxy to the point is the smallest
	double Getdxy(const pat::CompositeCandidate& c, const reco::Candidate::Point& p); //gets dxy when dz to the point =0
	const pat::CompositeCandidate makeChiCandidate(const pat::CompositeCandidate& dimuon, const reco::Conversion& photon);
	double Getdz(const reco::Conversion& c, const reco::Candidate::Point& p); //gets dz when dxy to the point is the smallest
	double Getdxy(const reco::Conversion& c, const reco::Candidate::Point& p); //gets dxy when dz to the point =0

	// gen functions
	template <typename particle_in, typename particle_type>
	int MatchGen(particle_in& genParticle, edm::Handle< std::vector <particle_type>>& collToBeMatched, double maxDeltaR, double maxDPtRel, int& nMatches_Out, double& rDelta_Out, double& ptDeltaRel_Out);
	template <typename particle_in, typename particle_type>
	int MatchGen(particle_in& genParticle, std::vector <particle_type>* collToBeMatched, double maxDeltaR, double maxDPtRel, int& nMatches_Out, double& rDelta_Out, double& ptDeltaRel_Out);
	template <typename T>  reco::LeafCandidate GetRecoCandidate(const T& a); //default, however momentum is stored differently for conversion and pats
	reco::LeafCandidate GetRecoCandidate(const reco::Conversion& conv); //momentum is stored differently for conversion
	reco::LeafCandidate GetRecoCandidate(const pat::CompositeCandidate& comp);
	reco::LeafCandidate GetRecoCandidate(const pat::Muon& comp);
	reco::LeafCandidate GetRecoCandidate(const reco::GenParticle& comp);

	//template <typename particle_type>
	//int MatchGen(particle_type& genParticle, double maxDeltaR, double maxDPtRel) { return 0; }

	//int MatchGen(const reco::Candidate& genParticle, double maxDeltaR, double maxDPtRel) { return 0; }

	// ----------member data ---------------------------
	std::string file_name;

	//edm::EDGetTokenT< edm::View <pat::Muon> > muon_label; //is a muon collection
	edm::EDGetTokenT<pat::MuonCollection>  muon_label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_label;
	//edm::EDGetTokenT<pat::CompositeCandidateCollection> photon_label;
	edm::EDGetTokenT<reco::ConversionCollection> conversion_label;
	//edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_label;
	edm::EDGetTokenT<reco::VertexCollection> primaryVertices_label;
	edm::EDGetTokenT<edm::TriggerResults> triggerResults_label;
	edm::EDGetTokenT<reco::Centrality> centrality_label;
	edm::EDGetTokenT<reco::TrackCollection> srcTracks_label;
	edm::EDGetTokenT<reco::GenParticleCollection> genParticles_label;

	bool flag_doMC;
	bool flag_saveExtraThings = false;  //saves all muons, and other not so necessary data. Should be false unless some debugging needs it
	bool flag_doKinematicRefit = true; // do kinematic refit based on github.com/alberto-sanchez/chi-analysis-miniaod/blob/master/src/OniaPhotonKinematicFit.cc from Alberto 
	
	const int PythCode_chic0 = 10441; //Pythia codes
	const int PythCode_chic1 = 20443;
	const int PythCode_chic2 = 445;
	double jpsi_mass = 3.097;

	// constants for events and PV
	const double cPVMatching_cutoff = 0.0; //If dimuon is within 2mm of PV[0], leave it. Based on dimuon.z - PV.z resolution. Has effect for peripheral events //Changed to 0 (no cut-off)


	// constants for muon cuts
	const double muon_maxDeltaR = 0.5;
	const double muon_maxDPtRel = 0.5;

	// constants for dimuon cuts
	const double jpsi_maxDeltaR = 0.5;
	const double jpsi_maxDPtRel = 0.5;

	// constants for conversion cuts
	const double conv_TkVtxCompSigmaCut = 5.0;
	const double conv_maxDeltaR = 0.2;
	const double conv_maxDPtRel = 1;

	const double conv_duplicateMaxDeltapT = 0.1; //what to consider a split track for purposes of removing those conversions (otherwise there are almost identical conversions in the selection)
	const double conv_duplicateMaxDeltaR = 0.05;

	std::string triggerName = "HLT_PAL1DoubleMuOpen_v1";
	std::string triggerFilter = "hltL1fL1sDoubleMuOpenBptxANDL1Filtered0";

	TTree* event_tree;

	//general
	long runNumber;
	long eventNumber;
	int nPrimVertices; 
	int muonPerEvent;
	int convPerTriggeredEvent;
	int dimuonPerEvent;
	int chiCandPerEvent;

	int ntracks_inEvent;
	double hfTowerSum_inEvent;
	HLTConfigProvider hltConfig;
	int Trig_Event_HLTDoubleMuOpen;


	// vertex
	std::vector <double> pvtx_z;
	std::vector <double> pvtx_zError;
	std::vector <double> pvtx_x;
	std::vector <double> pvtx_y;
	std::vector <double> pvtx_nTracksUncut;
	std::vector <double> pvtx_nTracks;
	std::vector <bool> pvtx_isFake;

	//muon info
	std::vector <bool> muonIsHLTDoubleMuOpen;
	std::vector <bool> muonIsHLTDoubleMuOpenFilter;
	std::vector <bool> muonIsGlobal;
	std::vector <bool> muonIsTracker;
	std::vector <bool> muonIsPF;
	std::vector <bool> muonIsSoft;
	std::vector <bool> muonIsTight;
	std::vector <bool> muonIsNotGlobalNorTracker;
	std::vector <bool> muonIDHas_TMOneStationTight;
	std::vector <int>  muon_pvtx_index;
	std::vector <double> muonInnerTrack_dxy;
	std::vector <double> muonInnerTrack_dz;
	std::vector <int> muonTrackerLayersWithMeasurement;
	std::vector <int> muonPixelLayersWithMeasurement;
	std::vector <bool> muonQuality_isHighPurity;
	std::vector <int> muon_charge;
	std::vector <double> muon_eta;
	std::vector <double> muon_pt;
	TClonesArray* muon_p4; //TLorentzVector
	std::vector <pat::Muon> patMuonStored;

	//muon MC  - mostly testing, not that useful
	std::vector <bool> muon_isMatchedMC;
	std::vector <double> muonGen_eta;
	std::vector <double> muonGen_pt;
	TClonesArray* muonGen_p4; //TLorentzVector
	std::vector <double> muonGen_rDelta;
	std::vector <double> muonGen_ptDelta;
	std::vector <double> muonGen_ptDeltaRel;


	//dimuon info

	TClonesArray*  dimuon_p4; //TLorentzVector
	std::vector <double> dimuon_eta;
	std::vector <double> dimuon_pt;
	std::vector <double> dimuon_charge; 
	TClonesArray* dimuon_vtx; //TVector3
	std::vector <int>  dimuon_pvtx_indexFromOniaMuMu;
    std::vector <int>  dimuon_pvtx_index;
	std::vector <double> dimuon_dz_dimuonvtx_pvtx;
	std::vector <double> dimuon_vtxProb;
	std::vector <pat::CompositeCandidate> dimuonStored;
	std::vector <int> dimuon_muon1_position; //stores position of first muon in muon collection (probably one with higher pT due to ordering in the collection)
	std::vector <int> dimuon_muon2_position; //stores position of second muon in muon collection 
	std::vector <double> dimuon_ctpv; 
	std::vector <double> dimuon_ctpvError; 
	


	//conversion info

	std::vector <int> convRaw_duplicityStatus; // 0: is not duplicate to any, 1: shares a track, but isn't split 2: shares a track, and is split, but is kept, 3: doesn't have 2 tracks, 4: shares the track, is split, and is removed
	std::vector <double> convRaw_splitDR;
	std::vector <double> convRaw_splitDpT;
	std::vector <int> convRaw_duplicityStatus_AV; // 0: is not duplicate to any, 1: shares a track, but is largest in prob 2: shares a track, and is not largest in prob, 3: doesn't have 2 tracks
	std::vector <double> convRaw_splitDR_AV;
	std::vector <double> convRaw_splitDpT_AV;
	std::vector <int> conv_positionRaw;
	std::vector <int> conv_tk1ValidHits;
	std::vector <int> conv_tk2ValidHits;
	std::vector <bool> convQuality_isHighPurity;
	std::vector <bool> convQuality_isGeneralTracksOnly;
	TClonesArray* conv_vtx; //TVector3
	std::vector <double> conv_vertexPositionRho;
	std::vector <double> conv_sigmaTkVtx1;
	std::vector <double> conv_sigmaTkVtx2;
	std::vector <bool> conv_tkVtxCompatibilityOK;
	std::vector <bool> conv_tkVtxCompatible_bestVertex;
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexA;
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexB;
	std::vector <bool> conv_tkVtxCompatibilityOK_test; //test ones
	std::vector <bool> conv_tkVtxCompatible_bestVertex_test; //test ones
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexA_test; //test
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexB_test; //test

	std::vector <int> conv_compatibleInnerHitsOK; //-1: less than 2 tracks, 0: not compatible, 1: yes
	std::vector <reco::HitPattern> conv_hitPat1;
	std::vector <reco::HitPattern> conv_hitPat2;
	std::vector <double> conv_vertexChi2Prob;
	std::vector <int>  conv_pvtx_index;
	std::vector <double> conv_zOfPriVtx; // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
	std::vector <double> conv_zOfPriVtxFromTracks;
	std::vector <double> conv_dzToClosestPriVtx;
	std::vector <double> conv_dxyPriVtx_Tr1;
	std::vector <double> conv_dxyPriVtx_Tr2;
	std::vector <double> conv_dxyPriVtxTimesCharge_Tr1;
	std::vector <double> conv_dxyPriVtxTimesCharge_Tr2;
	std::vector <double> conv_dxyError_Tr1;
	std::vector <double> conv_dxyError_Tr2;

	std::vector <int> conv_tk1NumOfDOF;
	std::vector <int> conv_tk2NumOfDOF;
	std::vector <double> conv_track1Chi2;
	std::vector <double> conv_track2Chi2;
	std::vector <double> conv_Tr1_pt;
	std::vector <double> conv_Tr2_pt;

	std::vector <double> conv_minDistanceOfApproach;
	TClonesArray*  conv_p4; //TLorentzVector
	std::vector <double> conv_eta;
	std::vector <double> conv_pt;


	//conv MC  - mostly testing, not that useful
	std::vector <bool> conv_isMatchedMC;
	std::vector <double> convGen_eta;
	std::vector <double> convGen_pt;
	TClonesArray* convGen_p4; //TLorentzVector
	std::vector <double> convGen_rDelta;
	std::vector <double> convGen_ptDelta;
	std::vector <double> convGen_ptDeltaRel;
	std::vector <int> convGen_motherCode;
	

	// MC general - main MC tree

	std::vector <bool> gen_isGoodChicDecay; // saves true if chic decay was -> J/psi gamma -> mumu(+gammas) gamma; false otherwise
	std::vector <int> gen_pdgId;
	std::vector <double> gen_chic_pt;
	std::vector <double> gen_chic_eta;
	TClonesArray* gen_chic_p4; //TLorentzVector
	std::vector <int> gen_chic_matchPosition;
	std::vector <int> gen_chic_nMatches;
	std::vector <double> gen_chic_rDelta; //in principle duplicates information
	std::vector <double> gen_chic_ptDeltaRel;//in principle duplicates information

	std::vector <double> gen_Jpsi_pt;
	std::vector <double> gen_Jpsi_eta;
	std::vector <int> gen_Jpsi_matchPosition;
	std::vector <int> gen_Jpsi_nMatches;
	std::vector <double> gen_Jpsi_rDelta; //in principle duplicates information
	std::vector <double> gen_Jpsi_ptDeltaRel;//in principle duplicates information
	TClonesArray* gen_Jpsi_p4; //TLorentzVector

	std::vector <int> gen_Jpsi_photon_n;
	std::vector <double> gen_Jpsi_photon_pt;
	TClonesArray* gen_Jpsi_photon_p4; //TLorentzVector

	std::vector <int> gen_muon_charge;
	std::vector <double> gen_muon_pt;
	std::vector <double> gen_muon_eta;
	std::vector <int> gen_muon_matchPosition;
	std::vector <int> gen_muon_nMatches;
	std::vector <double> gen_muon_rDelta; //in principle duplicates information
	std::vector <double> gen_muon_ptDeltaRel;//in principle duplicates information
	TClonesArray* gen_muon_p4; //TLorentzVector

	std::vector <double> gen_phot_pt;
	std::vector <double> gen_phot_eta;
	TClonesArray* gen_phot_p4; //TLorentzVector
	std::vector <int> gen_conv_matchPosition;
	std::vector <int> gen_conv_nMatches;
	std::vector <double> gen_conv_rDelta; //in principle duplicates information
	std::vector <double> gen_conv_ptDeltaRel;//in principle duplicates information


	// chi

	TClonesArray*  chi_p4; //TLorentzVector
	std::vector <double> chi_eta;
	std::vector <double> chi_pt;
	pat::CompositeCandidate chi_cand;
	std::vector <int> chi_daughterJpsi_position; //stores position of daughter Jpsi in dimuon collection
	std::vector <int> chi_daughterConv_position; //stores position of daughter photon (conversion)
	std::vector <double> chi_dzPhotToDimuonVtx; //z distance of photon to dimuon vertex when dxy is minimal
	std::vector <double> chi_dxyPhotToDimuonVtx; //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions
	std::vector <pat::CompositeCandidate> chiStored;
	std::vector <int> chi_kinematicRefitFlag; // -1 kinematic refit not done, 1 done: +2 needed extra refit for photon +4 something wrong with photon at the end +8 something wrong with the final fit 
	std::vector <int> chi_refit_origChicPosition; //stores position of the chic candidate for the refit (there will be gaps if refit fails) 
	std::vector <pat::CompositeCandidate> chi_refitStored;
	std::vector <double> chi_refit_vprob;
	std::vector <double> chi_refit_ctauPV;
	std::vector <double> chi_refit_ctauErrPV;
	std::vector <double> chi_refit_ctauPV3D;
	std::vector <double> chi_refit_pvtxFromPVwithMuons_x;
	std::vector <double> chi_refit_pvtxFromPVwithMuons_y;
	std::vector <double> chi_refit_pvtxFromPVwithMuons_z;


};

#endif 

