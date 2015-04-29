// -*- C++ -*-
//
// Package:    OniaPhotonRootupler
// Class:      OniaPhotonRootupler
// 
/**
 Description: Create a rootuple of the Onia Photon candidates

 Implementation:  Original work from Stefano Argiro, Alessandro Degano  and the Torino team
                  Adapted by Alberto Sanchez
*/

// system include files
#include <memory>

// user include files
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

class OniaPhotonRootupler:public edm::EDAnalyzer {
      public:
	explicit OniaPhotonRootupler(const edm::ParameterSet &);
	~OniaPhotonRootupler();

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

      private:
	virtual void beginJob();
	virtual void analyze(const edm::Event &, const edm::EventSetup &);
	virtual void endJob();

	virtual void beginRun(edm::Run const &, edm::EventSetup const &);
	virtual void endRun(edm::Run const &, edm::EventSetup const &);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
	virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);

	// ----------member data ---------------------------
	std::string file_name;

        edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_Label;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> ups_Label;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> refit1_Label;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> refit2_Label;
        edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
        edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

	bool isMC_;

	UInt_t run;
	UInt_t event;

	TLorentzVector chi_p4;
	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;
	TLorentzVector photon_p4;

	TLorentzVector rf1S_chi_p4;
	TLorentzVector rf1S_dimuon_p4;
	TLorentzVector rf1S_muonP_p4;
	TLorentzVector rf1S_muonN_p4;
	TLorentzVector rf1S_photon_p4;

	TLorentzVector rf2S_chi_p4;
	TLorentzVector rf2S_dimuon_p4;
	TLorentzVector rf2S_muonP_p4;
	TLorentzVector rf2S_muonN_p4;
	TLorentzVector rf2S_photon_p4;
        TVector3       primary_v;
        TVector3       secondary_v;
        TVector3       dimuon_v;
        TVector3       conversion_v;


	Double_t ele_lowerPt_pt;
	Double_t ele_higherPt_pt;
	Double_t ctpv;
	Double_t ctpv_error;
	Double_t pi0_abs_mass;
	Double_t psi1S_nsigma;
	Double_t psi2S_nsigma;
	Double_t conv_vertex;
	Double_t dz;
	Double_t numPrimaryVertices;
	Int_t trigger;

	Double_t probFit1S;
	Int_t rf1S_rank;
	Double_t probFit2S;

	TTree *chib_tree;

	TLorentzVector gen_chic_p4;
	Int_t chic_pdgId;
	TLorentzVector gen_jpsi_p4;
	//Int_t dimuon_pdgId;
	TLorentzVector gen_photon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;
        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const double pi0_mass = 0.1349766;
static const Double_t psi1SMass = 3.09691;
static const Double_t psi2SMass = 3.68610;

/*
// 2011 par
static const double Y_sig_par_A = 0.058;
static const double Y_sig_par_B = 0.047;
static const double Y_sig_par_C = 0.22;
*/

// 2012 par
static const double Y_sig_par_A = 62.62;
static const double Y_sig_par_B = 56.3;
static const double Y_sig_par_C = -20.77;

//
// constructors and destructor
//

OniaPhotonRootupler::OniaPhotonRootupler(const edm::ParameterSet & iConfig): 
chi_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("chi_cand"))),
ups_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("ups_cand"))),
refit1_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refit1S"))),
refit2_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refit2S"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))), 
isMC_(iConfig.getParameter < bool > ("isMC"))
{

    edm::Service < TFileService > fs;
    chib_tree = fs->make < TTree > ("chicTree", "Tree of chic");

    chib_tree->Branch("run", &run, "run/I");
    chib_tree->Branch("event", &event, "event/I");

    chib_tree->Branch("primary_v",    "TVector3", &primary_v);
    chib_tree->Branch("secondary_v",    "TVector3", &secondary_v);
    chib_tree->Branch("dimuon_v",     "TVector3", &dimuon_v);
    chib_tree->Branch("conversion_v", "TVector3", &conversion_v);

    chib_tree->Branch("chi_p4", "TLorentzVector", &chi_p4);
    chib_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
    chib_tree->Branch("muonP_p4", "TLorentzVector", &muonP_p4);
    chib_tree->Branch("muonN_p4", "TLorentzVector", &muonN_p4);
    chib_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);

    chib_tree->Branch("rf1S_chi_p4", "TLorentzVector", &rf1S_chi_p4);
    chib_tree->Branch("rf1S_dimuon_p4", "TLorentzVector", &rf1S_dimuon_p4);
    chib_tree->Branch("rf1S_muonP_p4", "TLorentzVector", &rf1S_muonP_p4);
    chib_tree->Branch("rf1S_muonN_p4", "TLorentzVector", &rf1S_muonN_p4);
    chib_tree->Branch("rf1S_photon_p4", "TLorentzVector", &rf1S_photon_p4);

    chib_tree->Branch("rf2S_chi_p4", "TLorentzVector", &rf2S_chi_p4);
    chib_tree->Branch("rf2S_dimuon_p4", "TLorentzVector", &rf2S_dimuon_p4);
    chib_tree->Branch("rf2S_muonP_p4", "TLorentzVector", &rf2S_muonP_p4);
    chib_tree->Branch("rf2S_muonN_p4", "TLorentzVector", &rf2S_muonN_p4);
    chib_tree->Branch("rf2S_photon_p4", "TLorentzVector", &rf2S_photon_p4);

    chib_tree->Branch("ele_lowerPt_pt", &ele_lowerPt_pt, "ele_lowerPt_pt/D");
    chib_tree->Branch("ele_higherPt_pt", &ele_higherPt_pt, "ele_higherPt_pt/D");
    chib_tree->Branch("ctpv", &ctpv, "ctpv/D");
    chib_tree->Branch("ctpv_error", &ctpv_error, "ctpv_error/D");
    chib_tree->Branch("pi0_abs_mass", &pi0_abs_mass, "pi0_abs_mass/D");
    chib_tree->Branch("psi1S_nsigma", &psi1S_nsigma, "psi1S_nsigma/D");
    chib_tree->Branch("psi2S_nsigma", &psi2S_nsigma, "psi2S_nsigma/D");

    chib_tree->Branch("conv_vertex", &conv_vertex, "conv_vertex/D");
    chib_tree->Branch("dz", &dz, "dz/D");
    chib_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/D");
    chib_tree->Branch("trigger", &trigger, "trigger/I");
    chib_tree->Branch("probFit1S", &probFit1S, "probFit1S/D");
    chib_tree->Branch("rf1S_rank", &rf1S_rank, "rf1S_rank/I");
    chib_tree->Branch("probFit2S", &probFit2S, "probFit2S/D");

    if (isMC_) {
       chib_tree->Branch("gen_chic_p4",      "TLorentzVector",  &gen_chic_p4);
       chib_tree->Branch("chic_pdgId",   &chic_pdgId,   "chic_pdgId/I");
       chib_tree->Branch("gen_jpsi_p4",    "TLorentzVector",  &gen_jpsi_p4);
       //chib_tree->Branch("dimuon_pdgId", &dimuon_pdgId, "dimuon_pdgId/I");
       chib_tree->Branch("gen_photon_p4",    "TLorentzVector",  &gen_photon_p4);
       chib_tree->Branch("gen_muonP_p4",     "TLorentzVector",  &gen_muonP_p4);
       chib_tree->Branch("gen_muonM_p4",     "TLorentzVector",  &gen_muonM_p4);
    }
    genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
    packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");

}

OniaPhotonRootupler::~OniaPhotonRootupler() {}

//
// member functions
//

//Check recursively if any ancestor of particle is the given one
bool OniaPhotonRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

// ------------ method called for each event  ------------
void OniaPhotonRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;

  Handle < std::vector < pat::CompositeCandidate > >chi_cand_handle;
  iEvent.getByToken(chi_Label, chi_cand_handle);

  Handle < std::vector < pat::CompositeCandidate > >ups_hand;
  iEvent.getByToken(ups_Label, ups_hand);

  Handle < std::vector < pat::CompositeCandidate > >refit1S_handle;
  iEvent.getByToken(refit1_Label, refit1S_handle);

  Handle < std::vector < pat::CompositeCandidate > >refit2S_handle;
  iEvent.getByToken(refit2_Label, refit2S_handle);

  Handle < std::vector < reco::Vertex > >primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle < edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken(triggerResults_Label, triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();

  pat::CompositeCandidate chi_cand;
  pat::CompositeCandidate refit1S;
  pat::CompositeCandidate refit2S;

  if (isMC_) {
   // Pruned particles are the one containing "important" stuff
   edm::Handle<std::vector<reco::GenParticle> > pruned;
   iEvent.getByToken(genCands_,pruned);
   // Packed particles are all the status 1, so usable to remake jets
   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
   edm::Handle<std::vector<pat::PackedGenParticle> > packed;
   iEvent.getByToken(packCands_,packed);

   //let's try to find all status1 originating directly from a Chi_c1 meson decay

   chic_pdgId = 0;
   int foundit = 0;
   for (size_t i=0; i<pruned->size();i++) {
      int p_id = abs((*pruned)[i].pdgId());
      if ( p_id == 20443 || p_id == 445 || p_id == 10443) {
         chic_pdgId = p_id;
         foundit++;
         const reco::Candidate * chic = &(*pruned)[i];
         gen_chic_p4.SetPtEtaPhiM(chic->pt(),chic->eta(),chic->phi(),chic->mass());

         const reco::Candidate * j1 = nullptr;
         const reco::Candidate * j2 = nullptr;

         for (size_t j=0; j<packed->size();j++) {
//get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
            const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
            if (motherInPrunedCollection != nullptr && isAncestor( chic , motherInPrunedCollection)) {
               const reco::Candidate * d = &(*packed)[j];
               int dauId = d->pdgId();

               if (dauId == 13 && motherInPrunedCollection->pdgId() == 443) {
                  gen_muonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
                  j1 = motherInPrunedCollection;
               }
               if (dauId == -13 && motherInPrunedCollection->pdgId() == 443) {
                  gen_muonP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
                  j2 = motherInPrunedCollection;
               }
               if (dauId == 22) {
                  gen_photon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
            }
            if ( foundit == 4 ) break;
         }
         if ( foundit == 4 ) {
            if ( j1 != nullptr && j2 != nullptr && j1 == j2 ) {
               gen_jpsi_p4.SetPtEtaPhiM(j1->pt(),j1->eta(),j1->phi(),j1->mass());
               break;
            }
            else {
              std::cout << "Mother of muons does not match (" << j1->pdgId() << "," << j2->pdgId() << ")" << std::endl;
              foundit = 0;
              chic_pdgId = 0;
            }
         } else {
           std::cout << "Found just " << foundit << " out of 4 particles" << std::endl;
           foundit = 0;
           chic_pdgId = 0;
         }
      }  // if ( p_id
   } // for (size

   if ( ! chic_pdgId ) { // sanity check
      std::cout << "Rootupler does not found the given decay " <<
                iEvent.id().run() << "," << iEvent.id().event() << std::endl;
   } 
}


   //grab Trigger informations
   // save it in variable trigger, trigger is an int between 0 and 15, in binary it is:
   // (pass 11)(pass 8)(pass 7)(pass 5)
   // es. 11 = pass 5, 7 and 11
   // es. 4 = pass only 8

   trigger = 0;
   if (triggerResults_handle.isValid()) {

      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);

      unsigned int NTRIGGERS = 9;
      string TriggersToTest[NTRIGGERS] = {"HLT_Dimuon10_Jpsi_Barrel","HLT_Dimuon16_Jpsi","HLT_Dimuon20_Jpsi","HLT_Dimuon8_Upsilon_Barrel","HLT_Dimuon13_Upsilon","HLT_Dimuon8_PsiPrime_Barrel","HLT_Dimuon13_PsiPrime","HLT_Mu16_TkMu0_dEta18_Onia","HLT_Mu25_TkMu0_dEta18_Onia"};
       
      for (unsigned int i = 0; i < NTRIGGERS; i++) {
         for (int version = 1; version < 5; version++) {
            stringstream ss;
            ss << TriggersToTest[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label().c_str());
            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
               trigger += (1<<i);
               break;
            }
         }
      }

    } else {
      std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
    } // if (trigger...

    bool bestCandidateOnly_ = true;
    // grabbing chi inforamtion
    if (chi_cand_handle.isValid() && chi_cand_handle->size() > 0) {

       unsigned int csize = chi_cand_handle->size();
       if (bestCandidateOnly_) csize = 1;

       for (unsigned int i = 0; i < csize; i++) {
	   chi_cand = chi_cand_handle->at(i);

	   chi_p4.SetPtEtaPhiM(chi_cand.pt(), chi_cand.eta(), chi_cand.phi(), chi_cand.mass());
	   dimuon_p4.SetPtEtaPhiM(chi_cand.daughter("dimuon")->pt(), chi_cand.daughter("dimuon")->eta(), 
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

	   Double_t ele1_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->
                                    userData<reco::Track>("track0"))->pt();
	   Double_t ele2_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->
                                    userData<reco::Track>("track1"))->pt();

	   if (ele1_pt > ele2_pt) {
	      ele_higherPt_pt = ele1_pt;
	      ele_lowerPt_pt = ele2_pt;
	   } else {
	      ele_higherPt_pt = ele2_pt;
	      ele_lowerPt_pt = ele1_pt;
	   }

	   ctpv = (dynamic_cast < pat::CompositeCandidate * >(chi_cand.daughter("dimuon")))->userFloat("ppdlPV");
	   ctpv_error = (dynamic_cast < pat::CompositeCandidate * >(chi_cand.daughter("dimuon")))->userFloat("ppdlErrPV");
	   pi0_abs_mass = 1.;//pi0_abs_values[0];

           reco::Candidate::Point vtx = chi_cand.daughter("photon")->vertex();
           conversion_v.SetXYZ(vtx.X(),vtx.Y(),vtx.Z());
	   conv_vertex = vtx.rho(); //chi_cand.daughter("photon")->vertex().rho();

           reco::Candidate::Point dimuon_vtx = chi_cand.daughter("dimuon")->vertex();
           dimuon_v.SetXYZ(dimuon_vtx.X(),dimuon_vtx.Y(),dimuon_vtx.Z());

           const reco::Vertex *opv = (dynamic_cast < pat::CompositeCandidate * >(chi_cand.daughter("dimuon")))->userData<reco::Vertex>("PVwithmuons");
           primary_v.SetXYZ(opv->x(),opv->y(),opv->z());

	   dz = chi_cand.userFloat("dz");

	   // 2012 parameterization
	   double sigma = Y_sig_par_A + Y_sig_par_B * pow(fabs(dimuon_p4.Rapidity()), 2) + 
                          Y_sig_par_C * pow(fabs(dimuon_p4.Rapidity()), 3);

	   psi1S_nsigma = fabs(dimuon_p4.M() - psi1SMass) / sigma;
	   psi2S_nsigma = fabs(dimuon_p4.M() - psi2SMass) / sigma;

           int j = -1;
           if (refit1S_handle.isValid() && i < refit1S_handle->size()) { j = (refit1S_handle->at(i)).userInt("Index"); }

	   if ( j >= 0  && (unsigned int) j == i ) {

	      refit1S = refit1S_handle->at(i);
              reco::Candidate::Point s_vtx = refit1S.vertex();
              secondary_v.SetXYZ(s_vtx.X(),s_vtx.Y(),s_vtx.Z());

	      rf1S_chi_p4.SetPtEtaPhiM(refit1S.pt(), refit1S.eta(), refit1S.phi(), refit1S.mass());
	      rf1S_dimuon_p4.SetPtEtaPhiM(refit1S.daughter("dimuon")->pt(), refit1S.daughter("dimuon")->eta(), 
                                          refit1S.daughter("dimuon")->phi(), refit1S.daughter("dimuon")->mass());
	      rf1S_photon_p4.SetPtEtaPhiM(refit1S.daughter("photon")->pt(), refit1S.daughter("photon")->eta(), 
                                          refit1S.daughter("photon")->phi(), refit1S.daughter("photon")->mass());

	      reco::Candidate::LorentzVector vP = refit1S.daughter("dimuon")->daughter("muon1")->p4();
	      reco::Candidate::LorentzVector vM = refit1S.daughter("dimuon")->daughter("muon2")->p4();

	      if (refit1S.daughter("dimuon")->daughter("muon1")->charge() < 0) {
		 vP = refit1S.daughter("dimuon")->daughter("muon2")->p4();
		 vM = refit1S.daughter("dimuon")->daughter("muon1")->p4();
	      }

	      rf1S_muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	      rf1S_muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	      probFit1S = refit1S.userFloat("vProb");

	      rf1S_rank = i;
	   } else {
              //if (j>=0) std::cout << "Something wrong in refit1S_handle " << i << " " << j << std::endl;
	      rf1S_chi_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf1S_dimuon_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf1S_photon_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf1S_muonP_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf1S_muonN_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      probFit1S = 0;
	      rf1S_rank = -1;
	   }	// if rf1S is valid

           j = -1;
           if (refit2S_handle.isValid() && i < refit2S_handle->size()) { j = (refit2S_handle->at(i)).userInt("Index"); }
           if ( j >= 0  && (unsigned int) j == i ) {
	      refit2S = refit2S_handle->at(i);

	      rf2S_chi_p4.SetPtEtaPhiM(refit2S.pt(), refit2S.eta(), refit2S.phi(), refit2S.mass());
	      rf2S_dimuon_p4.SetPtEtaPhiM(refit2S.daughter("dimuon")->pt(), refit2S.daughter("dimuon")->eta(), 
                                          refit2S.daughter("dimuon")->phi(), refit2S.daughter("dimuon")->mass());
	      rf2S_photon_p4.SetPtEtaPhiM(refit2S.daughter("photon")->pt(), refit2S.daughter("photon")->eta(), 
                                          refit2S.daughter("photon")->phi(), refit2S.daughter("photon")->mass());

	      reco::Candidate::LorentzVector vP = refit2S.daughter("dimuon")->daughter("muon1")->p4();
	      reco::Candidate::LorentzVector vM = refit2S.daughter("dimuon")->daughter("muon2")->p4();

	      if (refit2S.daughter("dimuon")->daughter("muon1")->charge() < 0) {
		 vP = refit2S.daughter("dimuon")->daughter("muon2")->p4();
		 vM = refit2S.daughter("dimuon")->daughter("muon1")->p4();
	      }

	      rf2S_muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	      rf2S_muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	      probFit2S = refit2S.userFloat("vProb");

           } else {
	      rf2S_chi_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf2S_dimuon_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf2S_photon_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf2S_muonP_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      rf2S_muonN_p4.SetPtEtaPhiM(0, 0, 0, 0);
	      probFit2S = 0;
	   }	// if rf2S valid

	   run = iEvent.id().run();
	   event = iEvent.id().event();
	   chib_tree->Fill();
	}		// for i on chi_cand_handle

	} else { std::cout << "no valid chi handle" << std::endl; }			// if chi handle valid
}

// ------------ method called once each job just before starting event loop  ------------
void OniaPhotonRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void OniaPhotonRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void OniaPhotonRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void OniaPhotonRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void OniaPhotonRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void OniaPhotonRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OniaPhotonRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaPhotonRootupler);
