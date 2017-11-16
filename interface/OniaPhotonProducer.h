#ifndef __OniaPhotonProducer_h_
#define __OniaPhotonProducer_h_

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


class OniaPhotonProducer : public edm::EDProducer {

 public:
  explicit OniaPhotonProducer(const edm::ParameterSet& ps);
 
 private:

  virtual void produce(edm::Event& event, const edm::EventSetup& esetup);  
  virtual void endJob();

  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> photon_Label;

  const pat::CompositeCandidate makeChiCandidate(const pat::CompositeCandidate&, 
						 const pat::CompositeCandidate&);
 
  float Getdz(const pat::CompositeCandidate&, const reco::Candidate::Point &);
  // check if the mass difference is in desired range
  bool cutDeltaMass(const pat::CompositeCandidate&,const pat::CompositeCandidate&); 

  bool cutdz(float dz){return dz<dzMax_; }

  bool pi0OnlineSwitch_;
  
  // delta mass range
  std::vector<double> deltaMass_;
  double dzMax_;

  // use only trigger-matched J/Psi or Upsilon   
  bool triggerMatch_;  
    
  int candidates;
  int delta_mass_fail;
  int dz_cut_fail;
  int pizero_fail;
};

#endif // __OniaPhotonProducer_h_
