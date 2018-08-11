#include <HeavyIonsAnalysis/ChiAnalysis/interface/ChiProducer.h>



ChiProducer::ChiProducer(const edm::ParameterSet& ps): //initialization list follows
  dimuon_label(consumes<pat::CompositeCandidateCollection>(ps.getParameter< edm::InputTag>("dimuon_cand"))),
  photon_label(consumes<pat::CompositeCandidateCollection>(ps.getParameter< edm::InputTag>("photon_cand")))
 // pi0OnlineSwitch_(ps.getParameter<bool>("pi0OnlineSwitch")),
 // deltaMass_(ps.getParameter<std::vector<double> >("deltaMass")),
//  dzMax_(ps.getParameter<double>("dzmax")),
//  triggerMatch_(ps.getParameter<bool>("triggerMatch"))*/
{
  produces<pat::CompositeCandidateCollection>("ChiCandidates");
  //candidates = 0;
  //delta_mass_fail = 0;
  //dz_cut_fail = 0;
  //pizero_fail = 0;
}
 
ChiProducer::~ChiProducer() {}
void ChiProducer::beginJob(const edm::EventSetup& iSetup) {}

void ChiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
	
	std::auto_ptr<pat::CompositeCandidateCollection> chiCandColl(new pat::CompositeCandidateCollection);

	edm::Handle<pat::CompositeCandidateCollection> dimuon_handle;
	iEvent.getByToken(dimuon_label, dimuon_handle);

	edm::Handle<pat::CompositeCandidateCollection> photon_handle;
	iEvent.getByToken(photon_label, photon_handle);

    // Note: since Dimuon cand are sorted by decreasing vertex probability then the first chi cand is the one associated with the "best" dimuon 
    for (pat::CompositeCandidateCollection::const_iterator  dimuonCand = dimuon_handle->begin(); dimuonCand!= dimuon_handle->end(); ++dimuonCand){

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
			double dz = fabs(Getdz(*photCand, dimuonCand->vertex())); // 
			//double dz = fabs( photCand->dz(dimuonCand->vertex()) );
			chiCand.addUserFloat("dz", (float)dz);

			//if (!cutdz(dz)) {
			//	dz_cut_fail++;
			//	continue;
			//}

			//int flags = (photCand->userInt("flags") % 32);
			//bool pi0_fail = flags & 8;
			//if (pi0OnlineSwitch_ && pi0_fail) {
			//	pizero_fail++;
			//	continue;
			//}

			chiCandColl->push_back(chiCand);
			candidates++;
		}
	}
	iEvent.put(chiCandColl, "ChiCandidates");
}

double ChiProducer::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {

  reco::Candidate::LorentzVector mom = c.p4();
  reco::Candidate::Point vtx = c.vertex();
  
  double dz = (vtx.Z()-p.Z()) - ((vtx.X()-p.X())*mom.X()+(vtx.Y()-p.Y())*mom.Y())/mom.Rho() * mom.Z()/mom.Rho();
  return dz;  
  
}

void ChiProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "Chi Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  //std::cout << "Delta mass fail: " << delta_mass_fail << std::endl;
  //std::cout << "Dz fail:         " << dz_cut_fail << std::endl;
  //std::cout << "Pi0 fail:        " << pizero_fail << std::endl;
  //std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " Chi candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}
  
const pat::CompositeCandidate ChiProducer::makeChiCandidate(const pat::CompositeCandidate& dimuon,
				  const pat::CompositeCandidate& photon){
  pat::CompositeCandidate chiCand;
  chiCand.addDaughter(dimuon,"dimuon");
  chiCand.addDaughter(photon,"photon");
  chiCand.setVertex(dimuon.vertex());
  reco::Candidate::LorentzVector vChic = dimuon.p4() + photon.p4();
  chiCand.setP4(vChic);
  return chiCand;
}
/*
// check if the mass difference is in desired range
bool ChiProducer::cutDeltaMass(const pat::CompositeCandidate& chiCand,
				   const pat::CompositeCandidate& dimuonCand){
  float deltam = chiCand.p4().M() - dimuonCand.p4().M();
  float m1     = deltaMass_[0];
  float m2     = deltaMass_[1];
  return (deltam > m1 && deltam < m2);
}
//*/
//define this as a plug-in
DEFINE_FWK_MODULE(ChiProducer);
