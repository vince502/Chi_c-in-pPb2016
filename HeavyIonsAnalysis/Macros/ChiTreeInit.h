#ifndef ChiTreeInit_h
#define ChiTreeInit_h

  // Declaration of branches that are saved for chi


#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <TClonesArray.h>
#include <vector>
#include <sstream>

//#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "FWCore/Common/interface/TriggerNames.h"
//#include "DataFormats/HeavyIonEvent/interface/Centrality.h"


// not really used:
//pat::Muon
//pat::CompositeCandidate
//reco::HitPattern

//general
long long int runNumber;
long long int eventNumber;
int nPrimVertices;
int muonPerEvent;
int convPerTriggeredEvent;
int dimuonPerEvent;
int chiCandPerEvent;

int ntracks_inEvent;
double hfTowerSum_inEvent;
int Trig_Event_HLTDoubleMuOpen;


// vertex
std::vector <double>* pvtx_z = 0;
std::vector <double>* pvtx_zError = 0;
std::vector <double>* pvtx_x = 0;
std::vector <double>* pvtx_y = 0;
std::vector <double>* pvtx_nTracks = 0;
std::vector <bool>* pvtx_isFake = 0;

//muon info
std::vector <bool>* muonIsHLTDoubleMuOpen = 0;
std::vector <bool>* muonIsHLTDoubleMuOpenFilter = 0;
std::vector <bool>* muonIsGlobal = 0;
std::vector <bool>* muonIsTracker = 0;
std::vector <bool>* muonIsPF = 0;
std::vector <bool>* muonIsSoft = 0;
std::vector <bool>* muonIsTight = 0;
std::vector <bool>* muonIsNotGlobalNorTracker = 0;
std::vector <bool>* muonIDHas_TMOneStationTight = 0;
std::vector <int>*  muon_pvtx_index = 0;
std::vector <double>* muonInnerTrack_dxy = 0;
std::vector <double>* muonInnerTrack_dz = 0;
std::vector <int>* muonTrackerLayersWithMeasurement = 0;
std::vector <int>* muonPixelLayersWithMeasurement = 0;
std::vector <bool>* muonQuality_isHighPurity = 0;
std::vector <int>* muon_charge = 0;
std::vector <double>* muon_eta = 0;
std::vector <double>* muon_pt = 0;
TClonesArray* muon_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
//std::vector <pat::Muon>* patMuonStored = 0;

//muon MC
std::vector <bool>* muon_isMatchedMC = 0;
std::vector <double>* muonGen_eta = 0;
std::vector <double>* muonGen_pt = 0;
TClonesArray* muonGen_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
std::vector <double>* muonGen_rDelta = 0;
std::vector <double>* muonGen_ptDelta = 0;
std::vector <double>* muonGen_ptDeltaRel = 0;


//dimuon info

TClonesArray*  dimuon_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
std::vector <double>* dimuon_eta = 0;
std::vector <double>* dimuon_pt = 0;
std::vector <double>* dimuon_charge = 0;
TClonesArray* dimuon_vtx = new TClonesArray("TVector3", 100); //TVector3
std::vector <int>*  dimuon_pvtx_indexFromOniaMuMu = 0;
std::vector <int>*  dimuon_pvtx_index = 0;
std::vector <double>* dimuon_dz_dimuonvtx_pvtx = 0;
std::vector <double>* dimuon_vtxProb = 0;
//std::vector <pat::CompositeCandidate>* dimuonStored = 0;
std::vector <int>* dimuon_muon1_position = 0; //stores position of first muon in muon collection (probably one with higher pT due to ordering in the collection)
std::vector <int>* dimuon_muon2_position = 0; //stores position of second muon in muon collection 
std::vector <double>* dimuon_ctpv = 0;
std::vector <double>* dimuon_ctpvError = 0;



//conversion info

std::vector <int>* conv_duplicityStatus = 0; // 0: is not duplicate to any, 1: shares a track, but is largest in prob 2: shares a track, and is not largest in prob, 3: doesn't have 2 tracks
std::vector <double>* conv_splitDR = 0;
std::vector <double>* conv_splitDpT = 0;
std::vector <int>* conv_tk1ValidHits = 0;
std::vector <int>* conv_tk2ValidHits = 0;
std::vector <bool>* convQuality_isHighPurity = 0;
std::vector <bool>* convQuality_isGeneralTracksOnly = 0;
TClonesArray* conv_vtx = new TClonesArray("TVector3", 100); //TVector3
std::vector <double>* conv_vertexPositionRho = 0;
std::vector <double>* conv_sigmaTkVtx1 = 0;
std::vector <double>* conv_sigmaTkVtx2 = 0;
std::vector <bool>* conv_tkVtxCompatibilityOK = 0;
std::vector <bool>* conv_tkVtxCompatible_bestVertex = 0;
std::vector <bool>* conv_tkVtxCompatible_secondBestVertexA = 0;
std::vector <bool>* conv_tkVtxCompatible_secondBestVertexB = 0;
//std::vector <bool>* conv_tkVtxCompatibilityOK_test = 0; //test ones
//std::vector <bool>* conv_tkVtxCompatible_bestVertex_test = 0; //test ones
//std::vector <bool>* conv_tkVtxCompatible_secondBestVertexA_test = 0; //test
//std::vector <bool>* conv_tkVtxCompatible_secondBestVertexB_test = 0; //test

std::vector <int>* conv_compatibleInnerHitsOK = 0; //-1: less than 2 tracks, 0: not compatible, 1: yes
//std::vector <reco::HitPattern>* conv_hitPat1 = 0;
//std::vector <reco::HitPattern>* conv_hitPat2 = 0;
//std::vector <bool>* conv_isCustomHighPurity = 0;//tbd - is just a sum of some other cuts, not creating at the time
std::vector <double>* conv_vertexChi2Prob = 0;
std::vector <int>*  conv_pvtx_index = 0;
std::vector <double>* conv_zOfPriVtx = 0; // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
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
std::vector <double>* conv_Tr1_pt = 0;
std::vector <double>* conv_Tr2_pt = 0;

std::vector <double>* conv_minDistanceOfApproach = 0;
TClonesArray*  conv_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
std::vector <double>* conv_eta = 0;
std::vector <double>* conv_pt = 0;


//conv MC
std::vector <bool>* conv_isMatchedMC = 0;
std::vector <double>* convGen_eta = 0;
std::vector <double>* convGen_pt = 0;
TClonesArray* convGen_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
std::vector <double>* convGen_rDelta = 0;
std::vector <double>* convGen_ptDelta = 0;
std::vector <double>* convGen_ptDeltaRel = 0;
std::vector <int>* convGen_motherCode = 0;


// MC general

std::vector <bool>* gen_isGoodChicDecay = 0; // saves true if chic decay was ->* J/psi gamma ->* mumu(+gammas) gamma = 0; false otherwise
std::vector <int>* gen_pdgId = 0;
std::vector <double>* gen_chic_pt = 0;
std::vector <double>* gen_chic_eta = 0;
TClonesArray* gen_chic_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
std::vector <int>* gen_chic_matchPosition = 0;
std::vector <int>* gen_chic_nMatches = 0;
std::vector <double>* gen_chic_rDelta = 0; //in principle duplicates information
std::vector <double>* gen_chic_ptDeltaRel = 0;//in principle duplicates information

std::vector <double>* gen_Jpsi_pt = 0;
std::vector <double>* gen_Jpsi_eta = 0;
std::vector <int>* gen_Jpsi_matchPosition = 0;
std::vector <int>* gen_Jpsi_nMatches = 0;
std::vector <double>* gen_Jpsi_rDelta = 0; //in principle duplicates information
std::vector <double>* gen_Jpsi_ptDeltaRel = 0;//in principle duplicates information
TClonesArray* gen_Jpsi_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector

std::vector <int>* gen_Jpsi_photon_n = 0;
std::vector <double>* gen_Jpsi_photon_pt = 0;
TClonesArray* gen_Jpsi_photon_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector

std::vector <int>* gen_muon_charge = 0;
std::vector <double>* gen_muon_pt = 0;
std::vector <double>* gen_muon_eta = 0;
std::vector <int>* gen_muon_matchPosition = 0;
std::vector <int>* gen_muon_nMatches = 0;
std::vector <double>* gen_muon_rDelta = 0; //in principle duplicates information
std::vector <double>* gen_muon_ptDeltaRel = 0;//in principle duplicates information
TClonesArray* gen_muon_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector

std::vector <double>* gen_phot_pt = 0;
std::vector <double>* gen_phot_eta = 0;
TClonesArray* gen_phot_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
std::vector <int>* gen_conv_matchPosition = 0;
std::vector <int>* gen_conv_nMatches = 0;
std::vector <double>* gen_conv_rDelta = 0; //in principle duplicates information
std::vector <double>* gen_conv_ptDeltaRel = 0;//in principle duplicates information


// chi

TClonesArray*  chi_p4 = new TClonesArray("TLorentzVector", 100); //TLorentzVector
std::vector <double>* chi_eta = 0;
std::vector <double>* chi_pt = 0;
std::vector <int>* chi_daughterJpsi_position = 0; //stores position of daughter Jpsi in dimuon collection
std::vector <int>* chi_daughterConv_position = 0; //stores position of daughter photon (conversion)
std::vector <double>* chi_dzPhotToDimuonVtx = 0; //z distance of photon to dimuon vertex when dxy is minimal
std::vector <double>* chi_dxyPhotToDimuonVtx = 0; //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions
//std::vector <pat::CompositeCandidate>* chiStored = 0;
std::vector <int>* chi_kinematicRefitFlag = 0; // -1 kinematic refit not done, 1 done: +2 needed extra refit for photon +4 something wrong with photon at the end +8 something wrong with the final fit 
std::vector <int>* chi_refit_origChicPosition = 0; //stores position of the chic candidate for the refit (there will be gaps if refit fails) 
std::vector <pat::CompositeCandidate>* chi_refitStored = 0;
std::vector <double>* chi_refit_vprob = 0;
std::vector <double>* chi_refit_ctauPV = 0;
std::vector <double>* chi_refit_ctauErrPV = 0;
std::vector <double>* chi_refit_ctauPV3D = 0;
std::vector <double>* chi_refit_pvtxFromPVwithMuons_x = 0;
std::vector <double>* chi_refit_pvtxFromPVwithMuons_y = 0;
std::vector <double>* chi_refit_pvtxFromPVwithMuons_z = 0;






/////////////////////////
/////// functions ///////
//////////////////////////

int LoadChiBranches(TTree* tree, bool isMC);


#endif 

