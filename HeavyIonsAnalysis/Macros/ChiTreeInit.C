#include "ChiTreeInit.h"


/////////////////////////
/////// functions ///////
//////////////////////////

int LoadChiBranches(TTree* tree, bool isMC) {
	if (!tree) {
		cout << "Problem with event tree";
		return 1;
	}



	//general
	tree->SetBranchAddress("runNumber", &runNumber);
	tree->SetBranchAddress("eventNumber", &eventNumber);
	tree->SetBranchAddress("nPrimVertices", &nPrimVertices);
	tree->SetBranchAddress("muonPerEvent_noCuts", &muonPerEvent);
	tree->SetBranchAddress("convPerEvent_noCuts", &convPerTriggeredEvent);
	tree->SetBranchAddress("dimuonPerEvent", &dimuonPerEvent);
	tree->SetBranchAddress("chiCandPerEvent", &chiCandPerEvent);

	tree->SetBranchAddress("ntracks_inEvent", &ntracks_inEvent);
	tree->SetBranchAddress("hfTowerSum_inEvent", &hfTowerSum_inEvent);
	tree->SetBranchAddress("Trig_Event_HLTDoubleMuOpen", &Trig_Event_HLTDoubleMuOpen);


	// vertex
	tree->SetBranchAddress("pvtx_z", &pvtx_z);
	tree->SetBranchAddress("pvtx_zError", &pvtx_zError);
	tree->SetBranchAddress("pvtx_x", &pvtx_x);
	tree->SetBranchAddress("pvtx_y", &pvtx_y);
	tree->SetBranchAddress("pvtx_nTracks", &pvtx_nTracks);
	tree->SetBranchAddress("pvtx_isFake", &pvtx_isFake);

	//muon info
	tree->SetBranchAddress("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
	tree->SetBranchAddress("muonIsHLTDoubleMuOpenFilter", &muonIsHLTDoubleMuOpenFilter);
	tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal);
	tree->SetBranchAddress("muonIsTracker", &muonIsTracker);
	tree->SetBranchAddress("muonIsPF", &muonIsPF);
	tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
	tree->SetBranchAddress("muonIsTight", &muonIsTight);
	tree->SetBranchAddress("muonIsNotGlobalNorTracker", &muonIsNotGlobalNorTracker);
	tree->SetBranchAddress("muonIDHas_TMOneStationTight", &muonIDHas_TMOneStationTight);
	tree->SetBranchAddress("muon_pvtx_index", &muon_pvtx_index);
	tree->SetBranchAddress("muonInnerTrack_dxy", &muonInnerTrack_dxy);
	tree->SetBranchAddress("muonInnerTrack_dz", &muonInnerTrack_dz);
	tree->SetBranchAddress("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);
	tree->SetBranchAddress("muonPixelLayersWithMeasurement", &muonPixelLayersWithMeasurement);
	tree->SetBranchAddress("muonQuality_isHighPurity", &muonQuality_isHighPurity);
	tree->SetBranchAddress("muon_charge", &muon_charge);
	tree->SetBranchAddress("muon_eta", &muon_eta);
	tree->SetBranchAddress("muon_pt", &muon_pt);
	tree->SetBranchAddress("muon_p4", &muon_p4); //TLorentzVector
	//tree->SetBranchAddress("patMuonStored", &patMuonStored);
	
	//muon MC
	if (isMC == true) {
		tree->SetBranchAddress("muon_isMatchedMC", &muon_isMatchedMC);
		tree->SetBranchAddress("muonGen_eta", &muonGen_eta);
		tree->SetBranchAddress("muonGen_pt", &muonGen_pt);
		tree->SetBranchAddress("muonGen_p4", &muonGen_p4); //TLorentzVector
		tree->SetBranchAddress("muonGen_rDelta", &muonGen_rDelta);
		tree->SetBranchAddress("muonGen_ptDelta", &muonGen_ptDelta);
		tree->SetBranchAddress("muonGen_ptDeltaRel", &muonGen_ptDeltaRel);
	}

	//dimuon info

	tree->SetBranchAddress("dimuon_p4", &dimuon_p4); //TLorentzVector
	tree->SetBranchAddress("dimuon_eta", &dimuon_eta);
	tree->SetBranchAddress("dimuon_pt", &dimuon_pt);
	tree->SetBranchAddress("dimuon_charge", &dimuon_charge);
	tree->SetBranchAddress("dimuon_vtx", &dimuon_vtx); //TVector3
	tree->SetBranchAddress("dimuon_pvtx_indexFromOniaMuMu", &dimuon_pvtx_indexFromOniaMuMu);
	tree->SetBranchAddress("dimuon_pvtx_index", &dimuon_pvtx_index);

	tree->SetBranchAddress("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);
	tree->SetBranchAddress("dimuon_vtxProb", &dimuon_vtxProb);
	//tree->SetBranchAddress("dimuonStored", &dimuonStored);
	tree->SetBranchAddress("dimuon_muon1_position", &dimuon_muon1_position); //stores position of first muon in muon collection (probably one with higher pT due to ordering in the collection)
	tree->SetBranchAddress("dimuon_muon2_position", &dimuon_muon2_position); //stores position of second muon in muon collection 
	tree->SetBranchAddress("dimuon_ctpv", &dimuon_ctpv);
	tree->SetBranchAddress("dimuon_ctpvError", &dimuon_ctpvError);



	//conversion info

	tree->SetBranchAddress("convRaw_duplicityStatus", &convRaw_duplicityStatus);
	tree->SetBranchAddress("convRaw_splitDR", &convRaw_splitDR);
	tree->SetBranchAddress("convRaw_splitDpT", &convRaw_splitDpT);
	tree->SetBranchAddress("conv_positionRaw", &conv_positionRaw);
 
	tree->SetBranchAddress("conv_tk1ValidHits", &conv_tk1ValidHits);
	tree->SetBranchAddress("conv_tk2ValidHits", &conv_tk2ValidHits);
	tree->SetBranchAddress("convQuality_isHighPurity", &convQuality_isHighPurity);
	tree->SetBranchAddress("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly);
	tree->SetBranchAddress("conv_vtx", &conv_vtx); //TVector3
	tree->SetBranchAddress("conv_vertexPositionRho", &conv_vertexPositionRho);
	tree->SetBranchAddress("conv_sigmaTkVtx1", &conv_sigmaTkVtx1);
	tree->SetBranchAddress("conv_sigmaTkVtx2", &conv_sigmaTkVtx2);
	tree->SetBranchAddress("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOK);
	tree->SetBranchAddress("conv_tkVtxCompatible_bestVertex", &conv_tkVtxCompatible_bestVertex);
	tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexA", &conv_tkVtxCompatible_secondBestVertexA);
	tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexB", &conv_tkVtxCompatible_secondBestVertexB);
	//tree->SetBranchAddress("conv_tkVtxCompatibilityOK_test", &conv_tkVtxCompatibilityOK_test); //test ones
	//tree->SetBranchAddress("conv_tkVtxCompatible_bestVertex_test", &conv_tkVtxCompatible_bestVertex_test); //test ones
	//tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexA_test", &conv_tkVtxCompatible_secondBestVertexA_test); //test
	//tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexB_test", &conv_tkVtxCompatible_secondBestVertexB_test); //test

	tree->SetBranchAddress("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK); //-1: less than 2 tracks, 0: not compatible, 1: yes
	//tree->SetBranchAddress("conv_hitPat1", &conv_hitPat1);
	//tree->SetBranchAddress("conv_hitPat2", &conv_hitPat2);
	//tree->SetBranchAddress("conv_isCustomHighPurity", &conv_isCustomHighPurity);//tbd - is just a sum of some other cuts, not creating at the time
	tree->SetBranchAddress("conv_vertexChi2Prob", &conv_vertexChi2Prob);
	tree->SetBranchAddress("conv_pvtx_index", &conv_pvtx_index);
	tree->SetBranchAddress("conv_zOfPriVtx", &conv_zOfPriVtx); // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
	tree->SetBranchAddress("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks);
	tree->SetBranchAddress("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx);
	tree->SetBranchAddress("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1);
	tree->SetBranchAddress("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2);
	tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1);
	tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2);
	tree->SetBranchAddress("conv_dxyError_Tr1", &conv_dxyError_Tr1);
	tree->SetBranchAddress("conv_dxyError_Tr2", &conv_dxyError_Tr2);

	tree->SetBranchAddress("conv_tk1NumOfDOF", &conv_tk1NumOfDOF);
	tree->SetBranchAddress("conv_tk2NumOfDOF", &conv_tk2NumOfDOF);
	tree->SetBranchAddress("conv_track1Chi2", &conv_track1Chi2);
	tree->SetBranchAddress("conv_track2Chi2", &conv_track2Chi2);
	tree->SetBranchAddress("conv_Tr1_pt", &conv_Tr1_pt);
	tree->SetBranchAddress("conv_Tr2_pt", &conv_Tr2_pt);

	tree->SetBranchAddress("conv_minDistanceOfApproach", &conv_minDistanceOfApproach);
	tree->SetBranchAddress("conv_p4", &conv_p4); //TLorentzVector
	tree->SetBranchAddress("conv_eta", &conv_eta);
	tree->SetBranchAddress("conv_pt", &conv_pt);


	//conv MC
	if (isMC == true) {
		tree->SetBranchAddress("conv_isMatchedMC", &conv_isMatchedMC);
		tree->SetBranchAddress("convGen_eta", &convGen_eta);
		tree->SetBranchAddress("convGen_pt", &convGen_pt);
		tree->SetBranchAddress("convGen_p4", &convGen_p4); //TLorentzVector
		tree->SetBranchAddress("convGen_rDelta", &convGen_rDelta);
		tree->SetBranchAddress("convGen_ptDelta", &convGen_ptDelta);
		tree->SetBranchAddress("convGen_ptDeltaRel", &convGen_ptDeltaRel);
		tree->SetBranchAddress("convGen_motherCode", &convGen_motherCode);
	}

	// MC general
	if (isMC == true) {
		tree->SetBranchAddress("gen_isGoodChicDecay", &gen_isGoodChicDecay);
		tree->SetBranchAddress("gen_pdgId", &gen_pdgId);
		tree->SetBranchAddress("gen_chic_pt", &gen_chic_pt);
		tree->SetBranchAddress("gen_chic_eta", &gen_chic_eta);
		tree->SetBranchAddress("gen_chic_p4", &gen_chic_p4); //TLorentzVector
		tree->SetBranchAddress("gen_chic_matchPosition", &gen_chic_matchPosition);
		tree->SetBranchAddress("gen_chic_nMatches", &gen_chic_nMatches);
		tree->SetBranchAddress("gen_chic_rDelta", &gen_chic_rDelta); //in principle duplicates information
		tree->SetBranchAddress("gen_chic_ptDeltaRel", &gen_chic_ptDeltaRel);//in principle duplicates information

		tree->SetBranchAddress("gen_Jpsi_pt", &gen_Jpsi_pt);
		tree->SetBranchAddress("gen_Jpsi_eta", &gen_Jpsi_eta);
		tree->SetBranchAddress("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
		tree->SetBranchAddress("gen_Jpsi_nMatches", &gen_Jpsi_nMatches);
		tree->SetBranchAddress("gen_Jpsi_rDelta", &gen_Jpsi_rDelta); //in principle duplicates information
		tree->SetBranchAddress("gen_Jpsi_ptDeltaRel", &gen_Jpsi_ptDeltaRel);//in principle duplicates information
		tree->SetBranchAddress("gen_Jpsi_p4", &gen_Jpsi_p4); //TLorentzVector

		tree->SetBranchAddress("gen_Jpsi_photon_n", &gen_Jpsi_photon_n);
		tree->SetBranchAddress("gen_Jpsi_photon_pt", &gen_Jpsi_photon_pt);
		tree->SetBranchAddress("gen_Jpsi_photon_p4", &gen_Jpsi_photon_p4); //TLorentzVector

		tree->SetBranchAddress("gen_muon_charge", &gen_muon_charge);
		tree->SetBranchAddress("gen_muon_pt", &gen_muon_pt);
		tree->SetBranchAddress("gen_muon_eta", &gen_muon_eta);
		tree->SetBranchAddress("gen_muon_matchPosition", &gen_muon_matchPosition);
		tree->SetBranchAddress("gen_muon_nMatches", &gen_muon_nMatches);
		tree->SetBranchAddress("gen_muon_rDelta", &gen_muon_rDelta); //in principle duplicates information
		tree->SetBranchAddress("gen_muon_ptDeltaRel", &gen_muon_ptDeltaRel);//in principle duplicates information
		tree->SetBranchAddress("gen_muon_p4", &gen_muon_p4); //TLorentzVector

		tree->SetBranchAddress("gen_phot_pt", &gen_phot_pt);
		tree->SetBranchAddress("gen_phot_eta", &gen_phot_eta);
		tree->SetBranchAddress("gen_phot_p4", &gen_phot_p4); //TLorentzVector
		tree->SetBranchAddress("gen_conv_matchPosition", &gen_conv_matchPosition);
		tree->SetBranchAddress("gen_conv_nMatches", &gen_conv_nMatches);
		tree->SetBranchAddress("gen_conv_rDelta", &gen_conv_rDelta); //in principle duplicates information
		tree->SetBranchAddress("gen_conv_ptDeltaRel", &gen_conv_ptDeltaRel);//in principle duplicates information
	}

	// chi

	tree->SetBranchAddress("chi_p4", &chi_p4); //TLorentzVector
	tree->SetBranchAddress("chi_eta", &chi_eta);
	tree->SetBranchAddress("chi_pt", &chi_pt);
	tree->SetBranchAddress("chi_daughterJpsi_position", &chi_daughterJpsi_position); //stores position of daughter Jpsi in dimuon collection
	tree->SetBranchAddress("chi_daughterConv_position", &chi_daughterConv_position); //stores position of daughter photon (conversion)
	tree->SetBranchAddress("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx); //z distance of photon to dimuon vertex when dxy is minimal
	tree->SetBranchAddress("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx); //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions
	//tree->SetBranchAddress("chiStored", &chiStored);
	tree->SetBranchAddress("chi_kinematicRefitFlag", &chi_kinematicRefitFlag); // -1 kinematic refit not done, 1 done: +2 needed extra refit for photon +4 something wrong with photon at the end +8 something wrong with the final fit 
	tree->SetBranchAddress("chi_refit_origChicPosition", &chi_refit_origChicPosition); //stores position of the chic candidate for the refit (there will be gaps if refit fails) 
	tree->SetBranchAddress("chi_refitStored", &chi_refitStored);
	tree->SetBranchAddress("chi_refit_vprob", &chi_refit_vprob);
	tree->SetBranchAddress("chi_refit_ctauPV", &chi_refit_ctauPV);
	tree->SetBranchAddress("chi_refit_ctauErrPV", &chi_refit_ctauErrPV);
	tree->SetBranchAddress("chi_refit_ctauPV3D", &chi_refit_ctauPV3D);
	tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_x", &chi_refit_pvtxFromPVwithMuons_x);
	tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_y", &chi_refit_pvtxFromPVwithMuons_y);
	tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_z", &chi_refit_pvtxFromPVwithMuons_z);


	return 0;
}

//muon  - acceptance

bool MuonAcceptanceLoose(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 1.1 && pt < 3.3) return false;
	if (fabs(eta) >= 1.1 && fabs(eta) < 2.1 && pt < 5.5 - 2 * fabs(eta)) return false;
	if (fabs(eta) >= 2.1 && pt < 1.3) return false;
	return true;
}
bool MuonAcceptanceTight(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 2.4 && pt < 3.3) return false;
	return true;
}

bool MuonAcceptance(double eta, double pt) {
	return MuonAcceptanceLoose(eta, pt);
}
// muon selection

bool MuonSelectionPass(int muonPos)  //uses variables loaded in main function
{
	if (muonIsSoft->at(muonPos) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(muonPos) != 1) return false;
	return true;
}

bool MuonSelectionPassMC(int muonPos)  //uses variables loaded in main function
{
	int matchPosition = gen_muon_matchPosition->at(muonPos);
	if (matchPosition < -0.5) return false; //not matched
	if (gen_muon_rDelta->at(muonPos)> muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(muonPos) > muon_maxDPtRel_analysis ) return false;//not matched
	if (muonIsSoft->at(matchPosition) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(matchPosition) != 1) return false;
	return true;
}


// conversion

bool PhotAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false; //2.5
	if (pt < 0.5) return false; // 0.2
	return true;
}



bool PhotSelectionPassTight(int photPos)  //uses variables loaded in main function //nominal by Alberto and Jhovanny (8.2021) (with updated conversion removal)
{
	//if (conv_duplicityStatus->at(photPos) != 0 && conv_duplicityStatus->at(photPos) != 1 && conv_duplicityStatus->at(photPos) != 2) return false;
	if (convQuality_isHighPurity->at(photPos) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	if (conv_sigmaTkVtx1->at(photPos) > 5) return false;
	if (conv_sigmaTkVtx2->at(photPos) > 5) return false;
	//if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (conv_vertexChi2Prob->at(photPos) <= 0.0005) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	//if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_track1Chi2->at(photPos) >= 10) return false;
	if (conv_track2Chi2->at(photPos) >= 10) return false;
	if (conv_minDistanceOfApproach->at(photPos) <= -0.25) return false; //removes -1000 values
	if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}

bool PhotSelectionPassMedium(int photPos)  //uses variables loaded in main function
{
	//if (conv_duplicityStatus->at(photPos) != 0 && conv_duplicityStatus->at(photPos) != 1 && conv_duplicityStatus->at(photPos) != 2) return false;
	//if (convQuality_isHighPurity->at(photPos) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	//if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	//if (conv_sigmaTkVtx1->at(photPos) > 10) return false;
	//if (conv_sigmaTkVtx2->at(photPos) > 10) return false;
	//if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (conv_vertexChi2Prob->at(photPos) <= 0.0005) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_track1Chi2->at(photPos) >= 10) return false;
	//if (conv_track2Chi2->at(photPos) >= 10) return false;
	if (conv_minDistanceOfApproach->at(photPos) <= -10) return false; //removes -1000 values
	//if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}

bool PhotSelectionPassLoose(int photPos)  //uses variables loaded in main function
{
	//if (conv_duplicityStatus->at(photPos) != 0 && conv_duplicityStatus->at(photPos) != 1 && conv_duplicityStatus->at(photPos) != 2) return false;
	//if (convQuality_isHighPurity->at(photPos) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	//if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	//if (conv_sigmaTkVtx1->at(photPos) > 10) return false;
	//if (conv_sigmaTkVtx2->at(photPos) > 10) return false;
	//if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	//if (conv_vertexChi2Prob->at(photPos) <= 0.001) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	//if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_track1Chi2->at(photPos) >= 10) return false;
	//if (conv_track2Chi2->at(photPos) >= 10) return false;
	if (conv_minDistanceOfApproach->at(photPos) <= -10) return false; //removes -1000 values
	//if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}

bool PhotSelectionPass(int photPos) // nominal version
{
	return PhotSelectionPassLoose(photPos);
}

bool PhotSelectionPassMC(int photPos)  //uses variables loaded in main function
{
	int matchPosition = gen_conv_matchPosition->at(photPos);
	if (matchPosition < -0.5) return false; //not matched
	if (gen_conv_rDelta->at(photPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(photPos) > conv_maxDPtRel_analysis) return false;//not matched
	return PhotSelectionPass(matchPosition);
}

bool PhotSelectionPassMCLoose(int photPos)  //uses variables loaded in main function
{
	int matchPosition = gen_conv_matchPosition->at(photPos);
	if (matchPosition < -0.5) return false; //not matched
	if (convQuality_isGeneralTracksOnly->at(matchPosition) != 1) return false;
	return true;
}

/////////////////////
//  D I M U O N /////
/////////////////////

// dimuon - acceptance

bool DimuonAcceptanceLoose(double rap, double pt)
{
	if (fabs(rap) > 2.4) return false;
	if (pt < 6.5) return false;
	if (pt > 30.0) return false;
	return true;
}
bool DimuonAcceptanceTight(double rap, double pt)
{
	if (fabs(rap) > 1.0) return false;  //2.4
	if (pt < 6.5) return false;
	if (pt > 25) return false;
	return true;
}

bool DimuonAcceptance(double rap, double pt) {
	return DimuonAcceptanceLoose(rap, pt);
}



// dimuon selection

bool DimuonSelectionPass(int dimuonPos)  //uses variables loaded in main function // contains acceptance cut for historic reasons
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	double ctauSigJpsi = dimuon_ctpv->at(dimuonPos) / dimuon_ctpvError->at(dimuonPos);
	if (ctauSigJpsi > 3 || ctauSigJpsi < -3) return false;
	double rap = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity();
	return DimuonAcceptance(rap, dimuon_pt->at(dimuonPos));
}

bool DimuonSelectionPassTight(int dimuonPos)  //uses variables loaded in main function // contains acceptance cut for historic reasons
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	double ctauSigJpsi = dimuon_ctpv->at(dimuonPos) / dimuon_ctpvError->at(dimuonPos);
	if (ctauSigJpsi > 3 || ctauSigJpsi < -3) return false;
	double rap = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity();
	return DimuonAcceptanceTight(rap, dimuon_pt->at(dimuonPos));
}

bool DimuonSelectionPassMC(int dimuonMCPos)
{
	int matchPosition = DimuonMCMatched(dimuonMCPos);
	if (matchPosition < -0.5) return false; //not matched
	return DimuonSelectionPass(matchPosition);
}

bool DimuonSelectionPassNoCharge(int dimuonPos)  //to allow selecting same sign bkg (but no sign cut directly)
{
	//if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	double ctauSigJpsi = dimuon_ctpv->at(dimuonPos) / dimuon_ctpvError->at(dimuonPos);
	if (ctauSigJpsi > 3 || ctauSigJpsi < -3) return false;
	double rap = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity();
	return DimuonAcceptance(rap, dimuon_pt->at(dimuonPos));
}

int DimuonMCMatched(int dimuonMCPos = 0)// check if the gen jpsi was matched to reco. Usually one Jpsi per event (thus index 0). New version - matching to the muons only: 
//-1 no muon match, -2 muons matched, but the dimuon doesn't exist (probably removed by dimuon preselection - confirmed for most)
{
	int dimuonPos = -2; //set to default (no match)

	int iMuon1MC = gen_muon_matchPosition->at(2 * dimuonMCPos);
	int iMuon2MC = gen_muon_matchPosition->at(2 * dimuonMCPos + 1);

	if (iMuon1MC < 0 || gen_muon_rDelta->at(2 * dimuonMCPos) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * dimuonMCPos) > muon_maxDPtRel_analysis) return -1; //muon not matched
	if (iMuon2MC < 0 || gen_muon_rDelta->at(2 * dimuonMCPos + 1) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * dimuonMCPos + 1) > muon_maxDPtRel_analysis) return -1; //muon not matched

	// all the final objects were matched, let's find the chic that is the proper match 


	for (int iJpsi = 0; iJpsi < dimuon_p4->GetEntriesFast(); iJpsi++) // Jpsi loop
	{
		if (dimuon_muon1_position->at(iJpsi) != iMuon1MC && (dimuon_muon1_position->at(iJpsi) != iMuon2MC)) { continue; } //first muon don't agree with either gen
		if (dimuon_muon2_position->at(iJpsi) != iMuon1MC && (dimuon_muon2_position->at(iJpsi) != iMuon2MC)) { continue; } //second muon don't agree with either gen
		if (dimuonPos < 0) { dimuonPos = iJpsi; } //we found the one matched
		else { cout << "Something wrong, two Jpsi shouldn't be both matched" << endl; } //we already found one, this is a crosscheck (once tested, add break above)
	}
	/*if (dimuonPos == -1) {
		cout << "Something odd, all the objects were matched, but we didn't find the matching dimuon. Muon 1 and 2 acceptance and selection cuts: " << MuonAcceptance(muon_eta->at(iMuon1MC), muon_pt->at(iMuon1MC))<< " " << MuonAcceptance(muon_eta->at(iMuon2MC), muon_pt->at(iMuon2MC)) << "  sel: " << MuonSelectionPassMC(2 * dimuonMCPos) << " " << MuonSelectionPassMC(2 * dimuonMCPos+1) << " Tracker: " << muonIsTracker->at(iMuon1MC)<< " " << muonIsTracker->at(iMuon2MC) << endl;
		cout << "n dimuons " << dimuon_p4->GetEntriesFast() << endl;
		if (dimuon_p4->GetEntriesFast() > 0)cout << "Reco muon pTs: " << muon_pt->at(dimuon_muon1_position->at(0)) << " " << muon_pt->at(dimuon_muon2_position->at(0)) << " etas: " << muon_eta->at(dimuon_muon1_position->at(0)) << " " << muon_eta->at(dimuon_muon2_position->at(0)) <<endl;
		cout << "Maybe different vertices?: " << muon_pvtx_index->at(iMuon1MC) << "  " << muon_pvtx_index->at(iMuon2MC) << endl;
		TLorentzVector* LVmuon1 = (TLorentzVector*)muon_p4->At(iMuon1MC);
		TLorentzVector* LVmuon2 = (TLorentzVector*)muon_p4->At(iMuon2MC);
		TLorentzVector* LVmuon1Gen = (TLorentzVector*)gen_muon_p4->At(2 * dimuonMCPos);
		TLorentzVector* LVmuon2Gen = (TLorentzVector*)gen_muon_p4->At(2 * dimuonMCPos +1);
		cout << "Not enough pt?: " << (*LVmuon1+*LVmuon2).Pt() << " rap: " << (*LVmuon1 + *LVmuon2).Rapidity() << " mass: " << (*LVmuon1 + *LVmuon2).M() <<  endl;
		cout << "Delta r: " << gen_muon_rDelta->at(2 * dimuonMCPos) << " " << gen_muon_rDelta->at(2 * dimuonMCPos + 1) << " and dpT: " << gen_muon_ptDeltaRel->at(2 * dimuonMCPos) << " " << gen_muon_ptDeltaRel->at(2 * dimuonMCPos + 1) << endl;
		cout << "Gen Mass: " << ((TLorentzVector*)gen_Jpsi_p4->At(dimuonMCPos))->M() << " and muon pT: " << gen_muon_pt->at(2 * dimuonMCPos) << " " << gen_muon_pt->at(2 * dimuonMCPos+1) << "  compared to reco muon mass: " << (*LVmuon1 + *LVmuon2).M() << " and pt " << muon_pt->at(iMuon1MC) << " " << muon_pt->at(iMuon2MC) <<endl;
		cout << "Gen muon eta: " << gen_muon_eta->at(2 * dimuonMCPos) << " " << gen_muon_eta->at(2 * dimuonMCPos + 1) << " and reco eta: " << muon_eta->at(iMuon1MC) << " " << muon_eta->at(iMuon2MC) << endl;
		cout << "Gen Mass: " << ((TLorentzVector*)gen_Jpsi_p4->At(dimuonMCPos))->M() << " and calculated from decay muons, gen information: " << (*LVmuon1Gen + *LVmuon2Gen).M() << endl;
		cout << "is good gen decay: " << gen_isGoodChicDecay->at(dimuonMCPos) << " n photons: " << gen_Jpsi_photon_n->at(dimuonMCPos) << endl;
		if (gen_Jpsi_photon_n->at(dimuonMCPos)>0){
			TLorentzVector* LVphot1Gen = (TLorentzVector*)gen_Jpsi_photon_p4->At(0);
			cout<< "Photon pT: " << gen_Jpsi_photon_pt->at(0) << ", now we add the photon to the mass: " << (*LVmuon1Gen + *LVmuon2Gen + *LVphot1Gen).M() << endl;
		}
		cout << endl;
		return -2;
	}//*/
	return dimuonPos; //either matched, or -2 if such reco doesn't exist
}

// overall dimuon pass

bool DimuonPassAllCuts(int dimuonPos)
{
	if (DimuonSelectionPass(dimuonPos) == false) return false;
	int muon1Pos = dimuon_muon1_position->at(dimuonPos);
	int muon2Pos = dimuon_muon2_position->at(dimuonPos);
	if (MuonSelectionPass(muon1Pos) == false) return false;
	if (MuonSelectionPass(muon2Pos) == false) return false;
	if (MuonAcceptance(muon_eta->at(muon1Pos), muon_pt->at(muon1Pos)) == false) return false;
	if (MuonAcceptance(muon_eta->at(muon2Pos), muon_pt->at(muon2Pos)) == false) return false;
	return true;
}


int DimuonPassAllCutsMC(int dimuonMCPos = 0)  // -1 if failed, else returns the position of good reco
{
	int dimuonPos = DimuonMCMatched(dimuonMCPos);
	if (dimuonPos < 0) { return -1; }
	else if (DimuonPassAllCuts(dimuonPos) == false) { return -1; }
	else return dimuonPos;
}




/////////////////////
////  C H I C   /////
/////////////////////


bool ChiPassAllCuts(int chiPos)
{
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (DimuonPassAllCuts(dimuonPos) == false) return false;

	// photon
	int convPos = chi_daughterConv_position->at(chiPos);
	if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) return false;

	if (PhotSelectionPass(convPos) == false) return false;
	return true;

//	//cuts or MVA for photons   //NO TMVA FOR NOW
//#ifdef UsesTMVA
//	convQuality_isHighPurityValue = convQuality_isHighPurity->at(convPos);
//	convQuality_isGeneralTracksOnlyValue = convQuality_isGeneralTracksOnly->at(convPos);
//	conv_vertexPositionRhoValue = conv_vertexPositionRho->at(convPos);
//	conv_sigmaTkVtx1Value = conv_sigmaTkVtx1->at(convPos);
//	conv_sigmaTkVtx2Value = conv_sigmaTkVtx2->at(convPos);
//	conv_tkVtxCompatibilityOKValue = conv_tkVtxCompatibilityOK->at(convPos);
//	conv_compatibleInnerHitsOKValue = conv_compatibleInnerHitsOK->at(convPos);
//	conv_vertexChi2ProbValue = conv_vertexChi2Prob->at(convPos);
//	conv_dzToClosestPriVtxValue = conv_dzToClosestPriVtx->at(convPos);
//	conv_dxyPriVtxTimesCharge_Tr1Value = conv_dxyPriVtxTimesCharge_Tr1->at(convPos);
//	conv_dxyPriVtxTimesCharge_Tr2Value = conv_dxyPriVtxTimesCharge_Tr2->at(convPos);
//	conv_tk1NumOfDOFValue = conv_tk1NumOfDOF->at(convPos);
//	conv_tk2NumOfDOFValue = conv_tk2NumOfDOF->at(convPos);
//	conv_minDistanceOfApproachValue = conv_minDistanceOfApproach->at(convPos);
//	conv_etaValue = conv_eta->at(convPos);
//	conv_ptValue = conv_pt->at(convPos);
//
//
//	double photMVA = TMWAreader->EvaluateMVA("BDT");
//	hphoton_MVA_response->Fill(photMVA);
//	//cout << "Values MVA: HP: " << convQuality_isHighPurityValue << "   TracksOnly: " << convQuality_isGeneralTracksOnlyValue << "    Vprob: " << conv_vertexChi2ProbValue << "    and the response: "<< photMVA << endl;
//	//if (photMVA < -0.2) continue;
//#endif
//


}


int ChiMCMatched(int chiMCPos = 0)// check if the gen chic was matched to reco. Usually one chic per event (thus index 0). New version - matching to the muons and conversions only
//-1 no matches, -2 conversions and muons matched, but the chic doesn't exist (probably removed by dimuon preselection - confirmed for most, but could maybe rarely be caused by chic selection)
{
	int chiPos = -2; //set to default (no match)
	
	// check if the matching to the gen objects happened (this we can get directly from the stored information, without checking the actual matched reco object)
	if (gen_conv_matchPosition->at(chiMCPos) < 0 || gen_conv_rDelta->at(chiMCPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(chiMCPos) > conv_maxDPtRel_analysis) return -1; //conv was not matched
	if (gen_muon_matchPosition->at(2 * chiMCPos) < 0 || gen_muon_rDelta->at(2 * chiMCPos) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos) > muon_maxDPtRel_analysis) return -1; //muon not matched
	if (gen_muon_matchPosition->at(2 * chiMCPos + 1) < 0 || gen_muon_rDelta->at(2 * chiMCPos + 1) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos + 1) > muon_maxDPtRel_analysis) return -1; //muon not matched

	// all the final objects were matched, let's find the chic that is the proper match 


	for (int iChi = 0; iChi < chi_p4->GetEntriesFast(); iChi++) //loop over chic to find the one matched
	{
		if (chi_daughterConv_position->at(iChi) != gen_conv_matchPosition->at(chiMCPos)) { continue; }
		int dimuonPos = chi_daughterJpsi_position->at(iChi);
		int iMuon1MC = gen_muon_matchPosition->at(2 * chiMCPos);
		int iMuon2MC = gen_muon_matchPosition->at(2 * chiMCPos + 1);

		if (dimuon_muon1_position->at(dimuonPos) != iMuon1MC && (dimuon_muon1_position->at(dimuonPos) != iMuon2MC)) { continue; } //first muon don't agree with either gen
		if (dimuon_muon2_position->at(dimuonPos) != iMuon1MC && (dimuon_muon2_position->at(dimuonPos) != iMuon2MC)) { continue; } //second muon don't agree with either gen
		if (chiPos < 0) { chiPos = iChi; } //we found the one matched
		else { cout << "Something wrong, two chic shouldn't be both matched" << endl; } //we already found one, this is a crosscheck (once tested, add break above)
	}

	return chiPos; //either matched, or -2 if such reco doesn't exist
}



int ChiPassAllCutsMC(int chiMCPos = 0)
{
	int chiPos = ChiMCMatched(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (chiPos < -0.5) { return -1; }// -1 not matched
	else if (ChiPassAllCuts(chiPos) == false) { return -1; }
	else return chiPos;

}






bool ChiIsMatchedAllDaughters(int chiPos, int chiMCPos = 0)// check if the gen chic was matched to reco, as well as J/psi and photon. Usually one chic per event (thus index 0)
// mostly obsolete
{
	int iChiMC = gen_chic_matchPosition->at(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (iChiMC < -0.5) { return false; }// -1 not matched
	if (chiPos != iChiMC) { return false; }
	if (chi_daughterConv_position->at(chiPos) != gen_conv_matchPosition->at(chiMCPos)) { return false; } //conversions don't agree
	if (gen_conv_matchPosition->at(chiMCPos)<0 || gen_conv_rDelta->at(chiMCPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(chiMCPos) > conv_maxDPtRel_analysis) return false;//not matched
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (dimuonPos != gen_Jpsi_matchPosition->at(chiMCPos)) { return false; } //dimuons don't agree
	if (gen_Jpsi_matchPosition->at(chiMCPos) < 0 || gen_Jpsi_rDelta->at(chiMCPos) > jpsi_maxDeltaR_analysis || gen_Jpsi_ptDeltaRel->at(chiMCPos) > jpsi_maxDPtRel_analysis) return false;//not matched
	int iMuon1MC = gen_muon_matchPosition->at(2 * chiMCPos);
	int iMuon2MC = gen_muon_matchPosition->at(2 * chiMCPos + 1);
	if (dimuon_muon1_position->at(dimuonPos) != iMuon1MC && (dimuon_muon1_position->at(dimuonPos) != iMuon2MC)) { return false; } //first muon don't agree with either gen
	if (dimuon_muon2_position->at(dimuonPos) != iMuon1MC && (dimuon_muon2_position->at(dimuonPos) != iMuon2MC)) { return false; } //second muon don't agree with either gen
	if (gen_muon_matchPosition->at(2 * chiMCPos) < 0 || gen_muon_rDelta->at(2 * chiMCPos) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos) > muon_maxDPtRel_analysis) return false;//not matched
	if (gen_muon_matchPosition->at(2 * chiMCPos +1) < 0 || gen_muon_rDelta->at(2 * chiMCPos +1) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos +1) > muon_maxDPtRel_analysis) return false;//not matched

	return true;
}











