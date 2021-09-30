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
	tree->SetBranchAddress("dimuon_pvtx_index", &dimuon_pvtx_index);
	tree->SetBranchAddress("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);
	tree->SetBranchAddress("dimuon_vtxProb", &dimuon_vtxProb);
	//tree->SetBranchAddress("dimuonStored", &dimuonStored);
	tree->SetBranchAddress("dimuon_muon1_position", &dimuon_muon1_position); //stores position of first muon in muon collection (probably one with higher pT due to ordering in the collection)
	tree->SetBranchAddress("dimuon_muon2_position", &dimuon_muon2_position); //stores position of second muon in muon collection 
	tree->SetBranchAddress("dimuon_ctpv", &dimuon_ctpv);
	tree->SetBranchAddress("dimuon_ctpvError", &dimuon_ctpvError);



	//conversion info

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



bool MuonAcceptance(double eta, double pt)
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
	return true;
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 2.4 && pt < 3.3) return false;
	return true;
}

bool PhotAcceptance(double eta, double pt)
{
	//if (fabs(eta) > 2.5) return false; //2.5
	if (pt < 0.5) return false; // 0.2
	return true;
}

bool DimuonSelectionPass(int dimuonPos)  //uses variables loaded in main function
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	if (dimuon_pt->at(dimuonPos) < 6) return false;
	return true;
}
bool DimuonSelectionPassMC(int dimuonPos)  //uses variables loaded in main function
{
	int matchPosition = gen_Jpsi_matchPosition->at(dimuonPos);
	if (matchPosition < -0.5) return false; //not matched
	if (dimuon_charge->at(matchPosition) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	if (dimuon_vtxProb->at(matchPosition) < 0.01) return false;
	return true;
}

bool DimuonSelectionPassTight(int dimuonPos, double rap)  //uses variables loaded in main function
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_ctpv->at(dimuonPos) > 10 || dimuon_ctpv->at(dimuonPos) < -10) continue;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	if (fabs(rap) > 1.0) return false;
	if (dimuon_pt->at(dimuonPos) > 25.0) return false;
	if (dimuon_pt->at(dimuonPos) < 6) return false;
	return true;
}

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
	if (muonIsSoft->at(matchPosition) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(matchPosition) != 1) return false;
	return true;
}


bool PhotSelectionPassTight(int photPos)  //uses variables loaded in main function //nominal by Alberto and Jhovanny (8.2021)
{
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
	return PhotSelectionPass(matchPosition);
}

//bool PhotSelectionPassMCLoose(int photPos)  //uses variables loaded in main function
//{
//	int matchPosition = gen_conv_matchPosition->at(photPos);
//	if (matchPosition < -0.5) return false; //not matched
//	if (convQuality_isGeneralTracksOnly->at(matchPosition) != 1) return false;
//	return true;
//}

bool ChiSelectionPassMC(int chiPos)  //uses variables loaded in main function
{
	return true;
}


bool ChiIsMatchedAllDaughters(int chiPos, int chiMCPos = 0)// check if the gen chic was matched to reco, as well as J/psi and photon. Usually one chic per event (thus index 0)
{
	int iChiMC = gen_chic_matchPosition->at(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (iChiMC < -0.5) { return false; }// -1 not matched
	if (chiPos != iChiMC) { return false; }
	if (chi_daughterConv_position->at(chiPos) != gen_conv_matchPosition->at(chiMCPos)) { return false; } //conversions don't agree
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (dimuonPos != gen_Jpsi_matchPosition->at(chiMCPos)) { return false; } //dimuons don't agree
	if (dimuon_muon1_position->at(dimuonPos) != gen_muon_matchPosition->at(2 * chiMCPos) && (dimuon_muon1_position->at(dimuonPos) != gen_muon_matchPosition->at(2 * chiMCPos + 1))) { return false; } //first muon don't agree with either gen
	if (dimuon_muon2_position->at(dimuonPos) != gen_muon_matchPosition->at(2 * chiMCPos) && (dimuon_muon2_position->at(dimuonPos) != gen_muon_matchPosition->at(2 * chiMCPos + 1))) { return false; } //second muon don't agree with either gen
	return true;
}