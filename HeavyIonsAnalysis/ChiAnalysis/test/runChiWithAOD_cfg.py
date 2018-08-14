#This example can be run over files from AOD, therefore we need to build some information in fly.
#
outFileName = 'Chi_c_pPb8TeV_test7.root'
inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/0249A3C5-A2B1-E611-8E3E-FA163ED701FA.root'

import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("ChiAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v15', '')

#process.GlobalTag.toGet = cms.VPSet(
#  cms.PSet(
#    record = cms.string("HeavyIonRcd"),
#    tag = cms.string("CentralityTable_HFtowersPlusTrunc200_EPOS8TeV_v80x01_mc"),
#    connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#    label = cms.untracked.string("HFtowersPlusTruncEpos")
#    )
#  )

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inFileNames))
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


# we will select muons 
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.ChiPATMuons = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
    muonSource = 'muons',
    #useParticleFlow = cms.bool (True),
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    embedPFCandidate    = False,
    embedCaloMETMuonCorrs = cms.bool(False),
    embedTcMETMuonCorrs   = cms.bool(False),
    embedPfEcalEnergy     = cms.bool(False),
    embedPickyMuon = False,
    embedTpfmsMuon = False,
    userIsolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(),     # no heavy isodeposits
    addGenMatch = False,          # no mc
)

# cuts on muons
process.ChiSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('ChiPATMuons'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 1.5))'
   ),
   filter = cms.bool(False)
)


#DIMUONS 

process.load("HeavyIonsAnalysis.HiOnia2MuMu.HiOnia2MuMuPAT_cfi")
process.HiOnia2MuMuPAT.muons=cms.InputTag('ChiSelectedMuons')
process.HiOnia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
process.HiOnia2MuMuPAT.higherPuritySelection = cms.string("isGlobalMuon") #O "isGlobalMuon"
process.HiOnia2MuMuPAT.lowerPuritySelection = cms.string("isTrackerMuon") #O "isGlobalMuon"
process.HiOnia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.HiOnia2MuMuPAT.dimuonSelection=cms.string("0.2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
process.HiOnia2MuMuPAT.addMCTruth = cms.bool(False)
process.HiOnia2MuMuPAT.addMuonlessPrimaryVertex = cms.bool(False)



# PHOTONS

#mc matching done by hand in rootupler - because it doesn't work well with OniaPhotonConversionProducer

# The low energy photons are reconstructed here.
import HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi
process.PhotonCandidates = HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi.PhotonCandidates.clone(
    conversions = 'allConversions',
    convAlgo    = 'undefined',
    convQuality = [''], #O: Changed ['highPurity','generalTracksOnly']
    primaryVertexTag = 'offlinePrimaryVertices',
    convSelection = 'conversionVertex.position.rho>0.0', #O: Changed 1.5
    wantTkVtxCompatibility = False,
    sigmaTkVtxComp = 50, #O: Changed 5
    wantCompatibleInnerHits = True,
    pfcandidates = 'particleFlow',
    pi0OnlineSwitch = False,
    TkMinNumOfDOF = 0, #O: Changed 3
    wantHighpurity = False,
    #test
    vertexChi2ProbCut = 0.0000,
    trackchi2Cut = 1000,
    minDistanceOfApproachMinCut = -100.25,
    minDistanceOfApproachMaxCut = 100.00,
    #O: Original
    #vertexChi2ProbCut = 0.0005,
    #trackchi2Cut = 10,
    #minDistanceOfApproachMinCut = -0.25,
    #minDistanceOfApproachMaxCut = 1.00,
)


# Chi parts
process.load('HeavyIonsAnalysis.ChiAnalysis.ChiAnalyzer_cfi')
process.ChiRootuple.muon_cand=cms.InputTag('ChiPATMuons')

process.analysisPath = cms.Path(
            process.ChiPATMuons *
            process.ChiSelectedMuons * 
            process.HiOnia2MuMuPAT *
            process.PhotonCandidates *
            process.ChiSequence
)



# muons without trigger info, alternatively in recent version of 80x you can use muon with trigger as well.
## we will select muons and create Onia2MuMu pairs.
#import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
#process.oniaPATMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
#    muonSource = 'muons',
#    embedTrack          = True,
#    embedCombinedMuon   = True,
#    embedStandAloneMuon = True,
#    embedPFCandidate    = False,
#    embedCaloMETMuonCorrs = cms.bool(False),
#    embedTcMETMuonCorrs   = cms.bool(False),
#    embedPfEcalEnergy     = cms.bool(False),
#    embedPickyMuon = False,
#    embedTpfmsMuon = False,
#    userIsolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
#    isoDeposits = cms.PSet(),     # no heavy isodeposits
#    addGenMatch = False,          # no mc
#)

#process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
#   src = cms.InputTag('oniaPATMuonsWithoutTrigger'),
#   cut = cms.string('muonID(\"TMOneStationTight\")'
#                    #' && abs(innerTrack.dxy) < 0.3'
#                    #' && abs(innerTrack.dz)  < 20.'
#                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
#                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
#                    ' && innerTrack.quality(\"highPurity\")'
#                    ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 1.5))'
#   ),
#   filter = cms.bool(True)
#)

#process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
#process.onia2MuMuPAT.muons=cms.InputTag('oniaSelectedMuons')
#process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
#process.onia2MuMuPAT.higherPuritySelection = cms.string("isGlobalMuon") #O "isGlobalMuon"
#process.onia2MuMuPAT.lowerPuritySelection = cms.string("isTrackerMuon") #O "isGlobalMuon"
#process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
#process.onia2MuMuPAT.dimuonSelection=cms.string("0.2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
#process.onia2MuMuPAT.addMCTruth = cms.bool(False)

## make photon candidate conversions for P-wave studies.
## The low energy photons are reconstructed here.
#import HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi
#process.oniaPhotonCandidates = HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi.PhotonCandidates.clone(
#    conversions = 'allConversions',
#    convAlgo    = 'undefined',
#    convQuality = [''], #O: Changed ['highPurity','generalTracksOnly']
#    primaryVertexTag = 'offlinePrimaryVertices',
#    convSelection = 'conversionVertex.position.rho>0.0', #O: Changed 1.5
#    wantTkVtxCompatibility = False,
#    sigmaTkVtxComp = 50, #O: Changed 5
#    wantCompatibleInnerHits = True,
#    pfcandidates = 'particleFlow',
#    pi0OnlineSwitch = False,
#    TkMinNumOfDOF = 0, #O: Changed 3
#    wantHighpurity = False,
#    #test
#    vertexChi2ProbCut = 0.0000,
#    trackchi2Cut = 1000,
#    minDistanceOfApproachMinCut = -100.25,
#    minDistanceOfApproachMaxCut = 100.00,
#    #O: Original
#    #vertexChi2ProbCut = 0.0005,
#    #trackchi2Cut = 10,
#    #minDistanceOfApproachMinCut = -0.25,
#    #minDistanceOfApproachMaxCut = 1.00,
#)

## Centrality and other
## process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
## from Configuration.AlCa.GlobalTag import GlobalTag
## process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v10', '')

#process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
#process.centralityBin.Centrality = cms.InputTag("pACentrality")
#process.centralityBin.centralityVariable = cms.string("HFtowersPlusTrunc")
#process.centralityBin.nonDefaultGlauberModel = cms.string("Epos")

## with all pieces at hand, we procced in the standard way.
#process.load('Ponia.OniaPhoton.OniaPhotonProducer_cfi')
#process.load('Ponia.OniaPhoton.OniaPhotonRootupler_cfi')

#tag_dimuon="onia2MuMuPAT"
#tag_chi_conv_prod="oniaPhotonCandidates"
#tag_chi_conv_lab="conversions"
#process.DiMuonCounter.src = cms.InputTag(tag_dimuon)
#process.PhotonCounter.src  = cms.InputTag(tag_chi_conv_prod,tag_chi_conv_lab)
#process.ChiProducer.conversions = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab)

#process.p = cms.Path(
#            process.oniaPATMuonsWithoutTrigger *
#            process.oniaSelectedMuons *
#            process.onia2MuMuPAT *
#            process.centralityBin *
#            process.oniaPhotonCandidates *
#            process.ChiSequence*process.rootuple
#)
#process.rootuple.isMC = cms.bool(False)
