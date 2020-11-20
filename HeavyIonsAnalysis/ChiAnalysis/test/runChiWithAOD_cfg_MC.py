#This example can be run over files from AOD, therefore we need to build some information in fly.
#
outFileName = 'Chi_c_pPb8TeV_MC.root'
inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/ChiCJpsiMuMu_Pythia8_8p16TeV_TuneCUETP8M1_RECO_test4.root'

import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("ChiAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_pA_v4', '')

#process.GlobalTag.toGet = cms.VPSet(
#  cms.PSet(
#    record = cms.string("HeavyIonRcd"),
#    tag = cms.string("CentralityTable_HFtowersPlusTrunc200_EPOS8TeV_v80x01_mc"),
#    connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#    label = cms.untracked.string("HFtowersPlusTruncEpos")
#    )
#  )

process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inFileNames))
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))



## ==== Filters ====  from T&P
### pPb Event Selection
process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')
process.primaryVertexFilterPA = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 50 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

### Trigger selection
process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
#process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_PAL1DoubleMu*_v*')
process.triggerResultsFilter.triggerConditions = cms.vstring('*')
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")
process.triggerResultsFilter.l1tResults = cms.InputTag("gtStage2Digis") #O: needs to be this
process.triggerResultsFilter.throw = False
### Filter sequence
process.fastFilter = cms.Sequence(process.hfCoincFilter + process.primaryVertexFilterPA + process.noScraping + process.triggerResultsFilter)


# MUONS

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
## with some customization  O: values copied over from the HI onia trees (or T&P - same)
process.muonL1Info.maxDeltaR = 0.3
process.muonL1Info.maxDeltaEta = 0.2
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR = 0.3
process.muonMatchHLTL1.maxDeltaEta = 0.2
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
## For trigger muons
#switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon

## For L1 muons
addHLTL1Passthrough(process)
useL1Stage2Candidates(process)
process.patTrigger.collections.remove("hltL1extraParticles")
process.patTrigger.collections.append("hltGmtStage2Digis:Muon")
process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGmtStage2Digis:Muon")')
process.muonMatchHLTL1.useStage2L1 = cms.bool(True)
process.muonMatchHLTL1.useMB2InOverlap = cms.bool(True)
process.muonMatchHLTL1.preselection = cms.string("")
appendL1MatchingAlgo(process)


# cuts on muons
process.ChiSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('patMuonsWithTrigger'),
   cut = cms.string(''
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
process.HiOnia2MuMuPAT.dimuonSelection=cms.string("0.2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 50")
process.HiOnia2MuMuPAT.addMCTruth = cms.bool(True)
process.HiOnia2MuMuPAT.addCommonVertex = cms.bool(True)
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
process.ChiRootuple.muon_cand=cms.InputTag('ChiSelectedMuons')
process.ChiRootuple.isMC=cms.bool(True)

process.analysisPath = cms.Path(
            process.fastFilter *
            process.patMuonsWithTriggerSequence *
            process.ChiSelectedMuons * 
            process.HiOnia2MuMuPAT *
            process.PhotonCandidates *
            process.ChiSequence
)

