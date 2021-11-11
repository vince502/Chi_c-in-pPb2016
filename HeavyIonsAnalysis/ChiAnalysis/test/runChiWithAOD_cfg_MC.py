#This example can be run over files from AOD, therefore we need to build some information in fly.
#
outFileName = 'Chi_c_pPb8TeV_MC.root'
#inFileNames = '/store/himc/pPb816Summer16DR/ZZ_PbP-EmbEPOS_8p16_Pythia8/AODSIM/PbPEmb_80X_mcRun2_pA_v4-v1/60000/3A305B07-17F1-E711-8B05-001E677926C0.root'
inFileNames = 'file:/afs/cern.ch/work/o/okukral/ChicMC/CMSSW_8_0_36/src/ChiCJpsiMuMu_Pythia8_8p16TeV_TuneCUETP8M1_Pbp_RECO.root' 
#inFileNames = 'file:/afs/cern.ch/work/o/okukral/ChicMC/CMSSW_8_0_30/src/EPOStest_RECO.root' 
#inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/ChiCJpsiMuMu_Pythia8_8p16TeV_TuneCUETP8M1_RECO_7.root' #v7 MC
#inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/ChiCJpsiMuMu_Pythia8_8p16TeV_TuneCUETP8M1_RECO_3.root' #v6 MC
#inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/ChiCJpsiMuMu_Pythia8_8p16TeV_TuneCUETP8M1_RECO_2.root' #v5 MC
#inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/MCTestFile_v4.root' #v4 MC
#inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/0249A3C5-A2B1-E611-8E3E-FA163ED701FA.root'

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


## cuts on muons
#process.ChiSelectedMuons = cms.EDFilter('PATMuonSelector',
#   src = cms.InputTag('patMuonsWithTrigger'),
#   cut = cms.string('muonID(\"TMOneStationTight\")'
#                    ' && abs(innerTrack.dxy) < 0.3'
#                    ' && abs(innerTrack.dz)  < 20.'
#                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
#                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
#                    ' && innerTrack.quality(\"highPurity\")'
#                    ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 0.7))'
#   ),
#   filter = cms.bool(False)
#)

# cuts on muons
process.ChiSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('patMuonsWithTrigger'),
   cut = cms.string(''#'muonID(\"TMOneStationTight\")'
                   # ' && abs(innerTrack.dxy) < 0.3'
                   # ' && abs(innerTrack.dz)  < 20.'
                   # ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                   # ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                   # ' && innerTrack.quality(\"highPurity\")'
                   # ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 0.7))'
   ),
   filter = cms.bool(False)
)


#DIMUONS 

process.load("HeavyIonsAnalysis.HiOnia2MuMu.HiOnia2MuMuPAT_cfi")
process.HiOnia2MuMuPAT.muons=cms.InputTag('ChiSelectedMuons')
process.HiOnia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
process.HiOnia2MuMuPAT.higherPuritySelection = cms.string("isTrackerMuon") #O "isGlobalMuon"
process.HiOnia2MuMuPAT.lowerPuritySelection = cms.string("isTrackerMuon") #O "isGlobalMuon"
process.HiOnia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.HiOnia2MuMuPAT.dimuonSelection=cms.string("2.0 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25 && pt>5.0")
process.HiOnia2MuMuPAT.addMCTruth = cms.bool(True)
process.HiOnia2MuMuPAT.addCommonVertex = cms.bool(True)
process.HiOnia2MuMuPAT.addMuonlessPrimaryVertex = cms.bool(False)
process.HiOnia2MuMuPAT.resolvePileUpAmbiguity = cms.bool(True)

# PHOTONS

#mc matching done by hand in rootupler - because it doesn't work well with OniaPhotonConversionProducer



# Chi parts
process.load('HeavyIonsAnalysis.ChiAnalysis.ChiAnalyzer_cfi')
process.ChiRootuple.muon_cand=cms.InputTag('ChiSelectedMuons')
process.ChiRootuple.isMC=cms.bool(True)

process.analysisPath = cms.Path(
            process.fastFilter *
            process.patMuonsWithTriggerSequence *
            process.ChiSelectedMuons * 
            process.HiOnia2MuMuPAT *
          #  process.PhotonCandidates *
            process.ChiSequence
)

