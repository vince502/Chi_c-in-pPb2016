# To be run on AOD, done by Ota Kukral
#

MC = False

outFileName = 'Chi_c_PbPb_5360GeV.root'
inFileNames = 'file:/eos/cms/store/group/phys_heavyions/dileptons/Data2023/MINIAOD/HIPhysicsRawPrime0/Run375064/7ed5766f-6b1d-415e-8916-e62825a6347f.root'
#inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/ChiCJpsiMuMu_Pythia8_8p16TeV_TuneCUETP8M1_RECO_3.root' #v6 MC

import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("ChiAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
# load the Geometry and Magnetic Field
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v4', '')

### For Centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
print('\n\033[31m~*~ USING PRELIMINARY CENTRALITY TABLE FOR 2023 PbPb DATA (Run 374810)~*~\033[0m\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374810"),
        connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374810.db"),
        label = cms.untracked.string("HFtowers")
        ),
    ])



process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inFileNames))
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Offline event filters
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')

# HLT trigger firing events
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHI.HLTPaths = ["HLT_HIL*SingleMu*_v*", "HLT_HIL*DoubleMu*_v*", "HLT_HIMinimumBiasHF1AND*_v*"]
process.hltHI.throw = False
process.hltHI.andOr = True
### Filter sequence
process.fastFilter = cms.Sequence(process.phfCoincFilter2Th4 * process.primaryVertexFilter * process.hltHI * process.clusterCompatibilityFilter)

# Make PAT Muons
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *

# with some customization
if MC:
     # Prune generated particles to muons and their parents
    process.genMuons = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "drop  *  ",                      # this is the default
            "++keep abs(pdgId) = 13"          # keep muons and their parents
        )
    )
    addMCinfo(process)
    # since we match inner tracks, keep the matching tight and make it one-to-one
    process.muonMatch.maxDeltaR = cms.double(0.05)
    process.muonMatch.resolveByMatchQuality = True
    process.muonMatch.matched = "genMuons"
    process.muonMatch.src = "muons"
changeTriggerProcessName(process, 'HLT')
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
addHLTL1Passthrough(process)


#This is NEEDED in PbPb, but should be REMOVED in pp
process.patTrigger.collections.append("hltIterL3MuonCandidatesPPOnAA")
process.patTrigger.collections.append("hltL2MuonCandidatesPPOnAA")
process.patTrigger.collections.append("hltGtStage2Digis:Muon")
process.muonMatchHLTL3.matchedCuts = cms.string('coll("hltIterL3MuonCandidatesPPOnAA")')
process.muonMatchHLTL2.matchedCuts = cms.string('coll("hltL2MuonCandidatesPPOnAA")')

from HLTrigger.Configuration.HLT_FULL_cff import fragment
process.hltESPSteppingHelixPropagatorAlong = fragment.hltESPSteppingHelixPropagatorAlong
process.hltESPSteppingHelixPropagatorOpposite = fragment.hltESPSteppingHelixPropagatorOpposite

process.muonL1Info.maxDeltaR = 0.3
process.muonL1Info.maxDeltaEta   = 0.2
process.muonL1Info.fallbackToME1 = True
process.muonL1Info.useStation2 = cms.bool(True)
process.muonL1Info.cosmicPropagationHypothesis = cms.bool(False)
process.muonL1Info.propagatorAlong = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorAlong')
process.muonL1Info.propagatorAny = cms.ESInputTag('', 'SteppingHelixPropagatorAny')
process.muonL1Info.propagatorOpposite = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorOpposite')
process.muonMatchHLTL1.maxDeltaR = 0.3
process.muonMatchHLTL1.maxDeltaEta   = 0.2
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL1.useStation2 = cms.bool(True)
process.muonMatchHLTL1.cosmicPropagationHypothesis = cms.bool(False)
process.muonMatchHLTL1.propagatorAlong = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorAlong')
process.muonMatchHLTL1.propagatorAny = cms.ESInputTag('', 'SteppingHelixPropagatorAny')
process.muonMatchHLTL1.propagatorOpposite = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorOpposite')
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0
## For trigger muons
#switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon

## For L1 muons
useL1Stage2Candidates(process)
process.patTrigger.collections.remove("hltL1extraParticles")
process.patTrigger.collections.append("hltGmtStage2Digis:Muon")
process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGmtStage2Digis:Muon")')
process.muonMatchHLTL1.useStage2L1 = cms.bool(True)
process.muonMatchHLTL1.useMB2InOverlap = cms.bool(True)
process.muonMatchHLTL1.preselection = cms.string("")
appendL1MatchingAlgo(process)
# Make a sequence
process.patMuonSequence = cms.Sequence(
        process.patMuonsWithTriggerSequence
)
if MC:
    process.patMuonSequence.insert(0, process.genMuons)


# cuts on muons
process.ChiSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('unpackedMuons'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 1.0))'
   ),
   filter = cms.bool(False)
)

#DIMUONS 

process.load("HeavyIonsAnalysis.HiOnia2MuMu.HiOnia2MuMuPAT_cfi")
process.HiOnia2MuMuPAT.muons=cms.InputTag('ChiSelectedMuons')
process.HiOnia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
#process.HiOnia2MuMuPAT.higherPuritySelection = cms.string("isGlobalMuon") #O "isGlobalMuon"
#process.HiOnia2MuMuPAT.lowerPuritySelection = cms.string("isTrackerMuon") #O "isGlobalMuon"
process.HiOnia2MuMuPAT.higherPuritySelection = cms.string("isTrackerMuon") 
process.HiOnia2MuMuPAT.lowerPuritySelection = cms.string("isTrackerMuon") 
process.HiOnia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.HiOnia2MuMuPAT.dimuonSelection=cms.string("2.0 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25 && pt>5.0")
process.HiOnia2MuMuPAT.addMCTruth = cms.bool(False)
process.HiOnia2MuMuPAT.addCommonVertex = cms.bool(True)
process.HiOnia2MuMuPAT.addMuonlessPrimaryVertex = cms.bool(False)
process.HiOnia2MuMuPAT.resolvePileUpAmbiguity = cms.bool(True)



## PHOTONS

##mc matching done by hand in rootupler - because it doesn't work well with OniaPhotonConversionProducer

## The low energy photons are reconstructed here.
import HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi
process.PhotonCandidates = HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi.PhotonCandidates.clone(
    conversions = 'reducedEgamma',
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
process.ChiRootuple.conversions_ch = cms.InputTag("reducedEgamma","reducedConversions")
process.ChiRootuple.centralityInfo =cms.InputTag("hiCentrality")




process.analysisPath = cms.Path(
            process.fastFilter *
            process.patMuonSequence *
            process.ChiSelectedMuons * 
            process.HiOnia2MuMuPAT *
          #  process.PhotonCandidates *
            process.ChiSequence
)
from HeavyIonsAnalysis.HiOnia2MuMu.onia2MuMuPAT_cff import changeToMiniAOD
changeToMiniAOD(process)
process.unpackedMuons.addPropToMuonSt = cms.bool(True)

process.options.numberOfThreads = 4


process.schedule  = cms.Schedule( process.analysisPath )

