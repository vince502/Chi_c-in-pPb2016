# To be run on AOD, done by Ota Kukral
#
outFileName = 'Chi_c_pPb8TeV_testNew8.root'
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
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_v19', '')

#process.GlobalTag.toGet = cms.VPSet(
#  cms.PSet(
#    record = cms.string("HeavyIonRcd"),
#    tag = cms.string("CentralityTable_HFtowersPlusTrunc200_EPOS8TeV_v80x01_mc"),
#    connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#    label = cms.untracked.string("HFtowersPlusTruncEpos")
#    )
#  )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inFileNames))
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


## we will select muons 
#import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
#process.ChiPATMuons = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
#    muonSource = 'muons',
#    #useParticleFlow = cms.bool (True),
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

## cuts on muons
#process.ChiSelectedMuons = cms.EDFilter('PATMuonSelector',
#   src = cms.InputTag('ChiPATMuons'),
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
process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_PAL1DoubleMu*_v*')
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")
process.triggerResultsFilter.l1tResults = cms.InputTag("gtStage2Digis") #O: needs to be this
process.triggerResultsFilter.throw = False
### Filter sequence
process.fastFilter = cms.Sequence(process.hfCoincFilter + process.primaryVertexFilterPA + process.noScraping + process.triggerResultsFilter)


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
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 0.7))'
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
process.HiOnia2MuMuPAT.dimuonSelection=cms.string("2.0 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
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
process.ChiRootuple.muon_cand=cms.InputTag('ChiSelectedMuons')




process.analysisPath = cms.Path(
            process.fastFilter *
            process.patMuonsWithTriggerSequence *
            process.ChiSelectedMuons * 
            process.HiOnia2MuMuPAT *
            process.PhotonCandidates *
            process.ChiSequence
)


#process.out = cms.OutputModule("PoolOutputModule",
#        fileName = cms.untracked.string('TestOut.root'),
#        outputCommands =  cms.untracked.vstring(
#            'drop *',
#            'keep patMuons_patMuonsWithTrigger_*_*',    # All PAT muons including matches to triggers       
#            'keep *_centralityBin_*_*',                            # PA Centrality
#            'keep *_hiCentrality_*_*',                             # PA Centrality
#            'keep *_pACentrality_*_*',                             # PA Centrality
#            'keep *_offlinePrimaryVertices_*_*',
#            'keep *_TriggerResults_*_*',
#            ),
#    )

#process.o= cms.EndPath(process.out)


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
