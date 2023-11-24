import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import *

def onia2MuMuPAT(process, GlobalTag, MC=False, HLT='HLT', Filter=True, useL1Stage2=False):
    # Setup the process
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
        # fileMode = cms.untracked.string('MERGE'),
    )
     
    # Drop the DQM stuff on input
    if hasattr(process, "source") and hasattr(process.source, "inputCommands"): 
        process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
    else:
        process.source = cms.Source("PoolSource",
                                    inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*"),
                                    fileNames = cms.untracked.vstring())

    # Prune generated particles to muons and their parents
    process.genMuons = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring("drop  *  ",                      # this is the default
            "++keep abs(pdgId) = 13"          # keep muons and their parents
        ))

    # Make PAT Muons
    process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
    from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addHLTL1Passthrough, useL1Stage2Candidates

    # with some customization
    if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it
        # one-to-one
        process.muonMatch.maxDeltaR = cms.double(0.05)
        process.muonMatch.resolveByMatchQuality = True
        process.muonMatch.matched = "genMuons"
    changeTriggerProcessName(process, HLT)
    switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the
                                          # same trigger muon
    addHLTL1Passthrough(process)
    
    if useL1Stage2:
        useL1Stage2Candidates(process)
        process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGmtStage2Digis:Muon")')
        process.muonMatchHLTL1.useMB2InOverlap = cms.bool(True)
        process.muonMatchHLTL1.useStage2L1 = cms.bool(True)
        process.muonMatchHLTL1.preselection = cms.string("")

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
    process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
    process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
    process.muonMatchHLTTrackMu.maxDeltaR = 0.1
    process.muonMatchHLTTrackMu.maxDPtRel = 10.0
    
    process.patTrigger.collections.append("hltHIL3MuonCandidates") 
    process.muonMatchHLTL3T.maxDeltaR = 0.1
    process.muonMatchHLTL3T.maxDPtRel = 10.0
    process.muonMatchHLTL3T.matchedCuts = cms.string('coll("hltHIL3MuonCandidates")') 
    process.muonMatchHLTL3T.resolveAmbiguities = False
      
    # Make a sequence
    process.patMuonSequence = cms.Sequence(process.hltOniaHI * process.genMuons * process.patMuonsWithTriggerSequence)
    if not MC:
        process.patMuonSequence.remove(process.genMuons)
      
    # Make dimuon candidates
    process.onia2MuMuPatGlbGlb = cms.EDProducer('HiOnia2MuMuPAT',
        muons                    = cms.InputTag("patMuonsWithTrigger"),
        beamSpotTag              = cms.InputTag("offlineBeamSpot"),
        primaryVertexTag         = cms.InputTag("offlinePrimaryVertices"),
        srcTracks                = cms.InputTag("generalTracks"),
        genParticles             = cms.InputTag("genParticles"),
        # At least one muon must pass this selection
        higherPuritySelection    = cms.string("((isGlobalMuon && isTrackerMuon) || (innerTrack.isNonnull && genParticleRef(0).isNonnull)) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35"),
        # BOTH muons must pass this selection
        lowerPuritySelection     = cms.string("((isGlobalMuon && isTrackerMuon) || (innerTrack.isNonnull && genParticleRef(0).isNonnull)) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35"),
        dimuonSelection          = cms.string(""), ## The dimuon must pass this selection before vertexing
        addCommonVertex          = cms.bool(True), ## Embed the full reco::Vertex out of the common vertex fit
        addMuonlessPrimaryVertex = cms.bool(False), ## Embed the primary vertex re-made from all the tracks except the two muons
        addMCTruth               = cms.bool(MC),      ## Add the common MC mother of the two muons, if any
        resolvePileUpAmbiguity   = cms.bool(True)   ## Order PVs by their vicinity to the J/psi vertex, not by sumPt
    )

    # check if there is at least one (inclusive) global+tracker di-muon
    process.onia2MuMuPatGlbGlbFilter = cms.EDFilter("CandViewCountFilter",
        src       = cms.InputTag('onia2MuMuPatGlbGlb'),
        minNumber = cms.uint32(1),)

    # the onia2MuMu path
    process.Onia2MuMuPAT = cms.Path(process.patMuonSequence * process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilter)

    process.outOnia2MuMu = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('onia2MuMuPAT.root'),
        outputCommands =  cms.untracked.vstring('drop *',                       
            'keep *_mergedtruth_*_*',                              # tracking particles and tracking vertices for
                                                                   # hit by hit matching
            'keep *_genParticles_*_*',                             # generated particles
            'keep *_genMuons_*_Onia2MuMuPAT',                      # generated muons and parents
            'keep patMuons_patMuonsWithTrigger_*_Onia2MuMuPAT',    # All PAT muons including matches to triggers
            'keep patCompositeCandidates_*__Onia2MuMuPAT',         # PAT di-muons
            'keep *_offlinePrimaryVertices_*_*',                   # Primary vertices: you want these to compute impact
                                                                   # parameters
            'keep *_offlineBeamSpot_*_*',                          # Beam spot: you want this for the same reason
            'keep edmTriggerResults_TriggerResults_*_*',           # HLT info, per path (cheap)
            'keep *_hltGmtStage2Digis_*_*',                        # Stage2 L1 Muon info
            'keep *_gmtStage2Digis_*_*',                           # Stage2 L1 Muon info
            'keep *_hltGtStage2Digis_*_*',                         # Stage2 L1 Muon info
            'keep *_gtStage2Digis_*_*',                            # Stage2 L1 Muon info
            'keep *_hltGtStage2ObjectMap_*_*',                     # Stage2 L1 Muon info
            'keep L1GlobalTriggerReadoutRecord_*_*_*',             # For HLT and L1 prescales (cheap)
            'keep L1GlobalTriggerRecord_*_*_*',                    # For HLT and L1 prescales (cheap)
            'keep L1GtTriggerMenu_*_*_*',                          # L1 prescales
            'keep *_pACentrality_*_*',                             # PA Centrality
            ),
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('Onia2MuMuPAT')) if Filter else cms.untracked.PSet())

    process.e = cms.EndPath(process.outOnia2MuMu)