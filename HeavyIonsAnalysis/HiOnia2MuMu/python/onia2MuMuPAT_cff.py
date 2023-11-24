import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import *

def onia2MuMuPAT(process, GlobalTag, MC=False, HLT='HLT', Filter=True, useL1Stage2=False, doTrimuons=False, DimuonTrk=False, flipJpsiDir=0):
    # Setup the process
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
    )

    # Make PAT Muons
    process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
    from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addHLTL1Passthrough, useL1Stage2Candidates

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
    changeTriggerProcessName(process, HLT)
    switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
    addHLTL1Passthrough(process)

    if useL1Stage2:
        useL1Stage2Candidates(process)
        process.patTrigger.collections.append("hltGtStage2Digis:Muon")
        process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGtStage2Digis:Muon")')
        process.muonMatchHLTL1.useMB2InOverlap = cms.bool(True)
        process.muonMatchHLTL1.useStage2L1 = cms.bool(True)
        process.muonMatchHLTL1.preselection = cms.string("")

    #This is NEEDED in PbPb, but should be REMOVED in pp
    process.patTrigger.collections.append("hltIterL3MuonCandidatesPPOnAA")
    process.patTrigger.collections.append("hltL2MuonCandidatesPPOnAA")
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

    # Make a sequence
    process.patMuonSequence = cms.Sequence(
            process.patMuonsWithTriggerSequence
    )
    if MC:
        process.patMuonSequence.insert(0, process.genMuons)

    # Make dimuon candidates
    from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cfi import onia2MuMuPAT
    process.onia2MuMuPatGlbGlb = onia2MuMuPAT.clone(
            doTrimuons        = cms.bool(doTrimuons), ## Make collections of trimuon candidates in addition to dimuons, and keep only events with >0 trimuons
            DimuonTrk         = cms.bool(DimuonTrk), ## Make collections of Jpsi+track candidates in addition to dimuons, and keep only events with >0 Jpsi+trk
            flipJpsiDirection = cms.int32(flipJpsiDir), ## flip the Jpsi direction, before combining it with a third muon
    )

    # check if there is at least one (inclusive) di-muon. BEWARE this can cause trouble in .root output if no event is selected by onia2MuMuPatGlbGlbFilter!
    process.onia2MuMuPatGlbGlbFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('onia2MuMuPatGlbGlb',''),
        minNumber = cms.uint32(1),
    )
    process.onia2MuMuPatGlbGlbFilterDimutrk = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('onia2MuMuPatGlbGlb','dimutrk'),
        minNumber = cms.uint32(1),
    )
    process.onia2MuMuPatGlbGlbFilterTrimu = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('onia2MuMuPatGlbGlb','trimuon'),
        minNumber = cms.uint32(1),
    )
    process.filter3mu = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('muons'),
        minNumber = cms.uint32(3),
    )
    process.pseudoDimuon = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string('muons@+ muons@-'),
        cut = cms.string('2.4 < mass < 3.7'),
    )
    process.pseudoDimuonFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('pseudoDimuon'),
        minNumber = cms.uint32(1),
    )
    process.pseudoDimuonFilterSequence = cms.Sequence(
        process.pseudoDimuon *
        process.pseudoDimuonFilter
    )

    # the onia2MuMu path
    process.Onia2MuMuPAT = cms.Path(
        process.patMuonSequence *
        process.onia2MuMuPatGlbGlb *
        process.onia2MuMuPatGlbGlbFilter
    )

    if DimuonTrk:
        process.Onia2MuMuPAT.replace(process.onia2MuMuPatGlbGlbFilter, process.onia2MuMuPatGlbGlbFilterDimutrk)
    if doTrimuons:
        process.Onia2MuMuPAT.replace(process.onia2MuMuPatGlbGlbFilter, process.onia2MuMuPatGlbGlbFilterTrimu)

    process.outOnia2MuMu = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('onia2MuMuPAT.root'),
        outputCommands =  cms.untracked.vstring(
            'drop *',
            'keep *_mergedtruth_*_*',                              # tracking particles and tracking vertices for hit by hit matching
            'keep *_genParticles_*_*',                             # generated particles
            'keep *_genMuons_*_Onia2MuMuPAT',                      # generated muons and parents
            'keep patMuons_patMuonsWithTrigger_*_Onia2MuMuPAT',    # All PAT muons including matches to triggers
            'keep patCompositeCandidates_*__Onia2MuMuPAT',         # PAT di-muons
            'keep patCompositeCandidates_*_trimuon_Onia2MuMuPAT',  # PAT trimuons
            'keep patCompositeCandidates_*_dimutrk_Onia2MuMuPAT',  # PAT dimuon+track candidates
            'keep *_dedxHarmonic2_*_*',                            # dE/dx estimator for tracks (to do PID when doDimuTrk)
            'keep *_offlineSlimmedPrimaryVertices_*_*',                   # Primary vertices: you want these to compute impact parameters
            'keep *_offlinePrimaryVertices_*_*',                   # Primary vertices: you want these to compute impact parameters
            'keep *_inclusiveSecondaryVertices_*_*',      # Secondary vertices: to check if non-prompt muons come from a common SV
            'keep *_inclusiveSecondaryVerticesLoose_*_*',      # Secondary vertices: to check if non-prompt muons come from a common SV
            'keep *_inclusiveCandidateSecondaryVertices_*_*',      # Secondary vertices: to check if non-prompt muons come from a common SV
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
            'keep *_centralityBin_*_*',                            # PA Centrality
            'keep *_hiCentrality_*_*',                             # PA Centrality
            'keep *_pACentrality_*_*',                             # PA Centrality
            ),
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('Onia2MuMuPAT') ) if Filter else cms.untracked.PSet()
    )


def changeToMiniAOD(process):

    if hasattr(process, "patMuonsWithTrigger"):
        from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import useExistingPATMuons
        useExistingPATMuons(process, newPatMuonTag=cms.InputTag("unpackedMuons"), addL1Info=False)

        process.patTriggerFull = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
            patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
            triggerResults              = cms.InputTag('TriggerResults::HLT'),
            unpackFilterLabels          = cms.bool(True)
        )
        process.load('HiAnalysis.HiOnia.unpackedTracksAndVertices_cfi')
        process.patMuonSequence.insert(0, process.unpackedTracksAndVertices)
        process.load('HiAnalysis.HiOnia.unpackedMuons_cfi')
        process.patMuonSequence.insert(1, process.unpackedMuons)

        process.outOnia2MuMu.outputCommands.append('keep *Vert*_unpackedTracksAndVertices_*_*')
        process.outOnia2MuMu.outputCommands.append('keep patMuons_unpackedMuons_*_*')
        process.outOnia2MuMu.outputCommands.append('drop patMuons_patMuonsWith*_*_*')

        if hasattr(process, "muonMatch"):
            from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import changeRecoMuonInput
            changeRecoMuonInput(process, recoMuonCollectionTag=cms.InputTag("unpackedMuonsWithGenMatch"), oldRecoMuonCollectionTag=cms.InputTag("unpackedMuons"))
            process.load('HiAnalysis.HiOnia.unpackedMuonsWithGenMatch_cfi')
            process.patMuonsWithTriggerSequence.insert(1, process.unpackedMuonsWithGenMatch)

            process.onia2MuMuPatGlbGlb.genParticles = "prunedGenParticles"
            process.genMuons.src = "prunedGenParticles"
            process.muonMatch.src = "unpackedMuons"

            process.outOnia2MuMu.outputCommands.append('keep *_prunedGenParticles_*_*')
            process.outOnia2MuMu.outputCommands.append('keep patMuons_unpackedMuonsWithGenMatch_*_*')
            process.outOnia2MuMu.outputCommands.append('drop patMuons_unpackedMuons_*_*')


    from HLTrigger.Configuration.CustomConfigs import massReplaceInputTag
    process = massReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
    process = massReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")
    process = massReplaceInputTag(process,"genParticles","prunedGenParticles")
