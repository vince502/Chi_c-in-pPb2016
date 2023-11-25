import FWCore.ParameterSet.Config as cms


ChiCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("ChiProducer","ChiCandidates"),
    minNumber = cms.uint32(1),
    )

ChiRootuple = cms.EDAnalyzer('ChiRootupler',
    muon_cand =  cms.InputTag("ChiSelectedMuons"), 
    dimuon_cand = cms.InputTag("HiOnia2MuMuPAT"),
    conversions_ch = cms.InputTag("allConversions"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    centralityInfo = cms.InputTag("pACentrality"),
    srcTracks = cms.InputTag("generalTracks"), # this is added to determine ntrack variable with the cuts, the one we are to use for centrality dependencies
    genParticlesTag = cms.InputTag("genParticles"),
    isMC = cms.bool(False)
    )

ChiSequence = cms.Sequence(
    #* ChiCounter 
    ChiRootuple
    )
