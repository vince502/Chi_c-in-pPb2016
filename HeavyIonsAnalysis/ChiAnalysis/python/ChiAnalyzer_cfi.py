import FWCore.ParameterSet.Config as cms


ChiCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("ChiProducer","ChiCandidates"),
    minNumber = cms.uint32(1),
    )

ChiRootuple = cms.EDAnalyzer('ChiRootupler',
    muon_cand =  cms.InputTag("ChiSelectedMuons"), 
    dimuon_cand = cms.InputTag("HiOnia2MuMuPAT"),
    photon_cand = cms.InputTag("PhotonCandidates","conversions"),
    conversions_ch = cms.InputTag("allConversions"),
    #chi_cand = cms.InputTag("ChiProd","ChiCandidates"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    centralityInfo = cms.InputTag("pACentrality"),
    genParticlesTag = cms.InputTag("genParticles"),
    isMC = cms.bool(False)
    )

ChiSequence = cms.Sequence(
    #* ChiCounter 
    ChiRootuple
    )
