import FWCore.ParameterSet.Config as cms

############################# CONFIGURATION PARAMETERS #############################
tag_dimuon = 'onia2MuMuPAT'  # Tag name of the dimuon collection as saved from process.dimuonProducer
cut_dimuon_Mass_low = 2.9
cut_dimuon_Mass_high = 3.3
cut_dimuon_Pt_min = 10.0
cut_dimuon_rapidity = 2.1
cut_dimuon_vprob = 0.01     # Minimum vertex probability for dimuon candidate
#
tag_chi_conv_prod = 'PhotonCandidates'
tag_chi_conv_lab = 'conversions'
pi0_online_switch = False
chi_deltaM_min = 0.0        # This two values define the minimum and the maximum values 
chi_deltaM_max = 2.0        # required for the QValue of the Chi candidate
chi_dzMax = 0.5
jpsi_mass = 3.0969          # in GeV
triggermatch_switch = False
############################# CONFIGURATION END ####################################

# Filter to analyzing only events with at least one dimuon
DimuonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag(tag_dimuon),
    minNumber = cms.uint32(1),
    filter    = cms.bool(False)
    )

# Filter to analyzing only events with at least one photon
PhotonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag(tag_chi_conv_prod,tag_chi_conv_lab),
    minNumber = cms.uint32(1),
    filter    = cms.bool(False)
    )

#Producer creates the chi_c candidates 
ChiProd = cms.EDProducer('ChiProducer',
    #conversions     = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab),
    #dimuons         = cms.InputTag(tag_dimuon),
    #pi0OnlineSwitch = cms.bool(pi0_online_switch),
    #deltaMass       = cms.vdouble(chi_deltaM_min, chi_deltaM_max),
    #dzmax           = cms.double(chi_dzMax),
    #triggerMatch    = cms.bool(triggermatch_switch),
    )

ChiCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("ChiProducer","ChiCandidates"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(False)
    )

ChiRootuple = cms.EDAnalyzer('ChiRootupler',
    muon_cand =  cms.InputTag("ChiPATMuons"), 
    dimuon_cand = cms.InputTag("onia2MuMuPAT"),
    photon_cand = cms.InputTag("PhotonCandidates","conversions"),
    conversions_ch = cms.InputTag("allConversions"),
    chi_cand = cms.InputTag("ChiProducer","ChiCandidates"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    genParticlesTag = cms.InputTag("genParticles"),
    isMC = cms.bool(False)
    )

ChiSequence = cms.Sequence(
    #DimuonCounter
    #* PhotonCounter
    ChiProd
    #* ChiCounter 
    * ChiRootuple
    )
