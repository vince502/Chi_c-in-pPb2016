import FWCore.ParameterSet.Config as cms

############################# CONFIGURATION PARAMETERS #############################
tag_dimuon = 'onia2MuMuPAT'  # Tag name of the dimuon collection as saved from process.dimuonProducer
cut_dimuon_Mass_low = 2.9
cut_dimuon_Mass_high = 3.3
cut_dimuon_Pt_min = 10.0 # O: Not used here
cut_dimuon_rapidity = 2.1 # O: Not used here
cut_dimuon_vprob = 1.00     # Minimum vertex probability for dimuon candidate # O: Not used here
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

# Speedup, analyzing only events with at least one Onia and Conversion Candidates
DiMuonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag(tag_dimuon),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
    )

PhotonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag(tag_chi_conv_prod,tag_chi_conv_lab),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
    )

ChiProducer = cms.EDProducer('OniaPhotonProducer',
    conversions     = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab),
    dimuons         = cms.InputTag(tag_dimuon),
    pi0OnlineSwitch = cms.bool(pi0_online_switch),
    deltaMass       = cms.vdouble(chi_deltaM_min, chi_deltaM_max),
    dzmax           = cms.double(chi_dzMax),
    triggerMatch    = cms.bool(triggermatch_switch),
    )

ChiCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("ChiProducer","ChiCandidates"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
    )

ChiKinFitter = cms.EDProducer('OniaPhotonKinematicFit',
    chi_cand     = cms.InputTag("ChiProducer","ChiCandidates"),
    upsilon_mass = cms.double(jpsi_mass),
    product_name = cms.string("ChiCandidatesFit")
    )

ChiSequence = cms.Sequence(
                            DiMuonCounter
                            * PhotonCounter
                            * ChiProducer
                            * ChiCounter
                            * ChiKinFitter
    )
