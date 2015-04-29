#This example can be run over files from AODSIM (that is MC), therefore we need to build some information in fly.
#
outFileName = 'chic-aodsim-output.root'
inFileNames = 'file:aodsim-input.root'

import FWCore.ParameterSet.Config as cms

process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inFileNames))
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

# muons without trigger info, alternatively in recent version of 80x you can use muon with trigger as well.
# we will select muons and create Onia2MuMu pairs.
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.oniaPATMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
    muonSource = 'muons',
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    embedPFCandidate    = False,
    embedCaloMETMuonCorrs = cms.bool(False),
    embedTcMETMuonCorrs   = cms.bool(False),
    embedPfEcalEnergy     = cms.bool(False),
    embedPickyMuon = False,
    embedTpfmsMuon = False,
    userIsolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(),     # no heavy isodeposits
    addGenMatch = False,          # no mc
)

process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('oniaPATMuonsWithoutTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 1.5))'
   ),
   filter = cms.bool(True)
)

process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('oniaSelectedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.dimuonSelection=cms.string("0.2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

# make photon candidate conversions for P-wave studies.
# The low energy photons are reconstructed here.
import HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi
process.oniaPhotonCandidates = HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi.PhotonCandidates.clone()

process.load('Ponia.OniaPhoton.OniaPhotonProducer_cfi')
process.load('Ponia.OniaPhoton.OniaPhotonRootupler_cfi')

tag_dimuon="onia2MuMuPAT"
tag_chi_conv_prod="oniaPhotonCandidates"
tag_chi_conv_lab="conversions"
process.DiMuonCounter.src = cms.InputTag(tag_dimuon)
process.PhotonCounter.src  = cms.InputTag(tag_chi_conv_prod,tag_chi_conv_lab)
process.ChiProducer.conversions = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab)

# reduce MC genParticles a la miniAOD
process.load('PhysicsTools.PatAlgos.slimming.genParticles_cff')
process.packedGenParticles.inputVertices = cms.InputTag('offlinePrimaryVertices')

process.p = cms.Path(
            process.prunedGenParticlesWithStatusOne*process.prunedGenParticles*process.packedGenParticles *
            process.oniaPATMuonsWithoutTrigger *
            process.oniaSelectedMuons *
            process.onia2MuMuPAT *
            process.oniaPhotonCandidates *
            process.ChiSequence*process.rootuple
)
process.rootuple.isMC = cms.bool(True)
