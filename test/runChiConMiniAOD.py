# This example can be run over files from modified MiniAOD, You should no be using this unless you know what you are doing.
# The modified MiniAOD has the conversion information already added. We are working to get that into official cmssw. 
#
outFileName = 'chic-miniaod-output.root'
inFileNames = 'file:miniaod-input.root'

import FWCore.ParameterSet.Config as cms

process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inFileNames))
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

# In MiniAOD, the PATMuons are already present. We just need to run Onia2MuMu, with a selection of muons.
process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuons'),
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
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.dimuonSelection=cms.string("0.2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

# The low energy photon, are already stored in the same format as BPHSkim

process.load('Ponia.OniaPhoton.OniaPhotonProducer_cfi')
process.load('Ponia.OniaPhoton.OniaPhotonRootupler_cfi')

tag_dimuon="onia2MuMuPAT"
tag_chi_conv_prod="oniaPhotonCandidates"
tag_chi_conv_lab="conversions"
process.DiMuonCounter.src = cms.InputTag(tag_dimuon)
process.PhotonCounter.src  = cms.InputTag(tag_chi_conv_prod,tag_chi_conv_lab)
process.ChiProducer.conversions = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab)

process.p = cms.Path(
            process.oniaSelectedMuons *
            process.onia2MuMuPAT *
            process.ChiSequence*process.rootuple
)
process.rootuple.isMC = cms.bool(False)
process.rootuple.primaryVertices = cms.InputTag('offlineSlimmedPrimaryVertices')
