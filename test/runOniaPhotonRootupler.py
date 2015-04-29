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

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:inputfile.root'))
process.TFileService = cms.Service("TFileService",fileName = cms.string('runOniaPhotonRootupler.root'))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.load('Ponia.OniaPhoton.OniaPhotonProducer_cfi')
process.load('Ponia.OniaPhoton.OniaPhotonRootupler_cfi')
process.p = cms.Path(process.ChiSequence*process.rootuple)
process.rootuple.isMC = cms.bool(False)
