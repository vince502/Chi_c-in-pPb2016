#This example needs to be run over files from the BPH official SKIM, which contains all relevant branches.
#
outFileName = 'chic-bph-output.root'
inFileNames = 'file:bphskim-input.root'

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
process.load('Ponia.OniaPhoton.OniaPhotonProducer_cfi')
process.load('Ponia.OniaPhoton.OniaPhotonRootupler_cfi')

tag_dimuon="onia2MuMuPAT"
tag_chi_conv_prod="oniaPhotonCandidates"
tag_chi_conv_lab="conversions"
process.DiMuonCounter.src = cms.InputTag(tag_dimuon)
process.PhotonCounter.src  = cms.InputTag(tag_chi_conv_prod,tag_chi_conv_lab)
process.ChiProducer.conversions = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab)

process.p = cms.Path(process.ChiSequence*process.rootuple)
process.rootuple.isMC = cms.bool(False)
