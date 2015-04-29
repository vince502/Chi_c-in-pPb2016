input_filename = 'BPHSkim_*.root'
ouput_filename = 'rootuple-chib.root'

import FWCore.ParameterSet.Config as cms
import os,fnmatch

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

dir = '/meson/data/store/group/phys_bphys/asanchez/MuOnia/BPHSkim-v1-Run2015D-16Dec2015-v1/160519_182152/0000/'
inputfiles = []
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, input_filename):
        filename = 'file:'+dir+str(file)
        inputfiles.append(filename)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inputfiles))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.load('Ponia.OniaPhoton.OniaPhotonProducer_cfi')
process.ChiKinFitter.upsilon_mass = cms.double(9.4603)  # upsilon(1s) for now
process.load('Ponia.OniaPhoton.OniaPhotonRootupler_cfi')
process.p = cms.Path(process.ChiSequence*process.rootuple)
process.rootuple.isMC = cms.bool(False)
