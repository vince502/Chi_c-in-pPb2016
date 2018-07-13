#This example can be run over files from AOD, therefore we need to build some information in fly.
#
outFileName = 'Chi_c_pPb8TeV_testMC2.root'
inFileNames = 'file:/afs/cern.ch/user/o/okukral/Work/ChicData/ChiCJpsiMuMu_Pythia8_8p16TeV_TuneCUETP8M1_RECO_90.root'

import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("ChiAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_pA_v4', '')

#process.GlobalTag.toGet = cms.VPSet(
#  cms.PSet(
#    record = cms.string("HeavyIonRcd"),
#    tag = cms.string("CentralityTable_HFtowersPlusTrunc200_EPOS8TeV_v80x01_mc"),
#    connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#    label = cms.untracked.string("HFtowersPlusTruncEpos")
#    )
#  )

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inFileNames))
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


# MUONS

#mc
import PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi
process.muonGenMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("muons"), # RECO objects to match  
    matched = cms.InputTag("genParticles"),   # mc-truth particle collection
    mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(1),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.5),  # Max deltaR for the match
    maxDPtRel = cms.double(0.5),  # Max deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True), # False = just match input in order; True = pick lowest deltaR pair first
)

# muons
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.ChiPATMuons = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
    muonSource = 'muons',
    #useParticleFlow = cms.bool (True),
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
    addGenMatch = cms.bool(True),   #Following lines for MC
    embedGenMatch    = cms.bool(True),                  #MC
    genParticleMatch = cms.InputTag("muonGenMatch")    #MC
)

# cuts on muons
process.ChiSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('ChiPATMuons'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 1.5))'
   ),
   filter = cms.bool(False)
)


#create dimuons 

process.load("HeavyIonsAnalysis.HiOnia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('ChiSelectedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
process.onia2MuMuPAT.higherPuritySelection = cms.string("isGlobalMuon") #O "isGlobalMuon"
process.onia2MuMuPAT.lowerPuritySelection = cms.string("isTrackerMuon") #O "isGlobalMuon"
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.dimuonSelection=cms.string("0.2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)


# make photon candidate conversions for P-wave studies.
# The low energy photons are reconstructed here.
import HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi
process.PhotonCandidates = HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi.PhotonCandidates.clone(
    conversions = 'allConversions',
    convAlgo    = 'undefined',
    convQuality = [''], #O: Changed ['highPurity','generalTracksOnly']
    primaryVertexTag = 'offlinePrimaryVertices',
    convSelection = 'conversionVertex.position.rho>0.0', #O: Changed 1.5
    wantTkVtxCompatibility = False,
    sigmaTkVtxComp = 50, #O: Changed 5
    wantCompatibleInnerHits = True,
    pfcandidates = 'particleFlow',
    pi0OnlineSwitch = False,
    TkMinNumOfDOF = 0, #O: Changed 3
    wantHighpurity = False,
    #test
    vertexChi2ProbCut = 0.0000,
    trackchi2Cut = 1000,
    minDistanceOfApproachMinCut = -100.25,
    minDistanceOfApproachMaxCut = 100.00,
    #O: Original
    #vertexChi2ProbCut = 0.0005,
    #trackchi2Cut = 10,
    #minDistanceOfApproachMinCut = -0.25,
    #minDistanceOfApproachMaxCut = 1.00,
)


# Chi parts
process.load('HeavyIonsAnalysis.ChiAnalysis.ChiAnalyzer_cfi')
process.ChiRootuple.muon_cand=cms.InputTag('ChiPATMuons')
process.ChiRootuple.isMC=cms.bool(True)

process.analysisPath = cms.Path(
            process.muonGenMatch *
            process.ChiPATMuons *
            process.ChiSelectedMuons * 
            process.onia2MuMuPAT *
            process.PhotonCandidates *
            process.ChiSequence
)

