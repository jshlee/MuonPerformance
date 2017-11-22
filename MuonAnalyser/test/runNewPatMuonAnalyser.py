import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("NewPatMuonAnalyser",eras.Phase2)
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

# Beware, in this area the wild character is not working!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_0_0_pre1/RelValTTbar_13/MINIAODSIM/94X_mc2017_realistic_v10-v1/10000/7CBB8BF1-9FCB-E711-8FBD-0CC47A78A3F8.root'
    ),
)

process.load('RecoMuon.MuonIsolation.muonIsolationPUPPI_cff')

process.puppiNewIso = process.muonIsolationMiniAODPUPPINoLeptons.clone()
process.trkNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))
#process.trkNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False), useOnlyTrack = cms.bool(True))
process.pfNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))

process.minipuppiNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(True))
process.minitrkNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))
process.minipfNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))

process.TFileService = cms.Service("TFileService",fileName = cms.string("out.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200))

process.NewPatMuonAnalyser = cms.EDAnalyzer("NewPatMuonAnalyser",
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    addPileupInfo = cms.InputTag("slimmedAddPileupInfo"),
    muons = cms.InputTag("slimmedMuons"),
    pruned = cms.InputTag("prunedGenParticles"),
    pfCands = cms.InputTag("packedPFCandidates"),
    miniIsoParams = cms.vdouble(0.05, 0.2, 10.0, 0.5, 0.0001, 0.01, 0.01, 0.01, 0.0),
    
    puppiNewIso_ch = cms.InputTag("puppiNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    puppiNewIso_nh = cms.InputTag("puppiNewIso","h0-DR030-ThresholdVeto000-ConeVeto001"),
    puppiNewIso_ph = cms.InputTag("puppiNewIso","gamma-DR030-ThresholdVeto000-ConeVeto001"),    
    trkNewIso = cms.InputTag("trkNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    pfNewIso_ch = cms.InputTag("pfNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    pfNewIso_nh = cms.InputTag("pfNewIso","h0-DR030-ThresholdVeto000-ConeVeto001"),
    pfNewIso_ph = cms.InputTag("pfNewIso","gamma-DR030-ThresholdVeto000-ConeVeto001"),
    #pfNewIso_pu = cms.InputTag("pfNewIso","pu-DR030-ThresholdVeto000-ConeVeto000"),
    minipuppiNewIso_ch = cms.InputTag("minipuppiNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    minipuppiNewIso_nh = cms.InputTag("minipuppiNewIso","h0-DR030-ThresholdVeto000-ConeVeto001"),
    minipuppiNewIso_ph = cms.InputTag("minipuppiNewIso","gamma-DR030-ThresholdVeto000-ConeVeto001"),
    minitrkNewIso = cms.InputTag("minitrkNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    minipfNewIso_ch = cms.InputTag("minipfNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    minipfNewIso_nh = cms.InputTag("minipfNewIso","h0-DR030-ThresholdVeto000-ConeVeto001"),
    minipfNewIso_ph = cms.InputTag("minipfNewIso","gamma-DR030-ThresholdVeto000-ConeVeto001"),
    #minipfNewIso_pu = cms.InputTag("minipfNewIso","pu-DR030-ThresholdVeto000-ConeVeto000"),

  )

process.p = cms.Path(
    process.puppiNewIso+process.trkNewIso+process.pfNewIso+
    process.minipuppiNewIso+process.minitrkNewIso+process.minipfNewIso+
    process.NewPatMuonAnalyser)
