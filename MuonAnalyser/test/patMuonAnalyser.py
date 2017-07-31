import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("PatMuonAnalyser",eras.Phase2C2)
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff')

# Beware, in this area the wild character is not working!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:FE139FC6-BF57-E711-8144-68B59972BFD8.root'
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/ZZTo4L_14TeV_powheg_pythia8/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v1/00000/FE139FC6-BF57-E711-8144-68B59972BFD8.root'
    ),
)

process.load('RecoMuon.MuonIsolation.muonIsolationPUPPI_cff')
process.muonIsolationMiniAOD = process.muonIsolationMiniAODPUPPI.clone(usePUPPI = cms.bool(False))

process.TFileService = cms.Service("TFileService",fileName = cms.string("out.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
    
process.PatMuonAnalyser = cms.EDAnalyzer("PatMuonAnalyser",
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    addPileupInfo = cms.InputTag("slimmedAddPileupInfo"),
    muons = cms.InputTag("slimmedMuons"),
    pruned = cms.InputTag("prunedGenParticles"),
    puppiNoLepIsolationChargedHadrons = cms.InputTag("muonIsolationMiniAODPUPPINoLeptons","h+-DR030-ThresholdVeto000-ConeVeto000"),
    puppiNoLepIsolationNeutralHadrons = cms.InputTag("muonIsolationMiniAODPUPPINoLeptons","h0-DR030-ThresholdVeto000-ConeVeto001"),
    puppiNoLepIsolationPhotons        = cms.InputTag("muonIsolationMiniAODPUPPINoLeptons","gamma-DR030-ThresholdVeto000-ConeVeto001"),    
    pfIsolationChargedHadrons = cms.InputTag("muonIsolationMiniAOD","h+-DR030-ThresholdVeto000-ConeVeto000"),
    pfIsolationNeutralHadrons = cms.InputTag("muonIsolationMiniAOD","h0-DR030-ThresholdVeto000-ConeVeto001"),
    pfIsolationPhotons        = cms.InputTag("muonIsolationMiniAOD","gamma-DR030-ThresholdVeto000-ConeVeto001"),    

)

process.p = cms.Path(process.muonIsolationMiniAOD+process.muonIsolationMiniAODPUPPINoLeptons+process.PatMuonAnalyser)
