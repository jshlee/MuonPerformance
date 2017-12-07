import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("NewPatMuonAnalyser",eras.Phase2)
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_0_0_pre1/RelValTTbar_13/MINIAODSIM/94X_mc2017_realistic_v10-v1/10000/7CBB8BF1-9FCB-E711-8FBD-0CC47A78A3F8.root',
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_0_0_pre1/RelValTTbar_13/MINIAODSIM/94X_mc2017_realistic_v10-v1/10000/C65FC0F1-9FCB-E711-9C97-0CC47A78A3F8.root'
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_0_0_pre1/RelValZMM_13/MINIAODSIM/PU25ns_94X_mc2017_realistic_v10-v1/10000/92D71AF9-EFCB-E711-BC66-0CC47A7C34B0.root',
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_0_0_pre1/RelValZMM_13/MINIAODSIM/PU25ns_94X_mc2017_realistic_v10-v1/10000/A00FB1F9-EFCB-E711-BD66-0025905A60B6.root'
        'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_0_0_pre1/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_94X_mc2017_realistic_v10-v1/10000/92B83DCE-5CCC-E711-86FA-0025905B85C0.root',
        'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_0_0_pre1/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_94X_mc2017_realistic_v10-v1/10000/B6FA61CB-5CCC-E711-AC5F-0025905A6118.root'
    ),
)

process.load('RecoMuon.MuonIsolation.muonIsolationPUPPI_cff')

process.puppiNewIso = process.muonIsolationMiniAODPUPPINoLeptons.clone()
process.trkNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))
#process.trkNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False), useOnlyTrack = cms.bool(True))
process.pfNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))

#process.minipuppiNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(True))
#process.minitrkNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))
#process.minipfNewIso = process.puppiNewIso.clone(usePUPPI = cms.bool(False))

process.TFileService = cms.Service("TFileService",fileName = cms.string("outnew.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.NewPatMuonAnalyser = cms.EDAnalyzer("NewPatMuonAnalyser",
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    addPileupInfo = cms.InputTag("slimmedAddPileupInfo"),
    muons = cms.InputTag("slimmedMuons"),
    pruned = cms.InputTag("prunedGenParticles"),
    
    tmvaWeightLabel   = cms.string('MuonPerformance/MuonAnalyser/src/TMVAClassification_BDT.weights.xml'),

    puppiNewIso_ch = cms.InputTag("puppiNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    puppiNewIso_nh = cms.InputTag("puppiNewIso","h0-DR030-ThresholdVeto000-ConeVeto001"),
    puppiNewIso_ph = cms.InputTag("puppiNewIso","gamma-DR030-ThresholdVeto000-ConeVeto001"),    
    trkNewIso = cms.InputTag("trkNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    pfNewIso_ch = cms.InputTag("pfNewIso","h+-DR030-ThresholdVeto000-ConeVeto000"),
    pfNewIso_nh = cms.InputTag("pfNewIso","h0-DR030-ThresholdVeto000-ConeVeto001"),
    pfNewIso_ph = cms.InputTag("pfNewIso","gamma-DR030-ThresholdVeto000-ConeVeto001"),
  )

process.p = cms.Path(
    process.puppiNewIso+process.trkNewIso+process.pfNewIso+
    #process.minipuppiNewIso+process.minitrkNewIso+process.minipfNewIso+
    process.NewPatMuonAnalyser)
