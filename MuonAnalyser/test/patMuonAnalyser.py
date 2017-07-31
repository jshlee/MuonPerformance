import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("PatMuonAnalyser",eras.Phase2C2)
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff')

# Beware, in this area the wild character is not working!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:0C9F61F7-5058-E711-9EDD-0CC47A0AD6AA.root'
        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/ZZTo4L_14TeV_powheg_pythia8/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v1/00000/FE139FC6-BF57-E711-8144-68B59972BFD8.root'
    ),
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("out.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
    
process.PatMuonAnalyser = cms.EDAnalyzer("PatMuonAnalyser",
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    addPileupInfo = cms.InputTag("slimmedAddPileupInfo"),
    muons = cms.InputTag("slimmedMuons"),
    pruned = cms.InputTag("prunedGenParticles"),
)


process.p = cms.Path(process.PatMuonAnalyser)
