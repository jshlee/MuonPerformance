import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("PatMuonAnalyser",eras.Phase2C2)
process.load('FWCore.MessageService.MessageLogger_cfi')

# Beware, in this area the wild character is not working!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:0C9F61F7-5058-E711-9EDD-0CC47A0AD6AA.root'
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
