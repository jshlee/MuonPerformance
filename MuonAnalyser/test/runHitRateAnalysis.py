import FWCore.ParameterSet.Config as cms
process = cms.Process('HitRateAnalysis')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
# /DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO
process.source.fileNames.append("file:0224B4AB-31AF-E711-BD7E-0CC47A4C8EC6.root")

#print "fileNames: ", process.source.fileNames
process.options = cms.untracked.PSet()

process.TFileService = cms.Service("TFileService",fileName = cms.string("histo.root"))

process.load('SimMuon.GEMDigitizer.muonGEMPadDigis_cfi')
process.load('SimMuon.GEMDigitizer.muonGEMPadDigiClusters_cfi')
process.HitRateAnalysis = cms.EDAnalyzer('HitAnalysis',
    gemDigiInput      = cms.InputTag("simMuonGEMDigis"),
    gemPadDigiInput   = cms.InputTag("simMuonGEMCSCPadDigis"),
    gemCoPadDigiInput = cms.InputTag("simCscTriggerPrimitiveDigis"),
)

process.p = cms.Path(process.simMuonGEMPadDigis+process.simMuonGEMPadDigiClusters+process.HitRateAnalysis)


## for debugging
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
