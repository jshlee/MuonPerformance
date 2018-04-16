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
# /RelValZMM_14/CMSSW_9_3_7-PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/GEN-SIM-DIGI-RAW
process.source.fileNames.append("file:0034A16F-1831-E811-BB69-5065F381B271.root")
process.source.inputCommands=cms.untracked.vstring('keep *', 'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT','drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT','drop l1tEMTFHit2016s_simEmtfDigis__HLT','drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT','drop l1tEMTFTrack2016s_simEmtfDigis__HLT')

process.options = cms.untracked.PSet()

process.TFileService = cms.Service("TFileService",fileName = cms.string("histo.root"))

process.load('SimMuon.GEMDigitizer.muonGEMPadDigis_cfi')
process.load('SimMuon.GEMDigitizer.muonGEMPadDigiClusters_cfi')
process.HitRateAnalysis = cms.EDAnalyzer('HitAnalysis',
    gemDigiInput      = cms.InputTag("simMuonGEMDigis"),
    gemPadDigiInput   = cms.InputTag("simMuonGEMPadDigis"),
    gemCoPadDigiInput = cms.InputTag("simCscTriggerPrimitiveDigis"),
)
process.p = cms.Path(#process.simMuonGEMPadDigis+process.simMuonGEMPadDigiClusters+
                         process.HitRateAnalysis)


## for debugging
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
