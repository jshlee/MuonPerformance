import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('SliceTestAnalysis',eras.Run2_2018,eras.run3_GEM)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#process.maxEvents.input = cms.untracked.int32(10)
# Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.skipEvents = cms.untracked.uint32(0)

### Files available as at 17.7.2018
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_000.root 319337
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_349.root 319347
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_975.root 319348
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_993.root 319349
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_1140.root 319449
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_1510.root 319450
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_1518.root 319456
# Error in <TFile::TFile>: file /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_1808.root does not exist
# Skipping /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_1808.root
# /xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3_1857.root 319459


#process.source.fileNames.append('/store/data/Run2018B/Cosmics/AOD/PromptReco-v1/000/317/428/00000/E4BC1D7B-3F6A-E811-9E05-FA163E57A064.root')
from glob import glob
process.source.fileNames.extend(
    ['/store/user/jlee/SingleMuon/Run2018C-v1/RECOv2/step3_019.root']
    # ['file:/xrootd/store/user/iawatson/SingleMuon/Run2018C-v1/wGEM/RECOv1/step3_{0:03d}.root'.format(i) for i in range(682, 682+1)]
    # ['file:'+f for f in glob('/xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3*.root')][:]
)
#process.source.fileNames.append('/store/user/jlee/SingleMuon/Run2017G-v1/FEVTEvent/step3_313.root')

#fname = 'singleMuon.txt'
#f = open(fname)
#for line in f:
#    process.source.fileNames.append(line)

process.options = cms.untracked.PSet()

process.TFileService = cms.Service("TFileService",fileName = cms.string("histo2018.root"))

process.SliceTestAnalysis = cms.EDAnalyzer('SliceTestAnalysis',
    process.MuonServiceProxy,
    gemRecHits = cms.InputTag("gemRecHits"),
    muons = cms.InputTag("muons"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
)
process.p = cms.Path(process.SliceTestAnalysis)
