import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("MuonAnalyser",eras.Phase2C1)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

#process.MessageLogger.categories.append("MuonAnalyser")
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.destinations = cms.untracked.vstring("cout","junk")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string("DEBUG"),
    default = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
    FwkReport = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    #MuonAnalyser   = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    #MuonAnalyser_Matching = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'file:/cms/scratch/jlee/upgradeMuonReco/reco.root'
      #'file:/cms/scratch/jlee/upgradeMuonReco/RelValZMM_13_reco.root'
      #'file:/cms/scratch/jlee/upgradeMuonReco/RelValTenMuExtendedE_0_200_reco.root'
      #'file:/xrootd/store/user/jlee/RelValTenMuExtendedE_0_200/crab_20161014_004717/161013_154753/0000/*.root' #pu0
      #'file:/xrootd/store/user/jlee/RelValTenMuExtendedE_0_200/crab_20161014_004624/161013_154703/0000/*.root' #pu200
      'file:/cms/scratch/jlee/muonHisto/TTbar_13TeV_TuneCUETP8M1_2023D1.root' #TTbar
      #'file:/cms/scratch/jlee/muonHisto/QCD_Pt_600_800_13TeV_TuneCUETP8M1_2023D1.root' #QCD
      #'file:/cms/scratch/jlee/muonHisto/TenMuExtendedE_0_200_pythia8_2023D1.root' #pu0 in ~doc/my/QCDandTTbar.txt
    ),
    skipBadFiles = cms.untracked.bool(True), 
)

"""
#run for entire sample
dir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/doc/'
filelst = open(dir+"pu0.txt", "r")
#filelst = open(dir+"pu200.txt", "r")
process.source.fileNames = filelst.readlines()
"""
process.TFileService = cms.Service("TFileService",fileName = cms.string("pu0/TTbar_pu0.root"))
#process.TFileService = cms.Service("TFileService",fileName = cms.string("out.root"))

process.load('SimMuon.MCTruth.muonAssociatorByHitsHelper_cfi')
process.muonAssociatorByHitsHelper.useGEMs = cms.bool(True)
process.muonAssociatorByHitsHelper.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis:Pixel")
process.muonAssociatorByHitsHelper.stripSimLinkSrc = cms.InputTag("simSiPixelDigis:Tracker")

from Validation.RecoMuon.selectors_cff import muonTPSet
process.MuonAnalyser = cms.EDAnalyzer("MuonAnalyser",
    primaryVertex = cms.InputTag('offlinePrimaryVertices'),
    simLabel = cms.InputTag("mix","MergedTrackTruth"),
    muonLabel = cms.InputTag("muons"),
    muAssocLabel = cms.InputTag("muonAssociatorByHitsHelper"),
    tpSelector = muonTPSet, 
)
process.MuonAnalyser.tpSelector.maxRapidity = cms.double(3)
process.MuonAnalyser.tpSelector.minRapidity = cms.double(-3)

process.p = cms.Path(process.muonAssociatorByHitsHelper+process.MuonAnalyser)





############## To make a output root file ###############

#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.output = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string(
#        'file:out_test.root'
#    ),
#    outputCommands = cms.untracked.vstring(
#        'keep  *_*_*_*',
#    ),
#)
#process.out_step     = cms.EndPath(process.output)
#process.p = cms.Path(process.gemMuonSel*process.muonAssociatorByHits*process.MuonAnalyser)
#process.p = cms.Path(process.gemMuonSel)
#process.out_step = cms.EndPath(process.output)
#process.schedule = cms.Schedule(
#    process.p,
#    process.out_step
#)
