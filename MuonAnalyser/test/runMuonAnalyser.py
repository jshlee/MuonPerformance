import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("MuonAnalyser",eras.Phase2C2_timing)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2023_realistic_v1', '')

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
      #'file:/cms/scratch/jlee/muonHisto/TTbar_13TeV_TuneCUETP8M1_2023D1.root' #TTbar
      #'file:/cms/scratch/jlee/muonHisto/QCD_Pt_600_800_13TeV_TuneCUETP8M1_2023D1.root' #QCD
      #'file:/cms/scratch/jlee/muonHisto/TenMuExtendedE_0_200_pythia8_2023D1.root' #pu0 in ~doc/my/QCDandTTbar.txt
      #'file:/xrootd/store/relval/CMSSW_9_0_0_pre2/RelValQCD_Pt-20toInf_MuEnrichedPt15_14TeV/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4Timing-v1/10000/442A9B3F-56C3-E611-A341-0CC47A4D75F4.root'
    ),
    skipBadFiles = cms.untracked.bool(True), 
)


#run for entire sample
dir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/doc/my/'
filelst = open(dir+"TenMuPU0.txt", "r")
#filelst = open(dir+"pu200.txt", "r")
process.source.fileNames = filelst.readlines()

#process.TFileService = cms.Service("TFileService",fileName = cms.string("pu0/TTbar_pu0.root"))
process.TFileService = cms.Service("TFileService",fileName = cms.string("TenMuPU0.root"))

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

process.load('CommonTools.PileupAlgos.Puppi_cff')
process.pfNoLepPUPPI = cms.EDFilter("PdgIdCandViewSelector",
                                    src = cms.InputTag("particleFlow"), 
                                    pdgId = cms.vint32( 1,2,22,111,130,310,2112,211,-211,321,-321,999211,2212,-2212 )
                                    )
process.puppiNoLep = process.puppi.clone()
process.puppiNoLep.candName = cms.InputTag('pfNoLepPUPPI') 
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
process.load("PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi")
process.load("PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi")

process.p = cms.Path(process.packedPFCandidates+process.muonAssociatorByHitsHelper+process.MuonAnalyser)
