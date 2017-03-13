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

"""
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
"""

# Beware, in this area the wild character is not working!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/pnfs/user/jlee/pfmuon/src/21211.0_TenMuExtendedE_0_200+TenMuExtendedE_0_200_pythia8_2023D4_GenSimHLBeamSpotFull+DigiFullTrigger_2023D4+RecoFullGlobal_2023D4+HARVESTFullGlobal_2023D4/step3.root'
      #'file:step3.root'
    ),
    skipBadFiles = cms.untracked.bool(True), 
)

#to run for entire sample
dir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/doc/9_0_0_pre2/TenMu_'
filelst = open(dir+"pu0.txt", "r")
#filelst = open(dir+"pu200.txt", "r")
#process.source.fileNames = filelst.readlines()

process.TFileService = cms.Service("TFileService",fileName = cms.string("out.root"))

process.load('SimMuon.MCTruth.muonAssociatorByHitsHelper_cfi')
process.muonAssociatorByHitsHelper.useGEMs = cms.bool(True)
process.muonAssociatorByHitsHelper.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis:Pixel")
process.muonAssociatorByHitsHelper.stripSimLinkSrc = cms.InputTag("simSiPixelDigis:Tracker")

from Validation.RecoMuon.selectors_cff import muonTPSet
process.MuonAnalyser = cms.EDAnalyzer("MuonAnalyser",
    primaryVertex     = cms.InputTag('offlinePrimaryVertices'),
    primaryVertex1D   = cms.InputTag('offlinePrimaryVertices1D'),
    primaryVertex1DBS = cms.InputTag('offlinePrimaryVertices1DWithBS'),
    primaryVertex4D   = cms.InputTag('offlinePrimaryVertices4D'),
    primaryVertex4DBS = cms.InputTag('offlinePrimaryVertices4DWithBS'),
    primaryVertexBS   = cms.InputTag('offlinePrimaryVerticesWithBS'),
    
    simLabel = cms.InputTag("mix","MergedTrackTruth"),
    simVertexCollection = cms.InputTag("g4SimHits"),
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
