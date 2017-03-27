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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
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
      #'/store/relval/CMSSW_9_0_0_pre5/RelValZMM_14/GEN-SIM-RECO/90X_upgrade2023_realistic_v4_D4T-v1/00000/0E75D0DF-0501-E711-8018-0025905B85F6.root',  
      #'file:/cms/home/jlee/scratch/gemSeeding/src/crab/step3.root'
      'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_000.root',
      'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_001.root',
      'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_002.root',
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_003.root',
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_004.root',
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_005.root',
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_006.root',
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_007.root',
      #'file:/xrootd/store/user/jlee/RelValZMM_14/crab_zmmD4PU140_me0seed/170319_133447/0000/step3_1.root'
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
    puppiIsolationChargedHadrons = cms.InputTag("muonIsolationPUPPI","h+-DR030-ThresholdVeto000-ConeVeto000"),
    puppiIsolationNeutralHadrons = cms.InputTag("muonIsolationPUPPI","h0-DR030-ThresholdVeto000-ConeVeto001"),
    puppiIsolationPhotons        = cms.InputTag("muonIsolationPUPPI","gamma-DR030-ThresholdVeto000-ConeVeto001"),
    puppiNoLepIsolationChargedHadrons = cms.InputTag("muonIsolationPUPPINoLep","h+-DR030-ThresholdVeto000-ConeVeto000"),
    puppiNoLepIsolationNeutralHadrons = cms.InputTag("muonIsolationPUPPINoLep","h0-DR030-ThresholdVeto000-ConeVeto001"),
    puppiNoLepIsolationPhotons        = cms.InputTag("muonIsolationPUPPINoLep","gamma-DR030-ThresholdVeto000-ConeVeto001"),    
)

process.MuonAnalyser.tpSelector.maxRapidity = cms.double(3.0)
process.MuonAnalyser.tpSelector.minRapidity = cms.double(-3.0)

process.load('CommonTools.PileupAlgos.Puppi_cff')
process.particleFlowNoLep = cms.EDFilter("PdgIdCandViewSelector",
                                    src = cms.InputTag("particleFlow"), 
                                    pdgId = cms.vint32( 1,2,22,111,130,310,2112,211,-211,321,-321,999211,2212,-2212 )
                                    )
process.puppiNoLep = process.puppi.clone(candName = cms.InputTag('particleFlowNoLep'))

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
process.load("PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi")
process.load("PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi")

IsoConeDefinitions = cms.VPSet(
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.3),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.0001),
                  isolateAgainst = cms.string('h+'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.3),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.01),
                  isolateAgainst = cms.string('h0'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.3),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.01),
                  isolateAgainst = cms.string('gamma'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),                  
)

process.muonIsolationPUPPI = cms.EDProducer( "CITKPFIsolationSumProducerForPUPPI",
                srcToIsolate = cms.InputTag("muons"),
                srcForIsolationCone = cms.InputTag('packedPFCandidates'),
                puppiValueMap = cms.InputTag(''),
                usePUPPINoLepton = cms.bool(False),
                isolationConeDefinitions = IsoConeDefinitions
)
process.muonIsolationPUPPINoLep = process.muonIsolationPUPPI.clone(usePUPPINoLepton = cms.bool(True))

process.p = cms.Path(process.muonAssociatorByHitsHelper
                         +process.particleFlowNoLep+process.puppiNoLep
                         +process.packedPFCandidates
                         +process.muonIsolationPUPPI+process.muonIsolationPUPPINoLep
                         +process.MuonAnalyser)
