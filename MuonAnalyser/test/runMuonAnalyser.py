import FWCore.ParameterSet.Config as cms
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process("MuonAnalyser",eras.Phase2C2)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))
run2 = False
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
      #'file:step3.root'
      '/store/relval/CMSSW_9_1_0_pre3/RelValZMM_14/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D13-v1/10000/00A9B1AB-642E-E711-9759-0025905B85CA.root',
      #'/store/relval/CMSSW_9_1_0_pre1/RelValZMM_14/GEN-SIM-RECO/90X_upgrade2023_realistic_v9_D4Timing-v1/00000/0C650B46-2719-E711-9B1C-0025905A60FE.root',  
      #'/store/relval/CMSSW_9_0_0_pre5/RelValZMM_14/GEN-SIM-RECO/90X_upgrade2023_realistic_v4_D4T-v1/00000/0E75D0DF-0501-E711-8018-0025905B85F6.root',  
      #'/store/relval/CMSSW_9_1_0_pre1/RelValZMM_13/GEN-SIM-RECO/PU25ns_90X_mcRun2_asymptotic_v5-v1/00000/40559096-AE10-E711-A8EE-0CC47A4C8ED8.root'
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_000.root',
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_001.root',
      #'file:/xrootd/store/user/jlee/CMSSW_9_0_0_pre5/me0seed/zmmD4/step3_002.root',
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
if not run2:
    process.muonAssociatorByHitsHelper.usePhase2Tracker = cms.bool(True)
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
    addPileupInfo = cms.InputTag("addPileupInfo"),
    muonLabel = cms.InputTag("muons"),
    muAssocLabel = cms.InputTag("muonAssociatorByHitsHelper"),
    tpSelector = muonTPSet,
    puppiIsolationChargedHadrons = cms.InputTag("muonIsolationPUPPI","h+-DR040-ThresholdVeto000-ConeVeto000"),
    puppiIsolationNeutralHadrons = cms.InputTag("muonIsolationPUPPI","h0-DR040-ThresholdVeto000-ConeVeto001"),
    puppiIsolationPhotons        = cms.InputTag("muonIsolationPUPPI","gamma-DR040-ThresholdVeto000-ConeVeto001"),
    puppiNoLepIsolationChargedHadrons = cms.InputTag("muonIsolationPUPPINoLep","h+-DR040-ThresholdVeto000-ConeVeto000"),
    puppiNoLepIsolationNeutralHadrons = cms.InputTag("muonIsolationPUPPINoLep","h0-DR040-ThresholdVeto000-ConeVeto001"),
    puppiNoLepIsolationPhotons        = cms.InputTag("muonIsolationPUPPINoLep","gamma-DR040-ThresholdVeto000-ConeVeto001"),    
)

process.MuonAnalyser.primaryVertex1D   = cms.InputTag('offlinePrimaryVertices')
process.MuonAnalyser.primaryVertex1DBS = cms.InputTag('offlinePrimaryVertices')
process.MuonAnalyser.primaryVertex4D   = cms.InputTag('offlinePrimaryVertices')
process.MuonAnalyser.primaryVertex4DBS = cms.InputTag('offlinePrimaryVertices')
process.MuonAnalyser.primaryVertexBS   = cms.InputTag('offlinePrimaryVertices')

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
                  coneSize = cms.double(0.4),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.0001),
                  isolateAgainst = cms.string('h+'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.4),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.01),
                  isolateAgainst = cms.string('h0'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.4),
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
                         +process.primaryVertexAssociation
                         +process.puppi
                         +process.particleFlowNoLep+process.puppiNoLep
                         +process.offlineSlimmedPrimaryVertices+process.packedPFCandidates
                         +process.muonIsolationPUPPI+process.muonIsolationPUPPINoLep
                         +process.MuonAnalyser)
