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
        #'file:/cms/home/jlee/scratch/test/gemSeeding/src/24811.0_TenMuExtendedE_0_200+TenMuExtendedE_0_200_pythia8_2023D13_GenSimHLBeamSpotFull+DigiFullTrigger_2023D13+RecoFullGlobal_2023D13+HARVESTFullGlobal_2023D13/step3.root'
        'file:/cms/home/jlee/scratch/test/baseline/src/26211.0_TenMuExtendedE_0_200+TenMuExtendedE_0_200_pythia8_2023D14_GenSimHLBeamSpotFull+DigiFull_2023D14+RecoFullGlobal_2023D14+HARVESTFullGlobal_2023D14/step3.root'
    ),
    skipBadFiles = cms.untracked.bool(True), 
)

#to run for entire sample
#dir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/doc/9_1_1/'
#filelst = open(dir+"zmm.txt", "r")
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
