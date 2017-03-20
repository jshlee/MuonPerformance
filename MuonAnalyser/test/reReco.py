# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --processName RERECO --conditions auto:phase2_realistic --era Phase2C2_timing --eventcontent FEVTDEBUGHLT --runUnscheduled -s RECO --datatier GEN-SIM-RECO --geometry Extended2023D4 -n -1 --filein /store/relval/CMSSW_9_0_0_pre4/RelValZMM_14/GEN-SIM-RECO/90X_upgrade2023_realistic_v3_2023D4Timing-v1/10000/04013655-A9EC-E611-9783-0CC47A4C8E98.root --fileout step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RERECO',eras.Phase2C2_timing)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_9_0_0_pre4/RelValZMM_14/GEN-SIM-RECO/90X_upgrade2023_realistic_v3_2023D4Timing-v1/10000/04013655-A9EC-E611-9783-0CC47A4C8E98.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
slimOutputCommands = cms.untracked.vstring('drop *_*_*_*',
'keep *_simSiPixelDigis_*_*',
'keep *_siPixelClusters_*_*',
'keep *_genParticles_*_*',
'keep *_electronGsfTracks_*_*',
'keep *_*Muon*_*_*',
'keep *_*Segments*_*_*',
'keep *_*GEM*_*_*',
'keep *_*ME0*_*_*',
'keep *_*DT*_*_*',
'keep *_*CSC*_*_*',
'keep *_*RPC*_*_*',
'keep *_*offlinePrimaryVertices*_*_*',
'keep *_*particleFlow*_*_*',
'keep TrackingParticles_*_*_*',
'keep TrackingVertexs_*_*_*',
'keep *Muons_*_*_*',
'keep SimVertexs_g4SimHits_*_*',
'keep SimTracks_g4SimHits_*_*',
'keep PSimHits_g4SimHits_TrackerHits*_*',
'keep *_*_*Muon*_*',
'keep *_dt1DRecHits_*_*',
'keep recoTrack*_*_*_*',
'keep *_ak4PFJets_*_*',
'keep *_inclusiveCandidateSecondaryVertices*_*_*',
'keep TrackingRecHits*_*_*_*',
'keep CrossingFrame*_*_*_*',
)

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(10485760),
    fileName = cms.untracked.string('step3.root'),
    outputCommands = slimOutputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Path and EndPath definitions
process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
