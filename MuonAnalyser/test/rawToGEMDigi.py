# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2018_realistic -n 10 --era Run2_2018 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@ExtraHLT+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry DB:Extended --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECOGEM',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('EventFilter.GEMRawToDigi.muonGEMDigis_cfi')
process.load('RecoLocalMuon.GEMRecHit.gemLocalReco_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:EA53671A-89CD-E711-85FF-02163E011A18.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:reco.root'),
    outputCommands = process.RECOEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)
process.RECOoutput.outputCommands.append('drop *')
process.RECOoutput.outputCommands.append('keep *_GEM*_*_*')
process.RECOoutput.outputCommands.append('keep *_muons_*_*')
process.RECOoutput.outputCommands.append('keep *_globalMuons_*_*')
process.RECOoutput.outputCommands.append('keep *_standAloneMuons_*_*')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.muonGEMDigis+process.gemLocalReco)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)
# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.RECOoutput_step)
