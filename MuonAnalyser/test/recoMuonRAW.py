import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('reRECO',eras.Run2_2018,eras.run3_GEM)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:B2AA1EF5-E0AA-E811-BD11-FA163ED8B27F.root'),
    fileNames = cms.untracked.vstring('/store/data/Run2018C/SingleMuon/RAW/v1/000/319/347/00000/5CCC3C8F-E982-E811-80E1-FA163EAE9F8D.root'),
    secondaryFileNames = cms.untracked.vstring()
)
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'cert_306824-306826_5TeV.txt').getVLuminosityBlockRange()
#print process.source.lumisToProcess

process.options = cms.untracked.PSet()
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.gemSkim = cms.EDFilter("GEMSkim",gemRecHits = cms.InputTag("gemRecHits"))
process.GEMRecHitSkim = cms.Path(process.gemSkim)

# Output definition
from EventFilter.ScalersRawToDigi.Scalers_EventContent_cff import *
from EventFilter.OnlineMetaDataRawToDigi.OnlineMetaData_EventContent_cff import *
from EventFilter.Utilities.Tcds_EventContent_cff import *
process.RecoMuonRECO.outputCommands.extend(EvtScalersAOD.outputCommands)
process.RecoMuonRECO.outputCommands.extend(OnlineMetaDataContent.outputCommands)
process.RecoMuonRECO.outputCommands.extend(TcdsEventContent.outputCommands)
process.RecoMuonRECO.outputCommands.append("keep *_muonGEMDigis_*_*")
process.RecoMuonRECO.outputCommands.append("keep *_gemRecHits_*_*")
process.RecoMuonRECO.outputCommands.append("keep *_csc2DRecHits_*_*")

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('AOD.root'),
    outputCommands = process.RecoMuonRECO.outputCommands,
    splitLevel = cms.untracked.int32(0),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('GEMRecHitSkim'))
)

# Other statements 101X_dataRun2_Prompt_v11
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')
# using local emap db
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        connect = cms.string('sqlite_fip:MuonPerformance/MuonAnalyser/data/GEMELMap.db'),
        record = cms.string('GEMELMapRcd'),
        tag = cms.string('GEMELMap_v3')
    ))

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.highlevelreco = cms.Sequence(process.egammaHighLevelRecoPrePF*
                             process.particleFlowReco*
                             #process.egammaHighLevelRecoPostPF*
                             process.muoncosmichighlevelreco*
                             process.muonshighlevelreco)
process.reconstruction = cms.Sequence(process.localreco*process.globalreco
                                          *process.highlevelreco
                                          )

process.reconstruction_step = cms.Path(process.reconstruction)
#process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)

process.muonGEMDigis.unPackStatusDigis = cms.bool(True)
# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,##process.L1Reco_step,
                                    process.reconstruction_step,
                                    ##process.GEMRecHitSkim,process.endjob_step,
                                    process.RECOoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2018 

#call to customisation function customisePostEra_Run2_2018 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2018(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
