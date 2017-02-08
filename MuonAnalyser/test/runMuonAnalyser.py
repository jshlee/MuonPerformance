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

# Beware, in this area the wild character is not working!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'file:/cms/scratch/jlee/upgradeMuonReco/reco.root'
      #'file:/cms/scratch/jlee/upgradeMuonReco/RelValZMM_13_reco.root'
      #'file:/cms/scratch/jlee/upgradeMuonReco/RelValTenMuExtendedE_0_200_reco.root'
      #'file:/cms/scratch/jlee/muonHisto/TTbar_13TeV_TuneCUETP8M1_2023D1.root' #TTbar
      #'file:/cms/scratch/jlee/muonHisto/QCD_Pt_600_800_13TeV_TuneCUETP8M1_2023D1.root' #QCD
      #'file:/cms/scratch/jlee/muonHisto/TenMuExtendedE_0_200_pythia8_2023D1.root' #pu0 in ~doc/my/QCDandTTbar.txt
      #'file:086EEBBD-67C2-E611-BD7D-0CC47A4D761A.root'
      #'file:/cms/scratch/quark2930/Work/muon_upgrade/CMSSW_9_0_0_pre4/src/step3.root'
      #'file:/xrootd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v1_2023D4TimingPU200-v1/10000/*.root'
      #'file:/xrootd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v1_2023D4TimingPU200-v1/10000/0042B3C2-54C4-E611-9D39-0025905B85DA.root', 
      #'file:/xrootd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v1_2023D4TimingPU200-v1/10000/F8C53CD1-40C4-E611-B456-0025905A60CA.root', 
      #'file:/xrootd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v1_2023D4TimingPU200-v1/10000/E00CBAFE-66C4-E611-8525-0CC47A7C346E.root'
      #'file:/xrootd/store/relval/CMSSW_9_0_0_pre2/RelValQCD_Pt-20toInf_MuEnrichedPt15_14TeV/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4Timing-v1/10000/442A9B3F-56C3-E611-A341-0CC47A4D75F4.root', 
      #'file:/xrootd/store/relval/CMSSW_9_0_0_pre2/RelValQCD_Pt-20toInf_MuEnrichedPt15_14TeV/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4Timing-v1/10000/C0EE80D3-55C3-E611-B660-0CC47A7C351E.root'
      'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4Timing-v1/10000/086EEBBD-67C2-E611-BD7D-0CC47A4D761A.root', 
      'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4Timing-v1/10000/1E75363B-69C2-E611-9F20-0025905A6084.root', 
      'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4Timing-v1/10000/A440BAB7-67C2-E611-BFC3-0CC47A4D76D2.root', 
      'root://cms-xrdr.sdfarm.kr:1094///xrd/store/relval/CMSSW_9_0_0_pre2/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4Timing-v1/10000/D6BD7163-69C2-E611-9D20-0CC47A4D75F4.root'
    ),
    skipBadFiles = cms.untracked.bool(True), 
)

"""
dir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/'
filelst = open(dir+"test/TenMuPU0.txt", "r")
#filelst = open(dir+"pu200.txt", "r")
#process.source.fileNames = filelst.readlines()
print filelst.readlines()

# If you want to use multiple source files, use the following codes.
#run for entire sample
dir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/'
filelst = open(dir+"test/TenMuPU0.txt", "r")
#filelst = open(dir+"pu200.txt", "r")
process.source.fileNames = filelst.readlines()
"""

#process.TFileService = cms.Service("TFileService",fileName = cms.string("pu0/TTbar_pu0.root"))
process.TFileService = cms.Service("TFileService",fileName = cms.string("out_1.root"))

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
