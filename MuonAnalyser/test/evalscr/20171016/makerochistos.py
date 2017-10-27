import ROOT
import sys
from MuonPerformance.MuonAnalyser.histoHelper import *


dicStrSig = {
  "TDR":  "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_171023_2_zmm_200.root", 
  "Run2": "/xrootd/store/user/quark2930/muon_upgrade/Run2/run2_171026_1_zmm.root", 
  "RelVal": "/xrootd/store/user/quark2930/muon_upgrade/relval_zmm_0_171025_1.root", 
}

dicStrTTH = {
  "TDR":  "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_171019_1_ttbar_200.root", 
  "Run2": "/xrootd/store/user/quark2930/muon_upgrade/Run2/run2_171019_1_ttH.root", 
}

dicStrBkg = {
  "TDR":  "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_171023_2_qcd_200.root", 
  "Run2": "/xrootd/store/user/quark2930/muon_upgrade/Run2/run2_171026_1_qcd.root", 
  "RelVal": "/xrootd/store/user/quark2930/muon_upgrade/relval_qcd_0_171025_1.root", 
}

dicVar = {
  #"Trk":        "muon_trkNewIso", 
  "Trk":         "muon_TrkIsolation03", 
  "PF":          "muon_pfNewIsoPt02_ch / muon.Pt()", 
  "PF2":          "muon_pfNewIsoPt02_ch / muon.Pt()", 
  "NH":          "muon_pfNewIsoPt04_nh / muon.Pt()", 
  "PFbeta":      "muon_PFIsolation03", 
  "PUPPI":       "muon_puppiNewIso", 
  #"Trk_Run2":   "muon_trkNewIso", 
  "Trk_Run2":    "muon_TrkIsolation03", 
  "PFbeta_Run2": "muon_pfNewIso", 
  "PF_Run2":     "muon_pfNewIsoPt02_ch / muon.Pt()", 
  "PF2_Run2":     "muon_pfNewIsoPt02_ch / muon.Pt()", 
  #"PF_Run2":     "muon_PFIsolation03", 
  "PUPPI_Run2":  "muon_puppiNewIso", 
  "Trk_RelVal":  "muon_TrkIsolation03", 
  "PF_RelVal":   "muon_PFIsolation03", 
  
  "Trk_mini":       "muon_minitrkNewIso", 
  "PF_mini":        "muon_minipfNewIso", 
  "PUPPI_mini":     "muon_minipuppiNewIso", 
  "Trk_miniRun2":   "muon_minitrkNewIso", 
  "PF_miniRun2":    "muon_minipfNewIso", 
  "PUPPI_miniRun2": "muon_minipuppiNewIso", 
  
  "Trk_compmini":       "muon_TrkIsolation03", 
  "PF_compmini":        "muon_pfNewIso", 
  "PUPPI_compmini":     "muon_puppiNewIso", 
  "Trk_compminiRun2":   "muon_trkNewIso", 
  "PF_compminiRun2":    "muon_pfNewIso", 
  "PUPPI_compminiRun2": "muon_puppiNewIso", 
}

strIdx = sys.argv[ 1 ]
strSigBkg = sys.argv[ 2 ]

fBinMax = 50.0
nNumBin = 50000

strCut = "muon.Pt() > 15 && 0.0 <= abs(muon.Eta()) && abs(muon.Eta()) < 2.4 && muon_isLoose"
strExtraCut = ""

if "PF" == strIdx: strExtraCut = " && muon_pfNewIsoPt04_nh / muon.Pt() < 0.1"
if "PF2" == strIdx: strExtraCut = " && muon_pfNewIsoPt04_nh / muon.Pt() < 1.2"
if "PF_Run2" == strIdx: strExtraCut = " && muon_pfNewIsoPt04_nh / muon.Pt() < 0.001"
if "PF2_Run2" == strIdx: strExtraCut = " && muon_pfNewIsoPt04_nh / muon.Pt() < 0.00001"
if "mini" in strIdx: strExtraCut = " && jet_HTNorm >= 500 && met_Norm >= 200"
#if "Run2" not in strIdx: strExtraCut = strExtraCut.replace("Norm", "PUPPI")

strSig = dicStrSig[ "TDR" if "Run2" not in strIdx else "Run2" ]
strBkg = dicStrBkg[ "TDR" if "Run2" not in strIdx else "Run2" ]
if "mini" in strIdx: strSig = dicStrTTH[ "TDR" if "Run2" not in strIdx else "Run2" ]

if "RelVal" in strIdx: 
 strSig = dicStrSig[ "RelVal" ]
 strBkg = dicStrBkg[ "RelVal" ]

strVar = dicVar[ strIdx ]

hRes = ""

strTreename = "PatMuonAnalyser/reco" if "RelVal" not in strIdx else "MuonAnalyser/reco"

if strSigBkg == "sig" or strSigBkg == "Sig": 
  print strIdx + " (" + strSigBkg + "; " + strVar + ") : " + strCut + " && muon_signal" + strExtraCut
  hRes = makeTH1(strSig, strTreename, "title", [nNumBin, 0, fBinMax], strVar, strCut + " && muon_signal" + strExtraCut)
else:
  print strIdx + " (" + strSigBkg + "; " + strVar + ") : " + strCut + strExtraCut
  hRes = makeTH1(strBkg, strTreename, "title", [nNumBin, 0, fBinMax], strVar, strCut + strExtraCut)

hRes.SaveAs("makeroc_%s_%s.root"%(strIdx, strSigBkg))


