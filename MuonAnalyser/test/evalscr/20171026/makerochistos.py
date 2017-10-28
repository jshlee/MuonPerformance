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
  "run2Trk": "muon_TrkIsolation03", 
  "run2PF":  "muon_pfNewIsoPt02_ch / muon.Pt()", 
}

strIdx = sys.argv[ 1 ]
strSigBkg = sys.argv[ 2 ]
strPrefix = "" if len(sys.argv) < 4 else sys.argv[ 3 ]

strVar = dicVar[ strIdx.split("_")[ 0 ] ]
fCut = int(strIdx.split("_")[ 1 ]) * 0.005

fBinMax = 50.0
nNumBin = 50000

strCut = "muon.Pt() > 15 && 0.0 <= abs(muon.Eta()) && abs(muon.Eta()) < 2.4 && muon_isLoose"
strExtraCut = ""

if "PF" in strIdx and len(sys.argv) < 4: 
  strExtraCut = " && muon_pfNewIsoPt04_nh / muon.Pt() < %0.8f"%(fCut)

strSig = dicStrSig[ "Run2" ]
strBkg = dicStrBkg[ "Run2" ]

hRes = ""

strTreename = "PatMuonAnalyser/reco" if "RelVal" not in strIdx else "MuonAnalyser/reco"

if strSigBkg == "sig" or strSigBkg == "Sig": 
  #print strIdx + " (" + strSigBkg + "; " + strVar + ") : " + strCut + " && muon_signal" + strExtraCut
  hRes = makeTH1(strSig, strTreename, "title", [nNumBin, 0, fBinMax], strVar, strCut + " && muon_signal" + strExtraCut)
  #print "Sig (%s%s) : %lf"%(strIdx, strPrefix, hRes.GetIntegral())
else:
  #print strIdx + " (" + strSigBkg + "; " + strVar + ") : " + strCut + strExtraCut
  hRes = makeTH1(strBkg, strTreename, "title", [nNumBin, 0, fBinMax], strVar, strCut + strExtraCut)

hRes.SaveAs("makeroc_%s_%s%s.root"%(strIdx, strSigBkg, strPrefix))


