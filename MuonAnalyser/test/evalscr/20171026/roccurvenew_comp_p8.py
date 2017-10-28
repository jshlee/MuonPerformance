import ROOT
import array, sys
from MuonPerformance.MuonAnalyser.histoHelper import *


def setMarkerStyle(h,color,style,width=2):
  h.SetMarkerColor(color)
  #h.SetMarkerStyle(style)
  h.SetMarkerStyle(25)
  h.SetMarkerSize(0.2)
  h.SetLineColor(color)
  h.SetLineWidth(width)
  h.SetLineStyle(style)


def drawSampleName(samplename, fX, fY, fSizeTex):
  tex2 = ROOT.TLatex()
  
  tex2.SetNDC()
  tex2.SetTextFont(62)
  tex2.SetTextSize(fSizeTex)
  
  for i, strLine in enumerate(samplename.split("\n")): 
    tex2.DrawLatex(fX, fY - i * 1.2 * fSizeTex, strLine)


strSigTDR = "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_171015_1_zmm_200.root"
strBkgTDR = "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_171015_1_qcd_200.root"
strSigRun2 = "/xrootd/store/user/quark2930/muon_upgrade/Run2/run2_171015_1_zmm.root"
strBkgRun2 = "/xrootd/store/user/quark2930/muon_upgrade/Run2/run2_171015_1_qcd.root"
strSigRelVal = "/xrootd/store/user/quark2930/muon_upgrade/relval_zmm_0_171025_1.root"
strBkgRelVal = "/xrootd/store/user/quark2930/muon_upgrade/relval_qcd_0_171025_1.root"

dicTitle = {
  "Trk": "I^{Track}, Phase2 <PU> = 200, #sqrt{s} = 14 TeV", 
  "PF": "I^{PF} (new), Phase2 <PU> = 200, #sqrt{s} = 14 TeV", 
  "PFbeta": "I^{PF} (#beta), Phase2 <PU> = 200, #sqrt{s} = 14 TeV", 
  "PUPPI": "I^{PUPPI}, Phase2 <PU> = 200, #sqrt{s} = 14 TeV", 
  "Trk_Run2": "I^{Track}, Run 2, <PU> = 27, #sqrt{s} = 13 TeV", 
  "PF_Run2": "I^{PF}, Run 2, <PU> = 27, #sqrt{s} = 13 TeV", 
  "PUPPI_Run2": "I^{PUPPI}, Run 2, <PU> = 27, #sqrt{s} = 13 TeV", 
  "Trk_RelVal": "I^{Track}, RelVal, <PU> = 0, #sqrt{s} = 14 TeV", 
  "PF_RelVal": "I^{PF}, RelVal, <PU> = 0, #sqrt{s} = 14 TeV", 
  
  "Trk_mini": "I^{mini Track}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "PF_mini": "I^{mini PF}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "PUPPI_mini": "I^{mini PUPPI}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "Trk_miniRun2": "I^{mini Track}, Run 2, #sqrt{s} = 13 TeV", 
  "PF_miniRun2": "I^{mini PF}, Run 2, #sqrt{s} = 13 TeV", 
  "PUPPI_miniRun2": "I^{mini PUPPI}, Run 2, #sqrt{s} = 13 TeV", 
  
  "Trk_compmini": "I^{Track}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "PF_compmini": "I^{PF} (new), <PU> = 200, #sqrt{s} = 14 TeV", 
  "PUPPI_compmini": "I^{PUPPI}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "Trk_compminiRun2": "I^{Track}, Run 2, #sqrt{s} = 13 TeV", 
  "PF_compminiRun2": "I^{PF}, Run 2, #sqrt{s} = 13 TeV", 
  "PUPPI_compminiRun2": "I^{PUPPI}, Run 2, #sqrt{s} = 13 TeV", 
}

"""
dicTitle = {
  "Trk": "I^{Track}, R = 3.0", 
  "PF": "I^{PF} (new), R = 3.0", 
  "PFbeta": "I^{PF} (#beta), R = 3.0", 
  "PUPPI": "I^{PUPPI}, R = 3.0", 
  "Trk_Run2": "I^{Track}, R = 3.0", 
  "PF_Run2": "I^{PF}, R = 3.0", 
  "PUPPI_Run2": "I^{PUPPI}, R = 3.0", 
  
  "Trk_mini": "I^{mini Track}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "PF_mini": "I^{mini PF}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "PUPPI_mini": "I^{mini PUPPI}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "Trk_miniRun2": "I^{mini Track}, Run 2, #sqrt{s} = 13 TeV", 
  "PF_miniRun2": "I^{mini PF}, Run 2, #sqrt{s} = 13 TeV", 
  "PUPPI_miniRun2": "I^{mini PUPPI}, Run 2, #sqrt{s} = 13 TeV", 
  
  "Trk_compmini": "I^{Track}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "PF_compmini": "I^{PF} (new), <PU> = 200, #sqrt{s} = 14 TeV", 
  "PUPPI_compmini": "I^{PUPPI}, <PU> = 200, #sqrt{s} = 14 TeV", 
  "Trk_compminiRun2": "I^{Track}, Run 2, #sqrt{s} = 13 TeV", 
  "PF_compminiRun2": "I^{PF}, Run 2, #sqrt{s} = 13 TeV", 
  "PUPPI_compminiRun2": "I^{PUPPI}, Run 2, #sqrt{s} = 13 TeV", 
}
"""

dicVar = {
  "Trk":        "muon_trkNewIso", 
  "PF":         "muon_pfNewIso", 
  "PUPPI":      "muon_puppiNewIso", 
  "Trk_Run2":   "muon_trkNewIso", 
  "PF_Run2":    "muon_pfNewIso", 
  "PUPPI_Run2": "muon_puppiNewIso", 
  
  "Trk_mini":       "muon_minitrkNewIso", 
  "PF_mini":        "muon_minipfNewIso", 
  "PUPPI_mini":     "muon_minipuppiNewIso", 
  "Trk_miniRun2":   "muon_minitrkNewIso", 
  "PF_miniRun2":    "muon_minipfNewIso", 
  "PUPPI_miniRun2": "muon_minipuppiNewIso", 
}

dicColor = {"Trk": 1, "PF": 2, "PFbeta": 8, "PUPPI": 3}

dicLineStyle = {"": 1, "Run2": 2, "mini": 3, "miniRun2": 39, "compmini": 1, "compminiRun2": 0}

arrLineColor = [1, 2, 2, 1, 2, 4, 8, 8, 20]
arrLineStyle = [1, 1, 1, 2, 2, 1, 1, 1, 1]

#leg = ROOT.TLegend(0.16, 0.66, 0.45, 0.46)
leg = ROOT.TLegend(0.45, 0.44, 0.84, 0.24)

fMin = 9999999999.0
fMax = 0.0

#arrCurveName = ["Trk", "PFbeta", "PF", "Trk_Run2", "PF_Run2", "Trk_RelVal", "PF_RelVal"]
#arrCurveName = ["Trk", "PFbeta", "PUPPI", "PF"]
#arrCurveName = ["Trk_Run2", "PF_Run2", "PUPPI_Run2"]
#arrCurveName = ["Trk", "PF", "PUPPI", "Trk_mini", "PF_mini", "PUPPI_mini"]
#arrCurveName = ["Trk_compminiRun2", "PF_compminiRun2", "PUPPI_compminiRun2", "Trk_miniRun2", "PF_miniRun2", "PUPPI_miniRun2"]
#arrCurveName = ["Trk", "PF"]
arrCurveName = ["run2Trk", sys.argv[ 1 ]]

fBkg = ROOT.TFile(strBkgTDR)
nNumEvtBkgTDR  = fBkg.Get("PatMuonAnalyser/nevents").Integral()
fBkg.Close()

fBkg = ROOT.TFile(strBkgRun2)
nNumEvtBkgRun2 = fBkg.Get("PatMuonAnalyser/nevents").Integral()
fBkg.Close()
#nNumEvtBkgRun2 = 1.0

fBkg = ROOT.TFile(strBkgRelVal)
nNumEvtBkgRelVal = fBkg.Get("MuonAnalyser/nevents").Integral()
fBkg.Close()

strCut = "muon.Pt() > 15 && 0.0 <= abs(muon.Eta()) && abs(muon.Eta()) < 2.4 && ( ( abs(muon.Eta()) <= 2.4 && muon_isLoose ) || ( abs(muon.Eta()) > 2.4 && muon_isME0MuonLoose ) )"

### dics (infos) for strIdxType

"""
dicExtraCut= {"Trk": "", "PF001": "0.1", "PF005": "0.05", "PF01": "0.1"}
dicTitle = {
    "Trk": "I^{Track}", 
    "PF001": "I^{PF}", 
}
dicColor = {"Trk": 1, "PF001": 2, "PF005": 4, "PF01": 8}
"""

canv = makeCanvas("canv1", False)
setMargins(canv, False)

leg = ROOT.TLegend(0.16, 0.72, 0.45, 0.46)

arrGraph = []

fMin = 9999999999.0
fMax = 0.0

nIdxLine = 0

dicListSigPoint = {}
dicListBkgPoint = {}

#for strCurve in dicTitle.keys():
for strCurve in arrCurveName:
  strExtraCut = ""
  #if strCurve != "Trk":
  #  strExtraCut = " && muon_pfNewIsoPt04_nh / muon.Pt() < " + dicExtraCut[ strCurve ]
  
  #strVar = "muon_pfNewIsoPt04_ch / muon.Pt()" if strCurve != "Trk" else "muon_TrkIsolation03"
  
  #hSig = makeTH1(strSig, "PatMuonAnalyser/reco", "title", [50000, 0, 50.0], strVar, strCut + " && muon_signal" + strExtraCut)
  #hBkg = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", [50000, 0, 50.0], strVar, strCut + strExtraCut)
   
  fileSig = ROOT.TFile("makeroc_" + strCurve + "_Sig.root")
  fileBkg = ROOT.TFile("makeroc_" + strCurve + "_Bkg.root")
  
  hSig = copy.deepcopy(fileSig.Get("temp"))
  hBkg = copy.deepcopy(fileBkg.Get("temp"))

  #nTotalSig = hSig.GetEntries() # another way, but same result
  #nTotalBkg = hBkg.GetEntries()
  nTotalSig = hSig.Integral(0, 50000 + 1)
  nTotalBkg = hBkg.Integral(0, 50000 + 1)

  listSigPoint = []
  listBkgPoint = []

  nValSigCurr = 0
  nValBkgCurr = 0

  for i in range(50000 + 2): 
    #nValSigCurr += hSig.GetBinContent(i) # another way, but same result
    #nValBkgCurr += hBkg.GetBinContent(i)
    nValSigCurr = hSig.Integral(0, i)
    nValBkgCurr = hBkg.Integral(0, i)
    
    nNumEvtBkg = nNumEvtBkgTDR if "Run2" not in strCurve else nNumEvtBkgRun2
    if "RelVal" in strCurve: nNumEvtBkg = nNumEvtBkgRelVal
    
    fSigPoint = nValSigCurr / nTotalSig
    fBkgPoint = nValBkgCurr / nNumEvtBkg # total number of events
    
    listSigPoint.append(fSigPoint)
    listBkgPoint.append(fBkgPoint)
    
    if fSigPoint >= 0.8 and fMin > fBkgPoint: 
      fMin = fBkgPoint
    
    if fSigPoint >= 0.8 and fMax < fBkgPoint: 
      fMax = fBkgPoint
  
  dicListSigPoint[ strCurve ] = listSigPoint
  dicListBkgPoint[ strCurve ] = listBkgPoint

fXMin = 0.8
fXMax = 1.02
fYMin = fMin * 0.8
fYMax = fMax * 1.2
fYMin = 0.024
fMin = 0.024 / 0.8

for strCurve in arrCurveName:
  listSigPointCurr = []
  listBkgPointCurr = []
  
  fSigBox = ( fXMax - fXMin ) / 100
  fBkgBox = ( fYMax - fYMin ) / 100
  
  fSigPointPrev = dicListSigPoint[ strCurve ][ 0 ]
  fBkgPointPrev = dicListBkgPoint[ strCurve ][ 0 ]
  
  listSigPointCurr.append(fSigPointPrev)
  listBkgPointCurr.append(fBkgPointPrev)
  
  for i in range(len(dicListSigPoint[ strCurve ])): 
    if i == 0: continue
    
    fSigPointCurr = dicListSigPoint[ strCurve ][ i ]
    fBkgPointCurr = dicListBkgPoint[ strCurve ][ i ]
    
    if abs(fSigPointCurr - fSigPointPrev) > fSigBox or abs(fBkgPointCurr - fBkgPointPrev) > fBkgBox:
      listSigPointCurr.append(fSigPointCurr)
      listBkgPointCurr.append(fBkgPointCurr)
      
      fSigPointPrev = fSigPointCurr
      fBkgPointPrev = fBkgPointCurr
  
  nNumPoints = len(listSigPointCurr)
  
  graphROC = ROOT.TGraph(nNumPoints, array.array("d", listSigPointCurr), array.array("d", listBkgPointCurr))
  
  #arrPartStrCurve = strCurve.split("_")
  #strCurveP1 = arrPartStrCurve[ 0 ]
  #strCurveP2 = "" if len(arrPartStrCurve) <= 1 else arrPartStrCurve[ 1 ]

  strTitle = dicTitle[ "Trk" ] if "Trk" in strCurve else dicTitle[ "PF" ]
  arrCurveName = ["Trk", sys.argv[ 1 ].split("_")[ 0 ].replace("run2", "")]
  
  #graphROC.SetTitle("%s"%(dicTitle[ strCurve ]))
  graphROC.SetTitle(strTitle)
  #setMarkerStyle(graphROC, dicColor[ strCurveP1 ], dicLineStyle[ strCurveP2 ])
  #setMarkerStyle(graphROC, dicColor[ strCurveP1 ] + dicLineStyle[ strCurveP2 ], 1)
  #setMarkerStyle(graphROC, arrLineColor[ nIdxLine ], arrLineStyle[ nIdxLine ])
  setMarkerStyle(graphROC, arrLineColor[ nIdxLine ], arrLineStyle[ nIdxLine ], 2 if strCurve != "PF" else 3)
  nIdxLine += 1
  graphROC.GetXaxis().SetLimits(fXMin, fXMax)
  
  arrGraph.append(graphROC)

#### The following has to run only once!!!

for i in range(len(arrGraph)): 
  if i == 0: 
    arrGraph[ i ].GetXaxis().SetTitle("Signal efficiency")
    arrGraph[ i ].GetYaxis().SetTitle("Average Loose Muon Bkg Multiplicity")
    arrGraph[ i ].GetXaxis().SetTitleSize(0.05)
    arrGraph[ i ].GetYaxis().SetTitleOffset(1.20)
    arrGraph[ i ].GetYaxis().SetTitleSize(0.05)

    arrGraph[ i ].GetXaxis().SetLabelSize(0.037)
    arrGraph[ i ].GetYaxis().SetLabelSize(0.037)

    arrGraph[ i ].SetMinimum(fMin * 0.8)
    arrGraph[ i ].SetMaximum(fMax * 1.2)

    arrGraph[ i ].Draw("")
  else:
    arrGraph[ i ].Draw("same")

  #### done of 'once'

  leg.AddEntry(arrGraph[ i ], arrGraph[ i ].GetTitle(), "pl")

###############################################################
### Finishing up
###############################################################

drawSampleName("Z/#gamma*#rightarrow#font[12]{#mu#mu} and QCD events, Loose Muon, p_{T} > 15 GeV, 0.0 #leq |#eta| < 2.4", 0.17, 0.80, 0.033)
#drawSampleName("Z/#gamma*#rightarrow#font[12]{#mu#mu} and QCD events, Loose Muon, #sqrt{s} = 14 TeV, Phase-2 <PU> = 200, \np_{T} > 15 GeV, 0.0 #leq |#eta| < 2.4", 0.17, 0.80, 0.030)
#drawSampleName("Z/#gamma*#rightarrow#font[12]{#mu#mu} and QCD events, Loose Muon, #sqrt{s} = 13 TeV, Run II, \np_{T} > 15 GeV, 0.0 #leq |#eta| < 2.4", 0.17, 0.80, 0.033)
#drawSampleName("Z/#gamma*#rightarrow#font[12]{#mu#mu} and QCD events, Loose Muon, #sqrt{s} = 13 TeV, Run II, \np_{T} > 15 GeV, 0.0 #leq |#eta| < 2.4, HT #geq 500, MET #geq 200", 0.17, 0.80, 0.033)
#drawSampleName("ttH and QCD events, Loose Muon, #sqrt{s} = 13 TeV, Run II, \np_{T} > 15 GeV, 0.0 #leq |#eta| < 2.4, HT #geq 500, MET #geq 200", 0.17, 0.80, 0.033)

leg.SetTextFont(62)
leg.SetTextSize(0.030)
leg.SetBorderSize(0)
leg.Draw()

#CMS_lumi setting
iPos = 0
iPeriod = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
CMS_lumi.extraText = "Simulation"
CMS_lumi.lumi_sqrtS = ""
CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

canv.Modified()
canv.Update()
strOutFile = "plots/20171026/run2rocs/ROCCurve_TDR_new_for_page8_20171026_%s.root"%(sys.argv[ 1 ])
canv.SaveAs(strOutFile)
canv.SaveAs(strOutFile.replace("root", "png"))
print strOutFile + " has been drawn"


