import ROOT
import array, sys
from MuonPerformance.MuonAnalyser.histoHelper import *


def setMarkerStyle(h,color,style):
  h.SetMarkerColor(color)
  #h.SetMarkerStyle(style)
  h.SetMarkerStyle(25)
  h.SetMarkerSize(0.2)
  h.SetLineColor(color)
  h.SetLineWidth(2)
  h.SetLineStyle(style)


def drawSampleName(samplename, fX, fY, fSizeTex):
  tex2 = ROOT.TLatex()
  
  tex2.SetNDC()
  tex2.SetTextFont(62)
  tex2.SetTextSize(fSizeTex)
  
  for i, strLine in enumerate(samplename.split("\n")): 
    tex2.DrawLatex(fX, fY - i * 1.2 * fSizeTex, strLine)


strSigTDR = "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170806_1_zmm_200.root"
strBkgTDR = "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170806_1_qcd_200.root"
strSigTDR = "/cms/scratch/quark2930/Work/muon_upgrade/validation_forTDR/src/MuonPerformance/MuonAnalyser/test/samples_tdr_zmm_200_170806_1/test.root"
strBkgTDR = "/cms/scratch/quark2930/Work/muon_upgrade/validation_forTDR/src/MuonPerformance/MuonAnalyser/test/samples_tdr_qcd_200_170806_1/test.root"
strSigRun2 = "/xrootd/store/user/quark2930/muon_upgrade/Run2/run2_170815_1_zmm.root"
strBkgRun2 = "/xrootd/store/user/quark2930/muon_upgrade/Run2/run2_170815_1_qcd.root"
strSigRun2 = "/cms/scratch/quark2930/Work/muon_upgrade/validation_forTDR/src/MuonPerformance/MuonAnalyser/test/samples_run2_zmm_200_170815_1/out_415.root"
strBkgRun2 = "/cms/scratch/quark2930/Work/muon_upgrade/validation_forTDR/src/MuonPerformance/MuonAnalyser/test/samples_run2_qcd_200_170815_1/out_39.root"

fBkg = ROOT.TFile(strBkgTDR)
nNumEvtBkgTDR = fBkg.Get("PatMuonAnalyser/nevents").Integral()
fBkg.Close()

fBkg = ROOT.TFile(strBkgRun2)
nNumEvtBkgRun2 = fBkg.Get("PatMuonAnalyser/nevents").Integral()
fBkg.Close()

strCut = "muon.Pt() > 15 && 0.0 <= abs(muon.Eta()) && abs(muon.Eta()) < 2.4 && muon_isLoose"

### dics (infos) for strIdxType

dicExtraCut= {"Trk": "", "PF001": "0.1", "TrkRun2": "", "PFRun2": ""}

dicVar = {
  "Trk":     "muon_TrkIsolation03", 
  "PF001":   "muon_pfNewIsoPt04_ch / muon.Pt()", 
  "TrkRun2": "muon_TrkIsolation03", 
  "PFRun2":  "muon_PFIsolation03"
}

dicTitle = {
  "Trk": "I^{Track}, #sqrt{s} = 14 TeV, Phase-2 <PU> = 200", 
  "PF001": "I^{PF}, #sqrt{s} = 14 TeV, Phase-2 <PU> = 200", 
  "TrkRun2": "I^{Track}, Run 2, #sqrt{s} = 13 TeV", 
  "PFRun2": "I^{PF}, Run 2, #sqrt{s} = 13 TeV", 
}

dicColor = {"Trk": 2, "PF001": 4, "TrkRun2": 2, "PFRun2": 4}

dicLineStyle = {"Trk": 1, "PF001": 1, "TrkRun2": 2, "PFRun2": 2}

canv = makeCanvas("canv1", False)
setMargins(canv, False)

leg = ROOT.TLegend(0.16, 0.66, 0.45, 0.46)

arrGraph = []

fMin = 9999999999.0
fMax = 0.0

arrCurveName = ["Trk", "PF001", "TrkRun2", "PFRun2"]
#arrCurveName = ["TrkRun2", "PFRun2"]

arrBin = [5000, 0, 50.0]

print "%i is TDR, %i is Run2"%(nNumEvtBkgTDR, nNumEvtBkgRun2)

#for strCurve in dicTitle.keys():
for strCurve in arrCurveName:
  strSig = strSigTDR if "Run2" not in strCurve else strSigRun2
  strBkg = strBkgTDR if "Run2" not in strCurve else strBkgRun2
  
  strExtraCut = ""
  if strCurve == "PF001":
    strExtraCut = " && muon_pfNewIsoPt04_nh / muon.Pt() < " + dicExtraCut[ strCurve ]
  
  #strVar = "muon_pfNewIsoPt04_ch / muon.Pt()" if strCurve == "Trk" else "muon_TrkIsolation03"
  strVar = dicVar[ strCurve ]
  
  print strCurve + " (sig) : " + strSig
  print strCurve + " (bkg) : " + strBkg
  print strCurve + " (cut) : " + strCut + strExtraCut
  
  hSig = makeTH1(strSig, "PatMuonAnalyser/reco", "title", arrBin, strVar, strCut + " && muon_signal" + strExtraCut)
  hBkg = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", arrBin, strVar, strCut + strExtraCut)

  #nTotalSig = hSig.GetEntries() # another way, but same result
  #nTotalBkg = hBkg.GetEntries()
  nTotalSig = hSig.Integral(0, arrBin[ 0 ] + 1)
  nTotalBkg = hBkg.Integral(0, arrBin[ 0 ] + 1)
  
  nNumEvtBkg = nNumEvtBkgTDR if "Run2" not in strCurve else nNumEvtBkgRun2
  hBkgDen = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", arrBin, strVar, "muon_isLoose")
  #hBkgDen = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", arrBin, strVar, "1 == 1")
  nNumMuonBkg = hBkgDen.Integral(0, arrBin[ 0 ] + 1)
  
  print "Survived : %lf (sig), %lf - %lf (bkg)"%(nTotalSig, nTotalBkg, nNumMuonBkg)

  listSigPoint = []
  listBkgPoint = []

  nValSigCurr = 0
  nValBkgCurr = 0

  for i in range(arrBin[ 0 ] + 2): 
    #nValSigCurr += hSig.GetBinContent(i) # another way, but same result
    #nValBkgCurr += hBkg.GetBinContent(i)
    nValSigCurr = hSig.Integral(0, i)
    nValBkgCurr = hBkg.Integral(0, i)
    
    fSigPoint = nValSigCurr / nTotalSig
    #fBkgPoint = nValBkgCurr / nNumEvtBkg # total number of events
    fBkgPoint = nValBkgCurr / nNumMuonBkg # for fraction
    
    listSigPoint.append(fSigPoint)
    listBkgPoint.append(fBkgPoint)
    
    if fSigPoint >= 0.8 and fMin > fBkgPoint: 
      fMin = fBkgPoint
    
    if fSigPoint >= 0.8 and fMax < fBkgPoint: 
      fMax = fBkgPoint

  graphROC = ROOT.TGraph(arrBin[ 0 ] + 2, array.array("d", listSigPoint), array.array("d", listBkgPoint))

  graphROC.SetTitle("%s"%(dicTitle[ strCurve ]))
  setMarkerStyle(graphROC, dicColor[ strCurve ], dicLineStyle[ strCurve ])
  graphROC.GetXaxis().SetLimits(0.8, 1.02)
  
  arrGraph.append(graphROC)
  
  print "The curve named \"" + strCurve + "\" has been drawn" + " %i"%nNumEvtBkg

#### The following has to run only once!!!

for i in range(len(arrGraph)): 
  if i == 0: 
    arrGraph[ i ].GetXaxis().SetTitle("Signal efficiency")
    arrGraph[ i ].GetYaxis().SetTitle("Loose muon background fraction")
    arrGraph[ i ].GetXaxis().SetTitleSize(0.05)
    arrGraph[ i ].GetYaxis().SetTitleOffset(1.20)
    arrGraph[ i ].GetYaxis().SetTitleSize(0.05)

    arrGraph[ i ].GetXaxis().SetLabelSize(0.037)
    arrGraph[ i ].GetYaxis().SetLabelSize(0.037)

    #arrGraph[ i ].SetMinimum(fMin * 0.8)
    #arrGraph[ i ].SetMaximum(fMax * 1.2)
    arrGraph[ i ].SetMinimum(0.001)
    arrGraph[ i ].SetMaximum(0.08)

    arrGraph[ i ].Draw("")
  else:
    arrGraph[ i ].Draw("same")

  #### done of 'once'

  leg.AddEntry(arrGraph[ i ], arrGraph[ i ].GetTitle(), "pl")

###############################################################
### Finishing up
###############################################################

drawSampleName("Z/#gamma*#rightarrow#font[12]{#mu#mu} and QCD events, Loose Muon, \np_{T} > 15 GeV, 0.0 < |#eta| < 2.4", 0.17, 0.80, 0.033)

leg.SetTextFont(62)
leg.SetTextSize(0.030)
leg.SetBorderSize(0)
leg.Draw()

#CMS_lumi setting
iPos = 0
iPeriod = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
CMS_lumi.extraText = "Phase-2 Simulation"
CMS_lumi.lumi_sqrtS = ""
CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

canv.Modified()
canv.Update()
strOutFile = "ROCCurve_TDR_new_with_Run2_20170814.root"
canv.SaveAs(strOutFile)
canv.SaveAs(strOutFile.replace(".root", ".eps"))
canv.SaveAs(strOutFile.replace(".root", ".png"))
print strOutFile + " (and eps, png) has been drawn"


