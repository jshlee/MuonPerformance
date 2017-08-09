import ROOT
import array
from MuonPerformance.MuonAnalyser.histoHelper import *


def setMarkerStyle(h,color,style):
  h.SetMarkerColor(color)
  h.SetMarkerStyle(style)
  h.SetMarkerSize(0.2)
  h.SetLineColor(color)
  h.SetLineWidth(2)


def drawSampleName(samplename, fX, fY, fSizeTex):
  tex2 = ROOT.TLatex()
  
  tex2.SetNDC()
  tex2.SetTextFont(62)
  tex2.SetTextSize(fSizeTex)
  
  for i, strLine in enumerate(samplename.split("\n")): 
    tex2.DrawLatex(fX, fY - i * 1.2 * fSizeTex, strLine)


strSig = "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170806_1_zmm_200.root"
strBkg = "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170806_1_qcd_200.root"

strCut = "muon.Pt() > 15 && 0.0 <= abs(muon.Eta()) && abs(muon.Eta()) < 2.4 && ( ( abs(muon.Eta()) < 2.4 && muon_isLoose ) || ( abs(muon.Eta()) > 2.4 && muon_isME0MuonLoose ) ) && abs(muon_poszPV0 - muon_poszMuon) < 0.5"

canv = makeCanvas("canv1", False)
setMargins(canv, False)

leg = ROOT.TLegend(0.16, 0.66, 0.45, 0.46)

###############################################################
### Track isolation
###############################################################

strVar = "muon_TrkIsolation03"
hSig = makeTH1(strSig, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)
hBkg = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)

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
  
  fSigPoint = nValSigCurr / nTotalSig
  fBkgPoint = nValBkgCurr / 2555040 # total number of events
  
  listSigPoint.append(fSigPoint)
  listBkgPoint.append(fBkgPoint)

graphROC1 = ROOT.TGraph(50000 + 2, array.array("d", listSigPoint), array.array("d", listBkgPoint))

graphROC1.SetTitle("Phase-2 <PU> = 200, I^{Track}, R = 0.3, #sqrt{s} = 14 TeV")
setMarkerStyle(graphROC1,  2, 25) # Color : Red
graphROC1.GetXaxis().SetLimits(0.0, 1.1)

#### The following has to run only once!!!

graphROC1.GetXaxis().SetTitle("Signal efficiency")
graphROC1.GetYaxis().SetTitle("Average Loose Muon Bkg Multiplicity")
graphROC1.GetXaxis().SetTitleSize(0.05)
graphROC1.GetYaxis().SetTitleOffset(1.20)
graphROC1.GetYaxis().SetTitleSize(0.05)

graphROC1.GetXaxis().SetLabelSize(0.037)
graphROC1.GetYaxis().SetLabelSize(0.037)

graphROC1.SetMinimum(0.0)
graphROC1.SetMaximum(0.72)

graphROC1.Draw("")

#### done of 'once'

leg.AddEntry(graphROC1, graphROC1.GetTitle(), "pl")

###############################################################
### PF isolation with no beta correction
###############################################################

strVar = "( muon_pfNewIsoPt06_ch + muon_pfNewIsoPt04_nh + muon_pfNewIso_ph ) / muon.Pt()"
hSig = makeTH1(strSig, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)
hBkg = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)

#nTotalSig = hSig.GetEntries()
#nTotalBkg = hBkg.GetEntries()
nTotalSig = hSig.Integral(0, 50000 + 1)
nTotalBkg = hBkg.Integral(0, 50000 + 1)

listSigPoint = []
listBkgPoint = []

nValSigCurr = 0
nValBkgCurr = 0

for i in range(50000 + 2): 
  #nValSigCurr += hSig.GetBinContent(i)
  #nValBkgCurr += hBkg.GetBinContent(i)
  nValSigCurr = hSig.Integral(0, i)
  nValBkgCurr = hBkg.Integral(0, i)
  
  fSigPoint = nValSigCurr / nTotalSig
  fBkgPoint = nValBkgCurr / 2555040 # total number of events
  
  listSigPoint.append(fSigPoint)
  listBkgPoint.append(fBkgPoint)

graphROC2 = ROOT.TGraph(50000 + 2, array.array("d", listSigPoint), array.array("d", listBkgPoint))

graphROC2.SetTitle("Phase-2 <PU> = 200, I^{PF} best, R = 0.3, #sqrt{s} = 14 TeV")
setMarkerStyle(graphROC2,  8, 25) # Color : Green
graphROC2.GetXaxis().SetLimits(0.0, 1.1)

graphROC2.Draw("same")
leg.AddEntry(graphROC2, graphROC2.GetTitle(), "pl")

###############################################################
### PF isolation with beta correction
###############################################################

# Because in this fomula max function does not work, I have to use that conditional operator
strVar = "( muon_pfNewIsoPt06_ch + ( muon_pfNewIsoPt04_nh + muon_pfNewIso_ph - 0.5 * muon_PFIso03PUPt >= 0.0 ? muon_pfNewIsoPt04_nh + muon_pfNewIso_ph - 0.5 * muon_PFIso03PUPt : 0.0 ) ) / muon.Pt()"
hSig = makeTH1(strSig, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)
hBkg = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)

#nTotalSig = hSig.GetEntries()
#nTotalBkg = hBkg.GetEntries()
nTotalSig = hSig.Integral(0, 50000 + 1)
nTotalBkg = hBkg.Integral(0, 50000 + 1)

listSigPoint = []
listBkgPoint = []

nValSigCurr = 0
nValBkgCurr = 0

for i in range(50000 + 2): 
  #nValSigCurr += hSig.GetBinContent(i)
  #nValBkgCurr += hBkg.GetBinContent(i)
  nValSigCurr = hSig.Integral(0, i)
  nValBkgCurr = hBkg.Integral(0, i)
  
  fSigPoint = nValSigCurr / nTotalSig
  fBkgPoint = nValBkgCurr / 2555040 # total number of events
  
  listSigPoint.append(fSigPoint)
  listBkgPoint.append(fBkgPoint)

graphROC3 = ROOT.TGraph(50000 + 2, array.array("d", listSigPoint), array.array("d", listBkgPoint))

graphROC3.SetTitle("Phase-2 <PU> = 200, I^{PF}+#beta best, R = 0.3, #sqrt{s} = 14 TeV")
setMarkerStyle(graphROC3,  1, 25) # Color : Black
graphROC3.GetXaxis().SetLimits(0.0, 1.1)

graphROC3.Draw("same")
leg.AddEntry(graphROC3, graphROC3.GetTitle(), "pl")

###############################################################
### PUPPI isolation
###############################################################

strVar = "muon_puppiNewIso"
hSig = makeTH1(strSig, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)
hBkg = makeTH1(strBkg, "PatMuonAnalyser/reco", "title", [50000, 0, 3.0], strVar, strCut)

#nTotalSig = hSig.GetEntries()
#nTotalBkg = hBkg.GetEntries()
nTotalSig = hSig.Integral(0, 50000 + 1)
nTotalBkg = hBkg.Integral(0, 50000 + 1)

listSigPoint = []
listBkgPoint = []

nValSigCurr = 0
nValBkgCurr = 0

for i in range(50000 + 2): 
  #nValSigCurr += hSig.GetBinContent(i)
  #nValBkgCurr += hBkg.GetBinContent(i)
  nValSigCurr = hSig.Integral(0, i)
  nValBkgCurr = hBkg.Integral(0, i)
  
  fSigPoint = nValSigCurr / nTotalSig
  fBkgPoint = nValBkgCurr / 2555040 # total number of events
  
  listSigPoint.append(fSigPoint)
  listBkgPoint.append(fBkgPoint)

graphROC4 = ROOT.TGraph(50000 + 2, array.array("d", listSigPoint), array.array("d", listBkgPoint))

graphROC4.SetTitle("Phase-2 <PU> = 200, I^{PUPPI}, R = 0.3, #sqrt{s} = 14 TeV")
setMarkerStyle(graphROC4,  4, 25) # Color : Blue
graphROC4.GetXaxis().SetLimits(0.0, 1.1)

graphROC4.Draw("same")
leg.AddEntry(graphROC4, graphROC4.GetTitle(), "pl")

###############################################################
### Finishing up
###############################################################

drawSampleName("Z/#gamma*#rightarrow#font[12]{#mu#mu} and QCD events, Loose Muon, \np_{T} > 15 GeV, 0.0 < |#eta| < 2.4, |z_{reco} - z_{sim}| < 0.5 cm", 0.17, 0.80, 0.033)

leg.SetTextFont(62)
leg.SetTextSize(0.033)
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
canv.SaveAs("ROCCurve_TDR_new.png")


