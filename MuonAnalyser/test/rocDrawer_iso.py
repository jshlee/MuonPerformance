import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def drawSampleName(samplename, fX, fY, fSizeTex):
  tex2 = ROOT.TLatex()
  
  tex2.SetNDC()
  tex2.SetTextFont(62)
  tex2.SetTextSize(fSizeTex)
  
  for i, strLine in enumerate(samplename.split("\n")): 
    tex2.DrawLatex(fX, fY - i * 1.2 * fSizeTex, strLine)

def getRoc(pfminPt, rangecut):
    print pfminPt
    x=array.array('f',[])
    y=array.array('f',[])

    h_zmm = makeTH1(datadir+'zmm.root', treename+"/gen", isocut, [10000,0,80], pfminPt, rangecut)
    h_qcd = makeTH1(datadir+'qcd.root', treename+"/reco", isocut, [10000,0,80], pfminPt, rangecut)

    for j in range(0,10001):
        sigeff = h_zmm.Integral(0,j)/h_zmm.Integral(0,10001)
        bkgrej = h_qcd.Integral(0,j)/float(nevents)
        x.append(sigeff)
        y.append(bkgrej)
        #print j, sigeff, bkgrej
    gr = ROOT.TGraph(len(x), x, y)
    gr.SetTitle(pfminPt)
    return gr

datadir = './'
treename = "NewPatMuonAnalyser"
rangecut = 'muon.Pt()>15&&abs(muon.Eta())<2.4&&muon_isLoose'
text = "Loose Muon, \np_{T} > 15 GeV, 0.0 < |#eta| < 2.4"

tfile = ROOT.TFile(datadir+'qcd.root')
nevents = tfile.Get("NewPatMuonAnalyser/nevents").Integral()

grs = []
for pfminPt in ["muon_miniIso_ch", "muon_miniIso_nh", "muon_miniIso_ph"]:
    grs.append(getRoc(pfminPt, rangecut))

h_init = ROOT.TH1F("","",1000,0.8,1.02)
h_init.SetMaximum(grs[1].GetHistogram().GetMaximum()*1.01)
h_init.GetXaxis().SetTitle("Signal Efficiency")
h_init.GetYaxis().SetTitle("Average Loose Muoon Bkg Multiplicity")
h_init.GetYaxis().SetLabelSize(0.038)
h_init.GetYaxis().SetTitleSize(0.048)

canv = makeCanvas("tmva", False)
setMargins(canv, False)
h_init.Draw()

drawSampleName(text, 0.17, 0.8, 0.033)
leg = ROOT.TLegend(0.15,0.48,0.5,0.68)

for i, gr in enumerate(grs):
    gr.SetLineWidth(2)
    gr.SetLineColor(i+1)
    gr.Draw("Lsame")
    leg.AddEntry(gr,gr.GetTitle(),"l")

leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(61)
leg.SetTextSize(0.03)
leg.Draw()

iPos = 0
iPeriod = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
CMS_lumi.extraText = "Simulation"
CMS_lumi.lumi_sqrtS = "14 TeV"
CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

canv.Modified()
canv.Update()
canv.SaveAs("roc_iso.png")

