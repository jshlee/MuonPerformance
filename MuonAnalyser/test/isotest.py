import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

datadir = '/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170806_1_' 
rangecut = 'muon.Pt()>15&&abs(muon.Eta())<2.8&&abs(muon.Eta())>2.4&&muon_isME0MuonLoose&&abs(muon_poszPV0-muon_poszMuon)<0.5'
#rangecut = 'muon.Pt()>15&&abs(muon.Eta())<2.8&&abs(muon.Eta())>2.4&&muon_isLoose&&abs(muon_poszPV0-muon_poszMuon)<0.5'

tfile = ROOT.TFile(datadir+'qcd_200.root')
nevents = tfile.Get("PatMuonAnalyser/nevents").Integral()
treename = "PatMuonAnalyser"

iso = "muon_pfNewIso"
#iso = "muon_puppiNewIso"
grs = []
for i, pfminPt in enumerate(["", "Pt02", "Pt04", "Pt06", "Pt08", "Pt10"]):
    print pfminPt
    x=array.array('f',[])
    y=array.array('f',[])
    isocut = iso+pfminPt+"_ch"
    h_zmm = makeTH1(datadir+'zmm_200.root', treename+"/gen", iso, [1000,0,35], isocut, rangecut)
    h_qcd = makeTH1(datadir+'qcd_200.root', treename+"/reco", iso, [1000,0,35], isocut, rangecut)
    for j in range(1,1000):
        sigeff = h_zmm.Integral(1,j)/h_zmm.Integral(0,1001)
        bkgrej = h_qcd.Integral(1,j)/float(nevents)
        x.append(sigeff)
        y.append(bkgrej)
        print j, sigeff, bkgrej
    gr = ROOT.TGraph(len(x), x, y)
    gr.SetLineColor(i+1)
    gr.SetLineWidth(2)
    if (i+1==5) or (i+1==6): gr.SetLineColor(i+1+3)
    gr.SetTitle(isocut)
    grs.append(gr)

plotrange = (0.8,1.02)
rMIN = plotrange[0]
rMAX = plotrange[1]
h_init = ROOT.TH1F("","",1000,rMIN,rMAX)
h_init.SetMaximum(grs[-1].GetHistogram().GetMaximum()*1.01)
h_init.SetMinimum(grs[-1].GetHistogram().GetMinimum()*0.8)
h_init.GetXaxis().SetTitle("Signal Efficiency")
h_init.GetYaxis().SetTitle("Average Loose Muoon Bkg Multiplicity")
h_init.GetYaxis().SetTitleOffset(1.2)
h_init.GetYaxis().SetLabelSize(0.03)
h_init.GetYaxis().SetTitleSize(0.045)

canv = makeCanvas("tmva", False)
setMargins(canv, False)
h_init.Draw()

leg = ROOT.TLegend(0.17,0.57,0.6,0.8)

for i, h in enumerate(grs):
    grs[i].Draw("Lsame")
    leg.AddEntry(grs[i],grs[i].GetTitle(),"l")

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
canv.SaveAs("isotest_%s.png"%iso)

