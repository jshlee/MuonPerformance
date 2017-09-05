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

etarange = sys.argv[1]
channel = 'ch'
iso = "muon_pfNewIso"
#iso = "muon_puppiNewIso" #"muon_pfNewIso"

datadir = '/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170806_1_'
treename = "PatMuonAnalyser"
tfile = ROOT.TFile(datadir+'qcd_200.root')
nevents = tfile.Get("PatMuonAnalyser/nevents").Integral()

grs = []
#for i, pfminPt in enumerate(["", "Pt02", "Pt04", "Pt06", "Pt08", "Pt10","TrkIsolation03"]):
for i, pfminPt in enumerate(["TrkIsolation03", "pfNewIsoPt04_ch"]):#pf
#for i, pfminPt in enumerate(["TrkIsolation03", "0.01", "0.05","0.1","0.2","0.4"]):#pf
#for i, pfminPt in enumerate(["0.001", "0.002", "0.005", "0.008", "0.01", "0.015","TrkIsolation03"]):#puppi
#for i, pfminPt in enumerate(["pfNewIsoPt04_ch", "0.1", "TrkIsolation03"]):
#for i, pfminPt in enumerate(["0.1", "puppiNewIsoPt04_ch", "TrkIsolation03"]):
    print pfminPt
    x=array.array('f',[])
    y=array.array('f',[])

    if etarange == 'high':
        rangecut = 'muon.Pt()>15&&abs(muon.Eta())<2.8&&abs(muon.Eta())>2.4&&muon_isME0MuonLoose'
        text = "Loose Muon, \np_{T} > 15 GeV, 2.4 < |#eta| < 2.8"
    if etarange == 'low':
        rangecut = 'muon.Pt()>15&&abs(muon.Eta())<2.4&&muon_isLoose'
        text = "Loose Muon, \np_{T} > 15 GeV, 0.0 < |#eta| < 2.4"

    if "TrkIsolation03" in pfminPt:
        isocut = "muon_TrkIsolation03/muon.Pt()"
    elif "puppi" in pfminPt:
        isocut = "muon_puppiNewIsoPt04_ch/muon.Pt()"
    else:
        isocut = "muon_pfNewIsoPt04_ch/muon.Pt()"
        #rangecut = rangecut+'&&(muon_pfNewIsoPt04_nh)/muon.Pt()<'+pfminPt
        
    h_zmm = makeTH1(datadir+'zmm_200.root', treename+"/gen", isocut, [10000,0,5], isocut, rangecut)
    h_qcd = makeTH1(datadir+'qcd_200.root', treename+"/reco", isocut, [10000,0,5], isocut, rangecut)
    if h_zmm.Integral(0,10001)==0: continue
    for j in range(0,10001):
        sigeff = h_zmm.Integral(0,j)/h_zmm.Integral(0,10001)
        bkgrej = h_qcd.Integral(0,j)/float(nevents)
        x.append(sigeff)
        y.append(bkgrej)
        #print j, sigeff, bkgrej
    gr = ROOT.TGraph(len(x), x, y)
    gr.SetLineColor(i+1)
    if (i+1==5) or (i+1==6): gr.SetLineColor(i+1+3)
    #gr.SetTitle(isocut)
    if "pf" in isocut: gr.SetTitle("I^{PF}, Pt 0.4")
    #if "pf" in isocut: gr.SetTitle("I^{PF}, NH rel < "+pfminPt)
    else : gr.SetTitle("I^{Track}")
    grs.append(gr)

h_init = ROOT.TH1F("","",1000,0.8,1.02)
#h_init.SetMaximum(grs[-1].GetHistogram().GetMaximum()*1.01)
h_init.SetMaximum(grs[-1].GetHistogram().GetMaximum())
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
#canv.SaveAs("isoroc_%s_%s_%s_nhph.png"%(iso, channel, etarange))
#canv.SaveAs("isoroc_%s.png"%(etarange))
canv.SaveAs("isoroc_optimized_%s.png"%(etarange))

