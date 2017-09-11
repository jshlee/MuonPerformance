import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def getMarker(filename, sig_tot, bkg_tot, cut, style):
    sig = getWeightedEntries(filename, treename, "muon_tmva_bdt", rangecut+"&&muon_signal&&%s"%cut)
    bkg = getWeightedEntries(filename, treename, "muon_tmva_bdt", rangecut+"&&!muon_signal&&%s"%cut)
    sigeff = sig/sig_tot
    bkgrej = (bkg_tot-bkg)/bkg_tot
    return ROOT.TMarker(sigeff, bkgrej, style)

def getRoc(h_sig, h_bkg):
    x=array.array('f',[])
    y=array.array('f',[])
    for k in range(101,0,-1):
        sigeff = h_sig.Integral(k,101)/h_sig.Integral(0,101)
        bkgrej = (h_bkg.Integral(0,101)-h_bkg.Integral(k,101))/h_bkg.Integral(0,101)
        x.append(sigeff)
        y.append(bkgrej)
        #print j, h_sig.GetBinCenter(j), sigeff, bkgrej, fake
    return ROOT.TGraph(len(x), x, y)

def getTmvaRoc(sample):
    tfile = ROOT.TFile("./tmva/TMVA_"+sample+".root")
    htmp = tfile.Get("dataset/Method_%s/%s/MVA_%s_rejBvsS"%(method,method,method))
    return copy.deepcopy(htmp)
    

datadir = '/xrootd/store/user/tt8888tt/muon/9_1_1_patch1/tmva/'
treename = "MuonAnalyser/reco"
rangecut = 'muon.Pt()>5&&abs(muon.Eta())<2.4'
method = "BDT"

samples= {"tenmu":[0.12,0.038], "ttbar":[0.08,0.02], "zmm":[0.03,-0.07]}
col = [[ROOT.kRed, ROOT.kOrange+7, ROOT.kPink+9], [ROOT.kBlue,ROOT.kViolet+7,ROOT.kAzure+7], [ROOT.kGreen,ROOT.kTeal+9,ROOT.kSpring+9]]

grs = []
hists = []
marks = []
for i, sample in enumerate(samples):
    for j, applied in enumerate(samples):
        samplename = sample+"_"+applied+"applied"
        if sample==applied: samplename=sample
        #else break
        filename = datadir+samplename+".root"
        print samplename

        h_sig = makeTH1(filename, treename, "muon_tmva_bdt", [100,-0.25,0.25], "muon_tmva_bdt", rangecut+"&&muon_signal")
        h_bkg = makeTH1(filename, treename, "muon_tmva_bdt", [100,-0.25,0.25], "muon_tmva_bdt", rangecut+"&&!muon_signal")

        gr = getRoc(h_sig, h_bkg) 
        gr.SetLineColor(col[i][j])
        gr.SetTitle("TMVA trained with "+applied.capitalize()+" on "+sample.capitalize())
        grs.append(gr)

        #tmva roc
        #h = getTmvaRoc(sample)
        #h.SetLineColor(col[i][j])
        #h.SetLineStyle(2)
        #h.SetLineWidth(3)
        #h.SetTitle(sample.split('/')[-1].split('.')[0]+" (TMVA)")
        #hists.append(h)

        ##id point
        #m = getMarker(filename, h_sig.Integral(0,101), h_bkg.Integral(0,101), "muon_tmva_bdt<%d"%samples[sample][0], 26)
        #m.SetMarkerSize(3)
        #m.SetMarkerColor(col[i][j])
        #marks.append(m)

h_init = ROOT.TH1F("","",1000,0.85,1.01)
h_init.SetMaximum(1.01)
h_init.SetMinimum(0.85)
h_init.GetXaxis().SetTitle("Signal Efficiency")
h_init.GetYaxis().SetTitle("Background Rejection")
h_init.GetYaxis().SetTitleOffset(1)

canv = makeCanvas("tmva", False)
canv.SetGrid()
setMargins(canv, False)
h_init.Draw()

leg = ROOT.TLegend(0.15,0.18,0.62,0.58)

for i, gr in enumerate(grs):
    gr.SetLineWidth(3)
    gr.Draw("Lsame")
    leg.AddEntry(gr,gr.GetTitle(),"l")

    hists[i].Draw("same")
    leg.AddEntry(hists[i],hists[i].GetTitle(),"l")


leg.SetTextFont(61)
leg.SetTextSize(0.028)
leg.Draw()

iPos = 0
iPeriod = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
CMS_lumi.extraText = "Simulation"
CMS_lumi.lumi_sqrtS = "14 TeV"
CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

canv.Modified()
canv.Update()
canv.SaveAs("tmva_roc.png")

