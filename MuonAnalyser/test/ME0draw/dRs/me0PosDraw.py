#!/usr/bin/env python

import ROOT, copy
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
from MuonPerformance.MuonAnalyser.histoHelper import *

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.3, 0.8, samplename)

#filenames = ["../../pu0/TTbar_pu0","../../pu0/QCD_pu0","../../pu0/pu0","../../pu200/pu200"]
filenames = ["../../ZMM_PU200_pre5/ZMM_PU200_pre5.root"]
gentree = "MuonAnalyser/gen"
recotree = "MuonAnalyser/reco"

binning = [ [41,0,4],      # delta x 
            [21,0,20],     # delta y
            [41,0,4],      # pull x
            [21,0,20],     # pull y
            [81,0,0.8],    # delta dxdz
            [21,0,2],      # delta dydz
            [21,0,1.0],    # delta phi
            [8,0,7] ]      # noRecHit

var = [ "muon_ME0deltaX",
        "muon_ME0deltaY",
        "muon_ME0pullX",
        "muon_ME0pullY",
        "muon_ME0deltaDXDZ",
        "muon_ME0deltaDYDZ",
        "muon_ME0dPhi",
        "muon_ME0noRecHit" ]

tightCut = [ 1, 4, 1.5, 3, 0.05, 0.3, 0.05, 3]
looseCut = [ 2, 8, 3, 8, 0.075, 0.4, 0.1, 2 ]

for filename in filenames:
    pileup = "pu200"
   
    for i,plotvar in enumerate(var):
        ### Canvas ###
        canvasname = ""+pileup+plotvar
        c = makeCanvas(canvasname, False)
        if "Err" in plotvar:
            l = ROOT.TLegend(.8, .8, .95, .90)
        else:
            l = ROOT.TLegend(.85, .85, .95, .95)
        c.SetLogy()
        #c.SetMinimum(0.1)
        tightIndicator = ROOT.TLine(tightCut[i],0,tightCut[i],1)
        tightIndicator.SetLineStyle(2)
        tightIndicator.SetLineColor(2)
        tightIndicator.SetLineWidth(3)
        looseIndicator = ROOT.TLine(looseCut[i],0,looseCut[i],1)
        looseIndicator.SetLineStyle(2)
        looseIndicator.SetLineColor(3)
        looseIndicator.SetLineWidth(3)

        histList = []
        #h_dummy = ROOT.TH1D("dummy", "dummy", binning[0], binning[1], binning[2])
        #histList.append(h_dummy)

        ### Make histo ###
        h_comp = makeTH1(filename, gentree, pileup, binning[i], plotvar, "muon_isME0Muon")
        scale_comp = 1 / (h_comp.Integral())
        h_comp.Scale(scale_comp)
        h_comp.GetXaxis().SetTitle("%s"%plotvar)
        h_comp.GetXaxis().SetTitleSize(0.05)
        h_comp.GetYaxis().SetTitle(pileup + " Normalized")
        h_comp.SetLineColor(4)
        h_comp.SetLineWidth(3)
        histList.append(h_comp)

        h_fake = makeTH1(filename, recotree, pileup, binning[i], plotvar, "muon_isME0Muon && !muon_signal")
        scale_fake = 1 / (h_fake.Integral())
        h_fake.Scale(scale_fake)
        h_fake.SetLineColor(65)
        h_fake.SetLineWidth(3)
        histList.append(h_fake)

        for j, h in enumerate(histList):
            if h == h_comp:
                l.AddEntry(h, "True", "f")
            if h == h_fake:
                l.AddEntry(h, "Fake", "f")
            h.Draw("histsame")
        
        l.SetTextSize(0.033)
        l.Draw()
        tightIndicator.Draw()
        looseIndicator.Draw()

        ### CMS_lumi setting ###
        iPos = 0
        iPeriod = 0
        if (iPos == 0):
            CMS_lumi.relPosX = 0.12
        CMS_lumi.extraText = "Work in progress"
        CMS_lumi.lumi_sqrtS = "14TeV"
        CMS_lumi.CMS_lumi(c, iPeriod, iPos)
        #drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")

        c.Modified()
        c.Update()
        c.SaveAs("%s_%s.png" % (pileup, plotvar))
