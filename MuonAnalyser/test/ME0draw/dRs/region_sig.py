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


filename = "../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root"
#filename = "../../out.root"
treename = "MuonAnalyser/reco"
binning = [ [21,0,5], 
            [21,0,0.02],
            [21,0,0.2],
            [21,0,20],
            [36,0,3.5] ]

var = ["fabs(muon_ME0deltaX)", "fabs(muon_ME0dPhi)", "fabs(muon_ME0dEta)", "fabs(muon_ME0deltaY)", "fabs(muon_ME0pullPhi)"]
regions = ["fabs(muon_ME0chamX)<5", "fabs(muon_ME0chamX)>20"]

for j, plotvar in enumerate(var):

    pileup = "PU200"   

    for i,region in enumerate(regions):
        ### Canvas ###
        canvasname = region+plotvar
        c = makeCanvas(canvasname, False)
        c.SetLogy()

        ### Make histo ###
        if "seg" and ">20" in region:
            h = makeTH1(filename, treename, pileup+plotvar, binning[j], plotvar,"muon_isME0Muon && muon_signal && muon_ME0segX!=100 && %s"%region)
        if "seg" and "<5" in region:
            h = makeTH1(filename, treename, pileup+plotvar, binning[j], plotvar,"muon_isME0Muon && muon_signal && %s"%region)
        if "cham" and ">20" in region:
            h = makeTH1(filename, treename, pileup+plotvar, binning[j], plotvar,"muon_isME0Muon && muon_signal && muon_ME0chamX!=100 && %s"%region)
        if "cham" and "<5" in region:
            h = makeTH1(filename, treename, pileup+plotvar, binning[j], plotvar,"muon_isME0Muon && muon_signal && %s"%region)
        h.Scale(1/h.Integral())
        h.SetTitle("Stats")
        h.GetXaxis().SetTitle("%s"%plotvar)
        h.GetXaxis().SetTitleSize(0.05)
        h.GetYaxis().SetTitle("Normalized")
        h.SetMaximum(1)
        h.SetMinimum(0.001)
        h.SetLineWidth(3)
        if i == 0 or i==2:
            h.SetLineColor(2)
        if i == 1 or i==3:
            h.SetLineColor(4)
        ROOT.gStyle.SetOptStat()
        h.Draw("hist")

        ### CMS_lumi setting ###
        iPos = 0
        iPeriod = 0
        if (iPos == 0):
            CMS_lumi.relPosX = 0.12
        CMS_lumi.extraText = "Work in progress"
        CMS_lumi.lumi_sqrtS = ""
        CMS_lumi.CMS_lumi(c, iPeriod, iPos)
        #drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")

        c.Modified()
        c.Update()
        c.SaveAs("%s_%s_signal_%s.png" % (pileup, plotvar,region))
