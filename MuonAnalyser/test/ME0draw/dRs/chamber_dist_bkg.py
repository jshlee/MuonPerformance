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


filenames = ["../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root", "../../ZMM_PU200_910_pre3_me0seed/ZMM_PU200_910_pre3_me0seed.root"]
#filename = "../../out.root"
treename = "MuonAnalyser/reco"
binning = [ [31,0,0.03], [31,0,0.3], [36,0,3.5],
            [31,0,0.03], [31,0,0.3], [36,0,3.5],
            [31,0,0.03], [31,0,0.3], [36,0,3.5] ]

var = ["fabs(muon_ME0dPhi)", "fabs(muon_ME0dEta)", "fabs(muon_ME0pullPhi)", 
       "fabs(muon_GE11dPhi)", "fabs(muon_GE11dEta)", "fabs(muon_GE11pullPhi)",
       "fabs(muon_GE21dPhi)", "fabs(muon_GE21dEta)", "fabs(muon_GE21pullPhi)" ]

for filename in filenames:
    if "relval" in filename:
        sample = "relval"
    if "me0seed" in filename:
        sample = "me0seed"

    for j, plotvar in enumerate(var):

        pileup = "PU200"   
        ### Canvas ###
        canvasname = sample+plotvar
        c = makeCanvas(canvasname, False)
        c.SetLogy()

        ### Make histo ###
        if "ME0" in plotvar:
            h = makeTH1(filename, treename, pileup+plotvar, binning[j], plotvar,"muon_isME0Muon && !muon_signal")
        if "GE11" in plotvar:
            h = makeTH1(filename, treename, pileup+plotvar, binning[j], plotvar,"muon_isGE11Muon && !muon_signal")
        if "GE21" in plotvar:
            h = makeTH1(filename, treename, pileup+plotvar, binning[j], plotvar,"muon_isGE21Muon && !muon_signal")
        
        h.Scale(1/h.Integral())
        h.SetTitle("Stats")
        h.GetXaxis().SetTitle("%s"%plotvar)
        h.GetXaxis().SetTitleSize(0.05)
        h.GetYaxis().SetTitle("Normalized")
        h.SetMaximum(1)
        h.SetMinimum(0.001)
        h.SetLineWidth(3)
        h.SetLineColor(3)
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
        c.SaveAs("%s_%s_background_%s.png" % (pileup, plotvar, sample))
