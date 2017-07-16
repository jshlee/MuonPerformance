#!/usr/bin/env python
import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(1.0)
    h.SetLineColor(color)
    h.SetLineWidth(2)

#filenames = ["../../ZMM_PU0_910_pre3_relval/ZMM_PU0_910_pre3_relval.root" , "../../ZMM_PU0_910_pre3_me0seed/ZMM_PU0_910_pre3_me0seed.root"]
filenames = ["../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root" , "../../ZMM_PU200_910_pre3_me0seed/ZMM_PU200_910_pre3_me0seed.root"]

sample = ["RelVal", "ME0 seed"]
gentree = "MuonAnalyser/gen"
binning = [10, 1.4, 3.4]

CutValues = ["GE11", "GE21"]
plotvar = "fabs(muon.Eta())"

for cutval in CutValues:
    ### Canvas ###
    canvasname = cutval+plotvar+"Eff"
    c = makeCanvas(canvasname, False)
    l = ROOT.TLegend(.75, .85, .95, .95)

    markers = [20, 21]
    if cutval == "GE11":
        colors = [4, 65]
    if cutval == "GE21":
        colors = [2, 95] 

    hEfflist = []
    
    for i, filename in enumerate(filenames):
        if "PU0" in filename:
            pileup = "pu0"
        if "PU200" in filename:
            pileup = "pu200"
        
        if cutval == "GE11":
            h_p = makeTH1(filename, gentree, "passed", binning, plotvar, "muon_isGE11Muon")
        elif cutval == "GE21": 
            h_p = makeTH1(filename, gentree, "passed", binning, plotvar, "muon_isGE21Muon")
        h_t = makeTH1(filename, gentree, "total", binning, plotvar, "")

        h_eff = ROOT.TEfficiency(h_p, h_t)
        h_eff.SetLineColor(colors[i])
        hEfflist.append(h_eff)

    h_init = ROOT.TH1F("","",10, 1.3, 3.5)
    h_init.SetMaximum(1.05)
    h_init.SetMinimum(0)
    h_init.GetXaxis().SetTitle("Muon |#eta|")
    h_init.GetYaxis().SetTitle("%s Efficiency"%cutval)

    h_init.Draw("l")
    for i, h in enumerate(hEfflist):
        l.AddEntry(h, sample[i], "lp")
        h.Draw("p0lsame")
        setMarkerStyle(h, colors[i], markers[i])
    l.SetTextSize(0.02)
    l.SetBorderSize(0)
    l.Draw()

    ### For Eff ###
    iPos = 0
    iPeriod = 0
    if (iPos == 0):
        CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    c.Modified()
    c.Update()
    c.SaveAs("%s_910pre3_Efficiency_comparison_%s.png" % (pileup,cutval))

