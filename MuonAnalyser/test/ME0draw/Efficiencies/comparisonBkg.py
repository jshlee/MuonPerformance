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

filenames = ["../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root" , "../../ZMM_PU200_910_pre3_me0seed/ZMM_PU200_910_pre3_me0seed.root"]
sample = ["RelVal", "ME0 seed"]
recotree = "MuonAnalyser/reco"
binning = [10, 1.4, 3.4]

CutValues = ["gem", "me0"]
plotvar = "fabs(muon.Eta())"

for cutval in CutValues:
    ### Canvas ###
    canvasname = cutval+plotvar+"Bkg"
    c = makeCanvas(canvasname, False)
    l = ROOT.TLegend(.75, .85, .95, .95)

    markers = [20, 21]
    if cutval == "gem":
        colors = [4, 65]
    if cutval == "me0":
        colors = [2, 95] 

    hBkglist = []
    
    for i, filename in enumerate(filenames):
        if "PU0" in filename:
            pileup = "pu0"
        if "PU200" in filename:
            pileup = "pu200"

        tfile = ROOT.TFile(filename)
        nevents = tfile.Get("MuonAnalyser/nevents").Integral()
        if cutval == "gem":
            h_br = makeTH1(filename, recotree, "bkgrate", binning, plotvar, "!muon_signal && muon_isGEMMuon")
        elif cutval == "me0": 
            h_br = makeTH1(filename, recotree, "bkgrate", binning, plotvar, "!muon_signal && muon_isME0Muon")

        h_br.Scale(1/nevents)
        h_br.GetXaxis().SetTitle("Muon |#eta|")
        h_br.GetXaxis().SetRangeUser(1.3, 3.5)
        h_br.GetYaxis().SetTitle("Background/event")
        h_br.SetLineColor(colors[i])
        hBkglist.append(h_br)

    for i, h in enumerate(hBkglist):
        l.AddEntry(h, sample[i], "lp")
        h.Draw("p0lsame")
        setMarkerStyle(h, colors[i], markers[i])
    l.SetTextSize(0.02)
    #l.SetBorderSize(0)
    l.Draw()

    ### For Bkg ###
    iPos = 0
    iPeriod = 0
    if (iPos == 0):
        CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    c.Modified()
    c.Update()
    c.SaveAs("%s_910pre3_Backgroundrate_comparison_%s.png" % (pileup,cutval))
