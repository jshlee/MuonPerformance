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

#filename = "../../ZMM_PU140_pre6/ZMM_PU140_pre6.root" 
#filename = "../../ZMM_PU0_pre6_test/ZMM_PU0_pre6_test.root" 
#filename = "../../ZMM_PU140_RelVal_me0seed/ZMM_PU140_RelVal_me0seed.root" 
#filename = "../../ZMM_PU200_pre5/ZMM_PU200_pre5.root" 
filename = "../../ZMM_PU200_910_pre3_me0seed/ZMM_PU200_910_pre3_me0seed.root"
treename = "MuonAnalyser/reco"
binning = [10, 1.4, 3.4]
#binning = [60, -3.4, 3.4]

if "PU0" in filename:
    pileup = "pu0"
elif "PU200" in filename:
    pileup = "PU200"
elif "PU140" in filename:
    pileup = "pu140"

if "me0seed" in filename:
    sampleType = "me0seed"
else:
    sampleType = "relVal"
"""
if "pre4" in filename:
    version = "pre4"
elif "pre5" in filename:
    version = "pre5"
"""
version = "910_pre3"

plotvar = "fabs(muon.Eta())"
#CutValues = ["deltaX", "deltaY", "pullX", "pullY", "deltaDXDZ", "deltaDYDZ", "dPhi", "noRecHit", "combined"]
CutValues = ["latest", "new"]
tfile = ROOT.TFile(filename)
nevents = tfile.Get("MuonAnalyser/nevents").Integral()

for cutval in CutValues:
    ### Canvas ###
    canvasname = cutval+plotvar
    c = makeCanvas(canvasname, False)
    c.SetLogy()
    l = ROOT.TLegend(.80, .77, .95, .95)

    if cutval == "deltaX":
        cuts = ["1", "2"]
        pcut = "muon_ME0deltaX < "
    if cutval == "deltaY":
        cuts = ["4", "8"]
        pcut = "muon_ME0deltaY < "
    if cutval == "pullX":
        cuts = ["1.5", "3"]
        pcut = "muon_ME0pullX < "
    if cutval == "pullY":
        cuts = ["3", "8"]
        pcut = "muon_ME0pullY < "
    if cutval == "deltaDXDZ":
        cuts = ["0.05", "0.075"]
        pcut = "muon_ME0deltaDXDZ < "
    if cutval == "deltaDYDZ":
        cuts = ["0.3", "0.4"]
        pcut = "muon_ME0deltaDYDZ < "
    if cutval == "dPhi":
        cuts = ["0.05", "0.1"]
        pcut = "muon_ME0dPhi < "
    if cutval == "noRecHit":
        cuts = ["3", "2"]
        pcut = "muon_ME0noRecHit > "
    else:
        cuts = ""
        pcut = ""

    hlist = []
    h_me0 = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon")
    #h_me0.SetMaximum(100)
    h_me0.GetXaxis().SetTitle("Muon |#eta|")
    h_me0.GetYaxis().SetTitle("Background/Events")
    h_me0.Scale(1/nevents)
    h_me0.SetFillColor(4)
    h_me0.SetFillStyle(3003)
    hlist.append(h_me0)

    #me0seed
    if sampleType == "me0seed":
        if cutval == "new":
            h_bkgLoose = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11dPhi)<0.016 && fabs(muon_GE11dEta)<0.12") ## Loose
            h_bkgLoose.Scale(1/nevents)
            hlist.append(h_bkgLoose)

            h_bkgTight = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11dPhi)<0.006 && fabs(muon_GE11dEta)<0.04") ## Tight
            h_bkgTight.Scale(1/nevents)
            hlist.append(h_bkgTight)

        if cutval == "latest":
            h_defBkgLoose = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11deltaX)<2.3 && fabs(muon_GE11deltaY)<19 && fabs(muon_GE11deltaDXDZ)<0.28 && fabs(muon_GE11deltaDYDZ)<0.43 && fabs(muon_GE11pullX)<2.3 && fabs(muon_GE11pullY)<7") # Loose
            h_defBkgLoose.Scale(1/nevents)
            hlist.append(h_defBkgLoose)

            h_defBkgTight = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11deltaX)<1 && fabs(muon_GE11deltaY)<8 && fabs(muon_GE11deltaDXDZ)<0.15 && fabs(muon_GE11deltaDYDZ)<0.42 && fabs(muon_GE11pullX)<1.2 && fabs(muon_GE11pullY)<3.5") # Tight
            h_defBkgTight.Scale(1/nevents)
            hlist.append(h_defBkgTight)

    #relval
    if sampleType == "relVal":
        if cutval == "new":
            h_bkgLoose = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11dPhi)<0.014 && fabs(muon_GE11dEta)<0.27") ## Loose
            h_bkgLoose.Scale(1/nevents)
            hlist.append(h_bkgLoose)

            h_bkgTight = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11dPhi)<0.006 && fabs(muon_GE11dEta)<0.1") ## Tight
            h_bkgTight.Scale(1/nevents)
            hlist.append(h_bkgTight)

        if cutval == "latest":
            h_defBkgLoose = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11deltaX)<2.5 && fabs(muon_GE11deltaY)<45 && fabs(muon_GE11deltaDXDZ)<0.5 && fabs(muon_GE11deltaDYDZ)<6 && fabs(muon_GE11pullX)<25 && fabs(muon_GE11pullY)<25") # Loose
            h_defBkgLoose.Scale(1/nevents)
            hlist.append(h_defBkgLoose)

            h_defBkgTight = makeTH1(filename, treename, "", binning, plotvar, "!muon_signal && muon_isGE11Muon && fabs(muon_GE11deltaX)<1.4 && fabs(muon_GE11deltaY)<17 && fabs(muon_GE11deltaDXDZ)<0.27 && fabs(muon_GE11deltaDYDZ)<1.8 && fabs(muon_GE11pullX)<1.6 && fabs(muon_GE11pullY)<10") # Tight
            h_defBkgTight.Scale(1/nevents)
            hlist.append(h_defBkgTight)

    for i, h in enumerate(hlist):
        if i == 1:
            l.AddEntry(h, "Loose", "f")
            h.SetFillColor(3)
            h.SetFillStyle(3005)
            h.Draw("histsame")
            #setMarkerStyle(h, 3, 22)
        if i == 2:
            l.AddEntry(h, "Tight", "f")
            h.SetFillColor(2)
            h.SetFillStyle(3004)
            h.Draw("histsame")
            #setMarkerStyle(h, 2, 20)
        if i == 0:
            l.AddEntry(h, "isGE11Muon", "f")
            h.Draw("histsame")
            setMarkerStyle(h, 4, 23)
    l.SetTextSize(0.02)
    #l.SetBorderSize(0)
    l.Draw()

    ### CMS_lumi setting ###
    iPos = 0
    iPeriod = 0
    if (iPos == 0):
        CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Simulation"
    #CMS_lumi.lumi_sqrtS = "13TeV"
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    c.Modified()
    c.Update()
    c.SaveAs("GE11_%s_OptBkg_%s_%s_%s.png" % (pileup,cutval,version,sampleType))
