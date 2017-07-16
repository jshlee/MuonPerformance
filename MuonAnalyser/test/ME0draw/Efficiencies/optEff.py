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
treename = "MuonAnalyser/gen"
binning = [8, 1.4, 2.8]
#binning = [60, -3.4, 3.4]

if "PU0" in filename:
    pileup = "pu0"
elif "200" in filename:
    pileup = "PU200"
elif "PU140" in filename:
    pileup = "pu140"

if "me0seed" in filename:
    relval = "me0seed"
else:
    relval = "relVal"
"""
if "pre4" in filename:
    version = "pre4"
elif "pre5" in filename:
    version = "pre5"
"""
version = "910_pre3"

plotvar = "fabs(muon.Eta())"
#CutValues = ["deltaX", "deltaY", "pullX", "pullY", "deltaDXDZ", "deltaDYDZ", "dPhi", "noRecHit", "combined", "default"]
CutValues = ["latest","new"]

for cutval in CutValues:
    ### Canvas ###
    canvasname = cutval+plotvar
    c = makeCanvas(canvasname, False)
    l = ROOT.TLegend(.75, .75, .95, .95)

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

    hlist = []
    h_background = makeTH1(filename, treename, "", binning, plotvar, "")

    if cutval == "latest":
        h_passedLoose = makeTH1(filename, treename, "", binning, plotvar, "muon_isGE11Muon && fabs(muon_GE11deltaX)<2.3 && fabs(muon_GE11deltaY)<19 && fabs(muon_GE11deltaDXDZ)<0.28 && fabs(muon_GE11deltaDYDZ)<0.43 && fabs(muon_GE11pullX)<2.3 && fabs(muon_GE11pullY)<7") # Loose
        h_effLoose = ROOT.TEfficiency(h_passedLoose, h_background)
        hlist.append(h_effLoose)

        h_passedTight = makeTH1(filename, treename, "", binning, plotvar, "muon_isGE11Muon && fabs(muon_GE11deltaX)<1 && fabs(muon_GE11deltaY)<8 && fabs(muon_GE11deltaDXDZ)<0.15 && fabs(muon_GE11deltaDYDZ)<0.42 && fabs(muon_GE11pullX)<1.2 && fabs(muon_GE11pullY)<3.5") # Tight
        h_effTight = ROOT.TEfficiency(h_passedTight, h_background)
        hlist.append(h_effTight)

    if cutval == "new":
        h_passedLoose = makeTH1(filename, treename, "", binning, plotvar, "muon_isGE11Muon && fabs(muon_GE11dPhi)<0.016 && fabs(muon_GE11dEta)<0.12") ## Loose
        h_effLoose = ROOT.TEfficiency(h_passedLoose, h_background)
        hlist.append(h_effLoose)

        h_passedTight = makeTH1(filename, treename, "", binning, plotvar, "muon_isGE11Muon && fabs(muon_GE11dPhi)<0.006 && fabs(muon_GE11dEta)<0.04") ## Tight
        h_effTight = ROOT.TEfficiency(h_passedTight, h_background)
        hlist.append(h_effTight)

    """
    if cutval == "default":
        h_defPassedTight = makeTH1(filename, treename, "", binning, plotvar, "muon_isME0Muon && (muon_ME0deltaX<3 || muon_ME0pullX<3) && (muon_ME0deltaY<3 || muon_ME0pullY<3) && muon_ME0dPhi<0.1") ## default Tight
        h_defEffTight = ROOT.TEfficiency(h_defPassedTight, h_background)
        hlist.append(h_defEffTight)

        h_defPassedLoose = makeTH1(filename, treename, "", binning, plotvar, "muon_isME0Muon && (muon_ME0deltaX<3 || muon_ME0pullX<3) &&  (muon_ME0deltaY<3 || muon_ME0pullY<3) && muon_ME0dPhi<0.5") ## default Loose
        h_defEffLoose = ROOT.TEfficiency(h_defPassedLoose, h_background)
        hlist.append(h_defEffLoose)
    """
    """
    else:
        for i, cut in enumerate(cuts):
            h_passed = makeTH1(filename, treename, cut, binning, plotvar, "muon_isME0Muon && %s %s" %(pcut, cut))
            h_eff = ROOT.TEfficiency(h_passed, h_background)
            if i == 0:
                h_eff.SetLineColor(2)
            if i == 1:
                h_eff.SetLineColor(3)
            hlist.append(h_eff)
    """ 
    h_me0 = makeTH1(filename, treename, "", binning, plotvar, "muon_isGE11Muon")
    h_me0Eff = ROOT.TEfficiency(h_me0, h_background)
    h_me0Eff.SetLineColor(4)
    hlist.append(h_me0Eff)

    h_init = ROOT.TH1F("", "", 8, 1.4, 2.8)
    h_init.SetMaximum(1.05)
    h_init.SetMinimum(0)
    h_init.GetXaxis().SetTitle("Muon |#eta|")
    h_init.GetYaxis().SetTitle("Efficiency")
    h_init.Draw("l")

    for i, h in enumerate(hlist):
        if i == 0:
            l.AddEntry(h, "Loose", "lp")
            h.Draw("p0lsame")
            setMarkerStyle(h, 3, 22)
        if i == 1:
            l.AddEntry(h, "Tight", "lp")
            h.Draw("p0lsame")
            setMarkerStyle(h, 2, 20)
        if i == 2:
            l.AddEntry(h, "isGEMMuon", "lp")
            h.Draw("p0lsame")
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
    c.SaveAs("GE11_%s_OptEff_%s_%s_%s.png" % (pileup,cutval,version,relval))
