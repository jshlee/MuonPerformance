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

filenames = ["../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root","../../ZMM_PU200_910_pre3_me0seed/ZMM_PU200_910_pre3_me0seed.root"]
treename = "MuonAnalyser/reco"
"""
binning = [ [41,0,4],      # delta x 
            [21,0,20],     # delta y
            [41,0,4],      # pull x
            [21,0,20],     # pull y
            [21,0,1.0],    # delta phi
            [21,0,1.0] ]   # delta eta

var = [ "deltaX",
        "deltaY",
        "pullX",
        "pullY",
        "dPhi",
        "dEta" ]
"""

var = [ "deltaX", "deltaY", "pullX", "pullY", "deltaDXDZ", "deltaDYDZ", "dPhi", "dEta", "pullPhi" ]

binning = [[41,0,40], [41,0,40], [41,0,40], [41,0,40], [21,0,1], [21,0,1], [21,0,1], [21,0,1], [36,0,3.5] ]

chambers = ["GE11", "GE21", "ME0"]

for chamber in chambers:
    for i,plotvar in enumerate(var):
        hlist = []
        ### Canvas ###
        canvasname = plotvar + chamber
        c = makeCanvas(canvasname, False)
        c.SetLogy()
        l = ROOT.TLegend(.70, .70, .95, .95)

        for filename in filenames:

            if "PU0" in filename:
                pileup = "pu0"
            elif "PU200" in filename:
                pileup = "pu200"

            version = "910pre3"

            """
            if "me0seed" in filename:
                sampleType = "me0seed"
            if "relval" in filename:
                sampleType = "relval"
            """
            if "GE" in chamber:
                h_true = makeTH1(filename, treename, filename+chamber, binning[i], "muon_"+chamber+plotvar, "muon_isGEMMuon && muon_signal")
                h_fake = makeTH1(filename, treename, filename+chamber, binning[i], "muon_"+chamber+plotvar, "muon_isGEMMuon && !muon_signal")
            else:
                h_true = makeTH1(filename, treename, filename+chamber, binning[i], "muon_"+chamber+plotvar, "muon_isME0Muon && muon_signal")
                h_fake = makeTH1(filename, treename, filename+chamber, binning[i], "muon_"+chamber+plotvar, "muon_isME0Muon && !muon_signal")
            hlist.append(h_true)
            hlist.append(h_fake)

        h_init = ROOT.TH1F( "", "", binning[i][0], binning[i][1], binning[i][2] )
        h_init.SetMaximum(1)
        h_init.GetXaxis().SetTitle("%s" % (chamber+plotvar))
        h_init.GetYaxis().SetTitle("%s Normalised"%pileup)
        h_init.Draw("l")

        for j, h in enumerate(hlist):
            h.Scale( 1/h.Integral() )
            if j == 0:  #relval true
                h.SetLineColor(2) 
                h.SetMarkerStyle(20)
                h.SetMarkerColor(2)
                l.AddEntry(h, "RelVal, True: %s"%chamber, "pl")
            if j == 1:  #relval fake
                h.SetLineColor(95)
                h.SetMarkerStyle(24)
                h.SetMarkerColor(95)
                l.AddEntry(h, "RelVal, Fake: %s"%chamber, "pl")
            if j == 2:  #me0seed true
                h.SetLineColor(4) 
                h.SetMarkerStyle(22)
                h.SetMarkerColor(4)
                l.AddEntry(h, "me0seed, True: %s"%chamber, "pl")
            if j == 3:  #me0seed fake
                h.SetLineColor(65)
                h.SetMarkerStyle(26)
                h.SetMarkerColor(65)
                l.AddEntry(h, "me0seed, Fake: %s"%chamber, "pl")
            h.Draw("p0lsame")

        l.SetTextSize(0.027)
        l.Draw()

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
        c.SaveAs("170613/%s_%s_%s.png" % (pileup, version, chamber+plotvar))
