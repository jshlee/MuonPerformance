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

filenames = ["../../ZMM_PU200_pre5/ZMM_PU200_pre5.root" ]
#filenames = ["../../ZMM_PU200_pre5/ZMM_PU200_pre5.root" ,"../../ZMM_PU200_pre5_me0seed/ZMM_PU200_pre5_me0seed.root"]
treename = "MuonAnalyser/reco"
binning = [10, 1.4, 3.4]

for filename in filenames:
    tfile = ROOT.TFile(filename)
    nevents = tfile.Get("MuonAnalyser/nevents").Integral()

    if "PU0" in filename:
        pileup = "pu0"
    elif "PU200" in filename:
        pileup = "pu200"
    elif "PU140" in filename:
        pileup = "pu140"

    if "me0seed" in filename:
        sampleType = "me0seed"
    else:
        sampleType = "RelVal"

    if "pre4" in filename:
        version = "pre4"
    elif "pre5" in filename:
        version = "pre5"

    plotvar = "abs(muon.Eta())"
    CutValues = ["me0"]
    #CutValues = ["dx", "dy", "dDx", "dDy", "hit", "me0"]
    for cutval in CutValues:
        ### Canvas ###
        canvasname = filename+cutval+plotvar
        c = makeCanvas(canvasname, False)
        l = ROOT.TLegend(.75, .80, .95, .95)

        hlist = []
        defaultCut = "!muon_signal && muon_isME0Muon"
        if cutval == "hit":
            cuts = ["2","3","4","5","6","7","8"]
            colors = [91, 81, 71, 61, 51, 41, 31]
            markers = [20, 33, 34, 21, 22, 29, 5]
            pcut = defaultCut + " && muon_noRecHitME0>"
            legentries = "ME0Hit > "

        if cutval == "dx":
            cuts = ["4.0", "3.0", "2.0", "1.0", "0.5"]
            colors = [1, 8, 30, 4, 6]
            markers = [20, 33, 34, 21, 22]
            pcut = defaultCut + " && muon_deltaXME0<"
            legentries = "No #Deltax/#sigma_{x}, #Deltax < "

        if cutval == "dy":
            cuts = ["20.0", "10.0", "5.0", "4.0", "3.0", "2.0", "1.0"]
            colors = [ 6, 2, 93, 41, 8, 30, 4]
            markers = [23, 24, 25, 26, 33, 34, 21]
            pcut = defaultCut + " && muon_deltaYME0<"
            legentries = "No #Deltax/#sigma_{x}, #Deltay < "

        if cutval == "dDx":
            cuts = ["0.6","0.5","0.4","0.3","0.2","0.1"]
            colors = [100, 90, 80, 70, 60, 50]
            markers = [20, 21, 22, 27, 28, 29]
            pcut = defaultCut + " && muon_deltaDXDZME0<"
            legentries = "No #Deltax/#sigma_{x}, #Deltadx < "

        if cutval == "dDy":
            cuts = ["1.4","1.2","1.0","0.8","0.6","0.4","0.2"]
            colors = [94, 84, 74, 64, 54, 44, 34]
            markers = [20, 21, 22, 27, 28, 29, 34]
            pcut = defaultCut + " && muon_deltaDYDZME0<"
            legentries = "No #Deltax/#sigma_{x}, #Deltady < "

        if cutval == "me0":
            cuts = [""]
            colors = [63]
            markers = [20]
            pcut = defaultCut
            legentries = "isME0Muon"

        if cutval == "me0gem":
            cuts = ["muon_isGEMMuon", "muon_isME0Muon"]
            colors = [2, 4]
            markers = [20, 21]
            pcut = "muon.Pt()>5 && "
            legentries = ""

        for i, cut in enumerate(cuts):
            h_bkg = makeTH1(filename, treename, cut, binning, plotvar, pcut+"%s" % cut)
            h_bkg.Scale(1/nevents)
            h_bkg.GetXaxis().SetTitle("Muon |#eta|")
            h_bkg.GetXaxis().SetRangeUser(1.3, 3.5)
            h_bkg.GetYaxis().SetTitle("Background/event")
            h_bkg.SetFillColor(colors[i])
            hlist.append(h_bkg)

        for i, h in enumerate(hlist):
            l.AddEntry(h, legentries+"%s" % cuts[i], "f")
            h.Draw("histsame")
            #setMarkerStyle(h, colors[i], markers[i])
        l.SetTextSize(0.025)
        l.SetBorderSize(0)
        if cutval != "me0":
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
        c.SaveAs("%s_Backgroundrate_%s_%s_%s.png" % (pileup,cutval,sampleType,version))
