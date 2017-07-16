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
binning = [ [41,0,4],      # delta x 
            [21,0,20],     # delta y
            [81,0,0.8],    # delta dxdz
            [21,0,2],      # delta dydz
            [25,0,2.5],    # delta phi
            [32,0,3.2]]    # delta Eta

var = [ "deltaX",
        "deltaY",
        "deltaDXDZ",
        "deltaDYDZ",
        "dPhi",
        "dEta" ]

for filename in filenames:
    if "PU0" in filename:
        pileup = "pu0"
    elif "PU200" in filename:
        pileup = "pu200"
    version = "910pre3"

    for i,plotvar in enumerate(var):
        ### Canvas ###
        canvasname = pileup+plotvar+version+filename
        c = makeCanvas(canvasname, False)
        c.SetLogy()
        l = ROOT.TLegend(.70, .70, .95, .95)
    
        h_GE11true = makeTH1(filename, treename, pileup, binning[i], "fabs(muon_GE11%s)"%plotvar, "muon.Pt()>5 && muon_isGEMMuon && muon_signal")
        h_GE11fake = makeTH1(filename, treename, pileup, binning[i], "fabs(muon_GE11%s)"%plotvar, "muon.Pt()>5 && muon_isGEMMuon && !muon_signal")
        h_GE21true = makeTH1(filename, treename, pileup, binning[i], "fabs(muon_GE21%s)"%plotvar, "muon.Pt()>5 && muon_isGEMMuon && muon_signal")
        h_GE21fake = makeTH1(filename, treename, pileup, binning[i], "fabs(muon_GE21%s)"%plotvar, "muon.Pt()>5 && muon_isGEMMuon && !muon_signal")
        h_ME0true  = makeTH1(filename, treename, pileup, binning[i], "fabs(muon_ME0%s)"%plotvar, "muon.Pt()>5 && muon_isME0Muon && muon_signal")
        h_ME0fake  = makeTH1(filename, treename, pileup, binning[i], "fabs(muon_ME0%s)"%plotvar, "muon.Pt()>5 && muon_isME0Muon && !muon_signal")
        histList = [ h_GE11true, h_GE11fake, h_GE21true, h_GE21fake, h_ME0true, h_ME0fake ]

        h_init = ROOT.TH1F( "", "", binning[i][0], binning[i][1], binning[i][2] )
        h_init.SetMaximum(1)
        h_init.GetXaxis().SetTitle("%s" % plotvar)
        h_init.GetYaxis().SetTitle("PU0 Normalised")
        h_init.Draw("l")

        for h in histList:
            h.Scale( 1/h.Integral() )
            if h == h_GE11true:
                h.SetLineColor(2) 
                h.SetMarkerStyle(20)
                h.SetMarkerColor(2)
                l.AddEntry(h, "True: GE11", "pl")
            if h == h_GE11fake:
                h.SetLineColor(95)
                h.SetMarkerStyle(24)
                h.SetMarkerColor(95)
                l.AddEntry(h, "Fake: GE11", "pl")
            if h == h_GE21true:
                h.SetLineColor(3) 
                h.SetMarkerStyle(21)
                h.SetMarkerColor(3)
                l.AddEntry(h, "True: GE21", "pl")
            if h == h_GE21fake:
                h.SetLineColor(87)
                h.SetMarkerStyle(25)
                h.SetMarkerColor(87)
                l.AddEntry(h, "Fake: GE21", "pl")
            if h == h_ME0true:
                h.SetLineColor(4)
                h.SetMarkerStyle(34)
                h.SetMarkerColor(4)
                l.AddEntry(h, "True: ME0", "pl")
            if h == h_ME0fake:
                h.SetLineColor(64)
                h.SetMarkerStyle(28)
                h.SetMarkerColor(64)
                l.AddEntry(h, "Fake: ME0", "pl")
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
        c.SaveAs("GEM-vs-ME0_%s_%s_%s.png" % (pileup, plotvar, version))
