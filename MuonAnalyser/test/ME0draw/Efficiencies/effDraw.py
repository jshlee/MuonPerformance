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

filename = "../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root"
tfile = ROOT.TFile(filename)
nevents = tfile.Get("MuonAnalyser/nevents").Integral()
gentree = "MuonAnalyser/gen"
recotree = "MuonAnalyser/reco"
binning = [10, 1.4, 3.4] #Eta range
plotvar = "fabs(muon.Eta())"

pileup = "PU200"


### Dictionaries ###
isME0Muon  = { "isMuon": "muon_isME0Muon" }
isME0MuonSelNew  = { "isMuon": "muon_isME0MuonSelNew" }

usrdic = [isME0Muon, isME0MuonSelNew]


for draw in ["Eff", "Bkg"]:
    ### Canvas ###
    canvasname = draw+plotvar
    c = makeCanvas(canvasname, False)
    l = ROOT.TLegend(.75, .75, .95, .95)

    hlist = []
    if draw == "Eff":
        h_Effbkg = makeTH1(filename, gentree, "Effbkg", binning, plotvar, "")
        for i,opt in enumerate(usrdic):
            h_eff = makeTH1(filename, gentree, "EffnoCut", binning, plotvar, opt["isMuon"])
            h_teff = ROOT.TEfficiency(h_eff, h_Effbkg)
            hlist.append(h_teff)

        h_init = ROOT.TH1F("", "", 10, 1.4, 3.4)
        h_init.SetMaximum(1.05)
        h_init.SetMinimum(0)
        h_init.GetXaxis().SetTitle("Muon |#eta|")
        h_init.GetYaxis().SetTitle("Efficiency")
        h_init.Draw("l")

    if draw == "Bkg":
        c.SetLogy()
        for i,opt in enumerate(usrdic):
            h_noCut = makeTH1(filename, recotree, "BkgNoCut", binning, plotvar, "!muon_signal && " + opt["isMuon"])
            h_noCut.Scale(1/nevents)
            hlist.append(h_noCut)

    for i, h in enumerate(hlist):
        if i == 0:
            if draw == "Eff":
                l.AddEntry(h, "isME0Muon", "lp")
                h.Draw("p0lsame")
                setMarkerStyle(h, 4, 23)
            if draw == "Bkg":
                l.AddEntry(h, "isME0Muon", "f")
                h.SetFillColor(4)
                h.SetFillStyle(3003)
                h.GetXaxis().SetTitle("Muon |#eta|")
                h.GetYaxis().SetTitle("Background/event")
                h.Draw("histsame")
        if i == 1:
            if draw == "Eff":
                l.AddEntry(h, "isME0MuonSelNew", "lp")
                h.Draw("p0lsame")
                setMarkerStyle(h, 3, 22)
            if draw == "Bkg":
                l.AddEntry(h, "isME0MuonSelNew", "f")
                h.SetFillColor(3)
                h.SetFillStyle(3005)
                h.Draw("histsame")
    l.SetTextSize(0.02)
    l.Draw()

    ### CMS_lumi setting ###
    iPos = 0
    iPeriod = 0
    if (iPos == 0):
        CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.lumi_sqrtS = "13TeV"
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    c.Modified()
    c.Update()
    c.SaveAs("OptPlots/%s_ME0_Opt%s_910pre3_SelNew.png" % (pileup,draw))
