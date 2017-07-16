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

filenames = ["../../pu200/pu200"]
binning = [10, -0.5, 9.5]

for i, filename in enumerate(filenames):
    treename = "MuonAnalyser/reco"

    if "pu0" in filename:
        pileup = "pu0"
    if "pu200" in filename:
        pileup = "pu200"

    h_list = []

    ### Canvas ###
    canvasname = filename
    c = makeCanvas(canvasname, False)
    #c.SetLogy()

    ### Make histo ###
    plotvar = "recoMuon_noRecHitME0"
    l = ROOT.TLegend(.75, .80, .95, .95)

    h_true = makeTH1(filename+".root", treename, pileup, binning, plotvar, 
                "recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && recoMuon_signal")
    scale_true = 1 / (h_true.Integral())
    h_true.Scale(scale_true)
    h_true.GetXaxis().SetTitle("ME0 #mu hits")
    h_true.GetYaxis().SetTitle(pileup + " Normalized")
    h_true.SetLineColor(4)
    h_list.append(h_true)

    h_fake = makeTH1(filename+".root", treename, pileup, binning, plotvar, 
                "recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && !recoMuon_signal")
    scale_fake = 1 / (h_fake.Integral())
    h_fake.Scale(scale_fake)
    h_fake.SetLineColor(3)
    h_list.append(h_fake)

    for i, h in enumerate(h_list):
        if h == h_true:
            l.AddEntry(h, "True", "l")
        if h == h_fake:
            l.AddEntry(h, "Fake", "l")
        h.Draw("histsame")
    l.Draw()

    ### CMS_lumi setting ###
    iPos = 0
    iPeriod = 0
    if (iPos == 0):
        CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Work in progress"
    CMS_lumi.lumi_sqrtS = "13TeV"
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)
    #drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")

    c.Modified()
    c.Update()
    c.SaveAs("%s_matchByHitME0mu.png" % pileup)
