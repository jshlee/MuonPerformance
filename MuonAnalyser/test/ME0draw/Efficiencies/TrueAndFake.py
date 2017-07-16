#!/usr/bin/env python
import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

filenames = ["../../pu0/QCD_pu0.root", 
             "../../pu0/TTbar_pu0.root", 
             "../../pu0/pu0.root", 
             "../../pu200/pu200.root"]
treename = "MuonAnalyser/reco"
binning = [[12,5,65], #pT
           [32,1.8,3.0], #eta
           [12,-3,3] #phi
           ]

for filename in filenames:
    if "QCD" in filename:
        pileup = "QCD_pu0"
    elif "TTbar" in filename:
        pileup = "TTbar_pu0"
    elif "pu200" in filename:
        pileup = "pu200"
    else:
        pileup = "pu0"

    for plotvar in (["recoMuon.Pt()", "abs(recoMuon.Eta())", "recoMuon.Phi()"]):
        if "Pt()" in plotvar:
            bin = binning[0]
        if "Eta()" in plotvar:
            bin = binning[1]
        if "Phi()" in plotvar:
            bin = binning[2]

        ### Canvas ###
        canvasname = ""+pileup+plotvar
        c = makeCanvas(canvasname, False)
        l = ROOT.TLegend(.85, .9, .99, .99)
        #c.SetLogy()

        hlist = []
        if "Eta()" in plotvar:
            h_true = makeTH1(filename, treename, pileup, bin, plotvar,
                             "recoMuon.Pt()>5 && recoMuon_isME0Muon && recoMuon_signal")
        else:
            h_true = makeTH1(filename, treename, pileup, bin, plotvar,
                             "recoMuon.Pt()>5 && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && recoMuon_isME0Muon && recoMuon_signal")
        scale_true = 1 / h_true.Integral()
        h_true.Scale(scale_true)
        h_true.GetXaxis().SetTitle(plotvar)
        h_true.GetYaxis().SetTitle(pileup + " Normalized")
        h_true.SetLineColor(4)
        hlist.append(h_true)

        if "Eta()" in plotvar:
            h_fake = makeTH1(filename, treename, pileup, bin, plotvar,
                             "recoMuon.Pt()>5 && recoMuon_isME0Muon && !recoMuon_signal")
        else:
            h_fake = makeTH1(filename, treename, pileup, bin, plotvar,
                             "recoMuon.Pt()>5 && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && recoMuon_isME0Muon && !recoMuon_signal")
            
        scale_fake = 1 / h_fake.Integral()
        h_fake.Scale(scale_fake)
        h_fake.SetLineColor(3)
        hlist.append(h_fake)

        for i, h in enumerate(hlist):
            if h == h_true:
                l.AddEntry(h, "True", "l")
            if h == h_fake:
                l.AddEntry(h, "Fake", "l")
            h.Draw("histsame")
        l.SetTextSize(0.027)
        l.Draw()

        ### CMS_lumi setting ###
        iPos = 0
        iPeriod = 0
        if (iPos == 0):
            CMS_lumi.relPosX = 0.12
        CMS_lumi.extraText = "work in progress"
        #CMS_lumi.lumi_sqrtS = "13TeV"
        CMS_lumi.CMS_lumi(c, iPeriod, iPos)

        c.Modified()
        c.Update()
        c.SaveAs("%s_%s.png" % (pileup, plotvar))
