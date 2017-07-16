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

filenames = ["../../ZMM_PU200_pre5/ZMM_PU200_pre5.root","../../ZMM_PU200_pre5_me0seed/ZMM_PU200_pre5_me0seed.root"]
treename = "MuonAnalyser/reco"
binning = [ [41,0,4],      # delta x 
            [21,0,20],     # delta y
            [41,0,4],    # pull x
            [21,0,20],     # pull y
            [21,0,1.0]]    # delta phi

var = [ "deltaX",
        "deltaY",
        "pullX",
        "pullY",
        "dPhi" ]
chambers = ["GE11", "GE21", "ME0"]

for chamber in chambers:
    for i,plotvar in enumerate(var):
        hlist = []
        ### Canvas ###
        canvasname = plotvar + chamber
        c = makeCanvas(canvasname, False)
        #c.SetLogy()
        l = ROOT.TLegend(.2, .75, .45, .95)

        for filename in filenames:
            if "PU0" in filename:
                pileup = "pu0"
            elif "PU200" in filename:
                pileup = "pu200"
            if "pre4" in filename:
                version = "pre4"
            elif "pre5" in filename:
                version = "pre5"
            if "me0seed" in filename:
                sampleType = "me0seed"
            else:
                sampleType = "RelVal"

            h_true = makeTH1(filename, treename, filename+chamber, binning[i], "muon_"+chamber+plotvar, "muon.Pt()>5 && muon_isGEMMuon && muon_signal")
            h_trueTemp = ROOT.TH1F("", "", binning[i][0], binning[i][1], binning[i][2] )
            contents = 0
            for j in range(binning[i][0]):
                contents += h_true.GetBinContent(j) 
                h_trueTemp.SetBinContent( j, contents )
            hlist.append(h_trueTemp)

            h_fake = makeTH1(filename, treename, filename+chamber, binning[i], "muon_"+chamber+plotvar, "muon.Pt()>5 && muon_isGEMMuon && !muon_signal")
            h_fakeTemp = ROOT.TH1F("", "", binning[i][0], binning[i][1], binning[i][2] )
            contents = 0
            for j in range(binning[i][0]):
                contents += h_fake.GetBinContent(j) 
                h_fakeTemp.SetBinContent( j, contents )
            hlist.append(h_fakeTemp)

        h_init = ROOT.TH1F( "", "", binning[i][0], binning[i][1], binning[i][2] )
        h_init.SetMaximum(.1)
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
        c.SaveAs("%s_%s_%s_%s.png" % (pileup, chamber+plotvar, version, sampleType))
    
        print "%s_%s_%s_%s.png complete..." % (pileup, chamber+plotvar, version, sampleType)
