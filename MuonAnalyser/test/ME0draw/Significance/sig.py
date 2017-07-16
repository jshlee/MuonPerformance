#!/usr/bin/env python

import ROOT, copy, math
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
from MuonPerformance.MuonAnalyser.histoHelper import *

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.3, 0.8, samplename)

filenames = ["../../ZMM_PU200_pre5/ZMM_PU200_pre5.root"]
gentree = "MuonAnalyser/gen"
recotree = "MuonAnalyser/reco"
binning = [ [41,0,4],      # delta x 
            [21,0,20],     # delta y
            [41,0,4],      # pull x
            [21,0,20],     # pull y
            [81,0,0.8],    # delta dxdz
            [21,0,2],      # delta dydz
            [21,0,1.0],    # delta phi
            [8,0,7] ]      # noRecHit
var = [ "deltaX",
        "deltaY",
        "pullX",
        "pullY",
        "deltaDXDZ",
        "deltaDYDZ",
        "dPhi",
        "noRecHit" ]

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

    for i,plotvar in enumerate(var):
        ### Canvas ###
        canvasname = pileup+plotvar+version+filename
        c = makeCanvas(canvasname, False)
        #c.SetLogy()
    
        binsize = binning[i][0]
        binmin = binning[i][1]
        binmax = binning[i][2]
        h_sig = ROOT.TH1F(filename+plotvar, filename+plotvar, binsize, binmin, binmax)
        for j in range( binsize ):
            cut = float(j)*binmax/(binsize-1.0)
            print cut
            if plotvar == "noRecHit":
                h_ME0sig  = makeTH1(filename, gentree, pileup, [10,1.4,3.4], "muon.Eta()", "muon.Pt()>5 && muon_isME0Muon && muon_ME0"+plotvar+">"+str(cut) )
                h_ME0bkg  = makeTH1(filename, recotree, pileup, [10,1.4,3.4], "muon.Eta()", "muon.Pt()>5 && muon_isME0Muon && !muon_signal && muon_ME0"+plotvar+">"+str(cut) )
            else:
                h_ME0sig  = makeTH1(filename, gentree, pileup, [10,1.4,3.4], "muon.Eta()", "muon.Pt()>5 && muon_isME0Muon && muon_ME0"+plotvar+"<"+str(cut) )
                h_ME0bkg  = makeTH1(filename, recotree, pileup, [10,1.4,3.4], "muon.Eta()", "muon.Pt()>5 && muon_isME0Muon && !muon_signal && muon_ME0"+plotvar+"<"+str(cut) )
            sig = h_ME0sig.Integral()
            bkg = h_ME0bkg.Integral()
            if sig+bkg == 0:
                h_sig.SetBinContent( j, 0 )
            else:
                h_sig.SetBinContent( j, sig/(sig+bkg) )

        h_sig.SetLineColor(2)
        h_sig.SetMarkerStyle(20)
        h_sig.SetMarkerColor(2)
        h_sig.GetXaxis().SetTitle("%s" % plotvar)
        h_sig.GetYaxis().SetTitle("%s "%pileup + "Significance")
        h_sig.Draw("p0 L")

        ### CMS_lumi setting ###
        iPos = 0
        iPeriod = 0
        if (iPos == 0):
            CMS_lumi.relPosX = 0.12
        CMS_lumi.extraText = "Simulation"
        CMS_lumi.lumi_sqrtS = "14TeV"
        CMS_lumi.CMS_lumi(c, iPeriod, iPos)
        #drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")

        c.Modified()
        c.Update()
        c.SaveAs("%s_significance_%s_%s_%s.png" % (pileup, plotvar, version, sampleType))
