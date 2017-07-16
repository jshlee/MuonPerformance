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

#filenames = ["../../pu0/TTbar_pu0","../../pu0/QCD_pu0","../../pu0/pu0","../../pu200/pu200"]
filenames = ["../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root"]
gentree = "MuonAnalyser/gen"

binning = [ [21,-1,1] ] 

var = [ "muon_ME0dPhi" ]

for filename in filenames:
    pileup = "pu200"
   
    for i,plotvar in enumerate(var):
        ### Canvas ###
        canvasname = ""+pileup+plotvar
        c = makeCanvas(canvasname, False)
        c.SetLogy()
        histList = []
        ### Make histo ###
        h = makeTH1(filename, gentree, pileup, binning[i], plotvar, "muon_isME0Muon")
        h.Draw("hist")       

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
        c.SaveAs("%s_%s.png" % (pileup, plotvar))
