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
filenames = ["../../ZMM_PU200_pre4/ZMM_pu200_pre4"]
treename = "MuonAnalyser/reco"
binning = [# [40,-30,10],    # edge x
           # [40,-30,10],    # edge y
           # [150,-50,100],  # distance
           # [10,0,10],      # distance Err
           # [70,-35,35],    # chamber x
           # [90,-45,45],    # chamber y
           # [100,0,50],     # chamber xErr
           # [100,0,50],     # chamber yErr
           # [20,-1,1],      # chamber dxdz
           # [20,-1,1],      # chamber dydz
           # [20,0,1],       # chamber dxdzErr
           # [20,0,1],       # chamber dydzErr
           # [70,-35,35],    # segment x
           # [90,-45,45],    # segment y
           # [100,0,50],     # segment xErr
           # [100,0,50],     # segment yErr
           # [20,-1,1],      # segment dxdz
           # [20,-1,1],      # segment dydz
           # [20,0,1],       # segment dxdzErr
           # [20,0,1],       # segment dydzErr
            [81,0,8],       # delta x
            [101,0,100],     # delta y
            [151,0,150],     # delta x div by seg xErr
            [151,0,150],     # delta y div by seg yErr
            [31,0,30],      # delta x div by cham xErr
            [151,0,150],     # delta y div by cham yErr
            [101,0,1],       # delta dxdz
            [51,0,2],       # delta dydz
            [51,0,250],     # delta dxdz div by seg xErr
            [21,0,20],      # delta dydz div by seg yErr
            [101,0,200],     # delta dxdz div by cham xErr
            [21,0,20] ]     # delta dydz div by cham yErr
            
            

var = [ #"recoMuon_edgeXME0",
        #"recoMuon_edgeYME0",
        #"recoMuon_distance",
        #"recoMuon_distErr",
        #"recoMuon_chamberMatchXME0",
        #"recoMuon_chamberMatchYME0",
        #"recoMuon_chamberMatchXErrME0",
        #"recoMuon_chamberMatchYErrME0",
        #"recoMuon_chamberMatchDXDZME0",
        #"recoMuon_chamberMatchDYDZME0",
        #"recoMuon_chamberMatchDXDZErrME0",
        #"recoMuon_chamberMatchDYDZErrME0",
        #"recoMuon_segmentMatchXME0",
        #"recoMuon_segmentMatchYME0",
        #"recoMuon_segmentMatchXErrME0",
        #"recoMuon_segmentMatchYErrME0",
        #"recoMuon_segmentMatchDXDZME0",
        #"recoMuon_segmentMatchDYDZME0",
        #"recoMuon_segmentMatchDXDZErrME0",
        #"recoMuon_segmentMatchDYDZErrME0",
        "recoMuon_deltaXME0",
        "recoMuon_deltaYME0",
        "recoMuon_deltaXDivBySegErrME0",
        "recoMuon_deltaYDivBySegErrME0",
        "recoMuon_deltaXDivByChamErrME0",
        "recoMuon_deltaYDivByChamErrME0",
        "recoMuon_deltaDXDZME0",
        "recoMuon_deltaDYDZME0",
        "recoMuon_deltaDXDZDivBySegErrME0",
        "recoMuon_deltaDYDZDivBySegErrME0",
        "recoMuon_deltaDXDZDivByChamErrME0",
        "recoMuon_deltaDYDZDivByChamErrME0" ]

for filename in filenames:
    """
    if "pu0" in filename:
        pileup = "pu0"
    if "pu200" in filename:
        pileup = "pu200"
    """
    if "TTbar" in filename:
        pileup = "TTbar_pu0"
    elif "QCD" in filename:
        pileup = "QCD_pu0"
    elif "pu200" in filename:
        pileup = "pu200"
    else:
        pileup = "pu0"

    for i,plotvar in enumerate(var):
        ### Canvas ###
        canvasname = ""+pileup+plotvar
        c = makeCanvas(canvasname, False)
        """
        if "Err" in plotvar:
            l = ROOT.TLegend(.8, .8, .95, .90)
        else:
            l = ROOT.TLegend(.2, .8, .35, .90)
        """
        l = ROOT.TLegend(.8, .8, .95, .90)
        c.SetLogy()
        #c.SetMinimum(0.1)

        histList = []
        #h_dummy = ROOT.TH1D("dummy", "dummy", binning[0], binning[1], binning[2])
        #histList.append(h_dummy)

        ### Make histo ###
        h_comp = makeTH1(filename+".root", treename, pileup, binning[i], plotvar, 
                    #"recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && recoMuon_signal")
                    "recoMuon.Pt()>5 && recoMuon_isME0Muon &&abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && recoMuon_signal")
        scale_comp = 1 / (h_comp.Integral())
        h_comp.Scale(scale_comp)
        h_comp.GetXaxis().SetTitle("%s"%plotvar)
        h_comp.GetXaxis().SetTitleSize(0.05)
        h_comp.GetYaxis().SetTitle(pileup + " Normalized")
        h_comp.SetLineColor(4)
        histList.append(h_comp)

        h_fake = makeTH1(filename+".root", treename, pileup, binning[i], plotvar, 
                    "recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && !recoMuon_signal")
        scale_fake = 1 / (h_fake.Integral())
        h_fake.Scale(scale_fake)
        h_fake.SetLineColor(3)
        histList.append(h_fake)

        sigma = 0
        cut_bin = 0
        for j, h in enumerate(histList):
            if h == h_comp:
                for k in range(0, binning[i][0]):
                    sigma += h.GetBinContent(k)
                    if sigma > 0.956: # 1 sigma of true cut
                        cut_bin= k * binning[i][2]/(binning[i][0]-1.)
                        break
                    else:
                        continue
                l.AddEntry(h, "True", "l")
                l.AddEntry(0, "2#sigma cut = %s"%cut_bin, "")
                l.AddEntry(0, "%0.4f" %sigma, "")
            if h == h_fake:
                l.AddEntry(h, "Fake", "l")
            h.Draw("histsame")
        indicator = ROOT.TLine(cut_bin,0,cut_bin,1)
        indicator.SetLineStyle(2)
        indicator.Draw()
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
        c.SaveAs("%s_%s.png" % (pileup, plotvar))
