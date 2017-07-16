import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
import numpy as np

print "Start"
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def overFlow(hist):
    nbins = hist.GetNbinsX()
    hist.SetBinContent(nbins, hist.GetBinContent(nbins)+hist.GetBinContent(nbins+1))
    hist.SetBinError(nbins, math.sqrt(hist.GetBinError(nbins)**2+hist.GetBinError(nbins+1)**2))
    #hist.SetBinContent(1, hist.GetBinContent(1)+hist.GetBinContent(0))
    #hist.SetBinError(1, math.sqrt(hist.GetBinError(1)**2+hist.GetBinError(0)**2))
    """
    nx = hist.GetNbinsX()+1
    x1 = hist.GetBinLowEdge(1)
    bw = hist.GetBinWidth(nx)
    x2 = hist.GetBinLowEdge(nx)+bw

    htmp = ROOT.TH1D("","",nx,x1,x2)
    for i in range(0,nx):
        htmp.Fill(htmp.GetBinCenter(i), hist.GetBinContent(i))
    htmp.SetEntries(hist.GetEntries())

    return copy.deepcopy(htmp)
    """

print "Initiating..."

#filenames = ["../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root", "../../ZMM_PU200_910_pre3_me0seed/ZMM_PU200_910_pre3_me0seed.root"]
filename = "../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root"
#filenames = ["../../ZMM_PU200_pre5_me0seed/ZMM_PU200_pre5_me0seed.root"]
genTree = "MuonAnalyser/gen"
recoTree = "MuonAnalyser/reco"

pileup = "PU200"
binning = [1500, 1.4, 3.4]
plotvar =  "fabs(muon.Eta())"
LT = ["loose", "tight"]
cutvals = ["New", "Latest"]

for lt in LT:
    ### canvas ###
    canvasname = lt
    c = makeCanvas(canvasname, False)
    l = ROOT.TLegend(.18, .20, .33, .30)

    print "Canvas was made..."
    print "Starting loop"
    print "type: %s" % lt

    mg = ROOT.TMultiGraph()
    mg.SetTitle("ROC of %s" % lt)

    for cutval in cutvals:
        graphBinX = []
        graphBinY = []
        if lt == "loose":
            if "New" in cutval:
                h_sigEff = makeTH1(filename, genTree, "", binning, plotvar, "muon_isME0Muon && fabs(muon_ME0dPhi)<0.015 && fabs(muon_ME0dEta)<0.05")
                h_bkgRej = makeTH1(filename, recoTree, "", binning, plotvar, "muon_isME0Muon && !muon_signal && fabs(muon_ME0dPhi)<0.015 && fabs(muon_ME0dEta)<0.05")

            else:
                h_sigEff = makeTH1(filename, genTree, "", binning, plotvar, "muon_isME0Muon && (fabs(muon_ME0deltaX)<1.6 || fabs(muon_ME0pullX)<1.4) && (fabs(muon_ME0deltaY)<6 || fabs(muon_ME0pullY)<8) && fabs(muon_ME0deltaDXDZ)<0.07 && fabs(muon_ME0deltaDYDZ)<0.46 && muon_ME0noRecHit > 2")
                h_sigBkg = makeTH1(filename, genTree, "", binning, plotvar, "!muon_signal && muon_isME0Muon && (fabs(muon_ME0deltaX)<1.6 || fabs(muon_ME0pullX)<1.4) && (fabs(muon_ME0deltaY)<6 || fabs(muon_ME0pullY)<8) && fabs(muon_ME0deltaDXDZ)<0.07 && fabs(muon_ME0deltaDYDZ)<0.46 && muon_ME0noRecHit > 2")
        if lt == "tight":
            if "Latest" in cutval:
                h_sigEff = makeTH1(filename, genTree, "", binning, plotvar, "muon_isME0Muon && fabs(muon_ME0dPhi)<0.004 && fabs(muon_ME0dEta)<0.035")
                h_bkgRej = makeTH1(filename, recoTree, "", binning, plotvar, "muon_isME0Muon && !muon_signal && fabs(muon_ME0dPhi)<0.004 && fabs(muon_ME0dEta)<0.035")

            else:
                h_sigEff = makeTH1(filename, genTree, "", binning, plotvar, "muon_isME0Muon && (fabs(muon_ME0deltaX)<0.7 || fabs(muon_ME0pullX)<1) && (fabs(muon_ME0deltaY)<4 || fabs(muon_ME0pullY)<3) && fabs(muon_ME0deltaDXDZ)<0.05 && fabs(muon_ME0deltaDYDZ)<0.42 && muon_ME0noRecHit > 3")
                h_sigBkg = makeTH1(filename, genTree, "", binning, plotvar, "!muon_signal && muon_isME0Muon && (fabs(muon_ME0deltaX)<0.7 || fabs(muon_ME0pullX)<1) && (fabs(muon_ME0deltaY)<4 || fabs(muon_ME0pullY)<3) && fabs(muon_ME0deltaDXDZ)<0.05 && fabs(muon_ME0deltaDYDZ)<0.42 && muon_ME0noRecHit > 3")

        sigInt = h_sigEff.Integral(0, binning[0]+1)
        bkgInt = h_bkgRej.Integral(0, binning[0]+1)

        for j in range(binning[0]+2):
            sigEntries = h_sigEff.Integral(0, j) / sigInt
            bkgEntries = 1.0 - h_bkgRej.Integral(0, j) / bkgInt
            graphBinX.append(sigEntries)
            graphBinY.append(bkgEntries)

        arrX = np.asarray(graphBinX)
        arrY = np.asarray(graphBinY)

        graphROC = ROOT.TGraph(binning[0]+2, arrX, arrY)
        graphROC.SetLineWidth(2)
        if "New" in cutval:
            graphROC.SetLineColor(2)
            l.AddEntry(graphROC, "New", "l")
            print "    Added New ROC graph"
        else:
            graphROC.SetLineColor(4)
            l.AddEntry(graphROC, "Latest", "l")
            print "    Added Latest ROC graph"
        mg.Add(graphROC)

    mg.Draw("AC")
    mg.GetXaxis().SetLimits(0, 1.05)
    mg.GetXaxis().SetTitle("Signal Efficiency")
    mg.GetYaxis().SetTitle("Background Rejection")
    mg.SetMinimum(0)
    mg.SetMaximum(1.05)
    l.Draw()
    print "    Graphs were drawn on a canvas..."

    ### CMS_lumi setting ###
    iPos = 0
    iPeriod = 0
    if (iPos == 0):
        CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "work in progress"
    CMS_lumi.lumi_sqrtS = "14TeV"
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    c.Modified()
    c.Update()
    c.SaveAs("NewLatestCompare_ROC_%s_910pre3_%s.png" % (pileup, lt))
    print "    Done"






        
