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

filenames = ["../../ZMM_PU200_910_pre3_relval/ZMM_PU200_910_pre3_relval.root", "../../ZMM_PU200_910_pre3_me0seed/ZMM_PU200_910_pre3_me0seed.root"]
#filenames = ["../../ZMM_PU200_pre5_me0seed/ZMM_PU200_pre5_me0seed.root"]
genTree = "MuonAnalyser/gen"
recoTree = "MuonAnalyser/reco"

pileup = "PU200"
"""
binning = [ [1500,0,20],
            [1500,0,50],
            [1500,0,20],
            [1500,0,50],
            [1500,0,2.5],
            [1500,0,1.5],
            [1500,0,20],
            [1500,0,50],
            [1500,0,20],
            [1500,0,50],
            [1500,0,2.5],
            [1500,0,1.5],
            [1500,0,20],
            [1500,0,50],
            [1500,0,20],
            [1500,0,50],
            [1500,0,2.5],
            [1500,0,1.5] ]

var = [ "fabs(muon_GE11deltaX)",
        "fabs(muon_GE11deltaY)",
        "fabs(muon_GE11pullX)",
        "fabs(muon_GE11pullY)",
        "fabs(muon_GE11dPhi)",
        "fabs(muon_GE11dEta)",
        "fabs(muon_GE21deltaX)",
        "fabs(muon_GE21deltaY)",
        "fabs(muon_GE21pullX)",
        "fabs(muon_GE21pullY)",
        "fabs(muon_GE21dPhi)",
        "fabs(muon_GE21dEta)",
        "fabs(muon_ME0deltaX)",
        "fabs(muon_ME0deltaY)",
        "fabs(muon_ME0pullX)",
        "fabs(muon_ME0pullY)",
        "fabs(muon_ME0dPhi)",
        "fabs(muon_ME0dEta)" ]
"""
var = [ "fabs(muon_ME0dPhi)",
        "fabs(muon_ME0dEta)",
        "fabs(muon_ME0pullPhi)",
        "fabs(muon_GE11dPhi)",
        "fabs(muon_GE11dEta)",
        "fabs(muon_GE11pullPhi)",
        "fabs(muon_GE21dPhi)",
        "fabs(muon_GE21dEta)",
        "fabs(muon_GE21pullPhi)" ]

binning = [ [1500,0,0.05],
            [1500,0,0.5],
            [1500,0,3.5],     
            [1500,0,0.05],
            [1500,0,0.5],
            [1500,0,3.5],     
            [1500,0,0.05],
            [1500,0,0.5],
            [1500,0,3.5] ]
    
for i, plotvar in enumerate(var):
    ### canvas ###
    canvasname = plotvar
    c = makeCanvas(canvasname, False)
    l = ROOT.TLegend(.18, .20, .33, .30)

    print "Canvas was made..."
    print "Starting loop"
    print "var: %s" % plotvar

    mg = ROOT.TMultiGraph()
    mg.SetTitle("ROC of %s" % plotvar)

    for filename in filenames:
        graphBinX = []
        graphBinY = []

        if "GE" in plotvar:
            h_sigEff = makeTH1(filename, genTree, "", binning[i], plotvar, "muon_isGEMMuon && %s!=100"%plotvar)
            h_bkgRej = makeTH1(filename, recoTree, "", binning[i], plotvar, "muon_isGEMMuon && !muon_signal && %s!=100"%plotvar)
        elif "ME0" in plotvar:
            h_sigEff = makeTH1(filename, genTree, "", binning[i], plotvar, "muon_isME0Muon && %s!=100"%plotvar)
            h_bkgRej = makeTH1(filename, recoTree, "", binning[i], plotvar, "muon_isME0Muon && !muon_signal && %s!=100"%plotvar)

        sigInt = h_sigEff.Integral(0, binning[i][0]+1)
        bkgInt = h_bkgRej.Integral(0, binning[i][0]+1)

        for j in range(binning[i][0]+2):
            sigEntries = h_sigEff.Integral(0, j) / sigInt
            bkgEntries = 1.0 - h_bkgRej.Integral(0, j) / bkgInt
            graphBinX.append(sigEntries)
            graphBinY.append(bkgEntries)

        arrX = np.asarray(graphBinX)
        arrY = np.asarray(graphBinY)

        graphROC = ROOT.TGraph(binning[i][0]+2, arrX, arrY)
        graphROC.SetLineWidth(2)
        if "me0seed" in filename:
            graphROC.SetLineColor(2)
            l.AddEntry(graphROC, "ME0seed", "l")
            print "    Added me0seed ROC graph"
        else:
            graphROC.SetLineColor(4)
            l.AddEntry(graphROC, "RelVal", "l")
            print "    Added RelVal ROC graph"
        mg.Add(graphROC)

    mg.Draw("AC")
    mg.SetTitle("%s"%plotvar)
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
    c.SaveAs("me0seedCompare_ROC_%s_910pre3_%s2.png" % (pileup, plotvar))
    print "    Done"






        
