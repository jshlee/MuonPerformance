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
ME0Nocut  = { "isMuon": "muon_isME0Muon" }
ME0Loose  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","1.6"],     "dY"   : ["fabs(muon_ME0deltaY)","6"], 
              "pX"    : ["fabs(muon_ME0pullX)","1.4"],      "pY"   : ["fabs(muon_ME0pullY)","8"], 
              "dDXDZ" : ["fabs(muon_ME0deltaDXDZ)","0.07"], "dDYDZ": ["fabs(muon_ME0deltaDYDZ)","0.46"], 
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.014"],     "dEta" : ["fabs(muon_ME0dEta)","0.042"] }

ME0Tight  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","0.7"],     "dY"   : ["fabs(muon_ME0deltaY)","4"], 
              "pX"    : ["fabs(muon_ME0pullX)","1.0"],      "pY"   : ["fabs(muon_ME0pullY)","3"], 
              "dDXDZ" : ["fabs(muon_ME0deltaDXDZ)","0.05"], "dDYDZ": ["fabs(muon_ME0deltaDYDZ)","0.42"], 
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.005"],     "dEta" : ["fabs(muon_ME0dEta)","0.034"] }

ME0Loose1  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","1.6"],     "dY"   : ["fabs(muon_ME0deltaY)","6"], 
              "pX"    : ["fabs(muon_ME0pullX)","1.4"],      "pY"   : ["fabs(muon_ME0pullY)","8"], 
              "dDXDZ" : ["fabs(muon_ME0deltaDXDZ)","0.07"], "dDYDZ": ["fabs(muon_ME0deltaDYDZ)","0.46"] }
ME0Tight1  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","0.7"],     "dY"   : ["fabs(muon_ME0deltaY)","4"], 
              "pX"    : ["fabs(muon_ME0pullX)","1.0"],      "pY"   : ["fabs(muon_ME0pullY)","3"], 
              "dDXDZ" : ["fabs(muon_ME0deltaDXDZ)","0.05"], "dDYDZ": ["fabs(muon_ME0deltaDYDZ)","0.42"] }

ME0Loose2  = { "isMuon": "muon_isME0Muon",
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.014"],     "dEta" : ["fabs(muon_ME0dEta)","0.042"] }
ME0Tight2  = { "isMuon": "muon_isME0Muon",
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.005"],     "dEta" : ["fabs(muon_ME0dEta)","0.034"] }

ME0Loose3  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","1.6"],     "dY"   : ["fabs(muon_ME0deltaY)","6"], 
              "pX"    : ["fabs(muon_ME0pullX)","1.4"],      "pY"   : ["fabs(muon_ME0pullY)","8"], 
              "dDXDZ" : ["fabs(muon_ME0deltaDXDZ)","0.15"], "dDYDZ": ["fabs(muon_ME0deltaDYDZ)","0.8"], 
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.014"],     "dEta" : ["fabs(muon_ME0dEta)","0.042"] }

ME0Tight3  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","1.6"],     "dY"   : ["fabs(muon_ME0deltaY)","6"], 
              "pX"    : ["fabs(muon_ME0pullX)","1.4"],      "pY"   : ["fabs(muon_ME0pullY)","8"], 
              "dDXDZ" : ["fabs(muon_ME0deltaDXDZ)","0.15"], "dDYDZ": ["fabs(muon_ME0deltaDYDZ)","0.8"], 
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.005"],     "dEta" : ["fabs(muon_ME0dEta)","0.034"] }

ME0SelNew = { "isMuon": "muon_isME0MuonSelNew" }

ME0BDT = { "isMuon"   : "muon_isME0Muon",
           "BDT"      : ["muon_tmva_bdt", "-0.1283"] }


ME0LooseDef  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","1.6"],     "dY"   : ["fabs(muon_ME0deltaY)","6"], 
              "pX"    : ["fabs(muon_ME0pullX)","1.4"],      "pY"   : ["fabs(muon_ME0pullY)","8"],
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.5"] }

ME0TightDef  = { "isMuon": "muon_isME0Muon",
              "dX"    : ["fabs(muon_ME0deltaX)","3"],     "dY"   : ["fabs(muon_ME0deltaY)","3"], 
              "pX"    : ["fabs(muon_ME0pullX)","3"],      "pY"   : ["fabs(muon_ME0pullY)","3"],
              "dPhi"  : ["fabs(muon_ME0dPhi)","0.1"] }

GE11Nocut = { "isMuon": "muon_isGE11Muon" }

GE11Loose = { "isMuon": "muon_isGE11Muon",
              "dX"    : ["fabs(muon_GE11deltaX)","2.5"],    "dY"   : ["fabs(muon_GE11deltaY)","45"], 
              "pX"    : ["fabs(muon_GE11pullX)","25"],      "pY"   : ["fabs(muon_GE11pullY)","25"], 
              "dDXDZ" : ["fabs(muon_GE11deltaDXDZ)","0.5"], "dDYDZ": ["fabs(muon_GE11deltaDYDZ)","6"], 
              "dPhi"  : ["fabs(muon_GE11dPhi)","0.014"],    "dEta" : ["fabs(muon_GE11dEta)","0.27"] }
GE11Tight = { "isMuon": "muon_isGE11Muon",
              "dX"    : ["fabs(muon_GE11deltaX)","1.4"],    "dY"   : ["fabs(muon_GE11deltaY)","17"], 
              "pX"    : ["fabs(muon_GE11pullX)","1.6"],     "pY"   : ["fabs(muon_GE11pullY)","10"], 
              "dDXDZ" : ["fabs(muon_GE11deltaDXDZ)","0.27"],"dDYDZ": ["fabs(muon_GE11deltaDYDZ)","1.8"], 
              "dPhi"  : ["fabs(muon_GE11dPhi)","0.006"],    "dEta" : ["fabs(muon_GE11dEta)","0.10"] }

GE11Loose1 = { "isMuon": "muon_isGE11Muon",
              "dX"    : ["fabs(muon_GE11deltaX)","2.5"],    "dY"   : ["fabs(muon_GE11deltaY)","45"], 
              "pX"    : ["fabs(muon_GE11pullX)","25"],      "pY"   : ["fabs(muon_GE11pullY)","25"], 
              "dDXDZ" : ["fabs(muon_GE11deltaDXDZ)","0.5"], "dDYDZ": ["fabs(muon_GE11deltaDYDZ)","6"] }
GE11Tight1 = { "isMuon": "muon_isGE11Muon",
              "dX"    : ["fabs(muon_GE11deltaX)","1.4"],    "dY"   : ["fabs(muon_GE11deltaY)","17"], 
              "pX"    : ["fabs(muon_GE11pullX)","1.6"],     "pY"   : ["fabs(muon_GE11pullY)","10"], 
              "dDXDZ" : ["fabs(muon_GE11deltaDXDZ)","0.27"],"dDYDZ": ["fabs(muon_GE11deltaDYDZ)","1.8"] }

GE11Loose2 = { "isMuon": "muon_isGE11Muon",
              "dPhi"  : ["fabs(muon_GE11dPhi)","0.014"],    "dEta" : ["fabs(muon_GE11dEta)","0.27"] }
GE11Tight2 = { "isMuon": "muon_isGE11Muon",
              "dPhi"  : ["fabs(muon_GE11dPhi)","0.006"],    "dEta" : ["fabs(muon_GE11dEta)","0.10"] }

usrdic = [ME0Nocut,ME0SelNew, ME0Loose3, ME0Tight3, ME0BDT]
#usrdic = [GE11Nocut, GE11Loose, GE11Tight]
#usrdic = [ME0Nocut, ME0Loose1, ME0Tight1]
#usrdic = [GE11Nocut, GE11Loose1, GE11Tight1]
#usrdic = [ME0Nocut, ME0Loose2, ME0Tight2]
#usrdic = [GE11Nocut, GE11Loose2, GE11Tight2]


for draw in ["Eff", "Bkg"]:
    ### Canvas ###
    canvasname = draw+plotvar
    c = makeCanvas(canvasname, False)
    l = ROOT.TLegend(.75, .75, .95, .95)

    hlist = []
    if draw == "Eff":
        h_Effbkg = makeTH1(filename, gentree, "Effbkg", binning, plotvar, "")
        for i,opt in enumerate(usrdic):
            if i == 0: #nocut
                h_noCut = makeTH1(filename, gentree, "EffnoCut", binning, plotvar, opt["isMuon"])
                h_EffnoCut = ROOT.TEfficiency(h_noCut, h_Effbkg)
                NocutEntries = h_noCut.Integral()
                hlist.append(h_EffnoCut)
            elif i == 1: #SelNew
                h_Effsig = makeTH1(filename, gentree, "Effsig", binning, plotvar, opt["isMuon"])
                h_EffSelNew = ROOT.TEfficiency(h_Effsig, h_Effbkg)
                SelNewEntries = h_Effsig.Integral()
                hlist.append(h_EffSelNew)
            elif i == 4: #TMVA BDT
                h_Effsig = makeTH1(filename, gentree, "Effsig", binning, plotvar, 
                                   opt["isMuon"]   +                         "&&" +
                                   opt["BDT"][0]   + ">" + opt["BDT"][1])#   + "&&" +
                                   #opt["BDT"][0]   + "!=1") 
                h_EffBDT = ROOT.TEfficiency(h_Effsig, h_Effbkg)
                BDTEntries = h_Effsig.Integral()
                hlist.append(h_EffBDT)
            
            else:
                h_Effsig = makeTH1(filename, gentree, "Effsig", binning, plotvar,
                                   opt["isMuon"]   +                         "&&" +
                                   opt["dX"][0]    + "<" + opt["dX"][1]    + "&&" +
                                   opt["dY"][0]    + "<" + opt["dY"][1]    + "&&" +
                                   opt["pX"][0]    + "<" + opt["pX"][1]    + "&&" +
                                   opt["pY"][0]    + "<" + opt["pY"][1]    + "&&" +
                                   opt["dDXDZ"][0] + "<" + opt["dDXDZ"][1] + "&&" +
                                   opt["dDYDZ"][0] + "<" + opt["dDYDZ"][1] + "&&" +
                                   opt["dPhi"][0]  + "<" + opt["dPhi"][1]  + "&&" +
                                   opt["dEta"][0]  + "<" + opt["dEta"][1])
                if i == 2: #Loose
                    h_EffLoose = ROOT.TEfficiency(h_Effsig, h_Effbkg)
                    hlist.append(h_EffLoose)
                if i == 3: #Tight
                    h_EffTight = ROOT.TEfficiency(h_Effsig, h_Effbkg)
                    hlist.append(h_EffTight)

        print "Total Eff of SelNew = %f" % (SelNewEntries/NocutEntries)
        print "Total Eff of BDT    = %f" % (BDTEntries/NocutEntries)
        
        h_init = ROOT.TH1F("", "", 10, 1.4, 3.4)
        h_init.SetMaximum(1.05)
        h_init.SetMinimum(0)
        h_init.GetXaxis().SetTitle("Muon |#eta|")
        h_init.GetYaxis().SetTitle("Efficiency")
        h_init.Draw("l")

    if draw == "Bkg":
        c.SetLogy()
        for i,opt in enumerate(usrdic):
            if i == 0: #nocut
                h_noCut = makeTH1(filename, recotree, "BkgNoCut", binning, plotvar, "!muon_signal && " + opt["isMuon"])
                h_noCut.Scale(1/nevents)
                hlist.append(h_noCut)
            elif i == 1: #SelNew
                h_SelNew = makeTH1(filename, recotree, "BkgNoCut", binning, plotvar, "!muon_signal && " + opt["isMuon"])
                h_SelNew.Scale(1/nevents)
                hlist.append(h_SelNew)
            elif i == 4: #TMVA BDT
                h_BDT = makeTH1(filename, recotree, "BkgNoCut", binning, plotvar, 
                                   "!muon_signal"  +                         "&&" +
                                   opt["isMuon"]   +                         "&&" +
                                   opt["BDT"][0]   + ">" + opt["BDT"][1])#   + "&&" +
                                   #opt["BDT"][0]   + "!=1") 
                h_BDT.Scale(1/nevents)
                hlist.append(h_BDT)
            else:
                h_Bkgrt = makeTH1(filename, recotree, "Bkgrate", binning, plotvar, 
                                   "!muon_signal"  +                         "&&" +
                                   opt["isMuon"]   +                         "&&" +
                                   opt["dX"][0]    + "<" + opt["dX"][1]    + "&&" +
                                   opt["dY"][0]    + "<" + opt["dY"][1]    + "&&" +
                                   opt["pX"][0]    + "<" + opt["pX"][1]    + "&&" +
                                   opt["pY"][0]    + "<" + opt["pY"][1]    + "&&" +
                                   opt["dDXDZ"][0] + "<" + opt["dDXDZ"][1] + "&&" +
                                   opt["dDYDZ"][0] + "<" + opt["dDYDZ"][1] + "&&" +
                                   opt["dPhi"][0]  + "<" + opt["dPhi"][1]  + "&&" +
                                   opt["dEta"][0]  + "<" + opt["dEta"][1])
                if i == 2: #Loose
                    h_BkgLoose = h_Bkgrt
                    h_BkgLoose.Scale(1/nevents)
                    hlist.append(h_BkgLoose)
                if i == 3: #Tight
                    h_BkgTight = h_Bkgrt
                    h_BkgTight.Scale(1/nevents)
                    hlist.append(h_BkgTight)

    for i, h in enumerate(hlist):
        if i == 0:
            if draw == "Eff":
                l.AddEntry(h, "No Cuts", "lp")
                h.Draw("p0lsame")
                setMarkerStyle(h, 4, 23)
            if draw == "Bkg":
                l.AddEntry(h, "No Cuts", "f")
                h.SetFillColor(4)
                h.SetFillStyle(3003)
                h.GetXaxis().SetTitle("Muon |#eta|")
                h.GetYaxis().SetTitle("Background/event")
                h.Draw("histsame")
        if i == 1:
            if draw == "Eff":
                l.AddEntry(h, "SelNew", "lp")
                h.Draw("p0lsame")
                setMarkerStyle(h, 6, 21)
            if draw == "Bkg":
                l.AddEntry(h, "SelNew", "f")
                h.SetFillColor(6)
                h.SetFillStyle(3003)
                h.Draw("histsame")
        """
        if i == 2:
            if draw == "Eff":
                l.AddEntry(h, "Loose", "lp")
                h.Draw("p0lsame")
                setMarkerStyle(h, 3, 22)
            if draw == "Bkg":
                l.AddEntry(h, "Loose", "f")
                h.SetFillColor(3)
                h.SetFillStyle(3005)
                h.Draw("histsame")
        if i == 3:
            if draw == "Eff":
                l.AddEntry(h, "Tight", "lp")
                h.Draw("p0lsame")
                setMarkerStyle(h, 2, 20)
            if draw == "Bkg":
                l.AddEntry(h, "Tight", "f")
                h.SetFillColor(2)
                h.SetFillStyle(3004)
                h.Draw("histsame")
        """
        if i == 4:
            if draw == "Eff":
                l.AddEntry(h, "BDT", "lp")
                h.Draw("p0lsame")
                setMarkerStyle(h, 7, 34)
            if draw == "Bkg":
                l.AddEntry(h, "BDT", "f")
                h.SetFillColor(7)
                h.SetFillStyle(3006)
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
    c.SaveAs("OptPlots/%s_ME0_Opt%s_910pre3_BDTvsSelNew.png" % (pileup,draw))
