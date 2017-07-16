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
    h.SetMarkerSize(1.5)
    h.SetLineColor(color)
    h.SetLineWidth(2)

def getEff(filename,treename,title,binning,plotvar,dencut,numcut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,dencut)
    h2 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    h2.Divide(h1)
    return copy.deepcopy(h2)

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.18, 0.8, samplename)

#datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
#filenames = ["../pu0_/pu0.root","../pu200_/pu200.root"]
filenames = ["../../TenMuPU200/pu200.root"]
binning_l = [[12,5,65],[32,1.8,3.0],[12,-3,3],[12,5,65],[32,1.8,3.0],[12,-3,3]]

for i, plotvar in enumerate(["genMuon.Pt()", "abs(genMuon.Eta())", "genMuon.Phi()", "recoMuon.Pt()", "abs(recoMuon.Eta())", "recoMuon.Phi()"]):
    hlist = []
    for j, filename in enumerate(filenames):        
        #Get histos
        #if "pu0" in filename:
        if "pu0" in filename:
            if "genMuon" in plotvar:
                if ("Pt()" in plotvar) or ("Phi()" in plotvar):
                    h_ph2pu0 = getEff(filename, "MuonAnalyser", "PU0", binning_l[i], plotvar,
                                      "gen/genMuon.Pt()>5 && abs(gen/genMuon.Eta())>2.0 && abs(gen/genMuon.Eta())<2.8",
                                      "gen/genMuon.Pt()>5 && gen/genMuon_isME0Muon && abs(gen/genMuon.Eta())>2.0 && abs(gen/genMuon.Eta())<2.8")
                    hlist.append(h_ph2pu0)

                elif "Eta()" in plotvar:
                    #h_ph2pu0 = getEff(filename, "MuonAnalyser/gen", "PU0", binning_l[i], plotvar, 
                    #                  "genMuon.Pt()>5",
                    #                  "genMuon.Pt()>5 && genMuon_isME0Muon")
                    h_ph2pu0 = getEff(filename, "MuonAnalyser/gen", "PU0", binning_l[i], plotvar, 
                                      #"genMuon.Pt()>5 && abs((genMuon.X())/(genMuon.Z()))<0.7",
                                      #"genMuon.Pt()>5 && abs((genMuon.X())/(genMuon.Z()))<0.7 && genMuon_isME0Muon")
                                      "genMuon.Pt()>5",
                                      "genMuon.Pt()>5 && genMuon_isME0Muon")
                    hlist.append(h_ph2pu0)

            if "recoMuon" in plotvar:
                recoCut = "!recoMuon_signal"
                if ("Pt()" in plotvar) or ("Phi()" in plotvar):
                    h_ph2pu0 = getEff(filename, "MuonAnalyser/reco", "PU0", binning_l[i], plotvar, 
                                      "recoMuon.Pt()>5 && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && recoMuon_signal",
                                      "recoMuon.Pt()>5 && !recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && recoMuon_signal")
                    hlist.append(h_ph2pu0)

                elif "Eta()" in plotvar:
                    h_ph2pu0 = getEff(filename, "MuonAnalyser/reco", "PU0", binning_l[i], plotvar, 
                                      "recoMuon.Pt()>5 && recoMuon_signal",
                                      "recoMuon.Pt()>5 && !recoMuon_isME0Muon && recoMuon_signal")
                    hlist.append(h_ph2pu0)

        if "pu200" in filename:
            recoCut = "recoMuon_deltaXME0<2.8 && recoMuon_deltaYME0<56 && recoMuon_deltaXDivBySegErrME0<75 && recoMuon_deltaYDivBySegErrME0<70 && recoMuon_deltaXDivByChamErrME0<15 && recoMuon_deltaYDivByChamErrME0<63 && recoMuon_deltaDXDZME0<0.25 && recoMuon_deltaDYDZME0<0.4 && recoMuon_deltaDXDZDivBySegErrME0<85 && recoMuon_deltaDYDZDivBySegErrME0<6 && recoMuon_deltaDXDZDivByChamErrME0<64 && recoMuon_deltaDYDZDivByChamErrME0<6"
            #genCut = "genMuon_deltaXME0<2.8 && genMuon_deltaYME0<56 && genMuon_deltaXErrME0<15 && genMuon_deltaYErrME0<63 && genMuon_deltaDXDZME0<0.25 && genMuon_deltaDYDZME0<0.4 && genMuon_deltaDXDZErrME0<64 && genMuon_deltaDYDZErrME0<6"
            non = "999999"
            #genCut = "genMuon_deltaXME0!=%s && genMuon_deltaYME0!=%s && genMuon_deltaXErrME0!=%s && genMuon_deltaYErrME0!=%s && genMuon_deltaDXDZME0!=%s && genMuon_deltaDYDZME0!=%s && genMuon_deltaDXDZErrME0!=%s && genMuon_deltaDYDZErrME0!=%s"%(non,non,non,non,non,non,non,non)
            genCut="genMuon_deltaXME0<100"
            if "genMuon" in plotvar:
                """
                if ("Pt()" in plotvar) or ("Phi()" in plotvar):
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/gen", "PU200", binning_l[i], plotvar,
                                      "genMuon.Pt()>5 && abs(genMuon.Eta())>2.0 && abs(genMuon.Eta())<2.8 && %s"%genCut,
                                      "genMuon.Pt()>5 && genMuon_isME0Muon && abs(genMuon.Eta())>2.0 && abs(genMuon.Eta())<2.8 && %s"%genCut)
                    hlist.append(h_ph2pu200)

                elif "Eta()" in plotvar:
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/gen", "PU200", binning_l[i], plotvar, 
                                      "genMuon.Pt()>5 && %s"%genCut,
                                      "genMuon.Pt()>5 && genMuon_isME0Muon && %s"%genCut)
                    hlist.append(h_ph2pu200)
                """
                if ("Pt()" in plotvar) or ("Phi()" in plotvar):
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/gen", "PU200", binning_l[i], plotvar,
                                      "genMuon.Pt()>5 && abs(genMuon.Eta())>2.0 && abs(genMuon.Eta())<2.8",
                                      "genMuon.Pt()>5 && genMuon_isME0Muon && abs(genMuon.Eta())>2.0 && abs(genMuon.Eta())<2.8")
                    hlist.append(h_ph2pu200)

                elif "Eta()" in plotvar:
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/gen", "PU200", binning_l[i], plotvar, 
                                      "genMuon.Pt()>5",
                                      "genMuon.Pt()>5 && genMuon_isME0Muon")
                    hlist.append(h_ph2pu200)
            """
            if "recoMuon" in plotvar:
                if ("Pt()" in plotvar) or ("Phi()" in plotvar):
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/reco", "PU200", binning_l[i], plotvar, 
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8",
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && !recoMuon_signal")
                    hlist.append(h_ph2pu200)

                elif "Eta()" in plotvar:
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/reco", "PU200", binning_l[i], plotvar, 
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon",
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon && !recoMuon_signal")
                    hlist.append(h_ph2pu200)
            """
            if "recoMuon" in plotvar:
                if ("Pt()" in plotvar) or ("Phi()" in plotvar):
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/reco", "PU200", binning_l[i], plotvar, 
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && %s"%recoCut,
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon && abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8 && !recoMuon_signal && %s"%recoCut)
                    hlist.append(h_ph2pu200)

                elif "Eta()" in plotvar:
                    h_ph2pu200 = getEff(filename, "MuonAnalyser/reco", "PU200", binning_l[i], plotvar, 
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon && %s"%recoCut,
                                      "recoMuon.Pt()>5 && recoMuon_isME0Muon && !recoMuon_signal && %s"%recoCut)
                    hlist.append(h_ph2pu200)

        
    #Set init histo
    nbins = h_ph2pu200.GetNbinsX()
    h_init = ROOT.TH1F("","",nbins,h_ph2pu200.GetBinLowEdge(1),h_ph2pu200.GetBinLowEdge(nbins+1))
    """
    nbins = h_ph2pu0.GetNbinsX()
    h_init = ROOT.TH1F("","",nbins,h_ph2pu0.GetBinLowEdge(1),h_ph2pu0.GetBinLowEdge(nbins+1))
    """

    #Set axis
    x_name = "#mu "
    if "Pt"  in plotvar: x_name = x_name+"p_{T}"
    if "Eta" in plotvar: x_name = x_name+"|#eta|"
    if "Phi" in plotvar: x_name = x_name+"#phi"

    y_name = "ME0 #mu "
    if "genMuon" in plotvar:
        h_init.SetMaximum(1.1)
        h_init.SetMinimum(0)
        y_name = y_name+"Efficiency"
    if "recoMuon" in plotvar:
        #h_init.SetMaximum(max(h.GetMaximum() for h in hlist)*2.5)
        h_init.SetMaximum(1.2)
        h_init.SetMinimum(0)
        h_init.GetYaxis().SetLabelSize(0.035)
        h_init.GetYaxis().SetTitleOffset(1.2)
        y_name = y_name+"Fake Rate"
    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.GetYaxis().SetTitleOffset(1)

    """
    ############ Plot design ##############
    #Plot style
    for h in hlist:
        if "genMuon" in plotvar: 
            if "pu0" in h.GetTitle():
                setMarkerStyle(h, 4, 20) #blue, circle
            if "pu200" in h.GetTitle():
                setMarkerStyle(h, 2, 34) #red, cross
        if "recoMuon" in plotvar:
            if "pu0" in h.GetTitle():
                setMarkerStyle(h, 7, 20) #skyblue, circle
            if "pu200" in h.GetTitle():
                setMarkerStyle(h, 6, 34) #pink, cross
    #setMarkerStyle(h_ph2pu140, 1, 34) #black, cross
    #setMarkerStyle(h_ph2pu200, 2, 21) #red, square
    #setMarkerStyle(h_run2, 3, 22) #green, triangle
    """

    #Set canvas
    canv = makeCanvas(plotvar, False)
    setMargins(canv, False)
    h_init.Draw()
    if "Eta()" in plotvar:
        drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")
    else:
        drawSampleName("#splitline{Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV}{2.0<|#eta|<2.8}")
        

    #Legend and drawing
    legTop = ROOT.TLegend(0.6,0.7,0.85,0.82)
    legBot = ROOT.TLegend(0.6,0.2,0.85,0.32)
    #if "genMuon"  in plotvar: leg = legBot
    #if "recoMuon" in plotvar: leg = legTop

    ############ Plot design ##############
    #Plot style
    for h in hlist:
        if "genMuon" in plotvar:
            leg = legBot 
            if "PU0" in h.GetTitle():
                h.Draw("same")
                setMarkerStyle(h, 4, 20) #blue, circle
                leg.AddEntry(h, h.GetTitle(), "p")
            if "PU200" in h.GetTitle():
                h.Draw("same")
                setMarkerStyle(h, 2, 34) #red, cross
                leg.AddEntry(h, h.GetTitle(), "p")
        if "recoMuon" in plotvar:
            leg = legTop
            if "PU0" in h.GetTitle():
                h.Draw("same")
                setMarkerStyle(h, 7, 20) #skyblue, circle
                leg.AddEntry(h, h.GetTitle(), "p")
            if "PU200" in h.GetTitle():
                h.Draw("same")
                setMarkerStyle(h, 6, 34) #pink, cross
                leg.AddEntry(h, h.GetTitle(), "p")
    #setMarkerStyle(h_ph2pu140, 1, 34) #black, cross
    #setMarkerStyle(h_ph2pu200, 2, 21) #red, square
    #setMarkerStyle(h_run2, 3, 22) #green, triangle

    #for h in hlist:
    #    h.Draw("same")
    #    leg.AddEntry(h,h.GetTitle(),"p")
    leg.SetTextFont(61)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.Draw()

    #CMS_lumi setting
    iPos = 0
    iPeriod = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Work in progress"
    CMS_lumi.lumi_sqrtS = "14 TeV"
    CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

    canv.Modified()
    canv.Update()

    if "gen" in plotvar:
        canv.SaveAs("%s_%s_%s.png" % ("pu200","efficiency",plotvar))
    if "reco" in plotvar:
        canv.SaveAs("%s_%s_%s.png" % ("pu200","fake",plotvar))
    """
    if "gen" in plotvar:
        canv.SaveAs("%s_%s_%s.png" % ("pu0","efficiency",plotvar))
    if "reco" in plotvar:
        canv.SaveAs("%s_%s_%s.png" % ("pu0","fake",plotvar))
    """
