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
    hEff = ROOT.TEfficiency(h2,h1)
    hEff.SetTitle(title)
    return copy.deepcopy(hEff)

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.18, 0.85, samplename)

def draw(h_init, y_name, hlists, name, text):
    #Plot style
    setMarkerStyle(hlists[0], 4, 20) #blue, circle
    setMarkerStyle(hlists[1], 1, 34) #black, cross
    setMarkerStyle(hlists[2], 2, 21) #red, square

    #Set canvas
    canv = makeCanvas(plotvar+name, False)
    setMargins(canv, False)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.Draw()
    drawSampleName(text)

    #Legend and drawing
    legTop = ROOT.TLegend(0.6,0.6,0.8,0.75)
    legBot = ROOT.TLegend(0.6,0.2,0.8,0.35)
    if "Rate" in y_name: leg = legTop
    else: leg = legBot
    for h in hlists:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"p")
    leg.SetTextFont(61)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw()

    #CMS_lumi setting
    iPos = 0
    iPeriod = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.lumi_sqrtS = "14 TeV"
    CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

    canv.Modified()
    canv.Update()
    canv.SaveAs("%s_%s_%s.png"%(plotvar,id,name))


datadir = '/xrootd/store/user/tt8888tt/muon/'
filenames = ["zmm.root", "zmm140.root", "zmm200.root"]
#filenames = ["backup2/zmm.root", "zmm.root", "zmm.root"]

id = sys.argv[1]
binning_l = [[20,5,105],[30,0,3],[30,-3,3],[30,-3,3],[40,60,260]]
rangecut = "muon.Pt()>5&&abs(muon.Eta())<2.4"

for i, plotvar in enumerate(["muon.Pt()", "abs(muon.Eta())", "muon.Eta()", "muon.Phi()", "nvertex"]):
    #Get histos
    hl_eff = []
    hl_eff.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "PU 0",   binning_l[i], plotvar,  rangecut, "%s&&muon_is%s"%(rangecut,id)))
    hl_eff.append(getEff(datadir+filenames[1], "MuonAnalyser/gen", "PU 140", binning_l[i], plotvar, rangecut, "%s&&muon_is%s"%(rangecut,id)))
    hl_eff.append(getEff(datadir+filenames[2], "MuonAnalyser/gen", "PU 200", binning_l[i], plotvar, rangecut, "%s&&muon_is%s"%(rangecut,id)))

    #Backgorund rate
    hl_bkg = []
    tfile = ROOT.TFile(datadir+filenames[0])
    nevents_pu0 = tfile.Get("MuonAnalyser/nevents").Integral()
    tfile = ROOT.TFile(datadir+filenames[1])
    nevents_pu140 = tfile.Get("MuonAnalyser/nevents").Integral()
    tfile = ROOT.TFile(datadir+filenames[2])
    nevents_pu200 = tfile.Get("MuonAnalyser/nevents").Integral()
    hl_bkg.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "PU 0",   binning_l[i], plotvar, "%s&&!muon_signal&&muon_is%s"%(rangecut,id)))
    hl_bkg.append(makeTH1(datadir+filenames[1], "MuonAnalyser/reco", "PU 140", binning_l[i], plotvar, "%s&&!muon_signal&&muon_is%s"%(rangecut,id)))
    hl_bkg.append(makeTH1(datadir+filenames[2], "MuonAnalyser/reco", "PU 200", binning_l[i], plotvar, "%s&&!muon_signal&&muon_is%s"%(rangecut,id)))
    hl_bkg[0].Scale(1/nevents_pu0)
    hl_bkg[1].Scale(1/nevents_pu140)
    hl_bkg[2].Scale(1/nevents_pu200)

    #Fake rate
    hl_fake = []
    hl_fake.append(getEff(datadir+filenames[0], "MuonAnalyser/reco", "PU0",   binning_l[i], plotvar, "%s&&muon_is%s"%(rangecut,id), "%s&&!muon_signal&&muon_is%s"%(rangecut,id)))
    hl_fake.append(getEff(datadir+filenames[1], "MuonAnalyser/reco", "PU140", binning_l[i], plotvar, "%s&&muon_is%s"%(rangecut,id), "%s&&!muon_signal&&muon_is%s"%(rangecut,id)))
    hl_fake.append(getEff(datadir+filenames[2], "MuonAnalyser/reco", "PU200", binning_l[i], plotvar, "%s&&muon_is%s"%(rangecut,id), "%s&&!muon_signal&&muon_is%s"%(rangecut,id)))

    #Set X axis name
    x_name = "Muon "
    if "Pt"  in plotvar: x_name = x_name+"p_{T}"
    elif "Eta" in plotvar: x_name = x_name+"|#eta|"
    elif "Phi" in plotvar: x_name = x_name+"#phi"
    else : x_name = "Number of vertex"
    
    #Set extra text
    samplename = "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}"
    text = samplename+", p_{T} > 5 GeV, |#eta| < 2.4"

    #Set init histo
    h_init = ROOT.TH1F("","",binning_l[i][0],binning_l[i][1],binning_l[i][2])

    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitleOffset(1)
    h_init2 = h_init.Clone()

    #Set Y axis name
    y_name = id+" "
    if "Custom" in id: y_name = "TightNoVtx Muon "

    h_init.SetMaximum(1.15)
    h_init.SetMinimum(0.4)
    draw(h_init, y_name+"Efficiency", hl_eff, "eff", text)

    h_init2.GetYaxis().SetLabelSize(0.035)
    h_init2.GetYaxis().SetTitleSize(0.050)
    h_init2.GetYaxis().SetTitleOffset(1.2)
    h_init3 = h_init2.Clone()
    h_init2.SetMaximum(max(h.GetMaximum() for h in hl_bkg)*2.5)
    draw(h_init2, y_name+"Background Rate", hl_bkg, "bkgrate", text)
    
    maxY = hl_fake[2].GetEfficiency(0)
    for j in range(binning_l[i][0]):
        x = hl_fake[2].GetEfficiency(j)
        if maxY > x: continue
        maxY = x
    h_init3.SetMaximum(maxY*2.5)
    draw(h_init3, y_name+"Fake Rate", hl_fake, "fakerate", text)
