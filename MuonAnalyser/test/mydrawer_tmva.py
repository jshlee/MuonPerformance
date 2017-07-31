import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(1.2)
    h.SetLineColor(color)
    h.SetLineWidth(2)

def getEff(filename,treename,title,binning,plotvar,dencut,numcut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,dencut)
    h2 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    hEff = ROOT.TEfficiency(h2,h1)
    hEff = delLargeErrorPoints(hEff,h1,h2,binning)
    hEff.SetTitle(title)
    return copy.deepcopy(hEff)

def delLargeErrorPoints(hEff, h1, h2, binning):
    maxY = hEff.GetEfficiency(1)
    for j in range(binning[0]):
        if maxY < hEff.GetEfficiency(j+1): maxY = hEff.GetEfficiency(j+1)
    for i in range(1,binning[0]+1):
        if hEff.GetEfficiencyErrorLow(i) > maxY/4.5:
            h1.SetBinContent(i,0)
            h2.SetBinContent(i,0)
    hEff = ROOT.TEfficiency(h2,h1)
    return copy.deepcopy(hEff)

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(62)
    tex2.SetTextSize(0.03)
    tex2.DrawLatex(0.25, 0.88, samplename)

def divByNevents(filedir, hist):
    tfile = ROOT.TFile(filedir)
    nevents = tfile.Get("MuonAnalyser/nevents").Integral()
    hist.Scale(1/nevents)

def draw(h_init, y_name, hlists, name, text, binning):
    #Plot style
    setMarkerStyle(hlists[0], 4, 20) #blue, circle
    setMarkerStyle(hlists[1], 2, 34) #black, cross

    #Set canvas
    #canv = makeCanvas(plotvar+name, False)
    #setMargins(canv, False)
    canv = ROOT.TCanvas()
    canv.SetGrid()
    h_init.GetYaxis().SetTitle(y_name)
    h_init.Draw()
    drawSampleName(text)

    #Legend and drawing
    leg = ROOT.TLegend(0.25,0.74,0.45,0.84)
    lines=[]
    for h in hlists:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"p")

    hlists[0].Draw("e1same")
    #for line in lines: line.Draw("same")
    leg.SetTextFont(62)
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw()

    #CMS_lumi setting
    iPos = 0
    iPeriod = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_sqrtS = "14 TeV"
    CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

    canv.Modified()
    canv.Update()
    canv.SaveAs("%s_%s_%s.png"%(plotvar,muonid,name))


datadir = './'
#datadir = '/xrootd/store/user/tt8888tt/muon/9_1_1/'
filename = "ttbartotal.root"
#filename = "zmmtotal.root"

muonid = sys.argv[1]
binning_l = [[10,5,105],[12,0,2.4],[12,-2.4,2.4],[15,-3,3],[20,60,260],[10,0,10],[20,60,260]]
#rangecut = "abs(muon.Eta())<2.4"
rangecut = "muon.Pt()>5&&abs(muon.Eta())<2.4"
    
#Set extra text
samplename = "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}"
text = samplename+", p_{T} > 5 GeV, |#eta| < 2.4"

if "Tight" in muonid:
    #origid = "isTight"
    origid = "isTightCustom"
    tmvaid = "tmva_bdt>-0.165" #0.11
elif muonid == "Loose":
    origid = "isLoose"
    tmvaid = "tmva_bdt>-0.195" #0.167

for i, plotvar in enumerate(["muon.Pt()", "abs(muon.Eta())", "muon.Eta()", "muon.Phi()", "nvertex", "pu_density/2", "pu_numInteractions"]):
    #Efficiency
    hl_eff = []
    hl_eff.append(getEff(datadir+filename, "MuonAnalyser/gen", "PhaseII ", binning_l[i], plotvar, rangecut, "%s&&muon_%s"%(rangecut,origid)))
    hl_eff.append(getEff(datadir+filename, "MuonAnalyser/gen", "PhaseII  (TMVA ID)", binning_l[i], plotvar, rangecut, "%s&&muon_%s"%(rangecut,tmvaid)))

    #Backgorund rate
    hl_bkg = []
    if "()" in plotvar:
        hl_bkg.append(makeTH1(datadir+filename, "MuonAnalyser/reco", "PhaseII ", binning_l[i], plotvar, "%s&&!muon_signal&&muon_%s"%(rangecut,origid)))
        hl_bkg.append(makeTH1(datadir+filename, "MuonAnalyser/reco", "PhaseII  (TMVA ID)", binning_l[i], plotvar, "%s&&!muon_signal&&muon_%s"%(rangecut,tmvaid)))
        for j in range(len(hl_bkg)):
            divByNevents(datadir+filename, hl_bkg[j])

    else:
        break
        hl_bkg = []
        hl_bkg.append(getEff(datadir+filename, "MuonAnalyser/reco", "PhaseII ", binning_l[i], plotvar, "muon_no==1", "%s&&!muon_signal&&muon_%s"%(rangecut,origid)))
        hl_bkg.append(getEff(datadir+filename, "MuonAnalyser/reco", "PhaseII  (TMVA ID)", binning_l[i], plotvar, "muon_no==1", "%s&&!muon_signal&&muon_%s"%(rangecut,tmvaid)))
    if "density" in plotvar: plotvar = plotvar.split('/')[0]

    print "id  ", origid, hl_bkg[0].Integral(-1,2)
    print "bdt ", tmvaid, hl_bkg[1].Integral(-1,2)

    #Set X axis name
    x_name = "Muon "
    if "Pt"  in plotvar: x_name = x_name+"p_{T}"
    elif "Eta" in plotvar:
        x_name = x_name+"#eta"
        if "abs" in plotvar: x_name = "|"+x_name+"|"
    elif "Phi" in plotvar: x_name = x_name+"#phi"
    elif "vertex" in plotvar : x_name = "Number of vertex"
    elif "density" in plotvar : x_name = "PU density (number of PU per mm)"
    elif "Interaction" in plotvar : x_name = "Number of interactions"

    #Set init histo
    h_init = ROOT.TH1F("","",binning_l[i][0],binning_l[i][1],binning_l[i][2])

    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitleOffset(1)
    h_init.GetXaxis().SetTitleOffset(1.1)
    h_init.GetYaxis().SetTitleOffset(1.5)
    h_init.GetXaxis().SetTitleSize(0.05)
    h_init.GetXaxis().SetLabelSize(0.038)
    h_init.GetYaxis().SetTitleSize(0.05)
    h_init.GetYaxis().SetLabelSize(0.038)
    h_init2 = h_init.Clone()

    #Set Y axis name
    y_name = muonid+" Muon "
    if "Custom" in muonid: y_name = "Tight Muon "

    h_init.SetMaximum(1.1)
    h_init.SetMinimum(0.8)
    draw(h_init, y_name+"Efficiency", hl_eff, "eff", text, binning_l[i])

    if "()" in plotvar: h_init2.SetMaximum(max(h.GetMaximum() for h in hl_bkg)*2)
    else:
        maxY = hl_bkg[1].GetEfficiency(0)
        for j in range(binning_l[i][0]):
            x = hl_bkg[1].GetEfficiency(j)
            if maxY < x: maxY = x
        h_init2.SetMaximum(maxY*2.2)
    draw(h_init2, "Average "+y_name+"Bkg Multiplicity", hl_bkg, "bkgrate", text, binning_l[i])

