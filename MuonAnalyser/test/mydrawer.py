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
    print filedir
    hist.Scale(1/nevents)

def draw(h_init, y_name, hlists, name, text):
    #Plot style
    setMarkerStyle(hlists[0], 4, 20) #blue, circle
    setMarkerStyle(hlists[1], 1, 34) #black, cross
    setMarkerStyle(hlists[2], 2, 21) #red, square
    #setMarkerStyle(hlists[3], 3, 31) #green, patrol

    #Set canvas
    #canv = makeCanvas(plotvar+name, False)
    #setMargins(canv, False)
    canv = ROOT.TCanvas()
    canv.SetGrid()
    h_init.GetYaxis().SetTitle(y_name)
    h_init.Draw()
    drawSampleName(text)

    #Legend and drawing
    leg = ROOT.TLegend(0.25,0.72,0.45,0.86)
    for h in hlists:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"p")
    hlists[0].Draw("e1same")
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
    canv.SaveAs("%s_tmvaLoose_%s.png"%(plotvar,name))


datadir = './'
#datadir = '/xrootd/store/user/tt8888tt/muon/9_1_1/'
filenames = ["zmm.root", "zmm140.root", "zmm200.root"]
#filenames = ["zmm.root", "zmm140.root", "zmm200.root", "run2.root"]
#filenames = ["relval_"+x for x in filenames]

idname = sys.argv[1]
muonid = "is"+idname

#tmva
if "Tight" in idname: muonid = "tmva_bdt>-0.1472"
elif "Loose" in idname: muonid = "tmva_bdt>-0.167"
else: print idname, "Tight" in idname, "ID name has to be Tight or Loose."

binning_l = [[20,5,105],[24,0,2.4],[24,-2.4,2.4],[30,-3,3],[40,60,260],[20,0,10],[40,60,260]]
rangecut = "muon.Pt()>5&&abs(muon.Eta())<2.4"
    
#Set extra text
samplename = "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}"
text = samplename+", p_{T} > 5 GeV, |#eta| < 2.4"

for i, plotvar in enumerate(["muon.Pt()", "abs(muon.Eta())", "muon.Eta()", "muon.Phi()", "nvertex", "pu_density/2", "pu_numInteractions"]):
    #Efficiency
    hl_eff = []
    hl_eff.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "PhaseII PU 0",   binning_l[i], plotvar, rangecut, "%s&&muon_%s"%(rangecut,muonid)))
    hl_eff.append(getEff(datadir+filenames[1], "MuonAnalyser/gen", "PhaseII PU 140", binning_l[i], plotvar, rangecut, "%s&&muon_%s"%(rangecut,muonid)))
    hl_eff.append(getEff(datadir+filenames[2], "MuonAnalyser/gen", "PhaseII PU 200", binning_l[i], plotvar, rangecut, "%s&&muon_%s"%(rangecut,muonid)))
    #hl_eff.append(getEff(datadir+filenames[3], "MuonAnalyser/gen", "Run 2", binning_l[i], plotvar, rangecut, "%s&&muon_is%s"%(rangecut,muonid)))

    #Backgorund rate
    hl_bkg = []
    if "()" in plotvar:
        hl_bkg.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "PhaseII PU 0",   binning_l[i], plotvar, "%s&&!muon_signal&&muon_%s"%(rangecut,muonid)))
        hl_bkg.append(makeTH1(datadir+filenames[1], "MuonAnalyser/reco", "PhaseII PU 140", binning_l[i], plotvar, "%s&&!muon_signal&&muon_%s"%(rangecut,muonid)))
        hl_bkg.append(makeTH1(datadir+filenames[2], "MuonAnalyser/reco", "PhaseII PU 200", binning_l[i], plotvar, "%s&&!muon_signal&&muon_%s"%(rangecut,muonid)))
        #hl_bkg.append(makeTH1(datadir+filenames[3], "MuonAnalyser/reco", "Run 2", binning_l[i], plotvar, "%s&&!muon_signal&&muon_is%s"%(rangecut,muonid)))
        for j in range(len(hl_bkg)):
            divByNevents(datadir+filenames[j], hl_bkg[j])
            avg = 0
            for k in range(1,binning_l[i][0]+1):
                avg += hl_bkg[j].GetBinContent(k)
            avg = avg/float(binning_l[i][0])
            print plotvar, avg

    else:
        hl_bkg = []
        hl_bkg.append(getEff(datadir+filenames[0], "MuonAnalyser/reco", "PhaseII PU0",   binning_l[i], plotvar, "muon_no==1", "%s&&!muon_signal&&muon_%s"%(rangecut,muonid)))
        hl_bkg.append(getEff(datadir+filenames[1], "MuonAnalyser/reco", "PhaseII PU140", binning_l[i], plotvar, "muon_no==1", "%s&&!muon_signal&&muon_%s"%(rangecut,muonid)))
        hl_bkg.append(getEff(datadir+filenames[2], "MuonAnalyser/reco", "PhaseII PU200", binning_l[i], plotvar, "muon_no==1", "%s&&!muon_signal&&muon_%s"%(rangecut,muonid)))
        #hl_bkg.append(getEff(datadir+filenames[3], "MuonAnalyser/reco", "Run2", binning_l[i], plotvar, "muon_no==1", "%s&&!muon_signal&&muon_is%s"%(rangecut,muonid)))
    if "density" in plotvar: plotvar = plotvar.split('/')[0]

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
    y_name = idname+" Muon "
    if "Custom" in idname: y_name = "Tight Muon "

    h_init.SetMaximum(1.1)
    h_init.SetMinimum(0.8)
    draw(h_init, y_name+"Efficiency", hl_eff, "eff", text)

    if "()" in plotvar: h_init2.SetMaximum(max(h.GetMaximum() for h in hl_bkg)*2)
    else:
        maxY = hl_bkg[2].GetEfficiency(0)
        for j in range(binning_l[i][0]):
            x = hl_bkg[2].GetEfficiency(j)
            if maxY < x: maxY = x
        h_init2.SetMaximum(maxY*2)
    draw(h_init2, "Average "+y_name+"Bkg Multiplicity", hl_bkg, "bkgrate", text)

