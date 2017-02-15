import ROOT, copy, os, sys
import array
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(0.2)
    h.SetLineColor(color)
    h.SetLineWidth(2)


def getROC(fileSig,fileBkg,treename,title,binning,plotvar,cut):
    hSig = makeTH1(fileSig,treename,title,binning,plotvar,cut)
    hBkg = makeTH1(fileBkg,treename,title,binning,plotvar,cut)
    
    arrSigEff = []
    arrBkgRej = []
    
    fMSig = hSig.Integral(0, binning[ 0 ] + 1)
    fMBkg = hBkg.Integral(0, binning[ 0 ] + 1)
    
    for i in range(binning[ 0 ] + 2): 
      # DO NOT confuse this : 0th bin has all underflows, while (binning[0] + 1)-th bin has all overflows
      # Well, we have no underflows, but... for habit for safety.
      nX = i + 0
      
      fSigEff = hSig.Integral(0, nX) / fMSig
      fBkgRej = 1.0 - hBkg.Integral(0, nX) / fMBkg
      
      arrSigEff.append(fSigEff)
      arrBkgRej.append(fBkgRej)
    
    arrX = array.array("d", arrSigEff)
    arrY = array.array("d", arrBkgRej)
    
    graphROC = ROOT.TGraph(binning[ 0 ] + 1, arrX, arrY)
    graphROC.SetTitle(title)
    
    return copy.deepcopy(graphROC)


def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.038)
    tex2.DrawLatex(0.15, 0.68, samplename)


datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
#datadir = "TenMuExtendedE_"
id = sys.argv[1]

binMain = [500, 0, 2.0]

arrPlotvar = [
    {"plotvar": "recoMuon_puppiIsoWithLep",    "title": "PUPPI - with lepton, R = 0.4", 
        "color": 4, "shape": 20}, # blue,  filled circle
    {"plotvar": "recoMuon_puppiIsoWithoutLep", "title": "PUPPI - without lepton, R = 0.4", 
        "color": 2, "shape": 21}, # red,   filled square
    {"plotvar": "recoMuon_puppiIsoCombined",   "title": "PUPPI - combined (ratio : 0.5)", 
        "color": 3, "shape": 34}, # green, filled cross
    
    {"plotvar": "recoMuon_TrkIsolation03",     "title": "Track Isolation, R = 0.3", 
        "color": 1, "shape": 24}, # black, unfilled circle
    {"plotvar": "recoMuon_PFIsolation04",      "title": "PF Isolation, R = 0.4", 
        "color": 6, "shape": 25}, # pink,  unfilled square
]

"""
arrSampleType = [
    {"title": "Phase II PU0",   "filename": "out_PU0.root",   "id": "Tight", "color": 4, "shape": 20}, # blue,  filled circle
    {"title": "Phase II PU200", "filename": "out_PU200.root", "id": "Tight", "color": 2, "shape": 21}, # red,   filled square
    {"title": "Phase II QCD",   "filename": "out_QCD.root",   "id": "Tight", "color": 1, "shape": 34}, # black, filled cross
    
    {"title": "Phase II PU0",   "filename": "out_PU0.root",   "id": "Loose", "color": 3, "shape": 24}, # green, unfilled circle
    {"title": "Phase II PU200", "filename": "out_PU200.root", "id": "Loose", "color": 6, "shape": 25}, # pink,  unfilled square
    {"title": "Phase II QCD",   "filename": "out_QCD.root",   "id": "Loose", "color": 7, "shape": 28}, # cyan,  unfilled cross
]

arrSampleType = [
    {"title": "Phase II PU0",   "filename": "puppi_PU0.root",   "id": "Tight", "color": 4, "shape": 20}, # blue,  filled circle
    #{"title": "Phase II PU200", "filename": "puppi_PU200.root", "id": "Tight", "color": 2, "shape": 21}, # red,   filled square
    {"title": "Phase II QCD",   "filename": "puppi_QCD.root",   "id": "Tight", "color": 1, "shape": 34}, # black, filled cross
    
    {"title": "Phase II PU0",   "filename": "puppi_PU0.root",   "id": "Loose", "color": 3, "shape": 24}, # green, unfilled circle
    #{"title": "Phase II PU200", "filename": "puppi_PU200.root", "id": "Loose", "color": 6, "shape": 25}, # pink,  unfilled square
    {"title": "Phase II QCD",   "filename": "puppi_QCD.root",   "id": "Loose", "color": 7, "shape": 28}, # cyan,  unfilled cross
]
"""
strSampleSig = "puppi_PU200.root"
strSampleBkg = "puppi_QCD.root"

strPUTitle = "PU 200"

#strCutRecNor = "recoMuon.Pt() > 5 && recoMuon_isMuon"
strCutRecNor = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4"


for dicPlotvar in arrPlotvar: 
    plotvar = dicPlotvar[ "plotvar" ]
    
    dicPlotvar[ "graph" ] = getROC(strSampleSig, strSampleBkg, "MuonAnalyser/reco", 
        dicPlotvar[ "title" ], binMain, plotvar, strCutRecNor)
    
    setMarkerStyle(dicPlotvar[ "graph" ], dicPlotvar[ "color" ], dicPlotvar[ "shape" ])
    
    dicPlotvar[ "graph" ].GetXaxis().SetLimits(0.0, 1.1)
    dicPlotvar[ "graph" ].SetMaximum(1.1)

#Set canvas
canv = makeCanvas("canv1", False)
setMargins(canv, False)

#Legend and drawing
leg = ROOT.TLegend(0.18,0.2,0.45,0.40)

x_name = "Signal Efficiency"
y_name = "Background Rejection"

for i, dicPlotvar in enumerate(arrPlotvar): 
    #if "Trk" not in dicPlotvar[ "plotvar" ]: continue
    #i = 0
    if i == 0: 
        dicPlotvar[ "graph" ].GetXaxis().SetTitle(x_name)
        dicPlotvar[ "graph" ].GetYaxis().SetTitle(y_name)
        dicPlotvar[ "graph" ].GetYaxis().SetTitleOffset(0.95)
        
        dicPlotvar[ "graph" ].Draw("")
    else :
        dicPlotvar[ "graph" ].Draw("same")
    leg.AddEntry(dicPlotvar[ "graph" ], dicPlotvar[ "graph" ].GetTitle(), "pl")

drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} (%s) and QCD events, \np_{T} > 5 GeV, |#eta| < 2.4"%strPUTitle)

leg.SetTextFont(61)
leg.SetTextSize(0.04)
leg.SetBorderSize(0)
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
canv.SaveAs("roccurves_%s.png"%id)
    

