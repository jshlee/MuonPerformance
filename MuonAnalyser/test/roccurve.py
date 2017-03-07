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


def getROC(fileSig,fileBkg,treename,title,binning,plotvar,cutSig,cutBkg):
    hSig = makeTH1(fileSig,treename,title,binning,plotvar,cutSig)
    hBkg = makeTH1(fileBkg,treename,title,binning,plotvar,cutBkg)
    
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
    
    graphROC = ROOT.TGraph(binning[ 0 ] + 2, arrX, arrY)
    graphROC.SetTitle(title)
    
    return copy.deepcopy(graphROC)


def drawSampleName(samplename):
    fX = 0.22
    fY = 0.58
    fSizeTex = 0.038
    
    tex2 = ROOT.TLatex()
    
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(fSizeTex)
    
    for i, strLine in enumerate(samplename.split("\n")): 
        tex2.DrawLatex(fX, fY - i * 1.1 * fSizeTex, strLine)


datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'

binMain = [1500, 0, 6.0]

arrPlotvar = [
    {"plotvar": "recoMuon_puppiIsoWithLep",    "title": "PUPPI - with lepton, R = 0.4", 
        "color": 4, "shape": 20}, # blue,  filled circle
    {"plotvar": "recoMuon_puppiIsoWithoutLep", "title": "PUPPI - without lepton, R = 0.4", 
        "color": 2, "shape": 21}, # red,   filled square
    {"plotvar": "recoMuon_puppiIsoCombined",   "title": "PUPPI - combined (ratio : 0.5)", 
        "color": 3, "shape": 34}, # green, filled cross
    
    {"plotvar": "recoMuon_PFIsolation04",      "title": "PF Isolation, R = 0.4", 
        "color": 6, "shape": 25}, # pink,  unfilled square
]
"""
arrPlotvar = [
    {"plotvar": "recoMuon_puppiIsoWithLep",    "title": "PUPPI - with lepton, R = 0.3", 
        "color": 4, "shape": 20}, # blue,  filled circle
    {"plotvar": "recoMuon_puppiIsoWithoutLep", "title": "PUPPI - without lepton, R = 0.3", 
        "color": 2, "shape": 21}, # red,   filled square
    {"plotvar": "recoMuon_puppiIsoCombined",   "title": "PUPPI - combined (ratio : 0.5)", 
        "color": 3, "shape": 34}, # green, filled cross
    
    #{"plotvar": "recoMuon_TrkIsolation03",     "title": "Track Isolation, R = 0.3", 
    #    "color": 1, "shape": 24}, # black, unfilled circle
    {"plotvar": "recoMuon_PFIsolation03",      "title": "PF Isolation, R = 0.3", 
        "color": 6, "shape": 25}, # pink,  unfilled square
]
"""
"""
arrPlotvar = [
    {"plotvar": "recoMuon_puppiIsoWithLep05",    "title": "PUPPI - with lepton, R = 0.5", 
        "color": 4, "shape": 20}, # blue,  filled circle
    {"plotvar": "recoMuon_puppiIsoWithoutLep05", "title": "PUPPI - without lepton, R = 0.5", 
        "color": 2, "shape": 21}, # red,   filled square
    {"plotvar": "recoMuon_puppiIsoCombined05",   "title": "PUPPI - combined (ratio : 0.5)", 
        "color": 3, "shape": 34}, # green, filled cross
    
    {"plotvar": "recoMuon_TrkIsolation05",     "title": "Track Isolation, R = 0.5", 
        "color": 1, "shape": 24}, # black, unfilled circle
]
"""

"""
dicSampleType = {
    "PU0":   {"title": "PU 0",   "file_signal": "puppi_PU0.root"}, 
    "PU200": {"title": "PU 200", "file_signal": "puppi_PU200.root"}, 
}
"""
strPathSamp = "/cms/scratch/quark2930/Work/muon_upgrade/samples/"
dicSampleType = {
    "PU0":   {"title": "PU 0",   
        "file_sig": strPathSamp + "run_ZMM_PU0_pre4_rereco01_ver01.root", 
        "file_bkg": strPathSamp + "run_QCD_PU0_pre4_rereco01_ver01.root"}, 
    "PU140": {"title": "PU 140", 
        "file_sig": strPathSamp + "run_ZMM_PU140_pre4_rereco01_ver01.root", 
        "file_bkg": strPathSamp + "run_QCD_PU140_pre4_rereco01_ver01.root"}, 
    "PU200": {"title": "PU 200", 
        "file_sig": strPathSamp + "run_ZMM_PU200_pre4_rereco01_ver01.root", 
        "file_bkg": strPathSamp + "run_QCD_PU200_pre4_rereco01_ver01.root"}, 
}

arrListID = ["Tight", "Loose"]


strTypePU = sys.argv[1]
id = sys.argv[2]

if strTypePU not in dicSampleType: 
    print "Error : " + strTypePU + " is not in option: "
    print dicSampleType.keys()
    exit(1)

if id not in arrListID: 
    print "Error : " + id + " is not in option: "
    print arrListID
    exit(1)

strPUTitle = dicSampleType[ strTypePU ][ "title" ]

strSampleSig = dicSampleType[ strTypePU ][ "file_sig" ]
strSampleBkg = dicSampleType[ strTypePU ][ "file_bkg" ]

strCutPT  = "15"
strCutEta = "2.4"

strDPV = "0.5"

#strCutRecNor = "recoMuon.Pt() > 5 && recoMuon_isMuon"
#strCutRecNor = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4 && recoMuon_is%s"%id
#strCutDef = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4"
#strCutDef = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4 && recoMuon_is%s"%id
strCutDef = "recoMuon.Pt() > %(pT)s && abs(recoMuon.Eta()) < %(Eta)s && recoMuon_is%(ID)s && abs(recoMuon_poszPV0 - recoMuon_poszSimPV) < %(dPV)s"%{"pT":strCutPT, "Eta": strCutEta, "ID":id, "dPV": strDPV}
strCutSig = strCutDef + " && recoMuon_signal"
strCutBkg = strCutDef


for dicPlotvar in arrPlotvar: 
    plotvar = dicPlotvar[ "plotvar" ]
    
    dicPlotvar[ "graph" ] = getROC(strSampleSig, strSampleBkg, "MuonAnalyser/reco", 
        dicPlotvar[ "title" ], binMain, plotvar, strCutSig, strCutBkg)
    
    setMarkerStyle(dicPlotvar[ "graph" ], dicPlotvar[ "color" ], dicPlotvar[ "shape" ])
    
    dicPlotvar[ "graph" ].GetXaxis().SetLimits(0.0, 1.1)
    dicPlotvar[ "graph" ].SetMaximum(1.1)
    
    print "%s is done"%dicPlotvar[ "title" ]

#Set canvas
canv = makeCanvas("canv1", False)
setMargins(canv, False)

#Legend and drawing
leg = ROOT.TLegend(0.18,0.2,0.45,0.40)

x_name = "Signal Efficiency"
y_name = "Background Rejection"

for i, dicPlotvar in enumerate(arrPlotvar): 
    if i == 0: 
        dicPlotvar[ "graph" ].GetXaxis().SetTitle(x_name)
        dicPlotvar[ "graph" ].GetYaxis().SetTitle(y_name)
        dicPlotvar[ "graph" ].GetYaxis().SetTitleOffset(0.95)
        
        dicPlotvar[ "graph" ].Draw("")
    else :
        dicPlotvar[ "graph" ].Draw("same")

    leg.AddEntry(dicPlotvar[ "graph" ], dicPlotvar[ "graph" ].GetTitle(), "pl")

drawSampleName(("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} (%(PU)s) and QCD events\n"
    "p_{T} > %(pT)s GeV, |#eta| < %(Eta)s, %(ID)s Muon")%{"PU": strPUTitle, "pT": strCutPT, "Eta": strCutEta, "ID": id, "dPV": strDPV})

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
canv.SaveAs("roccurves_%s_%s_dPV%s.png"%(strPUTitle.replace(" ", ""), id, strDPV))
    

