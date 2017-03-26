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


def getRefSigEff(fileSig,treename,cutSig,cutIso):
    binning = [3000,0,30]
    hSigDen = makeTH1(fileSig,treename,"RefSigEffDen",binning,"muon_PFIsolation04",cutSig)
    hSigNum = makeTH1(fileSig,treename,"RefSigEffNum",binning,"muon_PFIsolation04",cutSig+" && "+cutIso)
    
    fMSigDen = hSigDen.Integral(0, binning[ 0 ] + 1)
    fMSigNum = hSigNum.Integral(0, binning[ 0 ] + 1)
    
    return fMSigNum / fMSigDen


def getRefBkgRej(fileBkg,treename,cutBkg,cutIso):
    binning = [3000,0,30]
    hBkgDen = makeTH1(fileBkg,treename,"RefBkgRejDen",binning,"muon_PFIsolation04",cutBkg)
    hBkgNum = makeTH1(fileBkg,treename,"RefBkgRejNum",binning,"muon_PFIsolation04",cutBkg+" && "+cutIso)
    
    fMBkgDen = hBkgDen.Integral(0, binning[ 0 ] + 1)
    fMBkgNum = hBkgNum.Integral(0, binning[ 0 ] + 1)
    
    return 1.0 - fMBkgNum / fMBkgDen


def getROC(fileSig,fileBkg,treename,title,binning,plotvar,cutSig,cutBkg,dicRef):
    hSig = makeTH1(fileSig,treename,title,binning,plotvar,cutSig)
    hBkg = makeTH1(fileBkg,treename,title,binning,plotvar,cutBkg)
    
    arrSigEff = []
    arrBkgRej = []
    
    fMSig = hSig.Integral(0, binning[ 0 ] + 1)
    fMBkg = hBkg.Integral(0, binning[ 0 ] + 1)
    
    fSigEffPrev = 1.0
    fBkgRejPrev = 1.0
    
    nPerf = 0.0
    
    for i in range(binning[ 0 ] + 2): 
      # DO NOT confuse this : 0th bin has all underflows, while (binning[0] + 1)-th bin has all overflows
      # Well, we have no underflows, but... for habit for safety.
      nX = i + 0
      
      fSigEff = hSig.Integral(0, nX) / fMSig
      fBkgRej = 1.0 - hBkg.Integral(0, nX) / fMBkg
      
      if "sigeff" in dicRef: 
          if fSigEffPrev <= dicRef[ "sigeff" ] and dicRef[ "sigeff" ] <= fSigEff: 
            nPerf = i
      elif "bkgrej" in dicRef: 
          if fBkgRej <= dicRef[ "bkgrej" ] and dicRef[ "bkgrej" ] <= fBkgRejPrev: 
            nPerf = i
      
      arrSigEff.append(fSigEff)
      arrBkgRej.append(fBkgRej)
      
      fSigEffPrev = fSigEff
      fBkgRejPrev = fBkgRej
    
    arrX = array.array("d", arrSigEff)
    arrY = array.array("d", arrBkgRej)
    
    graphROC = ROOT.TGraph(binning[ 0 ] + 2, arrX, arrY)
    graphROC.SetTitle(title)
    
    return {"graph": copy.deepcopy(graphROC), 
        "performance": binning[1] + nPerf * ( binning[2] - binning[1] ) / binning[0]}


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
    {"name": "PUPPIWL", "plotvar": "muon_puppiIsoWithLep",    "title": "PUPPI - with lepton, R = 0.4", 
        "color": 4, "shape": 20}, # blue,  filled circle
    {"name": "PUPPINL", "plotvar": "muon_puppiIsoWithoutLep", "title": "PUPPI - without lepton, R = 0.4", 
        "color": 2, "shape": 21}, # red,   filled square
    #{"name": "PUPPICB", "plotvar": "muon_puppiIsoCombined",   "title": "PUPPI - combined (ratio : 0.5)", 
    #    "color": 3, "shape": 34}, # green, filled cross
    {"name": "PUPPINEWWL", "plotvar": "muon_puppiIso",      "title": "new PUPPI - with lepton, R = 0.4", 
        "color": 3, "shape": 34}, # green, filled cross
    {"name": "PUPPINEWNL", "plotvar": "muon_puppiIsoNoLep", "title": "new PUPPI - without lepton, R = 0.4", 
        "color": 1, "shape": 21}, # black, star
    
    {"name": "PF", "plotvar": "muon_PFIsolation04",      "title": "PF Isolation, R = 0.4", 
        "color": 6, "shape": 25}, # pink,  unfilled square
]

"""
dicSampleType = {
    "PU0":   {"title": "PU 0",   "file_signal": "puppi_PU0.root"}, 
    "PU200": {"title": "PU 200", "file_signal": "puppi_PU200.root"}, 
}
"""
#strPathSamp = "/cms/scratch/quark2930/Work/muon_upgrade/samples/"
#strPathSamp = "/cms/scratch/gwheo/muonPerf_900_pre6/src/MuonPerformance/MuonAnalyser/test/"
strPathSamp = "/xrootd/store/user/tt8888tt/muon/"
# rereco
"""
dicSampleType = {
    "PU0":   {"title": "PU 0",   
        "file_sig": strPathSamp + "run_ZMM_PU0_pre4_rereco02_ver01.root", 
        "file_bkg": strPathSamp + "run_QCD_PU0_pre4_rereco02_ver01.root"}, 
    "PU140": {"title": "PU 140", 
        "file_sig": strPathSamp + "run_ZMM_PU140_pre4_rereco02_ver01.root", 
        "file_bkg": strPathSamp + "run_QCD_PU140_pre4_rereco02_ver01.root"}, 
    "PU200": {"title": "PU 200", 
        "file_sig": strPathSamp + "run_ZMM_PU200_pre4_rereco02_ver01.root", 
        "file_bkg": strPathSamp + "run_QCD_PU200_pre4_rereco02_ver01.root"}, 
}
"""
dicSampleType = {
    "PU0":   {"title": "PU 0",   
        "file_sig": strPathSamp + "zmm.root", 
        "file_bkg": strPathSamp + "qcd.root"}, 
    "PU140": {"title": "PU 140", 
        "file_sig": strPathSamp + "zmm140.root", 
        "file_bkg": strPathSamp + "qcd140.root"}, 
    "PU200": {"title": "PU 200", 
        "file_sig": strPathSamp + "zmm200.root", 
        "file_bkg": strPathSamp + "qcd200.root"}, 
}
# RelVal
if sys.argv[ 2 ] == "3" or sys.argv[ 2 ] == "4": 
  dicSampleType = {
      "PU0":   {"title": "PU 0",   
          "file_sig": strPathSamp + "run_ZMM_PU0_pre4_ver02.root", 
          "file_bkg": strPathSamp + "run_QCD_PU0_pre4_ver02.root"}, 
      "PU140": {"title": "PU 140", 
          "file_sig": strPathSamp + "run_ZMM_PU140_pre4_ver02.root", 
          "file_bkg": strPathSamp + "run_QCD_PU140_pre4_ver02.root"}, 
      "PU200": {"title": "PU 200", 
          "file_sig": strPathSamp + "run_ZMM_PU200_pre4_ver02.root", 
          "file_bkg": strPathSamp + "run_QCD_PU200_pre4_ver02.root"}, 
  }

arrListID = ["Loose", "Tight", "LooseMod", "TightModNoIP"]


strTypePU = sys.argv[1]
id = sys.argv[2]

if strTypePU not in dicSampleType: 
    print "Error : " + strTypePU + " is not in option: "
    print dicSampleType.keys()
    exit(1)

if 1 == 0 and id not in arrListID: 
    print "Error : " + id + " is not in option: "
    print arrListID
    exit(1)

strPUTitle = dicSampleType[ strTypePU ][ "title" ]

strSampleSig = dicSampleType[ strTypePU ][ "file_sig" ]
strSampleBkg = dicSampleType[ strTypePU ][ "file_bkg" ]

strCutPT  = "15"
strCutEta = "2.4"

strDPV = "0.5"

#id = "Tight"
id = sys.argv[ 3 ]

strCutDef = ""

#strCutRecNor = "recoMuon.Pt() > 5 && recoMuon_isMuon"
#strCutRecNor = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4 && recoMuon_is%s"%id
#strCutDef = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4"
#strCutDef = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4 && recoMuon_is%s"%id
if id in arrListID: 
    strCutDef = "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s && muon_is%(ID)s && abs(muon_poszPV0 - muon_poszSimPV) < %(dPV)s"%{"pT":strCutPT, "Eta": strCutEta, "ID":id, "dPV": strDPV}
    if sys.argv[ 2 ] == "2" or sys.argv[ 2 ] == "4": 
        strCutDef = "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s && muon_is%(ID)s"%{"pT":strCutPT, "Eta": strCutEta, "ID":id, "dPV": strDPV}
else: 
    strCutDef = "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s && abs(muon_poszPV0 - muon_poszSimPV) < %(dPV)s"%{"pT":strCutPT, "Eta": strCutEta, "ID":id, "dPV": strDPV}
    if sys.argv[ 2 ] == "2" or sys.argv[ 2 ] == "4": 
        strCutDef = "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s"%{"pT":strCutPT, "Eta": strCutEta, "ID":id, "dPV": strDPV}
strCutSig = strCutDef + " && muon_signal"
strCutBkg = strCutDef


strCutIsoPF = ""
fSigEff = 0.0

if "Tight" in id: 
    strCutIsoPF = "muon_PFIsolation04 < 0.15"
    fSigEff = 0.95
else: 
    strCutIsoPF = "muon_PFIsolation04 < 0.25"
    fSigEff = 0.98

#fRefBkgRej = getRefBkgRej(strSampleBkg, "MuonAnalyser/reco", 
#    strCutDef, "recoMuon_PFIsolation04 < 0.15")
dicRef = {"sigeff": fSigEff}

dicResCut = {}

strDataType = ""
strVtxCut   = ""

if sys.argv[ 2 ] == "1": strDataType, strVtxCut = "rereco", "withvtxcut"
if sys.argv[ 2 ] == "2": strDataType, strVtxCut = "rereco", "novtxcut"
if sys.argv[ 2 ] == "3": strDataType, strVtxCut = "relval", "withvtxcut"
if sys.argv[ 2 ] == "4": strDataType, strVtxCut = "relval", "novtxcut"

strPrefix = strDataType + "_" + strVtxCut

for dicPlotvar in arrPlotvar: 
    plotvar = dicPlotvar[ "plotvar" ]
    
    dicRes = getROC(strSampleSig, strSampleBkg, "MuonAnalyser/reco", 
        dicPlotvar[ "title" ], binMain, plotvar, strCutSig, strCutBkg, dicRef)
    
    dicPlotvar[ "graph" ] = dicRes[ "graph" ]
    
    setMarkerStyle(dicPlotvar[ "graph" ], dicPlotvar[ "color" ], dicPlotvar[ "shape" ])
    
    dicPlotvar[ "graph" ].GetXaxis().SetLimits(0.0, 1.1)
    dicPlotvar[ "graph" ].SetMaximum(1.1)
    
    #print "%s; %s; %s; %s; %0.5f"%(id, strTypePU, sys.argv[ 2 ], dicPlotvar[ "title" ], dicRes[ "performance" ])
    print "%s; %s; %s; %s%s; %0.5f"%(strVtxCut, strDataType, dicPlotvar[ "name" ], id, strTypePU, dicRes[ "performance" ])
    
    #dicNewRes = {}
    
    #dicNewRes[ "" ]
    
    #dicResCut[ dicPlotvar[ "name" ] ] = dicNewRes

#json.dump(dicResCut, open("isocut_%s_%s_%s.json"%(strPUTitle.replace(" ", ""), id, strPrefix), "w"))

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

strAddTitleVtx = ""

if strVtxCut == "withvtxcut": 
    strAddTitleVtx = ", |z_{reco} - z_{sim}| < 0.5"

if id in arrListID:
    drawSampleName(("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} (%(PU)s) and QCD events\n"
        "p_{T} > %(pT)s GeV, |#eta| < %(Eta)s, %(ID)s Muon%(VtxCut)s")%{"PU": strPUTitle, "pT": strCutPT, "Eta": strCutEta, "ID": id, "VtxCut": strAddTitleVtx})
else:
    drawSampleName(("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} (%(PU)s) and QCD events\n"
        "p_{T} > %(pT)s GeV, |#eta| < %(Eta)s%(VtxCut)s")%{"PU": strPUTitle, "pT": strCutPT, "Eta": strCutEta, "ID": id, "VtxCut": strAddTitleVtx})

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
#canv.SaveAs("ROCCurves_%s_%s_dPV%s.png"%(strPUTitle.replace(" ", ""), id, strDPV))
canv.SaveAs("ROCCurves_%s_%s_%s.png"%(strPUTitle.replace(" ", ""), id, strPrefix))
    

