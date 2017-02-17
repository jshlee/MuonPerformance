import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(0.5)
    h.SetLineColor(color)
    h.SetLineWidth(2)


def getH1_Normalized(filename,treename,title,binning,plotvar,cut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,cut)
    
    normfactor = h1.GetEntries()
    if normfactor > 0.0: h1.Scale(1.0 / normfactor)
    
    h1.SetTitle(title)
    
    return copy.deepcopy(h1)


def getEff(filename,treename,title,binning,plotvar,dencut,numcut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,dencut)
    h2 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    hEff = ROOT.TEfficiency(h2,h1)
    hEff.SetTitle(title)
    return copy.deepcopy(hEff)


def drawSampleName(samplename):
    fX = 0.18
    fY = 0.85
    fSizeTex = 0.038
    
    tex2 = ROOT.TLatex()
    
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(fSizeTex)
    
    for i, strLine in enumerate(samplename.split("\n")): 
        tex2.DrawLatex(fX, fY - i * 1.1 * fSizeTex, strLine)


datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
#datadir = "TenMuExtendedE_"
isotype = sys.argv[1]

listIso = ["Trk", "PF"]

if isotype not in listIso:
    print "Error : Isolation condition (second arg) should be Trk or PF"
    exit(1)

fMaxPuppi = 0.5
nNumBinPuppi = 100

"""
arrPlotvar = [
    #{"plotvar": "genMuon.Pt()",        "binning": [10,5,105], "effrate": False, "id": "Tight", 
    #    "xtitle": "p_{T} (GeV)", "ytitle": "Efficiency"}, 
    #{"plotvar": "abs(genMuon.Eta())",  "binning": [15,0,2.4], "effrate": False, "id": "Tight", 
    #    "xtitle": "|#eta|",      "ytitle": "Efficiency"}, 
    #{"plotvar": "genMuon.Phi()",       "binning": [12,-3,3],  "effrate": False, "id": "Tight", 
    #    "xtitle": "#phi",        "ytitle": "Efficiency"}, 
]
"""

"""
    {"plotvar": "genMuon.Pt()",        "binning": [10,5,105], "effrate": True, "id": "Tight", 
        "xtitle": "p_{T} (GeV)", "ytitle": "Efficiency", "min": 0.6, "max": 1.1}, 
    {"plotvar": "abs(genMuon.Eta())",  "binning": [15,0,2.4], "effrate": True, "id": "Tight", 
        "xtitle": "|#eta|",      "ytitle": "Efficiency", "min": 0.6, "max": 1.1}, 
    {"plotvar": "genMuon.Phi()",       "binning": [12,-3,3],  "effrate": True, "id": "Tight", 
        "xtitle": "#phi",        "ytitle": "Efficiency", "min": 0.6, "max": 1.1}, 
"""
    
"""
    {"plotvar": "recoMuon.Pt()",       "binning": [10,5,105], "effrate": True, "id": "Tight", 
        "xtitle": "p_{T} (GeV)", "ytitle": "Fake Rate", "max": 1.2}, 
    {"plotvar": "abs(recoMuon.Eta())", "binning": [8,0,2.4],  "effrate": True, "id": "Tight", 
        "xtitle": "|#eta|",      "ytitle": "Fake Rate", "max": 1.2}, 
    {"plotvar": "recoMuon.Phi()",      "binning": [12,-3,3],  "effrate": True, "id": "Tight", 
        "xtitle": "#phi",        "ytitle": "Fake Rate", "max": 1.2}, 
    
    {"plotvar": "recoMuon.Pt()",       "binning": [10,5,105], "effrate": True, "id": "Loose", 
        "xtitle": "p_{T} (GeV)", "ytitle": "Fake Rate", "max": 1.2}, 
    {"plotvar": "abs(recoMuon.Eta())", "binning": [8,0,2.4],  "effrate": True, "id": "Loose", 
        "xtitle": "|#eta|",      "ytitle": "Fake Rate", "max": 1.2}, 
    {"plotvar": "recoMuon.Phi()",      "binning": [12,-3,3],  "effrate": True, "id": "Loose", 
        "xtitle": "#phi",        "ytitle": "Fake Rate", "max": 1.2}, 
"""
    
arrPlotvar = [
    {"plotvar": "recoMuon_TrkIsolation03", "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Tight", "type": "all", 
        "ylog": True, "xtitle": "Track Isolation", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
]
"""
    {"plotvar": "recoMuon_PFIsolation04",  "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Tight", "type": "all",
        "ylog": True, "xtitle": "PF Isolation", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    
    {"plotvar": "recoMuon_TrkIsolation03", "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Loose", "type": "all",
        "ylog": True, "xtitle": "Track Isolation", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    {"plotvar": "recoMuon_PFIsolation04",  "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Loose", "type": "all",
        "ylog": True, "xtitle": "PF Isolation", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    
    {"plotvar": "recoMuon_puppiIsoWithLep",    "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Tight", "type": "all",
        "ylog": True, "xtitle": "PUPPI Isolation with lepton", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    {"plotvar": "recoMuon_puppiIsoWithoutLep", "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Tight", "type": "all",
        "ylog": True, "xtitle": "PUPPI Isolation without lepton", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    {"plotvar": "recoMuon_puppiIsoCombined",   "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Tight", "type": "all",
        "ylog": True, "xtitle": "Combined PUPPI Isolation (ratio = 0.5)", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    
    {"plotvar": "recoMuon_puppiIsoWithLep",    "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Loose", "type": "all",
        "ylog": True, "xtitle": "PUPPI Isolation with lepton", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    {"plotvar": "recoMuon_puppiIsoWithoutLep", "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Loose", "type": "all",
        "ylog": True, "xtitle": "PUPPI Isolation without lepton", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
    {"plotvar": "recoMuon_puppiIsoCombined",   "binning": [nNumBinPuppi,0,fMaxPuppi], "id": "Loose", "type": "all",
        "ylog": True, "xtitle": "Combined PUPPI Isolation (ratio = 0.5)", "ytitle": "# of events (normalized)", "min": 0.00005, "max":1.0}, 
]
"""

strFileSampZMM0   = "puppi_ZMM_PU0_pre4_fixed01_ver2.root"
strFileSampZMM140 = "puppi_ZMM_PU140_pre4_fixed01.root"
strFileSampQCD0   = "puppi_QCD_PU0_pre4_fixed01_ver2.root"
strFileSampQCD140 = "puppi_QCD_PU140_pre4_fixed01.root"

"""
arrSampleType = [
    {"title": "Phase II PU0 - pre2",   "filename": "puppi_ZMM_PU0_pre2.root",   "id": "Tight", 
        "color": 4, "shape": 20}, # blue,  filled circle
    {"title": "Phase II PU200 - pre2", "filename": "puppi_ZMM_PU200_pre2.root", "id": "Tight", 
        "color": 2, "shape": 21}, # red,   filled square
    {"title": "Phase II QCD - pre2",   "filename": "puppi_QCD_PU0_pre2.root",   "id": "Tight",
        "color": 1, "shape": 34}, # black, filled cross
    
    {"title": "Phase II PU0 - pre2",   "filename": "puppi_ZMM_PU0_pre2.root",   "id": "Loose",
        "color": 3, "shape": 24}, # green, unfilled circle
    {"title": "Phase II PU200 - pre2", "filename": "puppi_ZMM_PU200_pre2.root", "id": "Loose", 
        "color": 6, "shape": 25}, # pink,  unfilled square
    {"title": "Phase II QCD - pre2",   "filename": "puppi_QCD_PU0_pre2.root",   "id": "Loose",
        "color": 7, "shape": 28}, # cyan,  unfilled cross
];
"""
arrSampleType = [
    {"title": "Phase II Signal PU0 - pre4",   "filename": strFileSampZMM0,   "id": "Tight", 
        "color": 4, "shape": 20, "extracut": "recoMuon_signal"}, # blue,  filled circle
    {"title": "Phase II Signal PU140 - pre4", "filename": strFileSampZMM140, "id": "Tight", 
        "color": 2, "shape": 21, "extracut": "recoMuon_signal"}, # red,   filled square
    {"title": "Phase II QCD PU0 - pre4",      "filename": strFileSampQCD0,   "id": "Tight",
        "color": 1, "shape": 34, "extracut": "recoMuon_pdgId != 0"}, # black, filled cross
    {"title": "Phase II QCD PU140 - pre4",    "filename": strFileSampQCD140, "id": "Tight",
        "color": 7, "shape": 28, "extracut": "recoMuon_pdgId != 0"}, # cyan,  unfilled cross
    
    {"title": "Phase II Signal PU0 - pre4",   "filename": strFileSampZMM0,   "id": "Loose",
        "color": 3, "shape": 24, "extracut": "recoMuon_signal"}, # green, unfilled circle
    {"title": "Phase II Signal PU140 - pre4", "filename": strFileSampZMM140, "id": "Loose", 
        "color": 6, "shape": 25, "extracut": "recoMuon_signal"}, # pink,  unfilled square
    {"title": "Phase II QCD PU0 - pre4",      "filename": strFileSampQCD0,   "id": "Loose",
        "color": 1, "shape": 34, "extracut": "recoMuon_pdgId != 0"}, # black, filled cross
    {"title": "Phase II QCD PU140 - pre4",    "filename": strFileSampQCD140, "id": "Loose",
        "color": 7, "shape": 28, "extracut": "recoMuon_pdgId != 0"}, # cyan,  unfilled cross
];

listIDCfg = {
    "Tight": {"isoCut": {"Trk": "TrkIsolation03 < 0.05", "PF": "PFIsolation04 < 0.15"}}, 
    "Loose": {"isoCut": {"Trk": "TrkIsolation03 < 0.10", "PF": "PFIsolation04 < 0.25"}}, 
}

strCutPT  = "5"
strCutEta = "2.4"

strCutGenNor = "genMuon.Pt() > %(pT)s && abs(genMuon.Eta()) < %(Eta)s"%{"pT":strCutPT, "Eta": strCutEta}
strCutGenDen = strCutGenNor
#strCutGenNum = strCutGenDen + " && genMuon_is%(ID)s && genMuon_%(iso)s" # ID is held
strCutGenNum = strCutGenDen + " && genMuon_isMuon && genMuon_%(iso)s"

#strCutRecNor = "recoMuon.Pt() > 5 && recoMuon_is%(ID)s"
#strCutRecNor = "recoMuon.Pt() > 5"
#strCutRecNor = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4 && recoMuon_TrkIsolation03 <= 0.00001"
strCutRecNor = "recoMuon.Pt() > %(pT)s && abs(recoMuon.Eta()) < %(Eta)s"%{"pT":strCutPT, "Eta": strCutEta} + " && recoMuon_is%(ID)s"
#strCutRecNor = "recoMuon.Pt() > 20 && abs(recoMuon.Eta()) < 2"
strCutRecDen = strCutRecNor + " && recoMuon_isMuon && recoMuon_%(iso)s"
strCutRecNum = strCutRecDen + " && !recoMuon_signal"

strCutGenIso = ""
strCutRecIso = ""


for i, dicPlotvar in enumerate(arrPlotvar):
    plotvar = dicPlotvar[ "plotvar" ]
    binCurr = dicPlotvar[ "binning" ]
    
    strCutNor = ""
    strCutDen = ""
    strCutNum = ""
    
    strTree = ""
    
    if "genMuon" in plotvar: 
        strCutNor = strCutGenNor
        strCutDen = strCutGenDen
        strCutNum = strCutGenNum
        
        strTree = "MuonAnalyser/gen"
    elif "recoMuon" in plotvar: 
        strCutNor = strCutRecNor
        strCutDen = strCutRecDen
        strCutNum = strCutRecNum
        
        strTree = "MuonAnalyser/reco"
    
    arrHist = []
    
    nMax = 0
    nMin = 10000000
    
    #Get histos
    for sampHead in arrSampleType: 
        if "type" not in dicPlotvar or dicPlotvar[ "type" ] != "all": 
            if "gen" in plotvar and "QCD" in sampHead[ "title" ]: 
                continue
            if "reco" in plotvar and "QCD" not in sampHead[ "title" ]: 
                continue

        strIDCurr = sampHead[ "id" ]
        
        # This is the isolation cut condition
        strCutIso = listIDCfg[ strIDCurr ][ "isoCut" ][ isotype ]
        
        strCutExtra = ""
        if "extracut" in sampHead: strCutExtra = " && " + sampHead[ "extracut" ]
        
        #if "Iso" not in plotvar: 
        if "effrate" in dicPlotvar and dicPlotvar[ "effrate" ]: 
            # Drawing efficiency / fake rate plot
            sampHead[ "hist" ] = getEff(datadir + sampHead[ "filename" ], strTree, 
                sampHead[ "title" ] + " - " + strIDCurr, binCurr, plotvar, 
                strCutDen%{"ID": strIDCurr, "iso": strCutIso} + strCutExtra, # cut for denominator
                strCutNum%{"ID": strIDCurr, "iso": strCutIso} + strCutExtra) # cut for nominator
        else:
            if strIDCurr != dicPlotvar[ "id" ]:
                continue
            
            strCut = strCutNor%{"ID": dicPlotvar[ "id" ], "iso": strCutIso} + strCutExtra # cut for normal
            print sampHead[ "title" ] + " - " + strCut
            
            # Drawing normal plot (for isolation values)
            sampHead[ "hist" ] = getH1_Normalized(datadir + sampHead[ "filename" ], strTree, 
                sampHead[ "title" ], binCurr, plotvar, strCut)
        
        setMarkerStyle(sampHead[ "hist" ], sampHead[ "color" ], sampHead[ "shape" ])
        arrHist.append(sampHead)
        
        if nMin > sampHead[ "hist" ].GetMinimum(): 
            nMin = sampHead[ "hist" ].GetMinimum()
        
        if nMax < sampHead[ "hist" ].GetMaximum(): 
            nMax = sampHead[ "hist" ].GetMaximum()
    
    #Set init histo
    h_init = ROOT.TH1F("", "", binCurr[ 0 ], binCurr[ 1 ], binCurr[ 2 ])

    #Set axis
    x_name = "Muon " + dicPlotvar[ "xtitle" ]
    y_name = "Muon " + dicPlotvar[ "ytitle" ]
    
    #h_init.SetMinimum(0.0)
    #h_init.SetMaximum(nMax * 0.9)
    #h_init.SetMaximum(1.0)
    
    if "min" in dicPlotvar: h_init.SetMinimum(dicPlotvar[ "min" ])
    if "max" in dicPlotvar: h_init.SetMaximum(dicPlotvar[ "max" ])
    
    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.GetYaxis().SetTitleOffset(1)

    name = "%s_%s"%(plotvar, isotype)
    if "id" in dicPlotvar: name = name + "_%s"%dicPlotvar[ "id" ]

    #Set canvas
    canv = makeCanvas(name, False)
    setMargins(canv, False)
    h_init.Draw()
    #drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")
    drawSampleName(("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} and QCD events, "
        "p_{T} > %(pT)s GeV, |#eta| < %(Eta)s, %(ID)s Muon")%{"pT": strCutPT, "Eta": strCutEta, "ID": dicPlotvar[ "id" ]})
    
    if "ylog" in dicPlotvar and dicPlotvar[ "ylog" ]:
        ROOT.gPad.SetLogy()
    else: 
        ROOT.gPad.SetLogy(0)

    #Legend and drawing
    legTop = ROOT.TLegend(0.5,0.65,0.85,0.80)
    legBot = ROOT.TLegend(0.5,0.2,0.85,0.35)
    if "genMuon"  in plotvar: leg = legBot
    if "recoMuon" in plotvar: leg = legTop
    
    for sampHead in arrHist: 
        sampHead[ "hist" ].Draw("e1same")
        leg.AddEntry(sampHead[ "hist" ], sampHead[ "hist" ].GetTitle(), "p")

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
    canv.SaveAs(name+".png")
    
    print "[%s - %s] has been drawn"%(plotvar, "All" if "id" not in dicPlotvar else dicPlotvar[ "id" ])

