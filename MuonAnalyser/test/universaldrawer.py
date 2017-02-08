import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


hA = ""


def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(1.5)
    h.SetLineColor(color)
    h.SetLineWidth(2)


def getH1_Normalized(filename,treename,title,binning,plotvar,cut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,cut)
    
    normfactor = h1.GetEntries()
    #h1.Scale(1.0 / normfactor)
    
    h1.SetTitle(title)
    
    return copy.deepcopy(h1)


def getEff(filename,treename,title,binning,plotvar,dencut,numcut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,dencut)
    h2 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    hEff = ROOT.TEfficiency(h2,h1)
    hEff.SetTitle(title)
    hA = copy.deepcopy(h1)
    return copy.deepcopy(hEff)


def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.18, 0.8, samplename)


datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
#datadir = "TenMuExtendedE_"
#id = sys.argv[1]
#isotype = sys.argv[2]
isotype = sys.argv[1]
#isotype = "Trk"

listID = ["Tight", "Loose"]
listIso = ["Trk", "PF"]

"""
if id not in listID:
    print "Error : ID (first arg) should be Tight or Loose"
    exit(1)
"""

if isotype not in listIso:
    print "Error : Isolation condition (second arg) should be Trk or PF"
    exit(1)

arrPlotvar = [
    {"plotvar": "genMuon.Pt()",        "binning": [10,5,105]}, 
    {"plotvar": "abs(genMuon.Eta())",  "binning": [15,0,2.4]}, 
    #{"plotvar": "genMuon.Phi()",       "binning": [12,-3,3]}, 
    
    #{"plotvar": "recoMuon.Pt()",       "binning": [10,5,105]}, 
    #{"plotvar": "abs(recoMuon.Eta())", "binning": [8,0,2.4]}, 
    #{"plotvar": "recoMuon.Phi()",      "binning": [12,-3,3]}, 
    
    #{"plotvar": "recoMuon_TrkIsolation03", "binning": [20,0,1], "id": "Tight", "ylog": True, "type": "all"}, 
    #{"plotvar": "recoMuon_PFIsolation04",  "binning": [20,0,1], "id": "Tight", "ylog": True, "type": "all"}, 
    
    #{"plotvar": "recoMuon_TrkIsolation03", "binning": [20,0,1], "id": "Loose", "ylog": True, "type": "all"}, 
    #{"plotvar": "recoMuon_PFIsolation04",  "binning": [20,0,1], "id": "Loose", "ylog": True, "type": "all"}, 
]

arrSampleType = [
    {"title": "Phase II PU0",   "filename": "out_PU0.root",   "id": "Tight", "color": 4, "shape": 20}, # blue,  circle
    {"title": "Phase II PU200", "filename": "out_PU200.root", "id": "Tight", "color": 2, "shape": 21}, # red,   square
    {"title": "Phase II QCD",   "filename": "out_QCD.root",   "id": "Tight", "color": 1, "shape": 34}, # black, cross
    
    {"title": "Phase II PU0",   "filename": "out_PU0.root",   "id": "Loose", "color": 4, "shape": 20}, # blue,  circle
    {"title": "Phase II PU200", "filename": "out_PU200.root", "id": "Loose", "color": 2, "shape": 21}, # red,   square
    {"title": "Phase II QCD",   "filename": "out_QCD.root",   "id": "Loose", "color": 1, "shape": 34}, # black, cross
];

listIDCfg = {
    "Tight": {"isoCut": {"Trk": "TrkIsolation03 < 0.05", "PF": "PFIsolation04 < 0.15"}}, 
    "Loose": {"isoCut": {"Trk": "TrkIsolation03 < 0.10", "PF": "PFIsolation04 < 0.25"}}, 
}

strCutGenNor = "genMuon.Pt() > 5"
strCutGenDen = strCutGenNor
strCutGenNum = strCutGenDen + " && genMuon_is%(ID)s && genMuon_%(iso)s"

strCutRecNor = "recoMuon.Pt() > 5"
strCutRecDen = strCutRecNor + " && recoMuon_is%(ID)s && recoMuon_%(iso)s"
strCutRecNum = strCutRecDen + " && !recoMuon_signal"

strCutGenIso = ""
strCutRecIso = ""


for i, dicPlotvar in enumerate(arrPlotvar):
    print "######################"
    print "plotvar : " + dicPlotvar[ "plotvar" ]
    print "######################"
    
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
    
    """
    if "Isolation" not in plotvar: 
        strCutDen = strCutDen + ( " && " + strCutIso if "genMuon" in plotvat else "" )
        strCutNum = strCutNum + " && " + strCutIso
    """
    
    arrHist = []
    
    #Get histos
    for sampHead in arrSampleType: 
        if "type" not in dicPlotvar or dicPlotvar[ "type" ] != "all": 
            if "gen" in plotvar and "QCD" in sampHead[ "title" ]: 
                continue
            if "reco" in plotvar and "PU" in sampHead[ "title" ]: 
                continue

        strIDCurr = sampHead[ "id" ]
        
        # This is the isolation cut condition
        strCutIso = listIDCfg[ strIDCurr ][ "isoCut" ][ isotype ]
        
        if "Isolation" not in plotvar: 
            # Drawing efficiency / fake rate plot
            sampHead[ "hist" ] = getH1_Normalized(datadir + sampHead[ "filename" ], strTree, 
                sampHead[ "title" ] + " - " + strIDCurr, binCurr, plotvar, 
                strCutDen%{"ID": strIDCurr, "iso": strCutIso}) # cut for normal
            """
            sampHead[ "hist" ] = getEff(datadir + sampHead[ "filename" ], strTree, 
                sampHead[ "title" ] + " - " + strIDCurr, binCurr, plotvar, 
                strCutDen%{"ID": strIDCurr, "iso": strCutIso}, # cut for denominator
                strCutNum%{"ID": strIDCurr, "iso": strCutIso}) # cut for nominator
            """
            
            print "  " + sampHead[ "title" ] + " - " + strIDCurr
            print "    (" + strCutDen%{"ID": strIDCurr, "iso": strCutIso} + ", " + strCutNum%{"ID": strIDCurr, "iso": strCutIso} + ")"
        else:
            if strIDCurr != dicPlotvar[ "id" ]:
                continue
            
            # Drawing normal plot (for isolation values)
            sampHead[ "hist" ] = getH1_Normalized(datadir + sampHead[ "filename" ], strTree, 
                sampHead[ "title" ], binCurr, plotvar, 
                strCutNor%{"ID": dicPlotvar[ "id" ], "iso": strCutIso}) # cut for normal
            
            print "  " + sampHead[ "title" ] + " - " + strIDCurr
            print "    (" + strCutNor%{"ID": dicPlotvar[ "id" ], "iso": strCutIso} + ")"
        
        setMarkerStyle(sampHead[ "hist" ], sampHead[ "color" ], sampHead[ "shape" ])
        arrHist.append(sampHead)
    
    #Set init histo
    h_init = ROOT.TH1F("", "", binCurr[ 0 ], binCurr[ 1 ], binCurr[ 2 ])

    #Set axis
    x_name = "Muon "
    if "Pt"  in plotvar: x_name = x_name+"p_{T}"
    if "Eta" in plotvar: x_name = x_name+"|#eta|"
    if "Phi" in plotvar: x_name = x_name+"#phi"

    y_name = "Muon "
    
    if "genMuon" in plotvar:
        #h_init.SetMaximum(1.1)
        #h_init.SetMinimum(0.6)
        y_name = y_name+"Efficiency"
    
    if "recoMuon" in plotvar:
        #h_init.SetMaximum(max(h.GetMaximum() for h in hlist)*2.5)
        h_init.GetYaxis().SetLabelSize(0.035)
        h_init.GetYaxis().SetTitleOffset(1.2)
        y_name = y_name+"Fake Rate"
    
    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.GetYaxis().SetTitleOffset(1)

    name = "%s_%s"%(plotvar, isotype)

    #Set canvas
    canv = makeCanvas(name, False)
    setMargins(canv, False)
    h_init.Draw()
    drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")
    
    if "ylog" in dicPlotvar and dicPlotvar[ "ylog" ]:
        ROOT.gPad.SetLogy()
    else: 
        ROOT.gPad.SetLogy(0)

    #Legend and drawing
    legTop = ROOT.TLegend(0.6,0.7,0.85,0.85)
    legBot = ROOT.TLegend(0.6,0.2,0.85,0.35)
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

