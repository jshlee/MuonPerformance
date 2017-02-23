import ROOT, copy, os, sys, json
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


def setMarkerStyle(h,color,style,size):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(size)
    h.SetLineColor(color)
    h.SetLineWidth(2)


def getH1_Normalized(filename,treename,title,binning,plotvar,cut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,cut)
    
    # Normalizing
    normfactor = h1.GetEntries()
    if normfactor > 0.0: h1.Scale(1.0 / normfactor)
    
    # Showing overflow
    nValLastBin  = h1.GetBinContent(binning[ 0 ])
    nValOverflow = h1.GetBinContent(binning[ 0 ] + 1)
    h1.SetBinContent(binning[ 0 ], nValLastBin + nValOverflow)
    h1.SetBinContent(binning[ 0 ] + 1, 0)
    
    # Setting title
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

if len(sys.argv) < 2: 
    print "Usage : python universaldrawer.py (JSON file)"
    exit(1)

dicMainCmd = json.load(open(sys.argv[ 1 ]))

arrPlotvar =    dicMainCmd[ "plotvars" ]
arrSampleType = dicMainCmd[ "samples" ]

listIDCfg = {
    "Tight": {"isoCut": {"Trk": "TrkIsolation03 < 0.05", "PF": "PFIsolation04 < 0.15"}}, 
    "Loose": {"isoCut": {"Trk": "TrkIsolation03 < 0.10", "PF": "PFIsolation04 < 0.25"}}, 
}

strCutPT  = "15"
strCutEta = "2"

strCutGenNor = "genMuon.Pt() > %(pT)s && abs(genMuon.Eta()) < %(Eta)s"%{"pT":strCutPT, "Eta": strCutEta}
strCutGenDen = strCutGenNor
#strCutGenNum = strCutGenDen + " && genMuon_is%(ID)s && genMuon_%(iso)s" # ID is held
strCutGenNum = strCutGenDen + " && genMuon_isMuon && genMuon_%(iso)s"

#strCutRecNor = "recoMuon.Pt() > 5 && recoMuon_is%(ID)s"
#strCutRecNor = "recoMuon.Pt() > 5"
#strCutRecNor = "recoMuon.Pt() > 5 && abs(recoMuon.Eta()) < 2.4 && recoMuon_TrkIsolation03 <= 0.00001"
strCutRecNor = "recoMuon.Pt() > %(pT)s && abs(recoMuon.Eta()) < %(Eta)s"%{"pT":strCutPT, "Eta": strCutEta}
#strCutRecNor = "recoMuon.Pt() > 20 && abs(recoMuon.Eta()) < 2"
strCutRecDen = strCutRecNor + " && recoMuon_isMuon && recoMuon_%(iso)s"
strCutRecNum = strCutRecDen + " && !recoMuon_signal"

strCutGenIso = ""
strCutRecIso = ""

dicCutConfig = {
    "pT": strCutPT, "Eta": strCutEta
}

if "cutconfig" in dicMainCmd: dicCutConfig = dicMainCmd[ "cutconfig" ]


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
    
    if "id" in dicMainCmd: dicCutConfig[ "ID" ] = dicMainCmd[ "id" ]
    
    #Get histos
    for sampHead in arrSampleType: 
        if "type" not in dicPlotvar or dicPlotvar[ "type" ] != "all": 
            if "gen" in plotvar and "QCD" in sampHead[ "title" ]: 
                continue
            if "reco" in plotvar and "QCD" not in sampHead[ "title" ]: 
                continue

        strIDCurr = sampHead[ "id" ]
        dicCutConfig[ "ID" ] = strIDCurr
        
        # This is the isolation cut condition
        isotype = dicPlotvar[ "isotype" ]
        dicMainCmd[ "iso" ] = listIDCfg[ strIDCurr ][ "isoCut" ][ isotype ]
        
        strCutExtra = ""
        if "extracut" in sampHead: strCutExtra = " && " + sampHead[ "extracut" ]
        
        #if "Iso" not in plotvar: 
        if "effrate" in dicPlotvar and dicPlotvar[ "effrate" ]: 
            # Drawing efficiency / fake rate plot
            sampHead[ "hist" ] = getEff(datadir + sampHead[ "filename" ], strTree, 
                sampHead[ "title" ] + " - " + strIDCurr, binCurr, plotvar, 
                strCutDen%dicCutConfig + strCutExtra, # cut for denominator
                strCutNum%dicCutConfig + strCutExtra) # cut for nominator
        else:
            if strIDCurr != dicPlotvar[ "id" ]:
                continue
            
            strCut = strCutNor%dicCutConfig + strCutExtra # cut for normal
            print sampHead[ "title" ] + " - " + strCut
            
            # Drawing normal plot (for isolation values)
            sampHead[ "hist" ] = getH1_Normalized(datadir + sampHead[ "filename" ], strTree, 
                sampHead[ "title" ], binCurr, plotvar, strCut)
        
        fSizeDot = 1.5
        if "size" in sampHead: fSizeDot = sampHead[ "size" ]
        
        setMarkerStyle(sampHead[ "hist" ], sampHead[ "color" ], sampHead[ "shape" ], fSizeDot)
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
    
    if "ylog" in dicPlotvar and dicPlotvar[ "ylog" ] != 0:
        ROOT.gPad.SetLogy()
    else: 
        ROOT.gPad.SetLogy(0)

    #Legend and drawing
    legTop = ROOT.TLegend(0.5,0.65,0.85,0.80)
    legBot = ROOT.TLegend(0.5,0.2,0.85,0.35)
    if "genMuon"  in plotvar: leg = legBot
    if "recoMuon" in plotvar: leg = legTop
    
    if "legend" in dicMainCmd: 
        fLegLeft   = 0.5
        fLegTop    = 0.65
        fLegRight  = 0.80
        fLegBottom = 0.5
        
        if "left"   in dicMainCmd[ "legend" ]: 
            fLegLeft   = dicMainCmd[ "legend" ][ "left" ]
        
        if "top"    in dicMainCmd[ "legend" ]: 
            fLegTop    = dicMainCmd[ "legend" ][ "top" ]
        
        if "right"  in dicMainCmd[ "legend" ]: 
            fLegRight  = dicMainCmd[ "legend" ][ "right" ]
        
        if "bottom" in dicMainCmd[ "legend" ]: 
            fLegBottom = dicMainCmd[ "legend" ][ "bottom" ]
        
        leg = ROOT.TLegend(fLegLeft, fLegTop, fLegRight, fLegBottom)
    
    for sampHead in arrHist: 
        strExtraOpt = ""
        if "extradrawopt" in sampHead: strExtraOpt = sampHead[ "extradrawopt" ]
        
        sampHead[ "hist" ].Draw(strExtraOpt + "e1same")
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
    if "filename" in dicPlotvar: 
        canv.SaveAs(dicPlotvar[ "filename" ])
    else: 
        canv.SaveAs(name+".png")
    
    print "[%s - %s] has been drawn"%(plotvar, "All" if "id" not in dicPlotvar else dicPlotvar[ "id" ])

