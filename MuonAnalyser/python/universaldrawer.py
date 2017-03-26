########################################################################################
## ## universaldrawer.py
##
## usage : python universaldrawer.py JSON_File
##         or import this into a python program and use drawPlotFromDict().
##         e.g., drawPlotFromDict(json.load(open([JSON file])))
##
## The following is the form of JSON file.
## You can see an example in 
##   test/jsonconf/jsonconf_Isovalue_PFIso04_Tight.json or test/drawtest.py.
## Perhaps, it is more helpful to read only this examples and 
##   skip decryption of the following descriptions.
## 
## Every item is marked by (L) or (D) or (V) or (C).
## If makred by (L), it should be a list.
## If makred by (D), it should be a dictinoary.
## If makred by (V), it should be a single value, such as a number or a string.
## If makred by (C), it should be a string which can contain format specifier.
##   The format specifier will be filled with values in cutconfig.
##   In Python, format specifiers can have a 'key', so the following is possible.
## 
##   "recoMuon.Pt() > %(pT)s && abs(recoMuon.Eta()) < %(Eta)s"%{"pT": "15", "Eta": "2.4"}
##     => "recoMuon.Pt() > 15 && abs(recoMuon.Eta()) < 2.4"
## 
##   In the value with this mark you can use this specifiers, 
##   where the followings are some possible keys.
## 
##   "pT"  : for pT cut.
##   "Eta" : for eta cut.
##   "ID"  : for ID cut. It MUST be "Tight" or "Loose".
## 
##   Of course, you can set up this setups as your taste; you can add any other keys
##   and even type of a specifier.
## 
## 
## (Required)
## 
## "tree" (V) : the tree which this program uses
## 
## "cutconfig" (D) : cut configuration values, such as pT and eta for cut are in here.
##   This values (and some in "vars") will be applied to values marked by (C), as mentioned.
## 
## "cut" (C) : Cut condition.
## 
## "binning" (L) : Binning of histogram.
##   For more information, see ROOT:TH1 or and so on.
## 
## "title" (C) : This string will be printed at the histogram and 
##   describes what this plot is.
## 
## "xtitle" (V) : title of x-axis
## 
## "ytitle" (V) : title of y-axis
## 
## "filename" (C) : the name of output image (or root) file
## 
## "vars" (L) : you can give the data in this list.
##   For more information, see below.
## 
## 
## (Optional)
## 
## "plotvar" (V) : What value do you want to draw?
##   This value can be set in "vars" if you want to draw several values.
## 
## "plottype" (V) : Type of the plot
##   This program has 3 types: 
##     "normal" (default; if "plottype" is not given, it will be applied)
##     "effrate" : See "effrate" key (if it is not given but "effrate" key is there, it turns on)
##     "bkgrate" : See "bkgrate" key (if it is not given but "bkgrate" key is there, it turns on)
##     "accumulated" : drawing accumulated plot
## 
## "normoff" (V) : If this key is there, plots are normalized. Note that it is turned on default!
##   (not for "effrate")
## 
## "ylog" (V) : Using log scale on y-axis
## 
## "min" (V) : you can set the minimum of the plot.
## 
## "max" (V) : you can set the maximum of the plot.
## 
## "effrate" (C) : This variable contains a cut condition of nominator 
##   for efficiency (or fake) rate plot.
## 
## "bkgrate" (C) : This variable contains a cut condition of nominator 
##   for background rate plot, where the denominator is the total number of events
## 
## "bkgdenominator" (D) : It contains the information for how to set denominator; 
##   it contains the name of histogram and the range
## 
## "legend" (D) : you can customize the position and size of the legend.
##   It should be given as a dictionary with keys "left", "top", "right", "bottom".
##   If one of these keys are not given, this program will use the default setting.
##   For more information, see ROOT::TLegend.
## 
## 
## The followings are for description of "vars".
##   Each in the list of "vars" is a dictionary containing the following keys and values.
## 
## (Required)
## 
## "filename" (V) : the root file containng samples.
## 
## "title" (V) : you can see this value in the legend.
## 
## "color" (V) : the color of point in the histogram.
## 
## "shape" (V) : the shape of point in the histogram.
## 
## 
## (Optional)
## 
## "plotvar" (V) : What value do you want to draw?
##   It is optional, but if there is no central "plotvar", it becomes required.
## 
## "cutconfig" (C) : Additional cut configuration values for this samples.
## 
## "cut" (C) : Additional cut condition for this samples.
## 
## "effrate" (C) : This variable contains a cut condition of nominator 
##   for efficiency (or fake or background) rate plot.
## 
## "extradrawopt" (V) : Extra option for drawing histogram.
##   For more information, see ROOT:TH1.
## 
## "size" (V) : the size of point in the histogram.
## 
########################################################################################


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


def getH1_Normalized(filename,treename,title,binning,plotvar,cut,normaloff=False):
    h1 = makeTH1(filename,treename,title,binning,plotvar,cut)
    
    # Normalizing
    if not normaloff: 
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


def getH1_Accumulated(filename,treename,title,binning,plotvar,cut,normaloff=False):
    h1 = makeTH1(filename,treename,title,binning,plotvar,cut)
    
    # Normalizing
    if not normaloff: 
      normfactor = h1.GetEntries()
      if normfactor > 0.0: h1.Scale(1.0 / normfactor)
    
    # Showing overflow
    nValLastBin  = h1.GetBinContent(binning[ 0 ])
    nValOverflow = h1.GetBinContent(binning[ 0 ] + 1)
    h1.SetBinContent(binning[ 0 ], nValLastBin + nValOverflow)
    h1.SetBinContent(binning[ 0 ] + 1, 0)
    
    # Accumulating
    for i in range(binning[ 0 ]): 
        fBinPrev = h1.GetBinContent(i + 1)
        fBinCurr = h1.GetBinContent(i + 2)
        
        h1.SetBinContent(i + 2, fBinPrev + fBinCurr)
    
    # Setting title
    h1.SetTitle(title)
    
    return copy.deepcopy(h1)


def getEff(filename,treename,title,binning,plotvar,dencut,numcut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,dencut)
    h2 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    hEff = ROOT.TEfficiency(h2,h1)
    hEff.SetTitle(title)
    return copy.deepcopy(hEff)


def getBkg(filename,treename,title,binning,plotvar,numcut,denconf):
    # Scaling by # of events
    f = ROOT.TFile(filename)
    hNorm = f.Get("MuonAnalyser/" + denconf[ "name" ].encode("ascii", "ignore"))
    
    nBinMin = hNorm.GetXaxis().FindBin(denconf[ "min" ])
    nBinMax = hNorm.GetXaxis().FindBin(denconf[ "max" ])
    normfactor = hNorm.Integral(nBinMin, nBinMax)
    
    h1 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    if normfactor > 0.0: h1.Scale(1.0 / normfactor)
    
    # Showing overflow
    nValLastBin  = h1.GetBinContent(binning[ 0 ])
    nValOverflow = h1.GetBinContent(binning[ 0 ] + 1)
    h1.SetBinContent(binning[ 0 ], nValLastBin + nValOverflow)
    h1.SetBinContent(binning[ 0 ] + 1, 0)
    
    # Setting title
    h1.SetTitle(title)
    
    return copy.deepcopy(h1)


def drawSampleName(samplename, fX, fY):
    #fX = 0.18
    #fY = 0.85
    fSizeTex = 0.038
    
    tex2 = ROOT.TLatex()
    
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(fSizeTex)
    
    for i, strLine in enumerate(samplename.split("\n")): 
        tex2.DrawLatex(fX, fY - i * 1.1 * fSizeTex, strLine)


def drawPlotFromDict(dicMainCmd): 
    #datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
    #datadir = "/cms/scratch/quark2930/Work/muon_upgrade/samples/"
    datadir = ""

    #if len(sys.argv) < 2: 
    #    print "Usage : python universaldrawer.py (JSON file)"
    #    exit(1)

    #dicMainCmd = json.load(open(sys.argv[ 1 ]))

    # Getting required variables
    try: 
        strTreePre = dicMainCmd[ "tree" ]
        strCut = dicMainCmd[ "cut" ]
        binCurr = dicMainCmd[ "binning" ]
        
        dicCutConfig = dicMainCmd[ "cutconfig" ]
        
        strHistTitle = dicMainCmd[ "title" ]
        x_name = "Muon " + dicMainCmd[ "xtitle" ]
        y_name = "Muon " + dicMainCmd[ "ytitle" ]
        
        strFilename = dicMainCmd[ "filename" ]
        
        arrVars = dicMainCmd[ "vars" ]

    except KeyError, strErr:
        print "Error: the JSON file does not contain a required key: %s"%strErr

    # Variables which can have a defalut value
    strPlotvar = "" # It must be determined either in "general" or "vars"
    strPlotType = "normal"
    nIsUsedCommonVar = 0

    nIsLogY = 0

    nIsUseMin = 0
    nIsUseMax = 0
    fMin =  1000000000
    fMax = -1000000000

    strCutNum = ""

    fLegLeft   = 0.50
    fLegTop    = 0.65
    fLegRight  = 0.85
    fLegBottom = 0.80
    
    fLegTextSize = 0.04
    
    fTitleX = 0.18
    fTitleY = 0.85
    
    bNormalOff = False
    
    dicBkgDen = {"name": "nevents", "min": 0, "max": 1}

    # Now setup the configuration
    if "plotvar" in dicMainCmd:
        strPlotvar = dicMainCmd[ "plotvar" ]
        nIsUsedCommonVar = 1
    
    if "plottype" in dicMainCmd: 
        listPossible = [
            "normal", "effrate", "bkgrate", "accumulated"
        ]
        
        if dicMainCmd[ "plottype" ] not in listPossible: 
            print "Error : \"%s\" is not a possible type."%dicMainCmd[ "plottype" ]
            sys.exit(1)
        
        strPlotType = dicMainCmd[ "plottype" ]

    if "ylog" in dicMainCmd:
        nIsLogY = 1

    if "min" in dicMainCmd: 
        fMin = dicMainCmd[ "min" ]
        nIsUseMin = 1

    if "max" in dicMainCmd: 
        fMax = dicMainCmd[ "max" ]
        nIsUseMax = 1

    if "effrate" in dicMainCmd: 
        strCutNum = " && " + dicMainCmd[ "effrate" ]
        strPlotType = "effrate"

    if "bkgrate" in dicMainCmd: 
        strCutNum = " && " + dicMainCmd[ "bkgrate" ]
        strPlotType = "bkgrate"
    
    if "bkgdenominator" in dicMainCmd: 
        dicBkgDen = dicMainCmd[ "bkgdenominator" ]
        strPlotType = "bkgrate"
    
    if "titlepos" in dicMainCmd: 
        fTitleX = dicMainCmd[ "titlepos" ][ 0 ]
        fTitleY = dicMainCmd[ "titlepos" ][ 1 ]

    if "legend" in dicMainCmd: 
        try: 
            fLegLeft   = dicMainCmd[ "legend" ][ "left" ]
            fLegTop    = dicMainCmd[ "legend" ][ "top" ]
            fLegRight  = dicMainCmd[ "legend" ][ "right" ]
            fLegBottom = dicMainCmd[ "legend" ][ "bottom" ]
        except KeyError, strErr:
            print "Error: Wrong legend configuration : %s is missing"%strErr
            
            fLegLeft   = 0.50
            fLegTop    = 0.65
            fLegRight  = 0.85
            fLegBottom = 0.80
        
        if "textsize" in dicMainCmd[ "legend" ]: 
            fLegTextSize = dicMainCmd[ "legend" ][ "textsize" ]
    
    if "normoff" in dicMainCmd: 
        bNormalOff = True

    # Now all plots get being drawn
    for varHead in arrVars: 
        if nIsUsedCommonVar == 0: 
            strPlotvar = varHead[ "plotvar" ]
        
        #strTree = ""
        strCutNumVar = strCutNum
        
        #if "genMuon" in strPlotvar or "genMuon" in strCut: 
        #    strTree = "MuonAnalyser/gen"
        #if "recoMuon" in strPlotvar or "recoMuon" in strCut: 
        #    strTree = "MuonAnalyser/reco"
        if strTreePre == "gen": # because strTreePre may be in unicode
            strTree = "MuonAnalyser/gen"
        if strTreePre == "reco": 
            strTree = "MuonAnalyser/reco"
        
        if "cutconfig" in varHead: 
            for strKey in varHead[ "cutconfig" ].keys(): 
                dicCutConfig[ strKey ] = varHead[ "cutconfig" ][ strKey ]
        
        strCutExtra = ""
        if "cut" in varHead: strCutExtra = " && " + varHead[ "cut" ]
        if "effrate" in varHead: strCutNumVar = strCutNumVar + " && " + varHead[ "effrate" ]
        if "bkgrate" in varHead: strCutNumVar = strCutNumVar + " && " + varHead[ "bkgrate" ]
        
        if strPlotType == "normal": 
            print "Cut : %s"%(( strCut + strCutExtra )%dicCutConfig)
            # Drawing normal plot (for isolation values)
            varHead[ "hist" ] = getH1_Normalized(datadir + varHead[ "filename" ], strTree, 
                varHead[ "title" ], binCurr, strPlotvar, 
                ( strCut + strCutExtra )%dicCutConfig, bNormalOff)
        elif strPlotType == "effrate": 
            print "Num : %s\nDen : %s"%(( strCut + strCutNumVar + strCutExtra )%dicCutConfig, ( strCut + strCutExtra )%dicCutConfig)
            # Drawing efficiency / fake rate plot
            varHead[ "hist" ] = getEff(datadir + varHead[ "filename" ], strTree, 
                varHead[ "title" ], binCurr, strPlotvar, 
                ( strCut +                strCutExtra )%dicCutConfig, # cut for denominator
                ( strCut + strCutNumVar + strCutExtra )%dicCutConfig) # cut for nominator
        elif strPlotType == "bkgrate": 
            print "Num : %s"%(( strCut + strCutNumVar + strCutExtra )%dicCutConfig)
            # Drawing normal plot (for isolation values)
            varHead[ "hist" ] = getBkg(datadir + varHead[ "filename" ], strTree, 
                varHead[ "title" ], binCurr, strPlotvar, 
                ( strCut + strCutNumVar + strCutExtra )%dicCutConfig, # cut for nominator
                dicBkgDen)
        elif strPlotType == "accumulated": 
            print "Cut : %s"%(( strCut + strCutExtra )%dicCutConfig)
            # Drawing normal plot (for isolation values)
            varHead[ "hist" ] = getH1_Accumulated(datadir + varHead[ "filename" ], strTree, 
                varHead[ "title" ], binCurr, strPlotvar, 
                ( strCut + strCutExtra )%dicCutConfig, bNormalOff)
        
        fSizeDot = 1.0
        if "size" in varHead: fSizeDot = varHead[ "size" ]
        
        setMarkerStyle(varHead[ "hist" ], varHead[ "color" ], varHead[ "shape" ], fSizeDot)
        
        if nIsUseMin == 0 and fMin > varHead[ "hist" ].GetMinimum(): 
            fMin = varHead[ "hist" ].GetMinimum()
        
        if nIsUseMax == 0 and fMax < varHead[ "hist" ].GetMaximum(): 
            fMax = varHead[ "hist" ].GetMaximum()
        
    # The remainings are for drawing the total plot
    h_init = ROOT.TH1F("", "", binCurr[ 0 ], binCurr[ 1 ], binCurr[ 2 ])

    if nIsUseMax == 0: 
        fMax = fMax * 1.15

    h_init.SetMinimum(fMin)
    h_init.SetMaximum(fMax)

    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.GetYaxis().SetTitleOffset(1)

    #Set canvas
    canv = makeCanvas("canvMain", False)
    setMargins(canv, False)
    h_init.Draw()
    
    strTmpID = dicCutConfig[ "ID" ]
    
    if "Tight" in dicCutConfig[ "ID" ]: 
        dicCutConfig[ "ID" ] = "Tight"
    
    if "Loose" in dicCutConfig[ "ID" ]: 
        dicCutConfig[ "ID" ] = "Loose"
    
    drawSampleName(strHistTitle%dicCutConfig, fTitleX, fTitleY)
    
    dicCutConfig[ "ID" ] = strTmpID

    if nIsLogY != 0: ROOT.gPad.SetLogy()

    #Legend and drawing
    leg = ROOT.TLegend(fLegLeft, fLegTop, fLegRight, fLegBottom)

    for varHead in arrVars: 
        strExtraOpt = ""
        if "extradrawopt" in varHead: strExtraOpt = varHead[ "extradrawopt" ]
        
        varHead[ "hist" ].Draw(strExtraOpt + "same")
        leg.AddEntry(varHead[ "hist" ], varHead[ "hist" ].GetTitle(), "pl")

    leg.SetTextFont(61)
    leg.SetTextSize(fLegTextSize)
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
    canv.SaveAs(strFilename%dicCutConfig)

    print "%s has been drawn"%(strFilename%dicCutConfig)


if __name__ == "__main__": 
    if len(sys.argv) < 2: 
        print "Usage : python universaldrawer.py JSON_file"
        print "        or import this into a python program and use drawPlotFromDict()."
        exit(1)

    drawPlotFromDict(json.load(open(sys.argv[ 1 ])))


