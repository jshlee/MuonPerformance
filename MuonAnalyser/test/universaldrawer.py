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


#datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
datadir = "/cms/scratch/quark2930/Work/muon_upgrade/samples/"

if len(sys.argv) < 2: 
    print "Usage : python universaldrawer.py (JSON file)"
    exit(1)

dicMainCmd = json.load(open(sys.argv[ 1 ]))

# Getting vital variables
try: 
    strCut = dicMainCmd[ "cut" ]
    binCurr = dicMainCmd[ "binning" ]
    
    dicCutConfig = dicMainCmd[ "cutconfig" ]
    
    strHistTitle = dicMainCmd[ "title" ]
    x_name = "Muon " + dicMainCmd[ "xtitle" ]
    y_name = "Muon " + dicMainCmd[ "ytitle" ]
    
    strFilename = dicMainCmd[ "filename" ]
    
    arrVars = dicMainCmd[ "vars" ]

except KeyError, strErr:
    print "Error: the JSON file does not contain a vital key: %s"%strErr

# Variables which can have a defalut value
strPlotvar = "" # It must be determined either in "general" or "vars"
nIsUsedCommonVar = 0

nIsLogY = 0

nIsUseMin = 0
nIsUseMax = 0
fMin =  1000000000
fMax = -1000000000

nIsEffRate = 0
strCutDen = ""

fLegLeft   = 0.50
fLegTop    = 0.65
fLegRight  = 0.85
fLegBottom = 0.80

# Now setup the configuration
if "plotvar" in dicMainCmd:
    strPlotvar = dicMainCmd[ "plotvar" ]
    nIsUsedCommonVar = 1

if "ylog" in dicMainCmd:
    nIsLogY = 1

if "min" in dicMainCmd: 
    fMin = dicMainCmd[ "min" ]
    nIsUseMin = 1

if "max" in dicMainCmd: 
    fMax = dicMainCmd[ "max" ]
    nIsUseMax = 1

if "effrate" in dicMainCmd: 
    strCutDen = dicMainCmd[ "effrate" ]
    nIsEffRate = 1

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

# Now all plots get being drawn
for varHead in arrVars: 
    if nIsUsedCommonVar == 0: 
        strPlotvar = varHead[ "plotvar" ]
    
    strTree = ""
    
    if "genMuon" in strPlotvar: 
        strTree = "MuonAnalyser/gen"
    if "recoMuon" in strPlotvar: 
        strTree = "MuonAnalyser/reco"
    
    if "cutconfig" in varHead: 
        for strKey in varHead[ "cutconfig" ].keys(): 
            dicCutConfig[ strKey ] = varHead[ "cutconfig" ][ strKey ]
    
    strCutExtra = ""
    if "cut" in varHead: strCutExtra = " && " + varHead[ "cut" ]
    
    if nIsEffRate != 0: 
        # Drawing efficiency / fake rate plot
        varHead[ "hist" ] = getEff(datadir + varHead[ "filename" ], strTree, 
            varHead[ "title" ], binCurr, strPlotvar, 
            ( strCut + strCutDen + strCutExtra )%dicCutConfig, # cut for denominator
            ( strCut +             strCutExtra )%dicCutConfig) # cut for nominator
    else:
        # Drawing normal plot (for isolation values)
        varHead[ "hist" ] = getH1_Normalized(datadir + varHead[ "filename" ], strTree, 
            varHead[ "title" ], binCurr, strPlotvar, 
            ( strCut + strCutExtra )%dicCutConfig)
    
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
drawSampleName(strHistTitle%dicCutConfig)

if nIsLogY != 0: ROOT.gPad.SetLogy()

#Legend and drawing
leg = ROOT.TLegend(fLegLeft, fLegTop, fLegRight, fLegBottom)

for varHead in arrVars: 
    strExtraOpt = ""
    if "extradrawopt" in varHead: strExtraOpt = varHead[ "extradrawopt" ]
    
    varHead[ "hist" ].Draw(strExtraOpt + "e1same")
    leg.AddEntry(varHead[ "hist" ], varHead[ "hist" ].GetTitle(), "p")

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
canv.SaveAs(strFilename%dicCutConfig)

print "%s has been drawn"%(strFilename%dicCutConfig)

