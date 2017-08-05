import ROOT, copy, array, json, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi_forTDR as CMS_lumi
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


def getRefBkgVal(fileBkg,treename,cutBkg,cutIso):
  binning = [3000,0,30]
  hBkgDen = makeTH1(fileBkg,treename,"RefBkgValDen",binning,"muon_PFIsolation04",cutBkg)
  hBkgNum = makeTH1(fileBkg,treename,"RefBkgValNum",binning,"muon_PFIsolation04",cutBkg+" && "+cutIso)
  
  fMBkgDen = hBkgDen.Integral(0, binning[ 0 ] + 1)
  fMBkgNum = hBkgNum.Integral(0, binning[ 0 ] + 1)
  
  return 1.0 - fMBkgNum / fMBkgDen


def getROC(fileSig,fileBkg,treename,title,binning,plotvar,cutSig,cutBkg,dicRef):
  hSig = makeTH1(fileSig,treename,title,binning,plotvar,cutSig)
  hBkg = makeTH1(fileBkg,treename,title,binning,plotvar,cutBkg)
  
  arrSigEff = []
  arrBkgVal = []
  
  fMSig = hSig.Integral(0, binning[ 0 ] + 1)
  fMBkg = hBkg.Integral(0, binning[ 0 ] + 1)
  
  fSigEffPrev = 1.0
  fBkgValPrev = 1.0
  
  nPerf = 0.0
  fSigEffChosen = 0.0
  fBkgValChosen = 0.0
  
  fBkgMax = 0.0
  
  for i in range(binning[ 0 ] + 2): 
    # DO NOT confuse this : 0th bin has all underflows, while (binning[0] + 1)-th bin has all overflows
    # Well, we have no underflows, but... for habit for safety.
    nX = i + 0
    
    fSigIntegral = hSig.Integral(0, nX)
    fBkgIntegral = hBkg.Integral(0, nX)
    
    fSigEff = fSigIntegral / fMSig
    fBkgVal = 1.0 - fBkgIntegral / fMBkg
    
    if "bkgcustom" in dicRef:
      dicBkgCfg = dicRef[ "bkgcustom" ]
      
      if dicBkgCfg[ "name" ] == "bkgyield":
        fBkgVal = fBkgIntegral / dicBkgCfg[ "numEvents" ]
      elif dicBkgCfg[ "name" ] == "bkgeff":
        fBkgVal = fBkgIntegral / fMBkg
    
    if "sigeff" in dicRef: 
      if fSigEffPrev <= dicRef[ "sigeff" ] and dicRef[ "sigeff" ] <= fSigEff: 
        nPerf = i
        fSigEffChosen = fSigEff
        fBkgValChosen = fBkgVal
    elif "bkgrej" in dicRef: 
      if fBkgVal <= dicRef[ "bkgrej" ] and dicRef[ "bkgrej" ] <= fBkgValPrev: 
        nPerf = i
        fSigEffChosen = fSigEff
        fBkgValChosen = fBkgVal
    
    arrSigEff.append(fSigEff)
    arrBkgVal.append(fBkgVal)
    
    fSigEffPrev = fSigEff
    fBkgValPrev = fBkgVal
    
    if fBkgMax < fBkgVal: 
      fBkgMax = fBkgVal
  
  arrX = array.array("d", arrSigEff)
  arrY = array.array("d", arrBkgVal)
  
  graphROC = ROOT.TGraph(binning[ 0 ] + 2, arrX, arrY)
  graphROC.SetTitle(title)
  
  return {"graph": copy.deepcopy(graphROC), 
    "performance": binning[1] + nPerf * ( binning[2] - binning[1] ) / binning[0], 
    "sigeff": fSigEffChosen, "bkgrej": fBkgValChosen, "bkgmax": fBkgMax}


def drawSampleName(samplename, fX, fY, fSizeTex):
  tex2 = ROOT.TLatex()
  
  tex2.SetNDC()
  tex2.SetTextFont(62)
  tex2.SetTextSize(fSizeTex)
  
  for i, strLine in enumerate(samplename.split("\n")): 
    tex2.DrawLatex(fX, fY - i * 1.2 * fSizeTex, strLine)


def myGetVarID(varID):
  if type(varID) is dict:
    return varID[ "var" ]
  else:
    return varID


def myGetLabelID(varID):
  if type(varID) is dict:
    return varID[ "label" ]
  else:
    return varID


def myGetDicForFormatVar(dicCutConfig):
  dicRes = copy.deepcopy(dicCutConfig)
  if "ID" in dicRes: dicRes[ "ID" ] = myGetVarID(dicRes[ "ID" ])
  return dicRes


def myGetDicForFormatLabel(dicCutConfig):
  dicRes = copy.deepcopy(dicCutConfig)
  if "ID" in dicRes: dicRes[ "ID" ] = myGetLabelID(dicRes[ "ID" ])
  return dicRes


def drawPlotFromDict(dicMainCmd): 
  ## Reading input dictionary
  ## (Required)
  try: 
    binMain = dicMainCmd[ "binning" ]
    arrPlotvar = dicMainCmd[ "vars" ]
    
    strCut = dicMainCmd[ "cut" ]
    dicCutConfig = dicMainCmd[ "cutconfig" ]
    
    strTitleMain = dicMainCmd[ "title" ]
    
    SampleSig = dicMainCmd[ "inputfilename_sig" ]
    SampleBkg = dicMainCmd[ "inputfilename_bkg" ]
    
    strPrintoutFront = dicMainCmd[ "printout_front" ]
    strFilenameOutput = dicMainCmd[ "filename" ]
  except KeyError, strErr:
    print "Error: the JSON file does not contain a required key: %s" % strErr
    exit(1)
  
  ## Reading input dictionary
  ## (Optional)
  strMainTree = "MuonAnalyser"
  
  fSigEff = 0.0 if "sigeff" not in dicMainCmd else dicMainCmd[ "sigeff" ]
  
  dicLegPos = {
    "left":   0.16, 
    "top":    0.40, 
    "right":  0.45, 
    "bottom": 0.20, 
    
    "fonttype": 42, 
    "fontsize": 0.04, 
  }
  
  fTitleX = 0.14
  fTitleY = 0.85
  fTitleSize = 0.03
  
  fXMin = 0.0
  fXMax = 0.0
  
  fYMin = 0.0
  fYMax = 1.05
  
  strLabelYAxis = "Background rejection"
  
  strNameTagCut    = "extracut_"
  #strNameTagVarCut = "extravarcut_"
  
  strExtraText = "Working Progress"
  
  for strItemForCut in dicMainCmd.keys(): 
    if strNameTagCut in strItemForCut: 
      dicCutConfig[ strItemForCut.replace(strNameTagCut, "") ] = dicMainCmd[ strItemForCut ]
    #elif strNameTagVarCut in strItemForCut: 
    #  dicCutExtraVarConfig[ strItemForCut.replace(strNameTagVarCut, "") ] = dicMainCmd[ strItemForCut ]
  
  if "maintree" in dicMainCmd:
    strMainTree = dicMainCmd[ "maintree" ].encode("ascii", "ignore")
  
  if "sigeff" not in dicMainCmd: 
    if "ID" in dicCutConfig: 
      if "Tight" in myGetVarID(dicCutConfig[ "ID" ]): 
        fSigEff = 0.95
      elif "Loose" in myGetVarID(dicCutConfig[ "ID" ]): 
        fSigEff = 0.98
  
  if "legend" in dicMainCmd: 
    for strKey in dicMainCmd[ "legend" ]: 
      dicLegPos[ strKey ] = dicMainCmd[ "legend" ][ strKey ]
  
  if "titlepos" in dicMainCmd: 
    if "x" in dicMainCmd[ "titlepos" ]: fTitleX = dicMainCmd[ "titlepos" ][ "x" ]
    if "y" in dicMainCmd[ "titlepos" ]: fTitleY = dicMainCmd[ "titlepos" ][ "y" ]
    if "size" in dicMainCmd[ "titlepos" ]: fTitleSize = dicMainCmd[ "titlepos" ][ "size" ]
  
  if "ymax" in dicMainCmd: 
    fYMax = dicMainCmd[ "ymax" ]
  
  if "xmin" in dicMainCmd: 
    fXMin = dicMainCmd[ "xmin" ]
  
  if "extraText" in dicMainCmd: 
    strExtraText = dicMainCmd[ "extraText" ]
  
  ## Prepare to draw
  
  strCutDef = strCut % myGetDicForFormatVar(dicCutConfig)
  strCutSig = strCutDef + " && muon_signal"
  strCutBkg = strCutDef
  
  dicRef = {"sigeff": fSigEff}
  
  if "bkgcustom" in dicMainCmd:
    dicRef[ "bkgcustom" ] = dicMainCmd[ "bkgcustom" ]
    if "label" in dicRef[ "bkgcustom" ]: 
      strLabelYAxis = dicRef[ "bkgcustom" ][ "label" ]
  
  fBkgMax = 0.0
  
  ## Drawing ROCs
  
  for dicPlotvar in arrPlotvar: 
    plotvar = dicPlotvar[ "plotvar" ]
    
    strSampleSig = SampleSig
    strSampleBkg = SampleBkg
    
    dicCutConfigCurr = copy.deepcopy(dicCutConfig)
    
    if "cutconfig" in dicPlotvar: 
      for strKeyCut in dicPlotvar[ "cutconfig" ]: 
        dicCutConfigCurr[ strKeyCut ] = dicPlotvar[ "cutconfig" ][ strKeyCut ]
    
    if "inputname" in dicPlotvar:
      strSampleSig = SampleSig[ dicPlotvar[ "inputname" ] ]
      strSampleBkg = SampleBkg[ dicPlotvar[ "inputname" ] ]
    
    if "bkgcustom" in dicPlotvar:
      dicRef[ "bkgcustom" ] = dicPlotvar[ "bkgcustom" ]
    
    dicRes = getROC(strSampleSig, strSampleBkg, strMainTree + "/reco", 
      dicPlotvar[ "title" ], binMain, plotvar, strCutSig, strCutBkg, dicRef)
    
    dicPlotvar[ "graph" ] = dicRes[ "graph" ]
    
    setMarkerStyle(dicPlotvar[ "graph" ], dicPlotvar[ "color" ], dicPlotvar[ "shape" ])
    
    dicPlotvar[ "graph" ].GetXaxis().SetLimits(fXMin, 1.1)
    
    if fBkgMax < dicRes[ "bkgmax" ]:
      fBkgMax = dicRes[ "bkgmax" ]
    
    print "%s; %s; %s%s; %0.3f; %0.3f"%(strPrintoutFront, dicPlotvar[ "name" ], myGetVarID(dicCutConfigCurr[ "ID" ]), dicCutConfigCurr[ "PU" ], dicRes[ "performance" ], dicRes[ "bkgrej" ])
  
  #Set canvas
  canv = makeCanvas("canv1", False)
  setMargins(canv, False)
  
  #Legend and drawing
  leg = ROOT.TLegend(dicLegPos[ "left" ], dicLegPos[ "top" ], dicLegPos[ "right" ], dicLegPos[ "bottom" ])
  
  x_name = "Signal efficiency"
  y_name = strLabelYAxis
  
  fYMax = fBkgMax * 1.2
  
  for i, dicPlotvar in enumerate(arrPlotvar): 
    if i == 0: 
      dicPlotvar[ "graph" ].GetXaxis().SetTitle(x_name)
      dicPlotvar[ "graph" ].GetYaxis().SetTitle(y_name)
      dicPlotvar[ "graph" ].GetXaxis().SetTitleSize(0.05)
      dicPlotvar[ "graph" ].GetYaxis().SetTitleOffset(1.20)
      dicPlotvar[ "graph" ].GetYaxis().SetTitleSize(0.05)
      
      dicPlotvar[ "graph" ].GetXaxis().SetLabelSize(0.037)
      dicPlotvar[ "graph" ].GetYaxis().SetLabelSize(0.037)
      
      dicPlotvar[ "graph" ].SetMaximum(fYMax)
      
      dicPlotvar[ "graph" ].Draw("")
    else :
      dicPlotvar[ "graph" ].Draw("same")
    
    strOptLeg = "pl" if "legopt" not in dicPlotvar else dicPlotvar[ "legopt" ]
    
    leg.AddEntry(dicPlotvar[ "graph" ], dicPlotvar[ "graph" ].GetTitle(), strOptLeg)
  
  drawSampleName(strTitleMain % myGetDicForFormatLabel(dicCutConfig), fTitleX, fTitleY, fTitleSize)
  
  leg.SetTextFont(dicLegPos[ "fonttype" ])
  leg.SetTextSize(dicLegPos[ "fontsize" ])
  leg.SetBorderSize(0)
  leg.Draw()
  
  #CMS_lumi setting
  iPos = 0
  iPeriod = 0
  if( iPos==0 ): CMS_lumi.relPosX = 0.12
  CMS_lumi.extraText = strExtraText
  CMS_lumi.lumi_sqrtS = ""
  CMS_lumi.CMS_lumi(canv, iPeriod, iPos)
  
  canv.Modified()
  canv.Update()
  canv.SaveAs(strFilenameOutput)


if __name__ == "__main__": 
  if len(sys.argv) < 2: 
    print "Usage : python roccurve_withjson.py JSON_file"
    print "        or import this into a python program and use drawPlotFromDict()."
    exit(1)

  drawPlotFromDict(json.load(open(sys.argv[ 1 ])))
    

