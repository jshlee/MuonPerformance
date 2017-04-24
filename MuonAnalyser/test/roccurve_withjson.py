import ROOT, copy, array, json, os, sys
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
  fSigEffChosen = 0.0
  fBkgRejChosen = 0.0
  
  for i in range(binning[ 0 ] + 2): 
    # DO NOT confuse this : 0th bin has all underflows, while (binning[0] + 1)-th bin has all overflows
    # Well, we have no underflows, but... for habit for safety.
    nX = i + 0
    
    fSigEff = hSig.Integral(0, nX) / fMSig
    fBkgRej = 1.0 - hBkg.Integral(0, nX) / fMBkg
    
    if "sigeff" in dicRef: 
      if fSigEffPrev <= dicRef[ "sigeff" ] and dicRef[ "sigeff" ] <= fSigEff: 
        nPerf = i
        fSigEffChosen = fSigEff
        fBkgRejChosen = fBkgRej
    elif "bkgrej" in dicRef: 
      if fBkgRej <= dicRef[ "bkgrej" ] and dicRef[ "bkgrej" ] <= fBkgRejPrev: 
        nPerf = i
        fSigEffChosen = fSigEff
        fBkgRejChosen = fBkgRej
    
    arrSigEff.append(fSigEff)
    arrBkgRej.append(fBkgRej)
    
    fSigEffPrev = fSigEff
    fBkgRejPrev = fBkgRej
  
  arrX = array.array("d", arrSigEff)
  arrY = array.array("d", arrBkgRej)
  
  graphROC = ROOT.TGraph(binning[ 0 ] + 2, arrX, arrY)
  graphROC.SetTitle(title)
  
  return {"graph": copy.deepcopy(graphROC), 
    "performance": binning[1] + nPerf * ( binning[2] - binning[1] ) / binning[0], 
    "sigeff": fSigEffChosen, "bkgrej": fBkgRejChosen}


def drawSampleName(samplename, fX, fY, fSizeTex):
  tex2 = ROOT.TLatex()
  
  tex2.SetNDC()
  tex2.SetTextFont(62)
  tex2.SetTextSize(fSizeTex)
  
  for i, strLine in enumerate(samplename.split("\n")): 
    tex2.DrawLatex(fX, fY - i * 1.1 * fSizeTex, strLine)


def drawPlotFromDict(dicMainCmd): 
  ## Reading input dictionary
  ## (Required)
  try: 
    binMain = dicMainCmd[ "binning" ]
    arrPlotvar = dicMainCmd[ "vars" ]
    
    strCut = dicMainCmd[ "cut" ]
    dicCutConfig = dicMainCmd[ "cutconfig" ]
    
    strTitleMain = dicMainCmd[ "title" ]
    
    strSampleSig = dicMainCmd[ "inputfilename_sig" ]
    strSampleBkg = dicMainCmd[ "inputfilename_bkg" ]
    
    strPrintoutFront = dicMainCmd[ "printout_front" ]
    strFilenameOutput = dicMainCmd[ "filename" ]
  except KeyError, strErr:
    print "Error: the JSON file does not contain a required key: %s" % strErr
    exit(1)
  
  ## Reading input dictionary
  ## (Optional)
  fSigEff = 0.0 if "sigeff" not in dicMainCmd else dicMainCmd[ "sigeff" ]
  
  dicLegPos = {
    "left":   0.20, 
    "top":    0.83, 
    "right":  0.45, 
    "bottom": 0.68, 
    
    "fonttype": 62, 
    "fontsize": 0.03, 
  }
  
  fTitleX = 0.14
  fTitleY = 0.85
  fTitleSize = 0.03
  
  if "sigeff" not in dicMainCmd: 
    if "ID" in dicCutConfig: 
      if "Tight" in dicCutConfig[ "ID" ]: 
        fSigEff = 0.95
      elif "Loose" in dicCutConfig[ "ID" ]: 
        fSigEff = 0.98
  
  if "legend" in dicMainCmd: 
    for strKey in dicMainCmd[ "legend" ]: 
      dicLegPos[ strKey ] = dicMainCmd[ "legend" ][ strKey ]
  
  if "titlepos" in dicMainCmd: 
    if "x" in dicMainCmd[ "titlepos" ]: fTitleX = dicMainCmd[ "titlepos" ][ "x" ]
    if "y" in dicMainCmd[ "titlepos" ]: fTitleY = dicMainCmd[ "titlepos" ][ "y" ]
    if "size" in dicMainCmd[ "titlepos" ]: fTitleSize = dicMainCmd[ "titlepos" ][ "size" ]
  
  ## Prepare to draw
  
  strCutDef = strCut % dicCutConfig
  strCutSig = strCutDef + " && muon_signal"
  strCutBkg = strCutDef
  
  dicRef = {"sigeff": fSigEff}
  
  ## Drawing ROCs
  
  for dicPlotvar in arrPlotvar: 
    plotvar = dicPlotvar[ "plotvar" ]
    
    dicRes = getROC(strSampleSig, strSampleBkg, "MuonAnalyser/reco", 
      dicPlotvar[ "title" ], binMain, plotvar, strCutSig, strCutBkg, dicRef)
    
    dicPlotvar[ "graph" ] = dicRes[ "graph" ]
    
    setMarkerStyle(dicPlotvar[ "graph" ], dicPlotvar[ "color" ], dicPlotvar[ "shape" ])
    
    dicPlotvar[ "graph" ].GetXaxis().SetLimits(0.0, 1.1)
    dicPlotvar[ "graph" ].SetMaximum(1.5)
    
    print "%s; %s; %s%s; %0.3f; %0.3f"%(strPrintoutFront, dicPlotvar[ "name" ], dicCutConfig[ "ID" ], dicCutConfig[ "PU" ], dicRes[ "performance" ], dicRes[ "bkgrej" ])
  
  #Set canvas
  canv = makeCanvas("canv1", False)
  setMargins(canv, False)
  
  #Legend and drawing
  leg = ROOT.TLegend(dicLegPos[ "left" ], dicLegPos[ "top" ], dicLegPos[ "right" ], dicLegPos[ "bottom" ])
  
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
    
    strOptLeg = "pl" if "legopt" not in dicPlotvar else dicPlotvar[ "legopt" ]
    
    leg.AddEntry(dicPlotvar[ "graph" ], dicPlotvar[ "graph" ].GetTitle(), strOptLeg)
  
  drawSampleName(strTitleMain % dicCutConfig, fTitleX, fTitleY, fTitleSize)
  
  leg.SetTextFont(dicLegPos[ "fonttype" ])
  leg.SetTextSize(dicLegPos[ "fontsize" ])
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
  canv.SaveAs(strFilenameOutput)


if __name__ == "__main__": 
  if len(sys.argv) < 2: 
    print "Usage : python roccurve_withjson.py JSON_file"
    print "        or import this into a python program and use drawPlotFromDict()."
    exit(1)

  drawPlotFromDict(json.load(open(sys.argv[ 1 ])))
    

