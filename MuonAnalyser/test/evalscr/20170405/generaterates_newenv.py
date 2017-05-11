###############################################################################
## 
## 20170405 : more eta (0 ~ 2.8) region
## 
###############################################################################


import sys, os, json


#strOTP = "OT_"
strOTP = ""

dicVars = {
  "rereco": {
    "gen": [
      {"title": "Phase II Signal PU0 - pre4", 
        "filename": "/xrootd/store/user/tt8888tt/muon/zmm.root", 
        "color": 4, 
        "shape": 20}, 
      {"title": "Phase II Signal PU140 - pre4", 
        "filename": "/xrootd/store/user/tt8888tt/muon/zmm140.root",
        "color": 1, 
        "shape": 34}, 
      {"title": "Phase II Signal PU200 - pre4", 
        "filename": "/xrootd/store/user/tt8888tt/muon/zmm200.root",
        "color": 2, 
        "shape": 21}, 
    ], 
    "reco": [
      {"title": "Phase II QCD PU0 - pre4", 
        "filename": "/xrootd/store/user/tt8888tt/muon/qcd.root", 
        "color": 4, 
        "shape": 20}, 
      {"title": "Phase II QCD PU140 - pre4", 
        "filename": "/xrootd/store/user/tt8888tt/muon/qcd140.root",
        "color": 1, 
        "shape": 34}, 
      {"title": "Phase II QCD PU200 - pre4", 
        "filename": "/xrootd/store/user/tt8888tt/muon/qcd200.root",
        "color": 2, 
        "shape": 21}
    ]
  }, 
}

dicIsoCutValAll = {}

fIsoCutVal = open("isovaluecutlist_moreeta.txt")

#for strLine:
while True:
  strLine = fIsoCutVal.readline()
  if not strLine: break
  
  arrIsoVal = strLine.split("; ")
  
  strVtxCut = arrIsoVal[ 0 ].encode("ascii", "ignore")
  strDataType = arrIsoVal[ 1 ].encode("ascii", "ignore")
  strIsoType = arrIsoVal[ 2 ].encode("ascii", "ignore")
  strIDPU = arrIsoVal[ 3 ].encode("ascii", "ignore")
  strIsoVal = arrIsoVal[ 4 ].encode("ascii", "ignore").split("\n")[0]
  
  if strVtxCut   not in dicIsoCutValAll: 
    dicIsoCutValAll[ strVtxCut ] = {}
  
  if strDataType not in dicIsoCutValAll[ strVtxCut ]: 
    dicIsoCutValAll[ strVtxCut ][ strDataType ] = {}
  
  if strIsoType  not in dicIsoCutValAll[ strVtxCut ][ strDataType ]: 
    dicIsoCutValAll[ strVtxCut ][ strDataType ][ strIsoType ] = {}
  
  dicIsoCutValAll[ strVtxCut ][ strDataType ][ strIsoType ][ strIDPU ] = strIsoVal

dicIsoCutPre = {
  "PF":         "muon_PFIsolation04 < ", 
  "PUPPIWL":    "muon_puppiIsoWithLep < ", 
  "PUPPINL":    "muon_puppiIsoWithoutLep < ", 
  "PUPPICB":    "muon_puppiIsoCombined < ", 
  "PUPPINEWWL": "muon_puppiIso < ", 
  "PUPPINEWNL": "muon_puppiIsoNoLep < ", 
}

dicIsoFile = {
  "PF":         "PF", 
  "PUPPIWL":    "PUPPIWithLep", 
  "PUPPINL":    "PUPPIWithoutLep", 
  "PUPPICB":    "PUPPICombined", 
  "PUPPINEWWL": "PUPPINew", 
  "PUPPINEWNL": "PUPPINewNoLep", 
}

dicIsoTitle = {
  "PF":         "%(ID)s Muon, cut by PF isolation", 
  "PUPPIWL":    "%(ID)s Muon, cut by PUPPI isolation with lepton", 
  "PUPPINL":    "%(ID)s Muon, cut by PUPPI isolation without lepton", 
  "PUPPICB":    "%(ID)s Muon, cut by combined PUPPI isolation", 
  "PUPPINEWWL": "%(ID)s Muon, cut by new PUPPI isolation with lepton", 
  "PUPPINEWNL": "%(ID)s Muon, cut by new PUPPI isolation without lepton", 
}

arrLegend = [
  {"tree": "gen",  "var": "muon.Pt()",       "legend": [0.15, 0.63, 0.50, 0.78]}, 
  {"tree": "reco", "var": "abs(muon.Eta())", "legend": [0.15, 0.20, 0.50, 0.40]}, 
]

arrX = [
  {"var": "muon.Pt()", "bin": [20, 5, 105],  "axis": "p_{T}", "prefix": "Pt"}, 
  {"var": "abs(muon.Eta())", "bin": [25, 0, 2.5],  "axis": "|#eta|", "prefix": "Eta"}, 
  {"var": "muon.Phi()", "bin": [32, 0, 3.2], "axis": "#phi",  "prefix": "Phi"}, 
  {"var": "nvertex", "bin": [60, 0, 240], "axis": "Numer of vertex", "prefix": "NV", "min": 0}, 
]

#for strRelRe in ["relval", "rereco"]: 
for strIsVtxCut in ["novtxcut", "withvtxcut"]: 
  for strRelRe in ["rereco"]: 
    for strID in ["Loose", "Tight", "LooseMod", "TightModNoIP"]: 
      for strTree in ["gen", "reco"]: 
        #for strIsoType in ["PF", "PUPPIWL", "PUPPINL", "PUPPICB"]: 
        for strIsoType in ["PF", "PUPPIWL", "PUPPINL", "PUPPINEWWL", "PUPPINEWNL"]: 
          for dicX in arrX: 
            dicIsoCutVal = dicIsoCutValAll[ strIsVtxCut ]
            
            fMin = 0.0
            fMax = 0.0
            strYAxis = ""
            legend = {}
            
            if strTree == "gen": 
              fMin = 0.30
              fMax = 1.20
              
              strYAxis = "Efficiency"
              
              legend = {
                "left":   0.50, 
                "top":    0.20, 
                "right":  0.85, 
                "bottom": 0.40
              }
            else:
              fMin = 0.00
              fMax = 0.35
              
              strYAxis = "Background Rate"
              
              legend = {
                "left":   0.50, 
                "top":    0.63, 
                "right":  0.85, 
                "bottom": 0.78
              }

            dicInput = {
              #"cut": "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s && muon_is%(ID)s && abs(muon_poszPV0 - muon_poszSimPV) < 0.5", 
              "cut": "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s", 
              "plotvar": "muon_PFIsolation04", 
              "binning": [100,0,0.5], 
              
              "cutconfig": {
                "pT":  "15", 
                "Eta": "2.8", 
                "ID":  "Loose"
              }, 
              
              "title": "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} and QCD events, p_{T} > %(pT)s GeV, |#eta| < %(Eta)s", 
              "ytitle": strYAxis, 
              
              "filename": "PF04ALL_Tight.png", 
              
              #"effrate": "", 
              
              #"ylog": 1, 
              "min": fMin, 
              "max": fMax, 
              
              "legend": legend, 
              
              "vars": ""
            }

            strDest = ""
            strKeyRate = ""

            if strTree == "gen": 
              strDest = "eff"
              strKeyRate = "effrate"
            else: 
              strDest = "bkgrate"
              strKeyRate = "bkgrate"

            dicOutput = {}
            
            for strKey in dicInput.keys(): 
              dicOutput[ strKey ] = dicInput[ strKey ]
            
            dicOutput[ "plotvar" ] = dicX[ "var" ]
            dicOutput[ "binning" ] = dicX[ "bin" ]
            
            dicOutput[ "tree" ] = strTree
            dicOutput[ "cutconfig" ][ "ID" ] = strID
            
            if strIsVtxCut == "withvtxcut": 
              dicOutput[ "cut" ] = dicOutput[ "cut" ] + " && abs(muon_poszPV0 - muon_poszSimPV) < 0.5"
            
            dicOutput[ strKeyRate ] = "muon_is%(ID)s"
            
            dicOutput[ "vars" ] = dicVars[ strRelRe ][ strTree ]
            
            for arrRoots in dicOutput[ "vars" ]: 
              strPU = ""
              
              if "PU0" in arrRoots[ "title" ]: 
                strPU = "PU0"
              elif "PU140" in arrRoots[ "title" ]: 
                strPU = "PU140"
              elif "PU200" in arrRoots[ "title" ]: 
                strPU = "PU200"
              
              arrRoots[ strKeyRate ] = dicIsoCutPre[ strIsoType ] + dicIsoCutVal[ strRelRe ][ strIsoType ][ strID + strPU ]
            
            if strTree == "reco": 
              dicOutput[ strKeyRate ] = dicOutput[ strKeyRate ] + " && !muon_signal"
            
            dicOutput[ "xtitle" ] = dicX[ "axis" ]
            
            dicOutput[ "title" ] = dicOutput[ "title" ] + "\n"
            if strIsVtxCut == "withvtxcut": 
              dicOutput[ "title" ] = dicOutput[ "title" ] + "|z_{reco} - z_{sim}| < 0.5, "
            dicOutput[ "title" ] = dicOutput[ "title" ] + dicIsoTitle[ strIsoType ]
            
            strDirVtx = "novtxcut"
            
            if strIsVtxCut == "withvtxcut": 
              strDirVtx = "withvtxcut"
            
            if strIsVtxCut == "withvtxcut" and strTree == "reco": 
              dicOutput[ "bkgdenominator" ] = {"name": "vertex reco vs sim", "min": -0.5, "max": 0.5}
            
            strDir = "20170405/effbkgplots_" + strIsVtxCut + "/" + strRelRe + "/" + strID + "/" + strDest
            strFilename = dicX[ "prefix" ] + "_" + strIsoType
            
            dicOutput[ "filename" ] = "plots/" + strDir + "/" + strFilename + ".png"
            strJSONFile = strDir + "/" + strFilename + ".json"
            
            for dicLegend in arrLegend: 
              if dicLegend[ "tree" ] == strTree and dicLegend[ "var" ] == dicOutput[ "plotvar" ]: 
                dicNewLegend = {}
                
                dicNewLegend[ "left" ]   = dicLegend[ "legend" ][ 0 ]
                dicNewLegend[ "top" ]    = dicLegend[ "legend" ][ 1 ]
                dicNewLegend[ "right" ]  = dicLegend[ "legend" ][ 2 ]
                dicNewLegend[ "bottom" ] = dicLegend[ "legend" ][ 3 ]
                
                dicOutput[ "legend" ] = dicNewLegend
            
            if "min" in dicX: 
              dicOutput[ "min" ] = dicX[ "min" ]
            
            if not os.path.isdir("jsonconf"): os.mkdir("jsonconf")
            if not os.path.isdir("plots"): os.mkdir("plots")
            
            strDirCurr = ""
            
            for strPart in strDir.split("/"): 
              strDirCurr = strDirCurr + "/" + strPart
              if not os.path.isdir("jsonconf" + strDirCurr): os.mkdir("jsonconf" + strDirCurr)
              if not os.path.isdir("plots" + strDirCurr): os.mkdir("plots" + strDirCurr)
            
            json.dump(dicOutput, open("jsonconf/" + strJSONFile, "w"))


