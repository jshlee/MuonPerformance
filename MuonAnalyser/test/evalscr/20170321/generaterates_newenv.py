import sys, json


#strOTP = "OT_"
strOTP = ""

dicVars = {
  "rereco": {
    "gen": [
      {"title": "Phase II Signal PU0 - pre4", 
        "filename": "/cms/scratch/gwheo/muonPerf_900_pre6/src/MuonPerformance/MuonAnalyser/test/ZMM_PU0_pre6/ZMM_PU0_pre6.root", 
        "color": 4, 
        "shape": 20}, 
      {"title": "Phase II Signal PU140 - pre4", 
        "filename": "/cms/scratch/gwheo/muonPerf_900_pre6/src/MuonPerformance/MuonAnalyser/test/ZMM_PU140_pre6/ZMM_PU140_pre6.root",
        "color": 1, 
        "shape": 34}, 
      {"title": "Phase II Signal PU200 - pre4", 
        "filename": "/cms/scratch/gwheo/muonPerf_900_pre6/src/MuonPerformance/MuonAnalyser/test/ZMM_PU200_pre6/ZMM_PU200_pre6.root",
        "color": 2, 
        "shape": 21}, 
    ], 
    "reco": [
      {"title": "Phase II QCD PU0 - pre4", 
        "filename": "/cms/scratch/gwheo/muonPerf_900_pre6/src/MuonPerformance/MuonAnalyser/test/QCD_PU0_pre6/QCD_PU0_pre6.root", 
        "color": 4, 
        "shape": 20}, 
      {"title": "Phase II QCD PU140 - pre4", 
        "filename": "/cms/scratch/gwheo/muonPerf_900_pre6/src/MuonPerformance/MuonAnalyser/test/QCD_PU140_pre6/QCD_PU140_pre6.root",
        "color": 1, 
        "shape": 34}, 
      {"title": "Phase II QCD PU200 - pre4", 
        "filename": "/cms/scratch/gwheo/muonPerf_900_pre6/src/MuonPerformance/MuonAnalyser/test/QCD_PU200_pre6/QCD_PU200_pre6.root",
        "color": 2, 
        "shape": 21}
    ]
  }, 
}

dicIsoCutValAll = {}

fIsoCutVal = open("isovaluecutlist.txt")

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

dicIsoCutValNoVtx = {
  "relval": {
    "PF":      {
      "LoosePU0": "0.316",
      "LoosePU140": "0.772",
      "LoosePU200": "1",
      "LooseModPU0": "0.324",
      "LooseModPU140": "1.08",
      "LooseModPU200": "1.424",
      "TightPU0": "0.16",
      "TightPU140": "0.472",
      "TightPU200": "0.64",
      "TightModNoIPPU0": "0.164",
      "TightModNoIPPU140": "0.78",
      "TightModNoIPPU200": "1.044",
    }, 
    "PUPPIWL": {
      "LoosePU0": "0.324",
      "LoosePU140": "0.636",
      "LoosePU200": "0.652",
      "LooseModPU0": "0.332",
      "LooseModPU140": "0.748",
      "LooseModPU200": "0.752",
      "TightPU0": "0.168",
      "TightPU140": "0.432",
      "TightPU200": "0.456",
      "TightModNoIPPU0": "0.172",
      "TightModNoIPPU140": "0.52",
      "TightModNoIPPU200": "0.516",
    }, 
    "PUPPINL": {
      "LoosePU0": "0.292",
      "LoosePU140": "0.336",
      "LoosePU200": "0.352",
      "LooseModPU0": "0.3",
      "LooseModPU140": "0.368",
      "LooseModPU200": "0.372",
      "TightPU0": "0.152",
      "TightPU140": "0.188",
      "TightPU200": "0.196",
      "TightModNoIPPU0": "0.152",
      "TightModNoIPPU140": "0.212",
      "TightModNoIPPU200": "0.22",
    },  
    "PUPPICB": {
      "LoosePU0": "0.308",
      "LoosePU140": "0.448",
      "LoosePU200": "0.464",
      "LooseModPU0": "0.316",
      "LooseModPU140": "0.516",
      "LooseModPU200": "0.516",
      "TightPU0": "0.16",
      "TightPU140": "0.288",
      "TightPU200": "0.304",
      "TightModNoIPPU0": "0.164",
      "TightModNoIPPU140": "0.34",
      "TightModNoIPPU200": "0.344",
    },  
  }, 
  "rereco": {
    "PF":      {
      "LoosePU0": "0.328",
      "LoosePU140": "1.06",
      "LoosePU200": "1.328",
      "LooseModPU0": "0.324",
      "LooseModPU140": "1.06",
      "LooseModPU200": "1.328",
      "TightPU0": "0.164",
      "TightPU140": "0.732",
      "TightPU200": "0.94",
      "TightModNoIPPU0": "0.164",
      "TightModNoIPPU140": "0.76",
      "TightModNoIPPU200": "0.968",
    }, 
    "PUPPIWL": {
      "LoosePU0": "0.332",
      "LoosePU140": "0.752",
      "LoosePU200": "0.756",
      "LooseModPU0": "0.332",
      "LooseModPU140": "0.752",
      "LooseModPU200": "0.756",
      "TightPU0": "0.172",
      "TightPU140": "0.528",
      "TightPU200": "0.548",
      "TightModNoIPPU0": "0.172",
      "TightModNoIPPU140": "0.524",
      "TightModNoIPPU200": "0.52",
    }, 
    "PUPPINL": {
      "LoosePU0": "0.3",
      "LoosePU140": "0.376",
      "LoosePU200": "0.384",
      "LooseModPU0": "0.3",
      "LooseModPU140": "0.376",
      "LooseModPU200": "0.384",
      "TightPU0": "0.152",
      "TightPU140": "0.212",
      "TightPU200": "0.22",
      "TightModNoIPPU0": "0.152",
      "TightModNoIPPU140": "0.216",
      "TightModNoIPPU200": "0.228",
    },  
    "PUPPICB": {
      "LoosePU0": "0.316",
      "LoosePU140": "0.528",
      "LoosePU200": "0.524",
      "LooseModPU0": "0.312",
      "LooseModPU140": "0.524",
      "LooseModPU200": "0.524",
      "TightPU0": "0.16",
      "TightPU140": "0.344",
      "TightPU200": "0.36",
      "TightModNoIPPU0": "0.16",
      "TightModNoIPPU140": "0.344",
      "TightModNoIPPU200": "0.352",
    }, 
  }
}

dicIsoCutValWithVtx = {
  "relval": {
    "PF":      {
      "LoosePU0": "0.316",
      "LoosePU140": "0.74",
      "LoosePU200": "0.968",
      "LooseModPU0": "0.324",
      "LooseModPU140": "1.032",
      "LooseModPU200": "1.352",
      "TightPU0": "0.16",
      "TightPU140": "0.472",
      "TightPU200": "0.64",
      "TightModNoIPPU0": "0.164",
      "TightModNoIPPU140": "0.756",
      "TightModNoIPPU200": "1.012",
    }, 
    "PUPPIWL": {
      "LoosePU0": "0.324",
      "LoosePU140": "0.628",
      "LoosePU200": "0.668",
      "LooseModPU0": "0.332",
      "LooseModPU140": "0.744",
      "LooseModPU200": "0.764",
      "TightPU0": "0.168",
      "TightPU140": "0.432",
      "TightPU200": "0.456",
      "TightModNoIPPU0": "0.172",
      "TightModNoIPPU140": "0.524",
      "TightModNoIPPU200": "0.54",
    }, 
    "PUPPINL": {
      "LoosePU0": "0.292",
      "LoosePU140": "0.312",
      "LoosePU200": "0.336",
      "LooseModPU0": "0.3",
      "LooseModPU140": "0.356",
      "LooseModPU200": "0.36",
      "TightPU0": "0.152",
      "TightPU140": "0.188",
      "TightPU200": "0.196",
      "TightModNoIPPU0": "0.152",
      "TightModNoIPPU140": "0.204",
      "TightModNoIPPU200": "0.212",
    },  
    "PUPPICB": {
      "LoosePU0": "0.308",
      "LoosePU140": "0.432",
      "LoosePU200": "0.468",
      "LooseModPU0": "0.316",
      "LooseModPU140": "0.508",
      "LooseModPU200": "0.52",
      "TightPU0": "0.16",
      "TightPU140": "0.288",
      "TightPU200": "0.304",
      "TightModNoIPPU0": "0.164",
      "TightModNoIPPU140": "0.34",
      "TightModNoIPPU200": "0.352",
    },  
  }, 
  "rereco": {
    "PF":      {
      "LoosePU0": "0.328",
      "LoosePU140": "1.012",
      "LoosePU200": "1.268",
      "LooseModPU0": "0.324",
      "LooseModPU140": "1.012",
      "LooseModPU200": "1.268",
      "TightPU0": "0.164",
      "TightPU140": "0.732",
      "TightPU200": "0.94",
      "TightModNoIPPU0": "0.164",
      "TightModNoIPPU140": "0.732",
      "TightModNoIPPU200": "0.940",
    }, 
    "PUPPIWL": {
      "LoosePU0": "0.332",
      "LoosePU140": "0.748",
      "LoosePU200": "0.772",
      "LooseModPU0": "0.332",
      "LooseModPU140": "0.748",
      "LooseModPU200": "0.772",
      "TightPU0": "0.172",
      "TightPU140": "0.528",
      "TightPU200": "0.548",
      "TightModNoIPPU0": "0.172",
      "TightModNoIPPU140": "0.524",
      "TightModNoIPPU200": "0.548",
    }, 
    "PUPPINL": {
      "LoosePU0": "0.3",
      "LoosePU140": "0.36",
      "LoosePU200": "0.372",
      "LooseModPU0": "0.3",
      "LooseModPU140": "0.36",
      "LooseModPU200": "0.372",
      "TightPU0": "0.152",
      "TightPU140": "0.212",
      "TightPU200": "0.22",
      "TightModNoIPPU0": "0.152",
      "TightModNoIPPU140": "0.212",
      "TightModNoIPPU200": "0.22",
    },  
    "PUPPICB": {
      "LoosePU0": "0.316",
      "LoosePU140": "0.516",
      "LoosePU200": "0.528",
      "LooseModPU0": "0.312",
      "LooseModPU140": "0.516",
      "LooseModPU200": "0.528",
      "TightPU0": "0.16",
      "TightPU140": "0.344",
      "TightPU200": "0.36",
      "TightModNoIPPU0": "0.16",
      "TightModNoIPPU140": "0.344",
      "TightModNoIPPU200": "0.36",
    }, 
  }
}

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
  {"var": "muon.Pt()",      "legend": [0.15, 0.63, 0.50, 0.78]}, 
  {"var": "abs(muon.Eta())", "legend": [0.15, 0.20, 0.50, 0.40]}, 
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
                "Eta": "2.4", 
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
            
            strDir = "20170321/effbkgplots_" + strIsVtxCut + "/" + strRelRe + "/" + strID + "/" + strDest
            strFilename = dicX[ "prefix" ] + "_" + strIsoType
            
            dicOutput[ "filename" ] = "plots/" + strDir + "/" + strFilename + ".png"
            strJSONFile = strDir + "/" + strFilename + ".json"
            
            for dicLegend in arrLegend: 
              if dicLegend[ "var" ] == dicOutput[ "plotvar" ]: 
                dicNewLegend = {}
                
                dicNewLegend[ "left" ]   = dicLegend[ "legend" ][ 0 ]
                dicNewLegend[ "top" ]    = dicLegend[ "legend" ][ 1 ]
                dicNewLegend[ "right" ]  = dicLegend[ "legend" ][ 2 ]
                dicNewLegend[ "bottom" ] = dicLegend[ "legend" ][ 3 ]
                
                dicOutput[ "legend" ] = dicNewLegend
            
            if "min" in dicX: 
              dicOutput[ "min" ] = dicX[ "min" ]
            
            if os.path.isdir("jsonconf"): os.mkdir("jsonconf")
            if os.path.isdir("plots"): os.mkdir("plots")
            
            strDirCurr = ""
            
            for strPart in strDir.split("/"): 
              strDirCurr = strDirCurr + "/" + strPart
              if os.path.isdir("jsonconf" + strDirCurr): os.mkdir("jsonconf" + strDirCurr)
              if os.path.isdir("plots" + strDirCurr): os.mkdir("plots" + strDirCurr)
            
            json.dump(dicOutput, open("jsonconf/" + strJSONFile, "w"))


