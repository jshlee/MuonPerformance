import json, copy, os


def myMakeDir(strDir):
  strDirCurr = ""
  
  for strPart in strDir.split("/"): 
    strDirCurr = ( strDirCurr + "/" + strPart ) if strDirCurr != "" else strPart
    if not os.path.isdir(strDirCurr): os.mkdir(strDirCurr)


dicMain = {
    "tree": "reco", 
    "cut": "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s && muon_is%(ID)s", 
    #"plotvar": "muon_PFIsolation04", 
    "binning": [20,0,1.0], 
    
    "cutconfig": {
        "pT":  "15", 
        "Eta": "2.4", 
        #"ID":  "Tight"
    }, 
    
    "title": "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} and QCD events, p_{T} > %(pT)s GeV, |#eta| < %(Eta)s", 
    #"xtitle": "PF Isolation R=0.4", 
    "ytitle": "# of events (normalized)", 
    
    #"filename": "plots/20170405/me0/PF04ALL_Tight_PU200_withvtxcut.png", 
    
    "ylog": 1, 
    "min": 0.00005, 
    "max": 10, 
    
    "vars": ""
}

################################################################################
## Configs for 
## samples used
################################################################################

strPathSamp = "/cms/scratch/quark2930/Work/muon_upgrade/samples/170411/"
strFilenameZMM = "run_ZMM_PU%s_pre4_rereco02_170411.root"
strFilenameQCD = "run_QCD_PU%s_pre4_rereco02_170411.root"

dicPUInName = {"PU0": "0", "PU140": "140", "PU200": "200"}

dicVars = {
  "PU0" : [
    {"title": "Phase II ZMM PU0 - pre4", 
      "filename": strPathSamp + strFilenameZMM % (dicPUInName[ "PU0" ]), 
      "color": 8, 
      "shape": 22, 
      "cut": "muon_signal"
    }, 
    {"title": "Phase II QCD PU0 - pre4", 
      "filename": strPathSamp + strFilenameQCD % (dicPUInName[ "PU0" ]), 
      "color": 2, 
      "shape": 20, 
    }
  ], 
  "PU140" : [
    {"title": "Phase II ZMM PU140 - pre4", 
      "filename": strPathSamp + strFilenameZMM % (dicPUInName[ "PU140" ]), 
      "color": 8, 
      "shape": 22, 
      "cut": "muon_signal"
    }, 
    {"title": "Phase II QCD PU140 - pre4", 
      "filename": strPathSamp + strFilenameQCD % (dicPUInName[ "PU140" ]), 
      "color": 2, 
      "shape": 20, 
    }
  ], 
  "PU200" : [
    {"title": "Phase II ZMM PU200 - pre4", 
      "filename": strPathSamp + strFilenameZMM % (dicPUInName[ "PU200" ]), 
      "color": 8, 
      "shape": 22, 
      "cut": "muon_signal"
    }, 
    {"title": "Phase II QCD PU200 - pre4", 
      "filename": strPathSamp + strFilenameQCD % (dicPUInName[ "PU200" ]), 
      "color": 2, 
      "shape": 20, 
    }
  ], 
}

################################################################################
## Configs for 
## ploting variables
################################################################################

dicPlotvar = {
  "PF04ALL": {
    "plotvar": "muon_PFIsolation04", 
    "xtitle": "PF Isolation R=0.4"
  }, 
  "PF04CH": {
    "plotvar": "muon_PFIso04ChargedHadronPt", 
    "xtitle": "PF Isolation R=0.4, Charged hadron"
  }, 
  "PF04NH": {
    "plotvar": "muon_PFIso04NeutralHadronEt", 
    "xtitle": "PF Isolation R=0.4, Neutral hadron"
  }, 
  "PF04PH": {
    "plotvar": "muon_PFIso04PhotonEt", 
    "xtitle": "PF Isolation R=0.4, Photon"
  }, 
  "PF04PU": {
    "plotvar": "muon_PFIso04PUPt", 
    "xtitle": "PF Isolation R=0.4, Pileup"
  }, 
  "PUPPIWL04ALL": {
    "plotvar": "muon_puppiIsoWithLep", 
    "xtitle": "PUPPI Isolation with lepton R=0.4"
  }, 
  "PUPPIWL04CH": {
    "plotvar": "muon_puppiIsoWithLep04ChargedHadron", 
    "xtitle": "PUPPI Isolation with lepton R=0.4, CH"
  }, 
  "PUPPIWL04NH": {
    "plotvar": "muon_puppiIsoWithLep04NeutralHadron", 
    "xtitle": "PUPPI Isolation with lepton R=0.4, NH"
  }, 
  "PUPPIWL04PH": {
    "plotvar": "muon_puppiIsoWithLep04Photon", 
    "xtitle": "PUPPI Isolation with lepton R=0.4, PH"
  }, 
  "PUPPINL04ALL": {
    "plotvar": "muon_puppiIsoWithoutLep", 
    "xtitle": "PUPPI Isolation without lepton R=0.4"
  }, 
  "PUPPINL04CH": {
    "plotvar": "muon_puppiIsoWithoutLep04ChargedHadron", 
    "xtitle": "PUPPI Isolation without lepton R=0.4, CH"
  }, 
  "PUPPINL04NH": {
    "plotvar": "muon_puppiIsoWithoutLep04NeutralHadron", 
    "xtitle": "PUPPI Isolation without lepton R=0.4, NH"
  }, 
  "PUPPINL04PH": {
    "plotvar": "muon_puppiIsoWithoutLep04Photon", 
    "xtitle": "PUPPI Isolation without lepton R=0.4, PH"
  }, 
}

################################################################################
## Configs for 
## vertex cutting
################################################################################

dicVtxCut = {
  "novtxcut": {
    "cut": "", 
    "title": "", 
  }, 
  "withvtxcut": {
    "cut": " && abs(muon_poszPV0 - muon_poszMuon) < 0.5", 
    "title": "\n|z_{reco} - z_{sim}| < 0.5", 
  }
}

################################################################################
## Configs for 
## Isolation value type -- abs or rel
################################################################################

dicAbsRel = {
  "abs": {
    "plotvar": "", 
    "xtitle": " (abs)", 
  }, 
  "rel": {
    "plotvar": " / muon.Pt()", 
    "xtitle": " (rel)", 
  }, 
}

################################################################################
## Looping
################################################################################

strPathJSON = "jsonconf/"
strPathPlot = "plots/"
strPathDate = "20170423/"

for strVar in dicVars.keys(): 
  for strPlotvar in dicPlotvar.keys(): 
    for strVtxCut in dicVtxCut.keys(): 
      for strAbsRel in dicAbsRel.keys(): 
        for strID in ["Loose", "Tight"]: 
          strPath = strPathDate + "isovalues_" + strVtxCut + "/" + strVar + "_" + strID + "/"
          strFilename = strPlotvar + ( "_" + strAbsRel if "ALL" not in strPlotvar else "" )
          
          dicCurr = copy.deepcopy(dicMain)
          
          ################################
          ## Putting infos into the dic, except for "filename"
          ################################
          
          for strKey, strVal in dicPlotvar[ strPlotvar ].items():
            dicCurr[ strKey ] = strVal
          
          for strKey, strVal in dicVtxCut[ strVtxCut ].items():
            dicCurr[ strKey ] += strVal
          
          if "ALL" not in strPlotvar: 
            for strKey, strVal in dicAbsRel[ strAbsRel ].items():
              dicCurr[ strKey ] += strVal
          
          dicCurr[ "vars" ] = dicVars[ strVar ]
          dicCurr[ "cutconfig" ][ "ID" ] = strID
          
          ################################
          ## Making dirs and output; the json file
          ################################
          
          myMakeDir(strPathJSON + strPath)
          myMakeDir(strPathPlot + strPath)
          
          dicCurr[ "filename" ] = strPathPlot + strPath + strFilename + ".png"
          json.dump(dicCurr, open(strPathJSON + strPath + strFilename + ".json", "w"), indent=2)


