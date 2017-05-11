import json, copy, os


def myMakeDir(strDir):
  strDirCurr = ""
  
  for strPart in strDir.split("/"): 
    strDirCurr = ( strDirCurr + "/" + strPart ) if strDirCurr != "" else strPart
    if not os.path.isdir(strDirCurr): os.mkdir(strDirCurr)


dicMain = {
  "cut": "muon.Pt() > %(pT)s && abs(muon.Eta()) < %(Eta)s && muon_is%(ID)s", # " && abs(muon_poszPV0 - muon_poszMuon) < 0.5"
  "binning": [1500,0,6.0], 
  
  "cutconfig": {
    "pT":  "15", 
    "Eta": "2.4", 
    #"ID":  "Tight", 
    #"PU":  "PU200"
  }, 
  
  #"inputfilename_sig": "/cms/scratch/quark2930/Work/muon_upgrade/samples/170411/run_ZMM_PU200_pre4_rereco02_170411.root", 
  #"inputfilename_bkg": "/cms/scratch/quark2930/Work/muon_upgrade/samples/170411/run_QCD_PU200_pre4_rereco02_170411.root", 
  
  "title": "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} and QCD events (%(PU)s), p_{T} > %(pT)s GeV, |#eta| < %(Eta)s, %(ID)s Muon", # "\n|z_{reco} - z_{sim}| < 0.5"
  
  "titlepos": {
    "x": 0.15, 
    "y": 0.85, 
  }, 
  
  #"printout_front": "withvtxcut; rereco", 
  #"filename": "ROCCurve.png", 
  
  "vars": []
}

################################################################################
## Configs for 
## samples used
################################################################################

strPathSamp = "/cms/scratch/quark2930/Work/muon_upgrade/samples/170411/"
strFilenameZMM = strPathSamp + "run_ZMM_PU%s_pre4_rereco02_170411.root"
strFilenameQCD = strPathSamp + "run_QCD_PU%s_pre4_rereco02_170411.root"

dicPUInName = {"PU0": "0", "PU140": "140", "PU200": "200"}
arrID = ["Loose", "Tight"]

################################################################################
## Configs for 
## ploting variables
################################################################################

dicPlotvar = {
  "PUPPIWL": [
    {
      "name": "PUPPIWLALL", 
      "plotvar": "muon_puppiIsoWithLep", 
      "title": "PhaseII 2023D4Timing PUPPI - with lepton, R = 0.4, ALL", 
      "color": 4, 
      "shape": 20
    }, 
    {
      "name": "PUPPIWLCH", 
      "plotvar": "muon_puppiIsoWithLep04ChargedHadron", 
      "title": "PhaseII 2023D4Timing PUPPI - with lepton, R = 0.4, CH", 
      "color": 2, 
      "shape": 21
    }, 
    {
      "name": "PUPPIWLNH", 
      "plotvar": "muon_puppiIsoWithLep04NeutralHadron", 
      "title": "PhaseII 2023D4Timing PUPPI - with lepton, R = 0.4, NH", 
      "color": 3, 
      "shape": 34
    }, 
    {
      "name": "PUPPIWLPH", 
      "plotvar": "muon_puppiIsoWithLep04Photon", 
      "title": "PhaseII 2023D4Timing PUPPI - with lepton, R = 0.4, PH", 
      "color": 6, 
      "shape": 25
    },
  ], 
  
  "PUPPINL": [
    {
      "name": "PUPPINLALL", 
      "plotvar": "muon_puppiIsoWithoutLep", 
      "title": "PhaseII 2023D4Timing PUPPI - without lepton, R = 0.4, ALL", 
      "color": 4, 
      "shape": 20
    }, 
    {
      "name": "PUPPINLCH", 
      "plotvar": "muon_puppiIsoWithoutLep04ChargedHadron", 
      "title": "PhaseII 2023D4Timing PUPPI - without lepton, R = 0.4, CH", 
      "color": 2, 
      "shape": 21
    }, 
    {
      "name": "PUPPINLNH", 
      "plotvar": "muon_puppiIsoWithoutLep04NeutralHadron", 
      "title": "PhaseII 2023D4Timing PUPPI - without lepton, R = 0.4, NH", 
      "color": 3, 
      "shape": 34
    }, 
    {
      "name": "PUPPINLPH", 
      "plotvar": "muon_puppiIsoWithoutLep04Photon", 
      "title": "PhaseII 2023D4Timing PUPPI - without lepton, R = 0.4, PH", 
      "color": 6, 
      "shape": 25
    },
  ], 
  
  "PF": [
    {
      "name": "PFALL", 
      "plotvar": "muon_PFIsolation04", 
      "title": "PhaseII 2023D4Timing PF, R = 0.4, ALL", 
      "color": 4, 
      "shape": 20
    }, 
    {
      "name": "PFCH", 
      "plotvar": "muon_PFIso04ChargedHadronPt", 
      "title": "PhaseII 2023D4Timing PF, R = 0.4, CH", 
      "color": 2, 
      "shape": 21
    }, 
    {
      "name": "PFNH", 
      "plotvar": "muon_PFIso04NeutralHadronEt", 
      "title": "PhaseII 2023D4Timing PF, R = 0.4, NH", 
      "color": 3, 
      "shape": 34
    }, 
    {
      "name": "PFPH", 
      "plotvar": "muon_PFIso04PhotonEt", 
      "title": "PhaseII 2023D4Timing PF, R = 0.4, PH", 
      "color": 6, 
      "shape": 25
    },
    {
      "name": "PFPU", 
      "plotvar": "muon_PFIso04PUPt", 
      "title": "PhaseII 2023D4Timing PF, R = 0.4, PU", 
      "color": 2, 
      "shape": 21
    },
  ], 
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
    "title": ", |z_{reco} - z_{sim}| < 0.5", 
  }
}

################################################################################
## Configs for 
## Isolation value type -- abs or rel
################################################################################

dicAbsRel = {
  "abs": {
    "plotvar": "", 
    "title": " (abs)", 
  }, 
  "rel": {
    "plotvar": " / muon.Pt()", 
    "title": " (rel)", 
  }, 
}

################################################################################
## Looping
################################################################################

strPathJSON = "jsonconf/"
strPathPlot = "plots/"
strPathDate = "20170423/"

for strPU in dicPUInName.keys(): 
  for strPlotvar in dicPlotvar.keys(): 
    for strVtxCut in dicVtxCut.keys(): 
      for strAbsRel in dicAbsRel.keys(): 
        for strID in arrID: 
          strPath = strPathDate + "roccurves_" + strVtxCut + "/" + strPU + "_" + strID + "/"
          strFilename = strPlotvar + ( "_" + strAbsRel if "ALL" not in strPlotvar else "" )
          
          dicCurr = copy.deepcopy(dicMain)
          
          ################################
          ## Putting infos into the dic, except for "filename"
          ################################
          
          dicCurr[ "vars" ] = copy.deepcopy(dicPlotvar[ strPlotvar ])
          dicCurr[ "printout_front" ] = strVtxCut + "; rereco"
          
          dicCurr[ "cutconfig" ][ "ID" ] = strID
          dicCurr[ "cutconfig" ][ "PU" ] = strPU
          
          dicCurr[ "inputfilename_sig" ] = strFilenameZMM % dicPUInName[ strPU ]
          dicCurr[ "inputfilename_bkg" ] = strFilenameQCD % dicPUInName[ strPU ]
          
          for strKey, strVal in dicVtxCut[ strVtxCut ].items():
            dicCurr[ strKey ] += strVal
          
          for dicVar in dicCurr[ "vars" ]: 
            if "ALL" not in dicVar[ "name" ]: 
              for strKey, strVal in dicAbsRel[ strAbsRel ].items():
                dicVar[ strKey ] += strVal
          
          ################################
          ## Making dirs and output; the json file
          ################################
          
          myMakeDir(strPathJSON + strPath)
          myMakeDir(strPathPlot + strPath)
          
          dicCurr[ "filename" ] = strPathPlot + strPath + strFilename + ".png"
          json.dump(dicCurr, open(strPathJSON + strPath + strFilename + ".json", "w"), indent=2)


