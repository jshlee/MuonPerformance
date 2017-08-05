import copy
import json


dicMain = {
  "maintree": "PatMuonAnalyser", 
  "cut": "muon.Pt() > %(pT)s && %(EtaL)s <= abs(muon.Eta()) && abs(muon.Eta()) < %(EtaH)s && ( ( abs(muon.Eta()) < 2.4 && muon_is%(ID)s ) || ( abs(muon.Eta()) > 2.4 && muon_isME0Muon%(ID)s ) ) && abs(muon_poszPV0 - muon_poszMuon) < 0.5", 
  "binning": [1500,0,6.0], 
  
  "cutconfig": {
    "pT":  "0.5", 
    "Eta": "2.4", 
    "PU":  "PU200"
  }, 
  
  #"extracut_ID": {"###bend": "ID", "Loose": "Loose", "Tight": "Tight"}, 
  "extracut_ID": "Loose", 
  "extracut_EtaL": {"###bend": "EtaRange", "00_09": "0.0", "09_16": "0.9", "16_24": "1.6", "00_24": "0.0", "24_28": "2.4"}, 
  "extracut_EtaH": {"###bend": "EtaRange", "00_09": "0.9", "09_16": "1.6", "16_24": "2.4", "00_24": "2.4", "24_28": "2.8"}, 
  
  "inputfilename_sig": "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170804_1_zmm_200.root", 
  "inputfilename_bkg": "/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170804_1_qcd_200.root", 
  
  "title": "Z/#gamma*#rightarrow#font[12]{#mu#mu} and QCD events, #sqrt{s} = 14 TeV, Phase-2 <PU> = 200, %(ID)s Muon, \np_{T} > %(pT)s GeV, %(EtaL)s < |#eta| < %(EtaH)s, |z_{reco} - z_{sim}| < 0.5 cm", 
  
  "titlepos": {
    "x": 0.17, 
    "y": 0.80, 
    "size": 0.033
  }, 
  
  "legend": {
    "top":    0.66, 
    "bottom": 0.46, 
    "fonttype": 62, 
    "fontsize": 0.033
  }, 
  
  "printout_front": "withvtxcut; rereco", 

  "bkgcustom": {
    "name": "bkgyield",
    "label": "Average Loose Muon Bkg Multiplicity",
    "numEvents": 2555040
  }, 
  
  "vars": [], 
  
  "extraText": "Phase-2 Simulation", 
  
  "###jsoninfo": {
      "jsonname": "jsonconf/20170805/roc_various_pt/jsonconf_roccurve_var_pt_Loose_%(Var)s_%(EtaRange)s_20170805.json", 
      "cmd": "python roccurve_withjson.py %s"
  }, 
  
  "filename": {"###savedic": "ROCCurve_TDR_withvtxcut_PU200_var_pt_Loose_%(Var)s_%(EtaRange)s_20170805.png"}
}

listVarMain = [
  {
    "name": "Def", 
    "plotvar": "", 
    "title": "%(Type)s, Pt 0, R = 0.3", 
    "color": 1, 
    "shape": 34
  }, 
  {
    "name": "Pt05", 
    "plotvar": "", 
    "title": "%(Type)s, Pt 5, R = 0.3", 
    "color": 8, 
    "shape": 25
  }, 
  {
    "name": "Pt10", 
    "plotvar": "", 
    "title": "%(Type)s, Pt 10, R = 0.3", 
    "color": 2, 
    "shape": 25
  }, 
  {
    "name": "Pt15", 
    "plotvar": "", 
    "title": "%(Type)s, Pt 15, R = 0.3", 
    "color": 4, 
    "shape": 25
  }, 
  {
    "name": "Pt20", 
    "plotvar": "", 
    "title": "%(Type)s, Pt 20, R = 0.3", 
    "color": 30,
    "shape": 25
  }
]

dicVarType = {
  "PF_ALL": "( muon_pfNewIso%(Pt)s_ch + muon_pfNewIso%(Pt)s_nh + muon_pfNewIso%(Pt)s_ph ) / muon.Pt()", 
  "PF_CH": "muon_pfNewIso%(Pt)s_ch / muon.Pt()", 
  "PF_NH": "muon_pfNewIso%(Pt)s_nh / muon.Pt()", 
  "PF_PH": "muon_pfNewIso%(Pt)s_ph / muon.Pt()", 
  "PUPPI_ALL": "muon_puppiNewIso%(Pt)s", 
  "PUPPI_CH": "muon_puppiNewIso%(Pt)s_ch / muon.Pt()", 
  "PUPPI_NH": "muon_puppiNewIso%(Pt)s_nh / muon.Pt()", 
  "PUPPI_PH": "muon_puppiNewIso%(Pt)s_ph / muon.Pt()", 
}

dicTitleType = {
  "PF_ALL": "I^{PF}", 
  "PF_CH": "I^{PF}, CH", 
  "PF_NH": "I^{PF}, NH", 
  "PF_PH": "I^{PF}, PH", 
  "PUPPI_ALL": "I^{PUPPI}", 
  "PUPPI_CH": "I^{PUPPI}, CH", 
  "PUPPI_NH": "I^{PUPPI}, NH", 
  "PUPPI_PH": "I^{PUPPI}, PH", 
}

dicVar = {"###bend": "Var", }

for strItemType in dicVarType.keys():
    listVarCurr = copy.deepcopy(listVarMain)
    
    for dicItemVar in listVarCurr:
      strPt = dicItemVar[ "name" ]
      if strPt == "Def": strPt = ""
      
      dicItemVar[ "plotvar" ] = dicVarType[ strItemType ]%{"Pt": strPt}
      dicItemVar[ "title" ]   = dicItemVar[ "title" ]%{"Type": dicTitleType[ strItemType ]}
    
    dicVar[ strItemType ] = listVarCurr

dicMain[ "vars" ] = dicVar

json.dump(dicMain, open("evalscr/20170805/template_roc_various_pt.json", "w"), indent=2)


