from MuonPerformance.MuonAnalyser.universaldrawer import drawPlotFromDict

dicInput = {
    "cut": "recoMuon.Pt() > %(pT)s && abs(recoMuon.Eta()) < %(Eta)s && recoMuon_is%(ID)s", 
    "plotvar": "recoMuon_PFIsolation04", 
    "binning": [100,0,0.5], 
    
    "cutconfig": {
        "pT":  "15", 
        "Eta": "2.4", 
        "ID":  "Tight"
    }, 
    
    "title": "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu} and QCD events, p_{T} > %(pT)s GeV, |#eta| < %(Eta)s, %(ID)s Muon", 
    "xtitle": "PF Isolation R=0.4", 
    "ytitle": "# of events (normalized)", 
    
    "filename": "PF04ALL_Tight.png", 
    
    "ylog": 1, 
    "min": 0.00005, 
    "max": 1.0, 
    
    "vars": [
        {"title": "Phase II Signal PU0 - pre4", 
            "filename": "/cms/scratch/quark2930/Work/muon_upgrade/samples/run_ZMM_PU0_pre4_fixed01_test01.root", 
            "color": 8, 
            "shape": 22, 
            "cut": "recoMuon_signal"}, 
        {"title": "Phase II Signal PU200 - pre4", 
            "filename": "/cms/scratch/quark2930/Work/muon_upgrade/samples/run_ZMM_PU200_pre4_fixed01_test01.root", 
            "color": 2, 
            "shape": 20, 
            "cut": "recoMuon_signal"}, 
        {"title": "Phase II QCD PU0 - pre4", 
            "filename": "/cms/scratch/quark2930/Work/muon_upgrade/samples/run_QCD_PU0_pre4_fixed01_test01.root", 
            "color": 1, 
            "shape": 34}, 
        {"title": "Phase II QCD PU200 - pre4", 
            "filename": "/cms/scratch/quark2930/Work/muon_upgrade/samples/run_QCD_PU200_pre4_fixed01_test01.root", 
            "color": 4, 
            "shape": 21}
    ]
}

drawPlotFromDict(dicInput)


