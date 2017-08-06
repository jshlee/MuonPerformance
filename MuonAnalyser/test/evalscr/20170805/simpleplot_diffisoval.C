{
    // "..._170804_1_qcd_..." can be replaced by "..._170804_1_zmm_..."
    TFile *f = new TFile("/xrootd/store/user/quark2930/muon_upgrade/TDRSpring2017/tdr_170804_1_qcd_200.root");
    TTree *tt = (TTree *)f->Get("PatMuonAnalyser/reco");
    
    // you can check other differences by changing that value
    const char *strVal = "(muon_puppiNewIsoPt10_nh - muon_puppiNewIso_nh)";
    int nNumBin = 100;
    float fMinBin = -5.0;
    float fMaxBin =  5.0;
    
    TH1 *h1 = new TH1D("nametmp1", "muon_puppiNewIsoPt10_nh - muon_puppiNewIso_nh (QCD)", nNumBin, fMinBin, fMaxBin);
    tt->Project("nametmp1", strVal, "muon.Pt() > 0.5 && 0.0 <= abs(muon.Eta()) && abs(muon.Eta()) < 2.4 && ( ( abs(muon.Eta()) < 2.4 && muon_isLoose ) || ( abs(muon.Eta()) > 2.4 && muon_isME0MuonLoose ) ) && abs(muon_poszPV0 - muon_poszMuon) < 0.5");
    
    h1->Draw("H");
}


