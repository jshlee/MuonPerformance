import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(1.5)
    h.SetLineColor(color)
    h.SetLineWidth(2)

def getEff(filename,treename,title,binning,plotvar,dencut,numcut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,dencut)
    h2 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    hEff = ROOT.TEfficiency(h2,h1)
    hEff.SetTitle(title)
    return copy.deepcopy(hEff)

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.18, 0.85, samplename)

#datadir = '/cms/scratch/quark2930/Work/muon_upgrade/CMSSW_9_0_0_pre4/src/MuonPerformance/MuonAnalyser/test/puppi_ZMM_'
#datadir = '/pnfs/user/tt8888tt/muon/'
datadir = '/pnfs/user/quark2/Work/muon_upgrade/samples/'
#datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
#datadir = "TenMuExtendedE_"

id = sys.argv[1]
binning_l = [[20,5,105],[24,0,2.4],[30,-3,3],[20,5,105],[24,0,2.4],[30,-3,3],[40,60,260]]

for i, plotvar in enumerate(["genMuon.Pt()", "abs(genMuon.Eta())", "genMuon.Phi()", "recoMuon.Pt()", "abs(recoMuon.Eta())", "recoMuon.Phi()", "nvertex"]):
    #Get histos
    if "genMuon" in plotvar:
        h_ph2pu0   = getEff(datadir+"run_ZMM_PU0_pre4_rereco01_ver01.root",   "MuonAnalyser/gen", "PU 0",   binning_l[i], plotvar, "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4", "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4&&genMuon_is%s"%id)
        h_ph2pu140 = getEff(datadir+"run_ZMM_PU140_pre4_rereco01_ver01.root", "MuonAnalyser/gen", "PU 140", binning_l[i], plotvar, "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4", "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4&&genMuon_is%s"%id)
        h_ph2pu200 = getEff(datadir+"run_ZMM_PU200_pre4_rereco01_ver01.root", "MuonAnalyser/gen", "PU 200", binning_l[i], plotvar, "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4", "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4&&genMuon_is%s"%id)
        hlist = [h_ph2pu0, h_ph2pu140, h_ph2pu200]

    #Backgorund rate
    """
    if "recoMuon" in plotvar:
        tfile = ROOT.TFile(datadir+"run_ZMM_PU0_pre4_rereco01_ver01.root")
        nevents_pu0 = tfile.Get("MuonAnalyser/nevents").Integral()
        tfile = ROOT.TFile(datadir+"run_ZMM_PU140_pre4_rereco01_ver01.root")
        nevents_pu140 = tfile.Get("MuonAnalyser/nevents").Integral()
        tfile = ROOT.TFile(datadir+"run_ZMM_PU200_pre4_rereco01_ver01.root")
        nevents_pu200 = tfile.Get("MuonAnalyser/nevents").Integral()
        h_ph2pu0   = makeTH1(datadir+"run_ZMM_PU0_pre4_rereco01_ver01.root",   "MuonAnalyser/reco", "PU 0",   binning_l[i], plotvar, "recoMuon.Pt()>5&&abs(recoMuon.Eta())<2.4&&!recoMuon_signal&&recoMuon_is%s"%id)
        h_ph2pu140 = makeTH1(datadir+"run_ZMM_PU140_pre4_rereco01_ver01.root", "MuonAnalyser/reco", "PU 140", binning_l[i], plotvar, "recoMuon.Pt()>5&&abs(recoMuon.Eta())<2.4&&!recoMuon_signal&&recoMuon_is%s"%id)
        h_ph2pu200 = makeTH1(datadir+"run_ZMM_PU200_pre4_rereco01_ver01.root", "MuonAnalyser/reco", "PU 200", binning_l[i], plotvar, "recoMuon.Pt()>5&&abs(recoMuon.Eta())<2.4&&!recoMuon_signal&&recoMuon_is%s"%id)
        h_ph2pu0.Scale(1/nevents_pu0)
        h_ph2pu140.Scale(1/nevents_pu140)
        h_ph2pu200.Scale(1/nevents_pu200)
        hlist = [h_ph2pu0, h_ph2pu140, h_ph2pu200]
    """
    #Fake rate
    if "recoMuon" in plotvar:
        h_ph2pu0   = getEff(datadir+"run_ZMM_PU0_pre4_rereco01_ver01.root", "MuonAnalyser/reco", "PhaseII PU0", binning_l[i], plotvar, "abs(recoMuon.Eta())<2.4&&recoMuon_is%s"%id, "abs(recoMuon.Eta())<2.4&&!recoMuon_signal&&recoMuon_is%s"%id)
        h_ph2pu140 = getEff(datadir+"run_ZMM_PU140_pre4_rereco01_ver01.root", "MuonAnalyser/reco", "PhaseII PU140", binning_l[i], plotvar, "abs(recoMuon.Eta())<2.4&&recoMuon_is%s"%id, "abs(recoMuon.Eta())<2.4&&!recoMuon_signal&&recoMuon_is%s"%id)
        h_ph2pu200 = getEff(datadir+"run_ZMM_PU200_pre4_rereco01_ver01.root", "MuonAnalyser/reco", "PhaseII PU200", binning_l[i], plotvar, "abs(recoMuon.Eta())<2.4&&recoMuon_is%s"%id, "abs(recoMuon.Eta())<2.4&&!recoMuon_signal&&recoMuon_is%s"%id)
        hlist = [h_ph2pu0, h_ph2pu140, h_ph2pu200]

    if "nvertex" in plotvar:
        h_ph2pu0   = getEff(datadir+"run_ZMM_PU0_pre4_rereco01_ver01.root",   "MuonAnalyser/gen", "PU 0",   binning_l[i], plotvar, "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4", "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4&&genMuon_is%s"%id)
        h_ph2pu140 = getEff(datadir+"run_ZMM_PU140_pre4_rereco01_ver01.root", "MuonAnalyser/gen", "PU 140", binning_l[i], plotvar, "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4", "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4&&genMuon_is%s"%id)
        h_ph2pu200 = getEff(datadir+"run_ZMM_PU200_pre4_rereco01_ver01.root", "MuonAnalyser/gen", "PU 200", binning_l[i], plotvar, "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4", "genMuon.Pt()>5&&abs(genMuon.Eta())<2.4&&genMuon_is%s"%id)
        hlist = [h_ph2pu0, h_ph2pu140, h_ph2pu200]
        

    #Set init histo
    h_init = ROOT.TH1F("","",binning_l[i][0],binning_l[i][1],binning_l[i][2])

    #Set axis
    x_name = "Muon "
    if "Pt"  in plotvar: x_name = x_name+"p_{T}"
    elif "Eta" in plotvar: x_name = x_name+"|#eta|"
    elif "Phi" in plotvar: x_name = x_name+"#phi"
    else : x_name = "Number of vertex"

    y_name = id+" Muon "
    if "Custom" in id:
        y_name = "TightNoVtx Muon "
    if "recoMuon" in plotvar:
        #"""
        maxY = hlist[1].GetEfficiency(0)
        for j in range(binning_l[i][0]):
            x = hlist[1].GetEfficiency(j)
            if maxY > x: continue
            maxY = x
        h_init.SetMaximum(maxY*3)
        #"""
        #h_init.SetMaximum(max(h.GetMaximum() for h in hlist)*2.5)
        h_init.GetYaxis().SetLabelSize(0.035)
        h_init.GetYaxis().SetTitleSize(0.050)
        h_init.GetYaxis().SetTitleOffset(1.2)
        y_name = y_name+"Fake Rate"
        #y_name = y_name+"Background Rate"
    else:
        h_init.SetMaximum(1.15)
        h_init.SetMinimum(0)
        y_name = y_name+"Efficiency"
    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.GetYaxis().SetTitleOffset(1)


    ############ Plot design ##############
    #Plot style
    setMarkerStyle(h_ph2pu0, 4, 20) #blue, circle
    setMarkerStyle(h_ph2pu140, 1, 34) #black, cross
    setMarkerStyle(h_ph2pu200, 2, 21) #red, square
    #setMarkerStyle(h_run2, 3, 22) #green, triangle

    name = "%s_%s"%(plotvar,id)

    #Set canvas
    canv = makeCanvas(name, False)
    setMargins(canv, False)
    h_init.Draw()
    drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, |#eta| < 2.4")
    #drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV, |#eta| < 2.4")

    #Legend and drawing
    legTop = ROOT.TLegend(0.6,0.6,0.8,0.75)
    legBot = ROOT.TLegend(0.6,0.2,0.8,0.35)
    if "genMuon"  in plotvar: leg = legBot
    elif "recoMuon" in plotvar: leg = legTop
    else: leg = legBot
    """
    hlist[1].Draw("e1same")
    hlist[2].Draw("e1same")
    leg.AddEntry(hlist[1],hlist[1].GetTitle(),"p")
    leg.AddEntry(hlist[2],hlist[2].GetTitle(),"p")
    """
    for h in hlist:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"p")
    leg.SetTextFont(61)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
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
    canv.SaveAs(name+".png")

