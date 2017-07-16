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
    h2.Divide(h1)
    return copy.deepcopy(h2)

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.DrawLatex(0.18, 0.8, samplename)

datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/'
#datadir = "TenMuExtendedE_"
id = "ME0Muon"
nme0=0
binning_l = [[10,5,105],[8,1.6,2.4],[12,-3,3],[10,5,105],[8,1.6,2.4],[12,-3,3]]

for i, plotvar in enumerate(["genMuon.Pt()", "abs(genMuon.Eta())", "genMuon.Phi()", "recoMuon.Pt()", "abs(recoMuon.Eta())", "recoMuon.Phi()"]):
    #Get histos
    if "genMuon" in plotvar:
        regionCut = "genMuon.Pt()>5&&abs(genMuon.Eta())>2.0&&abs(genMuon.Eta())<2.8"
        h_ph2pu0 = getEff(datadir+"pu0.root", "MuonAnalyser/gen", "PhaseII PU0", binning_l[i], plotvar, regionCut, "%s&&genMuon_is%s"%(regionCut,id))
        h_ph2pu200 = getEff(datadir+"pu200.root", "MuonAnalyser/gen", "PhaseII PU200", binning_l[i], plotvar, regionCut, "%s&&genMuon_is%s"%(regionCut,id)) #&&genMuon_noRecHitGEM>=%d"%(regionCut,id,ngem))
        hlist = [h_ph2pu0, h_ph2pu200]

    if "recoMuon" in plotvar:
        regionCut = "recoMuon.Pt()>5&&abs(recoMuon.Eta())>2.0&&abs(recoMuon.Eta())<2.8"
        h_ph2pu0 = getEff(datadir+"pu0.root", "MuonAnalyser/reco", "PhaseII PU0", binning_l[i], plotvar, "%s&&recoMuon_is%s"%(regionCut,id), "%s&&recoMuon_is%s&&!recoMuon_signal&&recoMuon_noRecHitME0>=%d"%(regionCut,id,nme0))
        h_ph2pu200 = getEff(datadir+"pu200.root", "MuonAnalyser/reco", "PhaseII PU200", binning_l[i], plotvar, "%s&&recoMuon_is%s"%(regionCut,id), "%s&&recoMuon_is%s&&!recoMuon_signal&&recoMuon_noRecHitME0   >=%d"%(regionCut,id,nme0))
        hlist = [h_ph2pu0, h_ph2pu200]

    #Set init histo
    h_init = ROOT.TH1F("","",binning_l[i][0],binning_l[i][1],binning_l[i][2])

    #Set axis
    x_name = "Muon "
    if "Pt"  in plotvar: x_name = x_name+"p_{T}"
    if "Eta" in plotvar: x_name = x_name+"|#eta|"
    if "Phi" in plotvar: x_name = x_name+"#phi"

    y_name = id+" Muon "
    if "genMuon" in plotvar:
        h_init.SetMaximum(1.3)
        h_init.SetMinimum(0.2)
        y_name = y_name+"Efficiency"
    if "recoMuon" in plotvar:
        h_init.SetMaximum(max(h.GetMaximum() for h in hlist)*2.5)
        h_init.GetYaxis().SetLabelSize(0.035)
        h_init.GetYaxis().SetTitleOffset(1.2)
        y_name = y_name+"Fake Rate"
    h_init.GetXaxis().SetTitle(x_name)
    h_init.GetYaxis().SetTitle(y_name)
    h_init.GetYaxis().SetTitleOffset(1)


    ############ Plot design ##############
    #Plot style
    setMarkerStyle(h_ph2pu0, 4, 20) #blue, circle
    setMarkerStyle(h_ph2pu200, 2, 21) #red, square
    #setMarkerStyle(h_run2, 3, 22) #green, triangle

    name = "%s_%s"%(plotvar,id)

    #Set canvas
    canv = makeCanvas(name, False)
    setMargins(canv, False)
    h_init.Draw()
    drawSampleName("Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}, p_{T} > 5 GeV")

    #Legend and drawing
    legTop = ROOT.TLegend(0.6,0.7,0.85,0.85)
    legBot = ROOT.TLegend(0.6,0.2,0.85,0.35)
    if "genMuon"  in plotvar: leg = legBot
    if "recoMuon" in plotvar: leg = legTop
    for h in hlist:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"p")
    leg.SetTextFont(61)
    leg.SetTextSize(0.04)
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
    canv.SaveAs("Efficiencies/"+name+"_nME0%d.png"%nme0)

