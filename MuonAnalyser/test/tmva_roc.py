import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def getMarker(datadir, sample, h_sig, h_bkg, cut, style, color):
    sig = getWeightedEntries(datadir+sample+".root", treename, "muon_tmva_bdt", rangecut+"&&muon_signal&&%s"%cut)/h_sig.Integral(0,101)
    bkg = (h_bkg.Integral(0,101)-getWeightedEntries(datadir+sample+".root", treename, "muon_tmva_bdt", rangecut+"&&!muon_signal&&%s"%cut))/h_bkg.Integral(0,101)
    m = ROOT.TMarker(sig, bkg, style)
    m.SetMarkerSize(3)
    m.SetMarkerColor(color)
    return copy.deepcopy(m)

datadir = '/xrootd/store/user/tt8888tt/muon/9_1_1_patch1/tmva/'
#datadir = os.environ["CMSSW_BASE"]+'/src/MuonPerformance/MuonAnalyser/test/results/'
treename = "MuonAnalyser/reco"
rangecut = 'muon.Pt()>5&&abs(muon.Eta())<2.4'
method = "BDT"

tmvacuts= {"tenmu":[0.12,0.038], "ttbar":[0.08,0.02], "zmm":[0.03,-0.07]}
titles=["TMVA trained with Tenmu on Tenmu",
        "TMVA trained with Tenmu on ZMM",
        "TMVA trained with Tenmu on TTbar",
        "TMVA trained with ZMM on ZMM",
        "TMVA trained with ZMM on Tenmu",
        "TMVA trained with ZMM on TTbar",
        "TMVA trained with TTbar on TTbar",
        "TMVA trained with TTbar on ZMM",
        "TMVA trained with TTbar on Tenmu"]

grs = []
hlist = []
markers = []
for i, sample in enumerate(["zmm"]):
#for i, sample in enumerate(["tenmu","zmm_tenmuapplied", "ttbar_tenmuapplied","zmm","tenmu_zmmapplied", "ttbar_zmmapplied","ttbar","zmm_ttbarapplied","tenmu_ttbarapplied"]):
#for i, sample in enumerate(["zmm", "ttbar","tenmu"]):
    tfile = ROOT.TFile(datadir+sample+".root")
    nevents = tfile.Get("MuonAnalyser/nevents").Integral()

    print sample
    h_sig = makeTH1(datadir+sample+".root", treename, "muon_tmva_bdt", [100,-0.25,0.25], "muon_tmva_bdt", rangecut+"&&muon_signal")
    h_bkg = makeTH1(datadir+sample+".root", treename, "muon_tmva_bdt", [100,-0.25,0.25], "muon_tmva_bdt", rangecut+"&&!muon_signal")

    # output result
    x=array.array('f',[])
    y=array.array('f',[])
    for j in range(101,0,-1):
        sigeff = h_sig.Integral(j,101)/h_sig.Integral(0,101)
        bkgrej = (h_bkg.Integral(0,101)-h_bkg.Integral(j,101))/h_bkg.Integral(0,101)
        x.append(sigeff)
        y.append(bkgrej)
        fake = h_bkg.Integral(j,101)/nevents
        print j, h_sig.GetBinCenter(j), sigeff, bkgrej, fake
        #if abs(h_sig.GetBinCenter(j)-0.05) <0.005: print "<<<<<<"
    gr = ROOT.TGraph(len(x), x, y)
    gr.SetLineColor(i+1)
    #gr.SetTitle(sample.split('.')[0])
    gr.SetTitle(titles[i])
    grs.append(gr)

    # tmva result
    #tfile = ROOT.TFile("./tmva/TMVA_"+sample+".root")
    #htmp = tfile.Get("dataset/Method_%s/%s/MVA_%s_rejBvsS"%(method,method,method))
    #htmp.SetLineColor(i+2)
    #htmp.SetTitle(sample.split('/')[-1].split('.')[0]+" (TMVA)")
    #hlist.append(copy.deepcopy(htmp))

    # id marker
    m_tightid = getMarker(datadir, sample, h_sig, h_bkg, "muon_isTightCustom", 22, i+2)
    m_looseid = getMarker(datadir, sample, h_sig, h_bkg, "muon_isLoose", 23, i+2)
    m_tighttmva = getMarker(datadir, sample, h_sig, h_bkg, "muon_tmva_bdt>%f"%tmvacuts[sample][0], 26, i+2)
    m_loosetmva = getMarker(datadir, sample, h_sig, h_bkg, "muon_tmva_bdt>%f"%tmvacuts[sample][1], 32, i+2)
    print m_tightid.GetX(), m_tighttmva.GetX(), m_looseid.GetX(), m_loosetmva.GetX()
    print m_tightid.GetY(), m_tighttmva.GetY(), m_looseid.GetY(), m_loosetmva.GetY()
    markers.append([m_tightid, m_tighttmva, m_looseid, m_loosetmva])
    marker_titles=["Tight ID", "Tight ID (TMVA)", "Loose ID", "Loose ID (TMVA)"]
    #markers.append([m_tightid,m_looseid)

h_init = ROOT.TH1F("","",1000,0.85,1.01)
h_init.SetMaximum(1.01)
h_init.SetMinimum(0.85)
h_init.GetXaxis().SetTitle("Signal Efficiency")
h_init.GetYaxis().SetTitle("Background Rejection")
h_init.GetYaxis().SetTitleOffset(1)

canv = makeCanvas("tmva", False)
canv.SetGrid()
setMargins(canv, False)
h_init.Draw()

leg = ROOT.TLegend(0.15,0.18,0.65,0.55)
#leg = ROOT.TLegend(0.15,0.18,0.5,0.48)
#for ms in markers:
#    for j, m in enumerate(ms):
#        m.Draw()
#        leg.AddEntry(m,marker_titles[j],"p")

for i, gr in enumerate(grs):
    gr.SetLineWidth(3)
    gr.Draw("Lsame")
    leg.AddEntry(gr,gr.GetTitle(),"l")

    #hlist[i].SetLineStyle(2)
    #hlist[i].SetLineWidth(3)
    #hlist[i].Draw("same")
    #leg.AddEntry(hlist[i],hlist[i].GetTitle(),"l")

#leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(61)
leg.SetTextSize(0.03)
leg.Draw()

iPos = 0
iPeriod = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
CMS_lumi.extraText = "Simulation"
CMS_lumi.lumi_sqrtS = "14 TeV"
CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

canv.Modified()
canv.Update()
canv.SaveAs("tmva_roc.png")
#canv.SaveAs("tmva_roc_%s.png"%sample)

