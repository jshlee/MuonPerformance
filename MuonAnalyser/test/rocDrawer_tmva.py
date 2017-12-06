import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def getMarker(filename, totSig, totBkg, cut, style):
    sig = getWeightedEntries(filename, treename, "muon_defaultMva", rangecut+"&&(muon_isSignalMuon)&&%s"%cut)
    bkg = getWeightedEntries(filename, treename, "muon_defaultMva", rangecut+"&&!(muon_isSignalMuon)&&%s"%cut)
    sigeff = sig/totSig
    bkgrej = (totBkg-bkg)/totBkg
    print "marker: ", sigeff, bkgrej
    return ROOT.TMarker(sigeff, bkgrej, style)

def getRocPoint(filename, totSig, totBkg, cut, style, color):
    m = getMarker(filename, totSig, totBkg, cut, style)
    m.SetMarkerSize(3)
    m.SetMarkerColor(color)
    return m

def getRoc(h_sig, h_bkg):
    x=array.array('f',[])
    y=array.array('f',[])
    for k in range(101,0,-1):
        sigeff = h_sig.Integral(k,101)/h_sig.Integral(0,101)
        bkgrej = (h_bkg.Integral(0,101)-h_bkg.Integral(k,101))/h_bkg.Integral(0,101)
        x.append(sigeff)
        y.append(bkgrej)
        print k, h_sig.GetBinCenter(k), sigeff, bkgrej
    return ROOT.TGraph(len(x), x, y)

datadir = './'
#datadir = '/xrootd/store/user/tt8888tt/muon/9_3_0_pre4/tmva/'
treename = "NewPatMuonAnalyser/reco"
rangecut = 'muon.Pt()>5&&abs(muon.Eta())<2.4'
method = "BDT"

samples= {"zmmnew":[0.25,0.4]}
#samples= {"tenmu":[0,-0.06], "ttbar":[-0.02,-0.12], "zmm":[0,-0.06]}
col = [[ROOT.kRed, ROOT.kOrange+7, ROOT.kPink+9], [ROOT.kBlue,ROOT.kViolet+7,ROOT.kAzure+7], [ROOT.kGreen,ROOT.kTeal+9,ROOT.kSpring+9]]

grs = []
marks_tmva = []
marks_id = []
auc = []
for i, sample in enumerate(samples):
    for j, applied in enumerate(samples):
        samplename = sample+"_"+applied+"applied"
        if sample==applied: samplename=sample
        #else: continue 
        filename = datadir+samplename+".root"
        print samplename

        h_sig = makeTH1(filename, treename, "muon_defaultMva", [100,-1,1], "muon_defaultMva", rangecut+"&&muon_isSignalMuon")
        h_bkg = makeTH1(filename, treename, "muon_defaultMva", [100,-1,1], "muon_defaultMva", rangecut+"&&!muon_isSignalMuon")

        gr = getRoc(h_sig, h_bkg) 
        gr.SetLineColor(col[i][j])
        gr.SetTitle("TMVA trained with "+applied.capitalize()+" on "+sample.capitalize())
        grs.append(gr)

        totSig = h_sig.Integral(0,101)
        totBkg = h_bkg.Integral(0,101)

        ##id point
        m1 = getRocPoint(filename, totSig, totBkg, "muon_defaultMva>%f"%samples[sample][0], 22, col[i][j])
        m2 = getRocPoint(filename, totSig, totBkg, "muon_defaultMva>%f"%samples[sample][1], 23, col[i][j])
        marks_tmva.append([m1, m2])

        m_tight = getRocPoint(filename, totSig, totBkg, "muon_isTightCustom", 26, col[i][j])
        m_loose = getRocPoint(filename, totSig, totBkg, "muon_isLoose", 32, col[i][j])
        marks_id.append([m_tight, m_loose])

        #get auc...?
        bdt = ROOT.vector('float')()
        mvaid = ROOT.vector('bool')()
        tfile = ROOT.TFile(filename)
        for t in tfile.Get("NewPatMuonAnalyser/reco"):
            if t.muon.Pt() < 5: continue
            if abs(t.muon.Eta()) > 2.4: continue
            bdt.push_back(t.muon_defaultMva)
            mvaid.push_back(t.muon_isSignalMuon)
        auc = ROOT.TMVA.ROCCurve(bdt, mvaid).GetROCIntegral()
        gr.SetTitle(gr.GetTitle()+" (%.4f)"%auc)
        print gr.GetTitle()


h_init = ROOT.TH1F("","",1000,0.7,1.001)
h_init.SetMaximum(1.001)
h_init.SetMinimum(0.7)
h_init.GetXaxis().SetTitle("Signal Efficiency")
h_init.GetYaxis().SetTitle("Background Rejection")
h_init.GetYaxis().SetTitleOffset(1)

canv = makeCanvas("tmva", False)
canv.SetGrid()
setMargins(canv, False)
h_init.Draw()

#leg = ROOT.TLegend(0.15,0.18,0.4,0.4)
leg = ROOT.TLegend(0.15,0.17,0.7,0.23)

for i, gr in enumerate(grs):
    gr.SetLineWidth(3)
    gr.Draw("Lsame")
    leg.AddEntry(gr,gr.GetTitle(),"l")

    marks_tmva[i][0].Draw("same")
    marks_tmva[i][1].Draw("same")
    marks_id[i][0].Draw("same")
    marks_id[i][1].Draw("same")

leg.SetTextFont(61)
leg.SetTextSize(0.026)
leg.Draw()

iPos = 0
iPeriod = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
CMS_lumi.extraText = "Simulation"
CMS_lumi.lumi_sqrtS = "14 TeV"
CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

canv.Modified()
canv.Update()
canv.SaveAs("roc_tmva.png")

