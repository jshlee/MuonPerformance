import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def getMarker(filename, sig_tot, bkg_tot, cut, style):
    sig = getWeightedEntries(filename, "MuonAnalyser/gen", "muon_tmva_bdt", rangecut+"&&(muon_signal)&&%s"%cut)
    bkg = getWeightedEntries(filename, "MuonAnalyser/reco", "muon_tmva_bdt", rangecut+"&&!(muon_signal)&&%s"%cut)
    sigeff = sig/sig_tot
    bkgrej = (bkg_tot-bkg)/bkg_tot
    return ROOT.TMarker(sigeff, bkgrej, style)

def getRoc(h_sig, h_bkg):
    x=array.array('f',[])
    y=array.array('f',[])
    for k in range(101,0,-1):
        sigeff = h_sig.Integral(k,101)/h_sig.Integral(0,101)
        bkgrej = (h_bkg.Integral(0,101)-h_bkg.Integral(k,101))/h_bkg.Integral(0,101)
        x.append(sigeff)
        y.append(bkgrej)
        #print k, h_sig.GetBinCenter(k), sigeff, bkgrej
    return ROOT.TGraph(len(x), x, y)

def getTmvaRoc(sample):
    tfile = ROOT.TFile("./tmva/TMVA_"+sample+".root")
    htmp = tfile.Get("dataset/Method_%s/%s/MVA_%s_rejBvsS"%(method,method,method))
    return copy.deepcopy(htmp)
    

datadir = '/xrootd/store/user/tt8888tt/muon/9_3_0_pre4/tmva/'
treename = "MuonAnalyser/reco"
rangecut = 'muon.Pt()>5&&abs(muon.Eta())<2.4'
method = "BDT"

samples= {"tenmu":[0,-0.06], "ttbar":[-0.02,-0.12], "zmm":[0,-0.06]}
col = [[ROOT.kRed, ROOT.kOrange+7, ROOT.kPink+9], [ROOT.kBlue,ROOT.kViolet+7,ROOT.kAzure+7], [ROOT.kGreen,ROOT.kTeal+9,ROOT.kSpring+9]]

grs = []
hists = []
marks = []
marks_id = []
auc = []
for i, sample in enumerate(samples):
    for j, applied in enumerate(samples):
        samplename = sample+"_"+applied+"applied"
        if sample==applied: samplename=sample
        #else: continue 
        filename = datadir+samplename+".root"
        print samplename

        h_sig = makeTH1(filename, treename, "muon_tmva_bdt", [100,-0.25,0.25], "muon_tmva_bdt", rangecut+"&&muon_signal")
        h_bkg = makeTH1(filename, treename, "muon_tmva_bdt", [100,-0.25,0.25], "muon_tmva_bdt", rangecut+"&&!muon_signal")

        gr = getRoc(h_sig, h_bkg) 
        gr.SetLineColor(col[i][j])
        #gr.SetTitle(sample.capitalize())
        gr.SetTitle("TMVA trained with "+applied.capitalize()+" on "+sample.capitalize())
        grs.append(gr)

        ##id point
        m1 = getMarker(filename, h_sig.Integral(0,101), h_bkg.Integral(0,101), "muon_tmva_bdt>%f"%samples[sample][0], 22)
        m2 = getMarker(filename, h_sig.Integral(0,101), h_bkg.Integral(0,101), "muon_tmva_bdt>%f"%samples[sample][1], 23)
        m1.SetMarkerSize(3)
        m2.SetMarkerSize(3)
        m1.SetMarkerColor(col[i][j])
        m2.SetMarkerColor(col[i][j])
        marks.append([m1, m2])

        m_tight = getMarker(filename, h_sig.Integral(0,101), h_bkg.Integral(0,101), "muon_isTightCustom", 26)
        m_loose = getMarker(filename, h_sig.Integral(0,101), h_bkg.Integral(0,101), "muon_isLoose", 32)
        m_tight.SetMarkerSize(3)
        m_loose.SetMarkerSize(3)
        m_tight.SetMarkerColor(col[i][j])
        m_loose.SetMarkerColor(col[i][j])
        marks_id.append([m_tight, m_loose])

        bdt = ROOT.vector('float')()
        mvaid = ROOT.vector('bool')()
        tfile = ROOT.TFile(filename)
        for t in tfile.Get("MuonAnalyser/reco"):
            if t.muon.Pt() < 5: continue
            if abs(t.muon.Eta()) > 2.4: continue
            bdt.push_back(t.muon_tmva_bdt)
            mvaid.push_back(t.muon_signal)
        auc = ROOT.TMVA.ROCCurve(bdt, mvaid).GetROCIntegral()
        gr.SetTitle(gr.GetTitle()+" (%.4f)"%auc)
        print gr.GetTitle()


h_init = ROOT.TH1F("","",1000,0.9,1.001)
h_init.SetMaximum(1.001)
h_init.SetMinimum(0.9)
h_init.GetXaxis().SetTitle("Signal Efficiency")
h_init.GetYaxis().SetTitle("Background Rejection")
h_init.GetYaxis().SetTitleOffset(1)

canv = makeCanvas("tmva", False)
canv.SetGrid()
setMargins(canv, False)
h_init.Draw()

#leg = ROOT.TLegend(0.15,0.18,0.4,0.4)
leg = ROOT.TLegend(0.15,0.17,0.65,0.55)

for i, gr in enumerate(grs):
    gr.SetLineWidth(3)
    gr.Draw("Lsame")
    leg.AddEntry(gr,gr.GetTitle(),"l")

    hists[i].Draw("same")
    leg.AddEntry(hists[i],hists[i].GetTitle(),"l")

    marks[i][0].Draw("same")
    marks[i][1].Draw("same")
    print marks[i][0].GetX(), marks[i][0].GetY()
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
canv.SaveAs("tmva_roc.png")

