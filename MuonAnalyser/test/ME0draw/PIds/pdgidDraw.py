#!/usr/bin/env python

import ROOT, copy, math 
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
from MuonPerformance.MuonAnalyser.histoHelper import *

pileup = "pu200"
filename = "%s.root" % pileup
pid = "abs(recoMuon_pdgId)"
treename = "MuonAnalyser/reco"
plotvar = ["recoMuon.Eta()","recoMuon.Pt()"]
binnings = [[10,2.0,2.8],[20,0,200]]
isME0 = ["&& recoMuon_isME0Muon &&","&&"]

f = open("PIds/pId_log.txt", 'w')

for i, plotvar in enumerate(plotvar):
    for me0opt in isME0:
        ### Canvas ###
        c = makeCanvas(plotvar+me0opt, False)
        #c.SetLogy()
        l = ROOT.TLegend(0.85,0.7, 0.97, 0.95)
        l.SetTextSize(0.04)
        ### Histo ###
        hs = ROOT.THStack("hs","")  #stacked histo
        regionCut = "abs(recoMuon.Eta())>2.0 && abs(recoMuon.Eta())<2.8"

        h_mu = makeTH1(filename, treename, "Muon", binnings[i], plotvar,
                       "recoMuon.Pt()>5 %s %s && !recoMuon_signal && %s==13" % (me0opt, regionCut, pid))
        h_mu.SetFillColor(4)
        l.AddEntry(h_mu,"#mu^{#pm}")
        hs.Add(h_mu)

        h_pi = makeTH1(filename, treename, "Pion", binnings[i], plotvar,
                       "recoMuon.Pt()>5 %s %s && !recoMuon_signal && %s==211" % (me0opt, regionCut, pid))
        h_pi.SetFillColor(7)
        l.AddEntry(h_pi,"#pi^{#pm}")
        hs.Add(h_pi)

        h_K = makeTH1(filename, treename, "Kaon", binnings[i], plotvar,
                       "recoMuon.Pt()>5 %s %s && !recoMuon_signal && %s==321" % (me0opt, regionCut, pid))
        h_K.SetFillColor(5)
        l.AddEntry(h_K,"K^{#pm}")
        hs.Add(h_K)

        h_p = makeTH1(filename, treename, "Proton", binnings[i], plotvar,
                       "recoMuon.Pt()>5 %s %s && !recoMuon_signal && %s==2212" % (me0opt, regionCut, pid))
        h_p.SetFillColor(3)
        l.AddEntry(h_p,"p")
        hs.Add(h_p)
        
        h_other = makeTH1(filename, treename, "Others", binnings[i], plotvar,
                       "recoMuon.Pt()>5 %s %s && !recoMuon_signal && %s!=13 && %s!=211 && %s!=321 && %s!=2212 && %s!=0" % (me0opt, regionCut, pid, pid, pid, pid, pid))
        h_other.SetFillColor(2)
        l.AddEntry(h_other,"others")
        hs.Add(h_other)

        hs.Draw("hist")
        if "Eta()" in plotvar:
            hs.GetXaxis().SetTitle("|#eta|")
        if "Pt()" in plotvar:
            hs.GetXaxis().SetTitle("p_{T}")
        hs.GetYaxis().SetTitle("%s Entries" % pileup)
        if "Eta()" in plotvar:
            hs.SetMinimum(0)
            hs.SetMaximum(5000)
        if "Pt()" in plotvar:
            hs.SetMinimum(1)
            c.SetLogy()
        """
        nbins = h_mu.GetNbinsX()
        h_dummy = ROOT.TH1F("","",nbins,h_mu.GetBinLowEdge(1),h_mu.GetBinLowEdge(nbins+1))
        if "Pt()" in plotvar:
            h_dummy.GetXaxis().SetTitle("p_{T}")
        if "Eta()" in plotvar:
            h_dummy.GetXaxis().SetTitle("|#eta|")
        h_dummy.GetYaxis().SetTitle("PU200 Entries")
        h_dummy.SetMinimum(1)
        h_dummy.Draw("histsame")
        """
        l.Draw()

        ### CMS_lumi setting ###
        iPos = 0
        iPeriod = 0
        if (iPos == 0):
            CMS_lumi.relPosX = 0.12
        CMS_lumi.extraText = "Work in progress"
        CMS_lumi.lumi_sqrtS = "14TeV"
        CMS_lumi.CMS_lumi(c, iPeriod, iPos)

        if "&& recoMuon_isME0Muon &&" == me0opt:
            me0option = "isME0"
        if "&&" == me0opt:
            me0option = "isAll"

        c.Modified()
        c.Update()
        c.SaveAs("PIds/%s_pId_%s_%s.png" % (pileup, me0option, plotvar))

        ### Make table ###
        status = "%s_pId_%s_%s\n" % (pileup, me0option, plotvar)
        f.write(status)
        f.write("mu     = %s\n" % h_mu.GetEntries())
        f.write("pion   = %s\n" % h_pi.GetEntries())
        f.write("kaon   = %s\n" % h_K.GetEntries())
        f.write("proton = %s\n" % h_p.GetEntries())
        f.write("others = %s\n\n" % h_other.GetEntries())

f.close()

