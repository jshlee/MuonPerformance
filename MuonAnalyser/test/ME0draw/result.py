import ROOT, copy

def rebin(h, binning):
    xaxis = h.GetXaxis()
    binSize1 = h.GetBinWidth(1)
    binSize2 = abs(binning[2]-binning[1])/float(binning[0])

    h.Rebin(int(binSize2/binSize1))
    xaxis.SetRangeUser(binning[1],binning[2])

def histinit(h, binning):
    h.SetStats(0)
    rebin(h, binning)

#binninglst = [[30,0,60],[12,-2.4,2.4],[6,-3,3]]
binninglst = [[20,5,105],[24,-2.4,2.4],[30,-3,3]]
filenames=["out"]
#filenames=["nopu","pu35","pu200"]
for filename in filenames:
    tfile = ROOT.TFile(""+filename+".root")
    nevents = tfile.Get("MuonEff/nevents").GetEntries()
    for i, plotvar in enumerate(['pt','eta','phi']):
        gen = tfile.Get("MuonEff/gen_%s"%(plotvar))
        histinit(gen,binninglst[i])
        gen.Sumw2()
        eff=gen.Clone()
        lst=[]
        fakelst=[]
        for id in ['Tight','Medium','Loose']:
        #for id in ['T','M','L']:
            reco = tfile.Get("MuonEff/gen%s_%s"%(id[0],plotvar))
            histinit(reco,binninglst[i])
            reco.Sumw2()
            eff.Divide(reco, gen, 1, 1, "B")
            eff.SetTitle(id)
            lst.append(copy.deepcopy(eff))

            fake = tfile.Get("MuonEff/fake%s_%s"%(id[0],plotvar))
            fake.Scale(1/nevents)
            histinit(fake,binninglst[i])
            fake.SetTitle(id)
            fakelst.append(copy.deepcopy(fake))
        """
        reco = tfile.Get("MuonEff/genMatched %s %s"%('tightbystep',plotvar))
        histinit(reco,binninglst[i])
        reco.Sumw2()
        eff.Divide(reco, gen, 1, 1, "B")
        eff.SetTitle("tight(w/o dz(PV) cut)")
        fixed = copy.deepcopy(eff)

        print "MuonEff/reco %sfake %s"%('tightbystep',plotvar)
        fake = tfile.Get("MuonEff/reco %sfake %s"%('tightbystep',plotvar))
        fake.Scale(1/nevents)
        histinit(fake,binninglst[i])
        fake.SetTitle("tight(w/o dz(PV) cut)")
        fixed_fake = copy.deepcopy(fake)
        """

        eff.Reset()
        eff.SetMaximum(1.05)
        eff.SetMinimum(0.4)
        eff.SetLineColor(0)
        eff.SetTitle("%s_%s"%(filename,plotvar))

        cn = ROOT.TCanvas()
        eff.Draw()
        leg = ROOT.TLegend(0.58,0.14,0.88,0.4)
        for i, h in enumerate(lst):
            h.SetFillColor(i+2)
            h.SetLineColor(0)
            h.SetFillStyle(3001)
            h.SetMarkerStyle(0)
            h.Draw("e5same")
            leg.AddEntry(h,h.GetTitle(),"f")

        """
        fixed.SetLineColor(0)
        fixed.SetFillColor(2)
        fixed.SetMarkerStyle(0)
        fixed.Draw("e5same")
        leg.AddEntry(fixed,fixed.GetTitle(),"f")
        """

        leg.SetBorderSize(0)
        leg.Draw()
        cn.Print("temp/%s_%s.png"%(filename,plotvar))

        cn2 = ROOT.TCanvas()
        leg = ROOT.TLegend(0.58,0.58,0.88,0.88)
        #leg.SetNColumns(2)
        fake = copy.deepcopy(fakelst[0])
        fake.SetLineColor(0)
        fake.SetTitle("fake %s %s"%(filename,plotvar))
        fake.SetMaximum(max(h.GetMaximum() for h in fakelst)*2.2)
        fake.Draw()
        for i, h in enumerate(fakelst):
            h.SetLineColor(i+2)
            h.SetLineWidth(2)
            h.Draw("same")
            leg.AddEntry(h,h.GetTitle(),"f")
        
        """
        fixed_fake.SetLineColor(2)
        fixed_fake.SetLineWidth(3)
        fixed_fake.SetLineStyle(2)
        fixed_fake.Draw("same")
        leg.AddEntry(fixed_fake,fixed_fake.GetTitle(),"f")
        """

        leg.SetBorderSize(0)
        leg.Draw()
        cn2.Print("temp/fake_%s_%s.png"%(filename,plotvar))


