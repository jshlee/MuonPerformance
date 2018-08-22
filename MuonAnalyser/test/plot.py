import ROOT, array
f = ROOT.TFile("sliceTest2018.root")
fv2 = ROOT.TFile("sliceTest2018v2.root")

tevent = f.SliceTestAnalysis.Get("Events")
tmuon=f.SliceTestAnalysis.Get("Muon")
thit=f.SliceTestAnalysis.Get("Hit")
cnv = ROOT.TCanvas()

"""
tmuon.Draw("run>>h")
h_run = ROOT.gROOT.FindObject("h")
runs=[]
for i in range(1,h_run.GetNbinsX()+1):
  run = int(h_run.GetBinCenter(i))
  if tmuon.Draw("nhits","run==%d"%run) != 0: runs.append(int(h_run.GetBinCenter(i)))
print runs
htmp = ROOT.TH1D("", "",len(runs),0,len(runs))

tex = ROOT.TLatex()
tex.SetNDC()
"""

idcut = "pt>5&&quality>=2"
for ch in [27,28,29,30]:
    for lay in [1,2]:
	print ch, lay
        tmuon.Draw("in_roll:in_strip/128>>in_strip(3,0,3,8,1,9)", "%s&&nvalidhits>0&&in_layer==%d&&in_chamber==%d"%(idcut,lay,ch))
        tmuon.Draw("in_roll:in_strip/128>>hit_strip(3,0,3,8,1,9)", "%s&&nvalidhits>0&&in_nearGemRoll==in_roll&&in_nearGemFirstStrip!=-9&&(in_nearGemFirstStrip/128)==(in_strip/128)&&in_layer==%d&&in_chamber==%d"%(idcut,lay,ch))
        h = ROOT.gROOT.FindObject("in_strip")
        h2 = ROOT.gROOT.FindObject("hit_strip")
        eff = ROOT.TEfficiency(h2,h)
        eff.SetTitle("ch %d lay %d; VFAT; Roll"%(ch, lay))
        eff.Draw()
        cnv.Update()
        eff.GetPaintedHistogram().SetMaximum(0.2)
        eff.GetPaintedHistogram().SetMinimum(0)
        eff.SetMarkerSize(1.3)
        ROOT.gStyle.SetPaintTextFormat(".4f")
        eff.Draw("coltextz")
        cnv.Print("eff_ch%d_lay%d.png"%(ch,lay))

	"""
        for i,run in enumerate(runs):
	    tot = tmuon.Draw("nhits","run==%d&&in_chamber==%d&&in_layer==%d"%(run,ch,lay))
	    passed = tmuon.Draw("nhits","nhits>0&&run==%d&&in_chamber==%d&&in_layer==%d"%(run,ch,lay))
	    h.GetXaxis().SetBinLabel(i+1,str(run))
	    if tot != 0: h.SetBinContent(i+1,passed/float(tot))
        h.Draw()
        h.SetStats(0)
        #h.SetMaximum(0.22)
        h.SetMaximum(0.3)
        h.SetMinimum(0)
        h.SetTitle("Efficiency per Run (ch %d lay%d);Run;Efficiency"%(ch,lay))
        tex.DrawLatex(0.15, 0.75, "SliceTest 2018C")
        cnv.Print("runeff_ch%d_lay%d.png"%(ch,lay))
        h.Reset()
	"""
