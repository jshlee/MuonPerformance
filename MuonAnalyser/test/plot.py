import ROOT

f = ROOT.TFile("sliceTest2018.root")
#f = ROOT.TFile("histo2018.root")
t = f.SliceTestAnalysis.Get("Muon")

#t.Draw("track_roll>>eta(8,1,9)", "onDet&&abs(track_roll)<10")
#t.Draw("track_roll>>eta2(8,1,9)", "!onDet&&abs(track_roll)<10")
#eta = ROOT.gROOT.FindObject("eta")
#eta2 = ROOT.gROOT.FindObject("eta2")
#
#t.Draw("track_chamber>>ch(36,0,36)", "onDet&&abs(track_chamber)>0")
#t.Draw("track_chamber>>ch2(36,0,36)", "!onDet&&abs(track_chamber)>0")
#ch = ROOT.gROOT.FindObject("ch")
#ch2 = ROOT.gROOT.FindObject("ch2")
#
cnv = ROOT.TCanvas()
#eff = ROOT.TEfficiency(eta,eta2)
#eff.Draw()
#eff.SetTitle("efficiency per etaPartition; etaPartition; Efficiency")
#cnv.Print("eff_eta.png")
#
#eff2 = ROOT.TEfficiency(ch,ch2)
#eff2.Draw()
#eff2.SetTitle("efficiency per chamber; chamber; Efficiency")
#cnv.Print("eff_ch.png")

t2 = f.SliceTestAnalysis.Get("Hit")
"""
t2.Draw("nStrips>>s(10,0,10)","etaPartition==1")
s = ROOT.gROOT.FindObject("s")
s.SetMaximum(s.GetMaximum()*3.5)
s.Draw()
s.SetStats(0)
s.SetTitle("nStrips; nStrips; number of hit")
leg = ROOT.TLegend(0.6,0.45,0.88,0.88)
for i in range(1,9):
    t2.Draw("nStrips>>nstrip%d"%i,"etaPartition==%d"%i,"same")
    strip = ROOT.gROOT.FindObject("nstrip%d"%i)
    strip.SetLineColor(i+1)
    leg.AddEntry(strip, "etaPartition %d"%i, "l")
leg.Draw()
cnv.Print("clusterSize.png")

h = f.SliceTestAnalysis.Get("res_x")
h.Draw()
cnv.Print("res_x.png")
h = f.SliceTestAnalysis.Get("res_y")
h.Draw()
cnv.Print("res_y.png")

h = f.SliceTestAnalysis.Get("pull_x")
h.Draw()
cnv.Print("pull_x.png")
h = f.SliceTestAnalysis.Get("pull_y")
h.Draw()
cnv.Print("pull_y.png")

t2.Draw("y:x>>xy", "muonQuality>=2")
xy = ROOT.gROOT.FindObject("xy")
xy.SetTitle("xy occupancy; global x; global y")
xy.Draw()
cnv.Print("xy.png")
t2.Draw("sqrt(x**2+y**2):z>>rz")
rz = ROOT.gROOT.FindObject("rz")
rz.SetTitle("rz occupancy; global r; global z")
rz.Draw()
cnv.Print("rz.png")
"""

t.Draw("resx>>rtot(10,-200,200)", "nvalidhits==3")
t.Draw("resx[0]>>r0(10,-200,200)", "nvalidhits==3&&resx[0]>resx[1]&&resx[0]>resx[2]")
t.Draw("resx[1]>>r1(10,-200,200)", "nvalidhits==3&&resx[1]>resx[0]&&resx[1]>resx[2]")
t.Draw("resx[2]>>r2(10,-200,200)", "nvalidhits==3&&resx[2]>resx[0]&&resx[2]>resx[1]")
rtot = ROOT.gROOT.FindObject("rtot")
r0 = ROOT.gROOT.FindObject("r0")
r1 = ROOT.gROOT.FindObject("r1")
r2 = ROOT.gROOT.FindObject("r2")
large = r0.Clone()
large.Add(r1)
large.Add(r2)
large.SetLineColor(2)
rtot.Add(large,-1)
rtot.Draw()
large.Draw("same")
cnv.Print("hit3.png")





