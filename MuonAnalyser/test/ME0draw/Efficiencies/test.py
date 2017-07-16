import ROOT, copy
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)

h = makeTH1("../../TenMuPU0/pu0.root", "MuonAnalyser/gen", "PU200", [100,0,100], "genMuon_deltaXME0", "")

x_name = "genMuon_deltaXME0"
y_name = "events"
h.GetXaxis().SetTitle(x_name)
h.GetYaxis().SetTitle(y_name)

canv = makeCanvas("genMuon_deltaXME0", False)
setMargins(canv, False)
h.Draw()
canv.Modified()
canv.Update()

canv.SaveAs("genMuon_deltaXME0.png")

