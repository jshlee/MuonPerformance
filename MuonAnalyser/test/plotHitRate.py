import ROOT
ROOT.gROOT.GetGlobal( "gStyle", 1 )
ROOT.gStyle.SetOptStat(0);

tfile = ROOT.TFile('histo.root')
tree = tfile.Get('HitRateAnalysis')

c1 = ROOT.TCanvas()
legend = ROOT.TLegend(0.7,0.7,0.9,0.9);
ch1 = [0,1,2,3,8,9]
ch2 = [4,5,6,7,10,11]
ch = ch1

for i in ch:    
    t1 = tree.Get("nPadPerGEB M%1i"%(i+1))
    c = i+1
    line = 1
    if c > 4:
        c = c-4
        line = 2
        
    t1.SetLineColor(c)
    #t1.SetLineStyle(line)
    t1.SetAxisRange(0,50)
    t1.SetLineWidth(2*line)
    t1.DrawNormalized("same")
    legend.AddEntry(t1,"M%1i"%(i+1),"f");

legend.Draw()
c1.SetLogy()
c1.Print("nPadPerGEB.png")


c2 = ROOT.TCanvas()
for i in ch:
    t1 = tree.Get("nPadPerGEB M%1i"%(i+1))
    #t2 = t1.Clone()
    t2 = ROOT.TH1D( "cumulative nPadPerGEB", "cumulative nPadPerGEB", 50, 0, 50)
    #t2.SetAxisRange(0,50)
    
    for j in range(51):
        x = t1.Integral(j, 100)
        #print j,x
        t2.SetBinContent(j,x)
        
    c = i+1
    line = 1
    if c > 4:
        c = c-4
        line = 2
        
    t2.SetLineColor(c)
    #t2.SetLineStyle(line)
    t2.SetLineWidth(2*line)
   
    t2.DrawNormalized("same")

legend.Draw()
c2.SetLogy()
c2.Print("cumulative_nPadPerGEB.png")
