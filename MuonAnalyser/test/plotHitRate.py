import ROOT
ROOT.gROOT.GetGlobal( "gStyle", 1 )
ROOT.gStyle.SetOptStat(0);

tfile = ROOT.TFile('histo.root')
tree = tfile.Get('HitRateAnalysis')

c1 = ROOT.TCanvas()
legend = ROOT.TLegend(0.6,0.6,0.9,0.9);
ch1 = [8,9,0,1,2,3]
ch2 = [10,11,4,5,6,7]
ch = ch2
xrange2 = 10
histoName = "nPadPerGEB"
histoName = "nPadGEB"
histoName = "nPadRangeGEB"

for i in ch:    
    t1 = tree.Get(histoName+" M%1i"%(i+1))
    
    c = i+1
    line = 1
    if c > 4:
        c = c-4
    if i > 7:
        line = 3
        if i == 8 or i ==10:
            c = 6
        if i == 9 or i ==11:
            c = 7
        
    t1.GetYaxis().SetTitle("Fraction of Chambers");
    t1.SetLineColor(c)
    t1.SetLineStyle(line)
    t1.SetAxisRange(0,xrange2)
    t1.SetLineWidth(2*line)
    t1.SetMinimum(0.1)
    #t1.DrawNormalized("same")
    t1.Draw("same")
    legend.AddEntry(t1,"M%1i"%(i+1),"f")

legend.Draw()
c1.SetLogy()
c1.Print(histoName+".png")


c2 = ROOT.TCanvas()
for i in ch:
    t1 = tree.Get(histoName+" M%1i"%(i+1))
    #t2 = t1.Clone()
    t2 = ROOT.TH1D( "cumulative nPadPerGEB", "cumulative nPadPerGEB", xrange2+1, 0, xrange2+1)
    #t2.SetAxisRange(0,50)
    
    for j in range(51):
        x = t1.Integral(j, 100)
        #print j,x
        t2.SetBinContent(j,x)
        
    c = i+1
    line = 1
    if c > 4:
        c = c-4
    if i > 7:
        line = 3
        if i == 8 or i ==10:
            c = 6
        if i == 9 or i ==11:
            c = 7
        
    t2.SetLineColor(c)
    t2.SetLineStyle(line)
    t2.SetLineWidth(2*line)
    t2.GetYaxis().SetTitle("Fraction of Chambers")
    xaxistitle = t1.GetXaxis().GetTitle()
    t2.GetXaxis().SetTitle("cumulative "+xaxistitle)
    t2.SetMinimum(0.1);
    t2.DrawNormalized("same")

legend.Draw()
c2.SetLogy()
c2.Print(histoName+"_cumulative.png")
