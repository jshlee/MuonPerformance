import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(1.2)
    h.SetLineColor(color)
    h.SetLineWidth(2)

def getEff(filename,treename,title,binning,plotvar,dencut,numcut):
    h1 = makeTH1(filename,treename,title,binning,plotvar,dencut)
    h2 = makeTH1(filename,treename,title,binning,plotvar,numcut)
    hEff = ROOT.TEfficiency(h2,h1)
    hEff = delLargeErrorPoints(hEff,h1,h2,binning)
    hEff.SetTitle(title)
    return copy.deepcopy(hEff)

def delLargeErrorPoints(hEff, h1, h2, binning):
    maxY = hEff.GetEfficiency(1)
    for j in range(binning[0]):
        if maxY < hEff.GetEfficiency(j+1): maxY = hEff.GetEfficiency(j+1)
    for i in range(1,binning[0]+1):
        if hEff.GetEfficiencyErrorLow(i) > maxY/4.5:
            h1.SetBinContent(i,0)
            h2.SetBinContent(i,0)
    hEff = ROOT.TEfficiency(h2,h1)
    return copy.deepcopy(hEff)

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(62)
    tex2.SetTextSize(0.03)
    tex2.DrawLatex(0.25, 0.88, samplename)

def divByNevents(filedir, hist):
    tfile = ROOT.TFile(filedir)
    nevents = tfile.Get("MuonAnalyser/nevents").Integral()
    print filedir
    hist.Scale(1/nevents)

def draw(h_init, y_name, hlists, name, text):
    #Plot style
    setMarkerStyle(hlists[0], 4, 20) #blue, circle
    setMarkerStyle(hlists[1], 2, 34) #black, cross
    setMarkerStyle(hlists[2], 3, 21) #red, square
    #setMarkerStyle(hlists[3], 6, 24) #red, square
    #setMarkerStyle(hlists[3], 3, 31) #green, patrol

    #Set canvas
    #canv = makeCanvas(plotvar+name, False)
    #setMargins(canv, False)
    canv = ROOT.TCanvas()
    canv.SetGrid()
    h_init.GetYaxis().SetTitle(y_name)
    #h_init.GetYaxis().SetMinimum(0)
    h_init.Draw()
    drawSampleName(text)

    #Legend and drawing
    leg = ROOT.TLegend(0.25,0.72,0.45,0.86)
    for h in hlists:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"p")
    hlists[0].Draw("e1same")
    leg.SetTextFont(62)
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw()

    #CMS_lumi setting
    iPos = 0
    iPeriod = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_sqrtS = "14 TeV"
    CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

    canv.Modified()
    canv.Update()
    canv.SaveAs("%s_%s.png"%(plotvar,name))


#datadir = './'
datadir = '../'
#filenames = ["zmm.root", "zmm140.root", "zmm200.root"]
#filenames = ["GEMMuon.root"]
filenames = ["tenmu.root"]
#filenames = ["zmm.root", "zmm140.root", "zmm200.root", "run2.root"]
#filenames = ["relval_"+x for x in filenames]

#binning_l = [[20,-10,10], [40, -20, 20], [16, -8, 8], [16, -8, 8], [32,-0.6, 1]]
#binning_l = [[16,-4,4], [40, -10, 10], [12, -3, 3], [20, -5, 5], [32,-0.6, 1]]
#binning_l = [[20,5,105],[8,1.6,2.4],[30,-3,3]]
binning_l = [[20,5,105],[8,1.6,2.4]]
#binning_l = [[8, 1.6, 2.4], [8, 1.6, 2.4], [8, 1.6, 2.4], [8, 1.6, 2.4], [8, 1.6, 2.4]]
#rangecut = "muon.Pt()>5&&abs(muon.Eta())<2.4"
#rangecut_GE11 = "muon.Pt()>5 && abs(muon.Eta())<2.2 && abs(muon.Eta())>1.6"
#rangecut_GE21 = "muon.Pt()>5 && abs(muon.Eta())<2.4 && abs(muon.Eta())>1.6"

rangecut_GE11Loose = "muon_isGE11Muon && abs(muon_GE11deltaX) < 8.5185 && abs(muon_GE11deltaDXDZ) < 0.6269 && abs(muon_GE11dPhi) < 0.05115 && abs(muon_GE11dEta) < 0.3845  && abs(muon_GE11pullPhi) < 2.83683333 && abs(muon_GE11pullX) < 4.2275 && abs(muon_GE11dPhi*muon.Pt()*muon_Charge) < 0.4377125"
rangecut_GE11Tight = "muon_isGE11Muon && abs(muon_GE11deltaX) < 2.6285 && abs(muon_GE11deltaY) < 18.79875 && abs(muon_GE11deltaDXDZ) < 0.2873 && abs(muon_GE11dPhi) < 0.01605 && abs(muon_GE11dEta) < 0.1089  && abs(muon_GE11pullPhi) < 1.8285 && abs(muon_GE11pullX) < 2.1565 && abs(muon_GE11pullY) < 9.8805 && abs(muon_GE11dPhi*muon.P()) < 0.4944375 && abs(muon_GE11dPhi*muon.Pt()*muon_Charge) < 0.1445375"
rangecut_GE21Loose = "muon_isGE21Muon && abs(muon_GE21deltaDXDZ) < 0.9903 && abs(muon_GE21dPhi) < 0.08545 && abs(muon_GE21dEta) < 0.4721  && abs(muon_GE21pullPhi) < 2.96802631579 && abs(muon_GE21pullX) < 4.8915 && abs(muon_GE21dPhi*muon.Pt()*muon_Charge) < 0.488495"
rangecut_GE21Tight = "muon_isGE21Muon && abs(muon_GE21deltaX) < 5.2645 && abs(muon_GE21deltaY) < 20.02125 && abs(muon_GE21deltaDXDZ) < 0.3923 && abs(muon_GE21dPhi) < 0.02995 && abs(muon_GE21dEta) < 0.1027  && abs(muon_GE21pullPhi) < 2.24802631579 && abs(muon_GE21pullX) < 2.5155 && abs(muon_GE21pullY) < 7.2555 && abs(muon_GE21dPhi*muon.P()) < 0.898375 && abs(muon_GE21dPhi*muon.Pt()*muon_Charge) < 0.207375"
#rangecut_GE11Tight = "abs(muon_GE11deltaX) < 8 && abs(muon_GE11deltaY) < 13 && abs(muon_GE11pullX) < 8 && abs(muon_GE11pullY) < 9 && muon_GE11deltaDYDZ > 0.2  && abs(muon_GE11pullPhi) < 1.6"
#rangecut_GE21Loose = "abs(muon_GE21deltaX) < 10 && abs(muon_GE21deltaY) < 15 && abs(muon_GE21pullX) < 10 && abs(muon_GE21pullY) < 12 && muon_GE21deltaDYDZ > 0.33 && abs(muon_GE21pullPhi) < 2.0"
#rangecut_GE21Tight = "abs(muon_GE21deltaX) < 8 && abs(muon_GE21deltaY) < 13 && abs(muon_GE21pullX) < 8 && abs(muon_GE21pullY) < 9 && muon_GE21deltaDYDZ > 0.43 && abs(muon_GE21pullPhi) < 1.6"
    
#Set extra text
#samplename = "Z/#gamma^{*}#rightarrow#font[12]{#mu#mu}"
samplename = "10#mu"
GE11_sig = "abs(muon.Eta()) > 1.6 && abs(muon.Eta())<2.2 && muon.Pt() > 5"
#GE11_sig = "muon_signal && abs(muon.Eta()) > 1.6 && abs(muon.Eta())<2.2 && muon.Pt() > 5"
GE11_bkg = "!muon_signal && abs(muon.Eta()) > 1.6 && abs(muon.Eta())<2.2 && muon.Pt() > 5 "
GE21_sig = "muon_signal && abs(muon.Eta()) > 1.6 && abs(muon.Eta())<2.4 && muon.Pt() > 5"
#GE21_sig = "muon_isGE21Muon && muon_signal"
GE21_bkg = "!muon_signal && abs(muon.Eta()) > 1.6 && abs(muon.Eta())<2.4 && muon.Pt() > 5 "
text = samplename

#for i, plotvar in enumerate(["deltaX", "deltaY", "pullX", "pullY", "dPhi"]):
#for i, plotvar in enumerate(["muon.Pt()", "abs(muon.Eta())", "muon.Phi()",]):
for i, plotvar in enumerate(["muon.Pt()", "abs(muon.Eta())"]):
    #Efficiency
    hl_eff_GE11 = []
    hl_eff_GE21 = []
    #hl_eff.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE11Muon_Loose", binning_l[i], "muon_GE11"+plotvar, GE11_sig, "%s && %s"%(rangecut_GE11Loose, GE11_sig)))
    #hl_eff.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE11Muon_Tight", binning_l[i], "muon_GE11"+plotvar, GE11_sig, "%s && %s"%(rangecut_GE11Tight, GE11_sig)))
    #hl_eff.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE21Muon_Loose", binning_l[i], "muon_GE21"+plotvar, GE21_sig, "%s && %s"%(rangecut_GE21Loose, GE21_sig)))
    #hl_eff.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE21Muon_Tight", binning_l[i], "muon_GE21"+plotvar, GE21_sig, "%s && %s"%(rangecut_GE21Tight, GE21_sig)))
    hl_eff_GE11.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE11Muon", binning_l[i], plotvar, GE11_sig, GE11_sig+"&&muon_isGE11Muon"))
    hl_eff_GE11.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE11Muon_Loose", binning_l[i], plotvar, GE11_sig, "%s && %s"%(rangecut_GE11Loose, GE11_sig)))
    hl_eff_GE11.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE11Muon_Tight", binning_l[i], plotvar, GE11_sig, "%s && %s"%(rangecut_GE11Tight, GE11_sig)))
    
    hl_eff_GE21.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE21Muon", binning_l[i], plotvar, GE21_sig, GE21_sig+"&&muon_isGE21Muon"))
    hl_eff_GE21.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE21Muon_Loose", binning_l[i], plotvar, GE21_sig, "%s && %s"%(rangecut_GE21Loose, GE21_sig)))
    hl_eff_GE21.append(getEff(datadir+filenames[0], "MuonAnalyser/gen", "GE21Muon_Tight", binning_l[i], plotvar, GE21_sig, "%s && %s"%(rangecut_GE21Tight, GE21_sig)))
    
    #Backgorund rate
    hl_bkg_GE11 = []
    hl_bkg_GE21 = []
    #hl_bkg.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE11Muon_Loose", binning_l[i], "muon_GE11"+plotvar, "%s && %s"%(rangecut_GE11Loose, GE11_bkg)))
    #hl_bkg.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE11Muon_Tight", binning_l[i], "muon_GE11"+plotvar, "%s && %s"%(rangecut_GE11Tight, GE11_bkg)))
    #hl_bkg.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE21Muon_Loose", binning_l[i], "muon_GE21"+plotvar, "%s && %s"%(rangecut_GE21Loose, GE21_bkg)))
    #hl_bkg.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE21Muon_Tight", binning_l[i], "muon_GE21"+plotvar, "%s && %s"%(rangecut_GE21Tight, GE21_bkg)))
    #hl_eff_GE11.append(getEff(datadir+filenames[0], "MuonAnalyser/reco", "GE11Muon", binning_l[i], plotvar, nevents, GE11_sig+"&&muon_isGE21Muon"))
    tfile = ROOT.TFile(datadir+filenames[0])
    nevents = tfile.Get("MuonAnalyser/nevents").Integral()

    hl_bkg_GE11.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE11Muon_bkg", binning_l[i], plotvar, GE11_bkg))
    hl_bkg_GE11.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE11Muon_bkg_Loose", binning_l[i], plotvar, "%s && %s"%(rangecut_GE11Loose, GE11_bkg)))
    hl_bkg_GE11.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE11Muon_bkg_Tight", binning_l[i], plotvar, "%s && %s"%(rangecut_GE11Tight, GE11_bkg)))
    
    
    for a in hl_bkg_GE11:
        a.Scale(1/nevents)
    
    hl_bkg_GE21.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE21Muon_bkg", binning_l[i], plotvar, GE21_bkg))
    hl_bkg_GE21.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE21Muon_bkg_Loose", binning_l[i], plotvar, "%s && %s"%(rangecut_GE21Loose, GE21_bkg)))
    hl_bkg_GE21.append(makeTH1(datadir+filenames[0], "MuonAnalyser/reco", "GE21Muon_bkg_Tight", binning_l[i], plotvar, "%s && %s"%(rangecut_GE21Tight, GE21_bkg)))
    for b in hl_bkg_GE21:
        b.Scale(1/nevents)
    '''for j in range(len(hl_bkg)):
        divByNevents(datadir+filenames[0], hl_bkg[j])
        avg = 0
        for k in range(1,binning_l[i][0]+1):
            avg += hl_bkg[j].GetBinContent(k)
        avg = avg/float(binning_l[i][0])
        print plotvar, avg
    '''
    #if "density" in plotvar: plotvar = plotvar.split('/')[0]

    #Set X axis name
    x_name = "Muon "
    if "deltaX"  in plotvar: x_name = "deltaX"
    elif "deltaY" in plotvar: x_name = "deltaY"
    elif "pullX" in plotvar: x_name = "pullX"
    elif "pullY" in plotvar : x_name = "pullY"
    elif "dPhi" in plotvar : x_name = "dPhi"
    elif "Pt" in plotvar : x_name = "p_{T}"
    elif "Eta" in plotvar : x_name = "|#eta|"
    
    #Set init histo
    h_init_GE11 = ROOT.TH1F("","",binning_l[i][0],binning_l[i][1],binning_l[i][2])
    h_init_GE21 = ROOT.TH1F("","",binning_l[i][0],binning_l[i][1],binning_l[i][2])

    h_init_GE11.GetXaxis().SetTitle(x_name)
    h_init_GE11.GetYaxis().SetTitleOffset(1)
    h_init_GE11.GetXaxis().SetTitleOffset(1.1)
    h_init_GE11.GetYaxis().SetTitleOffset(1.5)
    h_init_GE11.GetXaxis().SetTitleSize(0.05)
    h_init_GE11.GetXaxis().SetLabelSize(0.038)
    h_init_GE11.GetYaxis().SetTitleSize(0.05)
    h_init_GE11.GetYaxis().SetLabelSize(0.038)

    h_init_GE21.GetXaxis().SetTitle(x_name)
    h_init_GE21.GetYaxis().SetTitleOffset(1)
    h_init_GE21.GetXaxis().SetTitleOffset(1.1)
    h_init_GE21.GetYaxis().SetTitleOffset(1.5)
    h_init_GE21.GetXaxis().SetTitleSize(0.05)
    h_init_GE21.GetXaxis().SetLabelSize(0.038)
    h_init_GE21.GetYaxis().SetTitleSize(0.05)
    h_init_GE21.GetYaxis().SetLabelSize(0.038)

    h_init_GE11_bkg = h_init_GE11.Clone()
    h_init_GE21_bkg = h_init_GE11.Clone()
   
     #Set Y axis name
    y_name = " Muon "
    #if "Custom" in idname: y_name = "Tight Muon "

    h_init_GE11.SetMaximum(1.5)
    h_init_GE11.SetMinimum(0.0)
    h_init_GE21.SetMaximum(1.5)
    h_init_GE21.SetMinimum(0.0)
    draw(h_init_GE11, "Efficiency", hl_eff_GE11, "eff_GE11", text)
    draw(h_init_GE21, "Efficiency", hl_eff_GE21, "eff_GE21", text)

    h_init_GE11_bkg.GetXaxis().SetTitle(x_name)
    h_init_GE11_bkg.GetYaxis().SetTitleOffset(1)
    h_init_GE11_bkg.GetXaxis().SetTitleOffset(1.1)
    h_init_GE11_bkg.GetYaxis().SetTitleOffset(1.5)
    h_init_GE11_bkg.GetXaxis().SetTitleSize(0.05)
    h_init_GE11_bkg.GetXaxis().SetLabelSize(0.038)
    h_init_GE11_bkg.GetYaxis().SetTitleSize(0.05)
    h_init_GE11_bkg.GetYaxis().SetLabelSize(0.038)
    h_init_GE11_bkg.SetMaximum(0.5)
    h_init_GE11_bkg.SetMinimum(0.0)

    h_init_GE21_bkg.GetXaxis().SetTitle(x_name)
    h_init_GE21_bkg.GetYaxis().SetTitleOffset(1)
    h_init_GE21_bkg.GetXaxis().SetTitleOffset(1.1)
    h_init_GE21_bkg.GetYaxis().SetTitleOffset(1.5)
    h_init_GE21_bkg.GetXaxis().SetTitleSize(0.05)
    h_init_GE21_bkg.GetXaxis().SetLabelSize(0.038)
    h_init_GE21_bkg.GetYaxis().SetTitleSize(0.05)
    h_init_GE21_bkg.GetYaxis().SetLabelSize(0.038)
    h_init_GE21_bkg.SetMaximum(0.5)
    h_init_GE21_bkg.SetMinimum(0.0)
    
    #draw(h_init_GE11_bkg, "Bkg", hl_bkg_GE11, "bkgrate_GE11", text)
    #draw(h_init_GE21_bkg, "Bkg", hl_bkg_GE11, "bkgrate_GE21", text)

    if "()" in plotvar: 
        h_init_GE11_bkg.SetMaximum(max(h.GetMaximum() for h in hl_bkg_GE11)*1.5)
        h_init_GE21_bkg.SetMaximum(max(h.GetMaximum() for h in hl_bkg_GE21)*1.5)
         
        draw(h_init_GE11_bkg, "Backgrounds", hl_bkg_GE11, "bkg_GE11", text)
        draw(h_init_GE21_bkg, "Backgorunds", hl_bkg_GE21, "bkg_GE21", text)
    


