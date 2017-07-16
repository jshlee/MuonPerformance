#!/usr/bin/env python
import ROOT, copy, os, sys, getopt
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

"""
opt.py -f "filename" -c "chamber"
"""

#get input
try:
    opts, args = getopt.getopt(sys.argv[1:],"f:c:",["filename","chamber"])
except getopt.GetoptError:
    print 'Usage : ./opt.py -f <filename> -c <chamber>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./opt.py -f <filename> -c <chamber>'
        sys.exit()
    elif opt in ("-f", "--filename"):
        filename = arg
    elif opt in ("-c", "--chamber"):
        chamber = arg       

#filename = "../../ZMM_PU200_pre5/ZMM_PU200_pre5.root"
gentree  = "MuonAnalyser/gen"
recotree = "MuonAnalyser/reco"

binning = [8, 1.4, 2.8]
plotvar = "abs(muon.Eta())"

if chamber == "GE11":
    defaultCut = "muon_isGEMMuon && muon.Pt()>5"
    CutVal = ["muon_GE11noRecHit", "fabs(muon_GE11deltaX)", "fabs(muon_GE11deltaY)", "fabs(muon_GE11deltaDXDZ)", "fabs(muon_GE11deltaDYDZ)", "fabs(muon_GE11pullX)", "fabs(muon_GE11pullY)", "fabs(muon_GE11dPhi)","fabs(muon_GE11dEta)",]
if chamber == "GE21":
    defaultCut = "muon_isGEMMuon && muon.Pt()>5"
    CutVal = ["muon_GE21noRecHit", "fabs(muon_GE21deltaX)", "fabs(muon_GE21deltaY)", "fabs(muon_GE21deltaDXDZ)", "fabs(muon_GE21deltaDYDZ)", "fabs(muon_GE21pullX)", "fabs(muon_GE21pullY)", "fabs(muon_GE21dPhi)","fabs(muon_GE21dEta)"]
if chamber == "ME0":
    defaultCut = "muon_isME0Muon"
    CutVal = ["muon_ME0noRecHit", "fabs(muon_ME0deltaX)", "fabs(muon_ME0deltaY)", "fabs(muon_ME0deltaDXDZ)", "fabs(muon_ME0deltaDYDZ)", "fabs(muon_ME0pullX)", "fabs(muon_ME0pullY)", "fabs(muon_ME0dPhi)", "fabs(muon_ME0dEta)"]
"""
if chamber == "ME0":
    defaultCut = "muon_isME0Muon"
    CutVal = ["fabs(muon_ME0dPhi)", "fabs(muon_ME0dEta)"]
"""
f = open("optLog_%s.txt"%chamber,"w")


for cutval in CutVal:
    if "noRecHit" in cutval:
        cuts = ["1", "2", "3", "4", "5", "6", "7"]
    if "deltaX" in cutval:
        cuts = ["2.4","2.3","2.2","1.2","1.1","1.0","0.9","0.8","0.7"]
    if "deltaY" in cutval:
        cuts = ["20","19","9","8","7","6","5","4","3"]
    if "deltaDXDZ" in cutval:
        cuts = ["0.1","0.09","0.08","0.07","0.06"]
    if "deltaDYDZ" in cutval:
        cuts = ["0.50","0.49","0.48","0.47","0.46","0.45","0.44","0.43","0.435","0.42","0.415"]
    if "pullX" in cutval:
        cuts = ["2.6","2.5","2.4","2.3","1.4","1.3","1.2","1.1","1.0"]
    if "pullY" in cutval:
        cuts = ["8","7","6","3.6","3.5","3.4"]
    if "dPhi" in cutval:
        cuts = ["0.015","0.014","0.013","0.012","0.011","0.010","0.009","0.008","0.007","0.006","0.005","0.004","0.003","0.002"]
    if "dEta" in cutval:
        cuts = ["0.05","0.048","0.046","0.044","0.042","0.040","0.038","0.036","0.035","0.034","0.032","0.030"]
    f.write( "--------------- %s Start ---------------\n" % cutval)

    for cut in cuts:
        h_nocut = makeTH1(filename, gentree, "nocut", binning, plotvar, "%s && %s!=100"%(defaultCut,cutval))
        nocutEntries = h_nocut.Integral()
        if "noRecHit" in cutval:
            h_pass = makeTH1(filename, gentree, cutval, binning, plotvar, "%s && %s > %s && %s!=100"%(defaultCut,cutval, cut, cutval))
        else:
            h_pass = makeTH1(filename, gentree, cutval, binning, plotvar, "%s && %s < %s && %s!=100"%(defaultCut,cutval, cut, cutval))
        cutEntries = h_pass.Integral()
        if "noRecHit" in cutval:
            f.write( "    %s Efficiency of %s > %s :  %s\n" % (chamber, cutval, cut, cutEntries/nocutEntries) )
        else:
            f.write( "    %s Efficiency of %s < %s :  %s\n" % (chamber, cutval, cut, cutEntries/nocutEntries) )
    
    f.write( "--------------- %s Done ---------------\n\n" % cutval )
        
#h_tight = makeTH1(filename, gentree, "tight", binning, plotvar, "muon_isME0Muon && muon.Pt()>5 && muon_ME0noRecHit > 3 && muon_ME0deltaX < 1 && muon_ME0deltaY < 4 && muon_ME0deltaDXDZ < 0.05 && muon_ME0deltaDYDZ < 0.3 && muon_ME0pullX < 1.5 && muon_ME0pullY < 3 && muon_ME0dPhi < 0.05" )
#h_loose = makeTH1(filename, gentree, "loose", binning, plotvar, "muon_isME0Muon && muon.Pt()>5 && muon_ME0noRecHit > 2 && muon_ME0deltaX < 2 && muon_ME0deltaY < 8 && muon_ME0deltaDXDZ < 0.075 && muon_ME0deltaDYDZ < 0.4 && muon_ME0pullX < 3 && muon_ME0pullY < 8 && muon_ME0dPhi < 0.1" )

"""
h_tight = makeTH1(filename, gentree, "tight", binning, plotvar, "muon_isME0Muon && muon.Pt()>5 && (muon_ME0pullX < 1.5 || muon_ME0deltaX < 1) && (muon_ME0pullY < 3 || muon_ME0deltaY < 4) && muon_ME0dPhi < 0.05" )
h_loose = makeTH1(filename, gentree, "loose", binning, plotvar, "muon_isME0Muon && muon.Pt()>5 && (muon_ME0pullX < 3 || muon_ME0deltaX < 2) && (muon_ME0pullY < 8 || muon_ME0deltaY < 8)  && muon_ME0dPhi < 0.1" )

h_tight = makeTH1(filename, gentree, "tight", binning, plotvar, "muon_isME0Muon && muon.Pt()>5 && (muon_ME0pullX < 3 || muon_ME0deltaX < 3) && (muon_ME0pullY < 3 || muon_ME0deltaY < 3) && muon_ME0dPhi < 0.1" )
h_loose = makeTH1(filename, gentree, "loose", binning, plotvar, "muon_isME0Muon && muon.Pt()>5 && (muon_ME0pullX < 3 || muon_ME0deltaX < 3) && (muon_ME0pullY < 3 || muon_ME0deltaY < 3) && muon_ME0dPhi < 0.5" )

print "Total Tight Eff = %s" % (h_tight.Integral()/ME0Entries)
print "Total Loose Eff = %s" % (h_loose.Integral()/ME0Entries)
"""

f.close()
