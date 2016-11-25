import ROOT
import sys
from array import array
from math import sqrt
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        print "Usage: %s <dir1, dir2, output_dir, legend1, legend2>" % sys.argv[0]  
        sys.exit(1)

#python makeComparisonsPlots.py dir1 dir2 output_dir


dir1 = str(sys.argv[1])
dir2 = str(sys.argv[2])
output_dir = str(sys.argv[3])
legend1 = str(sys.argv[4])
legend2 = str(sys.argv[5])


inputfile1 = TFile(dir1+"/dumptreeHistos.root")
inputfile2 = TFile(dir2+"/dumptreeHistos.root")


histoNames = ["diLeptonDiJetMass", "diLeptonMass", "leadLeptonPt", "subLeadLeptonPt", "leadLeptonEta", "subLeadLeptonEta", "leadLeptonPhi", "subLeadLeptonPhi", "leadJetPt", "subLeadJetPt", "leadJetEta", "subLeadJetEta", "leadJetPhi", "subLeadJetPhi"]#, "diLeptonDiJetSumPt", "leadElePt", "leadMuonPt"]

for histoName in histoNames:
	plot2HistosAndRatio(histoName+"Hist", "", inputfile1, inputfile2, output_dir, "", False, False, 0, 0, False, 0, 0, legend1, legend2, False, False, False, False, False, False, False, title="", xTitle=histoName, yTitle="[a.u.]")
