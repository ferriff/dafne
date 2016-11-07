import ROOT
import sys
from array import array
from math import sqrt, floor
from ROOT import *
from DrawingAndComparisonFunctions import *



if len(sys.argv) == 1:
        print "Usage: %s <input_file, output_dir, output_fileName>" % sys.argv[0]  
        sys.exit(1)
#python Analysis/makePlotsfromMicroAOD.py output_numEvent1000.root Analysis/results/prova dumpTreeHistos


input_file = str(sys.argv[1])
output_dir = str(sys.argv[2])
output_fileName = str(sys.argv[3])


tree = TTree()
f = TFile(input_file)

dir = f.Get(input_file+":/DiLeptonDiJetDumper/trees")
dir.GetObject("cmsWR_13TeV_all",tree)


branchesToPlot = ["Pt","Eta","Phi","Mass"]

makeGenPlotFromTree(tree, output_dir, output_fileName, branchesToPlot, False, 0, 0)




