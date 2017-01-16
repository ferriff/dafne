import ROOT
import sys
from array import array
from math import sqrt, floor
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        print "Usage: %s <input_file, output_dir, output_fileName, WRmass>" % sys.argv[0]  
        sys.exit(1)
#python Analysis/makePlotsfromMicroAOD.py output_numEvent1000.root Analysis/results/prova dumpTreeHistos


input_file = str(sys.argv[1])
output_dir = str(sys.argv[2])
output_fileName = str(sys.argv[3])
WRmass = float(sys.argv[4])


tree = TTree()
f = TFile(input_file)

dir = f.Get(input_file+":/DiLeptonDiJetDumper/trees")
dir.GetObject("cmsWR_13TeV_all",tree)

branches = tree.GetListOfBranches() 
#print branches.GetEntries()

f_output = TFile(output_dir+"/"+output_fileName+".root","recreate")

for i in range(branches.GetEntries()):  
	branch = branches.At(i)
	branchName = branch.GetName()	
	histName = branchName + "Hist"

	if "Pt" in branchName :	#  and not "Sum" in branch Name: altrimenti prende anche diLeptonDiJetSumPt
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, 1500, 100, False)

	if "Eta" in branchName:
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, -2.5, 2.5, 100, False)

	if "Phi" in branchName:
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, -3.5, 3.5, 100, False)

	if "diLeptonMass" in branchName:
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, WRmass, 100, False) 

	if "diLeptonDiJetMass" in branchName:
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, WRmass*2, 100, False) 

	if "diLeptonDiJetSumPt" in branchName:
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, WRmass*2, 100, False) 
