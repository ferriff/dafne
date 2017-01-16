import ROOT
import sys
from array import array
from math import sqrt, floor
from ROOT import *
from DrawingAndComparisonFunctions import *


if len(sys.argv) == 1:
        # print "Usage: %s <input_file, output_dir, output_fileName, WRmass>" % sys.argv[0]  
        print "Usage: %s <input_file, output_dir, output_fileName>" % sys.argv[0]  
        sys.exit(1)
#python Analysis/makePlotsfromMiniTrees.py /afs/cern.ch/user/g/gnegro/work/NuAnalysis-flashgg/CMSSW_8_0_20/src/flashgg/MetaData/scripts/cmsWR2016-ReReco-miniTrees/output_WR-1600_ToLNu-800_ToEEJJ_13TeV-2016.root Analysis/miniTrees/WR-1600_ToLNu-800_ToEEJJ distributions_WR-1600_ToLNu-800_ToEEJJ -b


input_file = str(sys.argv[1])
output_dir = str(sys.argv[2])
output_fileName = str(sys.argv[3])
# WRmass = float(sys.argv[4])


tree = TTree()
f = TFile(input_file)

dir = f.Get(input_file+":/analysisTree")  
dir.GetObject("event",tree)    

branches = tree.GetListOfBranches() 
#print branches.GetEntries()


distributionsNames = ["nvtx", "_eta", "_phi", "_invMass", "_pt"]
HEEPvariablesNames = ["etaSC", "isEcalDriven", "dEtaIn", "dPhiIn", "hOverE", "full5x5_sigmaIetaIeta", "full5x5_E2x5_Over_E5x5", "full5x5_E1x5_Over_E5x5", "EmHadDepth1Iso", "innerLayerLostHits", "dxy"] #, "ptTracksIso"]



f_output = TFile(output_dir+"/"+output_fileName+".root","recreate")

for i in range(branches.GetEntries()):  
	branch = branches.At(i)
	branchName = branch.GetName()	
	histName = branchName + "Hist"

	if "_pt" in branchName and not "_ptTrack" in branchName:  
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, 1000, 100, False)

	if "diLeptonDiJet_invMass" in branchName:  
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, 4000, 100, False)

	if "diLepton_invMass" in branchName:  
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, 2000, 100, False)

	if "diJet_invMass" in branchName:  
		dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0, 2000, 100, False)

	for distributionsName in distributionsNames:
		if distributionsName in branchName:
			dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, False, 0, 0, 0, False)
		




f_output_HEEPvariable = TFile(output_dir+"/"+output_fileName+"_HEEPvariables.root","recreate")

for i in range(branches.GetEntries()):  
	branch = branches.At(i)
	branchName = branch.GetName()	
	histName = branchName + "Hist"

	for HEEPvariablesName in HEEPvariablesNames:
		if HEEPvariablesName in branchName:

			if "dEtaIn" in HEEPvariablesName or "dPhiIn" in HEEPvariablesName:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, -0.2, 0.2, 100, False)
			elif "hOverE" in HEEPvariablesName:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0., 2., 40, False)
			elif "full5x5_sigmaIetaIeta" in HEEPvariablesName:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0., 0.05, 100, False)
			elif "full5x5_E2x5_Over_E5x5" in HEEPvariablesName:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0.5, 1.1, 60, False)
			elif "full5x5_E1x5_Over_E5x5" in HEEPvariablesName:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0., 1.1, 100, False)
			elif "EmHadDepth1Iso" in HEEPvariablesName:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0., 200., 200, False)
			elif "dxy" in HEEPvariablesName:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, True, 0., 0.5, 100, False)
			else:
				dumpPlotFromTreeAndEditHisto(tree, branchName, histName, output_dir, False, 0, 0, 0, False)


#dumpPlotFromTreeAndEditHisto(chain, branchName, histName, output_name, zoomX, xmin, xmax, nBins, log):







