import ROOT
import sys
import tdrstyle,CMS_lumi
from array import array
from math import sqrt, floor
from ROOT import *


tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = ""
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"	



def dumpPlotFromTree(chain, branchName, output_name):
	canv = TCanvas(branchName, branchName, 600, 600)
	canv.cd()
	chain.Draw(branchName)
	# canv.SaveAs(output_name+"_"+branchName+".png","recreate")
	canv.SaveAs(output_name+"/"+branchName+".png","recreate")
	canv.Close()



def dumpPlotFromTreeAndEditHisto(chain, branchName, histName, output_name, zoomX, xmin, xmax, nBins, log):
	canv = TCanvas(branchName, branchName, 600, 600)
	canv.cd()
	
	if zoomX:
		chain.Draw(branchName+">>"+histName+"("+str(nBins)+", "+str(xmin)+", "+str(xmax)+")") #to specify binning x axis
	else:
		chain.Draw(branchName+">>"+histName)  #a new histogram called histName is created and kept in the current directory
	
	hist = TH1F()
	hist = gDirectory.Get(histName)  #to retrieve the histogram

	# hist.SetTitle(branchName)
	# hist.SetLineColor(1)
	# hist.SetLineWidth(3)

	hist.GetXaxis().SetTitle(branchName)
	hist.GetYaxis().SetTitle("Events")

	hist.Write()
	hist.Draw()

	if log:
		canv.SetLogy()  #visible only on .png plots

	canv.SaveAs(output_name+"/"+branchName+".png","recreate")
	canv.Close()



def makeGenPlotFromTree(chain, output_name, outputFile_name, branchesToPlot, zoomX, xmin, xmax, nBins, log):
	branches = chain.GetListOfBranches() 
	#print branches.GetEntries()

	f_output = TFile(output_name+"/"+outputFile_name+".root","recreate")

	for i in range(branches.GetEntries()):  
		branch = branches.At(i)
		branchName = branch.GetName()	
		histName = branchName + "Hist"

		for branchToPlot in branchesToPlot:
			if branchToPlot in branchName:	
				dumpPlotFromTreeAndEditHisto(chain, branchName, histName, output_name, zoomX, xmin, xmax, nBins, log)



def makeProfile(histo2DName, etaRegion, inputfile, output_name, suff, ymin, ymax, yTitle, pfxRange, pfxMin, pfxMax, zoomY, yRangeMin, yRangeMax, xRangeMin, xRangeMax, xTitle = "#eta", xbins = 100 , xmin = -4., xmax = 4., saveCanvas = False):

		histo_vsEta = TH2D()
		inputfile.GetObject(str(histo2DName)+str(etaRegion),histo_vsEta)

		histo_vsEta_profile = TProfile(str(histo2DName)+str(etaRegion)+"_profile", str(histo2DName)+str(etaRegion)+"_profile", xbins, xmin, xmax, ymin, ymax)


		if pfxRange:
			histo_vsEta_profile = histo_vsEta.ProfileX("_pfx", pfxMin, pfxMax)    
		else:
			histo_vsEta_profile = histo_vsEta.ProfileX() 			
	
		histo_vsEta_profile.GetXaxis().SetTitle(xTitle)
		histo_vsEta_profile.GetYaxis().SetTitle(yTitle)

		if zoomY:
			histo_vsEta_profile.GetYaxis().SetRangeUser(yRangeMin, yRangeMax) 

		histo_vsEta_profile.GetXaxis().SetRangeUser(xRangeMin, xRangeMax)


		if saveCanvas:
			c = TCanvas(str(histo2DName)+str(etaRegion)+"_profile",str(histo2DName)+str(etaRegion)+"_profile")
			c.cd()
			CMS_lumi.CMS_lumi(c, 4, 33)
			histo_vsEta_profile.Draw()   
			c.SaveAs(str(output_name)+str(suff)+"/"+str(histo2DName)+str(etaRegion)+"_profile"+str(sys.argv[4])+".png")

		histo_vsEta_profile.Write()



def plot2HistosAndRatio(histoName, etaRegion, inputfile1, inputfile2, output_name, suff, dataMC, zoomX, xRangeMin, xRangeMax, zoomY, yRangeMin, yRangeMax, leg1_name, leg2_name, leftLegends, profile, log, doRebin, doRebinVariableBinSize, doRestrictedIntegral, MCreweighted, title="", xTitle="", yTitle=""):

	if profile:
		histo1 = TProfile()
		histo2 = TProfile()
		inputfile1.GetObject(str(histoName)+str(etaRegion)+"_pfx",histo1)
		inputfile2.GetObject(str(histoName)+str(etaRegion)+"_pfx",histo2)
	else:
		histo1 = TH1D() #TH1F()
		histo2 = TH1D() #TH1F()
		inputfile1.GetObject(str(histoName)+str(etaRegion),histo1)
		inputfile2.GetObject(str(histoName)+str(etaRegion),histo2)

		histo1.Sumw2()
		histo2.Sumw2()

		if doRestrictedIntegral == "True" and zoomX:
			integral_histo1 = histo1.Integral(histo1.FindBin(xRangeMin), histo1.FindBin(xRangeMax))
			integral_histo2 = histo2.Integral(histo2.FindBin(xRangeMin), histo2.FindBin(xRangeMax))
		else:
			integral_histo1 = histo1.Integral()   
			integral_histo2 = histo2.Integral()		

		histo1.Scale(1/integral_histo1)
		histo2.Scale(1/integral_histo2)


	if profile:
		c1 = TCanvas(str(histoName)+str(etaRegion)+"_pfx"+" Comparison",str(histoName)+str(etaRegion)+"_pfx"+" Comparison")
	else:
		c1 = TCanvas(str(histoName)+str(etaRegion)+" Comparison",str(histoName)+str(etaRegion)+" Comparison")

	c1.Range(0,0,1,1)


	if doRebin: 
		print "doing rebinning"
		histo1.Rebin()
		histo2.Rebin()


	histo1_top = histo1.Clone()
	histo2_top = histo2.Clone()



	if doRebinVariableBinSize == "True":
		binLowEdges = []	
		binLowEdges.append(histo1.GetBinLowEdge(1))

		for i in range(2, histo1.GetNbinsX()+1):	
			#if (histo1.GetBinContent(i) > 1.e-4):
			if (histo1.GetBinContent(i) > 1.e-3):
				binLowEdge = histo1.GetBinLowEdge(i)
				binLowEdges.append(binLowEdge)
			else:
				continue

		binLowEdges.append(histo1.GetBinLowEdge(histo1.GetNbinsX())+histo1.GetBinWidth(histo1.GetNbinsX()))

		binLowEdgesArray = array('d', (i for i in binLowEdges)) 
   
		histo1_bottom = histo1.Rebin(len(binLowEdges)-1, "histo1 rebinned", binLowEdgesArray)
		histo2_bottom = histo2.Rebin(len(binLowEdges)-1, "histo2 rebinned", binLowEdgesArray)    

	else:
		histo1_bottom = histo1.Clone()
		histo2_bottom = histo2.Clone()
	

	errors_data = []
	for i in range(histo1_bottom.GetNbinsX()+1):	
		error_data = 0.
		if (abs(histo1_bottom.GetBinContent(i)) > 1.e-32):
			error_data = histo1_bottom.GetBinError(i) / histo1_bottom.GetBinContent(i)   #prendo errori di histo1 rebinnato prima di ratio
		errors_data.append(error_data)



	histo1_bottom.Divide(histo2_bottom)


	for i in range(histo1_bottom.GetNbinsX()+1):	
		histo1_bottom.SetBinError(i,errors_data[i] * histo1_bottom.GetBinContent(i))

	
	if profile:
		histo_MCerr = TProfile(histo1_bottom)
	else:
		histo_MCerr = TH1D(histo1_bottom) #TH1F(histo1_bottom)

	histo_MCerr.SetName("histo_ratio_MCerrors")
	for i in range(histo2_bottom.GetNbinsX()+1):	
		ey = 0.
		if (abs(histo2_bottom.GetBinContent(i)) > 1.e-32):
			ey = histo2_bottom.GetBinError(i) / histo2_bottom.GetBinContent(i)

		histo_MCerr.SetBinError(i,ey * histo1_bottom.GetBinContent(i))
		histo_MCerr.SetBinContent(i,1)


	yTitle2 = "ratio" #bottom plot y axis title

	if profile:
		defaultRatioYmin = 0.8
		defaultRatioYmax = 1.2
	else:
		defaultRatioYmin = 0.
		defaultRatioYmax = 2.


	#Bottom plot
	c1_1 = TPad("c1_1", "newpad",0.01,0.01,0.99,0.32)
	#c1_1.DrawFrame(xRangeMin, defaultRatioYmin, xRangeMax, defaultRatioYmax, ";"+xTitle+";"+yTitle2)
	c1_1.Draw()
	c1_1.cd()
	c1_1.SetTopMargin(0.01)
	c1_1.SetBottomMargin(0.3)
	c1_1.SetRightMargin(0.05) #0.1)
  	c1_1.SetFillStyle(0)
	#c1_1.SetGridy(5)


	histo1_bottom.SetMinimum(defaultRatioYmin)
	histo1_bottom.SetMaximum(defaultRatioYmax)
	histo1_bottom.GetYaxis().SetNdivisions(5)
	histo1_bottom.SetTitle(";"+xTitle+";"+yTitle2)
	histo1_bottom.GetXaxis().SetTitleSize(0.14)
	histo1_bottom.GetXaxis().SetLabelSize(0.14)
	histo1_bottom.GetYaxis().SetLabelSize(0.11)
	histo1_bottom.GetYaxis().SetTitleSize(0.14)
	histo1_bottom.GetYaxis().SetTitleOffset(0.4)#0.28)

	histo1_bottom.SetLineColor(1)
	histo1_bottom.SetMarkerStyle(8)
	histo1_bottom.SetMarkerSize(0.5)
	histo1_bottom.SetMarkerColor(1)
	histo1_bottom.Draw("E1P")

	if MCreweighted == "True":
		histo_MCerr.SetFillColor(4)
	else:
		histo_MCerr.SetFillColor(2)
	histo_MCerr.SetFillStyle(3001)
	#histo_MCerr.Draw("E2same")


	if zoomX:
		l = TLine(xRangeMin, 1., xRangeMax, 1.)
	else:
		l = TLine(histo1_bottom.GetXaxis().GetXmin(), 1., histo1_bottom.GetXaxis().GetXmax(), 1.)

	l.SetLineColor(1)
	l.Draw("same")




 	#Top Plot
	c1.cd()
	c1_2 = TPad("c1_2", "newpad",0.01,0.33,0.99,0.99)
	c1_2.Draw()
	c1_2.cd()
	c1_2.SetTopMargin(0.1)
	c1_2.SetBottomMargin(0.01)
	c1_2.SetRightMargin(0.05)#0.1)
	c1_2.SetFillStyle(0)


	histo1_top.SetLabelSize(0.0)   #non stampa label su asse x
	histo1_top.GetXaxis().SetTitleSize(0.00)
	histo1_top.GetYaxis().SetLabelSize(0.06)#0.07)
	histo1_top.GetYaxis().SetTitleSize(0.07)#0.08)
	histo1_top.GetYaxis().SetTitleOffset(1.) #0.76)
	histo1_top.SetTitle(title+";;"+yTitle)

	if dataMC:
		histo1_top.SetMarkerColor(1)
		histo1_top.SetLineColor(1)
		histo1_top.SetMarkerStyle(8)
		histo1_top.SetMarkerSize(0.5) 
		histo1_top.Draw("E1P") 
	else:
		histo1_top.SetLineColor(4)
		histo1_top.SetFillColor(kBlue-3)

		if profile:
			histo1_top.SetMarkerStyle(8)
			histo1_top.SetMarkerSize(0.5) 
			histo1_top.SetMarkerColor(4)  
			histo1_top.Draw("E1P") 
		else:
			histo1_top.SetFillStyle(3001)
			histo1_top.Draw("HIST") 


	histo2_top.SetLineColor(2)
	histo2_top.SetFillColor(kRed-3)

	if profile:
		histo2_top.SetMarkerStyle(8)
		histo2_top.SetMarkerSize(0.5) 
		histo2_top.SetMarkerColor(2)
		histo2_top.Draw("same")
	else:
		histo2_top.SetFillStyle(3001)
		histo2_top.Draw("HISTsame")



	if zoomY:		
		histo1_top.GetYaxis().SetRangeUser(yRangeMin,yRangeMax)

	if zoomX:		
		histo1_top.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)
		histo1_bottom.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)   #non lo prende per RebinVariableSize


	if profile:
		CMS_lumi.CMS_lumi(c1_2, 4, 11)
		c1.cd()
		leg1 = TLegend(0.6,0.85,0.8,0.9)
		leg2 = TLegend(0.6,0.8,0.8,0.85)

	else:
		if leftLegends:
			CMS_lumi.CMS_lumi(c1_2, 4, 11)
			c1.cd()
			leg1 = TLegend(0.2,0.75,0.3,0.8)
			leg2 = TLegend(0.2,0.7,0.3,0.75)
		else:
			CMS_lumi.CMS_lumi(c1_2, 4, 33)
			c1.cd()
			leg1 = TLegend(0.6,0.75,0.8,0.8)
			leg2 = TLegend(0.6,0.7,0.8,0.75)



	leg1.AddEntry(0, leg1_name, '')
	leg2.AddEntry(0, leg2_name, '')

	leg1.SetBorderSize(0)
	leg1.SetFillColor(0)
	leg1.SetTextSize(0.028) 

	if dataMC:
		leg1.SetTextColor(1)
	else: 
		leg1.SetTextColor(4)  
	leg1.Draw('same')

	leg2.SetBorderSize(0)
	leg2.SetFillColor(0)
	leg2.SetTextSize(0.028)
	leg2.SetTextColor(2)	
	leg2.Draw('same')


	if log:
		c1_2.SetLogy()


	if profile:
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_pfx"+"_"+str(output_name)+str(suff)+".png")
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_pfx"+"_"+str(output_name)+str(suff)+".root")
	else:
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+".png")
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+".root")
		# c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_"+str(output_name)+str(suff)+".png")
		# c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_"+str(output_name)+str(suff)+".root")

	c1.Close()


def plot3HistosAndRatioFirst2(histoName, etaRegion, inputfile1, inputfile2, inputfile3, output_name, suff, zoomX, xRangeMin, xRangeMax, zoomY, yRangeMin, yRangeMax, leg1_name, leg2_name, leg3_name, leftLegends, profile, log, doRebin, doRebinVariableBinSize, doRestrictedIntegral, MCreweighted, title="", xTitle="", yTitle=""):

	if profile:
		histo1 = TProfile()
		histo2 = TProfile()
		histo3 = TProfile()
		inputfile1.GetObject(str(histoName)+str(etaRegion)+"_pfx",histo1)
		inputfile2.GetObject(str(histoName)+str(etaRegion)+"_pfx",histo2)
		inputfile3.GetObject(str(histoName)+str(etaRegion)+"_pfx",histo3)
	else:
		histo1 = TH1F()
		histo2 = TH1F()
		histo3 = TH1F()
		inputfile1.GetObject(str(histoName)+str(etaRegion),histo1)
		inputfile2.GetObject(str(histoName)+str(etaRegion),histo2)
		inputfile3.GetObject(str(histoName)+str(etaRegion),histo3)

		histo1.Sumw2()
		histo2.Sumw2()
		histo3.Sumw2()

		if doRestrictedIntegral == "True" and zoomX:
			integral_histo1 = histo1.Integral(histo1.FindBin(xRangeMin), histo1.FindBin(xRangeMax))
			integral_histo2 = histo2.Integral(histo2.FindBin(xRangeMin), histo2.FindBin(xRangeMax))
			integral_histo3 = histo3.Integral(histo3.FindBin(xRangeMin), histo3.FindBin(xRangeMax))
		else:
			integral_histo1 = histo1.Integral()   
			integral_histo2 = histo2.Integral()		
			integral_histo3 = histo3.Integral()

		histo1.Scale(1/integral_histo1)
		histo2.Scale(1/integral_histo2)
		histo3.Scale(1/integral_histo3)


	if profile:
		c1 = TCanvas(str(histoName)+str(etaRegion)+"_pfx"+" Comparison",str(histoName)+str(etaRegion)+"_pfx"+" Comparison")
	else:
		c1 = TCanvas(str(histoName)+str(etaRegion)+" Comparison",str(histoName)+str(etaRegion)+" Comparison")

	c1.Range(0,0,1,1)


	if doRebin:
		histo1.Rebin()
		histo2.Rebin()
		histo3.Rebin()

	histo1_top = histo1.Clone()
	histo2_top = histo2.Clone()
	histo3_top = histo3.Clone()


	if doRebinVariableBinSize == "True":
		binLowEdges = []	
		binLowEdges.append(histo1.GetBinLowEdge(1))

		for i in range(2, histo1.GetNbinsX()+1):	
			#if (histo1.GetBinContent(i) > 1.e-4):
			if (histo1.GetBinContent(i) > 1.e-3):
				binLowEdge = histo1.GetBinLowEdge(i)
				binLowEdges.append(binLowEdge)
			else:
				continue

		binLowEdges.append(histo1.GetBinLowEdge(histo1.GetNbinsX())+histo1.GetBinWidth(histo1.GetNbinsX()))

		binLowEdgesArray = array('d', (i for i in binLowEdges)) 
   
		histo1_bottom = histo1.Rebin(len(binLowEdges)-1, "histo1 rebinned", binLowEdgesArray)
		histo2_bottom = histo2.Rebin(len(binLowEdges)-1, "histo2 rebinned", binLowEdgesArray)    

	else:
		histo1_bottom = histo1.Clone()
		histo2_bottom = histo2.Clone()
	

	errors_data = []
	for i in range(histo1_bottom.GetNbinsX()+1):	
		error_data = 0.
		if (abs(histo1_bottom.GetBinContent(i)) > 1.e-32):
			error_data = histo1_bottom.GetBinError(i) / histo1_bottom.GetBinContent(i)   #prendo errori di histo1 rebinnato prima di ratio
		errors_data.append(error_data)



	histo1_bottom.Divide(histo2_bottom)


	for i in range(histo1_bottom.GetNbinsX()+1):	
		histo1_bottom.SetBinError(i,errors_data[i] * histo1_bottom.GetBinContent(i))

	
	if profile:
		histo_MCerr = TProfile(histo1_bottom)
	else:
		histo_MCerr = TH1F(histo1_bottom)

	histo_MCerr.SetName("histo_ratio_MCerrors")
	for i in range(histo2_bottom.GetNbinsX()+1):	
		ey = 0.
		if (abs(histo2_bottom.GetBinContent(i)) > 1.e-32):
			ey = histo2_bottom.GetBinError(i) / histo2_bottom.GetBinContent(i)

		histo_MCerr.SetBinError(i,ey * histo1_bottom.GetBinContent(i))
		histo_MCerr.SetBinContent(i,1)


	yTitle2 = "ratio" #bottom plot y axis title

	if profile:
		defaultRatioYmin = 0.8
		defaultRatioYmax = 1.2
	else:
		defaultRatioYmin = 0.
		defaultRatioYmax = 2.


	#Bottom plot
	c1_1 = TPad("c1_1", "newpad",0.01,0.01,0.99,0.32)
	#c1_1.DrawFrame(xRangeMin, defaultRatioYmin, xRangeMax, defaultRatioYmax, ";"+xTitle+";"+yTitle2)
	c1_1.Draw()
	c1_1.cd()
	c1_1.SetTopMargin(0.01)
	c1_1.SetBottomMargin(0.3)
	c1_1.SetRightMargin(0.05) #0.1)
  	c1_1.SetFillStyle(0)
	#c1_1.SetGridy(5)


	histo1_bottom.SetMinimum(defaultRatioYmin)
	histo1_bottom.SetMaximum(defaultRatioYmax)
	histo1_bottom.GetYaxis().SetNdivisions(5)
	histo1_bottom.SetTitle(";"+xTitle+";"+yTitle2)
	histo1_bottom.GetXaxis().SetTitleSize(0.14)
	histo1_bottom.GetXaxis().SetLabelSize(0.14)
	histo1_bottom.GetYaxis().SetLabelSize(0.11)
	histo1_bottom.GetYaxis().SetTitleSize(0.14)
	histo1_bottom.GetYaxis().SetTitleOffset(0.4)#0.28)

	histo1_bottom.SetLineColor(1)
	histo1_bottom.SetMarkerStyle(8)
	histo1_bottom.SetMarkerSize(0.5)
	histo1_bottom.SetMarkerColor(1)
	histo1_bottom.Draw("E1P")

	if MCreweighted == "True":
		histo_MCerr.SetFillColor(4)
	else:
		histo_MCerr.SetFillColor(2)
	histo_MCerr.SetFillStyle(3001)
	#histo_MCerr.Draw("E2same")


	if zoomX:
		l = TLine(xRangeMin, 1., xRangeMax, 1.)
	else:
		l = TLine(histo1_bottom.GetXaxis().GetXmin(), 1., histo1_bottom.GetXaxis().GetXmax(), 1.)

	l.SetLineColor(1)
	l.Draw("same")




 	#Top Plot
	c1.cd()
	c1_2 = TPad("c1_2", "newpad",0.01,0.33,0.99,0.99)
	c1_2.Draw()
	c1_2.cd()
	c1_2.SetTopMargin(0.1)
	c1_2.SetBottomMargin(0.01)
	c1_2.SetRightMargin(0.05)#0.1)
	c1_2.SetFillStyle(0)


	histo1_top.SetLineColor(2)
	histo1_top.SetFillColor(kRed-3)

	histo1_top.SetLabelSize(0.0)   #non stampa label su asse x
	histo1_top.GetXaxis().SetTitleSize(0.00)
	histo1_top.GetYaxis().SetLabelSize(0.06)#0.07)
	histo1_top.GetYaxis().SetTitleSize(0.07)#0.08)
	histo1_top.GetYaxis().SetTitleOffset(1.) #0.76)
	histo1_top.SetTitle(title+";;"+yTitle)

	histo2_top.SetLineColor(4)
	histo2_top.SetFillColor(kBlue-3)

	histo3_top.SetLineColor(3)
	histo3_top.SetFillColor(kGreen-3)


	if profile:
		histo1_top.SetMarkerStyle(8)
		histo1_top.SetMarkerSize(0.5) 
		histo1_top.SetMarkerColor(2)  

		histo2_top.SetMarkerStyle(8)
		histo2_top.SetMarkerSize(0.5) 
		histo2_top.SetMarkerColor(4)

		histo3_top.SetMarkerStyle(8)
		histo3_top.SetMarkerSize(0.5) 
		histo3_top.SetMarkerColor(3)

		histo1_top.Draw("E1P") 
		histo2_top.Draw("same")
		histo3_top.Draw("E1Psame")		

	else:
		histo1_top.SetFillStyle(3001)
		histo1_top.Draw("HIST") 
		histo2_top.SetFillStyle(3001)
		histo2_top.Draw("HISTsame")
		histo3_top.SetFillStyle(3001)
		histo3_top.Draw("HISTsame")


	if zoomY:		
		histo1_top.GetYaxis().SetRangeUser(yRangeMin,yRangeMax)
		#histo2_top.GetYaxis().SetRangeUser(yRangeMin,yRangeMax)

	if zoomX:		
		histo1_top.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)
		#histo2_top.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)	
		histo1_bottom.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)   #non lo prende per RebinVariableSize


	if profile:
		CMS_lumi.CMS_lumi(c1_2, 4, 11)
		c1.cd()
		leg1 = TLegend(0.6,0.85,0.8,0.9)
		leg2 = TLegend(0.6,0.8,0.8,0.85)
		leg3 = TLegend(0.6,0.75,0.8,0.8)

	else:
		if leftLegends:
			CMS_lumi.CMS_lumi(c1_2, 4, 11)
			c1.cd()
			leg1 = TLegend(0.2,0.75,0.3,0.8)
			leg2 = TLegend(0.2,0.7,0.3,0.75)
			leg3 = TLegend(0.2,0.65,0.3,0.7)
		else:
			CMS_lumi.CMS_lumi(c1_2, 4, 33)
			c1.cd()
			leg1 = TLegend(0.6,0.75,0.8,0.8)
			leg2 = TLegend(0.6,0.7,0.8,0.75)
			leg3 = TLegend(0.6,0.65,0.8,0.7)

	leg1.AddEntry(0, leg1_name,'')
	leg2.AddEntry(0, leg2_name,'')
	leg3.AddEntry(0, leg3_name,'')

	leg1.SetBorderSize(0)
	leg1.SetFillColor(0)
	leg1.SetTextSize(0.028) 
	leg1.SetTextColor(2)
	leg1.Draw('same')

	leg2.SetBorderSize(0)
	leg2.SetFillColor(0)
	leg2.SetTextSize(0.028)
	leg2.SetTextColor(4)	
	leg2.Draw('same')

	leg3.SetBorderSize(0)
	leg3.SetFillColor(0)
	leg3.SetTextSize(0.028)
	leg3.SetTextColor(3)	
	leg3.Draw('same')


	if log:
		c1_2.SetLogy()


	if profile:
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_pfx"+"_"+str(output_name)+str(suff)+".png")
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_pfx"+"_"+str(output_name)+str(suff)+".root")
	else:
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_"+str(output_name)+str(suff)+".png")
		c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+"_"+str(output_name)+str(suff)+".root")

	c1.Close()




def plotDataAndMCHistos(histoName, etaRegion, inputfile1, fileList, output_name, suff, zoomX, xRangeMin, xRangeMax, zoomY, yRangeMin, yRangeMax, leg1_name, legendList, leftLegends, log, doRebin, doRestrictedIntegral, colorList, title="", xTitle="", yTitle=""):

	histoList = []

	histo1 = TH1D()
	inputfile1.GetObject(str(histoName)+str(etaRegion),histo1)
	histo1.Sumw2()

	if doRestrictedIntegral == "True" and zoomX:
	# if zoomX:
		integral_histo1 = histo1.Integral(histo1.FindBin(xRangeMin), histo1.FindBin(xRangeMax))
	else:
		integral_histo1 = histo1.Integral()   

	histo1.Scale(1/integral_histo1)


	for file in fileList: 
		histo = TH1D()
		file.GetObject(str(histoName)+str(etaRegion),histo)
		histo.Sumw2()

		if doRestrictedIntegral == "True" and zoomX:
		# if zoomX:
			integral_histo = histo.Integral(histo.FindBin(xRangeMin), histo.FindBin(xRangeMax))
		else:
			integral_histo = histo.Integral()   

		# histo.Scale(1/integral_histo)
		histo.Scale(1/integral_histo1)

		histoList.append(histo)


	c1 = TCanvas(str(histoName)+str(etaRegion)+" Comparison",str(histoName)+str(etaRegion)+" Comparison")
	c1.cd()

	if doRebin:
		histo1.Rebin()
		for histo in histoList:
			histo.Rebin()

	histo1.SetLineColor(1)
	histo1.SetMarkerColor(1)
	histo1.SetMarkerStyle(8)
	histo1.SetMarkerSize(0.5) 
	# histo1.GetXaxis().SetTitle(xTitle)
	# histo1.GetYaxis().SetTitle(yTitle)
	histo1.Draw("E1P") 

	for histo,color in zip(histoList, colorList):
		histo.SetLineColor(color)
		histo.SetFillColor(color)
		histo.Draw("HISTsame")

	if zoomY:		
		histo1.GetYaxis().SetRangeUser(yRangeMin,yRangeMax)

	if zoomX:		
		histo1.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)
	

	if leftLegends:
		CMS_lumi.CMS_lumi(c1, 4, 11)
		leg = TLegend(0.15,0.6,0.35,0.8)
	else:
		CMS_lumi.CMS_lumi(c1, 4, 33)
		leg = TLegend(0.65,0.6,0.85,0.8)

	leg.AddEntry(histo1, leg1_name,'P')
	for histo,leg_name in zip(histoList, legendList):
		leg.AddEntry(histo, leg_name,'F')

	leg.SetBorderSize(0)
	leg.SetFillColor(0)
	leg.SetTextSize(0.028) 
	leg.SetTextColor(1)
	leg.Draw('same')


	if log:
		c1.SetLogy()

	c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+".png")
	c1.SaveAs(str(output_name)+str(suff)+"/"+str(histoName)+str(etaRegion)+".root")

	c1.Close()
