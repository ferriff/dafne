import FWCore.ParameterSet.Config as cms


emptyBins = cms.PSet(
	variables = cms.vstring("1"),
	bins = cms.VPSet()
	)


scalesAndSmearingsPrefix = cms.string("EgammaAnalysis/ElectronTools/data/Winter_2016_reReco_v1_ele")


MCSmearHighR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet2D"),
		Label = cms.string("MCSmearHighR9EB"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the photon object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)

MCSmearLowR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet2D"),
		Label = cms.string("MCSmearLowR9EB"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the photon object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)

MCSmearHighR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),  
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet2D"),
		Label = cms.string("MCSmearHighR9EE"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the photon object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)

MCSmearLowR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronSmearStochasticEGMTool"),
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet2D"),
		Label = cms.string("MCSmearLowR9EE"),
		FirstParameterName = cms.string("Rho"),
		SecondParameterName = cms.string("Phi"),
		CorrectionFile = scalesAndSmearingsPrefix,
		NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
							secondVar = cms.vint32(0,0,1,-1)),
		OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		# has to match the labels embedded in the photon object as
		# defined e.g. in dafne/MicroAOD/python/flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		#           or in flashgg/MicroAOD/python/flashggRandomizedElectronProducer_cff.py (if at MicroAOD prod.)
		# RandomLabel = cms.string("rnd_g_E"), #for flashggRandomizedElectronProducer_cff.py
		RandomLabel = cms.string("smearE"),    #for flashggRandomizedElectronForDiLeptonDiJetProducer_cff.py
		Debug = cms.untracked.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		ApplyCentralValue = cms.bool(True)
		)


MCScaleHighR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet"),
		Label = cms.string("MCScaleHighR9EB"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)

MCScaleLowR9EB_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet"),
		Label = cms.string("MCScaleLowR9EB"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)<1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)

MCScaleHighR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet"),
		Label = cms.string("MCScaleHighR9EE"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)

MCScaleLowR9EE_EGM = cms.PSet( ElectronMethodName = cms.string("FlashggElectronScaleEGMTool"),
		MethodName = cms.string("FlashggElectronFromDiLeptonDiJet"),
		Label = cms.string("MCScaleLowR9EE"),
		NSigmas = cms.vint32(-1,1),
		OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)>=1.5"),
		BinList = emptyBins,
		CorrectionFile = scalesAndSmearingsPrefix,
		ApplyCentralValue = cms.bool(False),
		ExaggerateShiftUp = cms.bool(False),
		Debug = cms.untracked.bool(False)
		)




flashggDiLeptonDiJetSystematics = cms.EDProducer('FlashggDiLeptonDiJetSystematicProducer',
		src = cms.InputTag("flashggDiLeptonDiJet"),
		SystMethods2D = cms.VPSet(),
		# the number of syst methods matches the number of nuisance parameters
		# assumed for a given systematic uncertainty and is NOT required
		# to match 1-to-1 the number of bins above.
		SystMethods = cms.VPSet(
				MCSmearHighR9EB_EGM,
				MCSmearLowR9EB_EGM,
				MCSmearHighR9EE_EGM,
				MCSmearLowR9EE_EGM,
				MCScaleHighR9EB_EGM,
				MCScaleLowR9EB_EGM,
				MCScaleHighR9EE_EGM,
				MCScaleLowR9EE_EGM
		)
)