import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

import os
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v13')
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_2016_miniAODv2')
else:
    raise Exception,"The default setup for microAODstd.py does not support releases other than 76X and 80X"

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
process.RandomNumberGeneratorService.flashggRandomizedPhotons = cms.PSet(
          initialSeed = cms.untracked.uint32(16253245)
        )


#80x MC
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
	#"/store/mc/RunIISpring16MiniAODv2/GluGluHToGG_M-125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/083E0403-2825-E611-89C4-0CC47A6C1038.root"
	#signal 2015
	#"/store/mc/RunIISpring15MiniAODv2/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/02126566-9571-E511-8171-00259073E3FA.root"
	#"/store/mc/RunIISpring15MiniAODv2/WRToNuEToEEJJ_MW-1000_MNu-500_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/EE853593-5171-E511-B80F-00259073E3F0.root"
	# "/store/mc/RunIISpring15MiniAODv2/WRToNuEToEEJJ_MW-1200_MNu-600_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/402595D5-8A72-E511-B9C7-90E6BA5CAE1C.root"
	# "/store/mc/RunIISpring15MiniAODv2/WRToNuMuToMuMuJJ_MW-800_MNu-400_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/4E112AD7-5071-E511-8DB6-00259073E456.root"
	# "/store/mc/RunIISpring15MiniAODv2/WRToNuMuToMuMuJJ_MW-1000_MNu-500_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/A27E6B54-2072-E511-9065-0002C90F8088.root"
	# "/store/mc/RunIISpring15MiniAODv2/WRToNuMuToMuMuJJ_MW-1200_MNu-600_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/D8B78FA4-6175-E511-BBA8-003048C56FD8.root"
	#bkg
	"/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/04C2B31D-4441-E611-AF44-24BE05CE1E01.root"
	#"/store/mc/RunIISpring16MiniAODv2/WZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/1E02599D-861B-E611-9524-A4BADB22A4AE.root"
	#"/store/mc/RunIISpring16MiniAODv2/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/1C9FF7D0-E21A-E611-A2D2-B083FED42488.root"
	#"/store/mc/RunIISpring16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v4/00000/20923337-ED2B-E611-83FE-02163E013D06.root"
	#"/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/00709321-002A-E611-A59B-0CC47A74527A.root"
	#"/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/0251DBB7-201B-E611-8653-0CC47A4F1C2E.root"
	#"/store/mc/RunIISpring16MiniAODv2/DYToEE_NNPDF30_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/00503B81-D21A-E611-9222-008CFA00317C.root"
	#"/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/90000/1223ED32-0333-E611-B0A9-549F35AF44F0.root"
	#"/store/mc/RunIISpring16MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/00F10062-152C-E611-99FF-549F358EB789.root"
	))

#80x data
# process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
	#"/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/1E5ABF54-E019-E611-AAED-02163E01293F.root"
	#"/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/02D9C19F-571A-E611-AD8E-02163E013732.root"
	#"/store/data/Run2016B/MuonEG/MINIAOD/PromptReco-v2/000/273/158/00000/26281378-291A-E611-AE69-02163E011E9B.root"
	# ))


process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now
# process.MessageLogger.suppressWarning.extend(['SimpleMemoryCheck','MemoryCheck']) # this would have been better...

# Uncomment the following if you notice you have a memory leak
# This is a lightweight tool to digg further
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#                                        ignoreTotal = cms.untracked.int32(1),
#                                        monitorPssAndPrivate = cms.untracked.bool(True)
#                                       )
 
process.load("dafne/MicroAOD/flashggMicroAODSequence_DiLeptonDiJet_cff")  

# NEEDED FOR ANYTHING PRIOR TO reMiniAOD
#process.weightsCount.pileupInfo = "addPileupInfo"

# from flashgg.MicroAOD.flashggMicroAODOutputCommands_cff import microAODDefaultOutputCommand
# process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myMicroAODOutputFile_DiLeptonDiJet.root'),
#                                outputCommands = microAODDefaultOutputCommand
#                                )

from dafne.MicroAOD.flashggMicroAODOutputCommands_cff import microAODmassiveNuOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myMicroAODOutputFile_DiLeptonDiJet.root'),
                               outputCommands = microAODmassiveNuOutputCommand
                               )


# All jets are now handled in MicroAODCustomize.py
# Switch from PFCHS to PUPPI with puppi=1 argument (both if puppi=2)



process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.globalTightHalo2016Filter_cfi')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')

process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)


process.p = cms.Path(process.flashggMicroAODSequenceDiLeptonDiJet)
process.e = cms.EndPath(process.out)

# Uncomment these lines to run the example commissioning module and send its output to root
#process.commissioning = cms.EDAnalyzer('flashggCommissioning',
#                                       PhotonTag=cms.untracked.InputTag('flashggPhotons'),
#                                       DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
#                                       VertexTag=cms.untracked.InputTag('offlineSlimmedPrimaryVertices')
#)
#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("commissioningTree.root")
#)
#process.p *= process.commissioning


from flashgg.MicroAOD.MicroAODCustomize import customize
customize(process)

if "DY" in customize.datasetName or "SingleElectron" in customize.datasetName or "DoubleEG" in customize.datasetName:
  customize.customizeHLT(process)
