#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.source = cms.Source("PoolSource",
				fileNames=cms.untracked.vstring(
				# "file:myMicroAODOutputFile_DiLeptonDiJet.root"  
				# "root://node12.datagrid.cea.fr//store/user/gnegro/cmsWR/WR-ToLNu_GEN_SIM_13TeV-2016/WR-2400_ToLNu-1200_ToEEJJ_microAOD_13TeV-2016/161026_083911/0000/dafneMicroAOD_WR_10.root"
				"/store/user/gnegro/cmsWR/cmsWR2016/dafne/DoubleEG/cmsWR2016-dafne-v0-Run2016B-PromptReco-v2/161027_122605/0000/dafneMicroAOD_100.root"
				# "/store/user/gnegro/cmsWR/cmsWR2016/dafne/MuonEG/cmsWR2016-dafne-v0-Run2016B-PromptReco-v2/161027_122809/0000/dafneMicroAOD_10.root"
				#"file:root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/DoubleEG/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-Run2016B-PromptReco-v2/160707_143218/0000/myMicroAODOutputFile_938.root" 
				)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.TFileService = cms.Service("TFileService",
									fileName = cms.string("test.root")
)


from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
   	"HLT_DoubleEle33_CaloIdL_MW_v*",
   	"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
   	"HLT_Mu50_v*",
   	"HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v*",
   	"HLT_Ele27_WPTight_Gsf_v*",
   	"HLT_IsoMu24_v*",
   	"HLT_IsoMu27_v*"
	), 
	throw = False
)


#process.load("dafne.MicroAOD.flashggDiLeptonDiJet_cfi")  ##  import DiLeptonDiJet (producer)

process.load("dafne.Taggers.DiLeptonDiJetDumper_cfi")  ##  import DiLeptonDiJetDumper 

import flashgg.Taggers.dumperConfigTools as cfgTools

#process.DiLeptonDiJetDumper.src = "flashggDiLeptonDiJet"  
process.DiLeptonDiJetDumper.dumpTrees = True
process.DiLeptonDiJetDumper.dumpWorkspace = False
process.DiLeptonDiJetDumper.quietRooFit = True


# split tree, histogram and datasets by process
# process.DiLeptonDiJetDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"
process.DiLeptonDiJetDumper.nameTemplate = "cmsWR_$SQRTS_$LABEL_$SUBCAT" #"massiveNu_$SQRTS_$LABEL_$SUBCAT"


## do not split by process
## process.DiLeptonDiJetDumper.nameTemplate = "minitree_$SQRTS_$LABEL_$SUBCAT"

## define categories and associated objects to dump
cfgTools.addCategory(process.DiLeptonDiJetDumper,
										 "Reject",
										 "0",
										 # "abs(leadingPhoton.superCluster.eta)>=1.4442&&abs(leadingPhoton.superCluster.eta)<=1.566||abs(leadingPhoton.superCluster.eta)>=2.5"
										 # "||abs(subLeadingPhoton.superCluster.eta)>=1.4442 && abs(subLeadingPhoton.superCluster.eta)<=1.566||abs(subLeadingPhoton.superCluster.eta)>=2.5",
											-1 ## if nSubcat is -1 do not store anythings
										 )

# interesting categories 
cfgTools.addCategories(process.DiLeptonDiJetDumper,
											 ## categories definition
											 ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
											 [("all","1",0)
											 # [("EBHighR9","max(abs(leadingPhoton.superCluster.eta),abs(leadingPhoton.superCluster.eta))<1.4442"
											 #   "&& min(leadingPhoton.r9,subLeadingPhoton.r9)>0.94",0), ## EB high R9
											 #  ("EBLowR9","max(abs(leadingPhoton.superCluster.eta),abs(leadingPhoton.superCluster.eta))<1.4442",0), ## remaining EB is low R9
											 #  ("EEHighR9","min(leadingPhoton.r9,subLeadingPhoton.r9)>0.94",0), ## then EE high R9
											 #  ("EELowR9","1",0), ## evereything elese is EE low R9
												],
											 ## variables to be dumped in trees/datasets. Same variables for all categories
											 ## if different variables wanted for different categories, can add categorie one by one with cfgTools.addCategory
											 variables=[ 
														"leadElePt                   := ? isEEJJ ? leadingEle.pt : -999",
														"leadMuonPt                  := ? isMMJJ ? leadingMuon.pt : -999",
														"leadLeptonPt                :=leadingLeptonPt", 
														"subLeadLeptonPt             :=subLeadingLeptonPt",
														"leadLeptonEta               :=leadingLeptonEta",
														"subLeadLeptonEta            :=subLeadingLeptonEta", 
														"leadLeptonPhi               :=leadingLeptonPhi",
														"subLeadLeptonPhi            :=subLeadingLeptonPhi",  
														"leadJetPt                   :=leadingJet.pt",
														"subLeadJetPt                :=subLeadingJet.pt",
														"leadJetEta                  :=leadingJet.eta",
														"subLeadJetEta               :=subLeadingJet.eta",
														"leadJetPhi                  :=leadingJet.phi",
														"subLeadJetPhi               :=subLeadingJet.phi",
														"diLeptonDiJetSumPt          :=sumPt",
														"diLeptonDiJetMass           :=invMass",
														"diLeptonMass                :=diLeptonInvMass",
														"isEEJJ                      :=isEEJJ",
														"isEETT                      :=isEETT",
														"isMMJJ                      :=isMMJJ",
														"isMMTT                      :=isMMTT"
														],
											 ## histograms to be plotted. 
											 ## the variables need to be defined first
											 histograms=[
														# # "subLeadLeptonPt:leadLeptonPt>>ptSubVsLead(180,20,200:180,20,200)",
														# "leadLeptonPt>>leadLeptonPt(100, 0, 2500)",
														# "subLeadLeptonPt>>subLeadLeptonPt(100, 0, 2500)",
														# "leadJetPt>>leadJetPt(100, 0, 2500)",
														# "subLeadJetPt>>subLeadJetPt(100, 0, 2500)",
														# "leadLeptonEta>>leadLeptonEta(100, -4, 4)",
														# "subLeadLeptonEta>>subLeadLeptonEta(100, -4, 4)",
														# "leadJetEta>>leadJetEta(100, -4, 4)",
														# "subLeadJetEta>>subLeadJetEta(100, -4, 4)",
														# "leadLeptonPhi>>leadLeptonPhi(100, -3.2, 3.2)",
														# "subLeadLeptonPhi>>subLeadLeptonPhi(100, -3.2, 3.2)",
														# "leadJetPhi>>leadJetPhi(100, -3.2, 3.2)",
														# "subLeadJetPhi>>subLeadJetPhi(100, -3.2, 3.2)",
														# "diLeptonDiJetMass>>diLeptonDiJetMass(500, 0, 6000)",
														# "diLeptonMass>>diLeptonMass(500, 0, 6400)"
														]
											 )


# process.p1 = cms.Path(
# 		#process.flashggDiLeptonDiJet*process.DiLeptonDiJetDumper  #add producer if microAOD producted without DiLeptonDiJet
# 		process.DiLeptonDiJetDumper   #if microAOD producted with DiLeptonDiJet
# 		)


from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents",1000)
customize(process)

print "processType:", customize.processType, "processId:", customize.processId

if customize.processType == "data":
        print 'data'
        process.DiLeptonDiJetDumper.globalVariables.addTriggerBits = cms.PSet(
            tag = cms.InputTag("TriggerResults::HLT"),
            bits = cms.vstring(
            	"HLT_DoubleEle33_CaloIdL_MW",
			   	"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW",
			   	"HLT_Mu50",
			   	"HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL",
			   	"HLT_Ele27_WPTight_Gsf",
			   	"HLT_IsoMu24",
			   	"HLT_IsoMu27"
            )
        )
        process.p = cms.Path( process.hltHighLevel * process.DiLeptonDiJetDumper )
else:
        print 'NOT data'
        process.p = cms.Path( process.DiLeptonDiJetDumper )
        process.DiLeptonDiJetDumper.globalVariables.puReWeight = True