import FWCore.ParameterSet.Config as cms

flashggDiLeptonDiJetTagSorter = cms.EDProducer('FlashggDiLeptonDiJetTagSorter',
                                  DiLeptonDiJetTag = cms.InputTag('flashggDiLeptonDiJets'), 
                                  # Top of list is highest priority
                                  # Optionally can add category ranges if priority depends on category number
                                  TagPriorityRanges = cms.VPSet(
                                                                 cms.PSet(TagName = cms.InputTag('flashggDiEleDiJetTag')), 
                                                                 cms.PSet(TagName = cms.InputTag('flashggDiMuDiJetTag')),   
                                                                 cms.PSet(TagName = cms.InputTag('flashggDiEleDiTrackTag')),     
                                                                 cms.PSet(TagName = cms.InputTag('flashggDiMuDiTrackTag'))
                                                                ),
                                  MassCutUpper=cms.double(180.), #??
                                  MassCutLower=cms.double(100),  #??
                                  MinObjectWeightException = cms.double(0.1),
                                  MaxObjectWeightException = cms.double(10.),
                                  MinObjectWeightWarning = cms.double(0.5),
                                  MaxObjectWeightWarning = cms.double(2.),
                                  StoreOtherTagInfo = cms.bool(False),
                                  BlindedSelectionPrintout = cms.bool(False),
                                  Debug = cms.untracked.bool(False)
                                  )
