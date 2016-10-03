import FWCore.ParameterSet.Config as cms

flashggDiLeptonDiJetTagTester = cms.EDAnalyzer('FlashggDiLeptonDiJetTagTestAnalyzer',
                                  TagSorter = cms.InputTag('flashggDiLeptonDiJetTagSorter'),
                                  )