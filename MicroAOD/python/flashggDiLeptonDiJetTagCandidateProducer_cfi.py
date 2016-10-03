import FWCore.ParameterSet.Config as cms

flashggDiLeptonDiJetTagCandidateProducer = cms.EDProducer('FlashggDiLeptonDiJetTagCandidateProducer',
                                             TagSorter = cms.untracked.InputTag('flashggDiLeptonDiJetTagSorter')
                             )
