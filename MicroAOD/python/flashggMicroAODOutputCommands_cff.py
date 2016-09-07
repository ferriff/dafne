import FWCore.ParameterSet.Config as cms


microAODmassiveNuOutputCommand = cms.untracked.vstring("drop *",
													"keep recoVertexs_offlineSlimmedPrimaryVertices_*_*", 
													"keep *_reducedEgamma_reducedSuperClusters_*",
													"keep *_reducedEgamma_*PhotonCores_*",
													"keep *_slimmedMETs_*_*",
													"keep *_slimmedMETsNoHF_*_*",
													"keep *_*Rho*_*_*",
													"keep *_offlineBeamSpot_*_*",
													"keep *_TriggerResults_*_*",
													"keep *_eventCount_*_*",
													"keep *_weightsCount_*_*",
													"keep *_generator_*_*",
													"keep *_slimmedGenJets_*_*",                                            
													"keep *_slimmedAddPileupInfo_*_*", 
													"keep *GsfElectronCore*_*_*_*", 
													"keep recoGenParticles_prunedGenParticles_*_*", 
													"keep patPackedGenParticles_packedGenParticles_*_*", 
													"keep *_selectedPatTrigger_*_*",                                                     
													"keep *_slimmedMuons_*_*", 
													"keep *_slimmedElectrons_*_*",
													"keep *_slimmedJets_*_*",
													"keep *_flashggDiLeptonDiJet_*_*"
													)
