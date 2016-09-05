import FWCore.ParameterSet.Config as cms


microAODmassiveNuOutputCommand = cms.untracked.vstring("drop *",
													"keep *_flashgg*_*_*",
													"drop *_flashggVertexMap*_*_*",                                                      
													"drop *_flashggPuppi*_*_*", # this part drop all the tools used to build puppi jets
													"drop *_flashggPhotons_*_*", # Only keep the copies with random numbers added
													"drop *_flashggPrunedGenParticles_*_*",                                     
													"drop floatedmValueMap_electronMVAValueMapProducer_*_*",
													"drop intedmValueMap_electronMVAValueMapProducer_*_*",
													#
													""
													"keep patPackedCandidates_*_*_*", # for intermediate PFCHSLeg jet constituents
													"keep recoGenParticles_flashggPrunedGenParticles_*_*", # this line, and preceding, drop unneded association object
													"keep recoVertexs_offlineSlimmedPrimaryVertices_*_*", # leave out floatedmValueMap_offlineSlimmedPrimaryVertices__PAT
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
													"keep *_slimmedAddPileupInfo_*_*", # Was huge in old MiniAod - hopefully better now
													"keep *GsfElectronCore*_*_*_*", # needed by at least one Tag
													"keep recoGenParticles_prunedGenParticles_*_*", # MiniAOD important status non-1
													"keep patPackedGenParticles_packedGenParticles_*_*", # MiniAOD status 1
													"keep *_selectedPatTrigger_*_*",                                                     
													"keep *_flashggSelected*_*_*",  

													"keep *_flashggDiLeptonDiJet_*_*", 
													"keep *_slimmedMuons_*_*", 
													"keep *_slimmedElectrons_*_*",
													"keep *_slimmedJets_*_*",
													"keep *_flashgg*Jet*_*_*",
													)
