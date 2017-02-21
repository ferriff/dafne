import FWCore.ParameterSet.Config as cms

microAODmassiveNuOutputCommand = cms.untracked.vstring("drop *",       
													"keep *_flashgg*_*_*",
													"drop *_flashggVertexMap*_*_*", 
													## this part drop all the tools used to build puppi jets
													"drop *_flashggPuppi*_*_*",
													"drop *_flashggPhotons_*_*", # Only keep the copies with random numbers added
													#
													""
													"drop patPackedCandidates_*_*_*", # for intermediate PFCHSLeg jet constituents
													"drop *_flashggPrunedGenParticles_*_*",   
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
													"keep *_slimmedAddPileupInfo_*_*", 
													"keep *GsfElectronCore*_*_*_*", # needed by at least one Tag

													"keep *_flashggSelected*_*_*",
													# Drop intermediate collections in favor of selected/final collections
													"drop *_flashgg*Jet*_*_*",
													"drop *_flashggMuons_*_*",
													"drop *_flashggElectrons_*_*",

													"keep *_flashggFinalJets_*_*",
													"keep *_flashggFinalPuppiJets_*_*",
													"drop floatedmValueMap_electronMVAValueMapProducer_*_*",
													"drop intedmValueMap_electronMVAValueMapProducer_*_*",
													"keep *_selectedPatTrigger_*_*",
													"keep *_flashggDiLeptonDiJet_*_*",
													"keep *_flashggTriLeptons_*_*",
													"keep patPackedCandidates_packedPFCandidates_*_*"
													 )

