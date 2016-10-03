import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections #???

# bDiscriminator74X = cms.vdouble(0.605,0.890)
# bDiscriminator76X = cms.vdouble(0.460,0.800)

# flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
#                                      JetsTag = cms.InputTag("flashggFinalJets"),
#                                      NCollections = cms.uint32(maxJetCollections)
#                                      )

# UnpackedJetCollectionVInputTag = cms.VInputTag()
# for i in range(0,maxJetCollections):
#     UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))


flashggDiEleDiJetTag = cms.EDProducer("FlashggDiEleDiJetTagProducer",
                                      DiLeptonDiJetTag = cms.InputTag('flashggDiLeptonDiJets'), 
                                      inputTagJets= UnpackedJetCollectionVInputTag, #???
                                      # JetTag=cms.InputTag('slimmedJets'),                                     
                                      ElectronTag=cms.InputTag('flashggElectrons'),  
                                      VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      GenParticleTag=cms.InputTag('flashggPrunedGenParticles'),  #ce le ho?
                                      # SystLabel=cms.string(""),
                                      jetPtThreshold = cms.double(25.),
                                      jetEtaThreshold = cms.double(2.4),                          
                                      electronPtThreshold = cms.double(20),
                                      electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                      useStdLeptonID = cms.bool(True),
                                      useElectronMVARecipe = cms.bool(False),
                                      useElectronLooseID = cms.bool(True), 
                                      # electronIsoThreshold = cms.double(0.15),
  		                                elMiniIsoEBThreshold = cms.double(0.045),
                                      elMiniIsoEEThreshold = cms.double(0.08)#, 
                                      # electronNumOfHitsThreshold = cms.double(1),
                                      # TransverseImpactParam = cms.double(0.02),
                                      # LongitudinalImpactParam = cms.double(0.2),                                
                                      )


flashggDiMuDiJetTag = cms.EDProducer("FlashggDiMuDiJetTagProducer",
                                      DiLeptonDiJetTag = cms.InputTag('flashggDiLeptonDiJets'), 
                                      inputTagJets= UnpackedJetCollectionVInputTag, #???
                                      # JetTag=cms.InputTag('slimmedJets'),
                                      MuonTag=cms.InputTag('flashggMuons'),  
                                      VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      GenParticleTag=cms.InputTag('flashggPrunedGenParticles'), #ce le ho?
                                      # SystLabel=cms.string(""),
                                      jetPtThreshold = cms.double(25.), 
                                      jetEtaThreshold= cms.double(2.4),
                                      muonPtThreshold = cms.double(20),
                                      muonEtaThreshold = cms.double(2.4),
                                      useStdLeptonID = cms.bool(True),
                                      muPFIsoSumRelThreshold = cms.double(0.25), 
				                              muMiniIsoSumRelThreshold = cms.double(0.05)
                                      )


flashggDiEleDiTrackTag = cms.EDProducer("FlashggDiEleDiTrackTagProducer",
                                      DiLeptonDiJetTag = cms.InputTag('flashggDiLeptonDiJets'), 
                                      ElectronTag=cms.InputTag('flashggElectrons'),  
                                      TrackTag = cms.InputTag('packedPFCandidates'),
                                      VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      GenParticleTag=cms.InputTag('flashggPrunedGenParticles'),  #ce le ho?
                                      # SystLabel=cms.string(""),
                                      trackPtThreshold = cms.double(25.),
                                      trackEtaThreshold = cms.double(2.4),                          
                                      electronPtThreshold = cms.double(20),
                                      electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                      useStdLeptonID = cms.bool(True),
                                      useElectronMVARecipe = cms.bool(False),
                                      useElectronLooseID = cms.bool(True), 
                                      # electronIsoThreshold = cms.double(0.15),
                                      elMiniIsoEBThreshold = cms.double(0.045),
                                      elMiniIsoEEThreshold = cms.double(0.08)#, 
                                      # electronNumOfHitsThreshold = cms.double(1),
                                      # TransverseImpactParam = cms.double(0.02),
                                      # LongitudinalImpactParam = cms.double(0.2),                                
                                      )


flashggDiMuDiTrackTag = cms.EDProducer("FlashggDiMuDiTrackTagProducer",
                                      DiLeptonDiJetTag = cms.InputTag('flashggDiLeptonDiJets'), 
                                      MuonTag=cms.InputTag('flashggMuons'), 
                                      TrackTag = cms.InputTag('packedPFCandidates'), 
                                      VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      GenParticleTag=cms.InputTag('flashggPrunedGenParticles'), #ce le ho?
                                      # SystLabel=cms.string(""),
                                      trackPtThreshold = cms.double(25.),
                                      trackEtaThreshold = cms.double(2.4),  
                                      muonPtThreshold = cms.double(20),
                                      muonEtaThreshold = cms.double(2.4),
                                      useStdLeptonID = cms.bool(True),
                                      muPFIsoSumRelThreshold = cms.double(0.25), 
                                      muMiniIsoSumRelThreshold = cms.double(0.05)
                                      )


