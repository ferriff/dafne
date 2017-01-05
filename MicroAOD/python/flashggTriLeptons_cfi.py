import FWCore.ParameterSet.Config as cms

flashggTriLeptons = cms.EDProducer('FlashggTriLeptonsProducer',     
									ElectronTag=cms.InputTag('flashggElectrons'),  #'slimmedElectrons'), 
									MuonTag=cms.InputTag('flashggMuons'),          #'slimmedMuons'),
									VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'), 
									##Parameters  
                                    minElectronPt=cms.double(5.),
                                    maxElectronEta=cms.double(2.4),
									minMuonPt=cms.double(5.),
									maxMuonEta=cms.double(2.4)
                                  )