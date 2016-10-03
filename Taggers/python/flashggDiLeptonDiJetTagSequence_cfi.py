import FWCore.ParameterSet.Config as cms
# from flashgg.Taggers.flashggDiPhotonMVA_cfi import flashggDiPhotonMVA
# from flashgg.Taggers.flashggVBFMVA_cff import flashggVBFMVA,flashggVBFDiPhoDiJetMVA
from flashgg.Taggers.flashggTags_cff import *
# from flashgg.Taggers.flashggPreselectedDiPhotons_cfi import flashggPreselectedDiPhotons
from dafne.Taggers.flashggDiLeptonDiJetTagSorter_cfi import flashggDiLeptonDiJetTagSorter
# from flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi import flashggUpdatedIdMVADiPhotons

flashggDiLeptonDiJetTagSequence = cms.Sequence(#flashggUpdatedIdMVADiPhotons
                              #     * flashggPreselectedDiPhotons
		                        		  # * flashggDiPhotonMVA
                              #     * flashggUnpackedJets
                              #     * flashggVBFMVA
                              #     * flashggVBFDiPhoDiJetMVA
                                  #* 
                                    ( flashggDiEleDiJetTag
                                     + flashggDiMuDiJetTag
                                     + flashggDiEleDiTrackTag
                                     + flashggDiMuDiTrackTag                                      
				                            	)
                                    * flashggDiLeptonDiJetTagSorter
                                  )