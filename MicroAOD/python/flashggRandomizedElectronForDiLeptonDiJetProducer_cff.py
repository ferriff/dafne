
import FWCore.ParameterSet.Config as cms

# add the following lines to the module used to configure the RandomNumberGeneratorService
# process.RandomNumberGeneratorService.flashggRandomizedElectronForDiLeptonDiJets = cms.PSet(
#           initialSeed = cms.untracked.uint32(281765313)
#         # engineName = cms.untracked.string('TRandom3') # optional, default to HepJamesRandom if absent
#         )


flashggRandomizedElectronForDiLeptonDiJets = cms.EDProducer("FlashggRandomizedElectronForDiLeptonDiJetProducer",
                                    src = cms.InputTag("flashggDiLeptonDiJet"),
#                                    outputCollectionName = cms.string("flashggFinalRandomizedDiLeptonDiJets"),
                                    # labels of various gaussian random numbers with mean=0, sigma=1
                                    # to be associated with the photon object
                                    labels = cms.vstring("smearE")
                                    )
