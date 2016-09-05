import FWCore.ParameterSet.Config as cms

from dafne.Taggers.DiLeptonDiJetDumpConfig_cff import DiLeptonDiJetDumpConfig

DiLeptonDiJetDumper = cms.EDAnalyzer('CutBasedDiLeptonDiJetDumper',
                                **DiLeptonDiJetDumpConfig.parameters_()
                                )

