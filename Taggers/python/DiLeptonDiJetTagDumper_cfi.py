import FWCore.ParameterSet.Config as cms

from dafne.Taggers.DiLeptonDiJetTagsDumpConfig_cff import DiLeptonDiJetTagsDumpConfig

DiLeptonDiJetTagsDumper = cms.EDAnalyzer('DiLeptonDiJetTagDumper',
                            **DiLeptonDiJetTagsDumpConfig.parameters_()
                            )