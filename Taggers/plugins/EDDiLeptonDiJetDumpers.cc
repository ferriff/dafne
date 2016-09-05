#include "FWCore/Framework/interface/MakerMacros.h"
#include "dafne/Taggers/interface/DiLeptonDiJetDumpers.h"
#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"

typedef edm::AnalyzerWrapper<flashgg::DiLeptonDiJetDumper> DiLeptonDiJetDumper;
typedef edm::AnalyzerWrapper<flashgg::CutBasedDiLeptonDiJetDumper> CutBasedDiLeptonDiJetDumper;


DEFINE_FWK_MODULE( DiLeptonDiJetDumper );
DEFINE_FWK_MODULE( CutBasedDiLeptonDiJetDumper );
