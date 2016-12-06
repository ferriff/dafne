#include "FWCore/Framework/interface/MakerMacros.h"

#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"
#include "dafne/Validation/interface/miniTreeMaker.h"

typedef edm::AnalyzerWrapper<miniTreeMaker> EDminiTreeMaker;
DEFINE_FWK_MODULE(EDminiTreeMaker);