#ifndef flashgg_DiLeptonDiJetDumpers_h
#define flashgg_DiLeptonDiJetDumpers_h

#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "flashgg/Taggers/interface/CollectionDumper.h"
/// #include "PhysicsTools/UtilAlgos/interface/FWLiteAnalyzerWrapper.h"

namespace flashgg {

    typedef CollectionDumper<std::vector<DiLeptonDiJetCandidate> > DiLeptonDiJetDumper;
    typedef CollectionDumper<std::vector<DiLeptonDiJetCandidate>,
            DiLeptonDiJetCandidate,
            CutBasedClassifier<DiLeptonDiJetCandidate> > CutBasedDiLeptonDiJetDumper;

}

#endif // flashgg_DiLeptonDiJetDumpers_h


