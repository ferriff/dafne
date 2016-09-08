#ifndef flashgg_DiLeptonDiJetTagCandidateDumper_h
#define flashgg_DiLeptonDiJetTagCandidateDumper_h

#include "dafne/DataFormats/interface/DiLeptonDiJetTagCandidate.h"
#include "flashgg/Taggers/interface/CollectionDumper.h"

namespace flashgg {
    typedef CollectionDumper<std::vector<DiLeptonDiJetTagCandidate>, DiLeptonDiJetTagCandidate,
                             CutBasedClassifier<DiLeptonDiJetTagCandidate> > CutBasedDiLeptonDiJetTagCandidateDumper;
}

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4