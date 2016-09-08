#ifndef flashgg_DiLeptonDiJetTagsDumpers_h
#define flashgg_DiLeptonDiJetTagsDumpers_h

#include "flashgg/Taggers/interface/CollectionDumper.h"

#include "dafne/DataFormats/interface/DiEleDiJetTag.h"
#include "dafne/DataFormats/interface/DiEleDiTrackTag.h"
#include "dafne/DataFormats/interface/DiMuDiJetTag.h"
#include "dafne/DataFormats/interface/DiMuDiTrackTag.h"



namespace flashgg {
    typedef CollectionDumper<std::vector<DiEleDiJetTag>,
            DiEleDiJetTag,
            CutBasedClassifier<DiEleDiJetTag> > CutBasedDiEleDiJetTagDumper;

    typedef CollectionDumper<std::vector<DiEleDiTrackTag>,
            DiEleDiTrackTag,
            CutBasedClassifier<DiEleDiTrackTag> > CutBasedDiEleDiTrackTagDumper;

    typedef CollectionDumper<std::vector<DiMuDiJetTag>,
            DiMuDiJetTag,
            CutBasedClassifier<DiMuDiJetTag> > CutBasedDiMuDiJetTagDumper;

    typedef CollectionDumper<std::vector<DiMuDiTrackTag>,
            DiMuDiTrackTag,
            CutBasedClassifier<DiMuDiTrackTag> > CutBasedDiMuDiTrackTagDumper;
}

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
