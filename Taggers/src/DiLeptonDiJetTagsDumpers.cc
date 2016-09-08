#include "flashgg/Taggers/interface/PluggableAnalyzer.h"
#include "dafne/Taggers/interface/DiLeptonDiJetTagsDumpers.h"

namespace flashgg {
    namespace fwlite {

        PLUGGABLE_ANALYZER( CutBasedDiEleDiJetTagDumper );
        PLUGGABLE_ANALYZER( CutBasedDiEleDiTrackTagDumper );
        PLUGGABLE_ANALYZER( CutBasedDiMuDiJetTagDumper );
        PLUGGABLE_ANALYZER( CutBasedDiMuDiTrackTagDumper );    

    }
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
