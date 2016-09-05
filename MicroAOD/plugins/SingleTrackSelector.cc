#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelectorStream.h"
#include "CommonTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


typedef pat::PackedCandidate Track_t;


typedef SingleObjectSelector <
edm::View<Track_t>,
    StringCutObjectSelector<Track_t, true>,
    std::vector<Track_t>
    > TrackSelector;


DEFINE_FWK_MODULE( TrackSelector );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4