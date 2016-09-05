#ifndef flashgg_MicroAOD__FLASHggObjectSelector_h
#define flashgg_MicroAOD__FLASHggObjectSelector_h

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <vector>


typedef pat::PackedCandidate Track_t;


namespace flashgg {


    typedef SingleObjectSelector <        
    std::vector<Track_t>,                 
        StringCutObjectSelector<Track_t>  
        > FLASHggTrackSelector;

}

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
