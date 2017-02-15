#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "dafne/Systematics/interface/BaseSystMethod_DiLeptonDiJet.h"
// #include "flashgg/Systematics/interface/ObjectSystematicProducer.h"
#include "dafne/Systematics/interface/ObjectSystematicProducer_DiLeptonDiJet.h"

namespace flashgg {

    typedef ObjectSystematicProducer_DiLeptonDiJet<DiLeptonDiJetCandidate, int, std::vector> DiLeptonDiJetSystematicProducer;

}

typedef flashgg::DiLeptonDiJetSystematicProducer FlashggDiLeptonDiJetSystematicProducer;
DEFINE_FWK_MODULE( FlashggDiLeptonDiJetSystematicProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4