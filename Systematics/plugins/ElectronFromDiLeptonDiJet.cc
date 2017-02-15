#include "dafne/Systematics/interface/ElectronFromDiLeptonDiJetBase.h"

namespace flashgg {

    typedef ElectronFromDiLeptonDiJetBase<int> ElectronFromDiLeptonDiJet;

}

DEFINE_EDM_PLUGIN( FlashggSystematicDiLeptonDiJetMethodsFactory,
                   flashgg::ElectronFromDiLeptonDiJet,
                   "FlashggElectronFromDiLeptonDiJet" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4