#include "dafne/Systematics/interface/ElectronFromDiLeptonDiJetBase.h"

namespace flashgg {

    typedef ElectronFromDiLeptonDiJetBase<std::pair<int, int> > ElectronFromDiLeptonDiJet2D;

}

DEFINE_EDM_PLUGIN( FlashggSystematicDiLeptonDiJetMethodsFactory2D,
                   flashgg::ElectronFromDiLeptonDiJet2D,
                   "FlashggElectronFromDiLeptonDiJet2D" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4