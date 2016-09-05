#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/Merger.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"

typedef Merger<edm::OwnVector<flashgg::DiLeptonDiJetTagBase> > DiLeptoDiJetTagMerger;

DEFINE_FWK_MODULE( DiLeptoDiJetTagMerger );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
