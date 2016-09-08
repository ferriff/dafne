#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"

#include "dafne/Taggers/interface/DiLeptonDiJetTagsDumpers.h"


typedef edm::AnalyzerWrapper<flashgg::CutBasedDiEleDiJetTagDumper> CutBasedDiEleDiJetTagDumper;
typedef edm::AnalyzerWrapper<flashgg::CutBasedDiEleDiTrackTagDumper> CutBasedDiEleDiTrackTagDumper;
typedef edm::AnalyzerWrapper<flashgg::CutBasedDiMuDiJetTagDumper> CutBasedDiMuDiJetTagDumper;
typedef edm::AnalyzerWrapper<flashgg::CutBasedDiMuDiTrackTagDumper> CutBasedDiMuDiTrackTagDumper;


DEFINE_FWK_MODULE( CutBasedDiEleDiJetTagDumper );
DEFINE_FWK_MODULE( CutBasedDiEleDiTrackTagDumper );
DEFINE_FWK_MODULE( CutBasedDiMuDiJetTagDumper );
DEFINE_FWK_MODULE( CutBasedDiMuDiTrackTagDumper );


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4