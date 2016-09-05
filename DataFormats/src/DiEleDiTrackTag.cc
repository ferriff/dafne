#include "dafne/DataFormats/interface/DiEleDiTrackTag.h"
#include <algorithm>

using namespace flashgg;

DiEleDiTrackTag::DiEleDiTrackTag() : DiLeptonDiJetTagBase::DiLeptonDiJetTagBase()
{}

DiEleDiTrackTag::~DiEleDiTrackTag()
{}

// N.B. Other attributes are set using methods in header file
DiEleDiTrackTag::DiEleDiTrackTag( edm::Ptr<DiLeptonDiJetCandidate> diEleDiTrack ) : DiLeptonDiJetTagBase::DiLeptonDiJetTagBase( diEleDiTrack ) {}  


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4