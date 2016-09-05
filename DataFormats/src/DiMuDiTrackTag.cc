#include "dafne/DataFormats/interface/DiMuDiTrackTag.h"
#include <algorithm>

using namespace flashgg;

DiMuDiTrackTag::DiMuDiTrackTag() : DiLeptonDiJetTagBase::DiLeptonDiJetTagBase()
{}

DiMuDiTrackTag::~DiMuDiTrackTag()
{}

// N.B. Other attributes are set using methods in header file
DiMuDiTrackTag::DiMuDiTrackTag( edm::Ptr<DiLeptonDiJetCandidate> diMuDiTrack ) : DiLeptonDiJetTagBase::DiLeptonDiJetTagBase( diMuDiTrack ) {}  


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4