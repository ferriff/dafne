#include "dafne/DataFormats/interface/DiEleDiJetTag.h"
#include <algorithm>

using namespace flashgg;

DiEleDiJetTag::DiEleDiJetTag() : DiLeptonDiJetTagBase::DiLeptonDiJetTagBase()
{}

DiEleDiJetTag::~DiEleDiJetTag()
{}

// N.B. Other attributes are set using methods in header file
DiEleDiJetTag::DiEleDiJetTag( edm::Ptr<DiLeptonDiJetCandidate> diEleDiJet ) : DiLeptonDiJetTagBase::DiLeptonDiJetTagBase( diEleDiJet ) {}  


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4