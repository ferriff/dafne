#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

using namespace flashgg;

DiLeptonDiJetTagBase::DiLeptonDiJetTagBase()
{
    category_number_ = -1;
    // isGold_ = -1;
}

DiLeptonDiJetTagBase::DiLeptonDiJetTagBase( edm::Ptr<DiLeptonDiJetCandidate> diLeptonDiJet )
{
    category_number_ = -1;
    diLeptonDiJet_ = diLeptonDiJet;
    // isGold_ = -1;
}

DiLeptonDiJetTagBase::~DiLeptonDiJetTagBase()
{
}


bool DiLeptonDiJetTagBase::operator <( const DiLeptonDiJetTagBase &b ) const
{
    // For choosing which of two tags OF THE SAME TYPE are of higher priority
    // Comparison of different tags not currently supported - is it ever needed?
    // Overloading may be appropriate if different tags have different priorities

    if( categoryNumber() == b.categoryNumber() ) {
        return ( sumPt() < b.sumPt() );
    } else {
        return ( categoryNumber() < b.categoryNumber() );
    }
}

// void DiLeptonDiJetTagBase::setIsGold( int runNumber ) {
//     // Below is the subtraction between
//     // https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt
//     // and
//     // https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt
//     //
//     // They were last changed on 18 December and still the most recent (as of 28 January)
//     // I have confirmed that none of these runs are in the Gold at all, so checking runNumber without lumiSection suffices
//     // Of course the result may be wrong for events not even in the Silver JSON
//     // I am hard coding this to make dumping and comparing events faster, because this issue now needs intensive study
//     //
//     isGold_ = 1;
//     if ( runNumber == 256729 ) { isGold_ = 0; }
//     if ( runNumber == 256734 ) { isGold_ = 0; }
//     if ( (runNumber >= 257394) && (runNumber <= 257397) ) { isGold_ = 0; }
//     if ( (runNumber >= 257399) && (runNumber <=257400)) { isGold_ = 0; }
//     if ( runNumber == 257487 ) { isGold_ = 0; }
//     if ( runNumber == 257490 ) { isGold_ = 0; }
//     if ( (runNumber >= 257822) && (runNumber <=257823)) { isGold_ = 0; }
//     if ( runNumber == 258443 ) { isGold_ = 0; }
// }

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

