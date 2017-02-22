#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"


using namespace flashgg;

MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate() {}

MultiLeptonMultiJetCandidate::~MultiLeptonMultiJetCandidate() {}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Electron_ptr electron1, Electron_ptr electron2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	type_ = kEEJJ;
	vertex_ = vertex;

    // FIXME: order in pt()

    //if (electron1->pt() < electron2->pt()) std::swap(electron1, electron2);
	ptrEle_.push_back(electron1);
	ptrEle_.push_back(electron2);
	ptrJet_.push_back(jet1);
	ptrJet_.push_back(jet2);

    this->setP4(ptrEle_[0]->p4() + ptrEle_[1]->p4() + ptrJet_[0]->p4() + ptrJet_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Muon_ptr muon1, Muon_ptr muon2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	type_ = kMMJJ;
	vertex_ = vertex;

    // FIXME: order in pt()

	ptrMuon_.push_back(muon1);
	ptrMuon_.push_back(muon2);
	ptrJet_.push_back(jet1);
	ptrJet_.push_back(jet2);

    this->setP4(ptrMuon_[0]->p4() + ptrMuon_[1]->p4() + ptrJet_[0]->p4() + ptrJet_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Electron_ptr electron1, Electron_ptr electron2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	type_ = kEETT;
	vertex_ = vertex;

    // FIXME: order in pt()

	ptrEle_.push_back(electron1);
	ptrEle_.push_back(electron2);
	ptrTrack_.push_back(track1);
	ptrTrack_.push_back(track2);

    this->setP4(ptrEle_[0]->p4() + ptrEle_[1]->p4() + ptrJet_[0]->p4() + ptrJet_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Muon_ptr muon1, Muon_ptr muon2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	type_ = kMMTT;
	vertex_ = vertex;

	ptrMuon_.push_back(muon1);
	ptrMuon_.push_back(muon2);
	ptrTrack_.push_back(track1);
	ptrTrack_.push_back(track2);

    this->setP4(ptrMuon_[0]->p4() + ptrMuon_[1]->p4() + ptrJet_[0]->p4() + ptrJet_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Electron_ptr electron, Muon_ptr muon, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	type_ = kEMJJ;
	vertex_ = vertex;

    // FIXME: order in pt()

	ptrEle_.push_back(electron);
	ptrMuon_.push_back(muon);
	ptrJet_.push_back(jet1);
	ptrJet_.push_back(jet2);

    this->setP4(ptrEle_[0]->p4() + ptrMuon_[0]->p4() + ptrJet_[0]->p4() + ptrJet_[1]->p4());
}


void MultiLeptonMultiJetCandidate::embedElectrons()
{
    while (!ptrEle_.empty()) {
        ele_.push_back(*(ptrEle_.back()));
        ptrEle_.pop_back();
    }
}


std::vector<Electron_t> & MultiLeptonMultiJetCandidate::embeddedElectrons()
{
    if (!ele_.size()) embedElectrons();
    return ele_;
}


//const Electron_t * MultiLeptonMultiJetCandidate::electron1() const
//{
//	if (type_ == kEEJJ || type_ == kEETT) return dynamic_cast<const Electron_t *>( daughter( 0 ) );  //electrons_[0] 
//	else return nullptr;
//}
//
//
//const Electron_t *MultiLeptonMultiJetCandidate::electron2() const
//{
//	if (type_ == kEEJJ || type_ == kEETT) return dynamic_cast<const Electron_t *>( daughter( 1 ) );
//	else return nullptr;
//}
//
//const Electron_t *MultiLeptonMultiJetCandidate::electron() const
//{
//	if (type_ == kEMJJ) return dynamic_cast<const Electron_t *>( daughter( 0 ) ); //electron_ 
//	else return nullptr;
//}
//
//
//const Muon_t *MultiLeptonMultiJetCandidate::muon1() const
//{
//	if (type_ == kMMJJ || type_ == kMMTT) return dynamic_cast<const Muon_t *>( daughter( 0 ) );
//	else return nullptr;	
//}
//
//
//const Muon_t *MultiLeptonMultiJetCandidate::muon2() const
//{
//	if (type_ == kMMJJ || type_ == kMMTT) return dynamic_cast<const Muon_t *>( daughter( 1 ) );
//	else return nullptr;
//}
//
//const Muon_t *MultiLeptonMultiJetCandidate::muon() const
//{
//	if (type_ == kEMJJ) return dynamic_cast<const Muon_t *>( daughter( 1 ) ); //muon_ 
//	else return nullptr;	
//}
//
//
const Electron_t * MultiLeptonMultiJetCandidate::leadingEle() const
{
	if (ele_.size()) {
        return &ele_[0];
    } else if (ptrEle_.size()) {
        return &*ptrEle_[0];
	}
    return nullptr;
}

const Electron_t * MultiLeptonMultiJetCandidate::subLeadingEle() const
{
	if (ele_.size() > 1) {
        return &ele_[1];
    } else if (ptrEle_.size() > 1) {
        return &*ptrEle_[1];
	}
    return nullptr;
}


const Muon_t *MultiLeptonMultiJetCandidate::leadingMuon() const
{
	if (muon_.size()) {
        return &muon_[0];
    } else if (ptrMuon_.size()) {
        return &*ptrMuon_[0];
	}
    return nullptr;
}

const Muon_t *MultiLeptonMultiJetCandidate::subLeadingMuon() const
{
	if (muon_.size() > 1) {
        return &muon_[1];
    } else if (ptrMuon_.size() > 1) {
        return &*ptrMuon_[1];
	}
    return nullptr;
}


const Jet_t *MultiLeptonMultiJetCandidate::leadingJet() const
{
	if (jet_.size()) {
        return &jet_[0];
    } else if (ptrJet_.size()) {
        return &*ptrJet_[0];
	}
    return nullptr;
}


const Jet_t *MultiLeptonMultiJetCandidate::subLeadingJet() const
{
	if (jet_.size() > 1) {
        return &jet_[1];
    } else if (ptrJet_.size() > 1) {
        return &*ptrJet_[1];
	}
    return nullptr;
}


const Track_t * MultiLeptonMultiJetCandidate::leadingTrack() const
{
	if (track_.size()) {
        return &track_[0];
    } else if (ptrTrack_.size()) {
        return &*ptrTrack_[0];
	}
    return nullptr;
}


const Track_t * MultiLeptonMultiJetCandidate::subLeadingTrack() const
{
	if (track_.size() > 1) {
        return &track_[1];
    } else if (ptrTrack_.size() > 1) {
        return &*ptrTrack_[1];
	}
    return nullptr;
}


const reco::Candidate * MultiLeptonMultiJetCandidate::leadingLepton() const
{
    // FIXME: can be done much better
    const reco::Candidate * le = leadingEle();
    const reco::Candidate * lm = leadingMuon();
    if (le && lm) return le->pt() > lm->pt() ? le : lm;
    else if (le)  return le;
    else if (lm)  return lm;
    return nullptr;
}


const reco::Candidate * MultiLeptonMultiJetCandidate::subLeadingLepton() const
{
    // FIXME: can be done much better
    // (e.g. order a vector of candidates)
    const reco::Candidate * le = leadingEle();
    const reco::Candidate * lm = leadingMuon();
    const reco::Candidate * sl = 0;
    // if leptons and muons, get the subleading one
    if (le && lm) le->pt() < lm->pt() ? sl = le : sl = lm;
    const reco::Candidate * se = subLeadingEle();
    const reco::Candidate * sm = subLeadingMuon();
    const reco::Candidate * ssl = 0;
    // if at least one subleading present for electrons or muons...
    if (se && sm) se->pt() > sm->pt() ? ssl = se : ssl = sm;
    else if (se)  ssl = se;
    else if (sm)  ssl = sm;
    // ...compare it to the previous subleading and chose
    if (sl && ssl && sl->pt() < ssl->pt()) sl = ssl;
    return sl;
}


float MultiLeptonMultiJetCandidate::sumPt() const
{	
    float ptsum = 0.;

    for (auto & e : ptrEle_)   ptsum += e->pt();
    for (auto & e : ptrMuon_)  ptsum += e->pt();
    for (auto & e : ptrJet_)   ptsum += e->pt();
    for (auto & e : ptrTrack_) ptsum += e->pt();

    // assume that if embedded, the pointers are empty
    // so that there is no double counting...
    for (auto & e : ele_)   ptsum += e.pt();
    for (auto & e : muon_)  ptsum += e.pt();
    for (auto & e : jet_)   ptsum += e.pt();
    for (auto & e : track_) ptsum += e.pt();

    return ptsum;
}  


//float MultiLeptonMultiJetCandidate::diLeptonInvMass() const
//{ 
//	TLorentzVector l1, l2;  
//	l1.SetPxPyPzE( 0., 0., 0., 0. );
//	l2.SetPxPyPzE( 0., 0., 0., 0. );
//	if (type_ == kEEJJ || type_ == kEETT) {
//		l1.SetPxPyPzE( electron1()->px(), electron1()->py(), electron1()->pz(), electron1()->energy() );
//		l2.SetPxPyPzE( electron2()->px(), electron2()->py(), electron2()->pz(), electron2()->energy() );
//	} else if (type_ == kMMJJ || type_ == kMMTT) {
//		l1.SetPxPyPzE( muon1()->px(), muon1()->py(), muon1()->pz(), muon1()->energy() );
//		l2.SetPxPyPzE( muon2()->px(), muon2()->py(), muon2()->pz(), muon2()->energy() );		
//	} else if (type_ == kEMJJ) {
//		l1.SetPxPyPzE( electron()->px(), electron()->py(), electron()->pz(), electron()->energy() );
//		l2.SetPxPyPzE( muon()->px(), muon()->py(), muon()->pz(), muon()->energy() );	
//	}
//	return (l1+l2).M();
//}


//float MultiLeptonMultiJetCandidate::leadingLeptonPt() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (leadingEle()->pt());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (leadingMuon()->pt());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return electron()->pt(); 
//		else return muon()->pt(); 
//	} 
//	else return 0.;
//}
//
//float MultiLeptonMultiJetCandidate::subLeadingLeptonPt() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (subLeadingEle()->pt());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (subLeadingMuon()->pt());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return muon()->pt(); 
//		else return electron()->pt(); 
//	} 	
//	else return 0.;
//}
//
//float MultiLeptonMultiJetCandidate::leadingLeptonEta() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (leadingEle()->eta());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (leadingMuon()->eta());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return electron()->eta(); 
//		else return muon()->eta(); 
//	} 	
//	else return 0.;
//}
//
//float MultiLeptonMultiJetCandidate::subLeadingLeptonEta() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (subLeadingEle()->eta());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (subLeadingMuon()->eta());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return muon()->eta(); 
//		else return electron()->eta(); 
//	} 	
//	else return 0.;
//}
//
//float MultiLeptonMultiJetCandidate::leadingLeptonPhi() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (leadingEle()->phi());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (leadingMuon()->phi());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return electron()->phi(); 
//		else return muon()->phi(); 
//	} 
//	else return 0.;
//}
//
//float MultiLeptonMultiJetCandidate::subLeadingLeptonPhi() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (subLeadingEle()->phi());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (subLeadingMuon()->phi());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return muon()->phi(); 
//		else return electron()->phi(); 
//	} 	
//	else return 0.;
//}
//
//float MultiLeptonMultiJetCandidate::leadingLeptonCharge() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (leadingEle()->charge());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (leadingMuon()->charge());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return electron()->charge(); 
//		else return muon()->charge(); 
//	} 
//	else return 0.;
//}
//
//float MultiLeptonMultiJetCandidate::subLeadingLeptonCharge() const {
//	if (type_ == kEEJJ || type_ == kEETT) return (subLeadingEle()->charge());
//	else if (type_ == kMMJJ || type_ == kMMTT) return (subLeadingMuon()->charge());
//	else if (type_ == kEMJJ) {
//		if( electron()->pt() > muon()->pt() ) return muon()->charge(); 
//		else return electron()->charge(); 
//	} 	
//	else return 0.;
//}


bool MultiLeptonMultiJetCandidate::operator <( const MultiLeptonMultiJetCandidate &b ) const
{
	return ( sumPt() < b.sumPt() );
}

bool MultiLeptonMultiJetCandidate::operator >( const MultiLeptonMultiJetCandidate &b ) const
{
	return ( sumPt() > b.sumPt() );
}



// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
