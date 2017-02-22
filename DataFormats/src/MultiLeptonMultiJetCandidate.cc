#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"


using namespace flashgg;

MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate() {}

MultiLeptonMultiJetCandidate::~MultiLeptonMultiJetCandidate() {}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Electron_ptr electron1, Electron_ptr electron2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	type_ = kEEJJ;
	vertex_ = vertex;

    if (electron1->pt() < electron2->pt()) std::swap(electron1, electron2);
	ptrEle_.push_back(electron1);
	ptrEle_.push_back(electron2);
    if (jet1->pt() < jet2->pt()) std::swap(jet1, jet2);
	ptrJet_.push_back(jet1);
	ptrJet_.push_back(jet2);

    this->setP4(ptrEle_[0]->p4() + ptrEle_[1]->p4() + ptrJet_[0]->p4() + ptrJet_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Muon_ptr muon1, Muon_ptr muon2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	type_ = kMMJJ;
	vertex_ = vertex;

    if (muon1->pt() < muon2->pt()) std::swap(muon1, muon2);
	ptrMuon_.push_back(muon1);
	ptrMuon_.push_back(muon2);
    if (jet1->pt() < jet2->pt()) std::swap(jet1, jet2);
	ptrJet_.push_back(jet1);
	ptrJet_.push_back(jet2);

    this->setP4(ptrMuon_[0]->p4() + ptrMuon_[1]->p4() + ptrJet_[0]->p4() + ptrJet_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Electron_ptr electron1, Electron_ptr electron2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	type_ = kEETT;
	vertex_ = vertex;

    if (electron1->pt() < electron2->pt()) std::swap(electron1, electron2);
	ptrEle_.push_back(electron1);
	ptrEle_.push_back(electron2);
    if (track1->pt() < track2->pt()) std::swap(track1, track2);
	ptrTrack_.push_back(track1);
	ptrTrack_.push_back(track2);

    this->setP4(ptrEle_[0]->p4() + ptrEle_[1]->p4() + ptrTrack_[0]->p4() + ptrTrack_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Muon_ptr muon1, Muon_ptr muon2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	type_ = kMMTT;
	vertex_ = vertex;

    if (muon1->pt() < muon2->pt()) std::swap(muon1, muon2);
	ptrMuon_.push_back(muon1);
	ptrMuon_.push_back(muon2);
    if (track1->pt() < track2->pt()) std::swap(track1, track2);
	ptrTrack_.push_back(track1);
	ptrTrack_.push_back(track2);

    this->setP4(ptrMuon_[0]->p4() + ptrMuon_[1]->p4() + ptrTrack_[0]->p4() + ptrTrack_[1]->p4());
}


MultiLeptonMultiJetCandidate::MultiLeptonMultiJetCandidate( Electron_ptr electron, Muon_ptr muon, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	type_ = kEMJJ;
	vertex_ = vertex;

	ptrEle_.push_back(electron);
	ptrMuon_.push_back(muon);
    if (jet1->pt() < jet2->pt()) std::swap(jet1, jet2);
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

    // assume that, if candidate are embedded,
    // the pointers are empty so that there is
    // no double counting...
    for (auto & e : ele_)   ptsum += e.pt();
    for (auto & e : muon_)  ptsum += e.pt();
    for (auto & e : jet_)   ptsum += e.pt();
    for (auto & e : track_) ptsum += e.pt();

    return ptsum;
}  


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
