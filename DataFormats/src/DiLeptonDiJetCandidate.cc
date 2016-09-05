#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"


using namespace flashgg;

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate() {}

DiLeptonDiJetCandidate::~DiLeptonDiJetCandidate() {}


DiLeptonDiJetCandidate::DiLeptonDiJetCandidate( Electron_ptr electron1, Electron_ptr electron2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *electron1 );
	addDaughter( *electron2 );
	addDaughter( *jet1 );
	addDaughter( *jet2 );

	type_ = kEEJJ;
	electrons_[0] = electron1;
	electrons_[1] = electron2;
	jets_[0] = jet1;
	jets_[1] = jet2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate( Muon_ptr muon1, Muon_ptr muon2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *muon1 );
	addDaughter( *muon2 );
	addDaughter( *jet1 );
	addDaughter( *jet2 );

	type_ = kMMJJ;
	muons_[0] = muon1;
	muons_[1] = muon2;
	jets_[0] = jet1;
	jets_[1] = jet2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate( Electron_ptr electron1, Electron_ptr electron2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *electron1 );
	addDaughter( *electron2 );
	addDaughter( *track1 );
	addDaughter( *track2 );

	type_ = kEETT;
	electrons_[0] = electron1;
	electrons_[1] = electron2;
	tracks_[0] = track1;
	tracks_[1] = track2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate( Muon_ptr muon1, Muon_ptr muon2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *muon1 );
	addDaughter( *muon2 );
	addDaughter( *track1 );
	addDaughter( *track2 );

	type_ = kMMTT;
	muons_[0] = muon1;
	muons_[1] = muon2;
	tracks_[0] = track1;
	tracks_[1] = track2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}


DiLeptonDiJetCandidate::DiLeptonDiJetCandidate(const Electron_t &electron1, const Electron_t &electron2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( electron1 );
	addDaughter( electron2 );
	addDaughter( *jet1 );
	addDaughter( *jet2 );

	type_ = kEEJJ;
	electron1_ = electron1;
	electron2_ = electron2;
	jets_[0] = jet1;
	jets_[1] = jet2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate(const Muon_t &muon1, const Muon_t &muon2, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( muon1 );
	addDaughter( muon2 );
	addDaughter( *jet1 );
	addDaughter( *jet2 );

	type_ = kMMJJ;
	muon1_ = muon1;
	muon2_ = muon2;
	jets_[0] = jet1;
	jets_[1] = jet2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate(const Electron_t &electron1, const Electron_t &electron2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( electron1 );
	addDaughter( electron2 );
	addDaughter( *track1 );
	addDaughter( *track2 );

	type_ = kEETT;
	electron1_ = electron1;
	electron2_ = electron2;
	tracks_[0] = track1;
	tracks_[1] = track2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate(const Muon_t &muon1, const Muon_t &muon2, Track_ptr track1, Track_ptr track2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( muon1 );
	addDaughter( muon2 );
	addDaughter( *track1 );
	addDaughter( *track2 );

	type_ = kMMTT;
	muon1_ = muon1;
	muon2_ = muon2;
	tracks_[0] = track1;
	tracks_[1] = track2;

	AddFourMomenta addP4;  
	addP4.set( *this );
}



const Electron_t *DiLeptonDiJetCandidate::electron1() const
{
	if (type_ == kEEJJ || type_ == kEETT) return dynamic_cast<const Electron_t *>( daughter( 0 ) );  //electrons_[0] 
	else return nullptr;
}


const Electron_t *DiLeptonDiJetCandidate::electron2() const
{
	if (type_ == kEEJJ || type_ == kEETT) return dynamic_cast<const Electron_t *>( daughter( 1 ) );
	else return nullptr;
}


const Muon_t *DiLeptonDiJetCandidate::muon1() const
{
	if (type_ == kMMJJ || type_ == kMMTT) return dynamic_cast<const Muon_t *>( daughter( 0 ) );
	else return nullptr;	
}


const Muon_t *DiLeptonDiJetCandidate::muon2() const
{
	if (type_ == kMMJJ || type_ == kMMTT) return dynamic_cast<const Muon_t *>( daughter( 1 ) );
	else return nullptr;
}


const Electron_t *DiLeptonDiJetCandidate::leadingEle() const
{
	if (type_ == kEEJJ || type_ == kEETT) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Electron_t *>( daughter( 0 ) );
		} else {
			return dynamic_cast<const Electron_t *>( daughter( 1 ) );
		}
	} else return nullptr;
}

const Electron_t *DiLeptonDiJetCandidate::subLeadingEle() const
{
	if (type_ == kEEJJ || type_ == kEETT) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Electron_t *>( daughter( 1 ) );
		} else {
			return dynamic_cast<const Electron_t *>( daughter( 0 ) );
		}
	} else return nullptr;
}


const Muon_t *DiLeptonDiJetCandidate::leadingMuon() const
{
	if (type_ == kMMJJ || type_ == kMMTT) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Muon_t *>( daughter( 0 ) );
		} else {
			return dynamic_cast<const Muon_t *>( daughter( 1 ) );
		}
	} else return nullptr;
}

const Muon_t *DiLeptonDiJetCandidate::subLeadingMuon() const
{
	if (type_ == kMMJJ || type_ == kMMTT) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Muon_t *>( daughter( 1 ) );
		} else {
			return dynamic_cast<const Muon_t *>( daughter( 0 ) );
		}
	} else return nullptr;
}


const Jet_t *DiLeptonDiJetCandidate::leadingJet() const
{
	if (type_ == kEEJJ || type_ == kMMJJ) {
		if( daughter( 2 )->pt() > daughter( 3 )->pt() ) {
			return dynamic_cast<const Jet_t *>( daughter( 2 ) );
		} else {
			return dynamic_cast<const Jet_t *>( daughter( 3 ) );
		}
	} else return nullptr;
}


const Jet_t *DiLeptonDiJetCandidate::subLeadingJet() const
{
	if (type_ == kEEJJ || type_ == kMMJJ) {
		if( daughter( 2 )->pt() > daughter( 3 )->pt() ) {
			return dynamic_cast<const Jet_t *>( daughter( 3 ) );
		} else {
			return dynamic_cast<const Jet_t *>( daughter( 2 ) );
		}
	} else return nullptr;
}


const Track_t *DiLeptonDiJetCandidate::leadingTrack() const
{
	if (type_ == kEETT || type_ == kMMTT) {	
    	if( daughter( 2 )->pt() > daughter( 3 )->pt() ) {          //if( tracks_[0]->pt() > tracks_[1]->pt() ) {
        	return dynamic_cast<const Track_t *>( daughter( 2 ) ); //return tracks_[0];
	    } else {
    	    return dynamic_cast<const Track_t *>( daughter( 3 ) ); //return tracks_[1];
    	}
	} else return nullptr;
}


const Track_t *DiLeptonDiJetCandidate::subLeadingTrack() const
{
	if (type_ == kEETT || type_ == kMMTT) {	
	    if( daughter( 2 )->pt() > daughter( 3 )->pt() ) {          //if( tracks_[0]->pt() > tracks_[1]->pt() ) {
    	    return dynamic_cast<const Track_t *>( daughter( 3 ) ); //return tracks_[1];
	    } else {
    	    return dynamic_cast<const Track_t *>( daughter( 2 ) ); //return tracks_[0];
    	}
    } else return nullptr;
}


bool DiLeptonDiJetCandidate::operator <( const DiLeptonDiJetCandidate &b ) const
{
	return ( sumPt() < b.sumPt() );
}

bool DiLeptonDiJetCandidate::operator >( const DiLeptonDiJetCandidate &b ) const
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
