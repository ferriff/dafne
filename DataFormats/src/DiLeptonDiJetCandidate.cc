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

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate( Electron_ptr electron, Muon_ptr muon, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *electron );
	addDaughter( *muon );
	addDaughter( *jet1 );
	addDaughter( *jet2 );

	type_ = kEMJJ;
	electronPtr_ = electron; //ok?
	muonPtr_ = muon;
	jets_[0] = jet1;
	jets_[1] = jet2;

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

DiLeptonDiJetCandidate::DiLeptonDiJetCandidate(const Electron_t &electron, const Muon_t &muon, Jet_ptr jet1, Jet_ptr jet2, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( electron );
	addDaughter( muon );
	addDaughter( *jet1 );
	addDaughter( *jet2 );

	type_ = kEMJJ;
	electron_ = electron;
	muon_ = muon;
	jets_[0] = jet1;
	jets_[1] = jet2;

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

const Electron_t *DiLeptonDiJetCandidate::electron() const
{
	if (type_ == kEMJJ) return dynamic_cast<const Electron_t *>( daughter( 0 ) ); //electron_ 
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

const Muon_t *DiLeptonDiJetCandidate::muon() const
{
	if (type_ == kEMJJ) return dynamic_cast<const Muon_t *>( daughter( 1 ) ); //muon_ 
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


float DiLeptonDiJetCandidate::sumPt() const {	
	if (type_ == kEEJJ) return ( electron1()->pt() + electron2()->pt() + leadingJet()->pt() + subLeadingJet()->pt() );
	else if (type_ == kMMJJ) return ( muon1()->pt() + muon2()->pt() + leadingJet()->pt() + subLeadingJet()->pt() );
	else if (type_ == kEETT) return ( electron1()->pt() + electron2()->pt() + leadingTrack()->pt() + subLeadingTrack()->pt() );			
	else if (type_ == kMMTT) return ( muon1()->pt() + muon2()->pt() + leadingTrack()->pt() + subLeadingTrack()->pt() );
	else if (type_ == kEMJJ) return ( electron()->pt() + muon()->pt() + leadingJet()->pt() + subLeadingJet()->pt() );						
	else return 0.;
}  


float DiLeptonDiJetCandidate::invMass() const {  
	TLorentzVector l1, l2, j1, j2; 
	l1.SetPxPyPzE( 0., 0., 0., 0. );
	l2.SetPxPyPzE( 0., 0., 0., 0. );
	j1.SetPxPyPzE( 0., 0., 0., 0. );
	j2.SetPxPyPzE( 0., 0., 0., 0. );
	if (type_ == kEEJJ) {
		l1.SetPxPyPzE( electron1()->px(), electron1()->py(), electron1()->pz(), electron1()->energy() );
		l2.SetPxPyPzE( electron2()->px(), electron2()->py(), electron2()->pz(), electron2()->energy() );
		j1.SetPxPyPzE( leadingJet()->px(), leadingJet()->py(), leadingJet()->pz(), leadingJet()->energy() );
		j2.SetPxPyPzE( subLeadingJet()->px(), subLeadingJet()->py(), subLeadingJet()->pz(), subLeadingJet()->energy() );
	} else if (type_ == kMMJJ) {
		l1.SetPxPyPzE( muon1()->px(), muon1()->py(), muon1()->pz(), muon1()->energy() );
		l2.SetPxPyPzE( muon2()->px(), muon2()->py(), muon2()->pz(), muon2()->energy() );
		j1.SetPxPyPzE( leadingJet()->px(), leadingJet()->py(), leadingJet()->pz(), leadingJet()->energy() );
		j2.SetPxPyPzE( subLeadingJet()->px(), subLeadingJet()->py(), subLeadingJet()->pz(), subLeadingJet()->energy() );
	} else if (type_ == kEETT) {
		l1.SetPxPyPzE( electron1()->px(), electron1()->py(), electron1()->pz(), electron1()->energy() );
		l2.SetPxPyPzE( electron2()->px(), electron2()->py(), electron2()->pz(), electron2()->energy() );
		j1.SetPxPyPzE( leadingTrack()->px(), leadingTrack()->py(), leadingTrack()->pz(), leadingTrack()->energy() );
		j2.SetPxPyPzE( subLeadingTrack()->px(), subLeadingTrack()->py(), subLeadingTrack()->pz(), subLeadingTrack()->energy() );
	} else if (type_ == kMMTT) {
		l1.SetPxPyPzE( muon1()->px(), muon1()->py(), muon1()->pz(), muon1()->energy() );
		l2.SetPxPyPzE( muon2()->px(), muon2()->py(), muon2()->pz(), muon2()->energy() );
		j1.SetPxPyPzE( leadingTrack()->px(), leadingTrack()->py(), leadingTrack()->pz(), leadingTrack()->energy() );
		j2.SetPxPyPzE( subLeadingTrack()->px(), subLeadingTrack()->py(), subLeadingTrack()->pz(), subLeadingTrack()->energy() );
	} else if (type_ == kEMJJ) {
		l1.SetPxPyPzE( electron()->px(), electron()->py(), electron()->pz(), electron()->energy() );
		l2.SetPxPyPzE( muon()->px(), muon()->py(), muon()->pz(), muon()->energy() );
		j1.SetPxPyPzE( leadingJet()->px(), leadingJet()->py(), leadingJet()->pz(), leadingJet()->energy() );
		j2.SetPxPyPzE( subLeadingJet()->px(), subLeadingJet()->py(), subLeadingJet()->pz(), subLeadingJet()->energy() );
	}
	return (l1+l2+j1+j2).M();
}


float DiLeptonDiJetCandidate::diLeptonInvMass() const { 
	TLorentzVector l1, l2;  
	l1.SetPxPyPzE( 0., 0., 0., 0. );
	l2.SetPxPyPzE( 0., 0., 0., 0. );
	if (type_ == kEEJJ || type_ == kEETT) {
		l1.SetPxPyPzE( electron1()->px(), electron1()->py(), electron1()->pz(), electron1()->energy() );
		l2.SetPxPyPzE( electron2()->px(), electron2()->py(), electron2()->pz(), electron2()->energy() );
	} else if (type_ == kMMJJ || type_ == kMMTT) {
		l1.SetPxPyPzE( muon1()->px(), muon1()->py(), muon1()->pz(), muon1()->energy() );
		l2.SetPxPyPzE( muon2()->px(), muon2()->py(), muon2()->pz(), muon2()->energy() );		
	} else if (type_ == kEMJJ) {
		l1.SetPxPyPzE( electron()->px(), electron()->py(), electron()->pz(), electron()->energy() );
		l2.SetPxPyPzE( muon()->px(), muon()->py(), muon()->pz(), muon()->energy() );	
	}
	return (l1+l2).M();
}


float DiLeptonDiJetCandidate::leadingLeptonPt() const {
	if (type_ == kEEJJ || type_ == kEETT) return (leadingEle()->pt());
	else if (type_ == kMMJJ || type_ == kMMTT) return (leadingMuon()->pt());
	else if (type_ == kEMJJ) {
		if( electron()->pt() > muon()->pt() ) return electron()->pt(); 
		else return muon()->pt(); 
	} 
	else return 0.;
}

float DiLeptonDiJetCandidate::subLeadingLeptonPt() const {
	if (type_ == kEEJJ || type_ == kEETT) return (subLeadingEle()->pt());
	else if (type_ == kMMJJ || type_ == kMMTT) return (subLeadingMuon()->pt());
	else if (type_ == kEMJJ) {
		if( electron()->pt() > muon()->pt() ) return muon()->pt(); 
		else return electron()->pt(); 
	} 	
	else return 0.;
}

float DiLeptonDiJetCandidate::leadingLeptonEta() const {
	if (type_ == kEEJJ || type_ == kEETT) return (leadingEle()->eta());
	else if (type_ == kMMJJ || type_ == kMMTT) return (leadingMuon()->eta());
	else if (type_ == kEMJJ) {
		if( electron()->pt() > muon()->pt() ) return electron()->eta(); 
		else return muon()->eta(); 
	} 	
	else return 0.;
}

float DiLeptonDiJetCandidate::subLeadingLeptonEta() const {
	if (type_ == kEEJJ || type_ == kEETT) return (subLeadingEle()->eta());
	else if (type_ == kMMJJ || type_ == kMMTT) return (subLeadingMuon()->eta());
	else if (type_ == kEMJJ) {
		if( electron()->pt() > muon()->pt() ) return muon()->eta(); 
		else return electron()->eta(); 
	} 	
	else return 0.;
}

float DiLeptonDiJetCandidate::leadingLeptonPhi() const {
	if (type_ == kEEJJ || type_ == kEETT) return (leadingEle()->phi());
	else if (type_ == kMMJJ || type_ == kMMTT) return (leadingMuon()->phi());
	else if (type_ == kEMJJ) {
		if( electron()->pt() > muon()->pt() ) return electron()->phi(); 
		else return muon()->phi(); 
	} 
	else return 0.;
}

float DiLeptonDiJetCandidate::subLeadingLeptonPhi() const {
	if (type_ == kEEJJ || type_ == kEETT) return (subLeadingEle()->phi());
	else if (type_ == kMMJJ || type_ == kMMTT) return (subLeadingMuon()->phi());
	else if (type_ == kEMJJ) {
		if( electron()->pt() > muon()->pt() ) return muon()->phi(); 
		else return electron()->phi(); 
	} 	
	else return 0.;
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
