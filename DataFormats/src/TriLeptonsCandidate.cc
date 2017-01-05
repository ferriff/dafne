#include "dafne/DataFormats/interface/TriLeptonsCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"


using namespace flashgg;

TriLeptonsCandidate::TriLeptonsCandidate() {}

TriLeptonsCandidate::~TriLeptonsCandidate() {}


TriLeptonsCandidate::TriLeptonsCandidate( Electron_ptr electron1, Electron_ptr electron2, Electron_ptr electron3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *electron1 );
	addDaughter( *electron2 );
	addDaughter( *electron3 );

	type_ = kEEE;
	electrons_[0] = electron1;
	electrons_[1] = electron2;
	electrons_[2] = electron3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

TriLeptonsCandidate::TriLeptonsCandidate( Electron_ptr electron1, Electron_ptr electron2, Muon_ptr muon3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *electron1 );
	addDaughter( *electron2 );
	addDaughter( *muon3 );

	type_ = kEEM;
	electrons_[0] = electron1;
	electrons_[1] = electron2;
	muons_[2] = muon3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

TriLeptonsCandidate::TriLeptonsCandidate( Muon_ptr muon1, Muon_ptr muon2, Muon_ptr muon3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *muon1 );
	addDaughter( *muon2 );
	addDaughter( *muon3 );

	type_ = kMMM;
	muons_[0] = muon1;
	muons_[1] = muon2;
	muons_[2] = muon3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

TriLeptonsCandidate::TriLeptonsCandidate( Muon_ptr muon1, Muon_ptr muon2, Electron_ptr electron3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( *muon1 );
	addDaughter( *muon2 );
	addDaughter( *electron3 );

	type_ = kMME;
	muons_[0] = muon1;
	muons_[1] = muon2;
	electrons_[2] = electron3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}



TriLeptonsCandidate::TriLeptonsCandidate(const Electron_t &electron1, const Electron_t &electron2, const Electron_t &electron3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( electron1 );
	addDaughter( electron2 );
	addDaughter( electron3 );

	type_ = kEEE;
	electron1_ = electron1;
	electron2_ = electron2;
	electron3_ = electron3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

TriLeptonsCandidate::TriLeptonsCandidate(const Electron_t &electron1, const Electron_t &electron2, const Muon_t &muon3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( electron1 );
	addDaughter( electron2 );
	addDaughter( muon3 );

	type_ = kEEM;
	electron1_ = electron1;
	electron2_ = electron2;
	muon3_ = muon3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}

TriLeptonsCandidate::TriLeptonsCandidate(const Muon_t &muon1, const Muon_t &muon2, const Muon_t &muon3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( muon1 );
	addDaughter( muon2 );
	addDaughter( muon3 );

	type_ = kMMM;
	muon1_ = muon1;
	muon2_ = muon2;
	muon3_ = muon3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}


TriLeptonsCandidate::TriLeptonsCandidate(const Muon_t &muon1, const Muon_t &muon2, const Electron_t &electron3, Vertex_ptr vertex )
{
	vertex_ = vertex;

	addDaughter( muon1 );
	addDaughter( muon2 );
	addDaughter( electron3 );

	type_ = kMME;
	muon1_ = muon1;
	muon2_ = muon2;
	electron3_ = electron3;

	AddFourMomenta addP4;  
	addP4.set( *this );
}


const Electron_t *TriLeptonsCandidate::leadingEle() const
{
	if (type_ == kEEE) {
		if( (daughter( 0 )->pt() > daughter( 1 )->pt()) && (daughter( 0 )->pt() > daughter( 2 )->pt()) ) {
			return dynamic_cast<const Electron_t *>( daughter( 0 ) );
		} else if( (daughter( 1 )->pt() > daughter( 0 )->pt()) && (daughter( 1 )->pt() > daughter( 2 )->pt()) ) {
			return dynamic_cast<const Electron_t *>( daughter( 1 ) );
		} else if( (daughter( 2 )->pt() > daughter( 0 )->pt()) && (daughter( 2 )->pt() > daughter( 1 )->pt()) ) {
			return dynamic_cast<const Electron_t *>( daughter( 2 ) );
		} else return nullptr;
	} else if (type_ == kEEM) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Electron_t *>( daughter( 0 ) );
		} else {
			return dynamic_cast<const Electron_t *>( daughter( 1 ) );
		}
	} else return nullptr;
}

const Electron_t *TriLeptonsCandidate::subLeadingEle() const
{
	if (type_ == kEEE) {
		if( ((daughter( 0 )->pt() < daughter( 1 )->pt()) && (daughter( 0 )->pt() > daughter( 2 )->pt())) 
			|| ((daughter( 0 )->pt() > daughter( 1 )->pt()) && (daughter( 0 )->pt() < daughter( 2 )->pt())) ) {
			return dynamic_cast<const Electron_t *>( daughter( 0 ) );
		} else if( ((daughter( 1 )->pt() < daughter( 0 )->pt()) && (daughter( 1 )->pt() > daughter( 2 )->pt())) 
			|| ((daughter( 1 )->pt() > daughter( 0 )->pt()) && (daughter( 1 )->pt() < daughter( 2 )->pt())) ) {
			return dynamic_cast<const Electron_t *>( daughter( 1 ) );
		} else if( ((daughter( 2 )->pt() < daughter( 0 )->pt()) && (daughter( 2 )->pt() > daughter( 1 )->pt())) 
			|| ((daughter( 2 )->pt() > daughter( 0 )->pt()) && (daughter( 2 )->pt() < daughter( 1 )->pt())) ) {
			return dynamic_cast<const Electron_t *>( daughter( 2 ) );
		} else return nullptr;
	} else if (type_ == kEEM) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Electron_t *>( daughter( 1 ) );
		} else {
			return dynamic_cast<const Electron_t *>( daughter( 0 ) );
		}
	} else return nullptr;
}

const Electron_t *TriLeptonsCandidate::softEle() const
{
	if (type_ == kEEE) {
		if( (daughter( 0 )->pt() < daughter( 1 )->pt()) && (daughter( 0 )->pt() < daughter( 2 )->pt()) ) {
			return dynamic_cast<const Electron_t *>( daughter( 0 ) );
		} else if( (daughter( 1 )->pt() < daughter( 0 )->pt()) && (daughter( 1 )->pt() < daughter( 2 )->pt()) ) {
			return dynamic_cast<const Electron_t *>( daughter( 1 ) );
		} else if( (daughter( 2 )->pt() < daughter( 0 )->pt()) && (daughter( 2 )->pt() < daughter( 1 )->pt()) ) {
			return dynamic_cast<const Electron_t *>( daughter( 2 ) );
		} else return nullptr;
	} else return nullptr;
}

const Electron_t *TriLeptonsCandidate::electron() const
{
	if (type_ == kMME) return dynamic_cast<const Electron_t *>( daughter( 2 ) ); 
	else return nullptr;
}


const Muon_t *TriLeptonsCandidate::leadingMuon() const
{
	if (type_ == kMMM) {
		if( (daughter( 0 )->pt() > daughter( 1 )->pt()) && (daughter( 0 )->pt() > daughter( 2 )->pt()) ) {
			return dynamic_cast<const Muon_t *>( daughter( 0 ) );
		} else if( (daughter( 1 )->pt() > daughter( 0 )->pt()) && (daughter( 1 )->pt() > daughter( 2 )->pt()) ) {
			return dynamic_cast<const Muon_t *>( daughter( 1 ) );
		} else if( (daughter( 2 )->pt() > daughter( 0 )->pt()) && (daughter( 2 )->pt() > daughter( 1 )->pt()) ) {
			return dynamic_cast<const Muon_t *>( daughter( 2 ) );
		} else return nullptr;
	} else if (type_ == kMME) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Muon_t *>( daughter( 0 ) );
		} else {
			return dynamic_cast<const Muon_t *>( daughter( 1 ) );
		}
	} else return nullptr;
}

const Muon_t *TriLeptonsCandidate::subLeadingMuon() const
{
	if (type_ == kMMM) {
		if( ((daughter( 0 )->pt() < daughter( 1 )->pt()) && (daughter( 0 )->pt() > daughter( 2 )->pt())) 
			|| ((daughter( 0 )->pt() > daughter( 1 )->pt()) && (daughter( 0 )->pt() < daughter( 2 )->pt())) ) {
			return dynamic_cast<const Muon_t *>( daughter( 0 ) );
		} else if( ((daughter( 1 )->pt() < daughter( 0 )->pt()) && (daughter( 1 )->pt() > daughter( 2 )->pt())) 
			|| ((daughter( 1 )->pt() > daughter( 0 )->pt()) && (daughter( 1 )->pt() < daughter( 2 )->pt())) ) {
			return dynamic_cast<const Muon_t *>( daughter( 1 ) );
		} else if( ((daughter( 2 )->pt() < daughter( 0 )->pt()) && (daughter( 2 )->pt() > daughter( 1 )->pt())) 
			|| ((daughter( 2 )->pt() > daughter( 0 )->pt()) && (daughter( 2 )->pt() < daughter( 1 )->pt())) ) {
			return dynamic_cast<const Muon_t *>( daughter( 2 ) );
		} else return nullptr;
	} else if (type_ == kMME) {
		if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
			return dynamic_cast<const Muon_t *>( daughter( 1 ) );
		} else {
			return dynamic_cast<const Muon_t *>( daughter( 0 ) );
		}
	} else return nullptr;
}

const Muon_t *TriLeptonsCandidate::softMuon() const
{
	if (type_ == kMMM) {
		if( (daughter( 0 )->pt() < daughter( 1 )->pt()) && (daughter( 0 )->pt() < daughter( 2 )->pt()) ) {
			return dynamic_cast<const Muon_t *>( daughter( 0 ) );
		} else if( (daughter( 1 )->pt() < daughter( 0 )->pt()) && (daughter( 1 )->pt() < daughter( 2 )->pt()) ) {
			return dynamic_cast<const Muon_t *>( daughter( 1 ) );
		} else if( (daughter( 2 )->pt() < daughter( 0 )->pt()) && (daughter( 2 )->pt() < daughter( 1 )->pt()) ) {
			return dynamic_cast<const Muon_t *>( daughter( 2 ) );
		} else return nullptr;
	} else return nullptr;
}

const Muon_t *TriLeptonsCandidate::muon() const
{
	if (type_ == kEEM) return dynamic_cast<const Muon_t *>( daughter( 2 ) );  
	else return nullptr;	
}



float TriLeptonsCandidate::sumPt() const {	
	if (type_ == kEEE) return ( leadingEle()->pt() + subLeadingEle()->pt() + softEle()->pt() );
	else if (type_ == kEEM) return ( leadingEle()->pt() + subLeadingEle()->pt() + muon()->pt() );
	else if (type_ == kMMM) return ( leadingMuon()->pt() + subLeadingMuon()->pt() + softMuon()->pt() );	
	else if (type_ == kMME) return ( leadingMuon()->pt() + subLeadingMuon()->pt() + electron()->pt() ); 
	else return 0.;
}  


float TriLeptonsCandidate::invMass() const {  
	TLorentzVector l1, l2, l3; 
	l1.SetPxPyPzE( 0., 0., 0., 0. );
	l2.SetPxPyPzE( 0., 0., 0., 0. );
	l3.SetPxPyPzE( 0., 0., 0., 0. );

	if (type_ == kEEE) {
		l1.SetPxPyPzE( leadingEle()->px(), leadingEle()->py(), leadingEle()->pz(), leadingEle()->energy() );
		l2.SetPxPyPzE( subLeadingEle()->px(), subLeadingEle()->py(), subLeadingEle()->pz(), subLeadingEle()->energy() );
		l3.SetPxPyPzE( softEle()->px(), softEle()->py(), softEle()->pz(), softEle()->energy() );
	} else if (type_ == kEEM) {
		l1.SetPxPyPzE( leadingEle()->px(), leadingEle()->py(), leadingEle()->pz(), leadingEle()->energy() );
		l2.SetPxPyPzE( subLeadingEle()->px(), subLeadingEle()->py(), subLeadingEle()->pz(), subLeadingEle()->energy() );
		l3.SetPxPyPzE( muon()->px(), muon()->py(), muon()->pz(), muon()->energy() );
	} else if (type_ == kMMM) {
		l1.SetPxPyPzE( leadingMuon()->px(), leadingMuon()->py(), leadingMuon()->pz(), leadingMuon()->energy() );
		l2.SetPxPyPzE( subLeadingMuon()->px(), subLeadingMuon()->py(), subLeadingMuon()->pz(), subLeadingMuon()->energy() );
		l3.SetPxPyPzE( softMuon()->px(), softMuon()->py(), softMuon()->pz(), softMuon()->energy() );
	} else if (type_ == kMME) {
		l1.SetPxPyPzE( leadingMuon()->px(), leadingMuon()->py(), leadingMuon()->pz(), leadingMuon()->energy() );
		l2.SetPxPyPzE( subLeadingMuon()->px(), subLeadingMuon()->py(), subLeadingMuon()->pz(), subLeadingMuon()->energy() );
		l3.SetPxPyPzE( electron()->px(), electron()->py(), electron()->pz(), electron()->energy() );
	} 
	return (l1+l2+l3).M();
}


float TriLeptonsCandidate::leadingLeptonPt() const {
	if (type_ == kEEE) return (leadingEle()->pt());
	else if (type_ == kMMM) return (leadingMuon()->pt());
	else if (type_ == kEEM) {
		if( muon()->pt() > leadingEle()->pt() ) return muon()->pt(); 
		else return leadingEle()->pt(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() > leadingMuon()->pt() ) return electron()->pt();
		else return leadingMuon()->pt();
	}
	else return 0.;
}

float TriLeptonsCandidate::subLeadingLeptonPt() const {
	if (type_ == kEEE) return (subLeadingEle()->pt());
	else if (type_ == kMMM) return (subLeadingMuon()->pt());
	else if (type_ == kEEM) {
		if( (muon()->pt() < leadingEle()->pt()) && (muon()->pt() > subLeadingEle()->pt()) ) return muon()->pt(); 
		else if( muon()->pt() > leadingEle()->pt() ) return leadingEle()->pt(); 
		else if( muon()->pt() <  subLeadingEle()->pt() ) return subLeadingEle()->pt(); 
		else return 0.;
	} 
	else if (type_ == kMME) {
		if( (electron()->pt() < leadingMuon()->pt()) && (electron()->pt() > subLeadingMuon()->pt()) ) return electron()->pt(); 
		else if( electron()->pt() > leadingMuon()->pt()) return leadingMuon()->pt(); 
		else if( electron()->pt() < subLeadingMuon()->pt()) return subLeadingMuon()->pt(); 
		else return 0.;
	} 	
	else return 0.;
}

float TriLeptonsCandidate::softLeptonPt() const {
	if (type_ == kEEE) return (softEle()->pt());
	else if (type_ == kMMM) return (softMuon()->pt());
	else if (type_ == kEEM) {
		if( muon()->pt() < subLeadingEle()->pt() ) return muon()->pt(); 
		else return subLeadingEle()->pt(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() < subLeadingMuon()->pt() ) return electron()->pt(); 
		else return subLeadingMuon()->pt(); 
	} 	
	else return 0.;
}


float TriLeptonsCandidate::leadingLeptonEta() const {
	if (type_ == kEEE) return (leadingEle()->eta());
	else if (type_ == kMMM) return (leadingMuon()->eta());
	else if (type_ == kEEM) {
		if( muon()->pt() > leadingEle()->pt() ) return muon()->eta(); 
		else return leadingEle()->eta(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() > leadingMuon()->pt() ) return electron()->eta();
		else return leadingMuon()->eta();
	}
	else return 0.;
}

float TriLeptonsCandidate::subLeadingLeptonEta() const {
	if (type_ == kEEE) return (subLeadingEle()->eta());
	else if (type_ == kMMM) return (subLeadingMuon()->eta());
	else if (type_ == kEEM) {
		if( (muon()->pt() < leadingEle()->pt()) && (muon()->pt() > subLeadingEle()->pt()) ) return muon()->eta(); 
		else if( muon()->pt() > leadingEle()->pt() ) return leadingEle()->eta(); 
		else if( muon()->pt() <  subLeadingEle()->pt() ) return subLeadingEle()->eta(); 
		else return 0.;
	} 
	else if (type_ == kMME) {
		if( (electron()->pt() < leadingMuon()->pt()) && (electron()->pt() > subLeadingMuon()->pt()) ) return electron()->eta(); 
		else if( electron()->pt() > leadingMuon()->pt()) return leadingMuon()->eta(); 
		else if( electron()->pt() < subLeadingMuon()->pt()) return subLeadingMuon()->eta(); 
		else return 0.;
	} 	
	else return 0.;
}

float TriLeptonsCandidate::softLeptonEta() const {
	if (type_ == kEEE) return (softEle()->eta());
	else if (type_ == kMMM) return (softMuon()->eta());
	else if (type_ == kEEM) {
		if( muon()->pt() < subLeadingEle()->pt() ) return muon()->eta(); 
		else return subLeadingEle()->eta(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() < subLeadingMuon()->pt() ) return electron()->pt(); 
		else return subLeadingMuon()->pt(); 
	} 	
	else return 0.;
}


float TriLeptonsCandidate::leadingLeptonPhi() const {
	if (type_ == kEEE) return (leadingEle()->phi());
	else if (type_ == kMMM) return (leadingMuon()->phi());
	else if (type_ == kEEM) {
		if( muon()->pt() > leadingEle()->pt() ) return muon()->phi(); 
		else return leadingEle()->phi(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() > leadingMuon()->pt() ) return electron()->phi();
		else return leadingMuon()->phi();
	}
	else return 0.;
}

float TriLeptonsCandidate::subLeadingLeptonPhi() const {
	if (type_ == kEEE) return (subLeadingEle()->phi());
	else if (type_ == kMMM) return (subLeadingMuon()->phi());
	else if (type_ == kEEM) {
		if( (muon()->pt() < leadingEle()->pt()) && (muon()->pt() > subLeadingEle()->pt()) ) return muon()->phi(); 
		else if( muon()->pt() > leadingEle()->pt() ) return leadingEle()->phi(); 
		else if( muon()->pt() <  subLeadingEle()->pt() ) return subLeadingEle()->phi(); 
		else return 0.;
	} 
	else if (type_ == kMME) {
		if( (electron()->pt() < leadingMuon()->pt()) && (electron()->pt() > subLeadingMuon()->pt()) ) return electron()->phi(); 
		else if( electron()->pt() > leadingMuon()->pt()) return leadingMuon()->phi(); 
		else if( electron()->pt() < subLeadingMuon()->pt()) return subLeadingMuon()->phi(); 
		else return 0.;
	} 	
	else return 0.;
}

float TriLeptonsCandidate::softLeptonPhi() const {
	if (type_ == kEEE) return (softEle()->phi());
	else if (type_ == kMMM) return (softMuon()->phi());
	else if (type_ == kEEM) {
		if( muon()->pt() < subLeadingEle()->pt() ) return muon()->phi(); 
		else return subLeadingEle()->phi(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() < subLeadingMuon()->pt() ) return electron()->phi(); 
		else return subLeadingMuon()->phi(); 
	} 	
	else return 0.;
}


float TriLeptonsCandidate::leadingLeptonCharge() const {
	if (type_ == kEEE) return (leadingEle()->charge());
	else if (type_ == kMMM) return (leadingMuon()->charge());
	else if (type_ == kEEM) {
		if( muon()->pt() > leadingEle()->pt() ) return muon()->charge(); 
		else return leadingEle()->charge(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() > leadingMuon()->pt() ) return electron()->charge();
		else return leadingMuon()->charge();
	}
	else return 0.;
}

float TriLeptonsCandidate::subLeadingLeptonCharge() const {
	if (type_ == kEEE) return (subLeadingEle()->charge());
	else if (type_ == kMMM) return (subLeadingMuon()->charge());
	else if (type_ == kEEM) {
		if( (muon()->pt() < leadingEle()->pt()) && (muon()->pt() > subLeadingEle()->pt()) ) return muon()->charge(); 
		else if( muon()->pt() > leadingEle()->pt() ) return leadingEle()->charge(); 
		else if( muon()->pt() <  subLeadingEle()->pt() ) return subLeadingEle()->charge(); 
		else return 0.;
	} 
	else if (type_ == kMME) {
		if( (electron()->pt() < leadingMuon()->pt()) && (electron()->pt() > subLeadingMuon()->pt()) ) return electron()->charge(); 
		else if( electron()->pt() > leadingMuon()->pt()) return leadingMuon()->charge(); 
		else if( electron()->pt() < subLeadingMuon()->pt()) return subLeadingMuon()->charge(); 
		else return 0.;
	} 	
	else return 0.;
}

float TriLeptonsCandidate::softLeptonCharge() const {
	if (type_ == kEEE) return (softEle()->charge());
	else if (type_ == kMMM) return (softMuon()->charge());
	else if (type_ == kEEM) {
		if( muon()->pt() < subLeadingEle()->pt() ) return muon()->charge(); 
		else return subLeadingEle()->charge(); 
	} 
	else if (type_ == kMME) {
		if( electron()->pt() < subLeadingMuon()->pt() ) return electron()->charge(); 
		else return subLeadingMuon()->charge(); 
	} 	
	else return 0.;
}




bool TriLeptonsCandidate::operator <( const TriLeptonsCandidate &b ) const
{
	return ( sumPt() < b.sumPt() );
}

bool TriLeptonsCandidate::operator >( const TriLeptonsCandidate &b ) const
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
