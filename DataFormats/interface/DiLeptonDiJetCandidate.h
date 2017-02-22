#ifndef FLASHgg_DiLeptonDiJetCandidate_h
#define FLASHgg_DiLeptonDiJetCandidate_h

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/DataFormats/interface/WeightedObject.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TLorentzVector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"  
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h" 
#include "flashgg/DataFormats/interface/Jet.h"


using namespace std;

typedef flashgg::Electron Electron_t;  
typedef edm::Ptr<flashgg::Electron> Electron_ptr; 
// typedef pat::Electron Electron_t;  
// typedef edm::Ptr<pat::Electron> Electron_ptr; 

typedef flashgg::Muon Muon_t;
typedef edm::Ptr<flashgg::Muon> Muon_ptr;
// typedef pat::Muon Muon_t;
// typedef edm::Ptr<pat::Muon> Muon_ptr;

// typedef flashgg::Jet Jet_t;
// typedef edm::Ptr<flashgg::Jet> Jet_ptr;
typedef pat::Jet Jet_t;
typedef edm::Ptr<pat::Jet> Jet_ptr;

typedef pat::PackedCandidate Track_t;
typedef edm::Ptr<pat::PackedCandidate> Track_ptr;

typedef reco::Vertex Vertex_t;
typedef edm::Ptr<reco::Vertex> Vertex_ptr;



namespace flashgg {

	class DiLeptonDiJetCandidate : public WeightedObject, public reco::CompositeCandidate
	{
	public:
		DiLeptonDiJetCandidate();

		DiLeptonDiJetCandidate( Electron_ptr, Electron_ptr, const Jet_t &, const Jet_t &, Vertex_ptr );
		DiLeptonDiJetCandidate( Electron_ptr, Electron_ptr, Track_ptr, Track_ptr, Vertex_ptr );
		DiLeptonDiJetCandidate( Muon_ptr, Muon_ptr, const Jet_t &, const Jet_t &, Vertex_ptr );
		DiLeptonDiJetCandidate( Muon_ptr, Muon_ptr, Track_ptr, Track_ptr, Vertex_ptr );
		DiLeptonDiJetCandidate( Electron_ptr, Muon_ptr, const Jet_t &, const Jet_t &, Vertex_ptr );

		DiLeptonDiJetCandidate( const Electron_t &, const Electron_t &, Jet_ptr, Jet_ptr, Vertex_ptr );
		DiLeptonDiJetCandidate( const Electron_t &, const Electron_t &, Track_ptr, Track_ptr, Vertex_ptr );
		DiLeptonDiJetCandidate( const Muon_t &, const Muon_t &, Jet_ptr, Jet_ptr, Vertex_ptr );
		DiLeptonDiJetCandidate( const Muon_t &, const Muon_t &, Track_ptr, Track_ptr, Vertex_ptr );
		DiLeptonDiJetCandidate( const Electron_t &, const Muon_t &, Jet_ptr, Jet_ptr, Vertex_ptr );

		~DiLeptonDiJetCandidate();

		enum CandidateType_t { kEEJJ, kMMJJ, kEETT, kMMTT, kEMJJ };
		const CandidateType_t getType() const { return type_; }

		const Vertex_ptr vtx() const { return vertex_; }
        void setVtx( Vertex_ptr val ) { vertex_ = val; }

		const Electron_t *electron1() const; 
		const Electron_t *electron2() const; 
		const Electron_t *electron() const; 

		const Muon_t *muon1() const; 
		const Muon_t *muon2() const;
		const Muon_t *muon() const;

		const Electron_t *leadingEle() const; 
		const Electron_t *subLeadingEle() const; 

		const Muon_t *leadingMuon() const; 
		const Muon_t *subLeadingMuon() const;

		const Jet_t *leadingJet() const; 
		const Jet_t *subLeadingJet() const; 

		const Track_t *leadingTrack() const;
		const Track_t *subLeadingTrack() const; 


		void setLogSumPt2( float val ) { logsumpt2_ = val; }
		void setPtBal( float val ) { ptbal_ = val; }
		void setPtAsym( float val ) { ptasym_ = val; }

		void setVLogSumPt2( vector<float> vval ) { vlogsumpt2_ = vval; }
		void setVPtBal( vector<float> vval ) { vptbal_ = vval; }
		void setVPtAsym( vector<float> vval ) { vptasym_ = vval; }

		void setVertexIndex( int val ) { vertex_index_ = val; }
		void setJetCollectionIndex( unsigned int val ) { jetCollectionIndex_ = val; }
        void setGenPV( const Point genpv ) { genPV_ = genpv; }


		float logSumPt2() const { return logsumpt2_; }
		float ptBal() const { return ptbal_; }
		float ptAsym() const { return ptasym_; }

		float nVert() const { return nVert_; }
		unsigned int nVtxInfoSize() const { return ( vlogsumpt2_.size() ) ;}
		Vertex_ptr vertexPtr( unsigned int iVtx ) const  { return iVtx < vVtxPtr_.size() ? vVtxPtr_.at( iVtx ) : Vertex_ptr(); }

		float logSumPt2( unsigned int iVtx ) const { return ( iVtx < vlogsumpt2_.size() ) ? vlogsumpt2_.at( iVtx ) : -9999. ;}
		float ptBal( unsigned int iVtx ) const  { return iVtx < vptbal_.size() ? vptbal_.at( iVtx ) : -9999. ;}
		float ptAsym( unsigned int iVtx ) const  { return iVtx < vptasym_.size() ? vptasym_.at( iVtx ) : -9999. ; }

		int vertexIndex() const { return vertex_index_; }
		unsigned int jetCollectionIndex() const { return jetCollectionIndex_; }
        Point genPV() const { return genPV_; }

		float sumPt() const;
        float invMass() const; 
		float diLeptonInvMass() const; 

		float leadingLeptonPt() const; 
		float subLeadingLeptonPt() const; 
		float leadingLeptonEta() const;
		float subLeadingLeptonEta() const; 		
		float leadingLeptonPhi() const;
		float subLeadingLeptonPhi() const; 	
		float leadingLeptonCharge() const; 
		float subLeadingLeptonCharge() const; 


		bool operator <( const DiLeptonDiJetCandidate &b ) const;
		bool operator >( const DiLeptonDiJetCandidate &b ) const;

		bool isEEJJ() const {
			if (type_ == kEEJJ) return true;
			else return false;
		}

		bool isEETT() const {
			if (type_ == kEETT) return true;
			else return false;
		}

		bool isMMJJ() const {
			if (type_ == kMMJJ) return true;
			else return false;
		}

		bool isMMTT() const {
			if (type_ == kMMTT) return true;
			else return false;
		}

		bool isEMJJ() const {
			if (type_ == kEMJJ) return true;
			else return false;
		}

		DiLeptonDiJetCandidate *clone() const { return ( new DiLeptonDiJetCandidate( *this ) ); }


	private:
		
		CandidateType_t type_;

		Electron_t electron1_, electron2_, electron_;
		Muon_t muon1_, muon2_, muon_;
		Jet_t jet1_, jet2_;

		Electron_ptr electrons_[2], electronPtr_;
		Muon_ptr muons_[2], muonPtr_;
		Jet_ptr jets_[2];
		Track_ptr tracks_[2];  

		Vertex_ptr vertex_;

		float logsumpt2_;
		float ptbal_;
		float ptasym_;
		float nVert_;
		int vertex_index_;

		vector<float> vlogsumpt2_;
		vector<float> vptbal_;
		vector<float> vptasym_;
		vector< Vertex_ptr > vVtxPtr_;
		unsigned int jetCollectionIndex_; // index for which jet collection corresponds to the vertex choice in this diphoton
        Point genPV_;

	};


}


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

