#ifndef FLASHgg_TriLeptonsCandidate_h
#define FLASHgg_TriLeptonsCandidate_h

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/DataFormats/interface/WeightedObject.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TLorentzVector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"  
#include "DataFormats/PatCandidates/interface/Muon.h" 

#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h" 



using namespace std;

typedef flashgg::Electron Electron_t;  
typedef edm::Ptr<flashgg::Electron> Electron_ptr; 
// typedef pat::Electron Electron_t;  
// typedef edm::Ptr<pat::Electron> Electron_ptr; 

typedef flashgg::Muon Muon_t;
typedef edm::Ptr<flashgg::Muon> Muon_ptr;
// typedef pat::Muon Muon_t;
// typedef edm::Ptr<pat::Muon> Muon_ptr;

typedef reco::Vertex Vertex_t;
typedef edm::Ptr<reco::Vertex> Vertex_ptr;



namespace flashgg {

	class TriLeptonsCandidate : public WeightedObject, public reco::CompositeCandidate
	{
	public:
		TriLeptonsCandidate();

		TriLeptonsCandidate( Electron_ptr, Electron_ptr, Electron_ptr, Vertex_ptr );
		TriLeptonsCandidate( Electron_ptr, Electron_ptr, Muon_ptr, Vertex_ptr );
		TriLeptonsCandidate( Muon_ptr, Muon_ptr, Muon_ptr, Vertex_ptr );
		TriLeptonsCandidate( Muon_ptr, Muon_ptr, Electron_ptr, Vertex_ptr );

		TriLeptonsCandidate( const Electron_t &, const Electron_t &, const Electron_t &, Vertex_ptr );
		TriLeptonsCandidate( const Electron_t &, const Electron_t &, const Muon_t &, Vertex_ptr );
		TriLeptonsCandidate( const Muon_t &, const Muon_t &, const Muon_t &, Vertex_ptr );
		TriLeptonsCandidate( const Muon_t &, const Muon_t &, const Electron_t &, Vertex_ptr );

		~TriLeptonsCandidate();

		enum CandidateType_t { kEEE, kEEM, kMMM, kMME };
		const CandidateType_t getType() const { return type_; }

		const Vertex_ptr vtx() const { return vertex_; }
        void setVtx( Vertex_ptr val ) { vertex_ = val; }

		const Electron_t *leadingEle() const; 
		const Electron_t *subLeadingEle() const; 
		const Electron_t *softEle() const;
		const Electron_t *electron() const; 

		const Muon_t *leadingMuon() const; 
		const Muon_t *subLeadingMuon() const;
		const Muon_t *softMuon() const;
		const Muon_t *muon() const;

		void setLogSumPt2( float val ) { logsumpt2_ = val; }
		void setPtBal( float val ) { ptbal_ = val; }
		void setPtAsym( float val ) { ptasym_ = val; }

		void setVLogSumPt2( vector<float> vval ) { vlogsumpt2_ = vval; }
		void setVPtBal( vector<float> vval ) { vptbal_ = vval; }
		void setVPtAsym( vector<float> vval ) { vptasym_ = vval; }

		void setVertexIndex( int val ) { vertex_index_ = val; }
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
        Point genPV() const { return genPV_; }


		float sumPt() const;
        float invMass() const; 

		float leadingLeptonPt() const; 
		float subLeadingLeptonPt() const; 
		float softLeptonPt() const;

		float leadingLeptonEta() const;
		float subLeadingLeptonEta() const; 		
		float softLeptonEta() const;

		float leadingLeptonPhi() const;
		float subLeadingLeptonPhi() const; 	
		float softLeptonPhi() const;

		float leadingLeptonCharge() const; 
		float subLeadingLeptonCharge() const; 
		float softLeptonCharge() const;


		bool operator <( const TriLeptonsCandidate &b ) const;
		bool operator >( const TriLeptonsCandidate &b ) const;

		bool isEEE() const {
			if (type_ == kEEE) return true;
			else return false;
		}

		bool isEEM() const {
			if (type_ == kEEM) return true;
			else return false;
		}

		bool isMMM() const {
			if (type_ == kMMM) return true;
			else return false;
		}

		bool isMME() const {
			if (type_ == kMME) return true;
			else return false;
		}

		TriLeptonsCandidate *clone() const { return ( new TriLeptonsCandidate( *this ) ); }


	private:
		
		CandidateType_t type_;

		Electron_t electron1_, electron2_, electron3_;
		Muon_t muon1_, muon2_, muon3_;

		Electron_ptr electrons_[3], electronPtr_;
		Muon_ptr muons_[3], muonPtr_;  

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

