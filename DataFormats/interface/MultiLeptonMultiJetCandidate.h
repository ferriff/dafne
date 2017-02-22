#ifndef FLASHgg_MultiLeptonMultiJetCandidate_h
#define FLASHgg_MultiLeptonMultiJetCandidate_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/DataFormats/interface/WeightedObject.h"
//#include "DataFormats/Math/interface/Point3D.h"
//#include "TLorentzVector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"  
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h" 
#include "flashgg/DataFormats/interface/Jet.h"


namespace flashgg {

        typedef flashgg::Electron Electron_t;  
        typedef edm::Ptr<flashgg::Electron> Electron_ptr; 

        typedef flashgg::Muon Muon_t;
        typedef edm::Ptr<flashgg::Muon> Muon_ptr;

        typedef pat::Jet Jet_t;
        typedef edm::Ptr<pat::Jet> Jet_ptr;

        typedef pat::PackedCandidate Track_t;
        typedef edm::Ptr<pat::PackedCandidate> Track_ptr;

        typedef reco::Vertex Vertex_t;
        typedef edm::Ptr<reco::Vertex> Vertex_ptr;

	class MultiLeptonMultiJetCandidate : public reco::LeafCandidate, public WeightedObject
	{
	public:
		enum CandidateType_t { kEEJJ, kMMJJ, kEETT, kMMTT, kEMJJ };

		MultiLeptonMultiJetCandidate();

		MultiLeptonMultiJetCandidate( Electron_ptr, Electron_ptr, Jet_ptr, Jet_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Electron_ptr, Electron_ptr, Track_ptr, Track_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Muon_ptr, Muon_ptr, Jet_ptr, Jet_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Muon_ptr, Muon_ptr, Track_ptr, Track_ptr, Vertex_ptr );
		MultiLeptonMultiJetCandidate( Electron_ptr, Muon_ptr, Jet_ptr, Jet_ptr, Vertex_ptr );

		//MultiLeptonMultiJetCandidate( const Electron_t &, const Electron_t &, Jet_ptr, Jet_ptr, Vertex_ptr );
		//MultiLeptonMultiJetCandidate( const Electron_t &, const Electron_t &, Track_ptr, Track_ptr, Vertex_ptr );
		//MultiLeptonMultiJetCandidate( const Muon_t &, const Muon_t &, Jet_ptr, Jet_ptr, Vertex_ptr );
		//MultiLeptonMultiJetCandidate( const Muon_t &, const Muon_t &, Track_ptr, Track_ptr, Vertex_ptr );
		//MultiLeptonMultiJetCandidate( const Electron_t &, const Muon_t &, Jet_ptr, Jet_ptr, Vertex_ptr );

		~MultiLeptonMultiJetCandidate();

		const CandidateType_t type() const { return type_; }

		const Vertex_ptr vtx() const { return vertex_; }
        void setVtx( Vertex_ptr val ) { vertex_ = val; }

		const std::vector<Electron_ptr> & electrons() const { return ptrEle_; }
        void embedElectrons();
		std::vector<Electron_t> & embeddedElectrons();

		const reco::Candidate * leadingLepton() const; 
		const reco::Candidate * subLeadingLepton() const; 

		const Electron_t *leadingEle() const; 
		const Electron_t *subLeadingEle() const; 

		const Muon_t *leadingMuon() const; 
		const Muon_t *subLeadingMuon() const;

		const Jet_t *leadingJet() const; 
		const Jet_t *subLeadingJet() const; 

		const Track_t *leadingTrack() const;
		const Track_t *subLeadingTrack() const; 

        //Point genPV() const { return genPV_; }

		float sumPt() const;
		float leptonInvMass() const; 

		bool operator <( const MultiLeptonMultiJetCandidate &b ) const;
		bool operator >( const MultiLeptonMultiJetCandidate &b ) const;

		bool isEEJJ() const { return type_ == kEEJJ; }

		bool isEETT() const { return type_ == kEETT; }

		bool isMMJJ() const { return type_ == kMMJJ; }

		bool isMMTT() const { return type_ == kMMTT; }

		bool isEMJJ() const { return type_ == kEMJJ; }

		MultiLeptonMultiJetCandidate *clone() const { return ( new MultiLeptonMultiJetCandidate( *this ) ); }

	private:
		
		CandidateType_t type_;

        std::vector<Electron_ptr> ptrEle_;
        std::vector<Electron_t> ele_;

        std::vector<Muon_ptr> ptrMuon_;
        std::vector<Muon_t>   muon_;

        std::vector<Jet_ptr> ptrJet_;
        std::vector<Jet_t>   jet_;

        std::vector<Track_ptr> ptrTrack_;
        std::vector<Track_t>   track_;

		Vertex_ptr vertex_;
        //Point genPV_;
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

