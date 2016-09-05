#ifndef FLASHgg_DiMuDiTrackTag_h
#define FLASHgg_DiMuDiTrackTag_h

#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


using namespace std;

typedef flashgg::Muon Muon_t;
typedef edm::Ptr<flashgg::Muon> Muon_ptr;
// typedef pat::Muon Muon_t;
// typedef edm::Ptr<pat::Muon> Muon_ptr;

typedef pat::PackedCandidate Track_t;
typedef edm::Ptr<pat::PackedCandidate> Track_ptr;



namespace flashgg {
    class DiMuDiTrackTag: public DiLeptonDiJetTagBase
    {
    public:
        DiMuDiTrackTag();
        DiMuDiTrackTag( edm::Ptr<DiLeptonDiJetCandidate> );  

        ~DiMuDiTrackTag();

        DiMuDiTrackTag *clone() const override { return ( new DiMuDiTrackTag( *this ) ); }

        const vector<Muon_ptr > muons() const { return Muons_;}
        const vector<Track_ptr > tracks() const { return Tracks_;}

        void setMuons( vector<Muon_ptr > Muons ) {Muons_ = Muons;}
        void setTracks( vector<Track_ptr > Tracks )  { Tracks_ = Tracks;}

		DiLeptonDiJetTagBase::tag_t tagEnum() const override {return DiLeptonDiJetTagBase::kDiMuDiTrack; }

    private:
        vector<Muon_ptr > Muons_;
        vector<Track_ptr > Tracks_;
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