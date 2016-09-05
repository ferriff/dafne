#ifndef FLASHgg_DiEleDiTrackTag_h
#define FLASHgg_DiEleDiTrackTag_h

#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"  
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


using namespace std;

typedef flashgg::Electron Electron_t;  
typedef edm::Ptr<flashgg::Electron> Electron_ptr; 
// typedef pat::Electron Electron_t;  
// typedef edm::Ptr<pat::Electron> Electron_ptr; 

typedef pat::PackedCandidate Track_t;
typedef edm::Ptr<pat::PackedCandidate> Track_ptr;



namespace flashgg {
    class DiEleDiTrackTag: public DiLeptonDiJetTagBase
    {
    public:
        DiEleDiTrackTag();
        DiEleDiTrackTag( edm::Ptr<DiLeptonDiJetCandidate> );  

        ~DiEleDiTrackTag();

        DiEleDiTrackTag *clone() const override { return ( new DiEleDiTrackTag( *this ) ); }

        const vector<Electron_ptr > electrons() const {return Electrons_;}
        const vector<Track_ptr > tracks() const { return Tracks_;}

        void setElectrons( vector<Electron_ptr > Electrons ) {Electrons_ = Electrons;}
        void setTracks( vector<Track_ptr > Tracks )  { Tracks_ = Tracks;}

		DiLeptonDiJetTagBase::tag_t tagEnum() const override {return DiLeptonDiJetTagBase::kDiEleDiTrack; }

    private:
        vector<Electron_ptr > Electrons_;
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