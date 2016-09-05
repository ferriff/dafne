#ifndef FLASHgg_DiMuDiJetTag_h
#define FLASHgg_DiMuDiJetTag_h

#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/Jet.h"


using namespace std;

typedef flashgg::Muon Muon_t;
typedef edm::Ptr<flashgg::Muon> Muon_ptr;
// typedef pat::Muon Muon_t;
// typedef edm::Ptr<pat::Muon> Muon_ptr;

// typedef flashgg::Jet Jet_t;
// typedef edm::Ptr<flashgg::Jet> Jet_ptr;
typedef pat::Jet Jet_t;
typedef edm::Ptr<pat::Jet> Jet_ptr;



namespace flashgg {
    class DiMuDiJetTag: public DiLeptonDiJetTagBase
    {
    public:
        DiMuDiJetTag();
        DiMuDiJetTag( edm::Ptr<DiLeptonDiJetCandidate> );  

        ~DiMuDiJetTag();

        DiMuDiJetTag *clone() const override { return ( new DiMuDiJetTag( *this ) ); }

        const vector<Muon_ptr > muons() const { return Muons_;}
        const vector<Jet_ptr > jets() const { return Jets_;}

        void setMuons( vector<Muon_ptr > Muons ) {Muons_ = Muons;}
        void setJets( vector<Jet_ptr > Jets ) { Jets_ = Jets; }

        DiLeptonDiJetTagBase::tag_t tagEnum() const override {return DiLeptonDiJetTagBase::kDiMuDiJet; }

    private:
        vector<Muon_ptr > Muons_;        
        vector<Jet_ptr > Jets_;
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