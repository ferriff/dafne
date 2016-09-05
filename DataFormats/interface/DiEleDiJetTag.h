#ifndef FLASHgg_DiEleDiJetTag_h
#define FLASHgg_DiEleDiJetTag_h

#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"
#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"  
#include "DataFormats/PatCandidates/interface/Jet.h"


using namespace std;

typedef flashgg::Electron Electron_t;  
typedef edm::Ptr<flashgg::Electron> Electron_ptr; 
// typedef pat::Electron Electron_t;  
// typedef edm::Ptr<pat::Electron> Electron_ptr; 

// typedef flashgg::Jet Jet_t;
// typedef edm::Ptr<flashgg::Jet> Jet_ptr;
typedef pat::Jet Jet_t;
typedef edm::Ptr<pat::Jet> Jet_ptr;



namespace flashgg {
    class DiEleDiJetTag: public DiLeptonDiJetTagBase
    {
    public:
        DiEleDiJetTag();
        DiEleDiJetTag( edm::Ptr<DiLeptonDiJetCandidate> );  

        ~DiEleDiJetTag();

        DiEleDiJetTag *clone() const override { return ( new DiEleDiJetTag( *this ) ); }

        const vector<Electron_ptr > electrons() const {return Electrons_;}
        const vector<Jet_ptr > jets() const { return Jets_;}

        void setElectrons( vector<Electron_ptr > Electrons ) {Electrons_ = Electrons;}
        void setJets( vector<Jet_ptr > Jets ) { Jets_ = Jets; }

		DiLeptonDiJetTagBase::tag_t tagEnum() const override {return DiLeptonDiJetTagBase::kDiEleDiJet; }

    private:
        vector<Electron_ptr > Electrons_;
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