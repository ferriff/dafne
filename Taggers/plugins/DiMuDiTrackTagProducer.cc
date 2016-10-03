#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "dafne/DataFormats/interface/DiMuDiTrackTag.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>

using namespace std;
using namespace edm;


typedef flashgg::Muon Muon_t;
typedef edm::Ptr<flashgg::Muon> Muon_ptr;
// typedef pat::Muon Muon_t;
// typedef edm::Ptr<pat::Muon> Muon_ptr;

typedef pat::PackedCandidate Track_t;
typedef edm::Ptr<pat::PackedCandidate> Track_ptr;


namespace flashgg {

    class DiMuDiTrackTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        DiMuDiTrackTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        EDGetTokenT<View<DiLeptonDiJetCandidate> > DiLeptonDiJetToken_;
        EDGetTokenT<View<Muon_t> > muonToken_;
        EDGetTokenT<View<Track_t> > trackToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        // string systLabel_;

        //---thresholds---
        //---tracks
        double trackPtThreshold_;
        double trackEtaThreshold_;
        // //leptons
        double muonPtThreshold_;   
        double muonEtaThreshold_;
        bool useStdLeptonID_;     //??
        double muPFIsoSumRelThreshold_;
        double muMiniIsoSumRelThreshold_;

    };

    DiMuDiTrackTagProducer::DiMuDiTrackTagProducer( const ParameterSet &iConfig ) :
        DiLeptonDiJetToken_( consumes<View<flashgg::DiLeptonDiJetCandidate> >( iConfig.getParameter<InputTag> ( "DiLeptonDiJetTag" ) ) ),
        muonToken_( consumes<View<Muon_t> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        trackToken_( consumes<View<Track_t> >( iConfig.getParameter<InputTag> ( "TrackTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) )//,
        // systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
    {
        trackPtThreshold_ = iConfig.getParameter<double>( "trackPtThreshold" );
        trackEtaThreshold_ = iConfig.getParameter<double>( "trackEtaThreshold" );
        muonPtThreshold_ = iConfig.getParameter<double>( "muonPtThreshold" );
        muonEtaThreshold_ = iConfig.getParameter<double>( "muonEtaThreshold" );
        useStdLeptonID_=iConfig.getParameter<bool>( "useStdLeptonID" );  //??
        muPFIsoSumRelThreshold_ = iConfig.getParameter<double>( "muPFIsoSumRelThreshold" );
        muMiniIsoSumRelThreshold_ = iConfig.getParameter<double>( "muMiniIsoSumRelThreshold" );

        produces<vector<DiMuDiTrackTag> >();
        produces<vector<TagTruthBase> >();
    }


    void DiMuDiTrackTagProducer::produce( Event &evt, const EventSetup & )
    {
        Handle<View<flashgg::DiLeptonDiJetCandidate> > diLeptonDiJets;
        evt.getByToken( DiLeptonDiJetToken_, diLeptonDiJets );
        // const PtrVector<flashgg::DiLeptonDiJetCandidate>& diLeptonDiJetsPointers = diLeptonDiJets->ptrVector();

        Handle<View<Track_t> > theTracks;
        evt.getByToken( trackToken_, theTracks );

        Handle<View<Muon_t> > theMuons;
        evt.getByToken( muonToken_, theMuons );

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        Handle<View<reco::GenParticle> > genParticles;

        auto_ptr<vector<DiMuDiTrackTag> > dmdttags( new vector<DiMuDiTrackTag> );
        auto_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );


        Point Vtx;
        if( ! evt.isRealData() ) {
            evt.getByToken( genParticleToken_, genParticles );
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 ) {   //Higgs or gamma -> cosa metto? Wright?
                    Vtx = genParticles->ptrAt( genLoop )->vertex();
                    break;
                }
            }
        }

        edm::RefProd<vector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<vector<TagTruthBase> >();
        unsigned int idx = 0;        

        vector<Muon_ptr> goodMuons;
        if( !useStdLeptonID_) {
            goodMuons = selectAllMuonsSum16( theMuons->ptrs(), vertices->ptrs(), muonEtaThreshold_ , muonPtThreshold_, muMiniIsoSumRelThreshold_ );
        } else {
            goodMuons = selectAllMuons( theMuons->ptrs(), vertices->ptrs(), muonEtaThreshold_ , muonPtThreshold_, muPFIsoSumRelThreshold_ );
        }  // -> da controllare
        
        
        for( unsigned int diLeptonDiJetIndex = 0; diLeptonDiJetIndex < diLeptonDiJets->size(); diLeptonDiJetIndex++ ) {

            edm::Ptr<flashgg::DiLeptonDiJetCandidate> diLeptonDiJet = diLeptonDiJets->ptrAt( diLeptonDiJetIndex );

            if( goodMuons.size() < 2 && theTracks->size() < 2)  continue;  //o uguale a 2? (DiLeptonDiJet ne dovrebbe avere solo 2 memorizzati, ma qui prendo i GOOD electrons)
            // if( theMuons.size() < 2 && theTracks->size() < 2)  continue;  //or
            // if (theMuons.size() == 2 && theTracks->size() == 2) {}

            // DiMuDiTrackTag dmdttags_obj( diLeptonDiJet, mvares, JetVect, BJetVect );
            DiMuDiTrackTag dmdttags_obj( diLeptonDiJet );
            // devo prendere il diLeptonDiJet se ho 2mu e 2track


            // dmdttags_obj.setDiPhotonIndex( diLeptonDiJetIndex );  //non ce li ho
            // dmdttags_obj.setSystLabel( systLabel_ );
            // for( unsigned num = 0; num < JetVect.size(); num++ ) {
            //     dmdttags_obj.includeWeights( *JetVect[num] );
            // }
            // dmdttags_obj.includeWeights( *diLeptonDiJet ); 

            dmdttags->push_back( dmdttags_obj );


            if( ! evt.isRealData() ) {
                TagTruthBase truth_obj;
                truth_obj.setGenPV( Vtx );
                truths->push_back( truth_obj );
                dmdttags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, idx++ ) ) );
            }

            // count++;
        }

        evt.put( dmdttags );
        evt.put( truths );
        // cout << "tagged events = " << count << endl;
    }
}
typedef flashgg::DiMuDiTrackTagProducer FlashggDiMuDiTrackTagProducer;
DEFINE_FWK_MODULE( FlashggDiMuDiTrackTagProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
