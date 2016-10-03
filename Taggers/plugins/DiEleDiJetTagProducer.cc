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

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "dafne/DataFormats/interface/DiEleDiJetTag.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>

using namespace std;
using namespace edm;


typedef flashgg::Electron Electron_t;  
typedef edm::Ptr<flashgg::Electron> Electron_ptr; 
// typedef pat::Electron Electron_t;  
// typedef edm::Ptr<pat::Electron> Electron_ptr; 

// typedef flashgg::Jet Jet_t;
// typedef edm::Ptr<flashgg::Jet> Jet_ptr;
typedef pat::Jet Jet_t;
typedef edm::Ptr<pat::Jet> Jet_ptr;



namespace flashgg {

    class DiEleDiJetTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        DiEleDiJetTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        EDGetTokenT<View<DiLeptonDiJetCandidate> > DiLeptonDiJetToken_;
        vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<Electron_t> > electronToken_;
        // EDGetTokenT<View<Jet_t> > jetToken_; 
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        // string systLabel_;
        vector<edm::EDGetTokenT<View<Jet_t> > > tokenJets_;

        typedef std::vector<edm::Handle<edm::View<Jet_t> > > JetCollectionVector;
        //---thresholds---
        //---jets
        double jetPtThreshold_;
        double jetEtaThreshold_;
        // //leptons
        double electronPtThreshold_;   
        vector<double>  electronEtaThresholds_;
        bool useStdLeptonID_;  //??
        bool useElectronMVARecipe_;  
        bool useElectronLooseID_;   
        // double electronIsoThreshold_;  
        double elMiniIsoEBThreshold_;  
        double elMiniIsoEEThreshold_;  
        // double electronNumOfHitsThreshold_;   
        // double TransverseImpactParam_;
        // double LongitudinalImpactParam_;

    };

    DiEleDiJetTagProducer::DiEleDiJetTagProducer( const ParameterSet &iConfig ) :
        DiLeptonDiJetToken_( consumes<View<flashgg::DiLeptonDiJetCandidate> >( iConfig.getParameter<InputTag> ( "DiLeptonDiJetTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        electronToken_( consumes<View<Electron_t> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        // jetToken_( consumes<View<Jet_t> >( iConfig.getParameter<InputTag> ( "JetTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) )//,
        // systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
    {
        jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold" );
        jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold" );
        electronPtThreshold_ = iConfig.getParameter<double>( "electronPtThreshold" );
        electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds" );
        useStdLeptonID_=iConfig.getParameter<bool>( "useStdLeptonID" );  //??
        useElectronMVARecipe_=iConfig.getParameter<bool>( "useElectronMVARecipe" );  
        useElectronLooseID_=iConfig.getParameter<bool>("useElectronLooseID");  
        // electronIsoThreshold_ = iConfig.getParameter<double>( "electronIsoThreshold");  
        elMiniIsoEBThreshold_ = iConfig.getParameter<double>( "elMiniIsoEBThreshold" );  
        elMiniIsoEEThreshold_ = iConfig.getParameter<double>( "elMiniIsoEEThreshold" );  
        // electronNumOfHitsThreshold_ = iConfig.getParameter<double>( "electronNumOfHitsThreshold");  //??
        // TransverseImpactParam_ = iConfig.getParameter<double>( "TransverseImpactParam");
        // LongitudinalImpactParam_ = iConfig.getParameter<double>( "LongitudinalImpactParam");
        
        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<Jet_t> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<DiEleDiJetTag> >();
        produces<vector<TagTruthBase> >();
    }


    void DiEleDiJetTagProducer::produce( Event &evt, const EventSetup & )
    {
        //Handle<View<Jet_t> > theJets;
        //evt.getByToken( jetToken_, theJets );
        // const PtrVector<Jet_t>& jetPointers = theJets->ptrVector(); //questo o quello sotto??
        JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::DiLeptonDiJetCandidate> > diLeptonDiJets;
        evt.getByToken( DiLeptonDiJetToken_, diLeptonDiJets );
        // const PtrVector<flashgg::DiLeptonDiJetCandidate>& diLeptonDiJetsPointers = diLeptonDiJets->ptrVector();

        Handle<View<Electron_t> > theElectrons;
        evt.getByToken( electronToken_, theElectrons );

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        Handle<View<reco::GenParticle> > genParticles;

        auto_ptr<vector<DiEleDiJetTag> > dedjtags( new vector<DiEleDiJetTag> );
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

        vector<Electron_ptr> goodElectrons;
        if( !useStdLeptonID_ ) {
             // goodElectrons = selectAllElectrons( theElectrons->ptrs(), vertices->ptrs(), electronPtThreshold_, 
             //                                     TransverseImpactParam_, LongitudinalImpactParam_, nonTrigMVAThresholds_, nonTrigMVAEtaCuts_,
             //                                     electronIsoThreshold_, electronNumOfHitsThreshold_, electronEtaThresholds_ );
            goodElectrons = selectAllElectronsSum16( theElectrons->ptrs(), vertices->ptrs(), electronPtThreshold_, electronEtaThresholds_,
                                                     true, true, elMiniIsoEBThreshold_, elMiniIsoEEThreshold_);
        } else {
            goodElectrons = selectStdAllElectrons(theElectrons->ptrs(), vertices->ptrs(), electronPtThreshold_, electronEtaThresholds_,
                                                  useElectronMVARecipe_, useElectronLooseID_);
        }  // -> da controllare
        
        
        for( unsigned int diLeptonDiJetIndex = 0; diLeptonDiJetIndex < diLeptonDiJets->size(); diLeptonDiJetIndex++ ) {

            edm::Ptr<flashgg::DiLeptonDiJetCandidate> diLeptonDiJet = diLeptonDiJets->ptrAt( diLeptonDiJetIndex );

            unsigned int jetCollectionIndex = diLeptonDiJets->ptrAt( diLeptonDiJetIndex )->jetCollectionIndex();  //diLeptonDiJet->jetCollectionIndex();
            vector<Jet_ptr> JetVect;

            if( goodElectrons.size() < 2 && Jets[jetCollectionIndex]->size() < 2)  continue;  //o uguale a 2? (DiLeptonDiJet ne dovrebbe avere solo 2 memorizzati, ma qui prendo i GOOD electrons)

            for( unsigned int jetIndex = 0; jetIndex < Jets[jetCollectionIndex]->size() ; jetIndex++ ) {  //o solo fino a 2?
                Jet_ptr thejet = Jets[jetCollectionIndex]->ptrAt( jetIndex );
                if( fabs( thejet->eta() ) > jetEtaThreshold_ && thejet->pt() < jetPtThreshold_ ) { continue; }
                JetVect.push_back( thejet );   //mi serve??
            }


            // if( theElectrons.size() < 2 && Jets[jetCollectionIndex]->size() < 2)  continue;  //or
            // if (theElectrons.size() == 2 && Jets[jetCollectionIndex]->size() == 2) {}

            // DiEleDiJetTag dedjtags_obj( diLeptonDiJet, mvares, JetVect, BJetVect );
            DiEleDiJetTag dedjtags_obj( diLeptonDiJet );
            // devo prendere il diLeptonDiJet se ho 2ele e 2jet


            // dedjtags_obj.setDiPhotonIndex( diLeptonDiJetIndex );  //non ce li ho
            // dedjtags_obj.setSystLabel( systLabel_ );
            // for( unsigned num = 0; num < JetVect.size(); num++ ) {
            //     dedjtags_obj.includeWeights( *JetVect[num] );
            // }
            // dedjtags_obj.includeWeights( *diLeptonDiJet ); 

            dedjtags->push_back( dedjtags_obj );


            if( ! evt.isRealData() ) {
                TagTruthBase truth_obj;
                truth_obj.setGenPV( Vtx );
                truths->push_back( truth_obj );
                dedjtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, idx++ ) ) );
            }

            // count++;
        }

        evt.put( dedjtags );
        evt.put( truths );
        // cout << "tagged events = " << count << endl;
    }
}
typedef flashgg::DiEleDiJetTagProducer FlashggDiEleDiJetTagProducer;
DEFINE_FWK_MODULE( FlashggDiEleDiJetTagProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
