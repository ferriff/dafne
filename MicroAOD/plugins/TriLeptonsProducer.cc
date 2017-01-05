#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "dafne/DataFormats/interface/TriLeptonsCandidate.h"


using namespace edm;
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

	class TriLeptonsProducer : public EDProducer
	{

	public:
		TriLeptonsProducer( const ParameterSet & );

	private:
		void produce( Event &, const EventSetup & ) override;
		EDGetTokenT<View<Electron_t> > electronToken_;
		EDGetTokenT<View<Muon_t> > muonToken_;      	
		EDGetTokenT<View<Vertex_t> > vertexToken_;

		double minElePt_;
		double maxEleEta_;
		double minMuPt_;
		double maxMuEta_;

	};

	TriLeptonsProducer::TriLeptonsProducer( const ParameterSet &iConfig ) :
		electronToken_( consumes<View<Electron_t> >( iConfig.getParameter<InputTag> ( "ElectronTag" ) ) ),
		muonToken_( consumes<View<Muon_t> >( iConfig.getParameter<InputTag> ( "MuonTag" ) ) ),
		vertexToken_( consumes<View<Vertex_t> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) )

	{
		minElePt_ = iConfig.getParameter<double>( "minElectronPt" );
		maxEleEta_ = iConfig.getParameter<double>( "maxElectronEta" );		
		minMuPt_ = iConfig.getParameter<double>( "minMuonPt" );
		maxMuEta_ = iConfig.getParameter<double>( "maxMuonEta" );	

		produces<vector<flashgg::TriLeptonsCandidate> >();
	}

	void TriLeptonsProducer::produce( Event &evt, const EventSetup & )
	{
		Handle<View<Vertex_t> > primaryVertices;
		evt.getByToken( vertexToken_, primaryVertices );
		const vector<Vertex_ptr> &pvPointers = primaryVertices->ptrs();
		Vertex_ptr pvx = pvPointers[0]; //selected vertex 0

		Handle<View<Electron_t> > electrons;
		evt.getByToken( electronToken_, electrons );
		const vector<Electron_ptr > &electronPointers = electrons->ptrs();

		Handle<View<Muon_t> > muons;
		evt.getByToken( muonToken_, muons );
		const vector<Muon_ptr > &muonPointers = muons->ptrs();

		auto_ptr<vector<flashgg::TriLeptonsCandidate> > TriLeptonsColl( new vector<flashgg::TriLeptonsCandidate> );
		//    cout << "evt.id().event()= " << evt.id().event() << "\tevt.isRealData()= " << evt.isRealData() << "\tmuonPointers.size()= " << muonPointers.size() << "\tpvPointers.size()= " << pvPointers.size() << endl;



		if (electronPointers.size() > 1) {

			for( unsigned int i = 0 ; i < electronPointers.size() ; i++ ) {
				Electron_ptr electron1 = electronPointers[i];
				double pt_ele1 = electron1->pt();
				double eta_ele1 = electron1->eta();
				if( pt_ele1 < minElePt_ || fabs( eta_ele1 ) > maxEleEta_ ) { continue; }

				for( unsigned int j = i + 1 ; j < electronPointers.size() ; j++ ) {
					Electron_ptr electron2 = electronPointers[j];
					double pt_ele2 = electron2->pt();
					double eta_ele2 = electron2->eta();
					if( pt_ele2 < minElePt_ || fabs( eta_ele2 ) > maxEleEta_ ) { continue; }

					for( unsigned int l = j + 1 ; l < electronPointers.size() ; l++ ) {
						Electron_ptr electron3 = electronPointers[l];
						double pt_ele3 = electron3->pt();
						double eta_ele3 = electron3->eta();
						if( pt_ele3 < minElePt_ || fabs( eta_ele3 ) > maxEleEta_ ) { continue; }

						TriLeptonsCandidate TriEle( electron1, electron2, electron3, pvx);  //create TriLeptonsCandidate with 3ele
						TriEle.setVtx( pvx );  

						int ivtx = 0;
						for( unsigned int k = 0; k < primaryVertices->size() ; k++ ) {
							if( pvx == primaryVertices->ptrAt( k ) ) {
								ivtx = k;
								break;
							}
						}
						TriEle.setVertexIndex( ivtx );               				

						TriLeptonsColl->push_back( TriEle );  // store the TriEle into the collection
					}


					if (muonPointers.size() > 1) {
						for( unsigned int i = 0 ; i < muonPointers.size() ; i++ ) {
							Muon_ptr muon = muonPointers[i];
							double pt_mu = muon->pt();
							double eta_mu = muon->eta();
							if( pt_mu < minMuPt_ || fabs( eta_mu ) > maxMuEta_ ) { continue; }

							TriLeptonsCandidate DiEleMuon( electron1, electron2, muon, pvx);  //create TriLeptonsCandidate with 2ele and 1muon
							DiEleMuon.setVtx( pvx );

							int ivtx = 0;
							for( unsigned int k = 0; k < primaryVertices->size() ; k++ ) {
								if( pvx == primaryVertices->ptrAt( k ) ) {
									ivtx = k;
									break;
								}	
							}
							DiEleMuon.setVertexIndex( ivtx );				

							TriLeptonsColl->push_back( DiEleMuon );  // store the DiEleMuon into the collection
						}
					}

				}
			}
		} 



		if (muonPointers.size() > 1) {

			for( unsigned int i = 0 ; i < muonPointers.size() ; i++ ) {
				Muon_ptr muon1 = muonPointers[i];
				double pt_mu1 = muon1->pt();
				double eta_mu1 = muon1->eta();
				if( pt_mu1 < minMuPt_ || fabs( eta_mu1 ) > maxMuEta_ ) { continue; }
				
				for( unsigned int j = i + 1 ; j < muonPointers.size() ; j++ ) {
					Muon_ptr muon2 = muonPointers[j];
					double pt_mu2 = muon2->pt();
					double eta_mu2 = muon2->eta();
					if( pt_mu2 < minMuPt_ || fabs( eta_mu2 ) > maxMuEta_ ) { continue; }

					for( unsigned int l = j + 1 ; l < muonPointers.size() ; l++ ) {
						Muon_ptr muon3 = muonPointers[l];
						double pt_mu3 = muon3->pt();
						double eta_mu3 = muon3->eta();
						if( pt_mu3 < minMuPt_ || fabs( eta_mu3 ) > maxMuEta_ ) { continue; }

						TriLeptonsCandidate TriMuon( muon1, muon2, muon3, pvx);  //create TriLeptonsCandidate with 3muons
						TriMuon.setVtx( pvx );

						int ivtx = 0;
						for( unsigned int k = 0; k < primaryVertices->size() ; k++ ) {
							if( pvx == primaryVertices->ptrAt( k ) ) {
								ivtx = k;
								break;
							}	
						}
						TriMuon.setVertexIndex( ivtx );

						TriLeptonsColl->push_back( TriMuon );  // store the TriMuon into the collection
					}


					if (electronPointers.size() > 1) {
						for( unsigned int i = 0 ; i < electronPointers.size() ; i++ ) {
							Electron_ptr electron = electronPointers[i];
							double pt_ele = electron->pt();
							double eta_ele = electron->eta();
							if( pt_ele < minElePt_ || fabs( eta_ele ) > maxEleEta_ ) { continue; }

							TriLeptonsCandidate DiMuonEle( muon1, muon2, electron, pvx);  //create TriLeptonsCandidate with 2muons and 1electron
							DiMuonEle.setVtx( pvx );

							int ivtx = 0;
							for( unsigned int k = 0; k < primaryVertices->size() ; k++ ) {
								if( pvx == primaryVertices->ptrAt( k ) ) {
									ivtx = k;
									break;
								}	
							}
							DiMuonEle.setVertexIndex( ivtx );

							TriLeptonsColl->push_back( DiMuonEle );  // store the DiMuonEle into the collection
						}
					}

				}
			} 
		} 

		evt.put( TriLeptonsColl );

	}
}

typedef flashgg::TriLeptonsProducer FlashggTriLeptonsProducer;
DEFINE_FWK_MODULE( FlashggTriLeptonsProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4