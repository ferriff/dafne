#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "dafne/DataFormats/interface/MultiLeptonMultiJetCandidate.h"


using namespace edm;
using namespace std;

namespace flashgg {

	class MultiLeptonMultiJetProducer : public EDProducer
	{

	public:
		MultiLeptonMultiJetProducer( const ParameterSet & );

	private:
		void produce( Event &, const EventSetup & ) override;
		EDGetTokenT<View<Electron_t> > electronToken_;
		EDGetTokenT<View<Muon_t> > muonToken_;      
		EDGetTokenT<View<Jet_t> > jetToken_;  		
		EDGetTokenT<View<Track_t> > trackToken_;	
		EDGetTokenT<View<Vertex_t> > vertexToken_;

		double minElePt_;
		double maxEleEta_;
		double minMuPt_;
		double maxMuEta_;
		double minJetPt_;
		double maxJetEta_;
		double minTrackPt_;
		double maxTrackEta_;

	};

	MultiLeptonMultiJetProducer::MultiLeptonMultiJetProducer( const ParameterSet &iConfig ) :
		electronToken_( consumes<View<Electron_t> >( iConfig.getParameter<InputTag> ( "ElectronTag" ) ) ),
		muonToken_( consumes<View<Muon_t> >( iConfig.getParameter<InputTag> ( "MuonTag" ) ) ),
		jetToken_( consumes<View<Jet_t> >( iConfig.getParameter<InputTag> ( "JetTag" ) ) ),
		trackToken_( consumes<View<Track_t> >( iConfig.getParameter<InputTag> ( "TrackTag" ) ) ),
		vertexToken_( consumes<View<Vertex_t> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) )

	{
		minElePt_ = iConfig.getParameter<double>( "minElectronPt" );
		maxEleEta_ = iConfig.getParameter<double>( "maxElectronEta" );		
		minMuPt_ = iConfig.getParameter<double>( "minMuonPt" );
		maxMuEta_ = iConfig.getParameter<double>( "maxMuonEta" );
		minJetPt_ = iConfig.getParameter<double>( "minJetPt" );
		maxJetEta_ = iConfig.getParameter<double>( "maxJetEta" );
		minTrackPt_ = iConfig.getParameter<double>( "minTrackPt" );
		maxTrackEta_ = iConfig.getParameter<double>( "maxTrackEta" );		

		produces<vector<flashgg::MultiLeptonMultiJetCandidate> >();
	}

	void MultiLeptonMultiJetProducer::produce( Event &evt, const EventSetup & )
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

		Handle<View<Jet_t> > jets;
		evt.getByToken( jetToken_, jets );
		const vector<Jet_ptr > &jetPointers = jets->ptrs();

		Handle<View<Track_t> > tracks;
		evt.getByToken( trackToken_, tracks );
		const vector<Track_ptr > &trackPointers = tracks->ptrs();


		auto_ptr<vector<flashgg::MultiLeptonMultiJetCandidate> > multiLeptonMultiJetColl( new vector<flashgg::MultiLeptonMultiJetCandidate> );

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


					if (jetPointers.size() > 1) {
						for ( unsigned int l = 0 ; l < jetPointers.size() ; l++ ) {
							Jet_ptr jet1 = jetPointers[l];
							double pt_jet1 = jet1->pt();
							double eta_jet1 = jet1->eta();
							if( pt_jet1 < minJetPt_ || fabs( eta_jet1 ) > maxJetEta_ ) { continue; }
							for( unsigned int h = l + 1 ; h < jetPointers.size() ; h++ ) {		
								Jet_ptr jet2 = jetPointers[h];
								double pt_jet2 = jet2->pt();
								double eta_jet2 = jet2->eta();
								if( pt_jet2 < minJetPt_ || fabs( eta_jet2 ) > maxJetEta_ ) { continue; }

								MultiLeptonMultiJetCandidate DiEleDiJet( electron1, electron2, jet1, jet2, pvx);  //create MultiLeptonMultiJetCandidate with 2ele and 2jets

								DiEleDiJet.setVtx( pvx );  

								multiLeptonMultiJetColl->push_back( DiEleDiJet );  // store the DiEleDiJet into the collection
							}
						}
					}


					if (trackPointers.size() > 1) {
						for ( unsigned int l = 0 ; l < trackPointers.size() ; l++ ) {
							Track_ptr track1 = trackPointers[l];
							double pt_track1 = track1->pt();
							double eta_track1 = track1->eta();
							if( pt_track1 < minTrackPt_ || fabs( eta_track1 ) > maxTrackEta_ ) { continue; }
							for( unsigned int h = l + 1 ; h < trackPointers.size() ; h++ ) {		
								Track_ptr track2 = trackPointers[h];
								double pt_track2 = track2->pt();
								double eta_track2 = track2->eta();
								if( pt_track2 < minTrackPt_ || fabs( eta_track2 ) > maxTrackEta_ ) { continue; }

								MultiLeptonMultiJetCandidate DiEleDiTrack( electron1, electron2, track1, track2, pvx);  //create MultiLeptonMultiJetCandidate with 2ele and 2tracks

								DiEleDiTrack.setVtx( pvx );

								multiLeptonMultiJetColl->push_back( DiEleDiTrack );  // store the DiEleDiTrack into the collection
							}
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


					if (jetPointers.size() > 1) {
						for ( unsigned int l = 0 ; l < jetPointers.size() ; l++ ) {
							Jet_ptr jet1 = jetPointers[l];
							double pt_jet1 = jet1->pt();
							double eta_jet1 = jet1->eta();
							if( pt_jet1 < minJetPt_ || fabs( eta_jet1 ) > maxJetEta_ ) { continue; }
							for( unsigned int h = l + 1 ; h < jetPointers.size() ; h++ ) {		
								Jet_ptr jet2 = jetPointers[h];
								double pt_jet2 = jet2->pt();
								double eta_jet2 = jet2->eta();
								if( pt_jet2 < minJetPt_ || fabs( eta_jet2 ) > maxJetEta_ ) { continue; }

								MultiLeptonMultiJetCandidate DiMuonDiJet( muon1, muon2, jet1, jet2, pvx);  //create MultiLeptonMultiJetCandidate with 2muons and 2jets

								DiMuonDiJet.setVtx( pvx );

								multiLeptonMultiJetColl->push_back( DiMuonDiJet );  // store the DiMuonDiJet into the collection
							}
						}
					}


					if (trackPointers.size() > 1) {
						for ( unsigned int l = 0 ; l < trackPointers.size() ; l++ ) {
							Track_ptr track1 = trackPointers[l];
							double pt_track1 = track1->pt();
							double eta_track1 = track1->eta();
							if( pt_track1 < minTrackPt_ || fabs( eta_track1 ) > maxTrackEta_ ) { continue; }
							for( unsigned int h = l + 1 ; h < trackPointers.size() ; h++ ) {		
								Track_ptr track2 = trackPointers[h];
								double pt_track2 = track2->pt();
								double eta_track2 = track2->eta();
								if( pt_track2 < minTrackPt_ || fabs( eta_track2 ) > maxTrackEta_ ) { continue; }

								MultiLeptonMultiJetCandidate DiMuonDiTrack( muon1, muon2, track1, track2, pvx);  //create MultiLeptonMultiJetCandidate with 2muons and 2jets

								DiMuonDiTrack.setVtx( pvx );

								multiLeptonMultiJetColl->push_back( DiMuonDiTrack );  // store the DiMuonDiTrack into the collection
							}
						}
					}

				}
			} 
		} 



		if (electronPointers.size() > 1 && muonPointers.size() > 1 && jetPointers.size() > 1) {

			for( unsigned int i = 0 ; i < electronPointers.size() ; i++ ) {
				Electron_ptr electron = electronPointers[i];
				double pt_ele = electron->pt();
				double eta_ele = electron->eta();
				if( pt_ele < minElePt_ || fabs( eta_ele ) > maxEleEta_ ) { continue; } 

				for( unsigned int j = 0 ; j < muonPointers.size() ; j++ ) {
					Muon_ptr muon = muonPointers[j];
					double pt_mu = muon->pt();
					double eta_mu = muon->eta();
					if( pt_mu < minMuPt_ || fabs( eta_mu ) > maxMuEta_ ) { continue; }

					for ( unsigned int l = 0 ; l < jetPointers.size() ; l++ ) {
						Jet_ptr jet1 = jetPointers[l];
						double pt_jet1 = jet1->pt();
						double eta_jet1 = jet1->eta();
						if( pt_jet1 < minJetPt_ || fabs( eta_jet1 ) > maxJetEta_ ) { continue; }

						for( unsigned int h = l + 1 ; h < jetPointers.size() ; h++ ) {		
							Jet_ptr jet2 = jetPointers[h];
							double pt_jet2 = jet2->pt();
							double eta_jet2 = jet2->eta();
							if( pt_jet2 < minJetPt_ || fabs( eta_jet2 ) > maxJetEta_ ) { continue; }

							MultiLeptonMultiJetCandidate EleMuDiJet( electron, muon, jet1, jet2, pvx);  //create MultiLeptonMultiJetCandidate with 1ele, 1mu and 2jets

							EleMuDiJet.setVtx( pvx );  

							multiLeptonMultiJetColl->push_back( EleMuDiJet );  // store the EleMuDiJet into the collection
						}
					}
				}
			}
		}


		evt.put( multiLeptonMultiJetColl );

	}
}

typedef flashgg::MultiLeptonMultiJetProducer FlashggMultiLeptonMultiJetProducer;
DEFINE_FWK_MODULE( FlashggMultiLeptonMultiJetProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
