#include "dafne/Validation/interface/miniTreeMaker.h"


miniTreeMaker::miniTreeMaker( const ParameterSet &iConfig, TFileDirectory& fs, ConsumesCollector && cc ):
	BasicAnalyzer::BasicAnalyzer(iConfig, fs),
	genParticleToken_( cc.consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "genParticleTag" ) ) ),
	genInfoToken_(cc.consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "generatorInfo" ) ) ),
	PileUpToken_(cc.consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
	vertexToken_( cc.consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
	DiLeptonDiJetToken_( cc.consumes<View<flashgg::DiLeptonDiJetCandidate> >( iConfig.getParameter<InputTag> ( "DiLeptonDiJetTag" ) ) ),
	jetsToken_( cc.consumes<View<vector<flashgg::Jet> > >( iConfig.getParameter<InputTag> ( "JetsTag" ) ) ),
	genJetToken_( cc.consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
	electronToken_( cc.consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
	muonToken_( cc.consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
	triggerBitsToken_( cc.consumes<TriggerResults>( iConfig.getParameter<InputTag>( "triggerBits" ) ) ),
	rhoToken_(cc.consumes<double>(iConfig.getParameter <InputTag>("rhoFixedGridCollection" ) ) )
{
	bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
	lumiWeight_ = iConfig.getUntrackedParameter<double>( "lumiWeight", 1000. ); //pb                                                                                                                              
  
	globalVarsDumper_ = new GlobalVariablesDumper( iConfig.getParameter<ParameterSet>( "globalVariables" ), forward<ConsumesCollector>(cc) );
  
	eventTree = fs.make<TTree>( "event", "event" );
}



miniTreeMaker::~miniTreeMaker()
{
}



// ******************************************************************************************
void miniTreeMaker::beginJob()
{
	ngen = 0;
	ngenPre = 0;
	ndldj = 0;
	npre = 0;
	nEle = 0;
	nEleGood = 0;
	nElePassingHEEPid = 0;
	nmuons = 0;
	nmuonsGood = 0;

  // per-event tree
	eventTree->Branch( "run", &evInfo.run, "run/I" );
	eventTree->Branch( "event", &evInfo.event, "event/I" );
	eventTree->Branch( "lumi", &evInfo.lumi, "lumi/I" );
	eventTree->Branch( "rho", &evInfo.rho, "rho/F" );
	eventTree->Branch( "weight", &evInfo.weight, "weight/F" );
	eventTree->Branch( "puweight", &evInfo.puweight,"puweight/F");
	eventTree->Branch( "nvtx", &evInfo.nvtx, "nvtx/I" );
	eventTree->Branch( "npu", &evInfo.npu, "npu/I" );
	eventTree->Branch( "passEEJJhlt", &evInfo.passEEJJhlt, "passEEJJhlt/I" );
	eventTree->Branch( "passMMJJhlt", &evInfo.passMMJJhlt, "passMMJJhlt/I" );
	eventTree->Branch( "passEMJJhlt", &evInfo.passEMJJhlt, "passEMJJhlt/I" );
	eventTree->Branch( "passTandPEEhlt", &evInfo.passTandPEEhlt, "passTandPEEhlt/I" );
	eventTree->Branch( "passTandPMMhlt", &evInfo.passTandPMMhlt, "passTandPMMhlt/I" );

	eventTree->Branch( "vtx_x", &evInfo.vtx_x );
	eventTree->Branch( "vtx_y", &evInfo.vtx_y );
	eventTree->Branch( "vtx_z", &evInfo.vtx_z );

	eventTree->Branch( "ele_e", &evInfo.ele_e );
	eventTree->Branch( "ele_pt", &evInfo.ele_pt );
	eventTree->Branch( "ele_eta", &evInfo.ele_eta );
	eventTree->Branch( "ele_phi", &evInfo.ele_phi );
	eventTree->Branch( "ele_idmva", &evInfo.ele_idmva );
	eventTree->Branch( "ele_id", &evInfo.ele_id );
	eventTree->Branch( "ele_iso", &evInfo.ele_iso );
	eventTree->Branch( "ele_dz", &evInfo.ele_dz );
	eventTree->Branch( "ele_d0", &evInfo.ele_d0 );
	eventTree->Branch( "ele_passHEEPId", &evInfo.ele_passHEEPId );
	eventTree->Branch( "ele_passCutBasedEleId", &evInfo.ele_passCutBasedEleId );
	eventTree->Branch( "ele_isMatchedToGen", &evInfo.ele_isMatchedToGen );
	eventTree->Branch( "ele_charge", &evInfo.ele_charge );
	eventTree->Branch( "ele_etaSC", &evInfo.ele_etaSC );
	eventTree->Branch( "ele_isEcalDriven", &evInfo.ele_isEcalDriven );
	eventTree->Branch( "ele_dEtaIn", &evInfo.ele_dEtaIn );
	eventTree->Branch( "ele_dPhiIn", &evInfo.ele_dPhiIn );
	eventTree->Branch( "ele_hOverE", &evInfo.ele_hOverE );
	eventTree->Branch( "ele_full5x5_r9", &evInfo.ele_full5x5_r9 );
	eventTree->Branch( "ele_full5x5_sigmaIetaIeta", &evInfo.ele_full5x5_sigmaIetaIeta );
	eventTree->Branch( "ele_full5x5_E5x5", &evInfo.ele_full5x5_E5x5 );
	eventTree->Branch( "ele_full5x5_E1x5", &evInfo.ele_full5x5_E1x5 );
	eventTree->Branch( "ele_full5x5_E2x5", &evInfo.ele_full5x5_E2x5 );
	eventTree->Branch( "ele_full5x5_E2x5_Over_E5x5", &evInfo.ele_full5x5_E2x5_Over_E5x5 );
	eventTree->Branch( "ele_full5x5_E1x5_Over_E5x5", &evInfo.ele_full5x5_E1x5_Over_E5x5 );
	eventTree->Branch( "ele_EmHadDepth1Iso", &evInfo.ele_EmHadDepth1Iso );
	eventTree->Branch( "ele_ptTracksIso", &evInfo.ele_ptTracksIso );
	eventTree->Branch( "ele_innerLayerLostHits", &evInfo.ele_innerLayerLostHits );
	eventTree->Branch( "ele_dxy", &evInfo.ele_dxy );
	eventTree->Branch( "ele_eOverP", &evInfo.ele_eOverP );

	eventTree->Branch( "mu_e", &evInfo.mu_e );
	eventTree->Branch( "mu_pt", &evInfo.mu_pt );
	eventTree->Branch( "mu_eta", &evInfo.mu_eta );
	eventTree->Branch( "mu_phi", &evInfo.mu_phi );
	eventTree->Branch( "mu_iso", &evInfo.mu_iso );
	eventTree->Branch( "mu_isTight", &evInfo.mu_isTight );
	eventTree->Branch( "mu_isMedium", &evInfo.mu_isMedium );
	eventTree->Branch( "mu_isLoose", &evInfo.mu_isLoose );
	eventTree->Branch( "mu_isHighPt", &evInfo.mu_isHighPt );
	eventTree->Branch( "mu_isMatchedToGen", &evInfo.mu_isMatchedToGen );
	eventTree->Branch( "mu_charge", &evInfo.mu_charge );
	eventTree->Branch( "mu_dz", &evInfo.mu_dz );
	eventTree->Branch( "mu_dxy", &evInfo.mu_dxy );

	eventTree->Branch( "jet_e", &evInfo.jet_e );
	eventTree->Branch( "jet_pt", &evInfo.jet_pt );
	eventTree->Branch( "jet_eta", &evInfo.jet_eta );
	eventTree->Branch( "jet_phi", &evInfo.jet_phi );
	eventTree->Branch( "jet_bdiscriminant", &evInfo.jet_bdiscriminant );
	eventTree->Branch( "jet_partonFlavour", &evInfo.jet_partonFlavour );
	eventTree->Branch( "jet_hadronFlavour", &evInfo.jet_hadronFlavour );
	eventTree->Branch( "jet_isMatchedToGen", &evInfo.jet_isMatchedToGen );

	eventTree->Branch( "isEEJJ", &evInfo.isEEJJ ); 
	eventTree->Branch( "isEETT", &evInfo.isEETT ); 
	eventTree->Branch( "isMMJJ", &evInfo.isMMJJ ); 
	eventTree->Branch( "isMMTT", &evInfo.isMMTT ); 
	eventTree->Branch( "isEMJJ", &evInfo.isEMJJ ); 

	eventTree->Branch( "isSignalRegion", &evInfo.isSignalRegion );
	eventTree->Branch( "isLowMllCR", &evInfo.isLowMllCR );
	eventTree->Branch( "isLowMlljjCR", &evInfo.isLowMlljjCR );

	eventTree->Branch( "isBB", &evInfo.isBB );
	eventTree->Branch( "isEE", &evInfo.isEE );
	eventTree->Branch( "isEB", &evInfo.isEB );

	eventTree->Branch( "passPreselections", &evInfo.passPreselections );

	eventTree->Branch( "leadingLepton_e", &evInfo.leadingLepton_e ); 
	eventTree->Branch( "leadingLepton_pt", &evInfo.leadingLepton_pt ); 
	eventTree->Branch( "leadingLepton_eta", &evInfo.leadingLepton_eta ); 
	eventTree->Branch( "leadingLepton_phi", &evInfo.leadingLepton_phi ); 
	eventTree->Branch( "leadingLepton_charge", &evInfo.leadingLepton_charge ); 

	eventTree->Branch( "subLeadingLepton_e", &evInfo.subLeadingLepton_e ); 
	eventTree->Branch( "subLeadingLepton_pt", &evInfo.subLeadingLepton_pt ); 
	eventTree->Branch( "subLeadingLepton_eta", &evInfo.subLeadingLepton_eta ); 
	eventTree->Branch( "subLeadingLepton_phi", &evInfo.subLeadingLepton_phi ); 
	eventTree->Branch( "subLeadingLepton_charge", &evInfo.subLeadingLepton_charge ); 

	eventTree->Branch( "leadingJet_e", &evInfo.leadingJet_e );
	eventTree->Branch( "leadingJet_pt", &evInfo.leadingJet_pt ); 
	eventTree->Branch( "leadingJet_eta", &evInfo.leadingJet_eta );
	eventTree->Branch( "leadingJet_phi", &evInfo.leadingJet_phi ); 

	eventTree->Branch( "subLeadingJet_e", &evInfo.subLeadingJet_e ); 
	eventTree->Branch( "subLeadingJet_pt", &evInfo.subLeadingJet_pt ); 
	eventTree->Branch( "subLeadingJet_eta", &evInfo.subLeadingJet_eta ); 
	eventTree->Branch( "subLeadingJet_phi", &evInfo.subLeadingJet_phi ); 

	eventTree->Branch( "dRLeadLeptonLeadJet", &evInfo.dRLeadLeptonLeadJet );
	eventTree->Branch( "dRLeadLeptonSubLeadJet", &evInfo.dRLeadLeptonSubLeadJet );
	eventTree->Branch( "dRSubLeadLeptonLeadJet", &evInfo.dRSubLeadLeptonLeadJet );
	eventTree->Branch( "dRSubLeadLeptonSubLeadJet", &evInfo.dRSubLeadLeptonSubLeadJet );

	eventTree->Branch( "diLeptonDiJet_vtxIndex", &evInfo.diLeptonDiJet_vtxIndex );
	eventTree->Branch( "diLeptonDiJet_sumPt", &evInfo.diLeptonDiJet_sumPt );
	eventTree->Branch( "diLeptonDiJet_invMass", &evInfo.diLeptonDiJet_invMass ); 
	eventTree->Branch( "diLepton_invMass", &evInfo.diLepton_invMass ); 
	eventTree->Branch( "diLepton_pt", &evInfo.diLepton_pt ); 
	eventTree->Branch( "diJet_invMass", &evInfo.diJet_invMass ); 
	eventTree->Branch( "diJetLeadingLepton_invMass", &evInfo.diJetLeadingLepton_invMass ); 
	eventTree->Branch( "diJetSubLeadingLepton_invMass", &evInfo.diJetSubLeadingLepton_invMass ); 

	eventTree->Branch( "leadingEle_passHEEPId", &evInfo.leadingEle_passHEEPId );
	eventTree->Branch( "leadingEle_etaSC", &evInfo.leadingEle_etaSC );
	eventTree->Branch( "leadingEle_isEcalDriven", &evInfo.leadingEle_isEcalDriven );
	eventTree->Branch( "leadingEle_dEtaIn", &evInfo.leadingEle_dEtaIn );
	eventTree->Branch( "leadingEle_dPhiIn", &evInfo.leadingEle_dPhiIn );
	eventTree->Branch( "leadingEle_hOverE", &evInfo.leadingEle_hOverE );
	eventTree->Branch( "leadingEle_full5x5_r9", &evInfo.leadingEle_full5x5_r9 );
	eventTree->Branch( "leadingEle_full5x5_sigmaIetaIeta", &evInfo.leadingEle_full5x5_sigmaIetaIeta );
	eventTree->Branch( "leadingEle_full5x5_E5x5", &evInfo.leadingEle_full5x5_E5x5 );
	eventTree->Branch( "leadingEle_full5x5_E1x5", &evInfo.leadingEle_full5x5_E1x5 );
	eventTree->Branch( "leadingEle_full5x5_E2x5", &evInfo.leadingEle_full5x5_E2x5 );
	eventTree->Branch( "leadingEle_full5x5_E2x5_Over_E5x5", &evInfo.leadingEle_full5x5_E2x5_Over_E5x5 );
	eventTree->Branch( "leadingEle_full5x5_E1x5_Over_E5x5", &evInfo.leadingEle_full5x5_E1x5_Over_E5x5 );
	eventTree->Branch( "leadingEle_EmHadDepth1Iso", &evInfo.leadingEle_EmHadDepth1Iso );
	eventTree->Branch( "leadingEle_ptTracksIso", &evInfo.leadingEle_ptTracksIso );
	eventTree->Branch( "leadingEle_innerLayerLostHits", &evInfo.leadingEle_innerLayerLostHits );
	eventTree->Branch( "leadingEle_dxy", &evInfo.leadingEle_dxy );
	eventTree->Branch( "leadingEle_eOverP", &evInfo.leadingEle_eOverP );
	eventTree->Branch( "leadingEle_id", &evInfo.leadingEle_id );
	eventTree->Branch( "leadingEle_passCutBasedEleId", &evInfo.leadingEle_passCutBasedEleId );

	eventTree->Branch( "subLeadingEle_passHEEPId", &evInfo.subLeadingEle_passHEEPId );
	eventTree->Branch( "subLeadingEle_etaSC", &evInfo.subLeadingEle_etaSC );
	eventTree->Branch( "subLeadingEle_isEcalDriven", &evInfo.subLeadingEle_isEcalDriven );
	eventTree->Branch( "subLeadingEle_dEtaIn", &evInfo.subLeadingEle_dEtaIn );
	eventTree->Branch( "subLeadingEle_dPhiIn", &evInfo.subLeadingEle_dPhiIn );
	eventTree->Branch( "subLeadingEle_hOverE", &evInfo.subLeadingEle_hOverE );
	eventTree->Branch( "subLeadingEle_full5x5_r9", &evInfo.subLeadingEle_full5x5_r9 );
	eventTree->Branch( "subLeadingEle_full5x5_sigmaIetaIeta", &evInfo.subLeadingEle_full5x5_sigmaIetaIeta );
	eventTree->Branch( "subLeadingEle_full5x5_E5x5", &evInfo.subLeadingEle_full5x5_E5x5 );
	eventTree->Branch( "subLeadingEle_full5x5_E1x5", &evInfo.subLeadingEle_full5x5_E1x5 );
	eventTree->Branch( "subLeadingEle_full5x5_E2x5", &evInfo.subLeadingEle_full5x5_E2x5 );
	eventTree->Branch( "subLeadingEle_full5x5_E2x5_Over_E5x5", &evInfo.subLeadingEle_full5x5_E2x5_Over_E5x5 );
	eventTree->Branch( "subLeadingEle_full5x5_E1x5_Over_E5x5", &evInfo.subLeadingEle_full5x5_E1x5_Over_E5x5 );
	eventTree->Branch( "subLeadingEle_EmHadDepth1Iso", &evInfo.subLeadingEle_EmHadDepth1Iso );
	eventTree->Branch( "subLeadingEle_ptTracksIso", &evInfo.subLeadingEle_ptTracksIso );
	eventTree->Branch( "subLeadingEle_innerLayerLostHits", &evInfo.subLeadingEle_innerLayerLostHits );
	eventTree->Branch( "subLeadingEle_dxy", &evInfo.subLeadingEle_dxy );
	eventTree->Branch( "subLeadingEle_eOverP", &evInfo.subLeadingEle_eOverP );
	eventTree->Branch( "subLeadingEle_id", &evInfo.subLeadingEle_id );
	eventTree->Branch( "subLeadingEle_passCutBasedEleId", &evInfo.subLeadingEle_passCutBasedEleId );

	eventTree->Branch( "leadingMuon_isHighPt", &evInfo.leadingMuon_isHighPt );
	eventTree->Branch( "subLeadingMuon_isHighPt", &evInfo.subLeadingMuon_isHighPt );

	// cout << "inizialed branches" << endl;

}
// ******************************************************************************************




// ******************************************************************************************
// analyzer
//
void miniTreeMaker::analyze(const EventBase& evt)
{
	const Event *fullEvent = dynamic_cast<const Event *>(&evt);
	const Event &iEvent = (*fullEvent);  
  

	//-------------- access edm objects

	//--- only if MC
	Handle<View<reco::GenParticle> > genParticles;
	Handle<GenEventInfoProduct> genInfo;
	Handle<View< PileupSummaryInfo> > PileupInfos;
	Handle<View<reco::GenJet> > genJets;

	if ( !iEvent.isRealData() ) {
		// cout << "MC" << endl;
		iEvent.getByToken( genParticleToken_, genParticles );		
		iEvent.getByToken( genInfoToken_, genInfo );
		iEvent.getByToken( PileUpToken_, PileupInfos );
		iEvent.getByToken( genJetToken_, genJets );

		// for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		//    Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
		//    // cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
		// }    
		
		// for( unsigned int i = 0 ; i < genJets->size(); i++ ) {
		//    Ptr<reco::GenJet> genJet = genJets->ptrAt(i);
		// }  
	}

	Handle<View<reco::Vertex> > vertices;
	iEvent.getByToken( vertexToken_, vertices );
  
	Handle<View<flashgg::DiLeptonDiJetCandidate> > diLeptonDiJets;
	iEvent.getByToken( DiLeptonDiJetToken_, diLeptonDiJets );

	Handle<View<vector<flashgg::Jet> > > jets;
	iEvent.getByToken( jetsToken_, jets );
	// cout << "jets size " << jets->size() << endl;   

	Handle<View<flashgg::Electron> > electrons;
	iEvent.getByToken( electronToken_, electrons );
  
	Handle<View<flashgg::Muon> > muons;
	iEvent.getByToken( muonToken_, muons );
	
	Handle<TriggerResults> triggerBits;
	iEvent.getByToken( triggerBitsToken_, triggerBits );
  
	Handle<double> rhoHandle;
	iEvent.getByToken( rhoToken_, rhoHandle );
	double rho = *( rhoHandle.product() );




	//-------------- initialize tree
	initEventStructure();


	
	//-------------- check if event passes HLT
	const TriggerNames &triggerNames = iEvent.triggerNames( *triggerBits );

	bool passTrigger = false;

	for( unsigned index = 0; index < triggerNames.size(); ++index ) {

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_DoubleEle33_CaloIdL_MW") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passEEJJhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Mu50") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_TkMu50")) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passMMJJhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passEMJJhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Ele27_WPTight_Gsf") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passTandPEEhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_IsoMu24") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_IsoMu27") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passTandPMMhlt =  triggerBits->accept( index );
			passTrigger = true;
		}

	}


	globalVarsDumper_->fill( iEvent );

	evInfo.run = globalVarsDumper_->cache().run;
	evInfo.event = globalVarsDumper_->cache().event;		
	evInfo.lumi = globalVarsDumper_->cache().lumi;
	evInfo.rho = rho;

	// -- event weight (Lumi x cross section x gen weight)
	float w = 1.;
	// cout << "w = " << w << endl;
	if( ! iEvent.isRealData() ) {
		w = lumiWeight_;
		// cout << "lumiWeight = " << w << endl;
		if( genInfo.isValid() ) {
			const auto &weights = genInfo->weights();
			if( ! weights.empty() ) {
				w *= weights[0];
				// cout << "genWeight*lumiWeight = " << w << endl;
			}
		}

		if( globalVarsDumper_->puReWeight() ) {
			w *= globalVarsDumper_->cache().puweight;
		}

	}
	evInfo.weight = w;
	// cout << "weight = " << w << endl;


	// -- pileup weight
	if( globalVarsDumper_->puReWeight() ) {
		evInfo.puweight = globalVarsDumper_->cache().puweight;
	}
	// cout << "puweight = " << evInfo.puweight << endl;

	// -- number of reco vertices
	evInfo.nvtx = vertices->size() ;
	// cout << vertices->size() << " nvertices " << endl;

	// -- vertices
	for (UInt_t ivtx = 0 ; ivtx < vertices->size(); ivtx++){
		Ptr<reco::Vertex> vtx = vertices->ptrAt( ivtx );
		evInfo.vtx_x.push_back(vtx->x());
		evInfo.vtx_y.push_back(vtx->y());
		evInfo.vtx_z.push_back(vtx->z());
	}


	// -- number of pileup events
	float pu = 0.; 
	if( ! iEvent.isRealData() ) {
		for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
			Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
			if( pu_bunchcrossing == 0 ) {
				pu = PileupInfos->ptrAt( PVI )->getPU_NumInteractions();
			}
		}
	}
	evInfo.npu = pu;
	// cout << "pu " << pu << endl;


	// -- electrons
	for (UInt_t iele = 0 ; iele < electrons->size(); iele++){
		// cout << "enter electron loop" << endl;
		nEle++;

		Ptr<flashgg::Electron> electron = electrons->ptrAt( iele );
		if (fabs(electron->eta()) > 2.4) { continue; }
		// if (electron->pt() < electronPtThreshold_) { continue; }  
		if( electron->hasMatchedConversion() ) { continue; } // remove conversions
		nEleGood++;

		float isol = electronIsolation(electron, rho); 

		Ptr<reco::Vertex> ele_vtx = chooseElectronVertex( electron,  vertices->ptrs() );
		float dz = electron->gsfTrack()->dz( ele_vtx->position() );
		float d0 = electron->gsfTrack()->dxy( ele_vtx->position() ); 
	
		int passHEEPId = passHEEPIdCuts( electron, vertices->ptrs(), rho );
		if (passHEEPId) nElePassingHEEPid++;

		int passCutBasedId = passCutBasedEleId( electron, rho );

		int mcMatch = -1;
		if( ! iEvent.isRealData() ) mcMatch = electronMatchingToGen(electron, genParticles); 

		Ptr<reco::Vertex> best_vtx_ele = chooseBestVtx(vertices->ptrs(), electron);

		evInfo.ele_e.push_back(electron->energy());
		evInfo.ele_pt.push_back(electron->pt());
		evInfo.ele_eta.push_back(electron->eta());
		evInfo.ele_phi.push_back(electron->phi());
		evInfo.ele_idmva.push_back(electron->nonTrigMVA());
		evInfo.ele_id.push_back( getEleId(electron) );  
		evInfo.ele_iso.push_back(isol);
		evInfo.ele_dz.push_back(dz);
		evInfo.ele_d0.push_back(d0);
		evInfo.ele_passHEEPId.push_back(passHEEPId);
		evInfo.ele_passCutBasedEleId.push_back(passCutBasedId);
		evInfo.ele_isMatchedToGen.push_back(mcMatch);
		evInfo.ele_charge.push_back(electron->charge());
		evInfo.ele_etaSC.push_back(electron->superCluster()->eta());
		evInfo.ele_isEcalDriven.push_back(electron->ecalDrivenSeed());
		evInfo.ele_dEtaIn.push_back(electron->deltaEtaSuperClusterTrackAtVtx());  
		evInfo.ele_dPhiIn.push_back(electron->deltaPhiSuperClusterTrackAtVtx());
		evInfo.ele_hOverE.push_back(electron->hadronicOverEm());
		evInfo.ele_full5x5_r9.push_back(electron->full5x5_r9());  
		evInfo.ele_full5x5_sigmaIetaIeta.push_back(electron->full5x5_sigmaIetaIeta());
		evInfo.ele_full5x5_E5x5.push_back(electron->full5x5_e5x5());
		evInfo.ele_full5x5_E1x5.push_back(electron->full5x5_e1x5());
		evInfo.ele_full5x5_E2x5.push_back(electron->full5x5_e2x5Max());
		evInfo.ele_full5x5_E2x5_Over_E5x5.push_back((electron->full5x5_e2x5Max()) / (electron->full5x5_e5x5()));
		evInfo.ele_full5x5_E1x5_Over_E5x5.push_back((electron->full5x5_e1x5()) / (electron->full5x5_e5x5()));
		evInfo.ele_EmHadDepth1Iso.push_back(electron->dr03EcalRecHitSumEt()+electron->dr03HcalDepth1TowerSumEt());
		evInfo.ele_ptTracksIso.push_back(electron->dr03TkSumPt());  //TO MODIFY
		evInfo.ele_innerLayerLostHits.push_back(electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS));
		evInfo.ele_dxy.push_back(fabs(electron->gsfTrack()->dxy( best_vtx_ele->position())));
		evInfo.ele_eOverP.push_back(electron->eSuperClusterOverP());
	}       


	// -- muons
	for (UInt_t imu = 0 ; imu < muons->size(); imu++){
		// cout << "enter muon loop" << endl;
		nmuons++;

		Ptr<flashgg::Muon> muon = muons->ptrAt( imu );
		if (fabs(muon->eta()) > 2.4) { continue; }
		// if (muon->pt()  < muonPtThreshold_) { continue; }  
		nmuonsGood++;

		// muon ID and isolation: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
		float muPFCombRelIso = ( muon->pfIsolationR04().sumChargedHadronPt + max( 0.,muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - 0.5 * muon->pfIsolationR04().sumPUPt ) ) / ( muon->pt() );

		Ptr<reco::Vertex> muonVtx = chooseBestMuonVtx(vertices->ptrs(), muon);

		float dz = -999;
		float dxy = -999;
		if ( !(!muon->innerTrack()) ) {
			dz = muon->innerTrack()->dz( muonVtx->position() );
			dxy = muon->innerTrack()->dxy( muonVtx->position() ); 			
		}

		int mcMatch =  -1;
		if( ! iEvent.isRealData() ) mcMatch = muonMatchingToGen(muon, genParticles); 

		evInfo.mu_e.push_back(muon->energy());
		evInfo.mu_pt.push_back(muon->pt());
		evInfo.mu_eta.push_back(muon->eta());
		evInfo.mu_phi.push_back(muon->phi());
		evInfo.mu_iso.push_back(muPFCombRelIso);
		evInfo.mu_isTight.push_back(muon->isTightMuon( *muonVtx ));
		evInfo.mu_isMedium.push_back(muon->isMediumMuon( ));
		evInfo.mu_isLoose.push_back(muon->isLooseMuon( ));
		evInfo.mu_isHighPt.push_back(muon->isHighPtMuon( *muonVtx ));
		evInfo.mu_isMatchedToGen.push_back(mcMatch); 
		evInfo.mu_charge.push_back(muon->charge());
		evInfo.mu_dz.push_back( dz );
		evInfo.mu_dxy.push_back( fabs(dxy) );
	}


	// -- jets
	// for (UInt_t ijetVector = 0 ; ijetVector < jets->size(); ijetVector++){
		// Ptr<vector<flashgg::Jet> > jetVector = jets->ptrAt( ijetVector );
	if (jets->size() > 0) {
		Ptr<vector<flashgg::Jet> > jetVector = jets->ptrAt( 0 );
		// if (jetVector->size() > 0) cout << "jetVector size " << jetVector->size() << endl;

		for (UInt_t ijet = 0 ; ijet < jetVector->size(); ijet++){
			// cout << "enter jet loop" << endl;
			flashgg::Jet jet = jetVector->at( ijet );

			int isMatchedToGen =  -1;
			if( ! iEvent.isRealData() ) isMatchedToGen = jetMatchingToGen(jet, genJets); 

			evInfo.jet_e.push_back(jet.energy());
			evInfo.jet_pt.push_back(jet.pt());
			evInfo.jet_eta.push_back(jet.eta());
			evInfo.jet_phi.push_back(jet.phi());
			evInfo.jet_bdiscriminant.push_back(jet.bDiscriminator( bTag_ ));
			evInfo.jet_hadronFlavour.push_back(jet.hadronFlavour());
			evInfo.jet_partonFlavour.push_back(jet.partonFlavour());
			evInfo.jet_isMatchedToGen.push_back(isMatchedToGen);
		}
	}


	// -- diLeptonDiJets
	// cout << "dldj size " << diLeptonDiJets->size() << endl;

	for ( unsigned int idldj = 0; idldj < diLeptonDiJets->size(); idldj++){
		// cout << "enter dldj loop" << endl;
		Ptr<flashgg::DiLeptonDiJetCandidate> diLeptonDiJet = diLeptonDiJets->ptrAt( idldj );        

		if (! iEvent.isRealData() ) { 
			ngen++;
			if ( passDiLeptonDiJetPreselection(diLeptonDiJet) ) ngenPre++; 
		} else {
			ndldj++;
			if ( passDiLeptonDiJetPreselection(diLeptonDiJet) ) npre++; 
		}

	
		if (diLeptonDiJet->leadingLeptonPt() < 35 || diLeptonDiJet->subLeadingLeptonPt() < 35) continue;
		if (fabs(diLeptonDiJet->leadingLeptonEta()) > 2.4 || fabs(diLeptonDiJet->subLeadingLeptonEta()) > 2.4) continue;

		if (diLeptonDiJet->leadingJet()->pt() < 40 || diLeptonDiJet->subLeadingJet()->pt() < 40) continue;
		if (fabs(diLeptonDiJet->leadingJet()->eta()) > 2.4 || fabs(diLeptonDiJet->subLeadingJet()->eta()) > 2.4) continue;


		float leadingLeptonCharge    = -999.;
		float subLeadingLeptonCharge = -999.;
		float leadingLeptonE         = -999.;
		float subLeadingLeptonE      = -999.;

		bool leadElePassHEEPId         = false;
		float leadingEleEtaSC          = -999.;
		bool leadingEleIsEcalDriven    = false;
		float leadingEleDEtaIn         = -999.;
		float leadingEleDPhiIn         = -999.;
		float leadingEleHoE            = -999.;
		float leadingEleR9             = -999.;
		float leadingEleSigmaIetaIeta  = -999.;
		float leadingEleE5x5           = -999.;
		float leadingEleE1x5           = -999.; 
		float leadingEleE2x5           = -999.;
		float leadingEleE2x5_E5x5      = -999.;
		float leadingEleE1x5_E5x5      = -999.;
		float leadingEleEmHadDepth1Iso = -999.;
		float leadingElePtTracksIso    = -999.;
		float leadingEleMissingHits    = -999.;
		float leadingEleDxy            = -999.;
		float leadingEleEoverP         = -999.;
		bool leadingElePassCutBasedId  = false;
		unsigned leadingEleId          = 0;

		bool subLeadElePassHEEPId         = false;
		float subLeadingEleEtaSC          = -999.;
		bool subLeadingEleIsEcalDriven    = false;
		float subLeadingEleDEtaIn         = -999.;
		float subLeadingEleDPhiIn         = -999.;
		float subLeadingEleHoE            = -999.;
		float subLeadingEleR9             = -999.;
		float subLeadingEleSigmaIetaIeta  = -999.;
		float subLeadingEleE5x5           = -999.;
		float subLeadingEleE1x5           = -999.; 
		float subLeadingEleE2x5           = -999.;
		float subLeadingEleE2x5_E5x5      = -999.;
		float subLeadingEleE1x5_E5x5      = -999.;
		float subLeadingEleEmHadDepth1Iso = -999.;
		float subLeadingElePtTracksIso    = -999.;
		float subLeadingEleMissingHits    = -999.;
		float subLeadingEleDxy            = -999.;
		float subLeadingEleEoverP         = -999.;
		bool subLeadingElePassCutBasedId  = false;
		unsigned subLeadingEleId          = 0;

		bool leadingMuonIsHighPt = false;
		bool subLeadingMuonIsHighPt = false;


		if (diLeptonDiJet->isEEJJ() || diLeptonDiJet->isEETT()) {
			leadingLeptonCharge = diLeptonDiJet->leadingEle()->charge();
			subLeadingLeptonCharge = diLeptonDiJet->subLeadingEle()->charge();
			leadingLeptonE = diLeptonDiJet->leadingEle()->energy();
			subLeadingLeptonE = diLeptonDiJet->subLeadingEle()->energy();

			Ptr<reco::Vertex> best_vtx_leadEle = chooseBestVtx(vertices->ptrs(), diLeptonDiJet->leadingEle());
			Ptr<reco::Vertex> best_vtx_subLeadEle = chooseBestVtx(vertices->ptrs(), diLeptonDiJet->subLeadingEle());

			leadElePassHEEPId = passHEEPIdCuts( diLeptonDiJet->leadingEle(), vertices->ptrs(), rho );
			leadingEleEtaSC = diLeptonDiJet->leadingEle()->superCluster()->eta();
			leadingEleIsEcalDriven = diLeptonDiJet->leadingEle()->ecalDrivenSeed();
			leadingEleDEtaIn = diLeptonDiJet->leadingEle()->deltaEtaSuperClusterTrackAtVtx();  
			leadingEleDPhiIn = diLeptonDiJet->leadingEle()->deltaPhiSuperClusterTrackAtVtx();
			leadingEleHoE = diLeptonDiJet->leadingEle()->hadronicOverEm();
			leadingEleR9 = diLeptonDiJet->leadingEle()->full5x5_r9();
			leadingEleSigmaIetaIeta = diLeptonDiJet->leadingEle()->full5x5_sigmaIetaIeta();
			leadingEleE5x5 = diLeptonDiJet->leadingEle()->full5x5_e5x5();
			leadingEleE1x5 = diLeptonDiJet->leadingEle()->full5x5_e1x5();
			leadingEleE2x5 = diLeptonDiJet->leadingEle()->full5x5_e2x5Max();
			leadingEleE2x5_E5x5 = leadingEleE2x5 / leadingEleE5x5;
			leadingEleE1x5_E5x5 = leadingEleE1x5 / leadingEleE5x5;
			leadingEleEmHadDepth1Iso = diLeptonDiJet->leadingEle()->dr03EcalRecHitSumEt() + diLeptonDiJet->leadingEle()->dr03HcalDepth1TowerSumEt();
			leadingElePtTracksIso = diLeptonDiJet->leadingEle()->dr03TkSumPt();  //TO MODIFY
			leadingEleMissingHits = diLeptonDiJet->leadingEle()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
			leadingEleDxy = fabs(diLeptonDiJet->leadingEle()->gsfTrack()->dxy( best_vtx_leadEle->position()));
			leadingEleEoverP = diLeptonDiJet->leadingEle()->eSuperClusterOverP();
			leadingElePassCutBasedId = passCutBasedEleId( diLeptonDiJet->leadingEle(), rho );
			leadingEleId = getEleId(diLeptonDiJet->leadingEle());

			subLeadElePassHEEPId = passHEEPIdCuts( diLeptonDiJet->subLeadingEle(), vertices->ptrs(), rho );
			subLeadingEleEtaSC = diLeptonDiJet->subLeadingEle()->superCluster()->eta();
			subLeadingEleIsEcalDriven = diLeptonDiJet->subLeadingEle()->ecalDrivenSeed();
			subLeadingEleDEtaIn = diLeptonDiJet->subLeadingEle()->deltaEtaSuperClusterTrackAtVtx();  
			subLeadingEleDPhiIn = diLeptonDiJet->subLeadingEle()->deltaPhiSuperClusterTrackAtVtx();
			subLeadingEleHoE = diLeptonDiJet->subLeadingEle()->hadronicOverEm();
			subLeadingEleR9 = diLeptonDiJet->subLeadingEle()->full5x5_r9();
			subLeadingEleSigmaIetaIeta = diLeptonDiJet->subLeadingEle()->full5x5_sigmaIetaIeta();
			subLeadingEleE5x5 = diLeptonDiJet->subLeadingEle()->full5x5_e5x5();
			subLeadingEleE1x5 = diLeptonDiJet->subLeadingEle()->full5x5_e1x5();
			subLeadingEleE2x5 = diLeptonDiJet->subLeadingEle()->full5x5_e2x5Max();
			subLeadingEleE2x5_E5x5 = subLeadingEleE2x5 / subLeadingEleE5x5;
			subLeadingEleE1x5_E5x5 = subLeadingEleE1x5 / subLeadingEleE5x5;
			subLeadingEleEmHadDepth1Iso = diLeptonDiJet->subLeadingEle()->dr03EcalRecHitSumEt() + diLeptonDiJet->subLeadingEle()->dr03HcalDepth1TowerSumEt();
			subLeadingElePtTracksIso = diLeptonDiJet->subLeadingEle()->dr03TkSumPt();  //TO MODIFY
			subLeadingEleMissingHits = diLeptonDiJet->subLeadingEle()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
			subLeadingEleDxy = fabs(diLeptonDiJet->subLeadingEle()->gsfTrack()->dxy( best_vtx_subLeadEle->position()));
			subLeadingEleEoverP = diLeptonDiJet->subLeadingEle()->eSuperClusterOverP();
			subLeadingElePassCutBasedId = passCutBasedEleId( diLeptonDiJet->subLeadingEle(), rho );
			subLeadingEleId = getEleId(diLeptonDiJet->subLeadingEle());

		} else if (diLeptonDiJet->isMMJJ() || diLeptonDiJet->isMMTT()) {
			leadingLeptonCharge = diLeptonDiJet->leadingMuon()->charge();
			subLeadingLeptonCharge = diLeptonDiJet->subLeadingMuon()->charge();
			leadingLeptonE = diLeptonDiJet->leadingMuon()->energy();
			subLeadingLeptonE = diLeptonDiJet->subLeadingMuon()->energy();

			Ptr<reco::Vertex> leadingMuonVtx = chooseBestMuonVtx(vertices->ptrs(), diLeptonDiJet->leadingMuon());
			Ptr<reco::Vertex> subLeadingMuonVtx = chooseBestMuonVtx(vertices->ptrs(), diLeptonDiJet->subLeadingMuon());

			leadingMuonIsHighPt = diLeptonDiJet->leadingMuon()->isHighPtMuon( *leadingMuonVtx );
			subLeadingMuonIsHighPt = diLeptonDiJet->subLeadingMuon()->isHighPtMuon( *subLeadingMuonVtx );

		} else if (diLeptonDiJet->isEMJJ()) {
			if (diLeptonDiJet->electron()->pt() > diLeptonDiJet->muon()->pt()) {
				leadingLeptonCharge = diLeptonDiJet->electron()->charge();
				subLeadingLeptonCharge = diLeptonDiJet->muon()->charge();
				leadingLeptonE = diLeptonDiJet->electron()->energy();
				subLeadingLeptonE = diLeptonDiJet->muon()->energy();

				Ptr<reco::Vertex> best_vtx_leadEle = chooseBestVtx(vertices->ptrs(), diLeptonDiJet->electron());

				leadElePassHEEPId = passHEEPIdCuts( diLeptonDiJet->electron(), vertices->ptrs(), rho ); 
				leadingEleEtaSC = diLeptonDiJet->electron()->superCluster()->eta();
				leadingEleIsEcalDriven = diLeptonDiJet->electron()->ecalDrivenSeed();
				leadingEleDEtaIn = diLeptonDiJet->electron()->deltaEtaSuperClusterTrackAtVtx();  
				leadingEleDPhiIn = diLeptonDiJet->electron()->deltaPhiSuperClusterTrackAtVtx();
				leadingEleHoE = diLeptonDiJet->electron()->hadronicOverEm();
				leadingEleR9 = diLeptonDiJet->electron()->full5x5_r9();
				leadingEleSigmaIetaIeta = diLeptonDiJet->electron()->full5x5_sigmaIetaIeta();
				leadingEleE5x5 = diLeptonDiJet->electron()->full5x5_e5x5();
				leadingEleE1x5 = diLeptonDiJet->electron()->full5x5_e1x5();
				leadingEleE2x5 = diLeptonDiJet->electron()->full5x5_e2x5Max();
				leadingEleE2x5_E5x5 = leadingEleE2x5 / leadingEleE5x5;
				leadingEleE1x5_E5x5 = leadingEleE1x5 / leadingEleE5x5;
				leadingEleEmHadDepth1Iso = diLeptonDiJet->electron()->dr03EcalRecHitSumEt() + diLeptonDiJet->electron()->dr03HcalDepth1TowerSumEt();
				leadingElePtTracksIso = diLeptonDiJet->electron()->dr03TkSumPt();  //TO MODIFY
				leadingEleMissingHits = diLeptonDiJet->electron()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
				leadingEleDxy = fabs(diLeptonDiJet->electron()->gsfTrack()->dxy( best_vtx_leadEle->position()));
				leadingEleEoverP = diLeptonDiJet->electron()->eSuperClusterOverP();
				leadingElePassCutBasedId = passCutBasedEleId( diLeptonDiJet->electron(), rho );
				leadingEleId = getEleId(diLeptonDiJet->electron());

				Ptr<reco::Vertex> subLeadingMuonVtx = chooseBestMuonVtx(vertices->ptrs(), diLeptonDiJet->muon());
				subLeadingMuonIsHighPt = diLeptonDiJet->muon()->isHighPtMuon( *subLeadingMuonVtx );

			} else {
				leadingLeptonCharge = diLeptonDiJet->muon()->charge();
				subLeadingLeptonCharge = diLeptonDiJet->electron()->charge();		
				leadingLeptonE = diLeptonDiJet->muon()->energy();
				subLeadingLeptonE = diLeptonDiJet->electron()->energy();
				
				Ptr<reco::Vertex> best_vtx_subLeadEle = chooseBestVtx(vertices->ptrs(), diLeptonDiJet->electron());

				subLeadElePassHEEPId = passHEEPIdCuts( diLeptonDiJet->electron(), vertices->ptrs(), rho );
				subLeadingEleEtaSC = diLeptonDiJet->electron()->superCluster()->eta();
				subLeadingEleIsEcalDriven = diLeptonDiJet->electron()->ecalDrivenSeed();
				subLeadingEleDEtaIn = diLeptonDiJet->electron()->deltaEtaSuperClusterTrackAtVtx();  
				subLeadingEleDPhiIn = diLeptonDiJet->electron()->deltaPhiSuperClusterTrackAtVtx();
				subLeadingEleHoE = diLeptonDiJet->electron()->hadronicOverEm();
				subLeadingEleR9 = diLeptonDiJet->electron()->full5x5_r9();
				subLeadingEleSigmaIetaIeta = diLeptonDiJet->electron()->full5x5_sigmaIetaIeta();
				subLeadingEleE5x5 = diLeptonDiJet->electron()->full5x5_e5x5();
				subLeadingEleE1x5 = diLeptonDiJet->electron()->full5x5_e1x5();
				subLeadingEleE2x5 = diLeptonDiJet->electron()->full5x5_e2x5Max();
				subLeadingEleE2x5_E5x5 = subLeadingEleE2x5 / subLeadingEleE5x5;
				subLeadingEleE1x5_E5x5 = subLeadingEleE1x5 / subLeadingEleE5x5;
				subLeadingEleEmHadDepth1Iso = diLeptonDiJet->electron()->dr03EcalRecHitSumEt() + diLeptonDiJet->electron()->dr03HcalDepth1TowerSumEt();
				subLeadingElePtTracksIso = diLeptonDiJet->electron()->dr03TkSumPt();  //TO MODIFY
				subLeadingEleMissingHits = diLeptonDiJet->electron()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
				subLeadingEleDxy = fabs(diLeptonDiJet->electron()->gsfTrack()->dxy( best_vtx_subLeadEle->position()));
				subLeadingEleEoverP = diLeptonDiJet->electron()->eSuperClusterOverP();
				subLeadingElePassCutBasedId = passCutBasedEleId( diLeptonDiJet->electron(), rho );
				subLeadingEleId = getEleId(diLeptonDiJet->electron());

				Ptr<reco::Vertex> leadingMuonVtx = chooseBestMuonVtx(vertices->ptrs(), diLeptonDiJet->muon());
				leadingMuonIsHighPt = diLeptonDiJet->muon()->isHighPtMuon( *leadingMuonVtx );
			}
		}

		// if (diLeptonDiJet->isEEJJ()) cout << "isEEJJ" << endl;
		// if (diLeptonDiJet->isMMJJ()) cout << "isMMJJ" << endl;
		// if (diLeptonDiJet->isEETT()) cout << "isEETT" << endl;
		// if (diLeptonDiJet->isMMTT()) cout << "isMMTT" << endl;
		// if (diLeptonDiJet->isEMJJ()) cout << "isEMJJ" << endl; 
		// if (!(diLeptonDiJet->isEMJJ())) cout << "no EMJJ" << endl;

		evInfo.isEEJJ.push_back(diLeptonDiJet->isEEJJ());  
		evInfo.isEETT.push_back(diLeptonDiJet->isEETT());
		evInfo.isMMJJ.push_back(diLeptonDiJet->isMMJJ());
		evInfo.isMMTT.push_back(diLeptonDiJet->isMMTT());
		evInfo.isEMJJ.push_back(diLeptonDiJet->isEMJJ());  

		evInfo.isSignalRegion.push_back( isSignalRegion(diLeptonDiJet) );
		evInfo.isLowMllCR.push_back( isLowMllCR(diLeptonDiJet) );
		evInfo.isLowMlljjCR.push_back( isLowMlljjCR(diLeptonDiJet) ); 

		evInfo.isBB.push_back( isBB(diLeptonDiJet) );
		evInfo.isEE.push_back( isEE(diLeptonDiJet) );
		evInfo.isEB.push_back( isEB(diLeptonDiJet) );

		evInfo.passPreselections.push_back( passDiLeptonDiJetPreselection(diLeptonDiJet) );

		evInfo.leadingLepton_e.push_back(leadingLeptonE);
		evInfo.leadingLepton_pt.push_back(diLeptonDiJet->leadingLeptonPt());
		evInfo.leadingLepton_eta.push_back(diLeptonDiJet->leadingLeptonEta());
		evInfo.leadingLepton_phi.push_back(diLeptonDiJet->leadingLeptonPhi());
		evInfo.leadingLepton_charge.push_back(leadingLeptonCharge);

		evInfo.subLeadingLepton_e.push_back(subLeadingLeptonE);
		evInfo.subLeadingLepton_pt.push_back(diLeptonDiJet->subLeadingLeptonPt());
		evInfo.subLeadingLepton_eta.push_back(diLeptonDiJet->subLeadingLeptonEta());
		evInfo.subLeadingLepton_phi.push_back(diLeptonDiJet->subLeadingLeptonPhi());
		evInfo.subLeadingLepton_charge.push_back(subLeadingLeptonCharge);

		evInfo.leadingJet_e.push_back(diLeptonDiJet->leadingJet()->energy());
		evInfo.leadingJet_pt.push_back(diLeptonDiJet->leadingJet()->pt());
		evInfo.leadingJet_eta.push_back(diLeptonDiJet->leadingJet()->eta());
		evInfo.leadingJet_phi.push_back(diLeptonDiJet->leadingJet()->phi());

		evInfo.subLeadingJet_e.push_back(diLeptonDiJet->subLeadingJet()->energy());
		evInfo.subLeadingJet_pt.push_back(diLeptonDiJet->subLeadingJet()->pt());
		evInfo.subLeadingJet_eta.push_back(diLeptonDiJet->subLeadingJet()->eta());
		evInfo.subLeadingJet_phi.push_back(diLeptonDiJet->subLeadingJet()->phi());

		evInfo.dRLeadLeptonLeadJet.push_back(deltaR(diLeptonDiJet->leadingLeptonEta(), diLeptonDiJet->leadingLeptonPhi(), diLeptonDiJet->leadingJet()->eta(), diLeptonDiJet->leadingJet()->phi()));
		evInfo.dRLeadLeptonSubLeadJet.push_back(deltaR(diLeptonDiJet->leadingLeptonEta(), diLeptonDiJet->leadingLeptonPhi(), diLeptonDiJet->subLeadingJet()->eta(), diLeptonDiJet->subLeadingJet()->phi()));
		evInfo.dRSubLeadLeptonLeadJet.push_back(deltaR(diLeptonDiJet->subLeadingLeptonEta(), diLeptonDiJet->subLeadingLeptonPhi(), diLeptonDiJet->leadingJet()->eta(), diLeptonDiJet->leadingJet()->phi()));
		evInfo.dRSubLeadLeptonSubLeadJet.push_back(deltaR(diLeptonDiJet->subLeadingLeptonEta(), diLeptonDiJet->subLeadingLeptonPhi(), diLeptonDiJet->subLeadingJet()->eta(), diLeptonDiJet->subLeadingJet()->phi()));

		evInfo.diLeptonDiJet_vtxIndex.push_back(diLeptonDiJet->vertexIndex()); 
		evInfo.diLeptonDiJet_sumPt.push_back(diLeptonDiJet->sumPt());
		evInfo.diLeptonDiJet_invMass.push_back(diLeptonDiJet->invMass());
		evInfo.diLepton_invMass.push_back(diLeptonDiJet->diLeptonInvMass());
		evInfo.diLepton_pt.push_back( diLeptonPt(diLeptonDiJet) );
		evInfo.diJet_invMass.push_back( diJetInvMass(diLeptonDiJet) );
		evInfo.diJetLeadingLepton_invMass.push_back( diJetLeadingLeptonInvMass(diLeptonDiJet) );
		evInfo.diJetSubLeadingLepton_invMass.push_back( diJetSubLeadingLeptonInvMass(diLeptonDiJet) );

		evInfo.leadingEle_passHEEPId.push_back(leadElePassHEEPId);
		evInfo.leadingEle_etaSC.push_back(leadingEleEtaSC);
		evInfo.leadingEle_isEcalDriven.push_back(leadingEleIsEcalDriven);
		evInfo.leadingEle_dEtaIn.push_back(leadingEleDEtaIn);
		evInfo.leadingEle_dPhiIn.push_back(leadingEleDPhiIn);
		evInfo.leadingEle_hOverE.push_back(leadingEleHoE);
		evInfo.leadingEle_full5x5_r9.push_back(leadingEleR9); 
		evInfo.leadingEle_full5x5_sigmaIetaIeta.push_back(leadingEleSigmaIetaIeta);
		evInfo.leadingEle_full5x5_E5x5.push_back(leadingEleE5x5);
		evInfo.leadingEle_full5x5_E1x5.push_back(leadingEleE1x5);
		evInfo.leadingEle_full5x5_E2x5.push_back(leadingEleE2x5);
		evInfo.leadingEle_full5x5_E2x5_Over_E5x5.push_back(leadingEleE2x5_E5x5);
		evInfo.leadingEle_full5x5_E1x5_Over_E5x5.push_back(leadingEleE1x5_E5x5);
		evInfo.leadingEle_EmHadDepth1Iso.push_back(leadingEleEmHadDepth1Iso);
		evInfo.leadingEle_ptTracksIso.push_back(leadingElePtTracksIso);
		evInfo.leadingEle_innerLayerLostHits.push_back(leadingEleMissingHits);
		evInfo.leadingEle_dxy.push_back(leadingEleDxy);
		evInfo.leadingEle_eOverP.push_back(leadingEleEoverP);
		evInfo.leadingEle_id.push_back(leadingEleId);
		evInfo.leadingEle_passCutBasedEleId.push_back(leadingElePassCutBasedId);

		evInfo.subLeadingEle_passHEEPId.push_back(subLeadElePassHEEPId);		
		evInfo.subLeadingEle_etaSC.push_back(subLeadingEleEtaSC);		
		evInfo.subLeadingEle_isEcalDriven.push_back(subLeadingEleIsEcalDriven);
		evInfo.subLeadingEle_dEtaIn.push_back(subLeadingEleDEtaIn);
		evInfo.subLeadingEle_dPhiIn.push_back(subLeadingEleDPhiIn);
		evInfo.subLeadingEle_hOverE.push_back(subLeadingEleHoE);
		evInfo.subLeadingEle_full5x5_r9.push_back(subLeadingEleR9); 
		evInfo.subLeadingEle_full5x5_sigmaIetaIeta.push_back(subLeadingEleSigmaIetaIeta);
		evInfo.subLeadingEle_full5x5_E5x5.push_back(subLeadingEleE5x5);
		evInfo.subLeadingEle_full5x5_E1x5.push_back(subLeadingEleE1x5);
		evInfo.subLeadingEle_full5x5_E2x5.push_back(subLeadingEleE2x5);
		evInfo.subLeadingEle_full5x5_E2x5_Over_E5x5.push_back(subLeadingEleE2x5_E5x5);
		evInfo.subLeadingEle_full5x5_E1x5_Over_E5x5.push_back(subLeadingEleE1x5_E5x5);
		evInfo.subLeadingEle_EmHadDepth1Iso.push_back(subLeadingEleEmHadDepth1Iso);
		evInfo.subLeadingEle_ptTracksIso.push_back(subLeadingElePtTracksIso);
		evInfo.subLeadingEle_innerLayerLostHits.push_back(subLeadingEleMissingHits);
		evInfo.subLeadingEle_dxy.push_back(subLeadingEleDxy);
		evInfo.subLeadingEle_eOverP.push_back(subLeadingEleEoverP);
		evInfo.subLeadingEle_id.push_back(subLeadingEleId);
		evInfo.subLeadingEle_passCutBasedEleId.push_back(subLeadingElePassCutBasedId);

		evInfo.leadingMuon_isHighPt.push_back(leadingMuonIsHighPt);
		evInfo.subLeadingMuon_isHighPt.push_back(subLeadingMuonIsHighPt);


		// // -- jets
		// int ijetVector = diLeptonDiJet->vertexIndex();
		// // cout << "ijetVector = " << ijetVector << endl; //sempre = 0 -> se all'interno del loop sui DLDJ prende piu` volte lo stesso vettore di jet

		// Ptr<vector<flashgg::Jet> > jetVector = jets->ptrAt( ijetVector );
		// // if (jetVecotr->size() > 0) cout << "jetVector size " << jetVector->size() << endl;

		// for (UInt_t ijet = 0 ; ijet < jetVector->size(); ijet++){
		// 	flashgg::Jet jet = jetVector->at( ijet );

		// int isMatchedToGen =  -1;
		// if( ! iEvent.isRealData() ) isMatchedToGen = jetMatchingToGen(jet, genJets); 

		// 	evInfo.jet_e.push_back(jet.energy());
		// 	evInfo.jet_pt.push_back(jet.pt());
		// 	evInfo.jet_eta.push_back(jet.eta());
		// 	evInfo.jet_phi.push_back(jet.phi());
		// 	evInfo.jet_bdiscriminant.push_back(jet.bDiscriminator( bTag_ ));
		// 	evInfo.jet_hadronFlavour.push_back(jet.hadronFlavour());
		// 	evInfo.jet_partonFlavour.push_back(jet.partonFlavour());
		// 	evInfo.jet_isMatchedToGen.push_back(isMatchedToGen);
		// }

	}
	// cout << "exit dldj loop" << endl;

	// --- fill the tree  
	// eventTree->Fill();
	if ( passTrigger || (!iEvent.isRealData()) ) eventTree->Fill();  //per i data fillo il tree solo per gli eventi che passano i triggers che mi interessano 

	// cout << "fillo tree" << endl;

}
// ******************************************************************************************



// ******************************************************************************************
void miniTreeMaker::endJob() {
	cout << "Total number of generated dldj before preselection = "<< ngen << endl;
	cout << "Number of generated dldj after preselection        = "<< ngenPre << endl;
	cout << "Total number of dldj before preselection           = "<< ndldj << endl;
	cout << "Number of dldj after preselection                  = "<< npre << endl;
	cout << "Total number of electrons                          = "<< nEle << endl;
	cout << "Number of electrons after eta and conversions cuts = "<< nEleGood << endl;
	cout << "Number of electrons passing HEEP id                = "<< nElePassingHEEPid << endl;
	cout << "Total number of muons                              = "<< nmuons << endl;
	cout << "Number of muons after eta cuts                     = "<< nmuonsGood << endl;
 
} // end of endJob
// ******************************************************************************************



// ******************************************************************************************
void miniTreeMaker::initEventStructure() {
	// per-event tree:
	evInfo.run = -999;
	evInfo.event = -999.;
	evInfo.lumi = -999.;
	evInfo.rho = -999.;
	evInfo.weight = -999.;
	evInfo.puweight = -999.;
	evInfo.nvtx = -999;
	evInfo.npu = -999;
	evInfo.passEEJJhlt = -1;
	evInfo.passMMJJhlt = -1;
	evInfo.passEMJJhlt = -1;
	evInfo.passTandPEEhlt = -1;
	evInfo.passTandPMMhlt = -1;

	evInfo.vtx_x .clear();
	evInfo.vtx_y .clear();
	evInfo.vtx_z .clear();

	evInfo.ele_e .clear();
	evInfo.ele_pt .clear();
	evInfo.ele_eta .clear();
	evInfo.ele_phi .clear();
	evInfo.ele_idmva .clear();
	evInfo.ele_id . clear();
	evInfo.ele_iso .clear();
	evInfo.ele_dz .clear();
	evInfo.ele_d0 .clear();
	evInfo.ele_passHEEPId .clear();
	evInfo.ele_passCutBasedEleId .clear();
	evInfo.ele_isMatchedToGen .clear();
	evInfo.ele_charge .clear();
	evInfo.ele_etaSC .clear();
	evInfo.ele_isEcalDriven .clear();
	evInfo.ele_dEtaIn .clear(); 
	evInfo.ele_dPhiIn .clear();
	evInfo.ele_hOverE .clear();
	evInfo.ele_full5x5_r9 .clear();  
	evInfo.ele_full5x5_sigmaIetaIeta .clear();
	evInfo.ele_full5x5_E5x5 .clear();
	evInfo.ele_full5x5_E1x5 .clear();
	evInfo.ele_full5x5_E2x5 .clear();
	evInfo.ele_full5x5_E2x5_Over_E5x5 .clear();
	evInfo.ele_full5x5_E1x5_Over_E5x5 .clear();
	evInfo.ele_EmHadDepth1Iso .clear();
	evInfo.ele_ptTracksIso .clear();
	evInfo.ele_innerLayerLostHits .clear();
	evInfo.ele_dxy .clear();
	evInfo.ele_eOverP .clear();

	evInfo.mu_e .clear();
	evInfo.mu_pt .clear();
	evInfo.mu_eta .clear();
	evInfo.mu_phi .clear();
	evInfo.mu_iso .clear();
	evInfo.mu_isTight .clear();
	evInfo.mu_isMedium .clear();
	evInfo.mu_isLoose .clear();
	evInfo.mu_isHighPt .clear();
	evInfo.mu_isMatchedToGen .clear();
	evInfo.mu_charge .clear();
	evInfo.mu_dz .clear();
	evInfo.mu_dxy .clear();

	evInfo.jet_e .clear();
	evInfo.jet_pt .clear();
	evInfo.jet_eta .clear();
	evInfo.jet_phi .clear();
	evInfo.jet_bdiscriminant .clear();
	evInfo.jet_partonFlavour .clear();
	evInfo.jet_hadronFlavour .clear();
	evInfo.jet_isMatchedToGen .clear();

	
	evInfo.isEEJJ .clear();
	evInfo.isEETT .clear();
	evInfo.isMMJJ .clear();
	evInfo.isMMTT .clear();
	evInfo.isEMJJ .clear();

	evInfo.isSignalRegion .clear();
	evInfo.isLowMllCR .clear();
	evInfo.isLowMlljjCR .clear();

	evInfo.isBB .clear();
	evInfo.isEE .clear();
	evInfo.isEB .clear();

	evInfo.passPreselections .clear();

	evInfo.leadingLepton_e .clear(); 
	evInfo.leadingLepton_pt .clear();
	evInfo.leadingLepton_eta .clear(); 
	evInfo.leadingLepton_phi .clear();  
	evInfo.leadingLepton_charge .clear();

	evInfo.subLeadingLepton_e .clear();
	evInfo.subLeadingLepton_pt .clear();
	evInfo.subLeadingLepton_eta .clear(); 
	evInfo.subLeadingLepton_phi .clear();  
	evInfo.subLeadingLepton_charge .clear();

	evInfo.leadingJet_e .clear(); 
	evInfo.leadingJet_pt .clear();
	evInfo.leadingJet_eta .clear(); 
	evInfo.leadingJet_phi .clear();

	evInfo.subLeadingJet_e .clear();
	evInfo.subLeadingJet_pt .clear();
	evInfo.subLeadingJet_eta .clear();
	evInfo.subLeadingJet_phi .clear();  

	evInfo.dRLeadLeptonLeadJet .clear();
	evInfo.dRLeadLeptonSubLeadJet .clear();
	evInfo.dRSubLeadLeptonLeadJet .clear();
	evInfo.dRSubLeadLeptonSubLeadJet .clear();

	evInfo.diLeptonDiJet_vtxIndex .clear();
	evInfo.diLeptonDiJet_sumPt .clear();
	evInfo.diLeptonDiJet_invMass .clear();
	evInfo.diLepton_invMass .clear();
	evInfo.diLepton_pt .clear();
	evInfo.diJet_invMass .clear();
	evInfo.diJetLeadingLepton_invMass .clear();
	evInfo.diJetSubLeadingLepton_invMass .clear();

	evInfo.leadingEle_passHEEPId .clear();
	evInfo.leadingEle_etaSC .clear();
	evInfo.leadingEle_isEcalDriven .clear();
	evInfo.leadingEle_dEtaIn .clear();  
	evInfo.leadingEle_dPhiIn .clear();
	evInfo.leadingEle_hOverE .clear();
	evInfo.leadingEle_full5x5_r9 .clear();  
	evInfo.leadingEle_full5x5_sigmaIetaIeta .clear();
	evInfo.leadingEle_full5x5_E5x5 .clear();
	evInfo.leadingEle_full5x5_E1x5 .clear();
	evInfo.leadingEle_full5x5_E2x5 .clear();
	evInfo.leadingEle_full5x5_E2x5_Over_E5x5 .clear();
	evInfo.leadingEle_full5x5_E1x5_Over_E5x5 .clear();
	evInfo.leadingEle_EmHadDepth1Iso .clear();
	evInfo.leadingEle_ptTracksIso .clear();
	evInfo.leadingEle_innerLayerLostHits .clear();
	evInfo.leadingEle_dxy .clear();
	evInfo.leadingEle_eOverP .clear();
	evInfo.leadingEle_id .clear();
	evInfo.leadingEle_passCutBasedEleId .clear();

	evInfo.subLeadingEle_passHEEPId .clear();
	evInfo.subLeadingEle_etaSC .clear();
	evInfo.subLeadingEle_isEcalDriven .clear();
	evInfo.subLeadingEle_dEtaIn .clear();  
	evInfo.subLeadingEle_dPhiIn .clear();
	evInfo.subLeadingEle_hOverE .clear();
	evInfo.subLeadingEle_full5x5_r9 .clear();  
	evInfo.subLeadingEle_full5x5_sigmaIetaIeta .clear();
	evInfo.subLeadingEle_full5x5_E5x5 .clear();
	evInfo.subLeadingEle_full5x5_E1x5 .clear();
	evInfo.subLeadingEle_full5x5_E2x5 .clear();
	evInfo.subLeadingEle_full5x5_E2x5_Over_E5x5 .clear();
	evInfo.subLeadingEle_full5x5_E1x5_Over_E5x5 .clear();
	evInfo.subLeadingEle_EmHadDepth1Iso .clear();
	evInfo.subLeadingEle_ptTracksIso .clear();
	evInfo.subLeadingEle_innerLayerLostHits .clear();
	evInfo.subLeadingEle_dxy .clear();
	evInfo.subLeadingEle_eOverP .clear();
	evInfo.subLeadingEle_id .clear();
	evInfo.subLeadingEle_passCutBasedEleId .clear();

	evInfo.leadingMuon_isHighPt .clear();
	evInfo.subLeadingMuon_isHighPt .clear();

}
// ******************************************************************************************