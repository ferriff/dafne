#include "dafne/Validation/interface/miniTreeMaker.h"


miniTreeMaker::miniTreeMaker( const edm::ParameterSet &iConfig, TFileDirectory& fs, edm::ConsumesCollector && cc ):
	edm::BasicAnalyzer::BasicAnalyzer(iConfig, fs),
	genParticleToken_( cc.consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "genParticleTag" ) ) ),
	genInfoToken_(cc.consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "generatorInfo" ) ) ),
	PileUpToken_(cc.consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ),
	vertexToken_( cc.consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
	DiLeptonDiJetToken_( cc.consumes<View<flashgg::DiLeptonDiJetCandidate> >( iConfig.getParameter<InputTag> ( "DiLeptonDiJetTag" ) ) ),
	jetsToken_( cc.consumes<View<vector<flashgg::Jet> > >( iConfig.getParameter<InputTag> ( "JetsTag" ) ) ),
	genJetToken_( cc.consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
	electronToken_( cc.consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
	muonToken_( cc.consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
	triggerBitsToken_( cc.consumes<edm::TriggerResults>( iConfig.getParameter<InputTag>( "triggerBits" ) ) ),
	rhoToken_(cc.consumes<double>(iConfig.getParameter <edm::InputTag>("rhoFixedGridCollection" ) ) )
{
	bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
	lumiWeight_ = iConfig.getUntrackedParameter<double>( "lumiWeight", 1000. ); //pb                                                                                                                              
  
	globalVarsDumper_ = new GlobalVariablesDumper( iConfig.getParameter<edm::ParameterSet>( "globalVariables" ), std::forward<edm::ConsumesCollector>(cc) );
  
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

	eventTree->Branch( "weight", &evInfo.weight, "weight/F" );
	eventTree->Branch( "puweight", &evInfo.puweight,"puweight/F");

	eventTree->Branch( "nvtx", &evInfo.nvtx, "nvtx/I" );
	eventTree->Branch( "npu", &evInfo.npu, "npu/I" );

	eventTree->Branch( "passEEJJhlt", &evInfo.passEEJJhlt, "passEEJJhlt/I" );
	eventTree->Branch( "passMMJJhlt", &evInfo.passMMJJhlt, "passMMJJhlt/I" );
	eventTree->Branch( "passEMJJhlt", &evInfo.passEMJJhlt, "passEMJJhlt/I" );
	eventTree->Branch( "passTandPEEhlt", &evInfo.passTandPEEhlt, "passTandPEEhlt/I" );
	eventTree->Branch( "passTandPMMhlt", &evInfo.passTandPMMhlt, "passTandPMMhlt/I" );

	eventTree->Branch( "ele_e", &evInfo.ele_e);
	eventTree->Branch( "ele_pt", &evInfo.ele_pt);
	eventTree->Branch( "ele_eta", &evInfo.ele_eta);
	eventTree->Branch( "ele_phi", &evInfo.ele_phi);
	eventTree->Branch( "ele_idmva", &evInfo.ele_idmva);
	eventTree->Branch( "ele_iso", &evInfo.ele_iso);
	eventTree->Branch( "ele_dz", &evInfo.ele_dz);
	eventTree->Branch( "ele_d0", &evInfo.ele_d0);
	eventTree->Branch( "ele_passHEEPId", &evInfo.ele_passHEEPId);
	eventTree->Branch( "ele_isMatchedToGen", &evInfo.ele_isMatchedToGen);
	eventTree->Branch( "ele_charge", &evInfo.ele_charge);
	eventTree->Branch( "ele_etaSC", &evInfo.ele_etaSC);
	eventTree->Branch( "ele_isEcalDriven", &evInfo.ele_isEcalDriven);
	eventTree->Branch( "ele_dEtaIn", &evInfo.ele_dEtaIn);
	eventTree->Branch( "ele_dPhiIn", &evInfo.ele_dPhiIn);
	eventTree->Branch( "ele_hOverE", &evInfo.ele_hOverE);
	eventTree->Branch( "ele_full5x5_sigmaIetaIeta", &evInfo.ele_full5x5_sigmaIetaIeta);
	eventTree->Branch( "ele_full5x5_E5x5", &evInfo.ele_full5x5_E5x5);
	eventTree->Branch( "ele_full5x5_E1x5", &evInfo.ele_full5x5_E1x5);
	eventTree->Branch( "ele_full5x5_E2x5", &evInfo.ele_full5x5_E2x5);
	eventTree->Branch( "ele_full5x5_r9", &evInfo.ele_full5x5_r9);
	eventTree->Branch( "ele_full5x5_E2x5_Over_E5x5", &evInfo.ele_full5x5_E2x5_Over_E5x5);
	eventTree->Branch( "ele_full5x5_E1x5_Over_E5x5", &evInfo.ele_full5x5_E1x5_Over_E5x5);
	eventTree->Branch( "ele_innerLayerLostHits", &evInfo.ele_innerLayerLostHits);

	eventTree->Branch( "mu_e", &evInfo.mu_e);
	eventTree->Branch( "mu_pt", &evInfo.mu_pt);
	eventTree->Branch( "mu_eta", &evInfo.mu_eta);
	eventTree->Branch( "mu_phi", &evInfo.mu_phi);
	eventTree->Branch( "mu_iso", &evInfo.mu_iso);
	eventTree->Branch( "mu_isTight", &evInfo.mu_isTight);
	eventTree->Branch( "mu_isMedium", &evInfo.mu_isMedium);
	eventTree->Branch( "mu_isLoose", &evInfo.mu_isLoose);
	eventTree->Branch( "mu_isHighPt", &evInfo.mu_isHighPt);
	eventTree->Branch( "mu_isMatchedToGen", &evInfo.mu_isMatchedToGen);
	eventTree->Branch( "mu_charge", &evInfo.mu_charge);

	eventTree->Branch( "jet_e", &evInfo.jet_e);
	eventTree->Branch( "jet_pt", &evInfo.jet_pt);
	eventTree->Branch( "jet_eta", &evInfo.jet_eta);
	eventTree->Branch( "jet_phi", &evInfo.jet_phi);
	eventTree->Branch( "jet_bdiscriminant", &evInfo.jet_bdiscriminant);
	eventTree->Branch( "jet_partonFlavour", &evInfo.jet_partonFlavour);
	eventTree->Branch( "jet_hadronFlavour", &evInfo.jet_hadronFlavour);
	eventTree->Branch( "jet_isMatchedToGen", &evInfo.jet_isMatchedToGen);

	eventTree->Branch( "isEEJJ", &evInfo.isEEJJ); 
	eventTree->Branch( "isEETT", &evInfo.isEETT); 
	eventTree->Branch( "isMMJJ", &evInfo.isMMJJ); 
	eventTree->Branch( "isMMTT", &evInfo.isMMTT); 
	eventTree->Branch( "isEMJJ", &evInfo.isEMJJ); 

	eventTree->Branch( "isSignalRegion", &evInfo.isSignalRegion);
	eventTree->Branch( "isLowMllCR", &evInfo.isLowMllCR);
	eventTree->Branch( "isLowMlljjCR", &evInfo.isLowMlljjCR);

	eventTree->Branch( "isBB", &evInfo.isBB);
	eventTree->Branch( "isEE", &evInfo.isEE);
	eventTree->Branch( "isEB", &evInfo.isEB);

	eventTree->Branch( "passPreselections", &evInfo.passPreselections);

	eventTree->Branch( "leadingLepton_pt", &evInfo.leadingLepton_pt); 
	eventTree->Branch( "leadingLepton_eta", &evInfo.leadingLepton_eta); 
	eventTree->Branch( "leadingLepton_phi", &evInfo.leadingLepton_phi); 

	eventTree->Branch( "subLeadingLepton_pt", &evInfo.subLeadingLepton_pt); 
	eventTree->Branch( "subLeadingLepton_eta", &evInfo.subLeadingLepton_eta); 
	eventTree->Branch( "subLeadingLepton_phi", &evInfo.subLeadingLepton_phi); 

	eventTree->Branch( "leadingJet_pt", &evInfo.leadingJet_pt); 
	eventTree->Branch( "leadingJet_eta", &evInfo.leadingJet_eta);
	eventTree->Branch( "leadingJet_phi", &evInfo.leadingJet_phi); 

	eventTree->Branch( "subLeadingJet_pt", &evInfo.subLeadingJet_pt); 
	eventTree->Branch( "subLeadingJet_eta", &evInfo.subLeadingJet_eta); 
	eventTree->Branch( "subLeadingJet_phi", &evInfo.subLeadingJet_phi); 

	eventTree->Branch( "diLeptonDiJet_vtxIndex", &evInfo.diLeptonDiJet_vtxIndex);
	eventTree->Branch( "diLeptonDiJet_sumPt", &evInfo.diLeptonDiJet_sumPt);
	eventTree->Branch( "diLeptonDiJet_invMass", &evInfo.diLeptonDiJet_invMass); 
	eventTree->Branch( "diLepton_invMass", &evInfo.diLepton_invMass); 
	eventTree->Branch( "diJet_invMass", &evInfo.diJet_invMass); 

	eventTree->Branch( "leadingEle_passHEEPId", &evInfo.leadingEle_passHEEPId);
	eventTree->Branch( "leadingEle_etaSC", &evInfo.leadingEle_etaSC);
	eventTree->Branch( "leadingEle_isEcalDriven", &evInfo.leadingEle_isEcalDriven);
	eventTree->Branch( "leadingEle_dEtaIn", &evInfo.leadingEle_dEtaIn);
	eventTree->Branch( "leadingEle_dPhiIn", &evInfo.leadingEle_dPhiIn);
	eventTree->Branch( "leadingEle_hOverE", &evInfo.leadingEle_hOverE);
	eventTree->Branch( "leadingEle_full5x5_sigmaIetaIeta", &evInfo.leadingEle_full5x5_sigmaIetaIeta);
	eventTree->Branch( "leadingEle_full5x5_E5x5", &evInfo.leadingEle_full5x5_E5x5);
	eventTree->Branch( "leadingEle_full5x5_E1x5", &evInfo.leadingEle_full5x5_E1x5);
	eventTree->Branch( "leadingEle_full5x5_E2x5", &evInfo.leadingEle_full5x5_E2x5);
	eventTree->Branch( "leadingEle_full5x5_r9", &evInfo.leadingEle_full5x5_r9);
	eventTree->Branch( "leadingEle_full5x5_E2x5_Over_E5x5", &evInfo.leadingEle_full5x5_E2x5_Over_E5x5);
 	eventTree->Branch( "leadingEle_full5x5_E1x5_Over_E5x5", &evInfo.leadingEle_full5x5_E1x5_Over_E5x5);
	eventTree->Branch( "leadingEle_innerLayerLostHits", &evInfo.leadingEle_innerLayerLostHits);

	eventTree->Branch( "subLeadingEle_passHEEPId", &evInfo.subLeadingEle_passHEEPId);
	eventTree->Branch( "subLeadingEle_etaSC", &evInfo.subLeadingEle_etaSC);
	eventTree->Branch( "subLeadingEle_isEcalDriven", &evInfo.subLeadingEle_isEcalDriven);
	eventTree->Branch( "subLeadingEle_dEtaIn", &evInfo.subLeadingEle_dEtaIn);
	eventTree->Branch( "subLeadingEle_dPhiIn", &evInfo.subLeadingEle_dPhiIn);
	eventTree->Branch( "subLeadingEle_hOverE", &evInfo.subLeadingEle_hOverE);
	eventTree->Branch( "subLeadingEle_full5x5_sigmaIetaIeta", &evInfo.subLeadingEle_full5x5_sigmaIetaIeta);
	eventTree->Branch( "subLeadingEle_full5x5_E5x5", &evInfo.subLeadingEle_full5x5_E5x5);
	eventTree->Branch( "subLeadingEle_full5x5_E1x5", &evInfo.subLeadingEle_full5x5_E1x5);
	eventTree->Branch( "subLeadingEle_full5x5_E2x5", &evInfo.subLeadingEle_full5x5_E2x5);
	eventTree->Branch( "subLeadingEle_full5x5_r9", &evInfo.subLeadingEle_full5x5_r9);
	eventTree->Branch( "subLeadingEle_full5x5_E2x5_Over_E5x5", &evInfo.subLeadingEle_full5x5_E2x5_Over_E5x5);
 	eventTree->Branch( "subLeadingEle_full5x5_E1x5_Over_E5x5", &evInfo.subLeadingEle_full5x5_E1x5_Over_E5x5);
	eventTree->Branch( "subLeadingEle_innerLayerLostHits", &evInfo.subLeadingEle_innerLayerLostHits);

}
// ******************************************************************************************




// ******************************************************************************************
// analyzer
//
void miniTreeMaker::analyze(const edm::EventBase& evt)
{
	const edm::Event *fullEvent = dynamic_cast<const edm::Event *>(&evt);
	const edm::Event &iEvent = (*fullEvent);  
  

  	//-------------- access edm objects

	//--- only if MC
	Handle<View<reco::GenParticle> > genParticles;
	Handle<GenEventInfoProduct> genInfo;
	Handle<View< PileupSummaryInfo> > PileupInfos;
	Handle<View<reco::GenJet> > genJets;

	if ( !iEvent.isRealData() ) {
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
	
	Handle<edm::TriggerResults> triggerBits;
	iEvent.getByToken( triggerBitsToken_, triggerBits );
  
	Handle<double> rhoHandle;
	iEvent.getByToken( rhoToken_, rhoHandle );
	double rho = *( rhoHandle.product() );




	//-------------- initialize tree
	initEventStructure();


	
	//-------------- check if event passes HLT: "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*"  
	const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerBits );

	for( unsigned index = 0; index < triggerNames.size(); ++index ) {

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_DoubleEle33_CaloIdL_MW") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passEEJJhlt =  triggerBits->accept( index );
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Mu50") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passMMJJhlt =  triggerBits->accept( index );
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passEMJJhlt =  triggerBits->accept( index );
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Ele27_WPTight_Gsf") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passTandPEEhlt =  triggerBits->accept( index );
		}

		if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_IsoMu24") || (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_IsoMu27") ) {
			// cout << (triggerNames.triggerName( index )).c_str() <<endl;
			evInfo.passTandPMMhlt =  triggerBits->accept( index );
		}

	}


	globalVarsDumper_->fill( iEvent );

	evInfo.run = globalVarsDumper_->cache().run;
	evInfo.event = globalVarsDumper_->cache().event;		
	evInfo.lumi = globalVarsDumper_->cache().lumi;

	// -- event weight (Lumi x cross section x gen weight)
	float w = 1.;
	if( ! iEvent.isRealData() ) {
		w = lumiWeight_;
		if( genInfo.isValid() ) {
			const auto &weights = genInfo->weights();
			if( ! weights.empty() ) {
				w *= weights[0];
			}
		}
	}
	evInfo.weight = w;

	// -- pileup weight
	if( globalVarsDumper_->puReWeight() ) {
		evInfo.puweight = globalVarsDumper_->cache().puweight;
	}

	// -- number of reco vertices
	evInfo.nvtx = vertices->size() ;
	// cout << vertices->size() << "nvertices " << endl;

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
	// cout << "pu" << pu << endl;


	// -- electrons
	for (UInt_t iele = 0 ; iele < electrons->size(); iele++){
		nEle++;

		edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( iele );
		if (fabs(electron->eta()) > 2.4) { continue; }
		// if (electron->pt() < electronPtThreshold_) { continue; }  
		if( electron->hasMatchedConversion() ) { continue; } // remove conversions
		nEleGood++;

		float isol = electronIsolation(electron, rho); 

		Ptr<reco::Vertex> ele_vtx = chooseElectronVertex( electron,  vertices->ptrs() );
		float dz = electron->gsfTrack()->dz( ele_vtx->position() );
		float d0 = electron->gsfTrack()->dxy( ele_vtx->position() ); 
	
		int passHEEPId = passHEEPIdCuts( electron, vertices->ptrs() );
		if (passHEEPId) nElePassingHEEPid++;

		int mcMatch = -1;
		if( ! iEvent.isRealData() ) mcMatch = electronMatchingToGen(electron, genParticles); 

		evInfo.ele_e.push_back(electron->energy());
		evInfo.ele_pt.push_back(electron->pt());
		evInfo.ele_eta.push_back(electron->eta());
		evInfo.ele_phi.push_back(electron->phi());
		evInfo.ele_idmva.push_back(electron->nonTrigMVA());
		evInfo.ele_iso.push_back(isol);
		evInfo.ele_dz.push_back(dz);
		evInfo.ele_d0.push_back(d0);
		evInfo.ele_passHEEPId.push_back(passHEEPId);
		evInfo.ele_isMatchedToGen.push_back(mcMatch);
		evInfo.ele_charge.push_back(electron->charge());
		evInfo.ele_etaSC.push_back(electron->superCluster()->eta());
		// evInfo.ele_isEcalDriven.push_back(electron->gsfTrack()->ecalDriven());
		evInfo.ele_dEtaIn.push_back(electron->deltaEtaSuperClusterTrackAtVtx());  //ok ?
		evInfo.ele_dPhiIn.push_back(electron->deltaPhiSuperClusterTrackAtVtx());
		evInfo.ele_hOverE.push_back(electron->hcalOverEcal());
		evInfo.ele_full5x5_sigmaIetaIeta.push_back(electron->full5x5_sigmaIetaIeta());
		evInfo.ele_full5x5_E5x5.push_back(electron->full5x5_e5x5());
		evInfo.ele_full5x5_E1x5.push_back(electron->full5x5_e1x5());
		evInfo.ele_full5x5_E2x5.push_back(electron->full5x5_e2x5Max());
		evInfo.ele_full5x5_r9.push_back(electron->full5x5_r9());  //o r9() ?
		evInfo.ele_full5x5_E2x5_Over_E5x5.push_back((electron->full5x5_e2x5Max()) / (electron->full5x5_e5x5()));
		evInfo.ele_full5x5_E1x5_Over_E5x5.push_back((electron->full5x5_e1x5()) / (electron->full5x5_e5x5()));
		evInfo.ele_innerLayerLostHits.push_back(electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS));

	// double vtx_dz = 1000000.;
	// unsigned int min_dz_vtx = -1;
				
	// for( unsigned int vtxi = 0; vtxi < pvPointers.size(); vtxi++ ) {            
	// 	Ptr<reco::Vertex> vtx = pvPointers[vtxi];            
	// 	if( vtx_dz > fabs(electron->gsfTrack()->dz( vtx->position() )) ) {                
	// 		vtx_dz = fabs( electron->gsfTrack()->dz( vtx->position() ) );
	// 		min_dz_vtx = vtxi;
	// 	}
	// }
				
	// Ptr<reco::Vertex> best_vtx = pvPointers[min_dz_vtx];
	// float dXY = fabs( electron->gsfTrack()->dxy( best_vtx->position()) ) ;  

	}       


	// -- muons
	for (UInt_t imu = 0 ; imu < muons->size(); imu++){
		nmuons++;

		edm::Ptr<flashgg::Muon> muon = muons->ptrAt( imu );
		if (fabs(muon->eta()) > 2.4) { continue; }
		// if (muon->pt()  < muonPtThreshold_) { continue; }  
		nmuonsGood++;

		// muon ID and isolation: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
		float muPFCombRelIso = ( muon->pfIsolationR04().sumChargedHadronPt + max( 0.,muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - 0.5 * muon->pfIsolationR04().sumPUPt ) ) / ( muon->pt() );

		int vtxInd = 0;
		double dzmin = 9999;
		for( size_t ivtx = 0 ; ivtx < vertices->size(); ivtx++ ) {
			Ptr<reco::Vertex> vtx = vertices->ptrAt(ivtx);
			if( !muon->innerTrack() ) { continue; }
			if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {
				dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
				vtxInd = ivtx;
			}
		}
		Ptr<reco::Vertex> muonVtx = vertices->ptrAt(vtxInd);

		int mcMatch =  -1;
		if( ! iEvent.isRealData() ) mcMatch = muonMatchingToGen(muon, genParticles); 

		evInfo.mu_e.push_back(muon->energy());
		evInfo.mu_pt.push_back(muon->pt());
		evInfo.mu_eta.push_back(muon->eta());
		evInfo.mu_phi.push_back(muon->phi());
		evInfo.mu_iso.push_back(muPFCombRelIso);
		evInfo.mu_isTight.push_back(muon::isTightMuon( *muon, *muonVtx ));
		evInfo.mu_isMedium.push_back(muon::isMediumMuon( *muon ));
		evInfo.mu_isLoose.push_back(muon::isLooseMuon( *muon ));
		evInfo.mu_isHighPt.push_back(muon::isHighPtMuon( *muon, *muonVtx ));
		evInfo.mu_isMatchedToGen.push_back(mcMatch); 
		evInfo.mu_charge.push_back(muon->charge());

	}


	// -- jets
	// for (UInt_t ijetVector = 0 ; ijetVector < jets->size(); ijetVector++){
		// edm::Ptr<vector<flashgg::Jet> > jetVector = jets->ptrAt( ijetVector );
	edm::Ptr<vector<flashgg::Jet> > jetVector = jets->ptrAt( 0 );
	// if (jetVecotr->size() > 0) cout << "jetVector size " << jetVector->size() << endl;

	for (UInt_t ijet = 0 ; ijet < jetVector->size(); ijet++){
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



	// -- diLeptonDiJets
	// cout << "dldj size " << diLeptonDiJets->size() << endl;

	for ( unsigned int idldj = 0; idldj < diLeptonDiJets->size(); idldj++){
		edm::Ptr<flashgg::DiLeptonDiJetCandidate> diLeptonDiJet = diLeptonDiJets->ptrAt( idldj );        
		
		if (! iEvent.isRealData() ) { 
			ngen++;
			if ( passDiLeptonDiJetPreselection(diLeptonDiJet) ) ngenPre++; 
		} else {
			ndldj++;
			if ( passDiLeptonDiJetPreselection(diLeptonDiJet) ) npre++; 
		}
					
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

		evInfo.leadingLepton_pt.push_back(diLeptonDiJet->leadingLeptonPt());
		evInfo.leadingLepton_eta.push_back(diLeptonDiJet->leadingLeptonEta());
		evInfo.leadingLepton_phi.push_back(diLeptonDiJet->leadingLeptonPhi());

		evInfo.subLeadingLepton_pt.push_back(diLeptonDiJet->subLeadingLeptonPt());
		evInfo.subLeadingLepton_eta.push_back(diLeptonDiJet->subLeadingLeptonEta());
		evInfo.subLeadingLepton_phi.push_back(diLeptonDiJet->subLeadingLeptonPhi());

		evInfo.leadingJet_pt.push_back(diLeptonDiJet->leadingJet()->pt());
		evInfo.leadingJet_eta.push_back(diLeptonDiJet->leadingJet()->eta());
		evInfo.leadingJet_phi.push_back(diLeptonDiJet->leadingJet()->phi());

		evInfo.subLeadingJet_pt.push_back(diLeptonDiJet->subLeadingJet()->pt());
		evInfo.subLeadingJet_eta.push_back(diLeptonDiJet->subLeadingJet()->eta());
		evInfo.subLeadingJet_phi.push_back(diLeptonDiJet->subLeadingJet()->phi());

		evInfo.diLeptonDiJet_vtxIndex.push_back(diLeptonDiJet->vertexIndex()); 
		evInfo.diLeptonDiJet_sumPt.push_back(diLeptonDiJet->sumPt());
		evInfo.diLeptonDiJet_invMass.push_back(diLeptonDiJet->invMass());
		evInfo.diLepton_invMass.push_back(diLeptonDiJet->diLeptonInvMass());
		evInfo.diJet_invMass.push_back( diJetInvMass(diLeptonDiJet) );


		if (diLeptonDiJet->isEEJJ() || diLeptonDiJet->isEETT()) {
			int leadElePassHEEPId = passHEEPIdCuts( diLeptonDiJet->leadingEle(), vertices->ptrs() );
			int subLeadElePassHEEPId = passHEEPIdCuts( diLeptonDiJet->subLeadingEle(), vertices->ptrs() );

			evInfo.leadingEle_passHEEPId.push_back(leadElePassHEEPId);
			evInfo.leadingEle_etaSC.push_back(diLeptonDiJet->leadingEle()->superCluster()->eta());
			// evInfo.leadingEle_isEcalDriven.push_back(diLeptonDiJet->leadingEle()->gsfTrack()->ecalDriven());
			evInfo.leadingEle_dEtaIn.push_back(diLeptonDiJet->leadingEle()->deltaEtaSuperClusterTrackAtVtx());  //ok ?
			evInfo.leadingEle_dPhiIn.push_back(diLeptonDiJet->leadingEle()->deltaPhiSuperClusterTrackAtVtx());
			evInfo.leadingEle_hOverE.push_back(diLeptonDiJet->leadingEle()->hcalOverEcal());
			evInfo.leadingEle_full5x5_sigmaIetaIeta.push_back(diLeptonDiJet->leadingEle()->full5x5_sigmaIetaIeta());
			evInfo.leadingEle_full5x5_E5x5.push_back(diLeptonDiJet->leadingEle()->full5x5_e5x5());
			evInfo.leadingEle_full5x5_E1x5.push_back(diLeptonDiJet->leadingEle()->full5x5_e1x5());
			evInfo.leadingEle_full5x5_E2x5.push_back(diLeptonDiJet->leadingEle()->full5x5_e2x5Max());
			evInfo.leadingEle_full5x5_r9.push_back(diLeptonDiJet->leadingEle()->full5x5_r9());  //r9() 
			evInfo.leadingEle_full5x5_E2x5_Over_E5x5.push_back((diLeptonDiJet->leadingEle()->full5x5_e2x5Max()) / (diLeptonDiJet->leadingEle()->full5x5_e5x5()));
			evInfo.leadingEle_full5x5_E1x5_Over_E5x5.push_back((diLeptonDiJet->leadingEle()->full5x5_e1x5()) / (diLeptonDiJet->leadingEle()->full5x5_e5x5()));
			evInfo.leadingEle_innerLayerLostHits.push_back(diLeptonDiJet->leadingEle()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS));

			evInfo.subLeadingEle_passHEEPId.push_back(subLeadElePassHEEPId);		
			evInfo.subLeadingEle_etaSC.push_back(diLeptonDiJet->subLeadingEle()->superCluster()->eta());			
			// evInfo.subLeadingEle_isEcalDriven.push_back(diLeptonDiJet->subLeadingEle()->gsfTrack()->ecalDriven());
			evInfo.subLeadingEle_dEtaIn.push_back(diLeptonDiJet->subLeadingEle()->deltaEtaSuperClusterTrackAtVtx());  //ok ?
			evInfo.subLeadingEle_dPhiIn.push_back(diLeptonDiJet->subLeadingEle()->deltaPhiSuperClusterTrackAtVtx());
			evInfo.subLeadingEle_hOverE.push_back(diLeptonDiJet->subLeadingEle()->hcalOverEcal());
			evInfo.subLeadingEle_full5x5_sigmaIetaIeta.push_back(diLeptonDiJet->subLeadingEle()->full5x5_sigmaIetaIeta());
			evInfo.subLeadingEle_full5x5_E5x5.push_back(diLeptonDiJet->subLeadingEle()->full5x5_e5x5());
			evInfo.subLeadingEle_full5x5_E1x5.push_back(diLeptonDiJet->subLeadingEle()->full5x5_e1x5());
			evInfo.subLeadingEle_full5x5_E2x5.push_back(diLeptonDiJet->subLeadingEle()->full5x5_e2x5Max());
			evInfo.subLeadingEle_full5x5_r9.push_back(diLeptonDiJet->subLeadingEle()->full5x5_r9());  //r9() 
			evInfo.subLeadingEle_full5x5_E2x5_Over_E5x5.push_back((diLeptonDiJet->subLeadingEle()->full5x5_e2x5Max()) / (diLeptonDiJet->subLeadingEle()->full5x5_e5x5()));
			evInfo.subLeadingEle_full5x5_E1x5_Over_E5x5.push_back((diLeptonDiJet->subLeadingEle()->full5x5_e1x5()) / (diLeptonDiJet->subLeadingEle()->full5x5_e5x5()));
			evInfo.subLeadingEle_innerLayerLostHits.push_back(diLeptonDiJet->subLeadingEle()->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS));
		}
	
		// // -- jets
		// int ijetVector = diLeptonDiJet->vertexIndex();
		// // cout << "ijetVector = " << ijetVector << endl; //sempre = 0 -> se all'interno del loop sui DLDJ prende piu` volte lo stesso vettore di jet

		// edm::Ptr<vector<flashgg::Jet> > jetVector = jets->ptrAt( ijetVector );
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


	// --- fill the tree
	eventTree->Fill();
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

	evInfo.weight = -999.;
	evInfo.puweight = -999.;

	evInfo.nvtx = -999;
	evInfo.npu = -999;

	evInfo.passEEJJhlt = -1;
	evInfo.passMMJJhlt = -1;
	evInfo.passEMJJhlt = -1;
	evInfo.passTandPEEhlt = -1;
	evInfo.passTandPMMhlt = -1;


	evInfo.ele_e .clear();
	evInfo.ele_pt .clear();
	evInfo.ele_eta .clear();
	evInfo.ele_phi .clear();
	evInfo.ele_idmva .clear();
	evInfo.ele_iso .clear();
	evInfo.ele_dz .clear();
	evInfo.ele_d0 .clear();
	evInfo.ele_passHEEPId .clear();
	evInfo.ele_isMatchedToGen .clear();
	evInfo.ele_charge .clear();
	evInfo.ele_etaSC .clear();
	evInfo.ele_isEcalDriven .clear();
	evInfo.ele_dEtaIn .clear(); 
	evInfo.ele_dPhiIn .clear();
	evInfo.ele_hOverE .clear();
	evInfo.ele_full5x5_sigmaIetaIeta .clear();
	evInfo.ele_full5x5_E5x5 .clear();
	evInfo.ele_full5x5_E1x5 .clear();
	evInfo.ele_full5x5_E2x5 .clear();
	evInfo.ele_full5x5_r9 .clear();  
	evInfo.ele_full5x5_E2x5_Over_E5x5 .clear();
	evInfo.ele_full5x5_E1x5_Over_E5x5 .clear();
	evInfo.ele_innerLayerLostHits .clear();

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

	evInfo.leadingLepton_pt .clear();
	evInfo.leadingLepton_eta .clear(); 
	evInfo.leadingLepton_phi .clear();  

	evInfo.subLeadingLepton_pt .clear();
	evInfo.subLeadingLepton_eta .clear(); 
	evInfo.subLeadingLepton_phi .clear();  

	evInfo.leadingJet_pt .clear();
	evInfo.leadingJet_eta .clear(); 
	evInfo.leadingJet_phi .clear();

	evInfo.subLeadingJet_pt .clear();
	evInfo.subLeadingJet_eta .clear();
	evInfo.subLeadingJet_phi .clear();  

	evInfo.diLeptonDiJet_vtxIndex .clear();
	evInfo.diLeptonDiJet_sumPt .clear();
	evInfo.diLeptonDiJet_invMass .clear();
	evInfo.diLepton_invMass .clear();
	evInfo.diJet_invMass .clear();

	evInfo.leadingEle_passHEEPId .clear();
	evInfo.leadingEle_etaSC .clear();
	evInfo.leadingEle_isEcalDriven .clear();
	evInfo.leadingEle_dEtaIn .clear();  
	evInfo.leadingEle_dPhiIn .clear();
	evInfo.leadingEle_hOverE .clear();
	evInfo.leadingEle_full5x5_sigmaIetaIeta .clear();
	evInfo.leadingEle_full5x5_E5x5 .clear();
	evInfo.leadingEle_full5x5_E1x5 .clear();
	evInfo.leadingEle_full5x5_E2x5 .clear();
	evInfo.leadingEle_full5x5_r9 .clear();  
	evInfo.leadingEle_full5x5_E2x5_Over_E5x5 .clear();
	evInfo.leadingEle_full5x5_E1x5_Over_E5x5 .clear();
	evInfo.leadingEle_innerLayerLostHits .clear();

	evInfo.subLeadingEle_passHEEPId .clear();
	evInfo.subLeadingEle_etaSC .clear();
	evInfo.subLeadingEle_isEcalDriven .clear();
	evInfo.subLeadingEle_dEtaIn .clear();  
	evInfo.subLeadingEle_dPhiIn .clear();
	evInfo.subLeadingEle_hOverE .clear();
	evInfo.subLeadingEle_full5x5_sigmaIetaIeta .clear();
	evInfo.subLeadingEle_full5x5_E5x5 .clear();
	evInfo.subLeadingEle_full5x5_E1x5 .clear();
	evInfo.subLeadingEle_full5x5_E2x5 .clear();
	evInfo.subLeadingEle_full5x5_r9 .clear();  
	evInfo.subLeadingEle_full5x5_E2x5_Over_E5x5 .clear();
	evInfo.subLeadingEle_full5x5_E1x5_Over_E5x5 .clear();
	evInfo.subLeadingEle_innerLayerLostHits .clear();

}
// ******************************************************************************************