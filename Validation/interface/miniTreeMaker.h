#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "flashgg/Taggers/interface/GlobalVariablesDumper.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"



using namespace std;
using namespace edm;
using namespace flashgg;


// define the structures used to create tree branches and fill the trees
struct eventInfo {

	int run;
	int event;
	int lumi;
	float rho;
	float weight;
	float puweight;
	int nvtx;
	int npu;
	int passEEJJhlt; 
	int passMMJJhlt; 
	int passEMJJhlt; 
	int passTandPEEhlt; 
	int passTandPMMhlt; 

	//variabili per oggetti 
	vector<float> vtx_x;
	vector<float> vtx_y;
	vector<float> vtx_z;

	vector<float> ele_e;
	vector<float> ele_pt;
	vector<float> ele_eta;
	vector<float> ele_phi;
	vector<float> ele_idmva;
	vector<unsigned> ele_id;
	vector<float> ele_iso;
	vector<float> ele_dz;
	vector<float> ele_d0;
	vector<bool> ele_passHEEPId; 
	vector<bool> ele_passCutBasedEleId;
	vector<int> ele_isMatchedToGen;
	vector<int> ele_charge;
	vector<float> ele_etaSC;
	vector<bool> ele_isEcalDriven;
	vector<float> ele_dEtaIn; 
	vector<float> ele_dPhiIn;
	vector<float> ele_hOverE;
	vector<float> ele_full5x5_r9;	
	vector<float> ele_full5x5_sigmaIetaIeta;
	vector<float> ele_full5x5_E5x5;
	vector<float> ele_full5x5_E1x5;
	vector<float> ele_full5x5_E2x5;  
	vector<float> ele_full5x5_E2x5_Over_E5x5;
	vector<float> ele_full5x5_E1x5_Over_E5x5;
	vector<float> ele_EmHadDepth1Iso;
	vector<float> ele_ptTracksIso;
	vector<int> ele_innerLayerLostHits;
	vector<float> ele_dxy;
	vector<float> ele_eOverP;

	vector<float> mu_e;
	vector<float> mu_pt;
	vector<float> mu_eta;
	vector<float> mu_phi;
	vector<float> mu_iso;
	vector<bool> mu_isTight;
	vector<bool> mu_isMedium;
	vector<bool> mu_isLoose;
	vector<bool> mu_isHighPt;
	vector<int> mu_isMatchedToGen;
	vector<int> mu_charge;
	vector<float> mu_dz;
	vector<float> mu_dxy;

	vector<float> jet_e;
	vector<float> jet_pt;
	vector<float> jet_eta;
	vector<float> jet_phi;
	vector<float> jet_bdiscriminant;
	vector<int>   jet_hadronFlavour;
	vector<int>   jet_partonFlavour;
	vector<int>   jet_isMatchedToGen;
	//jetID

	//variabili di DLDJcandidate
	vector<bool> isEEJJ;
	vector<bool> isEETT;
	vector<bool> isMMJJ;
	vector<bool> isMMTT;
	vector<bool> isEMJJ;

	vector<bool> isSignalRegion;
	vector<bool> isLowMllCR;
	vector<bool> isLowMlljjCR;

	vector<bool> isBB;
	vector<bool> isEE;
	vector<bool> isEB;

	vector<bool> passPreselections;

	vector<float> leadingLepton_e;
	vector<float> leadingLepton_pt;
	vector<float> leadingLepton_eta;  
	vector<float> leadingLepton_phi; 
	vector<float> leadingLepton_charge;

	vector<float> subLeadingLepton_e;
	vector<float> subLeadingLepton_pt;
	vector<float> subLeadingLepton_eta;  
	vector<float> subLeadingLepton_phi;  
	vector<float> subLeadingLepton_charge;

	vector<float> leadingJet_e;
	vector<float> leadingJet_pt;
	vector<float> leadingJet_eta;  
	vector<float> leadingJet_phi;  

	vector<float> subLeadingJet_e;
	vector<float> subLeadingJet_pt;
	vector<float> subLeadingJet_eta;  
	vector<float> subLeadingJet_phi;  

	vector<float> dRLeadLeptonLeadJet;
	vector<float> dRLeadLeptonSubLeadJet;
	vector<float> dRSubLeadLeptonLeadJet;
	vector<float> dRSubLeadLeptonSubLeadJet;

	vector<int> diLeptonDiJet_vtxIndex;
	vector<float> diLeptonDiJet_sumPt;
	vector<float> diLeptonDiJet_invMass;
	vector<float> diLepton_invMass;
	vector<float> diLepton_pt;
	vector<float> diJet_invMass;
	vector<float> diJetLeadingLepton_invMass;
	vector<float> diJetSubLeadingLepton_invMass;

	vector<bool> leadingEle_passHEEPId;
	vector<float> leadingEle_etaSC;
	vector<bool> leadingEle_isEcalDriven;
	vector<float> leadingEle_dEtaIn;  
	vector<float> leadingEle_dPhiIn;
	vector<float> leadingEle_hOverE;
	vector<float> leadingEle_full5x5_r9;  
	vector<float> leadingEle_full5x5_sigmaIetaIeta;
	vector<float> leadingEle_full5x5_E5x5;
	vector<float> leadingEle_full5x5_E1x5;
	vector<float> leadingEle_full5x5_E2x5;
	vector<float> leadingEle_full5x5_E2x5_Over_E5x5;
	vector<float> leadingEle_full5x5_E1x5_Over_E5x5;
	vector<float> leadingEle_EmHadDepth1Iso;
	vector<float> leadingEle_ptTracksIso;
	vector<int> leadingEle_innerLayerLostHits;
	vector<float> leadingEle_dxy;
	vector<float> leadingEle_eOverP;
	vector<unsigned> leadingEle_id;
	vector<bool> leadingEle_passCutBasedEleId;

	vector<bool> subLeadingEle_passHEEPId;		
	vector<float> subLeadingEle_etaSC;			
	vector<bool> subLeadingEle_isEcalDriven;
	vector<float> subLeadingEle_dEtaIn;  
	vector<float> subLeadingEle_dPhiIn;
	vector<float> subLeadingEle_hOverE;
	vector<float> subLeadingEle_full5x5_r9; 
	vector<float> subLeadingEle_full5x5_sigmaIetaIeta;
	vector<float> subLeadingEle_full5x5_E5x5;
	vector<float> subLeadingEle_full5x5_E1x5;
	vector<float> subLeadingEle_full5x5_E2x5;
	vector<float> subLeadingEle_full5x5_E2x5_Over_E5x5;
	vector<float> subLeadingEle_full5x5_E1x5_Over_E5x5;
	vector<float> subLeadingEle_EmHadDepth1Iso;	
	vector<float> subLeadingEle_ptTracksIso;	
	vector<int> subLeadingEle_innerLayerLostHits;
	vector<float> subLeadingEle_dxy;	
	vector<float> subLeadingEle_eOverP;
	vector<unsigned> subLeadingEle_id;
	vector<bool> subLeadingEle_passCutBasedEleId;

	vector<bool> leadingMuon_isHighPt;
	vector<bool> subLeadingMuon_isHighPt;

};



// // ************************** 
int electronMatchingToGen(Ptr<flashgg::Electron> electron,  Handle<View<reco::GenParticle> > genParticles){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
		if ( fabs(gen->pdgId()) != 11 ) continue;
		if ( !(gen->isPromptFinalState())) continue;
		float dR = deltaR( electron->eta(), electron->phi(), gen->eta(), gen->phi() );
		if (dR < 0.1){ //??? 0.1 ok???
			mcmatch = 1;
		}
	}
	return (mcmatch);
}
// ******************************************************************************************



// *********************** 
int muonMatchingToGen(Ptr<flashgg::Muon> muon, Handle<View<reco::GenParticle> > genParticles){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
		Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
			if ( fabs(gen->pdgId()) != 13 ) continue;
			if ( !(gen)->isPromptFinalState()) continue;
			float dR = deltaR( muon->eta(), muon->phi(), gen->eta(), gen->phi() );
			//cout << " *** Found muon:  ***"<<endl;
			//cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
			//cout << "dR = " << dR <<endl;
			if (dR < 0.1){ //??? 0.1 ok???
				mcmatch = 1;
			}
	}
	return (mcmatch);
}
// ******************************************************************************************



// // ************************** 
int jetMatchingToGen(flashgg::Jet jet,  Handle<View<reco::GenJet> > genJets){
	int mcmatch = 0;
	for( unsigned int i = 0 ; i < genJets->size(); i++ ) {
		Ptr<reco::GenJet> genJet = genJets->ptrAt(i);
		float dR = deltaR( jet.eta(), jet.phi(), genJet->eta(), genJet->phi() );
		if (dR > 0.4) continue;
		mcmatch = 1;
	}
	return (mcmatch);
}
// ******************************************************************************************



// **************** 
Ptr<reco::Vertex> chooseBestVtx(const vector<Ptr<reco::Vertex> > &vertices, Ptr<flashgg::Electron> electron){
	double vtx_dz = 1000000.;
	unsigned int min_dz_vtx = -1;			
	for( unsigned int vtxi = 0; vtxi < vertices.size(); vtxi++ ) {            
		Ptr<reco::Vertex> vtx = vertices[vtxi];            
		if( vtx_dz > fabs(electron->gsfTrack()->dz( vtx->position() )) ) {                
			vtx_dz = fabs(electron->gsfTrack()->dz( vtx->position() ) );
			min_dz_vtx = vtxi;
		}
	}					
	return vertices[min_dz_vtx];
}
// **************** 


// **************** 
Ptr<reco::Vertex> chooseBestVtx(const vector<Ptr<reco::Vertex> > &vertices, const flashgg::Electron* electron){
	double vtx_dz = 1000000.;
	unsigned int min_dz_vtx = -1;			
	for( unsigned int vtxi = 0; vtxi < vertices.size(); vtxi++ ) {            
		Ptr<reco::Vertex> vtx = vertices[vtxi];            
		if( vtx_dz > fabs(electron->gsfTrack()->dz( vtx->position() )) ) {                
			vtx_dz = fabs(electron->gsfTrack()->dz( vtx->position() ) );
			min_dz_vtx = vtxi;
		}
	}					
	return vertices[min_dz_vtx];
}
// **************** 



// **************** 
Ptr<reco::Vertex> chooseBestMuonVtx(const vector<Ptr<reco::Vertex> > &vertices, Ptr<flashgg::Muon> muon){
	int vtxInd = 0;
	double dzmin = 9999;
	for( size_t ivtx = 0 ; ivtx < vertices.size(); ivtx++ ) {
		Ptr<reco::Vertex> vtx = vertices[ivtx];
		if( !muon->innerTrack() ) { continue; }
		if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {
			dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
			vtxInd = ivtx;
		}
	}
	return vertices[vtxInd];
}
// **************** 


// **************** 
Ptr<reco::Vertex> chooseBestMuonVtx(const vector<Ptr<reco::Vertex> > &vertices, const flashgg::Muon* muon){
	int vtxInd = 0;
	double dzmin = 9999;
	for( size_t ivtx = 0 ; ivtx < vertices.size(); ivtx++ ) {
		Ptr<reco::Vertex> vtx = vertices[ivtx];
		if( !muon->innerTrack() ) { continue; }
		if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {
			dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
			vtxInd = ivtx;
		}
	}
	return vertices[vtxInd];
}
// **************** 



// **************** 
bool passHEEPIdCuts(Ptr<flashgg::Electron> electron, const vector<Ptr<reco::Vertex> > &pvPointers, double rho){

	bool pass = false;

	bool isEB = (fabs(electron->superCluster()->eta())<1.4442);
	bool isEE = (fabs(electron->superCluster()->eta())>1.566 && fabs(electron->superCluster()->eta())<2.5);

	float et = electron->et();
	bool isEcalDriven = electron->ecalDrivenSeed();
	float dEtaInSeed = electron->deltaEtaSuperClusterTrackAtVtx();  
	float dPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();  
	float hOverE = electron->hadronicOverEm(); 
	float full5x5_sigmaIetaIeta = electron->full5x5_sigmaIetaIeta();  
	float full5x5_E5x5 = electron->full5x5_e5x5(); 
	float full5x5_E1x5 = electron->full5x5_e1x5();
	float full5x5_E2x5 = electron->full5x5_e2x5Max(); 
	float full5x5_E2x5_Over_E5x5 = full5x5_E2x5 / full5x5_E5x5;
	float full5x5_E1x5_Over_E5x5 = full5x5_E1x5 / full5x5_E5x5;

	float EmHadDepth1Iso = electron->dr03EcalRecHitSumEt()+electron->dr03HcalDepth1TowerSumEt();
	float EmHadDepth1Iso_EB = 0.28*rho + 2 + 0.03*et;
	float EmHadDepth1Iso_EE = 0.28*rho + 2.5; 
	if (et > 50) EmHadDepth1Iso_EE = 0.28*rho + 2.5 + 0.03*(et-50);

	float ptTracksIso = 1.;//electron->dr03TkSumPt();  //TO MODIFY
	int innerLayerLostHits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
				
	Ptr<reco::Vertex> best_vtx = chooseBestVtx(pvPointers, electron);
	float dXY = fabs( electron->gsfTrack()->dxy( best_vtx->position()) ) ;  
	

	if (isEB) {
		if (et > 35 
			&& isEcalDriven 
			&& fabs(dEtaInSeed) < 0.004 
			&& fabs(dPhiIn) < 0.06
			&& hOverE < (1.0/electron->superCluster()->energy() + 0.05) 
			&& (full5x5_E2x5_Over_E5x5 > 0.94 || full5x5_E1x5_Over_E5x5 > 0.83) 
			&& EmHadDepth1Iso < EmHadDepth1Iso_EB 
			&& ptTracksIso < 5 
			&& innerLayerLostHits <=1 
			&& fabs(dXY) < 0.02 
		) pass =  true;
	}

	if (isEE) {
		if (et > 35 
			&& isEcalDriven 
			&& fabs(dEtaInSeed) < 0.006 
			&& fabs(dPhiIn) < 0.06
			&& hOverE < (5.0/electron->superCluster()->energy() + 0.05) 
			&& full5x5_sigmaIetaIeta < 0.03 
			&& EmHadDepth1Iso < EmHadDepth1Iso_EE 
			&& ptTracksIso < 5 
			&& innerLayerLostHits <=1 
			&& fabs(dXY) < 0.05
		) pass = true; 
	}

	return pass;

}
// ******************************************************************************************


// **************** 
bool passHEEPIdCuts(const flashgg::Electron* electron, const vector<Ptr<reco::Vertex> > &pvPointers, double rho){

	bool pass = false;

	bool isEB = (fabs(electron->superCluster()->eta())<1.4442);
	bool isEE = (fabs(electron->superCluster()->eta())>1.566 && fabs(electron->superCluster()->eta())<2.5);

	float et = electron->et();
	bool isEcalDriven = electron->ecalDrivenSeed();
	float dEtaInSeed = electron->deltaEtaSuperClusterTrackAtVtx();  //con electron->deltaEtaSeedClusterTrackAtVtx(); Product not found
	float dPhiIn = electron->deltaPhiSuperClusterTrackAtVtx(); 
	float hOverE = electron-> hadronicOverEm(); 
	float full5x5_sigmaIetaIeta = electron->full5x5_sigmaIetaIeta();  
	float full5x5_E5x5 = electron->full5x5_e5x5(); 
	float full5x5_E1x5 = electron->full5x5_e1x5();
	float full5x5_E2x5 = electron->full5x5_e2x5Max(); 
	float full5x5_E2x5_Over_E5x5 = full5x5_E2x5 / full5x5_E5x5;
	float full5x5_E1x5_Over_E5x5 = full5x5_E1x5 / full5x5_E5x5;

	float EmHadDepth1Iso = electron->dr03EcalRecHitSumEt()+electron->dr03HcalDepth1TowerSumEt();
	float EmHadDepth1Iso_EB = 0.28*rho + 2 + 0.03*et;
	float EmHadDepth1Iso_EE = 0.28*rho + 2.5; 
	if (et > 50) EmHadDepth1Iso_EE = 0.28*rho + 2.5 + 0.03*(et-50);

	float ptTracksIso = 1.;//electron->dr03TkSumPt();  //TO MODIFY
	int innerLayerLostHits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
				
	Ptr<reco::Vertex> best_vtx = chooseBestVtx(pvPointers, electron);
	float dXY = fabs( electron->gsfTrack()->dxy( best_vtx->position()) ) ;  
	

	if (isEB) {
		if (et > 35 
			&& isEcalDriven 
			&& fabs(dEtaInSeed) < 0.004 
			&& fabs(dPhiIn) < 0.06
			&& hOverE < (1.0/electron->superCluster()->energy() + 0.05) 
			&& (full5x5_E2x5_Over_E5x5 > 0.94 || full5x5_E1x5_Over_E5x5 > 0.83) 
			&& EmHadDepth1Iso < EmHadDepth1Iso_EB 
			&& ptTracksIso < 5 
			&& innerLayerLostHits <=1 
			&& fabs(dXY) < 0.02 
		) pass =  true;
	}

	if (isEE) {
		if (et > 35 
			&& isEcalDriven 
			&& fabs(dEtaInSeed) < 0.006 
			&& fabs(dPhiIn) < 0.06
			&& hOverE < (5.0/electron->superCluster()->energy() + 0.05) 
			&& full5x5_sigmaIetaIeta < 0.03 
			&& EmHadDepth1Iso < EmHadDepth1Iso_EE 
			&& ptTracksIso < 5 
			&& innerLayerLostHits <=1 
			&& fabs(dXY) < 0.05
		) pass = true; 
	}

	return pass;

}
// ******************************************************************************************



// ***************************** 
unsigned getEleId(Ptr<flashgg::Electron> electron){
	unsigned eleId = 0;
	std::vector<std::pair<std::string,float> > idlist = electron->electronIDs();
	for (unsigned i  = 0 ; i < idlist.size(); ++i){
		// cout << idlist[i].first << ": " << idlist[i].second << endl;
		if(int(idlist[i].second)){  //se valore id diverso da zero
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-veto") eleId |= 0;
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-loose") eleId |= 1;
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-medium") eleId |= 2;
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-tight") eleId |= 3;
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-loose") eleId |= 4;
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-medium") eleId |= 5;	    		
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-tight") eleId |= 6;
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-veto") eleId |= 7;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-Trig-V1-wp80") eleId |= 8;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-Trig-V1-wp90") eleId |= 9;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-nonTrig-V1-wp80") eleId |= 10;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-nonTrig-V1-wp90") eleId |= 11;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose") eleId |= 12;
			if (idlist[i].first == "heepElectronID-HEEPV60") eleId |= 13;
			if (idlist[i].first == "heepElectronID-HEEPV70") eleId |= 14;
		}	
	}
	return eleId;
}
// ******************************************************************************************



// ***************************** 
unsigned getEleId(const flashgg::Electron* electron){
	unsigned eleId = 0;
	std::vector<std::pair<std::string,float> > idlist = electron->electronIDs();
	for (unsigned i  = 0 ; i < idlist.size(); ++i){
		// cout << idlist[i].first << ": " << idlist[i].second << endl;
		if(int(idlist[i].second)){  //se valore id diverso da zero
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-veto") eleId |= 0;
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-loose") eleId |= 1;
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-medium") eleId |= 2;
			if (idlist[i].first == "cutBasedElectronID-Summer16-80X-V1-tight") eleId |= 3;
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-loose") eleId |= 4;
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-medium") eleId |= 5;	    		
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-tight") eleId |= 6;
			if (idlist[i].first == "cutBasedElectronID-Spring15-25ns-V1-standalone-veto") eleId |= 7;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-Trig-V1-wp80") eleId |= 8;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-Trig-V1-wp90") eleId |= 9;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-nonTrig-V1-wp80") eleId |= 10;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-nonTrig-V1-wp90") eleId |= 11;
			if (idlist[i].first == "mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose") eleId |= 12;
			if (idlist[i].first == "heepElectronID-HEEPV60") eleId |= 13;
			if (idlist[i].first == "heepElectronID-HEEPV70") eleId |= 14;
		}	
	}
	return eleId;
}
// ******************************************************************************************



// ***************************** 
float electronIsolation(Ptr<flashgg::Electron> electron, double rho){
	// -- compute combined relative isolation: IsoCh + max( 0.0, IsoNh + IsoPh - PU ) )/pT, PU = rho * Aeff 
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	// //effective areas:  https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
	// effective areas:  https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt

	float Aeff = 0;
	float eta = fabs(electron->eta());
	if( eta <  1.0 )                  { Aeff = 0.1703; }
	if( eta >= 1.0   && eta < 1.479 ) { Aeff = 0.1715; }
	if( eta >= 1.479 && eta < 2.0 )   { Aeff = 0.1213; }
	if( eta >= 2.0   && eta < 2.2 )   { Aeff = 0.1230; }
	if( eta >= 2.2   && eta < 2.3 )   { Aeff = 0.1635; }
	if( eta >= 2.3   && eta < 2.4 )   { Aeff = 0.1937; }
	if( eta >= 2.4 )                  { Aeff = 0.2393; }

	//float iso = electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ); 
	reco::GsfElectron::PflowIsolationVariables pfIso = electron->pfIsolationVariables();
	float iso = pfIso.sumChargedHadronPt + max( pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * Aeff, 0. );

	//cout << electron->chargedHadronIso() << "  " <<  pfIso.sumChargedHadronPt << "   pt = " << electron->pt() << endl; 
	//cout << electron->neutralHadronIso() << "  " << pfIso.sumNeutralHadronEt << endl;
	//cout << electron->photonIso() << "  " << pfIso.sumPhotonEt <<endl;
	//cout << electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ) << "   "<< iso<< endl;

	return (iso/ electron->pt());	
}
// ******************************************************************************************



// ***************************** 
float electronIsolation(const flashgg::Electron* electron, double rho){
	// -- compute combined relative isolation: IsoCh + max( 0.0, IsoNh + IsoPh - PU ) )/pT, PU = rho * Aeff 
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	// //effective areas:  https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
	// effective areas:  https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt

	float Aeff = 0;
	float eta = fabs(electron->eta());
	if( eta <  1.0 )                  { Aeff = 0.1703; }
	if( eta >= 1.0   && eta < 1.479 ) { Aeff = 0.1715; }
	if( eta >= 1.479 && eta < 2.0 )   { Aeff = 0.1213; }
	if( eta >= 2.0   && eta < 2.2 )   { Aeff = 0.1230; }
	if( eta >= 2.2   && eta < 2.3 )   { Aeff = 0.1635; }
	if( eta >= 2.3   && eta < 2.4 )   { Aeff = 0.1937; }
	if( eta >= 2.4 )                  { Aeff = 0.2393; }

	//float iso = electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ); 
	reco::GsfElectron::PflowIsolationVariables pfIso = electron->pfIsolationVariables();
	float iso = pfIso.sumChargedHadronPt + max( pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * Aeff, 0. );

	//cout << electron->chargedHadronIso() << "  " <<  pfIso.sumChargedHadronPt << "   pt = " << electron->pt() << endl; 
	//cout << electron->neutralHadronIso() << "  " << pfIso.sumNeutralHadronEt << endl;
	//cout << electron->photonIso() << "  " << pfIso.sumPhotonEt <<endl;
	//cout << electron->chargedHadronIso() + max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ) << "   "<< iso<< endl;

	return (iso/ electron->pt());	
}
// ******************************************************************************************



// **************** 
bool passCutBasedEleId(Ptr<flashgg::Electron> electron, double rho){
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	// Medium 80X-tuned selection

	bool pass = false;

	bool isEB = (fabs(electron->superCluster()->eta())<=1.479);

	float full5x5_sigmaIetaIeta = electron->full5x5_sigmaIetaIeta();  
	float dEtaInSeed = electron->deltaEtaSuperClusterTrackAtVtx(); 
	float dPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();  
	float hOverE = electron->hcalOverEcal();
	float pfIso = electronIsolation(electron, rho);

	float ooEmooP =-999 ; 
	if( electron->ecalEnergy() == 0 ){
	  ooEmooP = 1e30;
	}else if( !std::isfinite(electron->ecalEnergy())){    
	  ooEmooP = 1e30;
	}else{
	  ooEmooP = fabs(1.0/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy() );
	}

	int missingInnerHits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
	bool passConversionVeto = !(electron->hasMatchedConversion());  //non prende ele con conversioni in miniTree


	if (isEB) {
		if (full5x5_sigmaIetaIeta < 0.00998 
			&& fabs(dEtaInSeed) < 0.00311 
			&& fabs(dPhiIn) < 0.103 
			&& hOverE < 0.253 
			&& pfIso < 0.0695 
			&& ooEmooP < 0.134 
			&& missingInnerHits <=1 
			&& passConversionVeto
		) pass =  true;
	} else {
		if (full5x5_sigmaIetaIeta < 0.0298 
			&& fabs(dEtaInSeed) < 0.00609 
			&& fabs(dPhiIn) < 0.045 
			&& hOverE < 0.0878 
			&& pfIso < 0.0821 
			&& ooEmooP < 0.13 
			&& missingInnerHits <=1 
			&& passConversionVeto
		) pass = true; 
	}

	return pass;

}
// ******************************************************************************************




// **************** 
bool passCutBasedEleId(const flashgg::Electron* electron, double rho){
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	// Medium 80X-tuned selection

	bool pass = false;

	bool isEB = (fabs(electron->superCluster()->eta())<=1.479);

	float full5x5_sigmaIetaIeta = electron->full5x5_sigmaIetaIeta();  
	float dEtaInSeed = electron->deltaEtaSuperClusterTrackAtVtx(); 
	float dPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();  
	float hOverE = electron->hcalOverEcal();
	float pfIso = electronIsolation(electron, rho);

	float ooEmooP =-999 ; 
	if( electron->ecalEnergy() == 0 ){
	  ooEmooP = 1e30;
	}else if( !std::isfinite(electron->ecalEnergy())){    
	  ooEmooP = 1e30;
	}else{
	  ooEmooP = fabs(1.0/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy() );
	}

	int missingInnerHits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
	bool passConversionVeto = !(electron->hasMatchedConversion());  //non prende ele con conversioni in miniTree


	if (isEB) {
		if (full5x5_sigmaIetaIeta < 0.00998 
			&& fabs(dEtaInSeed) < 0.00311 
			&& fabs(dPhiIn) < 0.103 
			&& hOverE < 0.253 
			&& pfIso < 0.0695 
			&& ooEmooP < 0.134 
			&& missingInnerHits <=1 
			&& passConversionVeto
		) pass =  true;
	} else {
		if (full5x5_sigmaIetaIeta < 0.0298 
			&& fabs(dEtaInSeed) < 0.00609 
			&& fabs(dPhiIn) < 0.045 
			&& hOverE < 0.0878 
			&& pfIso < 0.0821 
			&& ooEmooP < 0.13 
			&& missingInnerHits <=1 
			&& passConversionVeto
		) pass = true; 
	}

	return pass;

}
// ******************************************************************************************



// ************************* 
bool passDiLeptonDiJetPreselection(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){

	if (dldj->leadingLeptonPt() < 60.) return false;
	if (dldj->subLeadingLeptonPt() < 53.) return false;
  	if (dldj->leadingJet()->pt() < 40.) return false;
	if (dldj->subLeadingJet()->pt() < 40.) return false;

	if (fabs(dldj->leadingLeptonEta()) > 2.4) return false;
	if (fabs(dldj->subLeadingLeptonEta()) > 2.4) return false;
	if (fabs(dldj->leadingJet()->eta()) > 2.4) return false;
	if (fabs(dldj->subLeadingJet()->eta()) > 2.4) return false;


	float dRLeadLeptonLeadJet = deltaR(dldj->leadingLeptonEta(), dldj->leadingLeptonPhi(), dldj->leadingJet()->eta(), dldj->leadingJet()->phi());
	float dRLeadLeptonSubLeadJet = deltaR(dldj->leadingLeptonEta(), dldj->leadingLeptonPhi(), dldj->subLeadingJet()->eta(), dldj->subLeadingJet()->phi());
	float dRSubLeadLeptonLeadJet = deltaR(dldj->subLeadingLeptonEta(), dldj->subLeadingLeptonPhi(), dldj->leadingJet()->eta(), dldj->leadingJet()->phi());
	float dRSubLeadLeptonSubLeadJet = deltaR(dldj->subLeadingLeptonEta(), dldj->subLeadingLeptonPhi(), dldj->subLeadingJet()->eta(), dldj->subLeadingJet()->phi());

	if( dRLeadLeptonLeadJet < 0.4 )  return false;
	if( dRLeadLeptonSubLeadJet < 0.4 )  return false;
	if( dRSubLeadLeptonLeadJet < 0.4 )  return false;
	if( dRSubLeadLeptonSubLeadJet < 0.4 )  return false;

	return true;
}
// ******************************************************************************************



// **************** 
bool isBB(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	if ( (fabs(dldj->leadingLeptonEta())<1.4442) 
		&& (fabs(dldj->subLeadingLeptonEta())<1.4442) 
	) return true;

	return false;
}
// **************** 


// **************** 
bool isEE(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	if ( (fabs(dldj->leadingLeptonEta())>1.566 && fabs(dldj->leadingLeptonEta())<2.5)  
		&& (fabs(dldj->subLeadingLeptonEta())>1.566 && fabs(dldj->subLeadingLeptonEta())<2.5)
	) return true;
	
	return false;
}
// **************** 


// **************** 
bool isEB(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	if ( ((fabs(dldj->leadingLeptonEta())<1.4442) && (fabs(dldj->subLeadingLeptonEta())>1.566 && fabs(dldj->subLeadingLeptonEta())<2.5)) 
		|| ((fabs(dldj->leadingLeptonEta())>1.566 && fabs(dldj->leadingLeptonEta())<2.5) && (fabs(dldj->subLeadingLeptonEta())<1.4442))
	) return true;

	return false;
}
// **************** 


// **************** 
bool isSignalRegion(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	if (dldj->invMass() > 600. && dldj->diLeptonInvMass() > 200.) {
		return true;
	}
	return false;
}
// **************** 


// **************** 
bool isLowMllCR(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	if (dldj->diLeptonInvMass() < 200.) return true;
	return false;
}
// **************** 


// **************** 
bool isLowMlljjCR(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	if (dldj->invMass() < 600. && dldj->diLeptonInvMass() > 200.) return true;
	return false;
}
// **************** 



// **************** 
float diLeptonPt(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	TLorentzVector l1, l2;  
	l1.SetPxPyPzE( 0., 0., 0., 0. );
	l2.SetPxPyPzE( 0., 0., 0., 0. );
	if (dldj->isEEJJ() || dldj->isEETT() ) {
		l1.SetPxPyPzE( dldj->electron1()->px(), dldj->electron1()->py(), dldj->electron1()->pz(), dldj->electron1()->energy() );
		l2.SetPxPyPzE( dldj->electron2()->px(), dldj->electron2()->py(), dldj->electron2()->pz(), dldj->electron2()->energy() );
	} else if (dldj->isMMJJ() || dldj->isMMTT()) {
		l1.SetPxPyPzE( dldj->muon1()->px(), dldj->muon1()->py(), dldj->muon1()->pz(), dldj->muon1()->energy() );
		l2.SetPxPyPzE( dldj->muon2()->px(), dldj->muon2()->py(), dldj->muon2()->pz(), dldj->muon2()->energy() );		
	} else if (dldj->isEMJJ()) {
		l1.SetPxPyPzE( dldj->electron()->px(), dldj->electron()->py(), dldj->electron()->pz(), dldj->electron()->energy() );
		l2.SetPxPyPzE( dldj->muon()->px(), dldj->muon()->py(), dldj->muon()->pz(), dldj->muon()->energy() );	
	}
	return (l1+l2).Pt();
}
// **************** 



// **************** 
float diJetInvMass(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	TLorentzVector j1, j2;  
	j1.SetPxPyPzE( 0., 0., 0., 0. );
	j2.SetPxPyPzE( 0., 0., 0., 0. );
	if (dldj->isEEJJ() || dldj->isMMJJ() || dldj->isEMJJ()) {
		j1.SetPxPyPzE( dldj->leadingJet()->px(), dldj->leadingJet()->py(), dldj->leadingJet()->pz(), dldj->leadingJet()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingJet()->px(), dldj->subLeadingJet()->py(), dldj->subLeadingJet()->pz(), dldj->subLeadingJet()->energy() );
	} else if (dldj->isEETT() || dldj->isMMTT()) {
		j1.SetPxPyPzE( dldj->leadingTrack()->px(), dldj->leadingTrack()->py(), dldj->leadingTrack()->pz(), dldj->leadingTrack()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingTrack()->px(), dldj->subLeadingTrack()->py(), dldj->subLeadingTrack()->pz(), dldj->subLeadingTrack()->energy() );
	}
	return (j1+j2).M();
}
// **************** 



// **************** 
float diJetLeadingLeptonInvMass(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	TLorentzVector j1, j2, l;  
	j1.SetPxPyPzE( 0., 0., 0., 0. );
	j2.SetPxPyPzE( 0., 0., 0., 0. );
	l.SetPxPyPzE( 0., 0., 0., 0. );

	if (dldj->isEEJJ()) {
		j1.SetPxPyPzE( dldj->leadingJet()->px(), dldj->leadingJet()->py(), dldj->leadingJet()->pz(), dldj->leadingJet()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingJet()->px(), dldj->subLeadingJet()->py(), dldj->subLeadingJet()->pz(), dldj->subLeadingJet()->energy() );
		l.SetPxPyPzE( dldj->leadingEle()->px(), dldj->leadingEle()->py(), dldj->leadingEle()->pz(), dldj->leadingEle()->energy() );
	} else if (dldj->isMMJJ()) {
		j1.SetPxPyPzE( dldj->leadingJet()->px(), dldj->leadingJet()->py(), dldj->leadingJet()->pz(), dldj->leadingJet()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingJet()->px(), dldj->subLeadingJet()->py(), dldj->subLeadingJet()->pz(), dldj->subLeadingJet()->energy() );
		l.SetPxPyPzE( dldj->leadingMuon()->px(), dldj->leadingMuon()->py(), dldj->leadingMuon()->pz(), dldj->leadingMuon()->energy() );
	} else if (dldj->isEMJJ()) {
		j1.SetPxPyPzE( dldj->leadingJet()->px(), dldj->leadingJet()->py(), dldj->leadingJet()->pz(), dldj->leadingJet()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingJet()->px(), dldj->subLeadingJet()->py(), dldj->subLeadingJet()->pz(), dldj->subLeadingJet()->energy() );
		if( dldj->electron()->pt() > dldj->muon()->pt() ) l.SetPxPyPzE( dldj->electron()->px(), dldj->electron()->py(), dldj->electron()->pz(), dldj->electron()->energy() );
		else l.SetPxPyPzE( dldj->muon()->px(), dldj->muon()->py(), dldj->muon()->pz(), dldj->muon()->energy() );
	} else if (dldj->isEETT()) {
		j1.SetPxPyPzE( dldj->leadingTrack()->px(), dldj->leadingTrack()->py(), dldj->leadingTrack()->pz(), dldj->leadingTrack()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingTrack()->px(), dldj->subLeadingTrack()->py(), dldj->subLeadingTrack()->pz(), dldj->subLeadingTrack()->energy() );
		l.SetPxPyPzE( dldj->leadingEle()->px(), dldj->leadingEle()->py(), dldj->leadingEle()->pz(), dldj->leadingEle()->energy() );
	} else if (dldj->isMMTT()) {
		j1.SetPxPyPzE( dldj->leadingTrack()->px(), dldj->leadingTrack()->py(), dldj->leadingTrack()->pz(), dldj->leadingTrack()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingTrack()->px(), dldj->subLeadingTrack()->py(), dldj->subLeadingTrack()->pz(), dldj->subLeadingTrack()->energy() );
		l.SetPxPyPzE( dldj->leadingMuon()->px(), dldj->leadingMuon()->py(), dldj->leadingMuon()->pz(), dldj->leadingMuon()->energy() );
	}

	return (j1+j2+l).M();
}
// **************** 



// **************** 
float diJetSubLeadingLeptonInvMass(Ptr<flashgg::DiLeptonDiJetCandidate> dldj){
	TLorentzVector j1, j2, l;  
	j1.SetPxPyPzE( 0., 0., 0., 0. );
	j2.SetPxPyPzE( 0., 0., 0., 0. );
	l.SetPxPyPzE( 0., 0., 0., 0. );

	if (dldj->isEEJJ()) {
		j1.SetPxPyPzE( dldj->leadingJet()->px(), dldj->leadingJet()->py(), dldj->leadingJet()->pz(), dldj->leadingJet()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingJet()->px(), dldj->subLeadingJet()->py(), dldj->subLeadingJet()->pz(), dldj->subLeadingJet()->energy() );
		l.SetPxPyPzE( dldj->subLeadingEle()->px(), dldj->subLeadingEle()->py(), dldj->subLeadingEle()->pz(), dldj->subLeadingEle()->energy() );
	} else if (dldj->isMMJJ()) {
		j1.SetPxPyPzE( dldj->leadingJet()->px(), dldj->leadingJet()->py(), dldj->leadingJet()->pz(), dldj->leadingJet()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingJet()->px(), dldj->subLeadingJet()->py(), dldj->subLeadingJet()->pz(), dldj->subLeadingJet()->energy() );
		l.SetPxPyPzE( dldj->subLeadingMuon()->px(), dldj->subLeadingMuon()->py(), dldj->subLeadingMuon()->pz(), dldj->subLeadingMuon()->energy() );
	} else if (dldj->isEMJJ()) {
		j1.SetPxPyPzE( dldj->leadingJet()->px(), dldj->leadingJet()->py(), dldj->leadingJet()->pz(), dldj->leadingJet()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingJet()->px(), dldj->subLeadingJet()->py(), dldj->subLeadingJet()->pz(), dldj->subLeadingJet()->energy() );
		if( dldj->electron()->pt() < dldj->muon()->pt() ) l.SetPxPyPzE( dldj->electron()->px(), dldj->electron()->py(), dldj->electron()->pz(), dldj->electron()->energy() );
		else l.SetPxPyPzE( dldj->muon()->px(), dldj->muon()->py(), dldj->muon()->pz(), dldj->muon()->energy() );
	} else if (dldj->isEETT()) {
		j1.SetPxPyPzE( dldj->leadingTrack()->px(), dldj->leadingTrack()->py(), dldj->leadingTrack()->pz(), dldj->leadingTrack()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingTrack()->px(), dldj->subLeadingTrack()->py(), dldj->subLeadingTrack()->pz(), dldj->subLeadingTrack()->energy() );
		l.SetPxPyPzE( dldj->subLeadingEle()->px(), dldj->subLeadingEle()->py(), dldj->subLeadingEle()->pz(), dldj->subLeadingEle()->energy() );
	} else if (dldj->isMMTT()) {
		j1.SetPxPyPzE( dldj->leadingTrack()->px(), dldj->leadingTrack()->py(), dldj->leadingTrack()->pz(), dldj->leadingTrack()->energy() );
		j2.SetPxPyPzE( dldj->subLeadingTrack()->px(), dldj->subLeadingTrack()->py(), dldj->subLeadingTrack()->pz(), dldj->subLeadingTrack()->energy() );
		l.SetPxPyPzE( dldj->subLeadingMuon()->px(), dldj->subLeadingMuon()->py(), dldj->subLeadingMuon()->pz(), dldj->subLeadingMuon()->energy() );
	}

	return (j1+j2+l).M();
}
// **************** 




class miniTreeMaker : public BasicAnalyzer 
{
 public:
	miniTreeMaker( const ParameterSet & iConfig, TFileDirectory& fs, ConsumesCollector && cc);
	virtual ~miniTreeMaker();
	void beginJob();
	void analyze( const EventBase& event );
	void endJob();
 
 private:
	void initEventStructure();

	TTree *eventTree;
	eventInfo evInfo;
	int ngen;
	int ngenPre;
	int ndldj;
	int npre;
	int nEle;
	int nEleGood;
	int nElePassingHEEPid;
	int nmuons;
	int nmuonsGood;

	EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
	EDGetTokenT<GenEventInfoProduct> genInfoToken_;
	EDGetTokenT<View<PileupSummaryInfo> >  PileUpToken_;
	EDGetTokenT<View<reco::Vertex> > vertexToken_;
	EDGetTokenT<View<flashgg::DiLeptonDiJetCandidate> > DiLeptonDiJetToken_; 
	EDGetTokenT<View<vector<flashgg::Jet> > > jetsToken_;
	EDGetTokenT<View<reco::GenJet> > genJetToken_;
	EDGetTokenT<View<Electron> > electronToken_;
	EDGetTokenT<View<Muon> > muonToken_;
	EDGetTokenT<TriggerResults> triggerBitsToken_;
	EDGetTokenT<double> rhoToken_;
	vector<EDGetTokenT<View<flashgg::Jet> > > tokenJets_;

	typedef vector<Handle<View<flashgg::Jet> > > JetCollectionVector;

	string bTag_;
	double lumiWeight_;

	GlobalVariablesDumper *globalVarsDumper_;
};
// ******************************************************************************************                                                                                                                     


