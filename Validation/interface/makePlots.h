#ifndef _makePlots_h_
#define _makePlots_h_
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <vector>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include "dafne/Validation/interface/functions.h"


using namespace std;


class makePlots {
	public:
	// TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		TChain          *fChain; 
		Int_t           fCurrent; //!current Tree number in a TChain

	// Declaration of leaf types
		Int_t           run;
		Int_t           event;
		Int_t           lumi;
		Float_t         rho;
		Float_t         weight;
		Float_t         puweight;
		Int_t           nvtx;
		Int_t           npu;
		Int_t           passEEJJhlt;
		Int_t           passMMJJhlt;
		Int_t           passEMJJhlt;
		Int_t           passTandPEEhlt;
		Int_t           passTandPMMhlt;
		vector<float>   *vtx_x;
		vector<float>   *vtx_y;
		vector<float>   *vtx_z;     
		vector<float>   *ele_e;
		vector<float>   *ele_pt;
		vector<float>   *ele_eta;
		vector<float>   *ele_phi;
		vector<float>   *ele_idmva;
		vector<unsigned int> *ele_id;
		vector<float>   *ele_iso;
		vector<float>   *ele_dz;
		vector<float>   *ele_d0;
		vector<bool>    *ele_passHEEPId;
		vector<bool>    *ele_passCutBasedEleId;
		vector<int>     *ele_isMatchedToGen;
		vector<int>     *ele_charge;
		vector<float>   *ele_etaSC;
		vector<bool>    *ele_isEcalDriven;
		vector<float>   *ele_dEtaIn;
		vector<float>   *ele_dPhiIn;
		vector<float>   *ele_hOverE;
		vector<float>   *ele_full5x5_r9;
		vector<float>   *ele_full5x5_sigmaIetaIeta;
		vector<float>   *ele_full5x5_E5x5;
		vector<float>   *ele_full5x5_E1x5;
		vector<float>   *ele_full5x5_E2x5;
		vector<float>   *ele_full5x5_E2x5_Over_E5x5;
		vector<float>   *ele_full5x5_E1x5_Over_E5x5;
		vector<float>   *ele_EmHadDepth1Iso;
		vector<float>   *ele_ptTracksIso;
		vector<int>     *ele_innerLayerLostHits;
		vector<float>   *ele_dxy;
		vector<float>   *ele_eOverP;
		vector<float>   *mu_e;
		vector<float>   *mu_pt;
		vector<float>   *mu_eta;
		vector<float>   *mu_phi;
		vector<float>   *mu_iso;
		vector<bool>    *mu_isTight;
		vector<bool>    *mu_isMedium;
		vector<bool>    *mu_isLoose;
		vector<bool>    *mu_isHighPt;
		vector<int>     *mu_isMatchedToGen;
		vector<int>     *mu_charge;
		vector<float>   *mu_dz;
		vector<float>   *mu_dxy;
		vector<float>   *jet_e;
		vector<float>   *jet_pt;
		vector<float>   *jet_eta;
		vector<float>   *jet_phi;
		vector<float>   *jet_bdiscriminant;
		vector<int>     *jet_partonFlavour;
		vector<int>     *jet_hadronFlavour;
		vector<int>     *jet_isMatchedToGen;
		vector<bool>    *isEEJJ;
		vector<bool>    *isEETT;
		vector<bool>    *isMMJJ;
		vector<bool>    *isMMTT;
		vector<bool>    *isEMJJ;
		vector<bool>    *isSignalRegion;
		vector<bool>    *isLowMllCR;
		vector<bool>    *isLowMlljjCR;
		vector<bool>    *isBB;
		vector<bool>    *isEE;
		vector<bool>    *isEB;
		vector<bool>    *passPreselections;
		vector<float>   *leadingLepton_e;
		vector<float>   *leadingLepton_pt;
		vector<float>   *leadingLepton_eta;
		vector<float>   *leadingLepton_phi;
		vector<float>   *leadingLepton_charge;
		vector<float>   *subLeadingLepton_e;
		vector<float>   *subLeadingLepton_pt;
		vector<float>   *subLeadingLepton_eta;
		vector<float>   *subLeadingLepton_phi;
		vector<float>   *subLeadingLepton_charge;
		vector<float>   *leadingJet_e;
		vector<float>   *leadingJet_pt;
		vector<float>   *leadingJet_eta;
		vector<float>   *leadingJet_phi;
		vector<float>   *subLeadingJet_e;
		vector<float>   *subLeadingJet_pt;
		vector<float>   *subLeadingJet_eta;
		vector<float>   *subLeadingJet_phi;
		vector<float>   *dRLeadLeptonLeadJet;
		vector<float>   *dRLeadLeptonSubLeadJet;
		vector<float>   *dRSubLeadLeptonLeadJet;
		vector<float>   *dRSubLeadLeptonSubLeadJet;
		vector<int>     *diLeptonDiJet_vtxIndex;
		vector<float>   *diLeptonDiJet_sumPt;
		vector<float>   *diLeptonDiJet_invMass;
		vector<float>   *diLepton_invMass;
		vector<float>   *diJet_invMass;
		vector<float>   *diJetLeadingLepton_invMass;
		vector<float>   *diJetSubLeadingLepton_invMass;
		vector<bool>    *leadingEle_passHEEPId;
		vector<float>   *leadingEle_etaSC;
		vector<bool>    *leadingEle_isEcalDriven;
		vector<float>   *leadingEle_dEtaIn;
		vector<float>   *leadingEle_dPhiIn;
		vector<float>   *leadingEle_hOverE;
		vector<float>   *leadingEle_full5x5_r9;
		vector<float>   *leadingEle_full5x5_sigmaIetaIeta;
		vector<float>   *leadingEle_full5x5_E5x5;
		vector<float>   *leadingEle_full5x5_E1x5;
		vector<float>   *leadingEle_full5x5_E2x5;
		vector<float>   *leadingEle_full5x5_E2x5_Over_E5x5;
		vector<float>   *leadingEle_full5x5_E1x5_Over_E5x5;
		vector<float>   *leadingEle_EmHadDepth1Iso;
		vector<float>   *leadingEle_ptTracksIso;
		vector<int>     *leadingEle_innerLayerLostHits;
		vector<float>   *leadingEle_dxy;
		vector<float>   *leadingEle_eOverP;
		vector<unsigned int> *leadingEle_id;
		vector<bool>    *leadingEle_passCutBasedEleId;
		vector<bool>    *subLeadingEle_passHEEPId;
		vector<float>   *subLeadingEle_etaSC;
		vector<bool>    *subLeadingEle_isEcalDriven;
		vector<float>   *subLeadingEle_dEtaIn;
		vector<float>   *subLeadingEle_dPhiIn;
		vector<float>   *subLeadingEle_hOverE;
		vector<float>   *subLeadingEle_full5x5_r9;
		vector<float>   *subLeadingEle_full5x5_sigmaIetaIeta;
		vector<float>   *subLeadingEle_full5x5_E5x5;
		vector<float>   *subLeadingEle_full5x5_E1x5;
		vector<float>   *subLeadingEle_full5x5_E2x5;
		vector<float>   *subLeadingEle_full5x5_E2x5_Over_E5x5;
		vector<float>   *subLeadingEle_full5x5_E1x5_Over_E5x5;
		vector<float>   *subLeadingEle_EmHadDepth1Iso;
		vector<float>   *subLeadingEle_ptTracksIso;
		vector<int>     *subLeadingEle_innerLayerLostHits;
		vector<float>   *subLeadingEle_dxy; 
		vector<float>   *subLeadingEle_eOverP;
		vector<unsigned int> *subLeadingEle_id;
		vector<bool>    *subLeadingEle_passCutBasedEleId;
		vector<bool>    *leadingMuon_isHighPt;
		vector<bool>    *subLeadingMuon_isHighPt;

	// List of branches
		TBranch        *b_run;   //!
		TBranch        *b_event;   //!
		TBranch        *b_lumi;   //!
		TBranch        *b_rho;   //!
		TBranch        *b_weight;   //!
		TBranch        *b_puweight;   //!
		TBranch        *b_nvtx;   //!
		TBranch        *b_npu;   //!
		TBranch        *b_passEEJJhlt;   //!
		TBranch        *b_passMMJJhlt;   //!
		TBranch        *b_passEMJJhlt;   //!
		TBranch        *b_passTandPEEhlt;   //!
		TBranch        *b_passTandPMMhlt;   //!
		TBranch        *b_vtx_x;   //!
		TBranch        *b_vtx_y;   //!
		TBranch        *b_vtx_z;   //!
		TBranch        *b_ele_e;   //!
		TBranch        *b_ele_pt;   //!
		TBranch        *b_ele_eta;   //!
		TBranch        *b_ele_phi;   //!
		TBranch        *b_ele_idmva;   //!
		TBranch        *b_ele_id;   //!
		TBranch        *b_ele_iso;   //!
		TBranch        *b_ele_dz;   //!
		TBranch        *b_ele_d0;   //!
		TBranch        *b_ele_passHEEPId;   //!
		TBranch        *b_ele_passCutBasedEleId;   //!
		TBranch        *b_ele_isMatchedToGen;   //!
		TBranch        *b_ele_charge;   //!
		TBranch        *b_ele_etaSC;   //!
		TBranch        *b_ele_isEcalDriven;   //!
		TBranch        *b_ele_dEtaIn;   //!
		TBranch        *b_ele_dPhiIn;   //!
		TBranch        *b_ele_hOverE;   //!
		TBranch        *b_ele_full5x5_r9;   //!
		TBranch        *b_ele_full5x5_sigmaIetaIeta;   //!
		TBranch        *b_ele_full5x5_E5x5;   //!
		TBranch        *b_ele_full5x5_E1x5;   //!
		TBranch        *b_ele_full5x5_E2x5;   //!
		TBranch        *b_ele_full5x5_E2x5_Over_E5x5;   //!
		TBranch        *b_ele_full5x5_E1x5_Over_E5x5;   //!
		TBranch        *b_ele_EmHadDepth1Iso;   //!
		TBranch        *b_ele_ptTracksIso;   //!
		TBranch        *b_ele_innerLayerLostHits;   //!
		TBranch        *b_ele_dxy;   //!
		TBranch        *b_ele_eOverP;   //!
		TBranch        *b_mu_e;   //!
		TBranch        *b_mu_pt;   //!
		TBranch        *b_mu_eta;   //!
		TBranch        *b_mu_phi;   //!
		TBranch        *b_mu_iso;   //!
		TBranch        *b_mu_isTight;   //!
		TBranch        *b_mu_isMedium;   //!
		TBranch        *b_mu_isLoose;   //!
		TBranch        *b_mu_isHighPt;   //!
		TBranch        *b_mu_isMatchedToGen;   //!
		TBranch        *b_mu_charge;   //!
		TBranch        *b_mu_dz;   //!
		TBranch        *b_mu_dxy;   //!
		TBranch        *b_jet_e;   //!
		TBranch        *b_jet_pt;   //!
		TBranch        *b_jet_eta;   //!
		TBranch        *b_jet_phi;   //!
		TBranch        *b_jet_bdiscriminant;   //!
		TBranch        *b_jet_partonFlavour;   //!
		TBranch        *b_jet_hadronFlavour;   //!
		TBranch        *b_jet_isMatchedToGen;   //!
		TBranch        *b_isEEJJ;   //!
		TBranch        *b_isEETT;   //!
		TBranch        *b_isMMJJ;   //!
		TBranch        *b_isMMTT;   //!
		TBranch        *b_isEMJJ;   //!
		TBranch        *b_isSignalRegion;   //!
		TBranch        *b_isLowMllCR;   //!
		TBranch        *b_isLowMlljjCR;   //!
		TBranch        *b_isBB;   //!
		TBranch        *b_isEE;   //!
		TBranch        *b_isEB;   //!
		TBranch        *b_passPreselections;   //!
		TBranch        *b_leadingLepton_e;   //!
		TBranch        *b_leadingLepton_pt;   //!
		TBranch        *b_leadingLepton_eta;   //!
		TBranch        *b_leadingLepton_phi;   //!
		TBranch        *b_leadingLepton_charge;   //!
		TBranch        *b_subLeadingLepton_e;   //!
		TBranch        *b_subLeadingLepton_pt;   //!
		TBranch        *b_subLeadingLepton_eta;   //!
		TBranch        *b_subLeadingLepton_phi;   //!
		TBranch        *b_subLeadingLepton_charge;   //!
		TBranch        *b_leadingJet_e;   //!
		TBranch        *b_leadingJet_pt;   //!
		TBranch        *b_leadingJet_eta;   //!
		TBranch        *b_leadingJet_phi;   //!
		TBranch        *b_subLeadingJet_e;   //!
		TBranch        *b_subLeadingJet_pt;   //!
		TBranch        *b_subLeadingJet_eta;   //!
		TBranch        *b_subLeadingJet_phi;   //!
		TBranch        *b_dRLeadLeptonLeadJet;   //!
		TBranch        *b_dRLeadLeptonSubLeadJet;   //!
		TBranch        *b_dRSubLeadLeptonLeadJet;   //!
		TBranch        *b_dRSubLeadLeptonSubLeadJet;   //!
		TBranch        *b_diLeptonDiJet_vtxIndex;   //!
		TBranch        *b_diLeptonDiJet_sumPt;   //!
		TBranch        *b_diLeptonDiJet_invMass;   //!
		TBranch        *b_diLepton_invMass;   //!
		TBranch        *b_diJet_invMass;   //!
		TBranch        *b_diJetLeadingLepton_invMass;   //!
		TBranch        *b_diJetSubLeadingLepton_invMass;   //!
		TBranch        *b_leadingEle_passHEEPId;   //!
		TBranch        *b_leadingEle_etaSC;   //!
		TBranch        *b_leadingEle_isEcalDriven;   //!
		TBranch        *b_leadingEle_dEtaIn;   //!
		TBranch        *b_leadingEle_dPhiIn;   //!
		TBranch        *b_leadingEle_hOverE;   //!
		TBranch        *b_leadingEle_full5x5_r9;   //!
		TBranch        *b_leadingEle_full5x5_sigmaIetaIeta;   //!
		TBranch        *b_leadingEle_full5x5_E5x5;   //!
		TBranch        *b_leadingEle_full5x5_E1x5;   //!
		TBranch        *b_leadingEle_full5x5_E2x5;   //!
		TBranch        *b_leadingEle_full5x5_E2x5_Over_E5x5;   //!
		TBranch        *b_leadingEle_full5x5_E1x5_Over_E5x5;   //!
		TBranch        *b_leadingEle_EmHadDepth1Iso;   //!
		TBranch        *b_leadingEle_ptTracksIso;   //!
		TBranch        *b_leadingEle_innerLayerLostHits;   //!
		TBranch        *b_leadingEle_dxy;   //!
		TBranch        *b_leadingEle_eOverP;   //!
		TBranch        *b_leadingEle_id;   //!
		TBranch        *b_leadingEle_passCutBasedEleId;   //!
		TBranch        *b_subLeadingEle_passHEEPId;   //!
		TBranch        *b_subLeadingEle_etaSC;   //!
		TBranch        *b_subLeadingEle_isEcalDriven;   //!
		TBranch        *b_subLeadingEle_dEtaIn;   //!
		TBranch        *b_subLeadingEle_dPhiIn;   //!
		TBranch        *b_subLeadingEle_hOverE;   //!
		TBranch        *b_subLeadingEle_full5x5_r9;   //!
		TBranch        *b_subLeadingEle_full5x5_sigmaIetaIeta;   //!
		TBranch        *b_subLeadingEle_full5x5_E5x5;   //!
		TBranch        *b_subLeadingEle_full5x5_E1x5;   //!
		TBranch        *b_subLeadingEle_full5x5_E2x5;   //!
		TBranch        *b_subLeadingEle_full5x5_E2x5_Over_E5x5;   //!
		TBranch        *b_subLeadingEle_full5x5_E1x5_Over_E5x5;   //!
		TBranch        *b_subLeadingEle_EmHadDepth1Iso;   //!
		TBranch        *b_subLeadingEle_ptTracksIso;   //!
		TBranch        *b_subLeadingEle_innerLayerLostHits;   //!
		TBranch        *b_subLeadingEle_dxy;   //!
		TBranch        *b_subLeadingEle_eOverP;   //!
		TBranch        *b_subLeadingEle_id;   //!
		TBranch        *b_subLeadingEle_passCutBasedEleId;   //!
		TBranch        *b_leadingMuon_isHighPt;   //!
		TBranch        *b_subLeadingMuon_isHighPt;   //!

	//definition of histos
		TH1D *rho_histo;
		TH1D *nvtx_histo;

		TH1D *vtx_x_histo;
		TH1D *vtx_y_histo;
		TH1D *vtx_z_histo;

		TH1D *muon_dxy_histo;

		TH1D *pt_histo[9];
		TH1D *eta_histo[9];
		TH1D *phi_histo[9];

		TH1D *etaSC_histo[9][3];
		TH1D *isEcalDriven_histo[9][3];  
		TH1D *dEtaIn_histo[9][3];
		TH1D *dPhiIn_histo[9][3];
		TH1D *hOverE_histo[9][3];
		TH1D *full5x5_r9_histo[9][3];
		TH1D *full5x5_sigmaIetaIeta_histo[9][3];
		TH1D *full5x5_E5x5_histo[9][3];
		TH1D *full5x5_E1x5_histo[9][3];
		TH1D *full5x5_E2x5_histo[9][3];
		TH1D *full5x5_E2x5_Over_E5x5_histo[9][3];
		TH1D *full5x5_E1x5_Over_E5x5_histo[9][3];
		TH1D *EmHadDepth1Iso_histo[9][3];
		TH1D *ptTracksIso_histo[9][3];
		TH1D *innerLayerLostHits_histo[9][3];
		TH1D *dxy_histo[9][3];
		TH1D *eOverP[9][3];

		TH1D *mass_dldj_histo[4][5];
		TH1D *mass_dl_histo[4][5];
		TH1D *mass_dj_histo[4][5];
		TH1D *mass_djLl_histo[4][5];
		TH1D *mass_djSLl_histo[4][5];
		TH1D *Zmass_histo[4][5];



		makePlots(TString filename_, TString outputdir_, bool MC_,  bool MCpuReweighted_, bool signalEE_, bool signalMuMu_, bool eMuSideband_, bool TnPee_, bool TnPmumu_);

		void     Init();  
		Int_t    GetEntry(Long64_t entry);
		Long64_t LoadTree(Long64_t entry);

		TH1D* newTH1D(string, string, string, int, double, double);
		TH2D* newTH2D(string, string, string, string, int, double, double, int, double, double);
		TProfile* newTProfile(string, string, string, string, int, double, double, double, double);
		void SetHistos();

		void getElectrons(vector<eleStruct>& electrons);
		void getLeadingElectrons(vector<eleStruct>& leadingElectrons);
		void getSubLeadingElectrons(vector<eleStruct>& subLeadingElectrons);

		void getMuons(vector<muonStruct>& muons);
		void getLeadingMuons(vector<muonStruct>& leadingMuons);
		void getSubLeadingMuons(vector<muonStruct>& subLeadingMuons);

		void getJets(vector<jetStruct>& jets);
		void getLeadingJets(vector<jetStruct>& leadingJets);
		void getSubLeadingJets(vector<jetStruct>& subLeadingJets);

		void getDiLeptonDiJets(vector<diLeptonDiJetStruct>& diLeptonDiJets);

		void doEleDistributionsPlots(vector<eleStruct>& electrons, const int eleIdx, const int histoIdx);
		void doMuonDistributionsPlots(vector<muonStruct>& muons, const int muIdx, const int histoIdx);
		void doJetsDistributionsPlots(vector<jetStruct>& jets, const int jetIdx, const int histoIdx);
		
		void doElePlots(vector<eleStruct>& electrons, const int eleIdx, const int histoIdx, const int etaIdx);

		void doMassPlots(vector<diLeptonDiJetStruct>& diLeptonDiJets, const int dldjIdx, const int histoIdx, const int etaIdx);
		void doZmassPlots(vector<diLeptonDiJetStruct>& diLeptonDiJets, const int dldjIdx, const int histoIdx, const int etaIdx);

		void Loop();
		void saveHistosAndOutputFile(TString& outputdir);

		TString filename, outputdir; 
		bool MC, MCpuReweighted, signalEE, signalMuMu, eMuSideband, TnPee, TnPmumu;
		vector<TH1*> listOfHistograms;

		unsigned int nEvents=0, nEventsPassingTrigger=0;

		float w=1,lumiData=1;

		string suff = "";  

		string etaName[3] = {"", "_EB", "_EE"};
		string etaMassName[5] = {"", "_EB-EB", "_EE-EE", "_EB-EE", "_noEB-EB"};

		string eleName[9] = {"ele", "leadingEle", "subLeadingEle", "elePassingHeepId", "leadingElePassingHeepId", "subLeadingElePassingHeepId", "elePassingEleId", "leadingElePassingEleId", "subLeadingElePassingEleId"};
		string objName[9] = {"ele", "leadingEle", "subLeadingEle", "muons", "leadingMu", "subLeadingMu", "jet", "leadingJet", "subLeadingJet"};
		string dldjName[4] = {"dldj","dldjPreselected","dldjPassingHeepId","dldjPassingEleId"};
		string Zname[4] = {"toEE", "toMuMu", "toEEpassingHeepId", "toEEpassingEleId"};

};

#endif

