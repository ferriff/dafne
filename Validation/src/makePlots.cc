#include "dafne/Validation/interface/makePlots.h"
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <limits>
#include <iostream>
#include <sstream>
// #include "tdrstyle.C"
// #include "CMS_lumi.C"

using namespace std;



// **************** 
void makePlots::Init(){

	// Set object pointer
   	vtx_x = 0;
   	vtx_y = 0;
   	vtx_z = 0;
	ele_e = 0;
	ele_pt = 0;
	ele_eta = 0;
	ele_phi = 0;
	ele_idmva = 0;
	ele_id = 0;
	ele_iso = 0;
	ele_dz = 0;
	ele_d0 = 0;
	ele_passHEEPId = 0;
	ele_passCutBasedEleId = 0;
	ele_isMatchedToGen = 0;
	ele_charge = 0;
	ele_etaSC = 0;
	ele_isEcalDriven = 0;
	ele_dEtaIn = 0;
	ele_dPhiIn = 0;
	ele_hOverE = 0;
	ele_full5x5_r9 = 0;
	ele_full5x5_sigmaIetaIeta = 0;
	ele_full5x5_E5x5 = 0;
	ele_full5x5_E1x5 = 0;
	ele_full5x5_E2x5 = 0;
	ele_full5x5_E2x5_Over_E5x5 = 0;
	ele_full5x5_E1x5_Over_E5x5 = 0;
	ele_EmHadDepth1Iso = 0;
	ele_ptTracksIso = 0;
	ele_innerLayerLostHits = 0;
	ele_dxy = 0;
	ele_eOverP = 0;
	mu_e = 0;
	mu_pt = 0;
	mu_eta = 0;
	mu_phi = 0;
	mu_iso = 0;
	mu_isTight = 0;
	mu_isMedium = 0;
	mu_isLoose = 0;
	mu_isHighPt = 0;
	mu_isMatchedToGen = 0;
	mu_charge = 0;
	mu_dz = 0;
	mu_dxy = 0;
	jet_e = 0;
	jet_pt = 0;
	jet_eta = 0;
	jet_phi = 0;
	jet_bdiscriminant = 0;
	jet_partonFlavour = 0;
	jet_hadronFlavour = 0;
	jet_isMatchedToGen = 0;
	isEEJJ = 0;
	isEETT = 0;
	isMMJJ = 0;
	isMMTT = 0;
	isEMJJ = 0;
	isSignalRegion = 0;
	isLowMllCR = 0;
	isLowMlljjCR = 0;
	isBB = 0;
	isEE = 0;
	isEB = 0;
	passPreselections = 0;
	leadingLepton_e = 0;
	leadingLepton_pt = 0;
	leadingLepton_eta = 0;
	leadingLepton_phi = 0;
	leadingLepton_charge = 0;
	subLeadingLepton_e = 0;
	subLeadingLepton_pt = 0;
	subLeadingLepton_eta = 0;
	subLeadingLepton_phi = 0;
	subLeadingLepton_charge = 0;
	leadingJet_e = 0;
	leadingJet_pt = 0;
	leadingJet_eta = 0;
	leadingJet_phi = 0;
	subLeadingJet_e = 0;
	subLeadingJet_pt = 0;
	subLeadingJet_eta = 0;
	subLeadingJet_phi = 0;
	dRLeadLeptonLeadJet = 0;
	dRLeadLeptonSubLeadJet = 0;
	dRSubLeadLeptonLeadJet = 0;
	dRSubLeadLeptonSubLeadJet = 0;
	diLeptonDiJet_vtxIndex = 0;
	diLeptonDiJet_sumPt = 0;
	diLeptonDiJet_invMass = 0;
	diLepton_invMass = 0;
	diJet_invMass = 0;
	diJetLeadingLepton_invMass = 0;
	diJetSubLeadingLepton_invMass = 0;
	leadingEle_passHEEPId = 0;
	leadingEle_etaSC = 0;
	leadingEle_isEcalDriven = 0;
	leadingEle_dEtaIn = 0;
	leadingEle_dPhiIn = 0;
	leadingEle_hOverE = 0;
	leadingEle_full5x5_r9 = 0;
	leadingEle_full5x5_sigmaIetaIeta = 0;
	leadingEle_full5x5_E5x5 = 0;
	leadingEle_full5x5_E1x5 = 0;
	leadingEle_full5x5_E2x5 = 0;
	leadingEle_full5x5_E2x5_Over_E5x5 = 0;
	leadingEle_full5x5_E1x5_Over_E5x5 = 0;
	leadingEle_EmHadDepth1Iso = 0;
	leadingEle_ptTracksIso = 0;
	leadingEle_innerLayerLostHits = 0;
	leadingEle_dxy = 0;
	leadingEle_eOverP = 0;
	leadingEle_id = 0;
	leadingEle_passCutBasedEleId = 0;
	subLeadingEle_passHEEPId = 0;
	subLeadingEle_etaSC = 0;
	subLeadingEle_isEcalDriven = 0;
	subLeadingEle_dEtaIn = 0;
	subLeadingEle_dPhiIn = 0;
	subLeadingEle_hOverE = 0;
	subLeadingEle_full5x5_r9 = 0;
	subLeadingEle_full5x5_sigmaIetaIeta = 0;
	subLeadingEle_full5x5_E5x5 = 0;
	subLeadingEle_full5x5_E1x5 = 0;
	subLeadingEle_full5x5_E2x5 = 0;
	subLeadingEle_full5x5_E2x5_Over_E5x5 = 0;
	subLeadingEle_full5x5_E1x5_Over_E5x5 = 0;
	subLeadingEle_EmHadDepth1Iso = 0;
	subLeadingEle_ptTracksIso = 0;
	subLeadingEle_innerLayerLostHits = 0;
	subLeadingEle_dxy = 0;
	subLeadingEle_eOverP = 0;
	subLeadingEle_id = 0;
	subLeadingEle_passCutBasedEleId = 0;
	leadingMuon_isHighPt = 0;
   	subLeadingMuon_isHighPt = 0;

	// Set branch addresses and branch pointers
	// if (!tree) return;
	// fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("run", &run, &b_run);
	fChain->SetBranchAddress("event", &event, &b_event);
	fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
	fChain->SetBranchAddress("rho", &rho, &b_rho);
	fChain->SetBranchAddress("weight", &weight, &b_weight);
	fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
	fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
	fChain->SetBranchAddress("npu", &npu, &b_npu);
	fChain->SetBranchAddress("passEEJJhlt", &passEEJJhlt, &b_passEEJJhlt);
	fChain->SetBranchAddress("passMMJJhlt", &passMMJJhlt, &b_passMMJJhlt);
	fChain->SetBranchAddress("passEMJJhlt", &passEMJJhlt, &b_passEMJJhlt);
	fChain->SetBranchAddress("passTandPEEhlt", &passTandPEEhlt, &b_passTandPEEhlt);
	fChain->SetBranchAddress("passTandPMMhlt", &passTandPMMhlt, &b_passTandPMMhlt);
	fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
	fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
	fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
	fChain->SetBranchAddress("ele_e", &ele_e, &b_ele_e);
	fChain->SetBranchAddress("ele_pt", &ele_pt, &b_ele_pt);
	fChain->SetBranchAddress("ele_eta", &ele_eta, &b_ele_eta);
	fChain->SetBranchAddress("ele_phi", &ele_phi, &b_ele_phi);
	fChain->SetBranchAddress("ele_idmva", &ele_idmva, &b_ele_idmva);
	fChain->SetBranchAddress("ele_id", &ele_id, &b_ele_id);
	fChain->SetBranchAddress("ele_iso", &ele_iso, &b_ele_iso);
	fChain->SetBranchAddress("ele_dz", &ele_dz, &b_ele_dz);
	fChain->SetBranchAddress("ele_d0", &ele_d0, &b_ele_d0);
	fChain->SetBranchAddress("ele_passHEEPId", &ele_passHEEPId, &b_ele_passHEEPId);
	fChain->SetBranchAddress("ele_passCutBasedEleId", &ele_passCutBasedEleId, &b_ele_passCutBasedEleId);
	fChain->SetBranchAddress("ele_isMatchedToGen", &ele_isMatchedToGen, &b_ele_isMatchedToGen);
	fChain->SetBranchAddress("ele_charge", &ele_charge, &b_ele_charge);
	fChain->SetBranchAddress("ele_etaSC", &ele_etaSC, &b_ele_etaSC);
	fChain->SetBranchAddress("ele_isEcalDriven", &ele_isEcalDriven, &b_ele_isEcalDriven);
	fChain->SetBranchAddress("ele_dEtaIn", &ele_dEtaIn, &b_ele_dEtaIn);
	fChain->SetBranchAddress("ele_dPhiIn", &ele_dPhiIn, &b_ele_dPhiIn);
	fChain->SetBranchAddress("ele_hOverE", &ele_hOverE, &b_ele_hOverE);
	fChain->SetBranchAddress("ele_full5x5_r9", &ele_full5x5_r9, &b_ele_full5x5_r9);
	fChain->SetBranchAddress("ele_full5x5_sigmaIetaIeta", &ele_full5x5_sigmaIetaIeta, &b_ele_full5x5_sigmaIetaIeta);
	fChain->SetBranchAddress("ele_full5x5_E5x5", &ele_full5x5_E5x5, &b_ele_full5x5_E5x5);
	fChain->SetBranchAddress("ele_full5x5_E1x5", &ele_full5x5_E1x5, &b_ele_full5x5_E1x5);
	fChain->SetBranchAddress("ele_full5x5_E2x5", &ele_full5x5_E2x5, &b_ele_full5x5_E2x5);
	fChain->SetBranchAddress("ele_full5x5_E2x5_Over_E5x5", &ele_full5x5_E2x5_Over_E5x5, &b_ele_full5x5_E2x5_Over_E5x5);
	fChain->SetBranchAddress("ele_full5x5_E1x5_Over_E5x5", &ele_full5x5_E1x5_Over_E5x5, &b_ele_full5x5_E1x5_Over_E5x5);
	fChain->SetBranchAddress("ele_EmHadDepth1Iso", &ele_EmHadDepth1Iso, &b_ele_EmHadDepth1Iso);
	fChain->SetBranchAddress("ele_ptTracksIso", &ele_ptTracksIso, &b_ele_ptTracksIso);
	fChain->SetBranchAddress("ele_innerLayerLostHits", &ele_innerLayerLostHits, &b_ele_innerLayerLostHits);
	fChain->SetBranchAddress("ele_dxy", &ele_dxy, &b_ele_dxy);
	fChain->SetBranchAddress("ele_eOverP", &ele_eOverP, &b_ele_eOverP);
	fChain->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
	fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
	fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
	fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
	fChain->SetBranchAddress("mu_iso", &mu_iso, &b_mu_iso);
	fChain->SetBranchAddress("mu_isTight", &mu_isTight, &b_mu_isTight);
	fChain->SetBranchAddress("mu_isMedium", &mu_isMedium, &b_mu_isMedium);
	fChain->SetBranchAddress("mu_isLoose", &mu_isLoose, &b_mu_isLoose);
	fChain->SetBranchAddress("mu_isHighPt", &mu_isHighPt, &b_mu_isHighPt);
	fChain->SetBranchAddress("mu_isMatchedToGen", &mu_isMatchedToGen, &b_mu_isMatchedToGen);
	fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
	fChain->SetBranchAddress("mu_dz", &mu_dz, &b_mu_dz);
	fChain->SetBranchAddress("mu_dxy", &mu_dxy, &b_mu_dxy);
	fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
	fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
	fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
	fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
	fChain->SetBranchAddress("jet_bdiscriminant", &jet_bdiscriminant, &b_jet_bdiscriminant);
	fChain->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
	fChain->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
	fChain->SetBranchAddress("jet_isMatchedToGen", &jet_isMatchedToGen, &b_jet_isMatchedToGen);
	fChain->SetBranchAddress("isEEJJ", &isEEJJ, &b_isEEJJ);
	fChain->SetBranchAddress("isEETT", &isEETT, &b_isEETT);
	fChain->SetBranchAddress("isMMJJ", &isMMJJ, &b_isMMJJ);
	fChain->SetBranchAddress("isMMTT", &isMMTT, &b_isMMTT);
	fChain->SetBranchAddress("isEMJJ", &isEMJJ, &b_isEMJJ);
	fChain->SetBranchAddress("isSignalRegion", &isSignalRegion, &b_isSignalRegion);
	fChain->SetBranchAddress("isLowMllCR", &isLowMllCR, &b_isLowMllCR);
	fChain->SetBranchAddress("isLowMlljjCR", &isLowMlljjCR, &b_isLowMlljjCR);
	fChain->SetBranchAddress("isBB", &isBB, &b_isBB);
	fChain->SetBranchAddress("isEE", &isEE, &b_isEE);
	fChain->SetBranchAddress("isEB", &isEB, &b_isEB);
	fChain->SetBranchAddress("passPreselections", &passPreselections, &b_passPreselections);
	fChain->SetBranchAddress("leadingLepton_e", &leadingLepton_e, &b_leadingLepton_e);
	fChain->SetBranchAddress("leadingLepton_pt", &leadingLepton_pt, &b_leadingLepton_pt);
	fChain->SetBranchAddress("leadingLepton_eta", &leadingLepton_eta, &b_leadingLepton_eta);
	fChain->SetBranchAddress("leadingLepton_phi", &leadingLepton_phi, &b_leadingLepton_phi);
	fChain->SetBranchAddress("leadingLepton_charge", &leadingLepton_charge, &b_leadingLepton_charge);
	fChain->SetBranchAddress("subLeadingLepton_e", &subLeadingLepton_e, &b_subLeadingLepton_e);
	fChain->SetBranchAddress("subLeadingLepton_pt", &subLeadingLepton_pt, &b_subLeadingLepton_pt);
	fChain->SetBranchAddress("subLeadingLepton_eta", &subLeadingLepton_eta, &b_subLeadingLepton_eta);
	fChain->SetBranchAddress("subLeadingLepton_phi", &subLeadingLepton_phi, &b_subLeadingLepton_phi);
	fChain->SetBranchAddress("subLeadingLepton_charge", &subLeadingLepton_charge, &b_subLeadingLepton_charge);
	fChain->SetBranchAddress("leadingJet_e", &leadingJet_e, &b_leadingJet_e);
	fChain->SetBranchAddress("leadingJet_pt", &leadingJet_pt, &b_leadingJet_pt);
	fChain->SetBranchAddress("leadingJet_eta", &leadingJet_eta, &b_leadingJet_eta);
	fChain->SetBranchAddress("leadingJet_phi", &leadingJet_phi, &b_leadingJet_phi);
	fChain->SetBranchAddress("subLeadingJet_e", &subLeadingJet_e, &b_subLeadingJet_e);
	fChain->SetBranchAddress("subLeadingJet_pt", &subLeadingJet_pt, &b_subLeadingJet_pt);
	fChain->SetBranchAddress("subLeadingJet_eta", &subLeadingJet_eta, &b_subLeadingJet_eta);
	fChain->SetBranchAddress("subLeadingJet_phi", &subLeadingJet_phi, &b_subLeadingJet_phi);
	fChain->SetBranchAddress("dRLeadLeptonLeadJet", &dRLeadLeptonLeadJet, &b_dRLeadLeptonLeadJet);
	fChain->SetBranchAddress("dRLeadLeptonSubLeadJet", &dRLeadLeptonSubLeadJet, &b_dRLeadLeptonSubLeadJet);
	fChain->SetBranchAddress("dRSubLeadLeptonLeadJet", &dRSubLeadLeptonLeadJet, &b_dRSubLeadLeptonLeadJet);
	fChain->SetBranchAddress("dRSubLeadLeptonSubLeadJet", &dRSubLeadLeptonSubLeadJet, &b_dRSubLeadLeptonSubLeadJet);
	fChain->SetBranchAddress("diLeptonDiJet_vtxIndex", &diLeptonDiJet_vtxIndex, &b_diLeptonDiJet_vtxIndex);
	fChain->SetBranchAddress("diLeptonDiJet_sumPt", &diLeptonDiJet_sumPt, &b_diLeptonDiJet_sumPt);
	fChain->SetBranchAddress("diLeptonDiJet_invMass", &diLeptonDiJet_invMass, &b_diLeptonDiJet_invMass);
	fChain->SetBranchAddress("diLepton_invMass", &diLepton_invMass, &b_diLepton_invMass);
	fChain->SetBranchAddress("diJet_invMass", &diJet_invMass, &b_diJet_invMass);
	fChain->SetBranchAddress("diJetLeadingLepton_invMass", &diJetLeadingLepton_invMass, &b_diJetLeadingLepton_invMass);
	fChain->SetBranchAddress("diJetSubLeadingLepton_invMass", &diJetSubLeadingLepton_invMass, &b_diJetSubLeadingLepton_invMass);
	fChain->SetBranchAddress("leadingEle_passHEEPId", &leadingEle_passHEEPId, &b_leadingEle_passHEEPId);
	fChain->SetBranchAddress("leadingEle_etaSC", &leadingEle_etaSC, &b_leadingEle_etaSC);
	fChain->SetBranchAddress("leadingEle_isEcalDriven", &leadingEle_isEcalDriven, &b_leadingEle_isEcalDriven);
	fChain->SetBranchAddress("leadingEle_dEtaIn", &leadingEle_dEtaIn, &b_leadingEle_dEtaIn);
	fChain->SetBranchAddress("leadingEle_dPhiIn", &leadingEle_dPhiIn, &b_leadingEle_dPhiIn);
	fChain->SetBranchAddress("leadingEle_hOverE", &leadingEle_hOverE, &b_leadingEle_hOverE);
	fChain->SetBranchAddress("leadingEle_full5x5_r9", &leadingEle_full5x5_r9, &b_leadingEle_full5x5_r9);
	fChain->SetBranchAddress("leadingEle_full5x5_sigmaIetaIeta", &leadingEle_full5x5_sigmaIetaIeta, &b_leadingEle_full5x5_sigmaIetaIeta);
	fChain->SetBranchAddress("leadingEle_full5x5_E5x5", &leadingEle_full5x5_E5x5, &b_leadingEle_full5x5_E5x5);
	fChain->SetBranchAddress("leadingEle_full5x5_E1x5", &leadingEle_full5x5_E1x5, &b_leadingEle_full5x5_E1x5);
	fChain->SetBranchAddress("leadingEle_full5x5_E2x5", &leadingEle_full5x5_E2x5, &b_leadingEle_full5x5_E2x5);
	fChain->SetBranchAddress("leadingEle_full5x5_E2x5_Over_E5x5", &leadingEle_full5x5_E2x5_Over_E5x5, &b_leadingEle_full5x5_E2x5_Over_E5x5);
	fChain->SetBranchAddress("leadingEle_full5x5_E1x5_Over_E5x5", &leadingEle_full5x5_E1x5_Over_E5x5, &b_leadingEle_full5x5_E1x5_Over_E5x5);
	fChain->SetBranchAddress("leadingEle_EmHadDepth1Iso", &leadingEle_EmHadDepth1Iso, &b_leadingEle_EmHadDepth1Iso);
	fChain->SetBranchAddress("leadingEle_ptTracksIso", &leadingEle_ptTracksIso, &b_leadingEle_ptTracksIso);
	fChain->SetBranchAddress("leadingEle_innerLayerLostHits", &leadingEle_innerLayerLostHits, &b_leadingEle_innerLayerLostHits);
	fChain->SetBranchAddress("leadingEle_dxy", &leadingEle_dxy, &b_leadingEle_dxy);
	fChain->SetBranchAddress("leadingEle_eOverP", &leadingEle_eOverP, &b_leadingEle_eOverP);
	fChain->SetBranchAddress("leadingEle_id", &leadingEle_id, &b_leadingEle_id);
	fChain->SetBranchAddress("leadingEle_passCutBasedEleId", &leadingEle_passCutBasedEleId, &b_leadingEle_passCutBasedEleId);
	fChain->SetBranchAddress("subLeadingEle_passHEEPId", &subLeadingEle_passHEEPId, &b_subLeadingEle_passHEEPId);
	fChain->SetBranchAddress("subLeadingEle_etaSC", &subLeadingEle_etaSC, &b_subLeadingEle_etaSC);
	fChain->SetBranchAddress("subLeadingEle_isEcalDriven", &subLeadingEle_isEcalDriven, &b_subLeadingEle_isEcalDriven);
	fChain->SetBranchAddress("subLeadingEle_dEtaIn", &subLeadingEle_dEtaIn, &b_subLeadingEle_dEtaIn);
	fChain->SetBranchAddress("subLeadingEle_dPhiIn", &subLeadingEle_dPhiIn, &b_subLeadingEle_dPhiIn);
	fChain->SetBranchAddress("subLeadingEle_hOverE", &subLeadingEle_hOverE, &b_subLeadingEle_hOverE);
	fChain->SetBranchAddress("subLeadingEle_full5x5_r9", &subLeadingEle_full5x5_r9, &b_subLeadingEle_full5x5_r9);
	fChain->SetBranchAddress("subLeadingEle_full5x5_sigmaIetaIeta", &subLeadingEle_full5x5_sigmaIetaIeta, &b_subLeadingEle_full5x5_sigmaIetaIeta);
	fChain->SetBranchAddress("subLeadingEle_full5x5_E5x5", &subLeadingEle_full5x5_E5x5, &b_subLeadingEle_full5x5_E5x5);
	fChain->SetBranchAddress("subLeadingEle_full5x5_E1x5", &subLeadingEle_full5x5_E1x5, &b_subLeadingEle_full5x5_E1x5);
	fChain->SetBranchAddress("subLeadingEle_full5x5_E2x5", &subLeadingEle_full5x5_E2x5, &b_subLeadingEle_full5x5_E2x5);
	fChain->SetBranchAddress("subLeadingEle_full5x5_E2x5_Over_E5x5", &subLeadingEle_full5x5_E2x5_Over_E5x5, &b_subLeadingEle_full5x5_E2x5_Over_E5x5);
	fChain->SetBranchAddress("subLeadingEle_full5x5_E1x5_Over_E5x5", &subLeadingEle_full5x5_E1x5_Over_E5x5, &b_subLeadingEle_full5x5_E1x5_Over_E5x5);
	fChain->SetBranchAddress("subLeadingEle_EmHadDepth1Iso", &subLeadingEle_EmHadDepth1Iso, &b_subLeadingEle_EmHadDepth1Iso);
	fChain->SetBranchAddress("subLeadingEle_ptTracksIso", &subLeadingEle_ptTracksIso, &b_subLeadingEle_ptTracksIso);
	fChain->SetBranchAddress("subLeadingEle_innerLayerLostHits", &subLeadingEle_innerLayerLostHits, &b_subLeadingEle_innerLayerLostHits);
	fChain->SetBranchAddress("subLeadingEle_dxy", &subLeadingEle_dxy, &b_subLeadingEle_dxy);
	fChain->SetBranchAddress("subLeadingEle_eOverP", &subLeadingEle_eOverP, &b_subLeadingEle_eOverP);
	fChain->SetBranchAddress("subLeadingEle_id", &subLeadingEle_id, &b_subLeadingEle_id);
	fChain->SetBranchAddress("subLeadingEle_passCutBasedEleId", &subLeadingEle_passCutBasedEleId, &b_subLeadingEle_passCutBasedEleId);
	fChain->SetBranchAddress("leadingMuon_isHighPt", &leadingMuon_isHighPt, &b_leadingMuon_isHighPt);
   	fChain->SetBranchAddress("subLeadingMuon_isHighPt", &subLeadingMuon_isHighPt, &b_subLeadingMuon_isHighPt);

	cout << "Branches are properly initialized." << endl;
}
// ******************************************************************************************


// **************** 
Int_t makePlots::GetEntry(Long64_t entry){
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
// ******************************************************************************************


// **************** 
Long64_t makePlots::LoadTree(Long64_t entry){
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
	}
	return centry;
}
// ******************************************************************************************



// **************** 
TH1D* makePlots::newTH1D(string name, string title, string xTitle, int nBins, double xLow, double xUp){
    TH1D* hist = new TH1D(name.c_str(), title.c_str(), nBins, xLow, xUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
	hist->GetYaxis()->SetTitleOffset(1.9);
    hist->SetOption("HIST");           
    listOfHistograms.push_back(hist);
    return hist;
}
// ******************************************************************************************


// **************** 
TH2D* makePlots::newTH2D(string name, string title, string xTitle, string yTitle, int nBinsX, double xLow, double xUp, int nBinsY, double yLow, double yUp){
    TH2D* hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xLow, xUp, nBinsY, yLow, yUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetYaxis()->SetTitleOffset(1.9); 
    listOfHistograms.push_back(hist);
    return hist;
}
// ******************************************************************************************


// **************** 
TProfile* makePlots::newTProfile(string name, string title, string xTitle, string yTitle, int nBinsX, double xLow, double xUp, double yLow, double yUp){
	TProfile* hist = new TProfile(name.c_str(), title.c_str(), nBinsX, xLow, xUp, yLow, yUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
	hist->GetYaxis()->SetTitleOffset(1.9);
    listOfHistograms.push_back(hist);
    return hist;
}
// ******************************************************************************************


// **************** 
void makePlots::SetHistos(){

	rho_histo = newTH1D("rho_histo", "rho_histo", "rho", 100, 0., 100.);
	nvtx_histo = newTH1D("nvtx_histo", "nvtx_histo", "nvtx", 100, 0., 100.);
	vtx_x_histo = newTH1D("vtx_x_histo", "vtx_x_histo", "vtx_x", 100, -0.3, 0.3);
	vtx_y_histo = newTH1D("vtx_y_histo", "vtx_y_histo", "vtx_y", 100, -0.3, 0.3);
	vtx_z_histo = newTH1D("vtx_z_histo", "vtx_z_histo", "vtx_z", 100, -20., 20.);

	muon_dxy_histo = newTH1D("muon_dxy_histo", "muon_dxy_histo", "|dxy| muons", 100, 0., 0.5);

	for (unsigned short i(0); i < 9; i++) {
		pt_histo[i] = newTH1D("pt_"+objName[i]+"_histo", "pt_"+objName[i]+"_histo", "p_{T} [GeV] "+objName[i], 100, 0., 1000.);
		eta_histo[i] = newTH1D("eta_"+objName[i]+"_histo", "eta_"+objName[i]+"_histo", "#eta "+objName[i], 100, -2.5, 2.5);
		phi_histo[i] = newTH1D("phi_"+objName[i]+"_histo", "phi_"+objName[i]+"_histo", "#phi "+objName[i], 100, -3.5, 3.5);
		// d0_histo[i] = newTH1D("d0_"+leptonName[i]+"_histo", "d0_"+leptonName[i]+"_histo", "|d0| "+leptonName[i], 100, 0., 0.5);
	}

	for (unsigned short j(0); j < 9; j++) {
		for (unsigned short m(0); m < 3; m++) {
			etaSC_histo[j][m] = newTH1D("etaSC_"+eleName[j]+"_histo"+etaName[m], "etaSC_"+eleName[j]+"_histo"+etaName[m], "#eta_{SC} "+eleName[j], 100, -2.5, 2.5);
			isEcalDriven_histo[j][m] = newTH1D("isEcalDriven_"+eleName[j]+"_histo"+etaName[m], "isEcalDriven_"+eleName[j]+"_histo"+etaName[m], "isEcalDriven "+eleName[j], 2, 0., 2.);
			dEtaIn_histo[j][m] = newTH1D("dEtaIn_"+eleName[j]+"_histo"+etaName[m], "dEtaIn_"+eleName[j]+"_histo"+etaName[m], "dEtaIn "+eleName[j], 100, -0.2, 0.2);
			dPhiIn_histo[j][m] = newTH1D("dPhiIn_"+eleName[j]+"_histo"+etaName[m], "dPhiIn_"+eleName[j]+"_histo"+etaName[m], "dPhiIn "+eleName[j], 100, -0.2, 0.2);
			hOverE_histo[j][m] = newTH1D("hOverE_"+eleName[j]+"_histo"+etaName[m], "hOverE_"+eleName[j]+"_histo"+etaName[m], "H/E "+eleName[j], 40, 0., 2.);
			full5x5_r9_histo[j][m] = newTH1D("r9_"+eleName[j]+"_histo"+etaName[m], "r9_"+eleName[j]+"_histo"+etaName[m], "r9 "+eleName[j], 100, 0., 1.05);
			full5x5_sigmaIetaIeta_histo[j][m] = newTH1D("sigmaIetaIeta_"+eleName[j]+"_histo"+etaName[m], "sigmaIetaIeta_"+eleName[j]+"_histo"+etaName[m], "sigmaIetaIeta "+eleName[j], 100, 0., 0.05);
			// full5x5_E5x5_ele_histo[j][m];
			// full5x5_E1x5_ele_histo[j][m];
			// full5x5_E2x5_ele_histo[j][m];
			full5x5_E2x5_Over_E5x5_histo[j][m] = newTH1D("e2x5_e5x5_"+eleName[j]+"_histo"+etaName[m], "e2x5_e5x5_"+eleName[j]+"_histo"+etaName[m], "e2x5/e5x5 "+eleName[j], 60, 0.5, 1.1);
			full5x5_E1x5_Over_E5x5_histo[j][m] = newTH1D("e1x5_e5x5_"+eleName[j]+"_histo"+etaName[m], "e1x5_e5x5_"+eleName[j]+"_histo"+etaName[m], "e1x5/e5x5 "+eleName[j], 100, 0., 1.1);
			EmHadDepth1Iso_histo[j][m] = newTH1D("EmHadDepth1Iso_"+eleName[j]+"_histo"+etaName[m], "EmHadDepth1Iso_"+eleName[j]+"_histo"+etaName[m], "EmHadDepth1Iso "+eleName[j], 200, 0., 200.);
			// ptTracksIso_histo[j][m];
			innerLayerLostHits_histo[j][m] = newTH1D("missingHits_"+eleName[j]+"_histo"+etaName[m], "missingHits_"+eleName[j]+"_histo"+etaName[m], "missingHits "+eleName[j], 10, 0., 10.);
			dxy_histo[j][m] = newTH1D("dxy_"+eleName[j]+"_histo"+etaName[m], "dxy_"+eleName[j]+"_histo"+etaName[m], "|dxy| "+eleName[j], 100, 0., 0.5);
			eOverP[j][m] = newTH1D("eOverP_"+eleName[j]+"_histo"+etaName[m], "eOverP_"+eleName[j]+"_histo"+etaName[m], "E/p "+eleName[j], 100, 0., 100);
		}
	}

	for (unsigned short l(0); l < 4; l++) {
		for (unsigned short n(0); n < 5; n++) {
			mass_dldj_histo[l][n] = newTH1D("mass_dileptondijets_histo_"+dldjName[l]+etaMassName[n], "mass_dileptondijets_histo_"+dldjName[l]+etaMassName[n], "diLeptonDiJet mass "+dldjName[l], 100, 0., 4000.);
			mass_dl_histo[l][n] = newTH1D("mass_dileptons_histo_"+dldjName[l]+etaMassName[n], "mass_dileptons_histo_"+dldjName[l]+etaMassName[n], "diLepton mass "+dldjName[l], 100, 0., 2000.);
			mass_dj_histo[l][n] = newTH1D("mass_dijets_histo_"+dldjName[l]+etaMassName[n], "mass_dijets_histo_"+dldjName[l]+etaMassName[n], "diJet mass "+dldjName[l], 100, 0., 2000.);		
			mass_djLl_histo[l][n] = newTH1D("mass_dijetsLeadingLepton_histo_"+dldjName[l]+etaMassName[n], "mass_dijetsLeadingLepton_histo_"+dldjName[l]+etaMassName[n], "diJetLeadingLepton mass "+dldjName[l], 100, 0., 4000.);	
			mass_djSLl_histo[l][n] = newTH1D("mass_dijetsSubLeadingLepton_histo_"+dldjName[l]+etaMassName[n], "mass_dijetsSubLeadingLepton_histo_"+dldjName[l]+etaMassName[n], "diJetSubLeadingLepton mass "+dldjName[l], 100, 0., 4000.);
			Zmass_histo[l][n] = newTH1D("Z"+Zname[l]+"_mass_histo"+etaMassName[n], "Z"+Zname[l]+"_mass_histo"+etaMassName[n], "m(Z"+Zname[l]+") [GeV/c^{2}]", 200, 0, 200); 
		}
	}
}
// ******************************************************************************************



// **************** 
void makePlots::getElectrons(vector<eleStruct>& electrons){
	unsigned short nTotElectrons(ele_e->size());
	for (unsigned short i(0); i < nTotElectrons; i++){
		eleStruct ele(ele_pt->at(i), 
						ele_eta->at(i), 						
						ele_phi->at(i), 
						ele_e->at(i), 
						ele_idmva->at(i), 
						ele_id->at(i), 
						ele_iso->at(i), 
						ele_dz->at(i),
						ele_d0->at(i),
						ele_passHEEPId->at(i), 
						ele_passCutBasedEleId->at(i), 
						ele_isMatchedToGen->at(i),						
						ele_charge->at(i),
						ele_etaSC->at(i),
						ele_isEcalDriven->at(i),
						ele_dEtaIn->at(i),
						ele_dPhiIn->at(i),
						ele_hOverE->at(i),
						ele_full5x5_r9->at(i),	
						ele_full5x5_sigmaIetaIeta->at(i),
						ele_full5x5_E5x5->at(i),
						ele_full5x5_E1x5->at(i),
						ele_full5x5_E2x5->at(i), 
						ele_full5x5_E2x5_Over_E5x5->at(i),
						ele_full5x5_E1x5_Over_E5x5->at(i),
						ele_EmHadDepth1Iso->at(i),
						0, //ele_ptTracksIso->at(i),
						ele_innerLayerLostHits->at(i), 
						ele_dxy->at(i),
						ele_eOverP->at(i) 
					); 
		electrons.push_back(ele);
	} //End of loop over all the electrons
}
// ******************************************************************************************


// **************** 
void makePlots::getLeadingElectrons(vector<eleStruct>& leadingElectrons){
	unsigned short nTotLeadingElectrons(leadingEle_etaSC->size());
	// cout << "nTotLeadingElectrons = " << nTotLeadingElectrons << endl; 

	for (unsigned short i(0); i < nTotLeadingElectrons; i++){
		eleStruct ele(leadingLepton_pt->at(i),
						leadingLepton_eta->at(i),
						leadingLepton_phi->at(i),
						leadingLepton_e->at(i), 
						0,
						leadingEle_id->at(i),	
						0,
						0,
						0,
						leadingEle_passHEEPId->at(i),
						leadingEle_passCutBasedEleId->at(i),	
						0,
						leadingLepton_charge->at(i), 
						leadingEle_etaSC->at(i),
						leadingEle_isEcalDriven->at(i),
						leadingEle_dEtaIn->at(i),
						leadingEle_dPhiIn->at(i),
						leadingEle_hOverE->at(i),
						leadingEle_full5x5_r9->at(i),
						leadingEle_full5x5_sigmaIetaIeta->at(i),
						leadingEle_full5x5_E5x5->at(i),
						leadingEle_full5x5_E1x5->at(i),
						leadingEle_full5x5_E2x5->at(i),		
						leadingEle_full5x5_E2x5_Over_E5x5->at(i),
						leadingEle_full5x5_E1x5_Over_E5x5->at(i),
						leadingEle_EmHadDepth1Iso->at(i),
						0, //leadingEle_ptTracksIso->at(i),
						leadingEle_innerLayerLostHits->at(i),
						leadingEle_dxy->at(i),
						leadingEle_eOverP->at(i)						
					);
		leadingElectrons.push_back(ele);
	} //End of loop over all the leading electrons
}
// ******************************************************************************************


// **************** 
void makePlots::getSubLeadingElectrons(vector<eleStruct>& subLeadingElectrons){
	unsigned short nTotSubLeadingElectrons(subLeadingEle_etaSC->size());
	// cout << "nTotSubLeadingElectrons = " << nTotSubLeadingElectrons << endl; 

	for (unsigned short i(0); i < nTotSubLeadingElectrons; i++){
		eleStruct ele(subLeadingLepton_pt->at(i),
						subLeadingLepton_eta->at(i),
						subLeadingLepton_phi->at(i),
						subLeadingLepton_e->at(i), 
						0,
						subLeadingEle_id->at(i),	
						0,
						0,
						0,
						subLeadingEle_passHEEPId->at(i),
						subLeadingEle_passCutBasedEleId->at(i),	
						0,
						subLeadingLepton_charge->at(i), 
						subLeadingEle_etaSC->at(i),
						subLeadingEle_isEcalDriven->at(i),
						subLeadingEle_dEtaIn->at(i),
						subLeadingEle_dPhiIn->at(i),
						subLeadingEle_hOverE->at(i),
						subLeadingEle_full5x5_r9->at(i),
						subLeadingEle_full5x5_sigmaIetaIeta->at(i),
						subLeadingEle_full5x5_E5x5->at(i),
						subLeadingEle_full5x5_E1x5->at(i),
						subLeadingEle_full5x5_E2x5->at(i),		
						subLeadingEle_full5x5_E2x5_Over_E5x5->at(i),
						subLeadingEle_full5x5_E1x5_Over_E5x5->at(i),
						subLeadingEle_EmHadDepth1Iso->at(i),
						0, //subLeadingEle_ptTracksIso->at(i),
						subLeadingEle_innerLayerLostHits->at(i),
						subLeadingEle_dxy->at(i),
						subLeadingEle_eOverP->at(i)						
					);
		subLeadingElectrons.push_back(ele);
	} //End of loop over all the subLeading electrons
}
// ******************************************************************************************



// **************** 
void makePlots::getMuons(vector<muonStruct>& muons){
	//--- get the number of Muon candidates from the vector size ---
	unsigned short nTotMuons(mu_e->size());
	 
	for (unsigned short i(0); i < nTotMuons; i++) {
		muonStruct mu(mu_pt->at(i), 
				mu_eta->at(i), 
				mu_phi->at(i), 
				mu_e->at(i), 
				mu_iso->at(i),	
				mu_isHighPt->at(i),
				mu_isMatchedToGen->at(i),
				mu_charge->at(i),
				mu_dz->at(i),
				mu_dxy->at(i)				
			);					
		muons.push_back(mu);
	}//End of loop over all the muons
}
// ******************************************************************************************


// **************** 
void makePlots::getLeadingMuons(vector<muonStruct>& leadingMuons){
	//--- get the number of Muon candidates from the vector size ---
	unsigned short nTotLeadingMuons(leadingMuon_isHighPt->size());
	 
	for (unsigned short i(0); i < nTotLeadingMuons; i++) {
		muonStruct mu(leadingLepton_pt->at(i),
					leadingLepton_eta->at(i),
					leadingLepton_phi->at(i),
					leadingLepton_e->at(i),
					0,
					leadingMuon_isHighPt->at(i),
					0,
					leadingLepton_charge->at(i), 	
					0,
					0);					
		leadingMuons.push_back(mu);
	}//End of loop over all the leading muons
}
// ******************************************************************************************


// **************** 
void makePlots::getSubLeadingMuons(vector<muonStruct>& subLeadingMuons){
	//--- get the number of Muon candidates from the vector size ---
	unsigned short nTotSubLeadingMuons(subLeadingMuon_isHighPt->size());
	 
	for (unsigned short i(0); i < nTotSubLeadingMuons; i++) {
		muonStruct mu(subLeadingLepton_pt->at(i),
					subLeadingLepton_eta->at(i),
					subLeadingLepton_phi->at(i),
					subLeadingLepton_e->at(i),
					0,
					subLeadingMuon_isHighPt->at(i),
					0,
					subLeadingLepton_charge->at(i), 	
					0,
					0);					
		subLeadingMuons.push_back(mu);
	}//End of loop over all the subleading muons
}
// ******************************************************************************************



// **************** 
void makePlots::getJets(vector<jetStruct>& jets){
	unsigned short nTotJets(jet_e->size());

	for (unsigned short i(0); i < nTotJets; i++){
		jetStruct jet(jet_pt->at(i),
				jet_eta->at(i),
				jet_phi->at(i),
				jet_e->at(i),
				jet_isMatchedToGen->at(i)
				);
		jets.push_back(jet);
	} //End of loop over all the jets
}
// ******************************************************************************************


// **************** 
void makePlots::getLeadingJets(vector<jetStruct>& leadingJets){
	unsigned short nTotLeadingJets(leadingJet_pt->size());

	for (unsigned short i(0); i < nTotLeadingJets; i++){
		jetStruct jet(leadingJet_pt->at(i),
				leadingJet_eta->at(i),
				leadingJet_phi->at(i),
				leadingJet_e->at(i),
				0
				);
		leadingJets.push_back(jet);
	} //End of loop over all the leading jets
}
// ******************************************************************************************


// **************** 
void makePlots::getSubLeadingJets(vector<jetStruct>& subLeadingJets){
	unsigned short nTotSubLeadingJets(subLeadingJet_pt->size());

	for (unsigned short i(0); i < nTotSubLeadingJets; i++){
		jetStruct jet(subLeadingJet_pt->at(i),
				subLeadingJet_eta->at(i),
				subLeadingJet_phi->at(i),
				subLeadingJet_e->at(i),
				0
				);
		subLeadingJets.push_back(jet);
	} //End of loop over all the leading jets
}
// ******************************************************************************************



// **************** 
void makePlots::getDiLeptonDiJets(vector<diLeptonDiJetStruct>& diLeptonDiJets){
	unsigned short nTotDiLeptonDiJets(diLeptonDiJet_invMass->size());
	// cout << "nTotDiLeptonDiJets = " << nTotDiLeptonDiJets << endl;

	for (unsigned short i(0); i < nTotDiLeptonDiJets; i++){
		diLeptonDiJetStruct dldj(isEEJJ->at(i),
							isEETT->at(i),
							isMMJJ->at(i),
							isMMTT->at(i),
							isEMJJ->at(i),
							isSignalRegion->at(i),
							isLowMllCR->at(i),
							isLowMlljjCR->at(i),
							isBB->at(i),
							isEE->at(i),
							isEB->at(i),
							passPreselections->at(i),
							dRLeadLeptonLeadJet->at(i),
							dRLeadLeptonSubLeadJet->at(i),
							dRSubLeadLeptonLeadJet->at(i),
							dRSubLeadLeptonSubLeadJet->at(i),
							diLeptonDiJet_vtxIndex->at(i),
							diLeptonDiJet_sumPt->at(i),
							diLeptonDiJet_invMass->at(i),
							diLepton_invMass->at(i),
							diJet_invMass->at(i),
							diJetLeadingLepton_invMass->at(i),
							diJetSubLeadingLepton_invMass->at(i)
						);
		diLeptonDiJets.push_back(dldj);
	} //End of loop over all the diLeptonDiJets
}
// ******************************************************************************************



// **************** 
void makePlots::doEleDistributionsPlots(vector<eleStruct>& electrons, const int eleIdx, const int histoIdx) { 
	pt_histo[histoIdx]->Fill(electrons[eleIdx].v.Pt(),w);
	eta_histo[histoIdx]->Fill(electrons[eleIdx].v.Eta(),w);
	phi_histo[histoIdx]->Fill(electrons[eleIdx].v.Phi(),w); 
}
// ******************************************************************************************


// **************** 
void makePlots::doMuonDistributionsPlots(vector<muonStruct>& muons, const int muIdx, const int histoIdx) { 
	pt_histo[histoIdx]->Fill(muons[muIdx].v.Pt(),w);
	eta_histo[histoIdx]->Fill(muons[muIdx].v.Eta(),w);
	phi_histo[histoIdx]->Fill(muons[muIdx].v.Phi(),w); 
}
// ******************************************************************************************


// **************** 
void makePlots::doJetsDistributionsPlots(vector<jetStruct>& jets, const int jetIdx, const int histoIdx) { 
	pt_histo[histoIdx]->Fill(jets[jetIdx].v.Pt(),w);
	eta_histo[histoIdx]->Fill(jets[jetIdx].v.Eta(),w);
	phi_histo[histoIdx]->Fill(jets[jetIdx].v.Phi(),w); 
}
// ******************************************************************************************



// **************** 
void makePlots::doElePlots(vector<eleStruct>& electrons, const int eleIdx, const int histoIdx, const int etaIdx){
	etaSC_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].etaSC,w);
	isEcalDriven_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].isEcalDriven,w);
	dEtaIn_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].dEtaIn,w);
	dPhiIn_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].dPhiIn,w);
	hOverE_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].HoE,w);
	full5x5_r9_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].r9,w);
	full5x5_sigmaIetaIeta_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].sigmaIetaIeta,w);
	full5x5_E2x5_Over_E5x5_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].e2x5_e5x5,w);
	full5x5_E1x5_Over_E5x5_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].e1x5_e5x5,w);
	EmHadDepth1Iso_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].EmHadDepth1Iso,w);
	innerLayerLostHits_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].missingHits,w);
	dxy_histo[histoIdx][etaIdx]->Fill(electrons[eleIdx].dxy,w);
	eOverP[histoIdx][etaIdx]->Fill(electrons[eleIdx].eOverP,w);
}
// ******************************************************************************************



// **************** 
void makePlots::doMassPlots(vector<diLeptonDiJetStruct>& diLeptonDiJets, const int dldjIdx, const int histoIdx, const int etaIdx){
	mass_dldj_histo[histoIdx][etaIdx]->Fill(diLeptonDiJets[dldjIdx].diLeptonDiJet_invMass,w);
	mass_dl_histo[histoIdx][etaIdx]->Fill(diLeptonDiJets[dldjIdx].diLepton_invMass,w);
	mass_dj_histo[histoIdx][etaIdx]->Fill(diLeptonDiJets[dldjIdx].diJet_invMass,w);
	mass_djLl_histo[histoIdx][etaIdx]->Fill(diLeptonDiJets[dldjIdx].diJetLeadingLepton_invMass,w);
	mass_djSLl_histo[histoIdx][etaIdx]->Fill(diLeptonDiJets[dldjIdx].diJetSubLeadingLepton_invMass,w);
}
// ******************************************************************************************


// **************** 
void makePlots::doZmassPlots(vector<diLeptonDiJetStruct>& diLeptonDiJets, const int dldjIdx, const int histoIdx, const int etaIdx) { 	
	float mZ = diLeptonDiJets[dldjIdx].diLepton_invMass;
	Zmass_histo[histoIdx][etaIdx]->Fill(mZ,w);
}
// ******************************************************************************************






// **************** 
void makePlots::Loop(){
	//--- Initialize the tree branches ---
	Init();
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	cout << "nentries = " << nentries << endl;


	// for (Long64_t i=0; i<2153862; i+=1) { //for WJets
	// for (Long64_t i=0; i<10; i+=1) {
	for (Long64_t i=0; i<nentries; i+=1) {
		Long64_t ientry = LoadTree(i);
		if (ientry < 0) break;
		// cout << "evento " << ientry << endl;
		nEvents++;

		if(fChain->GetEntry(i) == 0){
			std::cerr << "Failed to read Tree entry " << i << "!\n";
			continue;
		}

		//Trigger
		bool passesTrigger = 1; //sempre vero per MC	
		if (!MC) {
			// cout << "data" << endl; 
			if (signalEE) passesTrigger = passEEJJhlt;   
			if (signalMuMu) passesTrigger = passMMJJhlt;
			if (eMuSideband) passesTrigger = passEMJJhlt;
			if (TnPee) passesTrigger = passTandPEEhlt;
			if (TnPmumu) passesTrigger = passTandPMMhlt;
		}

		if (!passesTrigger) continue;  //se evento non passa il trigger passo a quello successivo
		nEventsPassingTrigger++;

		w = weight;
		if (MC) w *= lumiData;
		if (MCpuReweighted) w *= puweight;
		// cout << "weight = " << w << endl;


		rho_histo->Fill(rho,w);
		nvtx_histo->Fill(nvtx,w);
		// cout << "filla nvtx" << endl;


		for (unsigned short i(0); i < nvtx; i++){
			vtx_x_histo->Fill(	vtx_x->at(i) );
			vtx_y_histo->Fill(	vtx_y->at(i) );
			vtx_z_histo->Fill(	vtx_z->at(i) );
		}


		//Collections
		vector<eleStruct> electrons, leadingElectrons, subLeadingElectrons;
		vector<diLeptonDiJetStruct> diLeptonDiJets;
		vector<muonStruct> muons, leadingMuons, subLeadingMuons;
		vector<jetStruct> jets, leadingJets, subLeadingJets;


		// cout << "get electrons" << endl;
		getElectrons(electrons);

		// cout << "get muons" << endl;
		getMuons(muons);

		// cout << "get jets" << endl;
		getJets(jets);

		// cout << "get dldj" << endl;
		getDiLeptonDiJets(diLeptonDiJets);

		// cout << "get leadingEle" << endl;
		getLeadingElectrons(leadingElectrons);

		// cout << "get subLeadingEle" << endl;
		getSubLeadingElectrons(subLeadingElectrons);

		// cout << "get leadingMuon" << endl;
		getLeadingMuons(leadingMuons);
		// cout << "get subLeadingMuon" << endl;
		getSubLeadingMuons(subLeadingMuons);

		// cout << "get leadingJet" << endl;
		getLeadingJets(leadingJets);
		// cout << "get subLeadingJet" << endl;
		getSubLeadingJets(subLeadingJets);


		unsigned short nElectrons = electrons.size();    
		// cout << "nElectrons = " << nElectrons << endl;    

		unsigned short nMuons = muons.size();
		unsigned short nJets = jets.size();

		unsigned short nDiLeptonDiJets = diLeptonDiJets.size();
		// cout << "nDiLeptonDiJets = " << nDiLeptonDiJets << endl;


	// --electrons
		for (unsigned short j(0); j < nElectrons; j++){

			doEleDistributionsPlots(electrons,j,0);
			doElePlots(electrons,j,0,0);

			bool inEB = isInEB(electrons[j].v.Eta());
			bool inEE = isInEE(electrons[j].v.Eta());

			if (inEB) doElePlots(electrons,j,0,1);
			if (inEE) doElePlots(electrons,j,0,2);

			if (electrons[j].v.Pt() < 35) continue;   

			if (electrons[j].passHEEPId) {
				doElePlots(electrons,j,3,0);
				if (inEB) doElePlots(electrons,j,3,1);
				if (inEE) doElePlots(electrons,j,3,2);
			}

			if (electrons[j].passCutBasedEleId) {
				doElePlots(electrons,j,6,0);
				if (inEB) doElePlots(electrons,j,6,1);
				if (inEE) doElePlots(electrons,j,6,2);
			}
		}


	// --muons
		// unsigned short nMuons(mu_e->size());
		for (unsigned short i(0); i < nMuons; i++){
			doMuonDistributionsPlots(muons,i,3);
			muon_dxy_histo->Fill( mu_dxy->at(i) );
		}


	// --jets
		for (unsigned short l(0); l < nJets; l++){
			doJetsDistributionsPlots(jets,l,6);
		}


	// --leadingEle
		for (unsigned short j(0); j < nDiLeptonDiJets; j++){
			// cout << "leadingEle " << j << endl;
			if ( diLeptonDiJets[j].isEEJJ || (diLeptonDiJets[j].isEMJJ && leadingElectrons[j].etaSC != -999) ) {  

				doEleDistributionsPlots(leadingElectrons,j,1);

				bool lInEB = isInEB(leadingElectrons[j].v.Eta());  
				bool lInEE = isInEE(leadingElectrons[j].v.Eta());

				doElePlots(leadingElectrons,j,1,0);
				if (lInEB) doElePlots(leadingElectrons,j,1,1);
				if (lInEE) doElePlots(leadingElectrons,j,1,2);

				if (leadingElectrons[j].v.Pt() < 35) continue;    

				if (leadingElectrons[j].passHEEPId) {
					doElePlots(leadingElectrons,j,4,0);
					if (lInEB) doElePlots(leadingElectrons,j,4,1);
					if (lInEE) doElePlots(leadingElectrons,j,4,2);
				}

				if (leadingElectrons[j].passCutBasedEleId) {
					doElePlots(leadingElectrons,j,7,0);
					if (lInEB) doElePlots(leadingElectrons,j,7,1);
					if (lInEE) doElePlots(leadingElectrons,j,7,2);
				}
			}
		}


	// --subLeadingEle
		for (unsigned short j(0); j < nDiLeptonDiJets; j++){
			// cout << "subLeadingEle " << j << endl;
			if ( diLeptonDiJets[j].isEEJJ || (diLeptonDiJets[j].isEMJJ && subLeadingElectrons[j].etaSC != -999) ) {  

				doEleDistributionsPlots(subLeadingElectrons,j,2);

				bool slInEB = isInEB(subLeadingElectrons[j].v.Eta());  
				bool slInEE = isInEE(subLeadingElectrons[j].v.Eta());

				doElePlots(subLeadingElectrons,j,2,0);
				if (slInEB) doElePlots(subLeadingElectrons,j,2,1);
				if (slInEE) doElePlots(subLeadingElectrons,j,2,2);

				if (subLeadingElectrons[j].v.Pt() < 35) continue;     
 
				if (subLeadingElectrons[j].passHEEPId) {
					doElePlots(subLeadingElectrons,j,5,0);
					if (slInEB) doElePlots(subLeadingElectrons,j,5,1);
					if (slInEE) doElePlots(subLeadingElectrons,j,5,2);
				}

				if (subLeadingElectrons[j].passCutBasedEleId) {
					doElePlots(subLeadingElectrons,j,8,0);
					if (slInEB) doElePlots(subLeadingElectrons,j,8,1);
					if (slInEE) doElePlots(subLeadingElectrons,j,8,2);
				}
			}
		}


	// --leadingMuon
		for (unsigned short j(0); j < nDiLeptonDiJets; j++){
			// cout << "leadingMuon " << j << endl;
			if ( diLeptonDiJets[j].isMMJJ || (diLeptonDiJets[j].isEMJJ && leadingMuons[j].isHighPt = -999) ) {  
				doMuonDistributionsPlots(leadingMuons,j,4);
			}
		}


	// --subLeadingMuon
		for (unsigned short j(0); j < nDiLeptonDiJets; j++){
			// cout << "subLeadingMuon " << j << endl;
			if ( diLeptonDiJets[j].isMMJJ || (diLeptonDiJets[j].isEMJJ && subLeadingMuons[j].isHighPt = -999) ) { 
				doMuonDistributionsPlots(subLeadingMuons,j,5);
			}
		}



	// --dldj
		int nEEJJ=0, nMMJJ=0;
		// bool signalRegion=0, flavourSidebandCR=0, lowMllCR=0, lowMlljjCR=0, TnP_CR=0;

		for (unsigned short j(0); j < nDiLeptonDiJets; j++){
			// cout << "dldj " << j << endl;

			// //-- definition of signal/contro regions
			// if (diLeptonDiJets[j].passPreselections) {
			// 	if (diLeptonDiJets[j].isEEJJ || diLeptonDiJets[j].isMMJJ) {
			// 		if (diLeptonDiJets[j].isSignalRegion) signalRegion=true;
			// 		if (diLeptonDiJets[j].isLowMllCR) lowMllCR=true;
			// 		if (diLeptonDiJets[j].isLowMlljjCR) lowMlljjCR=true;
			// 	}
			// 	if (diLeptonDiJets[j].isEMJJ && diLeptonDiJets[j].isSignalRegion) flavourSidebandCR=true;
			// }
			// if (TnPee || TnPmumu) TnP_CR=true; //eventi hanno passato il trigger (vedi sopra) 
			// //--
					
			bool isInEBEB = inEBEB(leadingLepton_eta->at(j), subLeadingLepton_eta->at(j));  //usare eta jets per djets mass ?
			bool isInEEEE = inEEEE(leadingLepton_eta->at(j), subLeadingLepton_eta->at(j));
			bool isInEBEE = inEBEE(leadingLepton_eta->at(j), subLeadingLepton_eta->at(j)); 

			if (diLeptonDiJets[j].isEEJJ || diLeptonDiJets[j].isMMJJ || diLeptonDiJets[j].isEMJJ) { 
				doJetsDistributionsPlots(leadingJets,j,7);
				doJetsDistributionsPlots(subLeadingJets,j,8);
			}

			doMassPlots(diLeptonDiJets,j,0,0);

			if (isInEBEB) doMassPlots(diLeptonDiJets,j,0,1);
			if (isInEEEE) doMassPlots(diLeptonDiJets,j,0,2);
			if (isInEBEE) doMassPlots(diLeptonDiJets,j,0,3);
			if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,j,0,4);

			if (diLeptonDiJets[j].passPreselections) {
				doMassPlots(diLeptonDiJets,j,1,0);
				if (isInEBEB) doMassPlots(diLeptonDiJets,j,1,1);
				if (isInEEEE) doMassPlots(diLeptonDiJets,j,1,2);
				if (isInEBEE) doMassPlots(diLeptonDiJets,j,1,3);
				if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,j,1,4);

				if (diLeptonDiJets[j].isEEJJ && leadingElectrons[j].passHEEPId && subLeadingElectrons[j].passHEEPId) {
					doMassPlots(diLeptonDiJets,j,2,0);
					if (isInEBEB) doMassPlots(diLeptonDiJets,j,2,1);
					if (isInEEEE) doMassPlots(diLeptonDiJets,j,2,2);
					if (isInEBEE) doMassPlots(diLeptonDiJets,j,2,3);
					if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,j,2,4);
				}

				if (diLeptonDiJets[j].isEEJJ && leadingElectrons[j].passCutBasedEleId && subLeadingElectrons[j].passCutBasedEleId) {
					doMassPlots(diLeptonDiJets,j,3,0);
					if (isInEBEB) doMassPlots(diLeptonDiJets,j,3,1);
					if (isInEEEE) doMassPlots(diLeptonDiJets,j,3,2);
					if (isInEBEE) doMassPlots(diLeptonDiJets,j,3,3);
					if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,j,3,4);
				}
			}


			if (diLeptonDiJets[j].isEEJJ) {  //ho 2 ele
				nEEJJ++;
				if (leadingElectrons[j].v.Pt() < 35 || subLeadingElectrons[j].v.Pt() < 35) continue;

				//Z->ee
				if (leadingElectrons[j].charge * subLeadingElectrons[j].charge < 0 ) {	// & TnP_CR		
					doZmassPlots(diLeptonDiJets,j,0,0);
					if (isInEBEB) doZmassPlots(diLeptonDiJets,j,0,1);
					if (isInEEEE) doZmassPlots(diLeptonDiJets,j,0,2);
					if (isInEBEE) doZmassPlots(diLeptonDiJets,j,0,3);
					if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,j,0,4);

					if (leadingElectrons[j].passHEEPId && subLeadingElectrons[j].passHEEPId) {
						doZmassPlots(diLeptonDiJets,j,2,0);
						if (isInEBEB) doZmassPlots(diLeptonDiJets,j,2,1);
						if (isInEEEE) doZmassPlots(diLeptonDiJets,j,2,2);
						if (isInEBEE) doZmassPlots(diLeptonDiJets,j,2,3);
						if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,j,2,4);
					}

					if (leadingElectrons[j].passCutBasedEleId && subLeadingElectrons[j].passCutBasedEleId) {
						doZmassPlots(diLeptonDiJets,j,3,0);
						if (isInEBEB) doZmassPlots(diLeptonDiJets,j,3,1);
						if (isInEEEE) doZmassPlots(diLeptonDiJets,j,3,2);
						if (isInEBEE) doZmassPlots(diLeptonDiJets,j,3,3);
						if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,j,3,4);
					}
				}
			}


			if (diLeptonDiJets[j].isMMJJ) {  //ho 2 muons
				nMMJJ++;
				if (leadingMuons[j].v.Pt() < 35 || subLeadingMuons[j].v.Pt() < 35) continue; 

				//Z->mumu
				if (leadingMuons[j].charge * subLeadingMuons[j].charge < 0 ) {  // & TnP_CR	
					doZmassPlots(diLeptonDiJets,j,1,0);
					if (isInEBEB) doZmassPlots(diLeptonDiJets,j,1,1);
					if (isInEEEE) doZmassPlots(diLeptonDiJets,j,1,2);
					if (isInEBEE) doZmassPlots(diLeptonDiJets,j,1,3);
					if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,j,1,4);

				}
			}

		}
		// cout << "nEEJJ = " << nEEJJ << endl;

	} //fine loop su eventi
	
}
// ******************************************************************************************




// **************** 
void makePlots::saveHistosAndOutputFile(TString& ouputdir){
	// setTDRStyle();
	// writeExtraText = true; 
	// extraText = "Preliminary"; 
	// lumi_13TeV = "0.6 fb^{-1}"; 

	unsigned short numbOfHistograms = listOfHistograms.size();

	TFile f(outputdir+suff+"/"+"distributions_histos"+suff+".root","recreate");
	for (unsigned short l(0); l < numbOfHistograms; l++){
		string histoName = listOfHistograms[l]->GetName();
		listOfHistograms[l]->Write(); 
		// if ( (histoName.find("nvtx") != string::npos) || (histoName.find("pt_") != string::npos) || 
		// 	(histoName.find("eta_") != string::npos) || (histoName.find("phi_") != string::npos) ||
		// 	(histoName.find("etaSC_") != string::npos) || (histoName.find("isEcalDriven_") != string::npos) ||
		// 	(histoName.find("dEtaIn_") != string::npos) || (histoName.find("dPhiIn_") != string::npos) ||
		// 	(histoName.find("hOverE_") != string::npos) || (histoName.find("sigmaIetaIeta_") != string::npos) ||
		// 	(histoName.find("e2x5_e5x5_") != string::npos) || (histoName.find("e1x5_e5x5_") != string::npos) ||
		// 	(histoName.find("EmHadDepth1Iso_") != string::npos) || (histoName.find("missingHits_") != string::npos) ||
		// 	(histoName.find("dxy_") != string::npos) || (histoName.find("mass_") != string::npos) ||
		// 	) {
		TCanvas* canv = new TCanvas();
		gPad->SetRightMargin(0.05);
		// CMS_lumi(canv, 4, 33);
		listOfHistograms[l]->Draw();
		canv->SaveAs(outputdir+suff+"/"+histoName+".png");
		delete canv;
		// }
	}

}
// ******************************************************************************************



makePlots::makePlots(TString filename_, TString outputdir_, bool MC_, bool MCpuReweighted_, bool signalEE_, bool signalMuMu_, bool eMuSideband_, bool TnPee_, bool TnPmumu_):
	filename(filename_), outputdir(outputdir_), MC(MC_), MCpuReweighted(MCpuReweighted_), signalEE(signalEE_), signalMuMu(signalMuMu_), eMuSideband(eMuSideband_), TnPee(TnPee_), TnPmumu(TnPmumu_)
{
	fChain = new TChain("", "");

	TFile *f = TFile::Open(filename);
	if(!f){
	    cerr << "Failed to open file " << filename << ".\n";
	} else {
	    cout << "Reading input file " << filename << "\n";
		TString treePath = filename + "/analysisTree/event";
		if (fChain) fChain->Add(treePath);
	}

	lumiData = 5.899; //RunB

	SetHistos();
	Loop();
	saveHistosAndOutputFile(outputdir);
}



void runMakePlots() {

	string inputDir = "root://node12.datagrid.cea.fr//dpm/datagrid.cea.fr/home/cms/trivcat/store/user/gnegro/miniTrees/";
	string outputDir = "distributions_lumiReweighted"; //distributions

	// makePlots("root://node12.datagrid.cea.fr//dpm/datagrid.cea.fr/home/cms/trivcat/store/user/gnegro/miniTrees/output_DoubleEG_Run2016B-ReReco-v3_miniTrees.root", "Analysis/miniTrees/DoubleEG_Run2016B-ReReco-v3/distributions", false, false, true, false, false, false, false);

	// makePlots("root://node12.datagrid.cea.fr//dpm/datagrid.cea.fr/home/cms/trivcat/store/user/gnegro/miniTrees/output_DYJetsToLL-amcatnloFXFX_miniTrees.root", "Analysis/miniTrees/DYJetsToLL-amcatnlo/distributions", true, true, false, false, false, false, false);

	// makePlots("root://node12.datagrid.cea.fr//dpm/datagrid.cea.fr/home/cms/trivcat/store/user/gnegro/miniTrees/output_WRToEEJJ_1600_miniTrees.root", "Analysis/miniTrees/WR-1600_ToLNu-800_ToEEJJ/distributions", true, false, true, false, false, false, false); //no ripesato per fare contronto con official MC
	// makePlots("root://node12.datagrid.cea.fr//dpm/datagrid.cea.fr/home/cms/trivcat/store/user/gnegro/miniTrees/output_WRToMuMuJJ_1600_miniTrees.root", "Analysis/miniTrees/WR-1600_ToLNu-800_ToMuMuJJ/distributions", true, true, false, true, false, false, false);

	// makePlots("root://node12.datagrid.cea.fr//dpm/datagrid.cea.fr/home/cms/trivcat/store/user/gnegro/miniTrees/output_officialWRToEEJJ_1600_miniTrees.root", "Analysis/miniTrees/WR-1600_ToLNu-800_ToEEJJ_official/distributions", true, false, true, false, false, false, false); //no ripesato per fare contronto con miei MC


	// makePlots(inputDir+"output_ZToEE_miniTrees.root", "Analysis/miniTrees/ZToEE/"+outputDir, true, true, false, false, false, false, false);

	// makePlots(inputDir+"output_ZToMuMu_miniTrees.root", "Analysis/miniTrees/ZToMuMu/"+outputDir, true, true, false, false, false, false, false);

	makePlots(inputDir+"output_ZZ_miniTrees.root", "Analysis/miniTrees/ZZ/"+outputDir, true, true, false, false, false, false, false);

	// makePlots(inputDir+"output_WZ_miniTrees.root", "Analysis/miniTrees/WZ/"+outputDir, true, true, false, false, false, false, false);

	// makePlots(inputDir+"output_WJetsToLNu_miniTrees.root", "Analysis/miniTrees/WJetsToLNu/"+outputDir, true, true, false, false, false, false, false);

	// makePlots(inputDir+"output_TTJets-v4_miniTrees.root", "Analysis/miniTrees/TTJets/"+outputDir, true, true, false, false, false, false, false);

}
