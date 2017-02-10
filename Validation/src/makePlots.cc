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
	diLepton_pt = 0;
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
	fChain->SetBranchAddress("diLepton_pt", &diLepton_pt, &b_diLepton_pt);
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

	nLeadingLeptons_histo = newTH1D("nLeadingLeptons_histo", "nLeadingLeptons_histo", "nLeadingLepton", 10, 0, 10); 
	nSubLeadingLeptons_histo = newTH1D("nSubLeadingLeptons_histo", "nSubLeadingLeptons_histo", "nSubLeadingLepton", 10, 0, 10); 

	for (unsigned short i(0); i < 11; i++) {
		pt_histo[i] = newTH1D("pt_"+objName[i]+"_histo", "pt_"+objName[i]+"_histo", "p_{T} [GeV] "+objName[i], 200, 0., 1000.);
		eta_histo[i] = newTH1D("eta_"+objName[i]+"_histo", "eta_"+objName[i]+"_histo", "#eta "+objName[i], 100, -2.5, 2.5);
		phi_histo[i] = newTH1D("phi_"+objName[i]+"_histo", "phi_"+objName[i]+"_histo", "#phi "+objName[i], 100, -3.5, 3.5);		
	}

	for (unsigned short j(0); j < 9; j++) {
		for (unsigned short m(0); m < 3; m++) {
			etaSC_histo[j][m] = newTH1D("etaSC_"+eleName[j]+etaName[m]+"_histo", "etaSC_"+eleName[j]+etaName[m]+"_histo", "#eta_{SC} "+eleName[j], 100, -2.5, 2.5);
			isEcalDriven_histo[j][m] = newTH1D("isEcalDriven_"+eleName[j]+etaName[m]+"_histo", "isEcalDriven_"+eleName[j]+etaName[m]+"_histo", "isEcalDriven "+eleName[j], 2, 0., 2.);
			dEtaIn_histo[j][m] = newTH1D("dEtaIn_"+eleName[j]+etaName[m]+"_histo", "dEtaIn_"+eleName[j]+etaName[m]+"_histo", "dEtaIn "+eleName[j], 100, -0.2, 0.2);
			dPhiIn_histo[j][m] = newTH1D("dPhiIn_"+eleName[j]+etaName[m]+"_histo", "dPhiIn_"+eleName[j]+etaName[m]+"_histo", "dPhiIn "+eleName[j], 100, -0.2, 0.2);
			hOverE_histo[j][m] = newTH1D("hOverE_"+eleName[j]+etaName[m]+"_histo", "hOverE_"+eleName[j]+etaName[m]+"_histo", "H/E "+eleName[j], 40, 0., 2.);
			full5x5_r9_histo[j][m] = newTH1D("r9_"+eleName[j]+etaName[m]+"_histo", "r9_"+eleName[j]+etaName[m]+"_histo", "r9 "+eleName[j], 100, 0., 1.05);
			full5x5_sigmaIetaIeta_histo[j][m] = newTH1D("sigmaIetaIeta_"+eleName[j]+etaName[m]+"_histo", "sigmaIetaIeta_"+eleName[j]+etaName[m]+"_histo", "sigmaIetaIeta "+eleName[j], 100, 0., 0.05);
			// full5x5_E5x5_ele_histo[j][m];
			// full5x5_E1x5_ele_histo[j][m];
			// full5x5_E2x5_ele_histo[j][m];
			full5x5_E2x5_Over_E5x5_histo[j][m] = newTH1D("e2x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e2x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e2x5/e5x5 "+eleName[j], 60, 0.5, 1.1);
			full5x5_E1x5_Over_E5x5_histo[j][m] = newTH1D("e1x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e1x5_e5x5_"+eleName[j]+etaName[m]+"_histo", "e1x5/e5x5 "+eleName[j], 100, 0., 1.1);
			EmHadDepth1Iso_histo[j][m] = newTH1D("EmHadDepth1Iso_"+eleName[j]+etaName[m]+"_histo", "EmHadDepth1Iso_"+eleName[j]+etaName[m]+"_histo", "EmHadDepth1Iso "+eleName[j], 200, 0., 200.);
			// ptTracksIso_histo[j][m];
			innerLayerLostHits_histo[j][m] = newTH1D("missingHits_"+eleName[j]+etaName[m]+"_histo", "missingHits_"+eleName[j]+etaName[m]+"_histo", "missingHits "+eleName[j], 10, 0., 10.);
			dxy_histo[j][m] = newTH1D("dxy_"+eleName[j]+etaName[m]+"_histo", "dxy_"+eleName[j]+etaName[m]+"_histo", "|dxy| "+eleName[j], 100, 0., 0.5);
			eOverP[j][m] = newTH1D("eOverP_"+eleName[j]+etaName[m]+"_histo", "eOverP_"+eleName[j]+etaName[m]+"_histo", "E/p "+eleName[j], 100, 0., 100);
		}
	}

	for (unsigned short n(0); n < 5; n++) {
		for (unsigned short l(0); l < 4; l++) {
			mass_dldj_histo[l][n] = newTH1D("mass_dileptondijets_"+dldjName[l]+etaMassName[n]+"_histo", "mass_dileptondijets_"+dldjName[l]+etaMassName[n]+"_histo", "m_{lljj} "+dldjName[l], 300, 0., 6000.);
			mass_dl_histo[l][n] = newTH1D("mass_dileptons_"+dldjName[l]+etaMassName[n]+"_histo", "mass_dileptons_"+dldjName[l]+etaMassName[n]+"_histo", "m_{ll} "+dldjName[l], 100, 0., 2000.);
			mass_dj_histo[l][n] = newTH1D("mass_dijets_"+dldjName[l]+etaMassName[n]+"_histo", "mass_dijets_"+dldjName[l]+etaMassName[n]+"_histo", "m_{jj} "+dldjName[l], 100, 0., 2000.);		
			mass_djLl_histo[l][n] = newTH1D("mass_dijetsLeadingLepton_"+dldjName[l]+etaMassName[n]+"_histo", "mass_dijetsLeadingLepton_"+dldjName[l]+etaMassName[n]+"_histo", "m_{jjl_{L}} "+dldjName[l], 200, 0., 4000.);	
			mass_djSLl_histo[l][n] = newTH1D("mass_dijetsSubLeadingLepton_"+dldjName[l]+etaMassName[n]+"_histo", "mass_dijetsSubLeadingLepton_"+dldjName[l]+etaMassName[n]+"_histo", "m_{jjl_{SL}} "+dldjName[l], 200, 0., 4000.);
		}

		for (unsigned short z(0); z < 5; z++) {
			Zmass_histo[z][n] = newTH1D("Z"+Zname[z]+"_mass"+etaMassName[n]+"_histo", "Z"+Zname[z]+"_mass"+etaMassName[n]+"_histo", "m(Z"+Zname[z]+") [GeV/c^{2}]", 200, 0, 200); 
		}
	}

	for (unsigned short l(0); l < 4; l++) {
		pt_dl_histo[l] = newTH1D("pt_dileptons_"+dldjName[l]+"_histo", "pt_dileptons_"+dldjName[l]+"_histo", "diLepton p_{T} [GeV]", 100, 0., 1000.);
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
							diLepton_pt->at(i), 
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

	if (etaIdx == 0) pt_dl_histo[histoIdx]->Fill(diLeptonDiJets[dldjIdx].diLepton_pt,w);

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

	// for (Long64_t i=0; i<300; i+=1) {
	for (Long64_t i=0; i<nentries; i+=1) {
		Long64_t ientry = LoadTree(i);
		if (ientry < 0) break;
		// cout << "evento " << ientry << endl;
		nEvents++;

		if(fChain->GetEntry(i) == 0){
			std::cerr << "Failed to read Tree entry " << i << "!\n";
			continue;
		}


		if (passEEJJhlt) nEventsPassingEEJJhlt++;
		if (passMMJJhlt) nEventsPassingMMJJhlt++;
		if (passEMJJhlt) nEventsPassingEMJJhlt++;
		if (passTandPEEhlt) nEventsPassingTandPEEhlt++;
		if (passTandPMMhlt) nEventsPassingTandPMMhlt++;


	// --Trigger
		// bool passTrigger = 1; //sempre vero per MC	
		// if (!MC) {
		// 	if (signalEE) passTrigger = passEEJJhlt;   
		// 	if (signalMuMu) passTrigger = passMMJJhlt;
		// 	if (eMuSideband) passTrigger = passEMJJhlt;
		// 	// if (TnPEE) passTrigger = passTandPEEhlt;
		// 	if (TnPEE) passTrigger = passEEJJhlt;
		// 	if (TnPMM) passTrigger = passTandPMMhlt;
		// }

		// if (!passTrigger) continue;  //se evento non passa il trigger passo a quello successivo
		nDataEventsPassingTrigger++;


		// cout << "sumWeights = " << sumWeights << endl;
		// sumWeights += weight;

		w = weight;
		// cout << "weight = " << w << endl;
		if (MCpuReweighted) w *= puweight;
		// cout << "puweight = " << puweight << endl;
		// cout << "weight*puweight = " << w << endl;


	// --Collections
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
		// cout << "nMuons = " << nMuons << endl;    
		unsigned short nJets = jets.size();
		// cout << "nJets = " << nJets << endl;   

		unsigned short nDiLeptonDiJets = diLeptonDiJets.size();
		// cout << "nDiLeptonDiJets = " << nDiLeptonDiJets << endl;

		// unsigned short nLeadingEle = leadingElectrons.size();
		// cout << "nLeadingEle = " << nLeadingEle << endl;
		// unsigned short nSubLeadingEle = subLeadingElectrons.size();
		// cout << "nSubLeadingEle = " << nSubLeadingEle << endl;

		// unsigned short nLeadingMuons = leadingMuons.size();
		// cout << "nLeadingMuons = " << nLeadingMuons << endl;
		// unsigned short nSubLeadingMuons = subLeadingMuons.size();
		// cout << "nSubLeadingMuons = " << nSubLeadingMuons << endl;

		// unsigned short nLeadingJets = leadingJets.size();
		// cout << "nLeadingJets = " << nLeadingJets << endl;
		// unsigned short nSubLeadingJets = subLeadingJets.size();
		// cout << "nSubLeadingJets = " << nSubLeadingJets << endl;



		// if (MC) {   //apply trigger pt cuts on MC
		int elePassingTrigger = 0, elePassingTandPTrigger = 0, muonsPassingTrigger = 0, muPassingEMuTrigger = 0, muPassingTandPTrigger = 0;

		for (unsigned short j(0); j < nElectrons; j++){
			if (electrons[j].v.Pt() > 33) elePassingTrigger++;
			if (electrons[j].v.Pt() > 27) elePassingTandPTrigger++;
		}

		for (unsigned short m(0); m < nMuons; m++){
			if (muons[m].v.Pt() > 50) muonsPassingTrigger++;
			if (muons[m].v.Pt() > 33) muPassingEMuTrigger++;
			if (muons[m].v.Pt() > 27) muPassingTandPTrigger++;
		}

		if (signalEE && elePassingTrigger < 2) continue;
		if (signalMuMu && muonsPassingTrigger < 2) continue;
		if (eMuSideband && (elePassingTrigger < 1 || muPassingEMuTrigger < 1)) continue;
		// if (TnPEE && elePassingTandPTrigger < 2) continue;
		if (TnPEE && elePassingTrigger < 2) continue;
		if (TnPMM && muPassingTandPTrigger < 2) continue;
		// }
		nMCEventsPassingTrigger++;


	// --vtx
		rho_histo->Fill(rho,w);
		nvtx_histo->Fill(nvtx,w);

		for (unsigned short i(0); i < nvtx; i++){
			vtx_x_histo->Fill(vtx_x->at(i),w);
			vtx_y_histo->Fill(vtx_y->at(i),w);
			vtx_z_histo->Fill(vtx_z->at(i),w);
		}


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
		for (unsigned short i(0); i < nMuons; i++){
			doMuonDistributionsPlots(muons,i,1);
			muon_dxy_histo->Fill( mu_dxy->at(i) );
		}


	// --jets
		for (unsigned short l(0); l < nJets; l++){
			doJetsDistributionsPlots(jets,l,2);
		}



	// --Zmass using electrons
		if (TnPEE) {   
			if (nElectrons < 2) continue;
			unsigned int lIdx=0, slIdx=1; //sono ordinati in ordine decrescente in pt

			if (electrons[lIdx].v.Pt() < 35 || electrons[lIdx].v.Pt() < 35) continue;
			nEventsWith2elePassingLoosePreselections++;

			if (electrons[lIdx].charge * electrons[slIdx].charge < 0) {
				nEventsWith2elePassingLoosePreselectionsAndCharge++;	

				if (electrons[lIdx].passCutBasedEleId && electrons[slIdx].passCutBasedEleId) {
					nEventsWith2elePassingLoosePreselectionsAndChargeAndEleId++;

					TLorentzVector lEle, slEle;
					lEle.SetPtEtaPhiE(electrons[lIdx].v.Pt(), electrons[lIdx].v.Eta(), electrons[lIdx].v.Phi(), electrons[lIdx].v.E());
					slEle.SetPtEtaPhiE(electrons[slIdx].v.Pt(), electrons[slIdx].v.Eta(), electrons[slIdx].v.Phi(), electrons[slIdx].v.E());
					float mZ = (lEle+slEle).M();

					Zmass_histo[0][0]->Fill(mZ,w);

					bool isInEBEB = inEBEB(electrons[lIdx].v.Eta(), electrons[slIdx].v.Eta()); 
					bool isInEEEE = inEEEE(electrons[lIdx].v.Eta(), electrons[slIdx].v.Eta());
					bool isInEBEE = inEBEE(electrons[lIdx].v.Eta(), electrons[slIdx].v.Eta()); 

					if (isInEBEB) Zmass_histo[0][1]->Fill(mZ,w);  
					if (isInEEEE) Zmass_histo[0][2]->Fill(mZ,w);  
					if (isInEBEE) Zmass_histo[0][3]->Fill(mZ,w);  
					if (isInEEEE || isInEBEE) Zmass_histo[0][4]->Fill(mZ,w);  
				}
			}
		}




	// --choose dldj candidate
		if (nDiLeptonDiJets < 1) continue;
		nEventsWithAtLeast1DLDJ++;

		int idxDLDJ = 0;

		if (nDiLeptonDiJets > 1) {
			for (unsigned short j(1); j < nDiLeptonDiJets; j++){
				// cout << "sumPt dldj[" << idxDLDJ << "] = " << diLeptonDiJets[idxDLDJ].diLeptonDiJet_sumPt << endl;
				// cout << "sumPt dldj[" << j << "] = " << diLeptonDiJets[j].diLeptonDiJet_sumPt << endl;				
				if (diLeptonDiJets[j].diLeptonDiJet_sumPt > diLeptonDiJets[idxDLDJ].diLeptonDiJet_sumPt) {
					idxDLDJ = j;
					// cout << "prendo idxDLDJ = " << j << endl;
				}
			}
		} 
		// cout << "idxDLDJ = " << idxDLDJ << endl;



	// --select right lepton pair
		if (signalEE || TnPEE) {
			if (diLeptonDiJets[idxDLDJ].isMMJJ || diLeptonDiJets[idxDLDJ].isEMJJ) continue;
		}
		if (signalMuMu || TnPMM) {
			if (diLeptonDiJets[idxDLDJ].isEEJJ || diLeptonDiJets[idxDLDJ].isEMJJ) continue;
		} 
		if (eMuSideband) {
			if (diLeptonDiJets[idxDLDJ].isEEJJ || diLeptonDiJets[idxDLDJ].isMMJJ) continue;
		}
		nEventsWithRightLeptonPair++;



	// --dldj
		bool isInEBEB = inEBEB(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));  //usare eta jets per djets mass ?
		bool isInEEEE = inEEEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ));
		bool isInEBEE = inEBEE(leadingLepton_eta->at(idxDLDJ), subLeadingLepton_eta->at(idxDLDJ)); 


		bool signalRegion=0, flavourSidebandCR=0, lowMlljjCR=0, lowMllCR=0;

		
		if (diLeptonDiJets[idxDLDJ].passPreselections) {
			nEventsWithDLDJpassingPreselections++;

		//-- definition of signal/control regions
			if (signalEE || signalMuMu) {
				if (diLeptonDiJets[idxDLDJ].isSignalRegion) {
					nEventsWithDLDJpassingPreselectionsInSignalRegion++;
					signalRegion=true;
				}
				if (diLeptonDiJets[idxDLDJ].isLowMlljjCR) {
					nEventsWithDLDJpassingPreselectionsInLowMlljjCR++;
					lowMlljjCR=true;
				}
				if (diLeptonDiJets[idxDLDJ].isLowMllCR) {
					nEventsWithDLDJpassingPreselectionsInLowMllCR++;
					lowMllCR=true;
				}
			}
			if (eMuSideband && diLeptonDiJets[idxDLDJ].isSignalRegion) {
				nEventsWithDLDJpassingPreselectionsInFlavourSidebandCR++;
				flavourSidebandCR=true;  
			}
		//--
		
			if ( (signalEE && leadingElectrons[idxDLDJ].passCutBasedEleId && subLeadingElectrons[idxDLDJ].passCutBasedEleId) || (signalMuMu && leadingMuons[idxDLDJ].isHighPt && subLeadingMuons[idxDLDJ].isHighPt) ) {
			// if ( (signalEE && leadingElectrons[idxDLDJ].passHEEPId && subLeadingElectrons[idxDLDJ].passHEEPId) || (signalMuMu && leadingMuons[idxDLDJ].isHighPt && subLeadingMuons[idxDLDJ].isHighPt) ) {
				nEventsWithDLDJpassingSelections++;

				if (signalRegion) {		
					nEventsWithDLDJpassingSelectionsInSignalRegion++;	
					doMassPlots(diLeptonDiJets,idxDLDJ,0,0);
					if (isInEBEB) doMassPlots(diLeptonDiJets,idxDLDJ,0,1);
					if (isInEEEE) doMassPlots(diLeptonDiJets,idxDLDJ,0,2);
					if (isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,0,3);
					if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,0,4);
				}
				if (lowMlljjCR) {
					nEventsWithDLDJpassingSelectionsInLowMlljjCR++;
					doMassPlots(diLeptonDiJets,idxDLDJ,1,0);
					if (isInEBEB) doMassPlots(diLeptonDiJets,idxDLDJ,1,1);
					if (isInEEEE) doMassPlots(diLeptonDiJets,idxDLDJ,1,2);
					if (isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,1,3);
					if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,1,4);
				}
				if (lowMllCR) {
					nEventsWithDLDJpassingSelectionsInLowMllCR++;
					doMassPlots(diLeptonDiJets,idxDLDJ,2,0);
					if (isInEBEB) doMassPlots(diLeptonDiJets,idxDLDJ,2,1);
					if (isInEEEE) doMassPlots(diLeptonDiJets,idxDLDJ,2,2);
					if (isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,2,3);
					if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,2,4);
				}
			}

			if ( flavourSidebandCR && ((leadingElectrons[idxDLDJ].passHEEPId && subLeadingMuons[idxDLDJ].isHighPt) || (subLeadingElectrons[idxDLDJ].passHEEPId && leadingMuons[idxDLDJ].isHighPt)) ) {
				nEventsWithDLDJpassingSelectionsInFlavourSidebandCR++;
				doMassPlots(diLeptonDiJets,idxDLDJ,3,0);
				if (isInEBEB) doMassPlots(diLeptonDiJets,idxDLDJ,3,1);
				if (isInEEEE) doMassPlots(diLeptonDiJets,idxDLDJ,3,2);
				if (isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,3,3);
				if (isInEEEE || isInEBEE) doMassPlots(diLeptonDiJets,idxDLDJ,3,4);
			}



			int nLeadingLeptons=0, nSubLeadingLeptons=0;

		// --leadingEle
			if ( signalEE || (eMuSideband && leadingElectrons[idxDLDJ].etaSC != -999) ) {
				if (leadingElectrons[idxDLDJ].v.Pt() < 60) continue;   //already in preselections
				nLeadingLeptons++;
				nLeadingEle++;

				// if (leadingElectrons[idxDLDJ].passHEEPId) {
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,3); 
				doEleDistributionsPlots(leadingElectrons,idxDLDJ,7);
				// }

				bool lInEB = isInEB(leadingElectrons[idxDLDJ].v.Eta());  
				bool lInEE = isInEE(leadingElectrons[idxDLDJ].v.Eta());

				doElePlots(leadingElectrons,idxDLDJ,1,0);
				if (lInEB) doElePlots(leadingElectrons,idxDLDJ,1,1);
				if (lInEE) doElePlots(leadingElectrons,idxDLDJ,1,2);
	   
				if (leadingElectrons[idxDLDJ].passHEEPId) {
					nLeadingElePassingHEEPId++;
					doElePlots(leadingElectrons,idxDLDJ,4,0);
					if (lInEB) doElePlots(leadingElectrons,idxDLDJ,4,1);
					if (lInEE) doElePlots(leadingElectrons,idxDLDJ,4,2);
				}

				if (leadingElectrons[idxDLDJ].passCutBasedEleId) {
					nLeadingElePassingEleId++;
					doElePlots(leadingElectrons,idxDLDJ,7,0);
					if (lInEB) doElePlots(leadingElectrons,idxDLDJ,7,1);
					if (lInEE) doElePlots(leadingElectrons,idxDLDJ,7,2);
				}
			}


		// --subLeadingEle
			if ( signalEE || (eMuSideband && subLeadingElectrons[idxDLDJ].etaSC != -999) ) {
				if (subLeadingElectrons[idxDLDJ].v.Pt() < 53) continue;   //already in preselections
				nSubLeadingLeptons++;
				nSubLeadingEle++;

				// if (subLeadingElectrons[idxDLDJ].passHEEPId) {
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,4); 
				doEleDistributionsPlots(subLeadingElectrons,idxDLDJ,8);
				// }

				bool slInEB = isInEB(subLeadingElectrons[idxDLDJ].v.Eta());  
				bool slInEE = isInEE(subLeadingElectrons[idxDLDJ].v.Eta());

				doElePlots(subLeadingElectrons,idxDLDJ,2,0);
				if (slInEB) doElePlots(subLeadingElectrons,idxDLDJ,2,1);
				if (slInEE) doElePlots(subLeadingElectrons,idxDLDJ,2,2); 
 
				if (subLeadingElectrons[idxDLDJ].passHEEPId) {
					nSubLeadingElePassingHEEPId++;
					doElePlots(subLeadingElectrons,idxDLDJ,5,0);
					if (slInEB) doElePlots(subLeadingElectrons,idxDLDJ,5,1);
					if (slInEE) doElePlots(subLeadingElectrons,idxDLDJ,5,2);
				}

				if (subLeadingElectrons[idxDLDJ].passCutBasedEleId) {
					nSubLeadingElePassingEleId++;
					doElePlots(subLeadingElectrons,idxDLDJ,8,0);
					if (slInEB) doElePlots(subLeadingElectrons,idxDLDJ,8,1);
					if (slInEE) doElePlots(subLeadingElectrons,idxDLDJ,8,2);
				}
			}


		// --leadingMuon
			if ( signalMuMu || (eMuSideband && subLeadingElectrons[idxDLDJ].etaSC != -999) ) {
				if (leadingMuons[idxDLDJ].v.Pt() < 60) continue;  //already in preselections
				nLeadingLeptons++;
				nLeadingMuon++;
				if (leadingMuons[idxDLDJ].isHighPt) { //da togliere?
					nLeadingMuonPassingHighPt++;
					doMuonDistributionsPlots(leadingMuons,idxDLDJ,3);  
					doMuonDistributionsPlots(leadingMuons,idxDLDJ,9);
				}
			}


		// --subLeadingMuon
			if ( signalMuMu || (eMuSideband && leadingElectrons[idxDLDJ].etaSC != -999) ) {
				if (subLeadingMuons[idxDLDJ].v.Pt() < 53) continue;  //already in preselections
				nSubLeadingLeptons++;
				nSubLeadingMuon++;
				if (subLeadingMuons[idxDLDJ].isHighPt) {  //da togliere?
					nSubLeadingMuonPassingHighPt++;
					doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,4);  
					doMuonDistributionsPlots(subLeadingMuons,idxDLDJ,10);
				}
			}


		// --leadingJet
			if (signalEE || signalMuMu || eMuSideband) {
				if (leadingJets[idxDLDJ].v.Pt() < 40) continue;  //already in preselections and miniTree	
				nLeadingJet++;
				doJetsDistributionsPlots(leadingJets,idxDLDJ,5);
			}


		// --subLeadingJet
			if (signalEE || signalMuMu || eMuSideband) {
				if (subLeadingJets[idxDLDJ].v.Pt() < 40) continue;   //already in preselections and miniTree
				nSubLeadingJet++;
				doJetsDistributionsPlots(subLeadingJets,idxDLDJ,6);
			}


			nLeadingLeptons_histo -> Fill(nLeadingLeptons);
			nSubLeadingLeptons_histo -> Fill(nSubLeadingLeptons);

		}//fine if su passPreselections



	// --Zmass
		if (TnPEE) {   //Z->ee 

			if (leadingElectrons[idxDLDJ].v.Pt() < 35 || subLeadingElectrons[idxDLDJ].v.Pt() < 35) continue;  //already in miniTree
			nEventsWithEEJJpassingLoosePreselections++; //=nEventsWithRightLeptonPair

			if (leadingElectrons[idxDLDJ].charge * subLeadingElectrons[idxDLDJ].charge < 0 ) {	
				nEventsWithEEJJpassingLoosePreselectionsAndCharge++;		
				// doZmassPlots(diLeptonDiJets,idxDLDJ,0,0);
				// if (isInEBEB) doZmassPlots(diLeptonDiJets,idxDLDJ,0,1);
				// if (isInEEEE) doZmassPlots(diLeptonDiJets,idxDLDJ,0,2);
				// if (isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,0,3);
				// if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,0,4);

				if (leadingElectrons[idxDLDJ].passHEEPId && subLeadingElectrons[idxDLDJ].passHEEPId) {
					nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId++;
					doZmassPlots(diLeptonDiJets,idxDLDJ,1,0);
					if (isInEBEB) doZmassPlots(diLeptonDiJets,idxDLDJ,1,1);
					if (isInEEEE) doZmassPlots(diLeptonDiJets,idxDLDJ,1,2);
					if (isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,1,3);
					if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,1,4);
				}

				if (leadingElectrons[idxDLDJ].passCutBasedEleId && subLeadingElectrons[idxDLDJ].passCutBasedEleId) {
					nEventsWithEEJJpassingLoosePreselectionsAndChargeAndEleId++;
					doZmassPlots(diLeptonDiJets,idxDLDJ,2,0);
					if (isInEBEB) doZmassPlots(diLeptonDiJets,idxDLDJ,2,1);
					if (isInEEEE) doZmassPlots(diLeptonDiJets,idxDLDJ,2,2);
					if (isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,2,3);
					if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,2,4);
				}
			}
		}

		if (TnPMM) {   //Z->mumu
				
			if (leadingMuons[idxDLDJ].v.Pt() < 35 || subLeadingMuons[idxDLDJ].v.Pt() < 35) continue;  //already in miniTree
			nEventsWithMMJJpassingLoosePreselections++; //=nEventsWithRightLeptonPair

			if (leadingMuons[idxDLDJ].charge * subLeadingMuons[idxDLDJ].charge < 0 ) { 
				nEventsWithMMJJpassingLoosePreselectionsAndCharge++;
				doZmassPlots(diLeptonDiJets,idxDLDJ,3,0);
				if (isInEBEB) doZmassPlots(diLeptonDiJets,idxDLDJ,3,1);
				if (isInEEEE) doZmassPlots(diLeptonDiJets,idxDLDJ,3,2);
				if (isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,3,3);
				if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,3,4);

				if (leadingMuons[idxDLDJ].isHighPt && subLeadingMuons[idxDLDJ].isHighPt) {
					nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt++;
					doZmassPlots(diLeptonDiJets,idxDLDJ,4,0);
					if (isInEBEB) doZmassPlots(diLeptonDiJets,idxDLDJ,4,1);
					if (isInEEEE) doZmassPlots(diLeptonDiJets,idxDLDJ,4,2);
					if (isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,4,3);
					if (isInEEEE || isInEBEE) doZmassPlots(diLeptonDiJets,idxDLDJ,4,4);
				}
			}
		}


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

	if (TnPEE || TnPMM) {
		TFile f(outputdir+suff+"/"+"Zmass_histos"+suff+".root","recreate");
		for (unsigned short l(0); l < numbOfHistograms; l++){
			string histoName = listOfHistograms[l]->GetName();
			if (histoName.find("Z") != string::npos) listOfHistograms[l]->Write(); 
		}

	} else {
		TFile f(outputdir+suff+"/"+"distributions_histos"+suff+".root","recreate");
		for (unsigned short l(0); l < numbOfHistograms; l++){
			string histoName = listOfHistograms[l]->GetName();
			if ( !(histoName.find("Z") != string::npos) ) listOfHistograms[l]->Write(); 

		// TCanvas* canv = new TCanvas();
		// gPad->SetRightMargin(0.05);
		// // CMS_lumi(canv, 4, 33);
		// listOfHistograms[l]->Draw();
		// canv->SaveAs(outputdir+suff+"/"+histoName+".png");
		// delete canv;
		// // }
		}
	}


	ofstream fOutput;
	if (TnPEE || TnPMM) fOutput.open(outputdir+suff+"/nEventsTnP"+suff+".dat");
	else fOutput.open(outputdir+suff+"/nEvents"+suff+".dat");
	fOutput << "nEvents = " << nEvents << " \n";
	fOutput << "nDataEventsPassingTrigger = " << nDataEventsPassingTrigger << " \n";
	fOutput << "nMCEventsPassingTrigger = " << nMCEventsPassingTrigger << " \n";
	fOutput << "nEventsPassingEEJJhlt = " << nEventsPassingEEJJhlt << " \n";
	fOutput << "nEventsPassingMMJJhlt = " << nEventsPassingMMJJhlt << " \n";
	fOutput << "nEventsPassingEMJJhlt = " << nEventsPassingEMJJhlt << " \n";
	fOutput << "nEventsPassingTandPEEhlt = " << nEventsPassingTandPEEhlt << " \n";
	fOutput << "nEventsPassingTandPMMhlt = " << nEventsPassingTandPMMhlt << " \n";
	fOutput << "nEventsWithAtLeast1DLDJ = " << nEventsWithAtLeast1DLDJ << " \n";
	fOutput << "nEventsWithRightLeptonPair = " << nEventsWithRightLeptonPair << " \n";
	fOutput << "nEventsWithDLDJpassingPreselections = " << nEventsWithDLDJpassingPreselections << " \n";
	fOutput << "nEventsWithDLDJpassingPreselectionsInSignalRegion = " << nEventsWithDLDJpassingPreselectionsInSignalRegion << " \n";
	fOutput << "nEventsWithDLDJpassingPreselectionsInLowMlljjCR = " << nEventsWithDLDJpassingPreselectionsInLowMlljjCR << " \n";
	fOutput << "nEventsWithDLDJpassingPreselectionsInLowMllCR = " << nEventsWithDLDJpassingPreselectionsInLowMllCR << " \n";
	fOutput << "nEventsWithDLDJpassingPreselectionsInFlavourSidebandCR = " << nEventsWithDLDJpassingPreselectionsInFlavourSidebandCR << " \n";
	fOutput << "nEventsWithDLDJpassingSelections = " << nEventsWithDLDJpassingSelections << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInSignalRegion = " << nEventsWithDLDJpassingSelectionsInSignalRegion << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInLowMlljjCR = " << nEventsWithDLDJpassingSelectionsInLowMlljjCR << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInLowMllCR = " << nEventsWithDLDJpassingSelectionsInLowMllCR << " \n";
	fOutput << "nEventsWithDLDJpassingSelectionsInFlavourSidebandCR = " << nEventsWithDLDJpassingSelectionsInFlavourSidebandCR << " \n";
	fOutput << "nLeadingEle = " << nLeadingEle << " \n";
	fOutput << "nLeadingElePassingHEEPId = " << nLeadingElePassingHEEPId << " \n";
	fOutput << "nLeadingElePassingEleId = " << nLeadingElePassingEleId << " \n";
	fOutput << "nSubLeadingEle = " << nSubLeadingEle << " \n";
	fOutput << "nSubLeadingElePassingHEEPId = " << nSubLeadingElePassingHEEPId << " \n";
	fOutput << "nSubLeadingElePassingEleId = " << nSubLeadingElePassingEleId << " \n";
	fOutput << "nLeadingMuon = " << nLeadingMuon << " \n";
	fOutput << "nLeadingMuonPassingHighPt = " << nLeadingMuonPassingHighPt << " \n";
	fOutput << "nSubLeadingMuon = " << nSubLeadingMuon << " \n";
	fOutput << "nSubLeadingMuonPassingHighPt = " << nSubLeadingMuonPassingHighPt << " \n";
	fOutput << "nLeadingJet = " << nLeadingJet << " \n";
	fOutput << "nSubLeadingJet = " << nSubLeadingJet << " \n";
	fOutput << "nEventsWithEEJJpassingLoosePreselections = " << nEventsWithEEJJpassingLoosePreselections << " \n";
	fOutput << "nEventsWithMMJJpassingLoosePreselections = " << nEventsWithMMJJpassingLoosePreselections << " \n";
	fOutput << "nEventsWithEEJJpassingLoosePreselectionsAndCharge = " << nEventsWithEEJJpassingLoosePreselectionsAndCharge << " \n";
	fOutput << "nEventsWithMMJJpassingLoosePreselectionsAndCharge = " << nEventsWithMMJJpassingLoosePreselectionsAndCharge << " \n";
	fOutput << "nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId = " << nEventsWithEEJJpassingLoosePreselectionsAndChargeAndHEEPId << " \n";
	fOutput << "nEventsWithEEJJpassingLoosePreselectionsAndChargeAndEleId = " << nEventsWithEEJJpassingLoosePreselectionsAndChargeAndEleId << " \n";
	fOutput << "nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt = " << nEventsWithMMJJpassingLoosePreselectionsAndChargeAndHighPt << " \n";
	fOutput << "nEventsWith2elePassingLoosePreselections = " << nEventsWith2elePassingLoosePreselections << " \n";
	fOutput << "nEventsWith2elePassingLoosePreselectionsAndCharge = " << nEventsWith2elePassingLoosePreselectionsAndCharge << " \n";
	fOutput << "nEventsWith2elePassingLoosePreselectionsAndChargeAndEleId = " << nEventsWith2elePassingLoosePreselectionsAndChargeAndEleId << " \n";
	// fOutput << " = " <<  << " \n";	
	fOutput.close();

}
// ******************************************************************************************



makePlots::makePlots(TString filename_, TString outputdir_, bool MC_, bool MCpuReweighted_, bool signalEE_, bool signalMuMu_, bool eMuSideband_, bool TnPEE_, bool TnPMM_):
	filename(filename_), outputdir(outputdir_), MC(MC_), MCpuReweighted(MCpuReweighted_), signalEE(signalEE_), signalMuMu(signalMuMu_), eMuSideband(eMuSideband_), TnPEE(TnPEE_), TnPMM(TnPMM_)
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

	SetHistos();
	Loop();
	saveHistosAndOutputFile(outputdir);
}



void runMakePlots() {

	string inputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_20/src/dafne/miniTrees-DoubleEG/";   
	string outputDir = "/home/gpfs/manip/mnt/cms/gnegro/CMSSW_8_0_20/src/dafne/DistributionsForDoubleEG/";

	bool MCpuReweighted = true;  // false quando giro MC bkg per confronto con SingleMuon (pureweighting gia` fatto in tree)

	bool signalEE = true;  //true for DoubleEG (if not TnP)
	bool signalMuMu = false;  //true for SingleMuon (if not TnP)
	bool eMuSideband = false; //true for MuonEG

	bool TnPEE = false; //true for DoubleEG 
	bool TnPMM = false; //true for SingleMuon 


	makePlots(inputDir+"cmsWR2016-DY/output_DYJetsToLL-amcatnloFXFX_miniTree.root", outputDir+"DYJetsToLL-amcatnlo", true, MCpuReweighted, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots(inputDir+"cmsWR2016-TT/output_TTJets_miniTree.root", outputDir+"TTJets", true, MCpuReweighted, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots(inputDir+"cmsWR2016-WJets/output_WJetsToLNu_miniTree.root", outputDir+"WJets", true, MCpuReweighted, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots(inputDir+"cmsWR2016-WZ/output_WZ_miniTree.root", outputDir+"WZ", true, MCpuReweighted, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);

	// makePlots(inputDir+"cmsWR2016-ZZ/output_ZZ_miniTree.root", outputDir+"ZZ", true, MCpuReweighted, signalEE, signalMuMu, eMuSideband, TnPEE, TnPMM);


	// makePlots(inputDir+"output_DoubleEG_miniTree.root", outputDir+"DoubleEG", false, false, signalEE, false, false, TnPEE, false);


}
