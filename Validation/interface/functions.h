#ifndef _functions_h_
#define _functions_h_
#include <iostream>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <TLorentzVector.h>

using namespace std;



struct eleStruct{
	eleStruct();
	eleStruct(float pt_, float eta_, float phi_, float en_, float idmva_, unsigned int id_, float iso_, 
		float dzvtx_, float d0vtx_, bool passHEEPId_, bool passCutBasedEleId_, int isMatchedToGen_, 
		int charge_, float etaSC_, bool isEcalDriven_, float dEtaIn_, float dPhiIn_, float HoE_, float r9_, 
		float sigmaIetaIeta_, float e5x5_, float e1x5_, float e2x5_, float e2x5_e5x5_, float e1x5_e5x5_, 
		float EmHadDepth1Iso_, float ptTracksIso_, int missingHits_, float dxy_, float eOverP_) {
		v.SetPtEtaPhiE(pt_, eta_, phi_, en_);
		idmva = idmva_;
		id = id_;
		iso = iso_;
		dzvtx = dzvtx_;
		d0vtx = d0vtx_;
		passHEEPId = passHEEPId_;
		passCutBasedEleId = passCutBasedEleId_;
		isMatchedToGen = isMatchedToGen_;
		charge = charge_;
		etaSC = etaSC_;
		isEcalDriven = isEcalDriven_;
		dEtaIn = dEtaIn_;
		dPhiIn = dPhiIn_;
		HoE = HoE_;		
		r9 = r9_;
		sigmaIetaIeta = sigmaIetaIeta_;
		e5x5 = e5x5_;
		e1x5 = e1x5_;
		e2x5 = e2x5_;
		e2x5_e5x5 = e2x5_e5x5_;
		e1x5_e5x5 = e1x5_e5x5_;
		EmHadDepth1Iso = EmHadDepth1Iso_;
		ptTracksIso = ptTracksIso_;
		missingHits = missingHits_;
		dxy = dxy_;
		eOverP = eOverP_;
	}
	TLorentzVector v;
	float idmva, iso, dzvtx, d0vtx, etaSC, dEtaIn, dPhiIn, HoE, r9, sigmaIetaIeta, e5x5, e1x5, e2x5, e2x5_e5x5, 
	e1x5_e5x5, EmHadDepth1Iso, ptTracksIso, dxy, eOverP; 
	bool passHEEPId, passCutBasedEleId, isEcalDriven;
	int isMatchedToGen, charge, missingHits;
	unsigned int id;
};



struct diLeptonDiJetStruct{
	diLeptonDiJetStruct();
	diLeptonDiJetStruct(bool isEEJJ_, bool isEETT_, bool isMMJJ_, bool isMMTT_, bool isEMJJ_, 
		bool isSignalRegion_, bool isLowMllCR_, bool isLowMlljjCR_, bool isBB_, bool isEE_, bool isEB_, 
		bool passPreselections_, float dRLeadLeptonLeadJet_, float dRLeadLeptonSubLeadJet_, 
		float dRSubLeadLeptonLeadJet_, float dRSubLeadLeptonSubLeadJet_, int diLeptonDiJet_vtxIndex_, 
		float diLeptonDiJet_sumPt_, float diLeptonDiJet_invMass_, float diLepton_invMass_, float diLepton_pt_, 
		float diJet_invMass_, float diJetLeadingLepton_invMass_, float diJetSubLeadingLepton_invMass_){
		isEEJJ = isEEJJ_;
		isEETT = isEETT_;
		isMMJJ = isMMJJ_;
		isMMTT = isMMTT_;
		isEMJJ = isEMJJ_;
		isSignalRegion = isSignalRegion_;
		isLowMllCR = isLowMllCR_;
		isLowMlljjCR = isLowMlljjCR_;
		isBB = isBB_;
		isEE = isEE_;
		isEB = isEB_;
		passPreselections = passPreselections_;
		dRLeadLeptonLeadJet = dRLeadLeptonLeadJet_;
		dRLeadLeptonSubLeadJet = dRLeadLeptonSubLeadJet_;
		dRSubLeadLeptonLeadJet = dRSubLeadLeptonLeadJet_;
		dRSubLeadLeptonSubLeadJet = dRSubLeadLeptonSubLeadJet_;
		diLeptonDiJet_vtxIndex = diLeptonDiJet_vtxIndex_;
		diLeptonDiJet_sumPt = diLeptonDiJet_sumPt_;
		diLeptonDiJet_invMass = diLeptonDiJet_invMass_;
		diLepton_invMass = diLepton_invMass_;
		diLepton_pt = diLepton_pt_;
		diJet_invMass = diJet_invMass_;
		diJetLeadingLepton_invMass = diJetLeadingLepton_invMass_;
		diJetSubLeadingLepton_invMass = diJetSubLeadingLepton_invMass_;		
	}
	bool isEEJJ, isEETT, isMMJJ, isMMTT, isEMJJ, isSignalRegion, isLowMllCR, isLowMlljjCR, isBB, isEE, isEB, 
	passPreselections; 
	float dRLeadLeptonLeadJet, dRLeadLeptonSubLeadJet, dRSubLeadLeptonLeadJet, dRSubLeadLeptonSubLeadJet, 
	diLeptonDiJet_sumPt, diLeptonDiJet_invMass, diLepton_invMass, diLepton_pt, diJet_invMass, 
	diJetLeadingLepton_invMass, diJetSubLeadingLepton_invMass;
	int diLeptonDiJet_vtxIndex; 
};


struct muonStruct{
	muonStruct();
	muonStruct(float pt_, float eta_, float phi_, float en_, float iso_, bool isHighPt_, int isMatchedToGen_, 
		int charge_, float dzvtx_, float dxy_) {
		v.SetPtEtaPhiE(pt_, eta_, phi_, en_);
		iso = iso_;
		isHighPt = isHighPt_;
		isMatchedToGen = isMatchedToGen_;
		charge = charge_;
		dzvtx = dzvtx_;
		dxy = dxy_;
	}
	TLorentzVector v;
	float iso, dzvtx, dxy; 
	bool isHighPt;
	int isMatchedToGen, charge;
};


struct jetStruct{
	jetStruct();
	jetStruct(double pt_, double eta_, double phi_, double en_, int isMatchedToGen_) {
		v.SetPtEtaPhiE(pt_, eta_, phi_,en_);
		isMatchedToGen = isMatchedToGen_;
	}
	TLorentzVector v;
	int isMatchedToGen;
};



bool isInEB(const float& eta){
	bool isInEB = fabs(eta) < 1.4442;
	return isInEB;
}


bool isInEE(const float& eta){
	bool isInEE = fabs(eta) > 1.566  &&  fabs(eta) < 2.5;
	return isInEE;
}



bool inEBEB(const float& eta1, const float& eta2){
	bool inEBEB = isInEB(eta1) && isInEB(eta2);
	return inEBEB;
}

bool inEEEE(const float& eta1, const float& eta2){
	bool inEEEE = isInEE(eta1) && isInEE(eta2);
	return inEEEE;
}

bool inEBEE(const float& eta1, const float& eta2){
	bool inEBEE = ((isInEB(eta1) && isInEE(eta2)) || (isInEE(eta1) && isInEB(eta2)));
	return inEBEE;
}



double phi0to2pi(double phi){
	double pi = 3.141592653589793238;
	while (phi >= 2.*pi) phi -= 2.*pi;
	while (phi < 0.) phi += 2.*pi;
	return phi;
}


double deltaPhi(TLorentzVector v1, TLorentzVector v2){
	// build the delta Phi angle between the two vectors
	// double pi = 3.141592653589793238;
	double dPhi = v1.Phi() - v2.Phi();
	return phi0to2pi(dPhi);
} 


double deltaR(TLorentzVector v1, TLorentzVector v2){
	double dEta = v1.Eta() - v2.Eta();
	double dPhi = deltaPhi(v1, v2);
	return sqrt(dEta * dEta + dPhi * dPhi);
}


float eff_sigma(std::vector<float> & v)
{
		size_t n = v.size();
		if (n < 2) return 0;
		std::sort(v.begin(), v.end());
		int s = floor(0.68269 * n);
		float d_min = v[s] - v[0];
		for (size_t i = s; i < n; ++i) {
				float d = v[i] - v[i - s];
				if (d < d_min) d_min = d;
		}
		return d_min / 2.;
}
#endif