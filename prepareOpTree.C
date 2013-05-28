// C++ headers
#include <iostream>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;

int main()
{
// /store/group/phys_higgs/Resonant_HH/trees/radion_tree_v04/Graviton_Radion-nm.root
//	TFile *infile = TFile::Open("root://eoscms//eos/cms/store/group/phys_higgs/Resonant_HH/trees/radion_tree_v04/Graviton_Radion-nm.root");
	TFile *infile = TFile::Open("histograms_CMS-HGG_1.root");
	TFile *outfile = new TFile("test.root", "RECREATE");
	TTree *intree = (TTree*)infile->Get("Radion_m300_8TeV_nm");
	TTree *outtree = new TTree("jets", "jets");

	cout << "#entries= " << intree->GetEntries() << endl;
	float gr_radion_p4_pt, gr_radion_p4_eta, gr_radion_p4_phi, gr_radion_p4_mass, gr_hgg_p4_pt, gr_hgg_p4_eta, gr_hgg_p4_phi, gr_hgg_p4_mass, gr_hbb_p4_pt, gr_hbb_p4_eta, gr_hbb_p4_phi, gr_hbb_p4_mass, gr_g1_p4_pt, gr_g1_p4_eta, gr_g1_p4_phi, gr_g1_p4_mass, gr_g2_p4_pt, gr_g2_p4_eta, gr_g2_p4_phi, gr_g2_p4_mass, gr_b1_p4_pt, gr_b1_p4_eta, gr_b1_p4_phi, gr_b1_p4_mass, gr_b2_p4_pt, gr_b2_p4_eta, gr_b2_p4_phi, gr_b2_p4_mass, gr_j1_p4_pt, gr_j1_p4_eta, gr_j1_p4_phi, gr_j1_p4_mass, gr_j2_p4_pt, gr_j2_p4_eta, gr_j2_p4_phi, gr_j2_p4_mass;
	float jet_DR_b1, jet_DR_b2, jet_DR_j1, jet_DR_j2;
	float j1_e, j1_pt, j1_phi, j1_eta, j1_beta, j1_betaStar, j1_betaStarClassic, j1_dR2Mean, j1_csvBtag, j1_csvMvaBtag, j1_jetProbBtag, j1_tcheBtag, j1_radionMatched, j1_ptD, j1_nSecondaryVertices, j1_secVtxPt, j1_secVtx3dL, j1_secVtx3deL, j1_emfrac, j1_hadfrac, j1_ntk, j1_nNeutrals, j1_nCharged, j1_axis1, j1_axis2, j1_pull, j1_Rchg, j1_Rneutral, j1_R;
	float j2_e, j2_pt, j2_phi, j2_eta, j2_beta, j2_betaStar, j2_betaStarClassic, j2_dR2Mean, j2_csvBtag, j2_csvMvaBtag, j2_jetProbBtag, j2_tcheBtag, j2_radionMatched, j2_ptD, j2_nSecondaryVertices, j2_secVtxPt, j2_secVtx3dL, j2_secVtx3deL, j2_emfrac, j2_hadfrac, j2_ntk, j2_nNeutrals, j2_nCharged, j2_axis1, j2_axis2, j2_pull, j2_Rchg, j2_Rneutral, j2_R;
	float j3_e, j3_pt, j3_phi, j3_eta, j3_beta, j3_betaStar, j3_betaStarClassic, j3_dR2Mean, j3_csvBtag, j3_csvMvaBtag, j3_jetProbBtag, j3_tcheBtag, j3_radionMatched, j3_ptD, j3_nSecondaryVertices, j3_secVtxPt, j3_secVtx3dL, j3_secVtx3deL, j3_emfrac, j3_hadfrac, j3_ntk, j3_nNeutrals, j3_nCharged, j3_axis1, j3_axis2, j3_pull, j3_Rchg, j3_Rneutral, j3_R;
	float j4_e, j4_pt, j4_phi, j4_eta, j4_beta, j4_betaStar, j4_betaStarClassic, j4_dR2Mean, j4_csvBtag, j4_csvMvaBtag, j4_jetProbBtag, j4_tcheBtag, j4_radionMatched, j4_ptD, j4_nSecondaryVertices, j4_secVtxPt, j4_secVtx3dL, j4_secVtx3deL, j4_emfrac, j4_hadfrac, j4_ntk, j4_nNeutrals, j4_nCharged, j4_axis1, j4_axis2, j4_pull, j4_Rchg, j4_Rneutral, j4_R;
	float jet_e, jet_pt, jet_phi, jet_eta, jet_beta, jet_betaStar, jet_betaStarClassic, jet_dR2Mean, jet_csvBtag, jet_csvMvaBtag, jet_jetProbBtag, jet_tcheBtag, jet_radionMatched, jet_ptD, jet_nSecondaryVertices, jet_secVtxPt, jet_secVtx3dL, jet_secVtx3deL, jet_emfrac, jet_hadfrac, jet_ntk, jet_nNeutrals, jet_nCharged, jet_axis1, jet_axis2, jet_pull, jet_Rchg, jet_Rneutral, jet_R;
	int jet_index1, jet_index2;

	intree->SetBranchAddress("gr_radion_p4_pt", &gr_radion_p4_pt);
	intree->SetBranchAddress("gr_radion_p4_eta", &gr_radion_p4_eta);
	intree->SetBranchAddress("gr_radion_p4_phi", &gr_radion_p4_phi);
	intree->SetBranchAddress("gr_radion_p4_mass", &gr_radion_p4_mass);
	intree->SetBranchAddress("gr_hgg_p4_pt", &gr_hgg_p4_pt);
	intree->SetBranchAddress("gr_hgg_p4_eta", &gr_hgg_p4_eta);
	intree->SetBranchAddress("gr_hgg_p4_phi", &gr_hgg_p4_phi);
	intree->SetBranchAddress("gr_hgg_p4_mass", &gr_hgg_p4_mass);
	intree->SetBranchAddress("gr_hbb_p4_pt", &gr_hbb_p4_pt);
	intree->SetBranchAddress("gr_hbb_p4_eta", &gr_hbb_p4_eta);
	intree->SetBranchAddress("gr_hbb_p4_phi", &gr_hbb_p4_phi);
	intree->SetBranchAddress("gr_hbb_p4_mass", &gr_hbb_p4_mass);
	intree->SetBranchAddress("gr_g1_p4_pt", &gr_g1_p4_pt);
	intree->SetBranchAddress("gr_g1_p4_eta", &gr_g1_p4_eta);
	intree->SetBranchAddress("gr_g1_p4_phi", &gr_g1_p4_phi);
	intree->SetBranchAddress("gr_g1_p4_mass", &gr_g1_p4_mass);
	intree->SetBranchAddress("gr_g2_p4_pt", &gr_g2_p4_pt);
	intree->SetBranchAddress("gr_g2_p4_eta", &gr_g2_p4_eta);
	intree->SetBranchAddress("gr_g2_p4_phi", &gr_g2_p4_phi);
	intree->SetBranchAddress("gr_g2_p4_mass", &gr_g2_p4_mass);
	intree->SetBranchAddress("gr_b1_p4_pt", &gr_b1_p4_pt);
	intree->SetBranchAddress("gr_b1_p4_eta", &gr_b1_p4_eta);
	intree->SetBranchAddress("gr_b1_p4_phi", &gr_b1_p4_phi);
	intree->SetBranchAddress("gr_b1_p4_mass", &gr_b1_p4_mass);
	intree->SetBranchAddress("gr_b2_p4_pt", &gr_b2_p4_pt);
	intree->SetBranchAddress("gr_b2_p4_eta", &gr_b2_p4_eta);
	intree->SetBranchAddress("gr_b2_p4_phi", &gr_b2_p4_phi);
	intree->SetBranchAddress("gr_b2_p4_mass", &gr_b2_p4_mass);
	intree->SetBranchAddress("gr_j1_p4_pt", &gr_j1_p4_pt);
	intree->SetBranchAddress("gr_j1_p4_eta", &gr_j1_p4_eta);
	intree->SetBranchAddress("gr_j1_p4_phi", &gr_j1_p4_phi);
	intree->SetBranchAddress("gr_j1_p4_mass", &gr_j1_p4_mass);
	intree->SetBranchAddress("gr_j2_p4_pt", &gr_j2_p4_pt);
	intree->SetBranchAddress("gr_j2_p4_eta", &gr_j2_p4_eta);
	intree->SetBranchAddress("gr_j2_p4_phi", &gr_j2_p4_phi);
	intree->SetBranchAddress("gr_j2_p4_mass", &gr_j2_p4_mass);

	intree->SetBranchAddress("j1_e", &j1_e);
	intree->SetBranchAddress("j1_pt", &j1_pt);
	intree->SetBranchAddress("j1_phi", &j1_phi);
	intree->SetBranchAddress("j1_eta", &j1_eta);
	intree->SetBranchAddress("j1_beta", &j1_beta);
	intree->SetBranchAddress("j1_betaStar", &j1_betaStar);
	intree->SetBranchAddress("j1_betaStarClassic", &j1_betaStarClassic);
	intree->SetBranchAddress("j1_dR2Mean", &j1_dR2Mean);
	intree->SetBranchAddress("j1_csvBtag", &j1_csvBtag);
	intree->SetBranchAddress("j1_csvMvaBtag", &j1_csvMvaBtag);
	intree->SetBranchAddress("j1_jetProbBtag", &j1_jetProbBtag);
	intree->SetBranchAddress("j1_tcheBtag", &j1_tcheBtag);
	intree->SetBranchAddress("j1_radionMatched", &j1_radionMatched);
	intree->SetBranchAddress("j1_ptD", &j1_ptD);
	intree->SetBranchAddress("j1_nSecondaryVertices", &j1_nSecondaryVertices);
	intree->SetBranchAddress("j1_secVtxPt", &j1_secVtxPt);
	intree->SetBranchAddress("j1_secVtx3dL", &j1_secVtx3dL);
	intree->SetBranchAddress("j1_secVtx3deL", &j1_secVtx3deL);
	intree->SetBranchAddress("j1_emfrac", &j1_emfrac);
	intree->SetBranchAddress("j1_hadfrac", &j1_hadfrac);
	intree->SetBranchAddress("j1_ntk", &j1_ntk);
	intree->SetBranchAddress("j1_nNeutrals", &j1_nNeutrals);
	intree->SetBranchAddress("j1_nCharged", &j1_nCharged);
	intree->SetBranchAddress("j1_axis1", &j1_axis1);
	intree->SetBranchAddress("j1_axis2", &j1_axis2);
	intree->SetBranchAddress("j1_pull", &j1_pull);
	intree->SetBranchAddress("j1_Rchg", &j1_Rchg);
	intree->SetBranchAddress("j1_Rneutral", &j1_Rneutral);
	intree->SetBranchAddress("j1_R", &j1_R);

	intree->SetBranchAddress("j2_e", &j2_e);
	intree->SetBranchAddress("j2_pt", &j2_pt);
	intree->SetBranchAddress("j2_phi", &j2_phi);
	intree->SetBranchAddress("j2_eta", &j2_eta);
	intree->SetBranchAddress("j2_beta", &j2_beta);
	intree->SetBranchAddress("j2_betaStar", &j2_betaStar);
	intree->SetBranchAddress("j2_betaStarClassic", &j2_betaStarClassic);
	intree->SetBranchAddress("j2_dR2Mean", &j2_dR2Mean);
	intree->SetBranchAddress("j2_csvBtag", &j2_csvBtag);
	intree->SetBranchAddress("j2_csvMvaBtag", &j2_csvMvaBtag);
	intree->SetBranchAddress("j2_jetProbBtag", &j2_jetProbBtag);
	intree->SetBranchAddress("j2_tcheBtag", &j2_tcheBtag);
	intree->SetBranchAddress("j2_radionMatched", &j2_radionMatched);
	intree->SetBranchAddress("j2_ptD", &j2_ptD);
	intree->SetBranchAddress("j2_nSecondaryVertices", &j2_nSecondaryVertices);
	intree->SetBranchAddress("j2_secVtxPt", &j2_secVtxPt);
	intree->SetBranchAddress("j2_secVtx3dL", &j2_secVtx3dL);
	intree->SetBranchAddress("j2_secVtx3deL", &j2_secVtx3deL);
	intree->SetBranchAddress("j2_emfrac", &j2_emfrac);
	intree->SetBranchAddress("j2_hadfrac", &j2_hadfrac);
	intree->SetBranchAddress("j2_ntk", &j2_ntk);
	intree->SetBranchAddress("j2_nNeutrals", &j2_nNeutrals);
	intree->SetBranchAddress("j2_nCharged", &j2_nCharged);
	intree->SetBranchAddress("j2_axis1", &j2_axis1);
	intree->SetBranchAddress("j2_axis2", &j2_axis2);
	intree->SetBranchAddress("j2_pull", &j2_pull);
	intree->SetBranchAddress("j2_Rchg", &j2_Rchg);
	intree->SetBranchAddress("j2_Rneutral", &j2_Rneutral);
	intree->SetBranchAddress("j2_R", &j2_R);

	intree->SetBranchAddress("j3_e", &j3_e);
	intree->SetBranchAddress("j3_pt", &j3_pt);
	intree->SetBranchAddress("j3_phi", &j3_phi);
	intree->SetBranchAddress("j3_eta", &j3_eta);
	intree->SetBranchAddress("j3_beta", &j3_beta);
	intree->SetBranchAddress("j3_betaStar", &j3_betaStar);
	intree->SetBranchAddress("j3_betaStarClassic", &j3_betaStarClassic);
	intree->SetBranchAddress("j3_dR2Mean", &j3_dR2Mean);
	intree->SetBranchAddress("j3_csvBtag", &j3_csvBtag);
	intree->SetBranchAddress("j3_csvMvaBtag", &j3_csvMvaBtag);
	intree->SetBranchAddress("j3_jetProbBtag", &j3_jetProbBtag);
	intree->SetBranchAddress("j3_tcheBtag", &j3_tcheBtag);
	intree->SetBranchAddress("j3_radionMatched", &j3_radionMatched);
	intree->SetBranchAddress("j3_ptD", &j3_ptD);
	intree->SetBranchAddress("j3_nSecondaryVertices", &j3_nSecondaryVertices);
	intree->SetBranchAddress("j3_secVtxPt", &j3_secVtxPt);
	intree->SetBranchAddress("j3_secVtx3dL", &j3_secVtx3dL);
	intree->SetBranchAddress("j3_secVtx3deL", &j3_secVtx3deL);
	intree->SetBranchAddress("j3_emfrac", &j3_emfrac);
	intree->SetBranchAddress("j3_hadfrac", &j3_hadfrac);
	intree->SetBranchAddress("j3_ntk", &j3_ntk);
	intree->SetBranchAddress("j3_nNeutrals", &j3_nNeutrals);
	intree->SetBranchAddress("j3_nCharged", &j3_nCharged);
	intree->SetBranchAddress("j3_axis1", &j3_axis1);
	intree->SetBranchAddress("j3_axis2", &j3_axis2);
	intree->SetBranchAddress("j3_pull", &j3_pull);
	intree->SetBranchAddress("j3_Rchg", &j3_Rchg);
	intree->SetBranchAddress("j3_Rneutral", &j3_Rneutral);
	intree->SetBranchAddress("j3_R", &j3_R);

	intree->SetBranchAddress("j4_e", &j4_e);
	intree->SetBranchAddress("j4_pt", &j4_pt);
	intree->SetBranchAddress("j4_phi", &j4_phi);
	intree->SetBranchAddress("j4_eta", &j4_eta);
	intree->SetBranchAddress("j4_beta", &j4_beta);
	intree->SetBranchAddress("j4_betaStar", &j4_betaStar);
	intree->SetBranchAddress("j4_betaStarClassic", &j4_betaStarClassic);
	intree->SetBranchAddress("j4_dR2Mean", &j4_dR2Mean);
	intree->SetBranchAddress("j4_csvBtag", &j4_csvBtag);
	intree->SetBranchAddress("j4_csvMvaBtag", &j4_csvMvaBtag);
	intree->SetBranchAddress("j4_jetProbBtag", &j4_jetProbBtag);
	intree->SetBranchAddress("j4_tcheBtag", &j4_tcheBtag);
	intree->SetBranchAddress("j4_radionMatched", &j4_radionMatched);
	intree->SetBranchAddress("j4_ptD", &j4_ptD);
	intree->SetBranchAddress("j4_nSecondaryVertices", &j4_nSecondaryVertices);
	intree->SetBranchAddress("j4_secVtxPt", &j4_secVtxPt);
	intree->SetBranchAddress("j4_secVtx3dL", &j4_secVtx3dL);
	intree->SetBranchAddress("j4_secVtx3deL", &j4_secVtx3deL);
	intree->SetBranchAddress("j4_emfrac", &j4_emfrac);
	intree->SetBranchAddress("j4_hadfrac", &j4_hadfrac);
	intree->SetBranchAddress("j4_ntk", &j4_ntk);
	intree->SetBranchAddress("j4_nNeutrals", &j4_nNeutrals);
	intree->SetBranchAddress("j4_nCharged", &j4_nCharged);
	intree->SetBranchAddress("j4_axis1", &j4_axis1);
	intree->SetBranchAddress("j4_axis2", &j4_axis2);
	intree->SetBranchAddress("j4_pull", &j4_pull);
	intree->SetBranchAddress("j4_Rchg", &j4_Rchg);
	intree->SetBranchAddress("j4_Rneutral", &j4_Rneutral);
	intree->SetBranchAddress("j4_R", &j4_R);

	outtree->Branch("jet_index1", &jet_index1, "jet_index1/i");
	outtree->Branch("jet_index2", &jet_index2, "jet_index2/i");
	outtree->Branch("jet_e", &jet_e, "jet_e/F");
	outtree->Branch("jet_pt", &jet_pt, "jet_pt/F");
	outtree->Branch("jet_phi", &jet_phi, "jet_phi/F");
	outtree->Branch("jet_eta", &jet_eta, "jet_eta/F");
	outtree->Branch("jet_beta", &jet_beta, "jet_beta/F");
	outtree->Branch("jet_betaStar", &jet_betaStar, "jet_betaStar/F");
	outtree->Branch("jet_betaStarClassic", &jet_betaStarClassic, "jet_betaStarClassic/F");
	outtree->Branch("jet_dR2Mean", &jet_dR2Mean, "jet_dR2Mean/F");
	outtree->Branch("jet_csvBtag", &jet_csvBtag, "jet_csvBtag/F");
	outtree->Branch("jet_csvMvaBtag", &jet_csvMvaBtag, "jet_csvMvaBtag/F");
	outtree->Branch("jet_jetProbBtag", &jet_jetProbBtag, "jet_jetProbBtag/F");
	outtree->Branch("jet_tcheBtag", &jet_tcheBtag, "jet_tcheBtag/F");
	outtree->Branch("jet_radionMatched", &jet_radionMatched, "jet_radionMatched/F");
	outtree->Branch("jet_ptD", &jet_ptD, "jet_ptD/F");
	outtree->Branch("jet_nSecondaryVertices", &jet_nSecondaryVertices, "jet_nSecondaryVertices/F");
	outtree->Branch("jet_secVtxPt", &jet_secVtxPt, "jet_secVtxPt/F");
	outtree->Branch("jet_secVtx3dL", &jet_secVtx3dL, "jet_secVtx3dL/F");
	outtree->Branch("jet_secVtx3deL", &jet_secVtx3deL, "jet_secVtx3deL/F");
	outtree->Branch("jet_emfrac", &jet_emfrac, "jet_emfrac/F");
	outtree->Branch("jet_hadfrac", &jet_hadfrac, "jet_hadfrac/F");
	outtree->Branch("jet_ntk", &jet_ntk, "jet_ntk/F");
	outtree->Branch("jet_nNeutrals", &jet_nNeutrals, "jet_nNeutrals/F");
	outtree->Branch("jet_nCharged", &jet_nCharged, "jet_nCharged/F");
	outtree->Branch("jet_axis1", &jet_axis1, "jet_axis1/F");
	outtree->Branch("jet_axis2", &jet_axis2, "jet_axis2/F");
	outtree->Branch("jet_pull", &jet_pull, "jet_pull/F");
	outtree->Branch("jet_Rchg", &jet_Rchg, "jet_Rchg/F");
	outtree->Branch("jet_Rneutral", &jet_Rneutral, "jet_Rneutral/F");
	outtree->Branch("jet_R", &jet_R, "jet_R/F");
	outtree->Branch("jet_DR_b1", &jet_DR_b1, "jet_DR_b1/F");
	outtree->Branch("jet_DR_b2", &jet_DR_b2, "jet_DR_b2/F");
	outtree->Branch("jet_DR_j1", &jet_DR_j1, "jet_DR_j1/F");
	outtree->Branch("jet_DR_j2", &jet_DR_j2, "jet_DR_j2/F");

	int npass = 0;
	int np[20] = {0};
	for(int ievt=0 ; ievt < intree->GetEntries() ; ievt++)
//	for(int ievt=0 ; ievt < 50 ; ievt++)
	{
		intree->GetEntry(ievt);
//		cout << "Matching: " << j1_radionMatched << "\t" << j2_radionMatched << "\t" << j3_radionMatched << "\t" << j4_radionMatched << endl;
		cout << gr_b1_p4_pt << "\t" << gr_b2_p4_pt << "\t" << gr_j1_p4_pt << "\t" << gr_j2_p4_pt << endl;
	
		// take only the subset of events where two genjets are matched to the signal
		TLorentzVector gen_b1, gen_b2, gen_j1, gen_j2;
		if( gr_b1_p4_pt < 1. ) continue;
		np[0]++;
		if( gr_b2_p4_pt < 1. ) continue;
		np[1]++;
		if( gr_j1_p4_pt < 1. ) continue;
		np[2]++;
		if( gr_j2_p4_pt < 1. ) continue;
		np[3]++;

		npass++;
/*
		cout << "gen_b1(pt, eta, phi, mass)= " << "(\t" << gr_b1_p4_pt << ",\t" << gr_b1_p4_eta << ",\t" << gr_b1_p4_phi << ",\t" << gr_b1_p4_mass << "\t)" << endl;;
		gen_b1.SetPtEtaPhiM(gr_b1_p4_pt, gr_b1_p4_eta, gr_b1_p4_phi, gr_b1_p4_mass);
		cout << "gen_b2" << endl;
		gen_b2.SetPtEtaPhiM(gr_b2_p4_pt, gr_b2_p4_eta, gr_b2_p4_phi, gr_b2_p4_mass);
		cout << "gen_j1" << endl;
		gen_j1.SetPtEtaPhiM(gr_j1_p4_pt, gr_j1_p4_eta, gr_j1_p4_phi, gr_j1_p4_mass);
		cout << "gen_j2" << endl;
		gen_j2.SetPtEtaPhiM(gr_j2_p4_pt, gr_j2_p4_eta, gr_j2_p4_phi, gr_j2_p4_mass);


		vector<float> DR_jet_b1;
		vector<float> DR_jet_b2;
		vector<float> DR_jet_j1;
		vector<float> DR_jet_j2;
		for(int ijet = 0 ; ijet < 4 ; ijet++ )
		{
			TLorentzVector jet;
			if(ijet == 0 && j1_e > 0.)
				jet.SetPtEtaPhiE(j1_pt, j1_eta, j1_phi, j1_e);
			if(ijet == 1 && j2_e > 0.)
				jet.SetPtEtaPhiE(j2_pt, j2_eta, j2_phi, j2_e);
			if(ijet == 2 && j3_e > 0.)
				jet.SetPtEtaPhiE(j3_pt, j3_eta, j3_phi, j3_e);
			if(ijet == 3 && j4_e > 0.)
				jet.SetPtEtaPhiE(j4_pt, j4_eta, j4_phi, j4_e);
			if(jet.Pt() < 1.0) continue;
			DR_jet_b1.push_back(jet.DeltaR(gen_b1));
			DR_jet_b2.push_back(jet.DeltaR(gen_b2));
			DR_jet_j1.push_back(jet.DeltaR(gen_j1));
			DR_jet_j2.push_back(jet.DeltaR(gen_j2));
		}
*/
		jet_e = j1_e;
		outtree->Fill();


/*
		// at least one match
		if( !(j1_radionMatched || j2_radionMatched || j3_radionMatched || j4_radionMatched)) continue;
		// at least a pair of matched
		if( !(
				(j1_radionMatched && j2_radionMatched)
				|| (j1_radionMatched && j3_radionMatched)
				|| (j1_radionMatched && j4_radionMatched)
				|| (j2_radionMatched && j3_radionMatched)
				|| (j2_radionMatched && j4_radionMatched)
				|| (j3_radionMatched && j4_radionMatched)
			) ) continue;

		int i1, i2;
		i1 = i2 = 0;
		if( j1_radionMatched )
			i1 = 1;
		if( j2_radionMatched )
		{
			if( i1 == 0 )
				i1 = 2;
			else
				i2 = 2;
		}
		if( j3_radionMatched && ((i1==0) || (i2==0)))
		{
			if( i1 == 0 )
				i1 = 3;
			else
				i2 = 3;
		}
		if( j4_radionMatched && ((i1==0) || (i2==0)))
		{
			if( i1 == 0 )
				i1 = 4;
			else
				i2 = 4;
		}

		jet_index1 = i1;
		jet_index1 = i2;
		outtree->Fill();
*/
	// treat only two jets per event: here is where the choice of the jets is made
// 

/*
		jet_index = 1;
		jet_e = j1_e;
		jet_pt = j1_pt;
		jet_phi = j1_phi;
		jet_eta = j1_eta;
		jet_beta = j1_beta;
		jet_betaStar = j1_betaStar;
		jet_betaStarClassic = j1_betaStarClassic;
		jet_dR2Mean = j1_dR2Mean;
		jet_csvBtag = j1_csvBtag;
		jet_csvMvaBtag = j1_csvMvaBtag;
		jet_jetProbBtag = j1_jetProbBtag;
		jet_tcheBtag = j1_tcheBtag;
		jet_radionMatched = j1_radionMatched;
		jet_ptD = j1_ptD;
		jet_nSecondaryVertices = j1_nSecondaryVertices;
		jet_secVtxPt = j1_secVtxPt;
		jet_secVtx3dL = j1_secVtx3dL;
		jet_secVtx3deL = j1_secVtx3deL;
		jet_emfrac = j1_emfrac;
		jet_hadfrac = j1_hadfrac;
		jet_ntk = j1_ntk;
		jet_nNeutrals = j1_nNeutrals;
		jet_nCharged = j1_nCharged;
		jet_axis1 = j1_axis1;
		jet_axis2 = j1_axis2;
		jet_pull = j1_pull;
		jet_Rchg = j1_Rchg;
		jet_Rneutral = j1_Rneutral;
		jet_R = j1_R;
		outtree->Fill();
		jet_index = 2;
		jet_e = j2_e;
		jet_pt = j2_pt;
		jet_phi = j2_phi;
		jet_eta = j2_eta;
		jet_beta = j2_beta;
		jet_betaStar = j2_betaStar;
		jet_betaStarClassic = j2_betaStarClassic;
		jet_dR2Mean = j2_dR2Mean;
		jet_csvBtag = j2_csvBtag;
		jet_csvMvaBtag = j2_csvMvaBtag;
		jet_jetProbBtag = j2_jetProbBtag;
		jet_tcheBtag = j2_tcheBtag;
		jet_radionMatched = j2_radionMatched;
		jet_ptD = j2_ptD;
		jet_nSecondaryVertices = j2_nSecondaryVertices;
		jet_secVtxPt = j2_secVtxPt;
		jet_secVtx3dL = j2_secVtx3dL;
		jet_secVtx3deL = j2_secVtx3deL;
		jet_emfrac = j2_emfrac;
		jet_hadfrac = j2_hadfrac;
		jet_ntk = j2_ntk;
		jet_nNeutrals = j2_nNeutrals;
		jet_nCharged = j2_nCharged;
		jet_axis1 = j2_axis1;
		jet_axis2 = j2_axis2;
		jet_pull = j2_pull;
		jet_Rchg = j2_Rchg;
		jet_Rneutral = j2_Rneutral;
		jet_R = j2_R;
		outtree->Fill();
		jet_index = 3;
		jet_e = j3_e;
		jet_pt = j3_pt;
		jet_phi = j3_phi;
		jet_eta = j3_eta;
		jet_beta = j3_beta;
		jet_betaStar = j3_betaStar;
		jet_betaStarClassic = j3_betaStarClassic;
		jet_dR2Mean = j3_dR2Mean;
		jet_csvBtag = j3_csvBtag;
		jet_csvMvaBtag = j3_csvMvaBtag;
		jet_jetProbBtag = j3_jetProbBtag;
		jet_tcheBtag = j3_tcheBtag;
		jet_radionMatched = j3_radionMatched;
		jet_ptD = j3_ptD;
		jet_nSecondaryVertices = j3_nSecondaryVertices;
		jet_secVtxPt = j3_secVtxPt;
		jet_secVtx3dL = j3_secVtx3dL;
		jet_secVtx3deL = j3_secVtx3deL;
		jet_emfrac = j3_emfrac;
		jet_hadfrac = j3_hadfrac;
		jet_ntk = j3_ntk;
		jet_nNeutrals = j3_nNeutrals;
		jet_nCharged = j3_nCharged;
		jet_axis1 = j3_axis1;
		jet_axis2 = j3_axis2;
		jet_pull = j3_pull;
		jet_Rchg = j3_Rchg;
		jet_Rneutral = j3_Rneutral;
		jet_R = j3_R;
		outtree->Fill();
		jet_index = 4;
		jet_e = j4_e;
		jet_pt = j4_pt;
		jet_phi = j4_phi;
		jet_eta = j4_eta;
		jet_beta = j4_beta;
		jet_betaStar = j4_betaStar;
		jet_betaStarClassic = j4_betaStarClassic;
		jet_dR2Mean = j4_dR2Mean;
		jet_csvBtag = j4_csvBtag;
		jet_csvMvaBtag = j4_csvMvaBtag;
		jet_jetProbBtag = j4_jetProbBtag;
		jet_tcheBtag = j4_tcheBtag;
		jet_radionMatched = j4_radionMatched;
		jet_ptD = j4_ptD;
		jet_nSecondaryVertices = j4_nSecondaryVertices;
		jet_secVtxPt = j4_secVtxPt;
		jet_secVtx3dL = j4_secVtx3dL;
		jet_secVtx3deL = j4_secVtx3deL;
		jet_emfrac = j4_emfrac;
		jet_hadfrac = j4_hadfrac;
		jet_ntk = j4_ntk;
		jet_nNeutrals = j4_nNeutrals;
		jet_nCharged = j4_nCharged;
		jet_axis1 = j4_axis1;
		jet_axis2 = j4_axis2;
		jet_pull = j4_pull;
		jet_Rchg = j4_Rchg;
		jet_Rneutral = j4_Rneutral;
		jet_R = j4_R;
		outtree->Fill();
*/
	}
	for(int i=0 ; i < 4 ; i++)
		cout << "#np[" << i << "]= " << np[i] << endl;
	cout << "#npass= " << npass << endl;


	outfile->cd();
	outtree->Write();
	outfile->Close();
	infile->Close();

	return 0;
}

