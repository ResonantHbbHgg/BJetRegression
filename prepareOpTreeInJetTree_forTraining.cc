// BJet regression: tree preparation for training from h2gglobe opTrees
// Original code by K. Kousouris
// O. Bondu (May 2013)
// C++ headers
#include <iostream>
#include <string>
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

int main(int argc, char *argv[])
{
	if(DEBUG) cout << "DEBUG: Initialisation: reading parameters" << endl;
	cout << "argc= " << argc << endl;
	for(int iarg = 0 ; iarg < argc; iarg++)
		cout << "argv[" << iarg << "]= " << argv[iarg] << endl;

	if( argc == 1 )
	{
		cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
		cerr << "WARNING: Syntax is " << argv[0] << " -i (inputfile) -it (inputtree) -o (outputfile) -ot (outputtree)" << endl;
	}

	string inputfile = "root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_08_tree_05/Graviton_Radion-nm.root";
	string inputtree = "Radion_m300_8TeV_nm";
	string outputfile = "jetTreeForTraining.root";
	string outputtree = "jets";

	for(int iarg=0 ; iarg < argc ; iarg++)
	{
		if(strcmp("-i", argv[iarg]) == 0 && argc >= iarg + 1)
			inputfile = argv[iarg+1];
		if(strcmp("-it", argv[iarg]) == 0 && argc >= iarg + 1)
			inputtree = argv[iarg+1];
		if(strcmp("-o", argv[iarg]) == 0 && argc >= iarg + 1)
			outputfile = argv[iarg+1];
		if(strcmp("-ot", argv[iarg]) == 0 && argc >= iarg + 1)
			outputtree = argv[iarg+1];
		if((strcmp("-h", argv[iarg]) == 0) || (strcmp("--help", argv[iarg]) == 0))
		{
			cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
			cerr << "WARNING: Syntax is " << argv[0] << " -i (inputfile) -it (inputtree) -o (outputfile) -ot (outputtree)" << endl;
			cerr << "inputfile= " << inputfile << endl;
			cerr << "inputtree= " << inputtree << endl;
			cerr << "outputfile= " << outputfile << endl;
			cerr << "outputtree= " << outputtree << endl;
			return 2;
		}
	}

	cout << "inputfile= " << inputfile << endl;
	cout << "inputtree= " << inputtree << endl;
	cout << "outputfile= " << outputfile << endl;
	cout << "outputtree= " << outputtree << endl;

	TFile *infile = TFile::Open(inputfile.c_str());
	TTree *intree = (TTree*)infile->Get(inputtree.c_str());
	TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
	TTree *outtree = new TTree(outputtree.c_str(), Form("%s for training", outputtree.c_str()));

	if(DEBUG) cout << "DEBUG: creating all the needed floats" << endl;

	float event;
	float ph1_pt, ph2_pt;
	float gr_radion_p4_pt, gr_radion_p4_eta, gr_radion_p4_phi, gr_radion_p4_mass, gr_hgg_p4_pt, gr_hgg_p4_eta, gr_hgg_p4_phi, gr_hgg_p4_mass, gr_hbb_p4_pt, gr_hbb_p4_eta, gr_hbb_p4_phi, gr_hbb_p4_mass, gr_g1_p4_pt, gr_g1_p4_eta, gr_g1_p4_phi, gr_g1_p4_mass, gr_g2_p4_pt, gr_g2_p4_eta, gr_g2_p4_phi, gr_g2_p4_mass, gr_b1_p4_pt, gr_b1_p4_eta, gr_b1_p4_phi, gr_b1_p4_mass, gr_b2_p4_pt, gr_b2_p4_eta, gr_b2_p4_phi, gr_b2_p4_mass, gr_j1_p4_pt, gr_j1_p4_eta, gr_j1_p4_phi, gr_j1_p4_mass, gr_j2_p4_pt, gr_j2_p4_eta, gr_j2_p4_phi, gr_j2_p4_mass;
	float met_corr_pfmet, met_corr_phi_pfmet, met_corr_eta_pfmet, met_corr_e_pfmet;
//	float met_pfmet, met_phi_pfmet, met_sumet_pfmet, met_mEtSig_pfmet, met_significance_pfmet, met_corrMet, met_corrMetPhi;
	float pu_n, nvtx, rho;
	float weight, evweight, pu_weight;
	float j1_e, j1_pt, j1_phi, j1_eta, j1_beta, j1_betaStar, j1_betaStarClassic, j1_dR2Mean, j1_csvBtag, j1_csvMvaBtag, j1_jetProbBtag, j1_tcheBtag, j1_ptD, j1_nSecondaryVertices, j1_secVtxPt, j1_secVtx3dL, j1_secVtx3deL, j1_emfrac, j1_hadfrac, j1_axis1, j1_axis2, j1_pull, /*j1_Rchg, j1_Rneutral, j1_R, j1_chargedMultiplicity, j1_neutralMultiplicity, j1_Chadfrac, j1_Nhadfrac, j1_Phofrac, j1_Mufrac, j1_Elefrac, j1_dPhiMet,*/ j1_radionMatched;
	int j1_nNeutrals, j1_nCharged, j1_ntk/*, j1_pfloose*/;
	float j2_e, j2_pt, j2_phi, j2_eta, j2_beta, j2_betaStar, j2_betaStarClassic, j2_dR2Mean, j2_csvBtag, j2_csvMvaBtag, j2_jetProbBtag, j2_tcheBtag, j2_ptD, j2_nSecondaryVertices, j2_secVtxPt, j2_secVtx3dL, j2_secVtx3deL, j2_emfrac, j2_hadfrac, j2_axis1, j2_axis2, j2_pull, /*j2_Rchg, j2_Rneutral, j2_R, j2_chargedMultiplicity, j2_neutralMultiplicity, j2_Chadfrac, j2_Nhadfrac, j2_Phofrac, j2_Mufrac, j2_Elefrac, j2_dPhiMet,*/ j2_radionMatched;
	int j2_nNeutrals, j2_nCharged, j2_ntk/*, j2_pfloose*/;
	float j3_e, j3_pt, j3_phi, j3_eta, j3_beta, j3_betaStar, j3_betaStarClassic, j3_dR2Mean, j3_csvBtag, j3_csvMvaBtag, j3_jetProbBtag, j3_tcheBtag, j3_ptD, j3_nSecondaryVertices, j3_secVtxPt, j3_secVtx3dL, j3_secVtx3deL, j3_emfrac, j3_hadfrac, j3_axis1, j3_axis2, j3_pull, /*j3_Rchg, j3_Rneutral, j3_R, j3_chargedMultiplicity, j3_neutralMultiplicity, j3_Chadfrac, j3_Nhadfrac, j3_Phofrac, j3_Mufrac, j3_Elefrac, j3_dPhiMet,*/ j3_radionMatched;
	int j3_nNeutrals, j3_nCharged, j3_ntk/*, j3_pfloose*/;
	float j4_e, j4_pt, j4_phi, j4_eta, j4_beta, j4_betaStar, j4_betaStarClassic, j4_dR2Mean, j4_csvBtag, j4_csvMvaBtag, j4_jetProbBtag, j4_tcheBtag, j4_ptD, j4_nSecondaryVertices, j4_secVtxPt, j4_secVtx3dL, j4_secVtx3deL, j4_emfrac, j4_hadfrac, j4_axis1, j4_axis2, j4_pull, /*j4_Rchg, j4_Rneutral, j4_R, j4_chargedMultiplicity, j4_neutralMultiplicity, j4_Chadfrac, j4_Nhadfrac, j4_Phofrac, j4_Mufrac, j4_Elefrac, j4_dPhiMet,*/ j4_radionMatched;
	int j4_nNeutrals, j4_nCharged, j4_ntk/*, j4_pfloose*/;
	float jet_e, jet_pt, jet_phi, jet_eta, jet_beta, jet_betaStar, jet_betaStarClassic, jet_dR2Mean, jet_csvBtag, jet_csvMvaBtag, jet_jetProbBtag, jet_tcheBtag, jet_ptD, jet_nSecondaryVertices, jet_secVtxPt, jet_secVtx3dL, jet_secVtx3deL, jet_emfrac, jet_hadfrac, jet_axis1, jet_axis2, jet_pull, /*jet_Rchg, jet_Rneutral, jet_R, jet_chargedMultiplicity, jet_neutralMultiplicity, jet_Chadfrac, jet_Nhadfrac, jet_Phofrac, jet_Mufrac, jet_Elefrac,*/ jet_dPhiMet, jet_radionMatched, jet_nConstituents;
	int jet_nNeutrals, jet_nCharged, jet_ntk/*, jet_pfloose*/;
	float jet_genDR, jet_genPt, jet_genE, jet_genR, jet_prtDR, jet_prtPt, jet_prtE, jet_prtR;
	int jet_index;
//	float ev_met_pfmet, ev_met_phi_pfmet, ev_met_sumet_pfmet, ev_met_mEtSig_pfmet, ev_met_significance_pfmet, ev_met_corrMet, ev_met_corrMetPhi;
	float ev_met_corr_pfmet, ev_met_corr_phi_pfmet, ev_met_corr_eta_pfmet, ev_met_corr_e_pfmet;
	float ev_pu_n, ev_nvtx, ev_rho;
	float ev_weight, ev_evweight, ev_pu_weight;

	int njets_passing_kLooseID;
	intree->SetBranchAddress("njets_passing_kLooseID", &njets_passing_kLooseID);
	int njets_passing_kLooseID_and_CSVM;
	intree->SetBranchAddress("njets_passing_kLooseID_and_CSVM", &njets_passing_kLooseID_and_CSVM);

/*	intree->SetBranchAddress("met_pfmet", &met_pfmet);
	intree->SetBranchAddress("met_phi_pfmet", &met_phi_pfmet);
	intree->SetBranchAddress("met_sumet_pfmet", &met_sumet_pfmet);
	intree->SetBranchAddress("met_mEtSig_pfmet", &met_mEtSig_pfmet);
	intree->SetBranchAddress("met_significance_pfmet", &met_significance_pfmet);
	intree->SetBranchAddress("met_corrMet", &met_corrMet);
	intree->SetBranchAddress("met_corrMetPhi", &met_corrMetPhi);
*/
	intree->SetBranchAddress("event", &event);
	intree->SetBranchAddress("ph1_pt", &ph1_pt);
	intree->SetBranchAddress("ph2_pt", &ph2_pt);
	intree->SetBranchAddress("met_corr_pfmet", &met_corr_pfmet);
	intree->SetBranchAddress("met_corr_phi_pfmet", &met_corr_phi_pfmet);
	intree->SetBranchAddress("met_corr_eta_pfmet", &met_corr_eta_pfmet);
	intree->SetBranchAddress("met_corr_e_pfmet", &met_corr_e_pfmet);
	intree->SetBranchAddress("pu_n", &pu_n);
	intree->SetBranchAddress("nvtx", &nvtx);
	intree->SetBranchAddress("rho", &rho);
	intree->SetBranchAddress("weight", &weight);
	intree->SetBranchAddress("evweight", &evweight);
	intree->SetBranchAddress("pu_weight", &pu_weight);
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
/*	intree->SetBranchAddress("j1_Rchg", &j1_Rchg);
	intree->SetBranchAddress("j1_Rneutral", &j1_Rneutral);
	intree->SetBranchAddress("j1_R", &j1_R);
	intree->SetBranchAddress("j1_chargedMultiplicity", &j1_chargedMultiplicity);
	intree->SetBranchAddress("j1_neutralMultiplicity", &j1_neutralMultiplicity);
	intree->SetBranchAddress("j1_Chadfrac", &j1_Chadfrac);
	intree->SetBranchAddress("j1_Nhadfrac", &j1_Nhadfrac);
	intree->SetBranchAddress("j1_Phofrac", &j1_Phofrac);
	intree->SetBranchAddress("j1_Mufrac", &j1_Mufrac);
	intree->SetBranchAddress("j1_Elefrac", &j1_Elefrac);
	intree->SetBranchAddress("j1_dPhiMet", &j1_dPhiMet);
	intree->SetBranchAddress("j1_pfloose", &j1_pfloose);
*/
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
/*	intree->SetBranchAddress("j2_Rchg", &j2_Rchg);
	intree->SetBranchAddress("j2_Rneutral", &j2_Rneutral);
	intree->SetBranchAddress("j2_R", &j2_R);
	intree->SetBranchAddress("j2_chargedMultiplicity", &j2_chargedMultiplicity);
	intree->SetBranchAddress("j2_neutralMultiplicity", &j2_neutralMultiplicity);
	intree->SetBranchAddress("j2_Chadfrac", &j2_Chadfrac);
	intree->SetBranchAddress("j2_Nhadfrac", &j2_Nhadfrac);
	intree->SetBranchAddress("j2_Phofrac", &j2_Phofrac);
	intree->SetBranchAddress("j2_Mufrac", &j2_Mufrac);
	intree->SetBranchAddress("j2_Elefrac", &j2_Elefrac);
	intree->SetBranchAddress("j2_dPhiMet", &j2_dPhiMet);
	intree->SetBranchAddress("j2_pfloose", &j2_pfloose);
*/
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
/*	intree->SetBranchAddress("j3_Rchg", &j3_Rchg);
	intree->SetBranchAddress("j3_Rneutral", &j3_Rneutral);
	intree->SetBranchAddress("j3_R", &j3_R);
	intree->SetBranchAddress("j3_chargedMultiplicity", &j3_chargedMultiplicity);
	intree->SetBranchAddress("j3_neutralMultiplicity", &j3_neutralMultiplicity);
	intree->SetBranchAddress("j3_Chadfrac", &j3_Chadfrac);
	intree->SetBranchAddress("j3_Nhadfrac", &j3_Nhadfrac);
	intree->SetBranchAddress("j3_Phofrac", &j3_Phofrac);
	intree->SetBranchAddress("j3_Mufrac", &j3_Mufrac);
	intree->SetBranchAddress("j3_Elefrac", &j3_Elefrac);
	intree->SetBranchAddress("j3_dPhiMet", &j3_dPhiMet);
	intree->SetBranchAddress("j3_pfloose", &j3_pfloose);
*/
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
/*	intree->SetBranchAddress("j4_Rchg", &j4_Rchg);
	intree->SetBranchAddress("j4_Rneutral", &j4_Rneutral);
	intree->SetBranchAddress("j4_R", &j4_R);
	intree->SetBranchAddress("j4_chargedMultiplicity", &j4_chargedMultiplicity);
	intree->SetBranchAddress("j4_neutralMultiplicity", &j4_neutralMultiplicity);
	intree->SetBranchAddress("j4_Chadfrac", &j4_Chadfrac);
	intree->SetBranchAddress("j4_Nhadfrac", &j4_Nhadfrac);
	intree->SetBranchAddress("j4_Phofrac", &j4_Phofrac);
	intree->SetBranchAddress("j4_Mufrac", &j4_Mufrac);
	intree->SetBranchAddress("j4_Elefrac", &j4_Elefrac);
	intree->SetBranchAddress("j4_dPhiMet", &j4_dPhiMet);
	intree->SetBranchAddress("j4_pfloose", &j4_pfloose);
*/
	outtree->Branch("event", &event, "event/F");
	outtree->Branch("ph1_pt", &ph1_pt, "ph1_pt/F");
	outtree->Branch("ph2_pt", &ph2_pt, "ph2_pt/F");
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
	outtree->Branch("jet_nConstituents", &jet_nConstituents, "jet_nConstituents/F");
	outtree->Branch("jet_ptD", &jet_ptD, "jet_ptD/F");
	outtree->Branch("jet_nSecondaryVertices", &jet_nSecondaryVertices, "jet_nSecondaryVertices/F");
	outtree->Branch("jet_secVtxPt", &jet_secVtxPt, "jet_secVtxPt/F");
	outtree->Branch("jet_secVtx3dL", &jet_secVtx3dL, "jet_secVtx3dL/F");
	outtree->Branch("jet_secVtx3deL", &jet_secVtx3deL, "jet_secVtx3deL/F");
	outtree->Branch("jet_emfrac", &jet_emfrac, "jet_emfrac/F");
	outtree->Branch("jet_hadfrac", &jet_hadfrac, "jet_hadfrac/F");
	outtree->Branch("jet_ntk", &jet_ntk, "jet_ntk/I");
	outtree->Branch("jet_nNeutrals", &jet_nNeutrals, "jet_nNeutrals/I");
	outtree->Branch("jet_nCharged", &jet_nCharged, "jet_nCharged/I");
	outtree->Branch("jet_axis1", &jet_axis1, "jet_axis1/F");
	outtree->Branch("jet_axis2", &jet_axis2, "jet_axis2/F");
	outtree->Branch("jet_pull", &jet_pull, "jet_pull/F");
/*	outtree->Branch("jet_Rchg", &jet_Rchg, "jet_Rchg/F");
	outtree->Branch("jet_Rneutral", &jet_Rneutral, "jet_Rneutral/F");
	outtree->Branch("jet_R", &jet_R, "jet_R/F");
	outtree->Branch("jet_chargedMultiplicity", &jet_chargedMultiplicity, "jet_chargedMultiplicity/F");
	outtree->Branch("jet_neutralMultiplicity", &jet_neutralMultiplicity, "jet_neutralMultiplicity/F");
	outtree->Branch("jet_Chadfrac", &jet_Chadfrac, "jet_Chadfrac/F");
	outtree->Branch("jet_Nhadfrac", &jet_Nhadfrac, "jet_Nhadfrac/F");
	outtree->Branch("jet_Phofrac", &jet_Phofrac, "jet_Phofrac/F");
	outtree->Branch("jet_Mufrac", &jet_Mufrac, "jet_Mufrac/F");
	outtree->Branch("jet_Elefrac", &jet_Elefrac, "jet_Elefrac/F");
	outtree->Branch("jet_pfloose", &jet_pfloose, "jet_pfloose/I");
*/
	outtree->Branch("jet_dPhiMet", &jet_dPhiMet, "jet_dPhiMet/F");
	outtree->Branch("jet_genDR", &jet_genDR, "jet_genDR/F");
	outtree->Branch("jet_genPt", &jet_genPt, "jet_genPt/F");
	outtree->Branch("jet_genE", &jet_genE, "jet_genE/F");
	outtree->Branch("jet_genR", &jet_genR, "jet_genR/F");
	outtree->Branch("jet_prtDR", &jet_prtDR, "jet_prtDR/F");
	outtree->Branch("jet_prtPt", &jet_prtPt, "jet_prtPt/F");
	outtree->Branch("jet_prtE", &jet_prtE, "jet_prtE/F");
	outtree->Branch("jet_prtR", &jet_prtR, "jet_prtR/F");
	outtree->Branch("jet_index", &jet_index, "jet_index/I");
	outtree->Branch("ev_met_corr_pfmet", &ev_met_corr_pfmet, "ev_met_corr_pfmet/F");
	outtree->Branch("ev_met_corr_phi_pfmet", &ev_met_corr_phi_pfmet, "ev_met_corr_phi_pfmet/F");
	outtree->Branch("ev_met_corr_eta_pfmet", &ev_met_corr_eta_pfmet, "ev_met_corr_eta_pfmet/F");
	outtree->Branch("ev_met_corr_e_pfmet", &ev_met_corr_e_pfmet, "ev_met_corr_e_pfmet/F");
/*	outtree->Branch("ev_met_pfmet", &ev_met_pfmet, "ev_met_pfmet/F");
	outtree->Branch("ev_met_phi_pfmet", &ev_met_phi_pfmet, "ev_met_phi_pfmet/F");
	outtree->Branch("ev_met_sumet_pfmet", &ev_met_sumet_pfmet, "ev_met_sumet_pfmet/F");
	outtree->Branch("ev_met_mEtSig_pfmet", &ev_met_mEtSig_pfmet, "ev_met_mEtSig_pfmet/F");
	outtree->Branch("ev_met_significance_pfmet", &ev_met_significance_pfmet, "ev_met_significance_pfmet/F");
	outtree->Branch("ev_met_corrMet", &ev_met_corrMet, "ev_met_corrMet/F");
	outtree->Branch("ev_met_corrMetPhi", &ev_met_corrMetPhi, "ev_met_corrMetPhi/F");
*/
	outtree->Branch("ev_pu_n", &ev_pu_n, "ev_pu_n/F");
	outtree->Branch("ev_nvtx", &ev_nvtx, "ev_nvtx/F");
	outtree->Branch("ev_rho", &ev_rho, "ev_rho/F");
	outtree->Branch("ev_weight", &ev_weight, "ev_weight/F");
	outtree->Branch("ev_evweight", &ev_evweight, "ev_evweight/F");
	outtree->Branch("ev_pu_weight", &ev_pu_weight, "ev_pu_weight/F");

	int np[20] = {0};
	int decade = 0;
	int totevents = intree->GetEntries();
	cout << "#entries= " << totevents << endl;
	// loop over events
	for(int ievt=0 ; ievt < totevents ; ievt++)
	{
		double progress = 10.0*ievt/(1.0*totevents);
		int k = TMath::FloorNint(progress);
		if (k > decade) cout<<10*k<<" %"<<endl;
		decade = k;

		intree->GetEntry(ievt);
	
		// take only the subset of events where at least one jet remains
		if( njets_passing_kLooseID < 1 )
			continue;
		np[0]++;

		// gen-level 4-momenta must be non-zero
		TLorentzVector gen_b1, gen_b2, gen_j1, gen_j2;
		if( gr_b1_p4_pt < .1 ) continue;
		if( gr_b2_p4_pt < .1 ) continue;
		if( gr_j1_p4_pt < .1 ) continue;
		if( gr_j2_p4_pt < .1 ) continue;
		np[1]++;

		gen_b1.SetPtEtaPhiM(gr_b1_p4_pt, gr_b1_p4_eta, gr_b1_p4_phi, gr_b1_p4_mass);
		gen_b2.SetPtEtaPhiM(gr_b2_p4_pt, gr_b2_p4_eta, gr_b2_p4_phi, gr_b2_p4_mass);
		gen_j1.SetPtEtaPhiM(gr_j1_p4_pt, gr_j1_p4_eta, gr_j1_p4_phi, gr_j1_p4_mass);
		gen_j2.SetPtEtaPhiM(gr_j2_p4_pt, gr_j2_p4_eta, gr_j2_p4_phi, gr_j2_p4_mass);
		float DR_jet_b1, DR_jet_b2;
		float DR_jet_j1, DR_jet_j2;
		TLorentzVector jet;
		TLorentzVector met;
		met.SetPtEtaPhiE(met_corr_pfmet, met_corr_eta_pfmet, met_corr_phi_pfmet, met_corr_e_pfmet);

		// loop over jets, store jet info + info on closest genjet / parton (no selection applied)
		for( int ijet = 0 ; ijet < min(njets_passing_kLooseID, 4); ijet ++ )
		{
			if( ijet == 0 )
			{
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
				jet_nConstituents = (float)(jet_nNeutrals + jet_nCharged);
				jet_axis1 = j1_axis1;
				jet_axis2 = j1_axis2;
				jet_pull = j1_pull;
/*				jet_Rchg = j1_Rchg;
				jet_Rneutral = j1_Rneutral;
				jet_R = j1_R;
				jet_chargedMultiplicity = j1_chargedMultiplicity;
				jet_neutralMultiplicity = j1_neutralMultiplicity;
				jet_Chadfrac = j1_Chadfrac;
				jet_Nhadfrac = j1_Nhadfrac;
				jet_Phofrac = j1_Phofrac;
				jet_Mufrac = j1_Mufrac;
				jet_Elefrac = j1_Elefrac;
				jet_dPhiMet = j1_dPhiMet;
				jet_pfloose = j1_pfloose;
*/
				jet_index = 1;
				jet.SetPtEtaPhiE(j1_pt, j1_eta, j1_phi, j1_e);
				jet_dPhiMet = jet.DeltaPhi(met);
			} // end if jet == 0

			if( ijet == 1 )
			{
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
				jet_nConstituents = (float)(jet_nNeutrals + jet_nCharged);
				jet_axis1 = j2_axis1;
				jet_axis2 = j2_axis2;
				jet_pull = j2_pull;
/*				jet_Rchg = j2_Rchg;
				jet_Rneutral = j2_Rneutral;
				jet_R = j2_R;
				jet_chargedMultiplicity = j2_chargedMultiplicity;
				jet_neutralMultiplicity = j2_neutralMultiplicity;
				jet_Chadfrac = j2_Chadfrac;
				jet_Nhadfrac = j2_Nhadfrac;
				jet_Phofrac = j2_Phofrac;
				jet_Mufrac = j2_Mufrac;
				jet_Elefrac = j2_Elefrac;
				jet_dPhiMet = j2_dPhiMet;
				jet_pfloose = j2_pfloose;
*/
				jet_index = 2;
				jet.SetPtEtaPhiE(j2_pt, j2_eta, j2_phi, j2_e);
				jet_dPhiMet = jet.DeltaPhi(met);
			} // end if jet == 1

			if( ijet == 2 )
			{
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
				jet_nConstituents = (float)(jet_nNeutrals + jet_nCharged);
				jet_axis1 = j3_axis1;
				jet_axis2 = j3_axis2;
				jet_pull = j3_pull;
/*				jet_Rchg = j3_Rchg;
				jet_Rneutral = j3_Rneutral;
				jet_R = j3_R;
				jet_chargedMultiplicity = j3_chargedMultiplicity;
				jet_neutralMultiplicity = j3_neutralMultiplicity;
				jet_Chadfrac = j3_Chadfrac;
				jet_Nhadfrac = j3_Nhadfrac;
				jet_Phofrac = j3_Phofrac;
				jet_Mufrac = j3_Mufrac;
				jet_Elefrac = j3_Elefrac;
				jet_dPhiMet = j3_dPhiMet;
				jet_pfloose = j3_pfloose;
*/
				jet_index = 3;
				jet.SetPtEtaPhiE(j3_pt, j3_eta, j3_phi, j3_e);
				jet_dPhiMet = jet.DeltaPhi(met);
			} // end if jet == 2

			if( ijet == 3 )
			{
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
				jet_nConstituents = (float)(jet_nNeutrals + jet_nCharged);
				jet_axis1 = j4_axis1;
				jet_axis2 = j4_axis2;
				jet_pull = j4_pull;
/*				jet_Rchg = j4_Rchg;
				jet_Rneutral = j4_Rneutral;
				jet_R = j4_R;
				jet_chargedMultiplicity = j4_chargedMultiplicity;
				jet_neutralMultiplicity = j4_neutralMultiplicity;
				jet_Chadfrac = j4_Chadfrac;
				jet_Nhadfrac = j4_Nhadfrac;
				jet_Phofrac = j4_Phofrac;
				jet_Mufrac = j4_Mufrac;
				jet_Elefrac = j4_Elefrac;
				jet_dPhiMet = j4_dPhiMet;
				jet_pfloose = j4_pfloose;
*/
				jet_index = 4;
				jet.SetPtEtaPhiE(j4_pt, j4_eta, j4_phi, j4_e);
				jet_dPhiMet = jet.DeltaPhi(met);
			} // end if jet == 3


			// generator level info
			DR_jet_b1 = jet.DeltaR(gen_b1);
			DR_jet_b2 = jet.DeltaR(gen_b2);
			DR_jet_j1 = jet.DeltaR(gen_j1);
			DR_jet_j2 = jet.DeltaR(gen_j2);
			jet_genDR = DR_jet_j1;
			jet_genPt = gen_j1.Pt();
			jet_genE = gen_j1.Energy();
			jet_genR = (float)gen_j1.Pt() / (float)jet_pt;
			if(DR_jet_j2 < DR_jet_j1)
			{
				jet_genDR = DR_jet_j2;
				jet_genPt = gen_j2.Pt();
				jet_genE = gen_j2.Energy();
				jet_genR = (float)gen_j2.Pt() / (float)jet_pt;
			}
			jet_prtDR = DR_jet_b1;
			jet_prtPt = gen_b1.Pt();
			jet_prtE = gen_b1.Energy();
			jet_prtR = (float)gen_b1.Pt() / (float)jet_pt;
			if(DR_jet_b2 < DR_jet_b1)
			{
				jet_prtDR = DR_jet_b2;
				jet_prtPt = gen_b2.Pt();
				jet_prtE = gen_b2.Energy();
				jet_prtR = (float)gen_b2.Pt() / (float)jet_pt;
			}
/*			ev_met_pfmet = met_pfmet;
			ev_met_phi_pfmet = met_phi_pfmet;
			ev_met_sumet_pfmet = met_sumet_pfmet;
			ev_met_mEtSig_pfmet = met_mEtSig_pfmet;
			ev_met_significance_pfmet = met_significance_pfmet;
			ev_met_corrMet = met_corrMet;
			ev_met_corrMetPhi = met_corrMetPhi;
*/
			ev_met_corr_pfmet = met_corr_pfmet;
			ev_met_corr_phi_pfmet = met_corr_phi_pfmet;
			ev_met_corr_eta_pfmet = met_corr_eta_pfmet;
			ev_met_corr_e_pfmet = met_corr_e_pfmet;
			ev_pu_n = pu_n;
			ev_nvtx = nvtx;
			ev_rho = rho;
			ev_weight = weight;
			ev_evweight = evweight;
			ev_pu_weight = pu_weight;
			np[2]++;
			outtree->Fill();
		} // end of loop over jets

	} // end of loop over events


	for(int i=0 ; i < 3 ; i++)
		cout << "#np[" << i << "]= " << np[i] << endl;

	outfile->cd();
	outtree->Write();
	outfile->Close();
	infile->Close();

	return 0;
}

