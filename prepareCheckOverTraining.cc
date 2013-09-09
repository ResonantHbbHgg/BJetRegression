// ROOT HEADERS
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TClonesArray.h"
// TMVA HEADERS
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
// C++ HEADERS
#include <string>
#include <iostream>
#include <sstream>
// DEFINES
#define DEBUG 0
#define USEHT 1
// namespaces
using namespace std;
using namespace TMVA;

int main(int argc, char *argv[])
{
	if(DEBUG) cout << "DEBUG: Initialisation: reading parameters" << endl;
	cout << "argc= " << argc << endl;
	for(int iarg = 0 ; iarg < argc; iarg++)
		cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
		string syntax = Form("WARNING: Syntax is %s -i (inputfile) -it (inputtree) -o (outputfile) -rf (regressionfile) -n (numberOfSplit) -j (treeSplit)", argv[0]);

	if( argc == 1 )
	{
		cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
		cerr << syntax << endl;
	}

	int nTrainingTrees = 1;
	string inputfile = "2013-08-07_jetTreeForTraining_m300.root";
	string inputtree = "jets";
	float w1 = 1.0;
	string inputfile2 = "jetTreeForTraining_m500.root";
	string inputtree2 = "jets";
	float w2 = 1.0;
	string inputfile3 = "jetTreeForTraining_m700.root";
	string inputtree3 = "jets";
	float w3 = 1.0;
	string inputfile4 = "jetTreeForTraining_m1000.root";
	string inputtree4 = "jets";
	float w4 = 1.0;
	string outputfile = "jetTreeForOverTrainingChecks_m300.root";
	string regressionfile = "weights/test_BDT.weights.xml";
	int n = 1;
	int j = 0;

	for(int iarg=0 ; iarg < argc ; iarg++)
	{
		if(strcmp("-i", argv[iarg]) == 0 && argc >= iarg + 1)
			inputfile = argv[iarg+1];
		if(strcmp("-i2", argv[iarg]) == 0 && argc >= iarg + 1)
			inputfile2 = argv[iarg+1];
		if(strcmp("-i3", argv[iarg]) == 0 && argc >= iarg + 1)
			inputfile3 = argv[iarg+1];
		if(strcmp("-i4", argv[iarg]) == 0 && argc >= iarg + 1)
			inputfile4 = argv[iarg+1];
		if(strcmp("-it", argv[iarg]) == 0 && argc >= iarg + 1)
			inputtree = argv[iarg+1];
		if(strcmp("-it2", argv[iarg]) == 0 && argc >= iarg + 1)
			inputtree2 = argv[iarg+1];
		if(strcmp("-it3", argv[iarg]) == 0 && argc >= iarg + 1)
			inputtree3 = argv[iarg+1];
		if(strcmp("-it4", argv[iarg]) == 0 && argc >= iarg + 1)
			inputtree4 = argv[iarg+1];
		if(strcmp("-w1", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> w1; }
		if(strcmp("-w2", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> w2; }
		if(strcmp("-w3", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> w3; }
		if(strcmp("-w4", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> w4; }
		if(strcmp("-n", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> n; }
		if(strcmp("-j", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> j; }
		if(strcmp("-o", argv[iarg]) == 0 && argc >= iarg + 1)
			outputfile = argv[iarg+1];
		if(strcmp("-rf", argv[iarg]) == 0 && argc >= iarg + 1)
			regressionfile = argv[iarg+1];
		if(strcmp("-ntt", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> nTrainingTrees; }
		if(strcmp("--help", argv[iarg]) == 0 || strcmp("-h", argv[iarg]) == 0)
		{
			cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
			cerr << syntax << endl;
			cerr << "inputfile= " << inputfile << endl;
			cerr << "inputtree= " << inputtree << endl;
			cerr << "outputfile= " << outputfile << endl;
			cerr << "regressionfile= " << regressionfile << endl;
			return 2;
		}
	}

	cout << "inputfile= " << inputfile << endl;
	cout << "inputtree= " << inputtree << endl;
	cout << "outputfile= " << outputfile << endl;
	cout << "regressionfile= " << regressionfile << endl;

	TFile *infile = TFile::Open(inputfile.c_str());
	TTree *intree = (TTree*)infile->Get(inputtree.c_str());
/*
	TFile *infile2 = TFile::Open(inputfile2.c_str());
	TTree *intree2 = (TTree*)infile2->Get(inputtree2.c_str());
	TFile *infile3 = TFile::Open(inputfile3.c_str());
	TTree *intree3 = (TTree*)infile3->Get(inputtree3.c_str());
	TFile *infile4 = TFile::Open(inputfile4.c_str());
	TTree *intree4 = (TTree*)infile4->Get(inputtree4.c_str());
*/
	TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
	TTree *outtree = new TTree(inputtree.c_str(), inputtree.c_str());

	float event, jet_pt, jet_eta, jet_emfrac, jet_nConstituents, jet_hadfrac, jet_secVtxPt, jet_secVtx3dL, ev_met_corr_pfmet, jet_dPhiMet, regjet_pt, jet_genDR, jet_csvBtag;
	float ph1_pt, ph2_pt, ev_rho, HT_gg;
	int lot;
	intree->SetBranchAddress("event", &event);
	intree->SetBranchAddress("jet_pt", &jet_pt);
	intree->SetBranchAddress("jet_eta", &jet_eta);
	intree->SetBranchAddress("jet_emfrac", &jet_emfrac);
	intree->SetBranchAddress("jet_nConstituents", &jet_nConstituents);
	intree->SetBranchAddress("jet_hadfrac", &jet_hadfrac);
	intree->SetBranchAddress("jet_secVtxPt", &jet_secVtxPt);
	intree->SetBranchAddress("jet_secVtx3dL", &jet_secVtx3dL);
	intree->SetBranchAddress("ev_met_corr_pfmet", &ev_met_corr_pfmet);
	intree->SetBranchAddress("jet_dPhiMet", &jet_dPhiMet);
	intree->SetBranchAddress("jet_genDR", &jet_genDR);
	intree->SetBranchAddress("jet_csvBtag", &jet_csvBtag);
	intree->SetBranchAddress("ph1_pt", &ph1_pt);
	intree->SetBranchAddress("ph2_pt", &ph2_pt);
	intree->SetBranchAddress("ev_rho", &ev_rho);
	outtree->Branch("event", &event, "event/F");
	outtree->Branch("jet_pt", &jet_pt, "jet_pt/F");
	outtree->Branch("jet_eta", &jet_eta, "jet_eta/F");
	outtree->Branch("jet_emfrac", &jet_emfrac, "jet_emfrac/F");
	outtree->Branch("jet_nConstituents", &jet_nConstituents, "jet_nConstituents/F");
	outtree->Branch("jet_hadfrac", &jet_hadfrac, "jet_hadfrac/F");
	outtree->Branch("jet_secVtxPt", &jet_secVtxPt, "jet_secVtxPt/F");
	outtree->Branch("jet_secVtx3dL", &jet_secVtx3dL, "jet_secVtx3dL/F");
	outtree->Branch("ev_met_corr_pfmet", &ev_met_corr_pfmet, "ev_met_corr_pfmet/F");
	outtree->Branch("jet_dPhiMet", &jet_dPhiMet, "jet_dPhiMet/F");
	outtree->Branch("jet_genDR", &jet_genDR, "jet_genDR/F");
	outtree->Branch("jet_csvBtag", &jet_csvBtag, "jet_csvBtag/F");
	outtree->Branch("ph1_pt", &ph1_pt, "ph1_pt/F");
	outtree->Branch("ph2_pt", &ph2_pt, "ph2_pt/F");
	outtree->Branch("HT_gg", &HT_gg, "HT_gg/F");
	outtree->Branch("ev_rho", &ev_rho, "ev_rho/F");
	outtree->Branch("regjet_pt", &regjet_pt, "regjet_pt/F");
	outtree->Branch("lot", &lot, "lot/I");

	TMVA::Reader* readerRegres_0 = new TMVA::Reader( "!Color:!Silent" );
	readerRegres_0->AddVariable( "jet_pt", &jet_pt);
	readerRegres_0->AddVariable( "jet_eta", &jet_eta);
	readerRegres_0->AddVariable( "jet_emfrac", &jet_emfrac);
	readerRegres_0->AddVariable( "jet_nConstituents", &jet_nConstituents);
	readerRegres_0->AddVariable( "jet_hadfrac", &jet_hadfrac);
	readerRegres_0->AddVariable( "jet_secVtxPt", &jet_secVtxPt);
	readerRegres_0->AddVariable( "jet_secVtx3dL", &jet_secVtx3dL);
	readerRegres_0->AddVariable( "ev_met_corr_pfmet", &ev_met_corr_pfmet);
	readerRegres_0->AddVariable( "jet_dPhiMet", &jet_dPhiMet);
// Adding variables
	readerRegres_0->AddVariable( "ev_rho", &ev_rho);
	if( USEHT ) readerRegres_0->AddVariable( "ph1_pt+ph2_pt", &HT_gg);
//	readerRegres_0->BookMVA("BDT", "weights/2013-08-08_test_n3_j0_BDT.weights.xml");
	readerRegres_0->BookMVA("BDT", regressionfile.c_str());
/*
	TMVA::Reader* readerRegres_1 = new TMVA::Reader( "!Color:!Silent" );
	readerRegres_1->AddVariable( "jet_pt", &jet_pt);
	readerRegres_1->AddVariable( "jet_eta", &jet_eta);
	readerRegres_1->AddVariable( "jet_emfrac", &jet_emfrac);
	readerRegres_1->AddVariable( "jet_nConstituents", &jet_nConstituents);
	readerRegres_1->AddVariable( "jet_hadfrac", &jet_hadfrac);
	readerRegres_1->AddVariable( "jet_secVtxPt", &jet_secVtxPt);
	readerRegres_1->AddVariable( "jet_secVtx3dL", &jet_secVtx3dL);
	readerRegres_1->AddVariable( "ev_met_corr_pfmet", &ev_met_corr_pfmet);
	readerRegres_1->AddVariable( "jet_dPhiMet", &jet_dPhiMet);
	readerRegres_1->BookMVA("BDT", "weights/2013-08-08_test_n3_j1_BDT.weights.xml");

	TMVA::Reader* readerRegres_2 = new TMVA::Reader( "!Color:!Silent" );
	readerRegres_2->AddVariable( "jet_pt", &jet_pt);
	readerRegres_2->AddVariable( "jet_eta", &jet_eta);
	readerRegres_2->AddVariable( "jet_emfrac", &jet_emfrac);
	readerRegres_2->AddVariable( "jet_nConstituents", &jet_nConstituents);
	readerRegres_2->AddVariable( "jet_hadfrac", &jet_hadfrac);
	readerRegres_2->AddVariable( "jet_secVtxPt", &jet_secVtxPt);
	readerRegres_2->AddVariable( "jet_secVtx3dL", &jet_secVtx3dL);
	readerRegres_2->AddVariable( "ev_met_corr_pfmet", &ev_met_corr_pfmet);
	readerRegres_2->AddVariable( "jet_dPhiMet", &jet_dPhiMet);
	readerRegres_2->BookMVA("BDT", "weights/2013-08-08_test_n3_j2_BDT.weights.xml");
*/
	int totevents = intree->GetEntries();
	cout << "totevents= " << totevents << endl;
	for(int ievt = 0 ; ievt < totevents ; ievt++)
	{
		intree->GetEntry(ievt);
		lot = (int)event % 3;
//		if(jet_genDR >= 0.4) continue;
//		if(jet_csvBtag <= 0.) continue;
//		if( (int)event % 3 == 0)
//			regjet_pt = readerRegres_2->EvaluateMVA("BDT");	
//		if( (int)event % 3 == 1)
			HT_gg = ph1_pt + ph2_pt;
			regjet_pt = readerRegres_0->EvaluateMVA("BDT");	
//		if( (int)event % 3 == 2)
//			regjet_pt = readerRegres_0->EvaluateMVA("BDT");	
		if(DEBUG) cout << "ievt= " << ievt << "\tevent= " << event << "\tinput= " << jet_pt << "\toutput= " << regjet_pt << endl; // Olivier regression
		outtree->Fill();
	}

	outfile->cd();
	outtree->Write();
	outfile->Close();
	infile->Close();

 	return 0;
}
