// C++ headers
#include <iostream>
#include <string>
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
// RooFit headers
// local files
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;

int main(int argc, char *argv[])
{
	cout << "argc= " << argc << endl;
	for(int iarg = 0 ; iarg < argc; iarg++)
		cout << "argv[" << iarg << "]= " << argv[iarg] << endl;

	if( argc == 1 )
	{
		cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
		cerr << "WARNING: Syntax is " << argv[0] << " -i (inputfile) -it (inputtree) -o (outputfile) -ot (outputtree)" << endl;
	}
	
	string inputfile = "Data_m300_StandardFullSelection_v2.root";
	string inputtree = "Data";
	string outputfile = "Data_m300_test_minimal.root";
	string outputtree = "TCVARS";

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
	TTree *outtree = new TTree(outputtree.c_str(), Form("%s minimal", outputtree.c_str()));


	float mgg, mjj, mtot;
	int cut_based_ct, njets_kRadionID_and_CSVM, selection_cut_level;
	float evWeight;
	intree->SetBranchAddress("gg_mass", &mgg);
	intree->SetBranchAddress("jj_mass", &mjj);
	intree->SetBranchAddress("ggjj_mass", &mtot);
	intree->SetBranchAddress("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM);
	intree->SetBranchAddress("selection_cut_level", &selection_cut_level);
	intree->SetBranchAddress("evweight", &evWeight);

	outtree->Branch("mgg", &mgg, "mgg/F");
	outtree->Branch("mjj", &mjj, "mjj/F");
	outtree->Branch("mtot", &mtot, "mtot/F");
	outtree->Branch("cut_based_ct", &cut_based_ct, "cut_based_ct/I");
	outtree->Branch("evWeight", &evWeight, "evWeight/F");

	for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
	{
		intree->GetEntry(ievt);
		// FIXME: add possible further selection
		if( njets_kRadionID_and_CSVM >= 2 ) cut_based_ct = 0;
		if( njets_kRadionID_and_CSVM == 1 ) cut_based_ct = 1;
		outtree->Fill();
	}
  outfile->cd();
  outtree->Write();
  outfile->Close();
  infile->Close();

	return 0;
}
