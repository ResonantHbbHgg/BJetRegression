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
	int cut_based_ct;
	float evWeight;

	float gg_mass, jj_mass, ggjj_mass;
	int njets_kRadionID_and_CSVM;
	float weight;
	intree->SetBranchAddress("gg_mass", &gg_mass);
	intree->SetBranchAddress("jj_mass", &jj_mass);
	intree->SetBranchAddress("ggjj_mass", &ggjj_mass);
	intree->SetBranchAddress("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM);
	intree->SetBranchAddress("weight", &weight);

	outtree->Branch("jj_mass", &jj_mass, "jj_mass/F");
	outtree->Branch("ggjj_mass", &ggjj_mass, "ggjj_mass/F");
	outtree->Branch("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM, "njets_kRadionID_and_CSVM/I");

	for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
	{
		intree->GetEntry(ievt);
		outtree->Fill();
	}
  outfile->cd();
  outtree->Write();
  outfile->Close();
/*
//	TFile *outfilereg = new TFile("simple_reg_parton.root", "RECREATE");
	TFile *outfilereg = new TFile("simple_reg_genjet_globeinputs.root", "RECREATE");
	TTree *outtreereg = new TTree("Radion_m300_8TeV_nm", "Radion_m300_8TeV_nm reg");
	outtreereg->Branch("jj_mass", &regjj_mass, "jj_mass/F");
	outtreereg->Branch("ggjj_mass", &regggjj_mass, "ggjj_mass/F");
	outtreereg->Branch("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM, "njets_kRadionID_and_CSVM/I");
	for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
	{
		intree->GetEntry(ievt);
		outtreereg->Fill();
	}

	outfilereg->cd();
  outtreereg->Write();
  outfilereg->Close();

//	TFile *outfileregMLP = new TFile("simple_regMLP_parton.root", "RECREATE");
	TFile *outfileregMLP = new TFile("simple_regMLP_genjet_globeinputs.root", "RECREATE");
	TTree *outtreeregMLP = new TTree("Radion_m300_8TeV_nm", "Radion_m300_8TeV_nm regMLP");
	outtreeregMLP->Branch("jj_mass", &regMLPjj_mass, "jj_mass/F");
	outtreeregMLP->Branch("ggjj_mass", &regMLPggjj_mass, "ggjj_mass/F");
	outtreeregMLP->Branch("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM, "njets_kRadionID_and_CSVM/I");
	for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
	{
		intree->GetEntry(ievt);
		outtreeregMLP->Fill();
	}

	outfileregMLP->cd();
  outtreeregMLP->Write();
  outfileregMLP->Close();
*/
  infile->Close();

	return 0;
}
