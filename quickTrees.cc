// C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
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
		cerr << "WARNING: Syntax is " << argv[0] << " -i (inputfile) -it (inputtree) -o (outputfile) -ot (outputtree) -wj (whichJet) -cutLevel (cutLevel) -fs (fitStrategy) -m (mass)" << endl;
	}
	
	string inputfile = "Data_m300_StandardFullSelection_v2.root";
	string inputtree = "Data";
	string outputfile = "Data_m300_test_minimal.root";
	string outputtree = "TCVARS";
	string whichJet = "";
	string fitStrategy = "mgg";
	int cutLevel = 0;
	int mass = 300;

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
		if(strcmp("-wj", argv[iarg]) == 0 && argc >= iarg + 1)
			whichJet = argv[iarg+1];
		if(strcmp("-fs", argv[iarg]) == 0 && argc >= iarg + 1)
			fitStrategy = argv[iarg+1];
		if(strcmp("-cutLevel", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> cutLevel; }
		if(strcmp("-m", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> mass; }
		if((strcmp("-h", argv[iarg]) == 0) || (strcmp("--help", argv[iarg]) == 0))
		{
			cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
			cerr << "WARNING: Syntax is " << argv[0] << " -i (inputfile) -it (inputtree) -o (outputfile) -ot (outputtree) -wj (whichJet) -cutLevel (cutLevel) -fs (fitStrategy) -m (mass)" << endl;
			cerr << "inputfile= " << inputfile << endl;
			cerr << "inputtree= " << inputtree << endl;
			cerr << "outputfile= " << outputfile << endl;
			cerr << "outputtree= " << outputtree << endl;
			cerr << "whichJet= " << whichJet << endl;
			cerr << "cutLevel= " << cutLevel << endl;
			cerr << "fitStrategy= " << fitStrategy << endl;
			cerr << "mass= " << mass << endl;
			return 2;
		}
	}

	if(strcmp(whichJet.c_str(), "base") == 0) whichJet="";

	cout << "inputfile= " << inputfile << endl;
	cout << "inputtree= " << inputtree << endl;
	cout << "outputfile= " << outputfile << endl;
	cout << "outputtree= " << outputtree << endl;
	cout << "whichJet= " << whichJet << endl;
	cout << "cutLevel= " << cutLevel << endl;
	cout << "fitStrategy= " << fitStrategy << endl;
	cout << "mass= " << mass << endl;

	TFile *infile = TFile::Open(inputfile.c_str());
	TTree *intree = (TTree*)infile->Get(inputtree.c_str());
	TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
	TTree *outtree = new TTree(outputtree.c_str(), Form("%s minimal", outputtree.c_str()));


	float mgg, mjj, mtot;
	int cut_based_ct, njets_kRadionID_and_CSVM, selection_cut_level;
	float evWeight;
	float regcosthetastar, minDRgregkinj;
	int njets_kLooseID;
	intree->SetBranchAddress("gg_mass", &mgg);
	intree->SetBranchAddress(Form("%sjj_mass", whichJet.c_str()), &mjj);
	intree->SetBranchAddress(Form("%sggjj_mass", whichJet.c_str()), &mtot);
	intree->SetBranchAddress("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM);
	intree->SetBranchAddress("selection_cut_level", &selection_cut_level);
	intree->SetBranchAddress("weight", &evWeight);
	intree->SetBranchAddress("regcosthetastar", &regcosthetastar);
	intree->SetBranchAddress("minDRgregkinj", &minDRgregkinj);
	intree->SetBranchAddress("njets_kLooseID", &njets_kLooseID);
	

	outtree->Branch("mgg", &mgg, "mgg/F");
	outtree->Branch("mjj", &mjj, "mjj/F");
	outtree->Branch("mtot", &mtot, "mtot/F");
	outtree->Branch("cut_based_ct", &cut_based_ct, "cut_based_ct/I");
	outtree->Branch("evWeight", &evWeight, "evWeight/F");

//	cout << "strcmp mgg= " << (strcmp("mgg", fitStrategy.c_str()) ) << endl;
//	cout << "strcmp mggjj= " << (strcmp("mggjj", fitStrategy.c_str()) ) << endl;

	for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
	{
		intree->GetEntry(ievt);
// EXTRA CUTS
//		if( selection_cut_level < cutLevel ) continue; // hard-coded in the trees, out of date wrt to the rest of the cuts
	if( cutLevel > 5)
	{
		if(fabs(regcosthetastar) >= .9) continue;
		if(minDRgregkinj <= 1.) continue;
		if( njets_kLooseID >= 4 ) continue;
	}

// FITTING THE MGG SPECTRUM
		if( strcmp("mgg", fitStrategy.c_str()) == 0 )
		{
			if( njets_kRadionID_and_CSVM >= 2 )
			{
				if( (strcmp("", whichJet.c_str()) == 0) || (strcmp("kin", whichJet.c_str()) == 0) )
					if( mjj < 95. || mjj > 175. ) continue;
				if( (strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
					if( mjj < 90. || mjj > 150. ) continue;
			}
			if( njets_kRadionID_and_CSVM == 1 )
			{
				if( (strcmp("", whichJet.c_str()) == 0) || (strcmp("kin", whichJet.c_str()) == 0) )
					if( mjj < 100. || mjj > 160. ) continue;
				if( (strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
					if( mjj < 95. || mjj > 140. ) continue;
			}
			if( mass == 300 )
			{
				if( (strcmp("", whichJet.c_str()) == 0) || (strcmp("kin", whichJet.c_str()) == 0) )
				{
					if( njets_kRadionID_and_CSVM >= 2 && (mtot < 255. || mtot > 320.) ) continue;
					if( njets_kRadionID_and_CSVM == 1 && (mtot < 260. || mtot > 335.) ) continue;
				}
				if( (strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
				{
					if( njets_kRadionID_and_CSVM >= 2 && (mtot < 260. || mtot > 335.) ) continue;
					if( njets_kRadionID_and_CSVM == 1 && (mtot < 265. || mtot > 345.) ) continue;
				}
			}
			if( mass == 500 && (mtot < 465. || mtot > 535.) ) continue;
			if( mass == 700 && (mtot < 660. || mtot > 740.) ) continue;
			if( mass == 1000 && (mtot < 955. || mtot > 1055.) ) continue;
		}

// FITTING THE MGGJJ SPECTRUM
		if( strcmp("mggjj", fitStrategy.c_str()) == 0 )
		{
			if( mgg < 120. || mgg > 130. ) continue;
			if( njets_kRadionID_and_CSVM >= 2 )
			{
				if( (strcmp("", whichJet.c_str()) == 0) || (strcmp("kin", whichJet.c_str()) == 0) )
					if( mjj < 95. || mjj > 150. ) continue;
				if( (strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
					if( mjj < 100. || mjj > 160. ) continue;
			}
			if( njets_kRadionID_and_CSVM == 1 ) 
				if( mjj < 90. || mjj > 170. ) continue;
		}
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
