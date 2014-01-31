// macro to compare final trees to check if selection is correct
// Olivier Bondu, January 2014
// C++ headers
#include <iostream>
#include <string>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
// Analysis headers
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;

TCanvas c1;

void display(string label, float a, float b)
{
	float r;
	if( a != 0 ) r = (a - b) / a;
	else r = 0;
	float r100 = r * 100;
	cout << label << ":\tdiff= " << r100 << " %\t(" << a << " - " << b << ") / " << a << endl;
	return;
}

float getMean(TTree *t, string var, string cut)
{
	c1.Clear();
	t->Draw(Form("%s>>htemp",var.c_str()), cut.c_str());
	float mean = ((TH1F*)gDirectory->Get("htemp"))->GetMean();
	((TH1F*)gDirectory->Get("htemp"))->Delete();
	return mean;
}

float getRMS(TTree *t, string var, string cut)
{
	c1.Clear();
	t->Draw(Form("%s>>htemp",var.c_str()), cut.c_str());
	float rms = ((TH1F*)gDirectory->Get("htemp"))->GetRMS();
	((TH1F*)gDirectory->Get("htemp"))->Delete();
	return rms;
}

int main()
{
// Radion m300 mgg
//	string file1 = "/afs/cern.ch/user/c/crovelli/public/4Alexandra/trees/v22/finalizedTrees_Radion_V07__fitToGG__noKinFit/Radion_m300.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/v22_fitToMgg_noKinFit/Radion_m300_8TeV_m300.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/minimum_test.root"; 
// Data m300 mgg
//	string file1 = "/afs/cern.ch/user/c/crovelli/public/4Alexandra/trees/v22/finalizedTrees_Radion_V07__fitToGG__noKinFit/data__selez300.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/minimum_test_data.root"; 
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/v22_fitToMgg_noKinFit/Data_m300.root";

// Radion m650 mggjj
//	string file1 = "/afs/cern.ch/user/c/crovelli/public/4Alexandra/trees/v22/finalizedTrees_Radion_V07__fitToGGJJ__withKinFit/Radion_m650.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/minimum_test_650.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/v22_fitToMggjj_withKinFit/Radion_m650_8TeV_m650.root";
// Data m650 mggjj
//	string file1 = "/afs/cern.ch/user/c/crovelli/public/4Alexandra/trees/v22/finalizedTrees_Radion_V07__fitToGGJJ__withKinFit/data.root";
//	string file1 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/minimum_test_PLOUF_data.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/v22_fitToMggjj_withKinFit/Data_m650.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/minimum_test_650_data.root";

// Radion m650 mggjj wokinfit
//	string file1 = "/afs/cern.ch/user/c/crovelli/public/4Alexandra/trees/v22/finalizedTrees_Radion_V07__fitToGGJJ__noKinFit/Radion_m650.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/v22_fitToMggjj_noKinFit/Radion_m650_8TeV_m650.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/minimum_test_650_wokinfit.root";
// Data m650 mggjj
	string file1 = "/afs/cern.ch/user/c/crovelli/public/4Alexandra/trees/v22/finalizedTrees_Radion_V07__fitToGGJJ__noKinFit/data.root";
	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/v22_fitToMggjj_noKinFit/Data_m650.root";
//	string file2 = "/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/BJetRegression/minimum_test_650_data_wokinfit.root";

	string label1 = "Chiara";
	string label2 = "Olivier";
	TFile *f1 = TFile::Open(file1.c_str());
	TFile *f2 = TFile::Open(file2.c_str());
	TTree *t1 = (TTree*)f1->Get("TCVARS");
	TTree *t2 = (TTree*)f2->Get("TCVARS");

	cout << "Results will be (" << label1 << " - " << label2 << ") / " << label1 << endl;
	cout << "##### ##### YIELDS ##### #####" << endl;
	display("entries tot", t1->GetEntries(), t2->GetEntries());
	display("entries 2btag", t1->GetEntries("cut_based_ct == 0"), t2->GetEntries("cut_based_ct == 0"));
	display("entries 1btag", t1->GetEntries("cut_based_ct == 1"), t2->GetEntries("cut_based_ct == 1"));
	vector<string> name;
	vector<string> cut;
	name.push_back("tot");
	cut.push_back("");
	name.push_back("1btag");
	cut.push_back("cut_based_ct == 0");
	name.push_back("2btag");
	cut.push_back("cut_based_ct == 1");
	for( unsigned int icat = 0 ; icat < name.size() ; icat++)
	{
		cout << "##### ##### MEAN " << name[icat] << " ##### #####" << endl;
		display("mgg", getMean(t1, "mgg", cut[icat]), getMean(t2, "mgg", cut[icat]));
		display("mjj", getMean(t1, "mjj", cut[icat]), getMean(t2, "mjj", cut[icat]));
		display("mjj_wkinfit", getMean(t1, "mjj_wkinfit", cut[icat]), getMean(t2, "mjj_wkinfit", cut[icat]));
		display("mtot", getMean(t1, "mtot", cut[icat]), getMean(t2, "mtot", cut[icat]));
		display("mtot_wokinfit", getMean(t1, "mtot_wokinfit", cut[icat]), getMean(t2, "mtot_wokinfit", cut[icat]));
		display("cut_based_ct", getMean(t1, "cut_based_ct", cut[icat]), getMean(t2, "cut_based_ct", cut[icat]));
		display("evWeight", getMean(t1, "evWeight", cut[icat]), getMean(t2, "evWeight", cut[icat]));
		display("weightBtagSF", getMean(t1, "weightBtagSF", cut[icat]), getMean(t2, "weightBtagSF", cut[icat]));
		display("weightBtagSFerrUp", getMean(t1, "weightBtagSFerrUp", cut[icat]), getMean(t2, "weightBtagSFerrUp", cut[icat]));
		display("weightBtagSFerrDown", getMean(t1, "weightBtagSFerrDown", cut[icat]), getMean(t2, "weightBtagSFerrDown", cut[icat]));
		cout << "##### ##### RMS " << name[icat] << " ##### #####" << endl;
		display("mgg", getRMS(t1, "mgg", cut[icat]), getRMS(t2, "mgg", cut[icat]));
		display("mjj", getRMS(t1, "mjj", cut[icat]), getRMS(t2, "mjj", cut[icat]));
		display("mjj_wkinfit", getRMS(t1, "mjj_wkinfit", cut[icat]), getRMS(t2, "mjj_wkinfit", cut[icat]));
		display("mtot", getRMS(t1, "mtot", cut[icat]), getRMS(t2, "mtot", cut[icat]));
		display("mtot_wokinfit", getRMS(t1, "mtot_wokinfit", cut[icat]), getRMS(t2, "mtot_wokinfit", cut[icat]));
		display("cut_based_ct", getRMS(t1, "cut_based_ct", cut[icat]), getRMS(t2, "cut_based_ct", cut[icat]));
		display("evWeight", getRMS(t1, "evWeight", cut[icat]), getRMS(t2, "evWeight", cut[icat]));
		display("weightBtagSF", getRMS(t1, "weightBtagSF", cut[icat]), getRMS(t2, "weightBtagSF", cut[icat]));
		display("weightBtagSFerrUp", getRMS(t1, "weightBtagSFerrUp", cut[icat]), getRMS(t2, "weightBtagSFerrUp", cut[icat]));
		display("weightBtagSFerrDown", getRMS(t1, "weightBtagSFerrDown", cut[icat]), getRMS(t2, "weightBtagSFerrDown", cut[icat]));
	}


	return 0;
}
