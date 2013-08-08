// ROOT HEADERS
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TClonesArray.h"
// TMVA HEADERS
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
// C++ HEADERS
#include <string>
#include <iostream>
#include <sstream>
// DEFINES
#define DEBUG 0
// namespaces
using namespace std;
using namespace TMVA;

int main(int argc, char *argv[])
{
	if(DEBUG) cout << "DEBUG: Initialisation: reading parameters" << endl;
	cout << "argc= " << argc << endl;
	for(int iarg = 0 ; iarg < argc; iarg++)
		cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
		string syntax = Form("WARNING: Syntax is %s -i (inputfile) -it (inputtree) -o (outputfile) -ox (outputxml) -n (numberOfSplit) -j (treeSplit)", argv[0]);

	if( argc == 1 )
	{
		cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
		cerr << syntax << endl;
	}

	int nTrainingTrees = 1;
	string inputfile = "jetTreeForTraining_m300.root";
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
	string outputfile = "regression_test.root";
	string outputxml = "test";
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
		if(strcmp("-ox", argv[iarg]) == 0 && argc >= iarg + 1)
			outputxml = argv[iarg+1];
		if(strcmp("-ntt", argv[iarg]) == 0 && argc >= iarg + 1)
			{ std::stringstream ss ( argv[iarg+1] ); ss >> nTrainingTrees; }
		if(strcmp("--help", argv[iarg]) == 0 || strcmp("-h", argv[iarg]) == 0)
		{
			cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
			cerr << syntax << endl;
			cerr << "inputfile= " << inputfile << endl;
			cerr << "inputtree= " << inputtree << endl;
			cerr << "outputfile= " << outputfile << endl;
			cerr << "outputxml= " << outputxml << endl;
			return 2;
		}
	}

	cout << "inputfile= " << inputfile << endl;
	cout << "inputtree= " << inputtree << endl;
	cout << "outputfile= " << outputfile << endl;
	cout << "outputxml= " << outputxml << endl;

	TFile *infile = TFile::Open(inputfile.c_str());
	TTree *intree = (TTree*)infile->Get(inputtree.c_str());
	TFile *infile2 = TFile::Open(inputfile2.c_str());
	TTree *intree2 = (TTree*)infile2->Get(inputtree2.c_str());
	TFile *infile3 = TFile::Open(inputfile3.c_str());
	TTree *intree3 = (TTree*)infile3->Get(inputtree3.c_str());
	TFile *infile4 = TFile::Open(inputfile4.c_str());
	TTree *intree4 = (TTree*)infile4->Get(inputtree4.c_str());
	TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");


//	TTree *outtree = new TTree(outputxml.c_str(), Form("%s reduced", outputxml.c_str()));
//  TFile *infile = TFile::Open("jetTreeForTraining.root");
//  TTree *intree  = (TTree*)infile->Get("jets");
  
//  TFile *outfile = new TFile("regressionParton2TMVA.root","RECREATE");
//  TFile *outfile = new TFile("regressionGen2TMVA_globeinputs.root","RECREATE");
  
//  TMVA::Factory* factory = new TMVA::Factory("factoryJetRegParton2",outfile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");
  TMVA::Factory* factory = new TMVA::Factory(outputxml.c_str(),outfile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");

	int nTrain = 0;
	int nTest = 0;
	int nEval = 0;

	if( n == 1 )
	{
  	factory->AddRegressionTree(intree, w1);
  	if( nTrainingTrees > 1) factory->AddRegressionTree(intree2, w2);
  	if( nTrainingTrees > 2) factory->AddRegressionTree(intree3, w3);
  	if( nTrainingTrees > 3) factory->AddRegressionTree(intree4, w4);
	} else {
		for(int i = 0 ; i < n ; i++)
		{
			if(i == 0)
			{
				factory->AddRegressionTree( intree->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w1, TMVA::Types::kTraining );
				nTrain = intree->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n))->GetEntries();
				cout << "nTrain= " << nTrain << endl; 
				if( nTrainingTrees > 1)
					factory->AddRegressionTree( intree2->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w2, TMVA::Types::kTraining );
				if( nTrainingTrees > 2)
					factory->AddRegressionTree( intree3->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w3, TMVA::Types::kTraining );
				if( nTrainingTrees > 3)
					factory->AddRegressionTree( intree4->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w4, TMVA::Types::kTraining );
			} else if(i == 1) {
				factory->AddRegressionTree( intree->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w1, TMVA::Types::kTesting );
				nTest = intree->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n))->GetEntries();
				cout << "nTest= " << nTest << endl;
				if( nTrainingTrees > 1)
					factory->AddRegressionTree( intree2->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w2, TMVA::Types::kTesting );
				if( nTrainingTrees > 2)
					factory->AddRegressionTree( intree3->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w3, TMVA::Types::kTesting );
				if( nTrainingTrees > 3)
					factory->AddRegressionTree( intree4->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w4, TMVA::Types::kTesting );
/*			} else {
				factory->AddRegressionTree( intree->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w1, TMVA::Types::kValidation );
				cout << "nEval= " << intree->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n))->GetEntries() << endl;
				if( nTrainingTrees > 1)
					factory->AddRegressionTree( intree2->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w2, TMVA::Types::kValidation );
				if( nTrainingTrees > 2)
					factory->AddRegressionTree( intree3->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w3, TMVA::Types::kValidation );
				if( nTrainingTrees > 3)
					factory->AddRegressionTree( intree4->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n)), w4, TMVA::Types::kValidation );
*/
			}
		}
	}
//		intree->CopyTree(Form("(event %% %i == (%i %% %i)) && jet_genDR<0.4 && jet_csvBtag > 0.", n, j+i, n));


  factory->AddVariable("jet_pt"						, "p_{T}^{j}", "GeV",'F');
  factory->AddVariable("jet_eta"					, "#eta^{j}", "",'F');
	factory->AddVariable("jet_emfrac"				, "#epsilon_{EM}^{j}", "", 'F');
	factory->AddVariable("jet_nConstituents", "n_{const}^{j}", "", 'I');
	factory->AddVariable("jet_hadfrac"			, "#epsilon_{HAD}^{j}", "", 'F');
  factory->AddVariable("jet_secVtxPt"			, "2^{nd}vtx_{p_{T}}", "", 'F');
  factory->AddVariable("jet_secVtx3dL"		, "2^{nd}vtx_{3dL}", "", 'F');
  factory->AddVariable("ev_met_corr_pfmet", "MET", "", 'F');
  factory->AddVariable("jet_dPhiMet"			, "#Delta #phi(j, MET)", "",'F');
//	factory->AddSpectator("ev_weight", 'F');
//	factory->SetWeightExpression("ev_weight");
//  factory->AddVariable("jet_csvBtag"   ,"CSV output", "",'F');
//  factory->AddVariable("ev_rho"       , "#rho", "GeV",'F');

//  factory->AddVariable("jet_Chadfrac"    ,'F');
//  factory->AddVariable("jet_Phofrac"    ,'F');
//  factory->AddVariable("jet_Nhadfrac"    ,'F');
//  factory->AddVariable("jet_Elefrac"    ,'F');
//  factory->AddVariable("jet_Mufrac"    ,'F');
//  factory->AddVariable("jet_ptD"    ,'F');
//  factory->AddVariable("jet_secVtx3deL",'F');
  //factory->AddVariable("jetE"      ,'F');
  
//  factory->AddTarget("jet_prtPt");
  factory->AddTarget("jet_genPt");
    
//  TCut preselectionCut("jet_prtDR<0.6 && jet_csvBtag > 0.");
  TCut preselectionCut("jet_genDR<0.4 && jet_csvBtag > 0.");
	unsigned int nentries = intree->GetEntries(preselectionCut);
	cout << "nentries= " << nentries << endl;
	cout << "nTrain + nTest = " << nTrain + nTest << endl;

	if(DEBUG) cout << "prepare training" << endl;
//  factory->PrepareTrainingAndTestTree(preselectionCut,"nTrain_Regression=10000:nTest_Regression=10000");
  factory->PrepareTrainingAndTestTree(preselectionCut,"SplitMode=Block:nTrain_Regression=0:nTest_Regression=0");
	if(DEBUG) cout << "book method" << endl;
//  factory->BookMethod(TMVA::Types::kMLP,"MLP","NCycles=700:HiddenLayers=N,N-1:TestRate=5:TrainingMethod=BFGS:VarTRansform=Norm");
//  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=200:nCuts=25"); // default
//  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=200:nCuts=-1:PruneStrength=-1:PruneMethod=CostComplexity"); 
//  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=1000:nCuts=25:MaxDepth=4"); // TMVA manual, page 110: Currently it looks as if in TMVA, better results for the whole forest are often achieved when pruning is not applied, but rather the maximal tree depth is set to a relatively small value (3 or 4) already during the tree building phase.
  factory->BookMethod(TMVA::Types::kBDT,"BDT",Form("NTrees=%i:nCuts=25:MaxDepth=4", 500)); 
  if(DEBUG) cout << "train" << endl; 
  factory->TrainAllMethods();
  if(DEBUG) cout << "test" << endl; 
  factory->TestAllMethods();
  if(DEBUG) cout << "evaluate" << endl; 
  factory->EvaluateAllMethods(); 

	return 0;
}
