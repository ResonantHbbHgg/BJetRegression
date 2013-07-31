// ROOT HEADERS
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
// TMVA HEADERS
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
// C++ HEADERS
#include <string>
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

	if( argc == 1 )
	{
		cerr << "WARNING: Arguments should be passed ! Default arguments will be used" << endl;
		cerr << "WARNING: Syntax is " << argv[0] << " -i (inputfile) -it (inputtree) -o (outputfile) -ox (outputxml)" << endl;
	}

	string inputfile = "jetTreeForTraining_m300.root";
	string inputtree = "jets";
	string inputfile2 = "jetTreeForTraining_m500.root";
	string inputtree2 = "jets";
	string outputfile = "regression_test.root";
	string outputxml = "test";

	for(int iarg=0 ; iarg < argc ; iarg++)
	{
		if(strcmp("-i", argv[iarg]) == 0 && argc >= iarg + 1)
			inputfile = argv[iarg+1];
		if(strcmp("-i2", argv[iarg]) == 0 && argc >= iarg + 1)
			inputfile2 = argv[iarg+1];
		if(strcmp("-it", argv[iarg]) == 0 && argc >= iarg + 1)
			inputtree = argv[iarg+1];
		if(strcmp("-it2", argv[iarg]) == 0 && argc >= iarg + 1)
			inputtree2 = argv[iarg+1];
		if(strcmp("-o", argv[iarg]) == 0 && argc >= iarg + 1)
			outputfile = argv[iarg+1];
		if(strcmp("-ox", argv[iarg]) == 0 && argc >= iarg + 1)
			outputxml = argv[iarg+1];
		if(strcmp("--help", argv[iarg]) == 0 || strcmp("-h", argv[iarg]) == 0)
		{
			cerr << "WARNING: Syntax is " << argv[0] << " -i (inputfile) -it (inputtree) -o (outputfile) -ox (outputxml)" << endl;
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
	TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
//	TTree *outtree = new TTree(outputxml.c_str(), Form("%s reduced", outputxml.c_str()));
//  TFile *infile = TFile::Open("jetTreeForTraining.root");
//  TTree *intree  = (TTree*)infile->Get("jets");
  
//  TFile *outfile = new TFile("regressionParton2TMVA.root","RECREATE");
//  TFile *outfile = new TFile("regressionGen2TMVA_globeinputs.root","RECREATE");
  
//  TMVA::Factory* factory = new TMVA::Factory("factoryJetRegParton2",outfile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");
  TMVA::Factory* factory = new TMVA::Factory(outputxml.c_str(),outfile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");

  factory->AddRegressionTree(intree, 1.0);
//  factory->AddRegressionTree(intree2, 100.0);
  
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

  factory->PrepareTrainingAndTestTree(preselectionCut,"nTrain_Regression=10000:nTest_Regression=10000");
//  factory->BookMethod(TMVA::Types::kMLP,"MLP","NCycles=700:HiddenLayers=N,N-1:TestRate=5:TrainingMethod=BFGS:VarTRansform=Norm");
  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=200;nCuts=25"); 
    
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 

	return 0;
}
