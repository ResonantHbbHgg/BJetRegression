#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
using namespace TMVA;

int main()
{
  TFile *inf = TFile::Open("test.root");
  TTree *tr  = (TTree*)inf->Get("jets");
  
//  TFile *outf = new TFile("regressionParton2TMVA.root","RECREATE");
  TFile *outf = new TFile("regressionGen2TMVA.root","RECREATE");
  
//  TMVA::Factory* factory = new TMVA::Factory("factoryJetRegParton2",outf,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");
  TMVA::Factory* factory = new TMVA::Factory("factoryJetRegGen2",outf,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");

  factory->AddRegressionTree(tr);
  
  factory->AddVariable("jet_csvBtag"   ,'F');
  factory->AddVariable("jet_pt"     ,'F');
  factory->AddVariable("jet_eta"    ,'F');
  factory->AddVariable("jet_dPhiMet" ,'F');
  factory->AddVariable("jet_Chadfrac"    ,'F');
  factory->AddVariable("jet_Phofrac"    ,'F');
  factory->AddVariable("jet_Nhadfrac"    ,'F');
  factory->AddVariable("jet_Elefrac"    ,'F');
  factory->AddVariable("jet_Mufrac"    ,'F');
  factory->AddVariable("jet_ptD"    ,'F');
  factory->AddVariable("jet_secVtxPt"  ,'F');
  factory->AddVariable("jet_secVtx3dL" ,'F');
  factory->AddVariable("jet_secVtx3deL",'F');
  factory->AddVariable("ev_met_corrMet"       ,'F');
  factory->AddVariable("ev_rho"       ,'F');
  //factory->AddVariable("jetE"      ,'F');
  
//  factory->AddTarget("jet_prtPt");
  factory->AddTarget("jet_genPt");
    
//  TCut preselectionCut("jet_prtDR<0.25");
  TCut preselectionCut("jet_genDR<0.25");

  factory->PrepareTrainingAndTestTree(preselectionCut,"nTrain_Regression=14000:nTest_Regression=14000");
  factory->BookMethod(TMVA::Types::kMLP,"MLP","NCycles=700:HiddenLayers=N,N-1:TestRate=5:TrainingMethod=BFGS:VarTRansform=Norm");
  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=200;nCuts=25"); 
    
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 

	return 0;
}
