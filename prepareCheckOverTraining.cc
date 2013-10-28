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
#include <boost/program_options.hpp>
// DEFINES
#define DEBUG 0
#define USEHT 1
// namespaces
using namespace std;
using namespace TMVA;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	// declare arguments
 	int numberOfRegressionFiles;
 	int nTrainingTrees;
	string inputfile;
	string inputtree;
	float w1;
	string inputfile2;
	string inputtree2;
	float w2;
	string inputfile3;
	string inputtree3;
	float w3;
	string inputfile4;
	string inputtree4;
	float w4;
	string outputfile;
	string regressionFolder;
	int n;
	int j;

 // print out passed arguments
  copy(argv, argv + argc, ostream_iterator<char*>(cout, " ")); cout << endl;
  // argument parsing
  try
  {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("nTrainingTrees", po::value<int>(&nTrainingTrees)->default_value(1), "number of training trees")
      (",n", po::value<int>(&numberOfRegressionFiles)->default_value(2), "number of split for circular training")
      (",j", po::value<int>(&numberOfRegressionFiles)->default_value(0), "index of split for circular checking")
      ("inputfile,i", po::value<string>(&inputfile)->default_value("2013-08-07_jetTreeForTraining_m300.root"), "input file 1")
      ("inputtree,t", po::value<string>(&inputtree)->default_value("jets"), "input tree 1")
      ("inputfile2", po::value<string>(&inputfile2)->default_value("jetTreeForTraining_m500.root"), "input file 2")
      ("inputtree2", po::value<string>(&inputtree2)->default_value("jets"), "input tree 2")
      ("inputfile3", po::value<string>(&inputfile3)->default_value("jetTreeForTraining_m700.root"), "input file 3")
      ("inputtree3", po::value<string>(&inputtree3)->default_value("jets"), "input tree 3")
      ("inputfile4", po::value<string>(&inputfile4)->default_value("jetTreeForTraining_m1000.root"), "input file 4")
      ("inputtree4", po::value<string>(&inputtree4)->default_value("jets"), "input tree 4")
      ("outputfile,o", po::value<string>(&outputfile)->default_value("jetTreeForOverTrainingChecks_m300.root"), "output file")
      ("regressionFolder", po::value<string>(&regressionFolder)->default_value("/afs/cern.ch/user/h/hebda/public/"), "regression folder")
      ("numberOfRegressionFiles,r", po::value<int>(&numberOfRegressionFiles)->default_value(2), "number of regression files")
      ("weight1", po::value<float>(&w1)->default_value(1.), "weight for inputfile 1")
      ("weight2", po::value<float>(&w2)->default_value(1.), "weight for inputfile 2")
      ("weight3", po::value<float>(&w3)->default_value(1.), "weight for inputfile 3")
      ("weight4", po::value<float>(&w4)->default_value(1.), "weight for inputfile 4")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
      cout << desc << "\n";
      return 1;
    }
  } catch(exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  } catch(...) {
    cerr << "Exception of unknown type!\n";
  }
  // end of argument parsing
  //################################################

  if(DEBUG) cout << "End of argument parsing" << endl;
	string regressionfile = regressionFolder + "weights/test_BDT.weights.xml";
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
