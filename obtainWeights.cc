// Get the pt weights to apply the control sample in order to get data-like shapes
// Original code by C. Rovelli (May 2013)
// C++ headers
#include <string>
#include <boost/program_options.hpp>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	// declare arguments
	string inputfile_data;
	string inputfile_CS;
	string inputtree_data;
	string inputtree_CS;
	string outputfile;
	// print out passed arguments
	copy(argv, argv + argc, ostream_iterator<char*>(cout, " ")); cout << endl;
	// argument parsing
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
		("help,h", "produce help message")
		("inputfile_data", po::value<string>(&inputfile_data)->default_value("2014-02-17_selection_noRegression_noMassCut_v10/Data_noRegression_noMassCut_v10.root"), "input data file")
		("inputfile_CS", po::value<string>(&inputfile_CS)->default_value("2014-02-17_selection_noRegression_noMassCut_v10/Data_noRegression_noMassCut_controlSample_v10.root"), "input control sample file")
		("inputtree_data", po::value<string>(&inputtree_data)->default_value("Data"), "input data tree")
		("inputtree_CS", po::value<string>(&inputtree_CS)->default_value("Data"), "input control sample tree")
		("outputfile", po::value<string>(&outputfile)->default_value("scales_2D_pt_data_4GeVbinning.root"), "output file")
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
	// ################################################
	
	TCanvas c1;
	if(DEBUG) cout << "setup files and trees" << endl;
  TFile* dataFile=TFile::Open(inputfile_data.c_str());
  TTree* tree_data=(TTree*)dataFile->Get(inputtree_data.c_str());

	if(DEBUG) cout << "inputfile_CS= " << inputfile_CS << "\tinputtree_CS= " << inputtree_CS << endl;
  TFile* csFile=TFile::Open(inputfile_CS.c_str());
  TTree* tree_cs=(TTree*)csFile->Get(inputtree_CS.c_str());
  
	if(DEBUG) cout << "setup and fill 2d histograms" << endl;
  // pt gamma 
	TH2F *h2D_pt_data = new TH2F("h2D_pt_data", "h2D_pt_data", 35,20.,160.,35,20.,160.);
	TH2F *h2D_pt_data_ub = new TH2F("h2D_pt_data_ub", "h2D_pt_data_ub", 35,20.,160.,35,20.,160.);
	TH2F *h2D_pt_cs = new TH2F("h2D_pt_cs", "h2D_pt_cs", 35,20.,160.,35,20.,160.);
	if(DEBUG) cout << "draw first histo" << endl;
  tree_data->Draw("pho1_pt:pho2_pt>>h2D_pt_data","evweight*(gg_mass>100 && gg_mass<180)*(gg_mass<115 || gg_mass>135)","colz");
  tree_data->Draw("pho1_pt:pho2_pt>>h2D_pt_data_ub","evweight*(gg_mass>100 && gg_mass<180)","colz");
	if(DEBUG) cout << "draw second histo" << endl;
  tree_cs->Draw("pho1_pt:pho2_pt>>h2D_pt_cs","evweight*(gg_mass>100 && gg_mass<180)","colz");
	if(DEBUG) cout << "Now compute integrals" << endl;
  float integral_2D_pt_data = h2D_pt_data->Integral();
  float integral_2D_pt_data_ub = h2D_pt_data_ub->Integral();
	cout << "integral_2D_pt_data= " << integral_2D_pt_data << endl;
	cout << "integral_2D_pt_data_ub= " << integral_2D_pt_data_ub << endl;
//  h2D_pt_data->Scale(1./integral_2D_pt_data);
  float integral_2D_pt_cs = h2D_pt_cs->Integral();
	cout << "integral_2D_pt_cs= " << integral_2D_pt_cs << endl;
//  h2D_pt_cs->Scale(1./integral_2D_pt_cs);
  h2D_pt_cs->Scale(integral_2D_pt_data_ub / integral_2D_pt_cs); // rescale normalization of the CS to match the one of unblinded data
  TH2F* cs_norm_pt = h2D_pt_cs;
  h2D_pt_data->Divide(cs_norm_pt);
  h2D_pt_data->GetXaxis()->SetTitle("p_{T} Sublead Photon");
  h2D_pt_data->GetYaxis()->SetTitle("p_{T} Lead Photon");
  h2D_pt_data->SaveAs(outputfile.c_str());
}
