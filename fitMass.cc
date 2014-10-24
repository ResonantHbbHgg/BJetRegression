// Fit the mjj and mggjj spectra to get the resolution out
// O. Bondu (May 2013)
// C++ headers
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <boost/program_options.hpp>
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TEntryList.h>
#include <TF1.h>
// RooFit headers
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooBernstein.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
// local files
#include "CMSStyle.C"
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;
using namespace RooFit;
namespace po = boost::program_options;

void plotParameters(RooArgList *r2_cat0_param, TCanvas *c, int canvasDivision, RooPlot* frame0, bool isThereSeveralFits, int iclass, string datasetName, float mvaInf, float mvaSup, int precision);
pair<float,float> sigmaEff(TTree* tree, string variable, string cut = "");
float getFWHM(RooAddPdf *pdf, RooRealVar *var);

int main (int argc, char *argv[])
{
	// declare arguments
	int compareWithRegression;
	int massHypothesis;
	int CRYSTALBALL;
	string inputfile;
	string inputtree;
	string inputregfile;
	string inputregtree;
	string outputfile;
	// print out passed arguments
	copy(argv, argv + argc, ostream_iterator<char*>(cout, " ")); cout << endl;
	// argument parsing itself
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message")
			("inputfile", po::value<string>(&inputfile)->default_value("2014-02-17_selection_noRegression_noMassCut_v10/Radion_m300_8TeV_noRegression_noMassCut_v10.root"), "input file")
			("inputtree", po::value<string>(&inputtree)->default_value("Radion_m300_8TeV"), "input tree")
			("inputregfile", po::value<string>(&inputregfile)->default_value("2014-02-17_selection_noRegression_noMassCut_v10/Radion_m300_8TeV_noRegression_noMassCut_v10.root"), "input file")
			("inputregtree", po::value<string>(&inputregtree)->default_value("Radion_m300_8TeV"), "input tree")
			("massHypothesis", po::value<int>(&massHypothesis)->default_value(300), "mass hypothesis")
			("compareWithRegression", po::value<int>(&compareWithRegression)->default_value(0), "switch on/off comparison with regression")
			("outputfile", po::value<string>(&outputfile)->default_value("resolution_Radion_m300_8TeV.txt"), "output file to store fit results")
			("crystalBall", po::value<int>(&CRYSTALBALL)->default_value(1), "perform CrystalBall fit")
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

	// Initialize plot style
	gROOT->Reset();
	gROOT->ProcessLine(".x setTDRStyle.C");
	CMSStyle();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	if( !compareWithRegression )
	{
		// Initialize stuff
		TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
		TFile *infile = TFile::Open(inputfile.c_str());
		TTree *intree = (TTree*)infile->Get(inputtree.c_str());
		ofstream outfile;
		outfile.open(outputfile.c_str());
//		outfile << "Spectrum\tCategory\tMethod\tmu\tsigma\tres" << endl;

		cout << "Performing the fits" << endl;
		// Observables
		float min_jj = 70.;
		float max_jj = 250.;
		float min_ggjj = massHypothesis - 200.;
		float max_ggjj = massHypothesis + 200.;
		RooRealVar jj_mass("jj_mass", "m_{jj}", min_jj, max_jj, "GeV");
		jj_mass.setBins(45);
		RooRealVar ggjj_mass("ggjj_mass", "m_{jj#gamma#gamma}", min_ggjj, max_ggjj, "GeV");
		ggjj_mass.setBins(50);
		RooRealVar njets_kRadionID_and_CSVM("njets_kRadionID_and_CSVM", "njets_kRadionID_and_CSVM", 0, 10);
	
		vector<string> categoryCut;
		categoryCut.clear();
		vector<string> categoryName;
		categoryName.clear();
	
//		categoryCut.push_back(Form("jj_mass > %f && jj_mass < %f && ggjj_mass > %f && ggjj_mass < %f", min_jj, max_jj, min_ggjj, max_ggjj));
//		categoryName.push_back("allcat");
		categoryCut.push_back(Form("njets_kRadionID_and_CSVM < 1.5 && jj_mass > %f && jj_mass < %f && ggjj_mass > %f && ggjj_mass < %f", min_jj, max_jj, min_ggjj, max_ggjj));
		categoryName.push_back("cat1");
		categoryCut.push_back(Form("njets_kRadionID_and_CSVM > 1.5 && jj_mass > %f && jj_mass < %f && ggjj_mass > %f && ggjj_mass < %f", min_jj, max_jj, min_ggjj, max_ggjj));
		categoryName.push_back("cat0");
		
	
		for(int icat = 0 ; icat < (int)categoryCut.size() ; icat++)
		{
			// datatset definition depend on category
			RooDataSet full_dataset("radion", "radion", intree, RooArgList(jj_mass, ggjj_mass, njets_kRadionID_and_CSVM));
			RooDataSet dataset = *((RooDataSet*)full_dataset.reduce(categoryCut[icat].c_str()));	
	
			float mean_jj = dataset.mean(jj_mass);
			float rms_jj = dataset.rmsVar(jj_mass)->getVal();
		
			float mean_ggjj = dataset.mean(ggjj_mass);
			float rms_ggjj = dataset.rmsVar(ggjj_mass)->getVal();
	
			pair<float, float> ms_jj = make_pair(0., 0.);
			pair<float, float> ms_gg = make_pair(0., 0.);
			pair<float, float> ms_ggjj = make_pair(0., 0.);
	
			ms_jj = sigmaEff(intree, "jj_mass", categoryCut[icat]);
			ms_gg = sigmaEff(intree, "gg_mass", categoryCut[icat]);
			ms_ggjj = sigmaEff(intree, "ggjj_mass", categoryCut[icat]);

			cout << "ncor: " << "\tmean= " << ms_jj.first << "\tsigma= " << ms_jj.second << "\tres= " << ms_jj.second / ms_jj.first * 100. << endl;
			cout << "ncor: " << "\tmean= " << ms_ggjj.first << "\tsigma= " << ms_ggjj.second << "\tres= " << ms_ggjj.second / ms_ggjj.first * 100. << endl;

/*		
			outfile << "mjj\t" << setprecision (2) << fixed << categoryName[icat] << "\t" << "sigmaEff" 
				<< "\t" << ms_jj.first << "\t" << ms_jj.second << "\t" << ms_jj.second / ms_jj.first * 100.
				<< endl;
			outfile << "mggjj\t" << setprecision (2) << fixed << categoryName[icat] << "\t" << "sigmaEff" 
				<< "\t" << ms_ggjj.first << "\t" << ms_ggjj.second << "\t" << ms_ggjj.second / ms_ggjj.first * 100.
				<< endl;
*/
			outfile << "mass= " << massHypothesis
				<< "\t" << "var= " << "mgg"
				<< "\t" << "cat= " << categoryName[icat] 
				<< "\t" << "method= " << "sigmaEff"
				<< "\t" << "width= " << ms_gg.second
				<< "\t" << "mean= " << 125.
				<< endl;
				outfile << "mass= " << massHypothesis
				<< "\t" << "var= " << "mjj"
				<< "\t" << "cat= " << categoryName[icat] 
				<< "\t" << "method= " << "sigmaEff"
				<< "\t" << "width= " << ms_jj.second
				<< "\t" << "mean= " << 125.
				<< endl;
				outfile << "mass= " << massHypothesis
				<< "\t" << "var= " << "mtot"
				<< "\t" << "cat= " << categoryName[icat] 
				<< "\t" << "method= " << "sigmaEff"
				<< "\t" << "width= " << ms_ggjj.second
				<< "\t" << "mean= " << massHypothesis
				<< endl;

			if( CRYSTALBALL )
			{
				// Fit is in three steps:
				// - gauss fit of the peak
				// - gauss + pol3 fit of the full range
				// - CB + pol3 fit of the full range
				
				// Define pol3
				RooRealVar a0_jj("a0_jj", "a0_jj", 0.5001, 0., 1.);
				RooRealVar a1_jj("a1_jj", "a1_jj", 0.5001, 0., 1.);
				RooRealVar a2_jj("a2_jj", "a2_jj", 0.5001, 0., 1.);
				RooRealVar a3_jj("a3_jj", "a3_jj", 0.5001, 0., 1.);
				RooBernstein pol3_jj("pol3_jj", "pol3_jj", jj_mass, RooArgList(a0_jj, a1_jj, a2_jj, a3_jj));
				RooRealVar f0_jj("f0_jj", "f0_jj", 0.2, 0.001, .5);
			
			
				// Define gauss and CB (common parameters)
				RooRealVar mu_CrystalBall_jj("mu_CrystalBall_jj", "#mu", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
				RooRealVar sigma_CrystalBall_jj("sigma_CrystalBall_jj", "#sigma", rms_jj, .01*rms_jj, 5*rms_jj, "GeV");
				RooRealVar alpha_CrystalBall_jj("alpha_CrystalBall_jj", "#alpha", 1., 0., 30., "GeV");
				RooRealVar n_CrystalBall_jj("n_CrystalBall_jj", "n", 10., 0., 30., "GeV");
				RooCBShape CrystalBall_jj("CrystalBall_jj", "CrystalBall_jj", jj_mass, mu_CrystalBall_jj, sigma_CrystalBall_jj, alpha_CrystalBall_jj, n_CrystalBall_jj);
				RooGaussian gauss_jj("gauss_jj", "gauss_jj", jj_mass, mu_CrystalBall_jj, sigma_CrystalBall_jj);
			
			
				// Models definitions
				RooAddPdf first_jj("first_jj", "first_jj", pol3_jj, gauss_jj, f0_jj);
				RooAddPdf model_jj("model_jj", "model_jj", pol3_jj, CrystalBall_jj, f0_jj);
			
			
				RooPlot * jj_frame = jj_mass.frame();
				// plot data
				dataset.plotOn(jj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				// first: gaussian fit
				gauss_jj.fitTo(dataset, Save(), Range(mean_jj-1.5*rms_jj, mean_jj+1.5*rms_jj));
				// then, pol3 fit
				pol3_jj.fitTo(dataset, Save());
				// then pol + gauss fit
				first_jj.fitTo(dataset, Save());
			//	first_jj.plotOn(jj_frame, LineColor(kBlue), LineWidth(2));
				// finally, CB + pol fit
				RooFitResult *f = model_jj.fitTo(dataset, Save());
				model_jj.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
				RooArgList p = f->floatParsFinal();
				jj_frame->Draw();
				double res = sigma_CrystalBall_jj.getVal()/mu_CrystalBall_jj.getVal()*100.;
				// printing out plot parameters + creating plots
				plotParameters( &p, c1, 0, jj_frame, true, 1, Form("CrystalBall+Pol3 (%s)", categoryName[icat].c_str()), 98, 99, 2);
				c1->Print(Form("pdf/mjj_m%d_CrystalBall_%s.pdf", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("root/mjj_m%d_CrystalBall_%s.root", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("gif/mjj_m%d_CrystalBall_%s.gif", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("png/mjj_m%d_CrystalBall_%s.png", massHypothesis, categoryName[icat].c_str()));
				c1->Clear();

				outfile << "mass= " << massHypothesis
				<< "\t" << "var= " << "mjj"
				<< "\t" << "cat= " << categoryName[icat] 
				<< "\t" << "method= " << "CrystalBall"
				<< "\t" << "width= " << sigma_CrystalBall_jj.getVal()
				<< "\t" << "mean= " << 125.
				<< endl;

				outfile << "mass= " << massHypothesis
				<< "\t" << "var= " << "mjj"
				<< "\t" << "cat= " << categoryName[icat] 
				<< "\t" << "method= " << "FWHM"
				<< "\t" << "width= " << getFWHM(&model_jj, &jj_mass)
				<< "\t" << "mean= " << 125.
				<< endl;

			// M_jjgg
				// Define pol3
				RooRealVar a0_ggjj("a0_ggjj", "a0_ggjj", 0.5001, 0., 1.);
				RooRealVar a1_ggjj("a1_ggjj", "a1_ggjj", 0.5001, 0., 1.);
				RooRealVar a2_ggjj("a2_ggjj", "a2_ggjj", 0.5001, 0., 1.);
				RooRealVar a3_ggjj("a3_ggjj", "a3_ggjj", 0.5001, 0., 1.);
				RooBernstein pol3_ggjj("pol3_ggjj", "pol3_ggjj", ggjj_mass, RooArgList(a0_ggjj, a1_ggjj, a2_ggjj, a3_ggjj));
				RooRealVar f0_ggjj("f0_ggjj", "f0_ggjj", 0.2, 0.001, .5);
			
				// Define gauss and CB (common parameters)
				RooRealVar mu_CrystalBall_ggjj("mu_CrystalBall_ggjj", "#mu", mean_ggjj, mean_ggjj-rms_ggjj, mean_ggjj+rms_ggjj, "GeV");
				RooRealVar sigma_CrystalBall_ggjj("sigma_CrystalBall_ggjj", "#sigma", rms_ggjj, .01*rms_ggjj, 5*rms_ggjj, "GeV");
				RooRealVar alpha_CrystalBall_ggjj("alpha_CrystalBall_ggjj", "#alpha", 1., 0., 30., "GeV");
				RooRealVar n_CrystalBall_ggjj("n_CrystalBall_ggjj", "n", 20., 0., 30., "GeV");
				RooCBShape CrystalBall_ggjj("CrystalBall_ggjj", "CrystalBall_ggjj", ggjj_mass, mu_CrystalBall_ggjj, sigma_CrystalBall_ggjj, alpha_CrystalBall_ggjj, n_CrystalBall_ggjj);
				RooGaussian gauss_ggjj("gauss_ggjj", "gauss_ggjj", ggjj_mass, mu_CrystalBall_ggjj, sigma_CrystalBall_ggjj);
			
				// Models definitions
				RooAddPdf first_ggjj("first_ggjj", "first_ggjj", pol3_ggjj, gauss_ggjj, f0_ggjj);
				RooAddPdf model_ggjj("model_ggjj", "model_ggjj", pol3_ggjj, CrystalBall_ggjj, f0_ggjj);
			
				RooPlot * ggjj_frame = ggjj_mass.frame();
				// plot data
				dataset.plotOn(ggjj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				// first: gaussian fit
				gauss_ggjj.fitTo(dataset, Save(), Range(mean_ggjj-1.5*rms_ggjj, mean_ggjj+1.5*rms_ggjj));
				// then, pol3 fit
				pol3_ggjj.fitTo(dataset, Save());
				// then pol + gauss fit
				first_ggjj.fitTo(dataset, Save());
			//	first_ggjj.plotOn(ggjj_frame, LineColor(kBlue), LineWidth(2));
				// finally, CB + pol fit
				RooFitResult *fggjj = model_ggjj.fitTo(dataset, Save());
				model_ggjj.plotOn(ggjj_frame, LineColor(kGreen+3), LineWidth(2));
				RooArgList pggjj = fggjj->floatParsFinal();
				ggjj_frame->Draw();
				// printing out plot parameters + creating plots
				plotParameters( &pggjj, c1, 0, ggjj_frame, true, 1, Form("CrystalBall+Pol3 (%s)", categoryName[icat].c_str()), 98, 99, 2);
				c1->Print(Form("pdf/mggjj_m%d_CrystalBall_%s.pdf", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("root/mggjj_m%d_CrystalBall_%s.root", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("gif/mggjj_m%d_CrystalBall_%s.gif", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("png/mggjj_m%d_CrystalBall_%s.png", massHypothesis, categoryName[icat].c_str()));
				c1->Clear();
				double res_ = sigma_CrystalBall_ggjj.getVal()/mu_CrystalBall_ggjj.getVal()*100.;

				outfile << "mass= " << massHypothesis
				<< "\t" << "var= " << "mtot"
				<< "\t" << "cat= " << categoryName[icat] 
				<< "\t" << "method= " << "CrystalBall"
				<< "\t" << "width= " << sigma_CrystalBall_ggjj.getVal()
				<< "\t" << "mean= " << massHypothesis
				<< endl;

				outfile << "mass= " << massHypothesis
				<< "\t" << "var= " << "mtot"
				<< "\t" << "cat= " << categoryName[icat] 
				<< "\t" << "method= " << "FWHM"
				<< "\t" << "width= " << getFWHM(&model_ggjj, &ggjj_mass)
				<< "\t" << "mean= " << massHypothesis
				<< endl;


	
			} // end of if CRYSTALBALL
	
	
		} // end of loop over categories
	
		// clean stuff
		outfile.close();
		delete c1;
		c1 = 0;
		infile->Close();

	} else {
	//	TFile *infile = TFile::Open("simple_genjet.root");
	//	TFile *infile = TFile::Open("simple_parton.root");
	//	TFile *infilereg = TFile::Open("simple_reg_genjet.root");
	//	TFile *infilereg = TFile::Open("simple_reg_parton.root");
	//	TFile *infile = TFile::Open("simple_genjet_globeinputs.root");
	//	TTree *intree = (TTree*)infile->Get("Radion_m300_8TeV_nm");
	//	TFile *infilereg = TFile::Open("simple_reg_genjet_globeinputs.root");
	//	TTree *intreereg = (TTree*)infilereg->Get("Radion_m300_8TeV_nm");
		TFile *infile = TFile::Open(inputfile.c_str());
		TTree *intree = (TTree*)infile->Get(inputtree.c_str());
		TFile *infilereg = TFile::Open(inputregfile.c_str());
		TTree *intreereg = (TTree*)infilereg->Get(inputregtree.c_str());
//		ofstream outfile;
//		outfile.open(outputfile.c_str());
//		TFile *infile = TFile::Open("2013-12-12_selection_noRegression_noMassCut_v02/Radion_m300_8TeV_noRegression_noMassCut_v02.root");
//		TTree *intree = (TTree*)infile->Get("Radion_m300_8TeV");
//		TFile *infilereg = TFile::Open("2013-12-12_selection_noRegression_noMassCut_v02/Radion_m350_8TeV_noRegression_noMassCut_v02.root");
//		TTree *intreereg = (TTree*)infilereg->Get("Radion_m350_8TeV");
		ofstream outfile_mjj, outfile_mggjj;
		outfile_mjj.open(Form("performanceSummary_m%d_mjj.txt", massHypothesis));
		outfile_mggjj.open(Form("performanceSummary_m%d_mggjj.txt", massHypothesis));
		outfile_mjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
		outfile_mggjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
		
		// Observables
		float min_jj = 70.;
		float max_jj = 250.;
		float min_ggjj = massHypothesis - 200.;
		float max_ggjj = massHypothesis + 200.;
		RooRealVar jj_mass("jj_mass", "m_{jj}", min_jj, max_jj, "GeV");
		jj_mass.setBins(45);
		RooRealVar ggjj_mass("ggjj_mass", "m_{jj#gamma#gamma}", min_ggjj, max_ggjj, "GeV");
		ggjj_mass.setBins(50);
		RooRealVar njets_kRadionID_and_CSVM("njets_kRadionID_and_CSVM", "njets_kRadionID_and_CSVM", 0, 10);
		TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	
		bool GAUSS = false;
		bool CRYSTALBALL = true;
		bool VOIGT = false;
		bool SIMVOIGT = false;
	
		vector<string> categoryCut;
		categoryCut.clear();
		vector<string> categoryName;
		categoryName.clear();
	
		categoryCut.push_back(Form("jj_mass > %f && jj_mass < %f && ggjj_mass > %f && ggjj_mass < %f", min_jj, max_jj, min_ggjj, max_ggjj));
		categoryName.push_back("allcat");
		categoryCut.push_back(Form("njets_kRadionID_and_CSVM < 1.5 && jj_mass > %f && jj_mass < %f && ggjj_mass > %f && ggjj_mass < %f", min_jj, max_jj, min_ggjj, max_ggjj));
		categoryName.push_back("cat1");
		categoryCut.push_back(Form("njets_kRadionID_and_CSVM > 1.5 && jj_mass > %f && jj_mass < %f && ggjj_mass > %f && ggjj_mass < %f", min_jj, max_jj, min_ggjj, max_ggjj));
		categoryName.push_back("cat0");
		
	
		for(int icat = 0 ; icat < (int)categoryCut.size() ; icat++)
		{
			// datatset definition depend on category
			RooDataSet full_dataset("radion", "radion", intree, RooArgList(jj_mass, ggjj_mass, njets_kRadionID_and_CSVM));
			RooDataSet full_datasetreg("radion", "radion", intreereg, RooArgList(jj_mass, ggjj_mass, njets_kRadionID_and_CSVM));
			RooDataSet dataset = *((RooDataSet*)full_dataset.reduce(categoryCut[icat].c_str()));	
			RooDataSet datasetreg = *((RooDataSet*)full_datasetreg.reduce(categoryCut[icat].c_str()));	
	
			float mean_jj = dataset.mean(jj_mass);
			float mean_regjj = datasetreg.mean(jj_mass);
			float rms_jj = dataset.rmsVar(jj_mass)->getVal();
			float rms_regjj = datasetreg.rmsVar(jj_mass)->getVal();
		
			float mean_ggjj = dataset.mean(ggjj_mass);
			float mean_regggjj = datasetreg.mean(ggjj_mass);
			float rms_ggjj = dataset.rmsVar(ggjj_mass)->getVal();
			float rms_regggjj = datasetreg.rmsVar(ggjj_mass)->getVal();
	
			pair<float, float> ms_jj = make_pair(0., 0.);
			pair<float, float> ms_regjj = make_pair(0., 0.);
			pair<float, float> ms_ggjj = make_pair(0., 0.);
			pair<float, float> ms_regggjj = make_pair(0., 0.);
	
			ms_jj = sigmaEff(intree, "jj_mass", categoryCut[icat]);
			ms_regjj = sigmaEff(intreereg, "jj_mass", categoryCut[icat]);
			ms_ggjj = sigmaEff(intree, "ggjj_mass", categoryCut[icat]);
			ms_regggjj = sigmaEff(intreereg, "ggjj_mass", categoryCut[icat]);
	
			cout << "ncor: " << "\tmean= " << ms_jj.first << "\tsigma= " << ms_jj.second << "\tres= " << ms_jj.second / ms_jj.first * 100. << endl;
			cout << "cor: " << "\tmean= " << ms_regjj.first << "\tsigma= " << ms_regjj.second << "\tres= " << ms_regjj.second / ms_regjj.first * 100. << endl;
			cout << "improvement= " << - (ms_regjj.second / ms_regjj.first - ms_jj.second / ms_jj.first) / (ms_jj.second / ms_jj.first) * 100. << endl;
		
	//		outfile_mjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
			outfile_mjj << setprecision (2) << fixed << categoryName[icat] << "\t" << "sigmaEff" 
				<< "\t" << ms_jj.first << "\t" << ms_jj.second << "\t" << ms_jj.second / ms_jj.first * 100.
				<< "\t" << ms_regjj.first << "\t" << ms_regjj.second << "\t" << ms_regjj.second / ms_regjj.first * 100.
				<< "\t" << - (ms_regjj.second / ms_regjj.first - ms_jj.second / ms_jj.first) / (ms_jj.second / ms_jj.first) * 100.
				<< endl;
	
			outfile_mggjj << setprecision (2) << fixed << categoryName[icat] << "\t" << "sigmaEff" 
				<< "\t" << ms_ggjj.first << "\t" << ms_ggjj.second << "\t" << ms_ggjj.second / ms_ggjj.first * 100.
				<< "\t" << ms_regggjj.first << "\t" << ms_regggjj.second << "\t" << ms_regggjj.second / ms_regggjj.first * 100.
				<< "\t" << - (ms_regggjj.second / ms_regggjj.first - ms_ggjj.second / ms_ggjj.first) / (ms_ggjj.second / ms_ggjj.first) * 100.
				<< endl;
	
			cout << "sigma improvement= " << (ms_jj.second - ms_regjj.second) / ms_jj.second * 100. << endl; 
				// fit parameters
			if(GAUSS)
			{
				RooRealVar mu_gauss("mu_gauss", "mean", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
				RooRealVar sigma_gauss("sigma_gauss", "sigma", rms_jj, .01*rms_jj, 5.*rms_jj, "GeV");
				RooGaussian gauss("gauss", "gauss", jj_mass, mu_gauss, sigma_gauss);
			
				RooRealVar mu_gaussreg("mu_gaussreg", "mean (reg)", mean_regjj, mean_regjj-rms_regjj, mean_regjj+rms_regjj, "GeV");
				RooRealVar sigma_gaussreg("sigma_gaussreg", "sigma (reg)", rms_regjj, .01*rms_regjj, 5.*rms_regjj, "GeV");
				RooGaussian gaussreg("gaussreg", "gaussreg", jj_mass, mu_gaussreg, sigma_gaussreg);
			
				RooRealVar mu_gauss_("mu_gauss_", "mean", mean_ggjj, mean_ggjj-rms_ggjj, mean_ggjj+rms_ggjj, "GeV");
				RooRealVar sigma_gauss_("sigma_gauss_", "sigma", rms_ggjj, .01*rms_ggjj, 5.*rms_ggjj, "GeV");
				RooGaussian gauss_("gauss_", "gauss_", ggjj_mass, mu_gauss_, sigma_gauss_);
			
				RooRealVar mu_gauss_reg("mu_gauss_reg", "mean (reg)", mean_regggjj, mean_regggjj-rms_regggjj, mean_regggjj+rms_regggjj, "GeV");
				RooRealVar sigma_gauss_reg("sigma_gauss_reg", "sigma (reg)", rms_regggjj, .01*rms_regggjj, 5.*rms_regggjj, "GeV");
				RooGaussian gauss_reg("gauss_reg", "gauss_reg", ggjj_mass, mu_gauss_reg, sigma_gauss_reg);
			
				RooPlot * jj_frame = jj_mass.frame();
				dataset.plotOn(jj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
			//	RooFitResult *f = gauss.fitTo(dataset, Save(), Range(mean_jj-0.9*rms_jj, mean_jj+1.5*rms_jj));
				RooFitResult *f = gauss.fitTo(dataset, Save(), Range(mean_jj-1.5*rms_jj, mean_jj+1.5*rms_jj));
				gauss.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
				RooFitResult * freg = gaussreg.fitTo(datasetreg, Save(), Range(mean_regjj-1.5*rms_regjj, mean_regjj+1.5*rms_regjj));
				gaussreg.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList p = f->floatParsFinal();
				RooArgList preg = freg->floatParsFinal();	
				jj_frame->Draw();
				float res = sigma_gauss.getVal() / mu_gauss.getVal() * 100.;
				float resreg = sigma_gaussreg.getVal() / mu_gaussreg.getVal() * 100.;
				plotParameters( &p, c1, 0, jj_frame, true, 1, "Gaussian", 98, 99, 2);
				plotParameters( &preg, c1, 0, jj_frame, true, 3, "test", res, resreg, 2);
				c1->Print(Form("pdf/mjj_gauss_%s.pdf", categoryName[icat].c_str()));
				c1->Print(Form("root/mjj_gauss_%s.root", categoryName[icat].c_str()));
				c1->Print(Form("gif/mjj_gauss_%s.gif", categoryName[icat].c_str()));
				c1->Print(Form("png/mjj_gauss_%s.png", categoryName[icat].c_str()));
				c1->Clear();
	//		outfile_mjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
				outfile_mjj << setprecision (2) << fixed << categoryName[icat] << "\tGauss"
					<< "\t" << mu_gauss.getVal() << "\t" << sigma_gauss.getVal() << "\t" << sigma_gauss.getVal() / mu_gauss.getVal() * 100.
					<< "\t" << mu_gaussreg.getVal() << "\t" << sigma_gaussreg.getVal() << "\t" << sigma_gaussreg.getVal() / mu_gaussreg.getVal() * 100.
					<< "\t" << - (resreg - res) / res * 100.
					<< endl;
			
				RooPlot * ggjj_frame = ggjj_mass.frame();
				dataset.plotOn(ggjj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(ggjj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
				RooFitResult *f_ = gauss_.fitTo(dataset, Save(), Range(mean_ggjj-1.5*rms_ggjj, mean_ggjj+1.5*rms_ggjj));
				gauss_.plotOn(ggjj_frame, LineColor(kGreen+3), LineWidth(2));
				RooFitResult * freg_ = gauss_reg.fitTo(datasetreg, Save(), Range(mean_regggjj-1.5*rms_regggjj, mean_regggjj+1.5*rms_regggjj));
				gauss_reg.plotOn(ggjj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList p_ = f_->floatParsFinal();
				RooArgList preg_ = freg_->floatParsFinal();	
				ggjj_frame->Draw();
				float res_ = sigma_gauss_.getVal() / mu_gauss_.getVal() * 100.;
				float resreg_ = sigma_gauss_reg.getVal() / mu_gauss_reg.getVal() * 100.;
				plotParameters( &p_, c1, 0, ggjj_frame, true, 1, "Gaussian", 98, 99, 2);
				plotParameters( &preg_, c1, 0, ggjj_frame, true, 3, "test", res_, resreg_, 2);
				c1->Print(Form("pdf/mggjj_gauss_%s.pdf", categoryName[icat].c_str()));
				c1->Print(Form("root/mggjj_gauss_%s.root", categoryName[icat].c_str()));
				c1->Print(Form("gif/mggjj_gauss_%s.gif", categoryName[icat].c_str()));
				c1->Print(Form("png/mggjj_gauss_%s.png", categoryName[icat].c_str()));
				c1->Clear();
			
				cout << "res= " << res << "\tresreg= " << resreg << "\timprov= " << fabs(res-resreg)/res*100. << endl;
				cout << "res_= " << res_ << "\tresreg_= " << resreg_ << "\timprov= " << fabs(res_-resreg_)/res_*100. << endl;
	//		outfile_mjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
				outfile_mggjj << setprecision (2) << fixed << categoryName[icat] << "\tGauss"
					<< "\t" << mu_gauss_.getVal() << "\t" << sigma_gauss_.getVal() << "\t" << sigma_gauss_.getVal() / mu_gauss_.getVal() * 100.
					<< "\t" << mu_gauss_reg.getVal() << "\t" << sigma_gauss_reg.getVal() << "\t" << sigma_gauss_reg.getVal() / mu_gauss_reg.getVal() * 100.
					<< "\t" << - (resreg_ - res_) / res_ * 100.
					<< endl;
				} // end of if GAUSS
	
			if( CRYSTALBALL )
			{
				// Fit is in three steps:
				// - gauss fit of the peak
				// - gauss + pol3 fit of the full range
				// - CB + pol3 fit of the full range
				
				// Define pol3
				RooRealVar a0_jj("a0_jj", "a0_jj", 0.5001, 0., 1.);
				RooRealVar a1_jj("a1_jj", "a1_jj", 0.5001, 0., 1.);
				RooRealVar a2_jj("a2_jj", "a2_jj", 0.5001, 0., 1.);
				RooRealVar a3_jj("a3_jj", "a3_jj", 0.5001, 0., 1.);
				RooBernstein pol3_jj("pol3_jj", "pol3_jj", jj_mass, RooArgList(a0_jj, a1_jj, a2_jj, a3_jj));
				RooRealVar f0_jj("f0_jj", "f0_jj", 0.2, 0.001, .5);
			
				RooRealVar a0_regjj("a0_regjj", "a0_regjj", 0.5001, 0., 1.);
				RooRealVar a1_regjj("a1_regjj", "a1_regjj", 0.5001, 0., 1.);
				RooRealVar a2_regjj("a2_regjj", "a2_regjj", 0.5001, 0., 1.);
				RooRealVar a3_regjj("a3_regjj", "a3_regjj", 0.5001, 0., 1.);
				RooBernstein pol3_regjj("pol3_regjj", "pol3_regjj", jj_mass, RooArgList(a0_regjj, a1_regjj, a2_regjj, a3_regjj));
				RooRealVar f0_regjj("f0_regjj", "f0_regjj", 0.2, 0.001, .5);
			
				// Define gauss and CB (common parameters)
				RooRealVar mu_CrystalBall_jj("mu_CrystalBall_jj", "#mu", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
				RooRealVar sigma_CrystalBall_jj("sigma_CrystalBall_jj", "#sigma", rms_jj, .01*rms_jj, 5*rms_jj, "GeV");
				RooRealVar alpha_CrystalBall_jj("alpha_CrystalBall_jj", "#alpha", 1., 0., 30., "GeV");
				RooRealVar n_CrystalBall_jj("n_CrystalBall_jj", "n", 10., 0., 30., "GeV");
				RooCBShape CrystalBall_jj("CrystalBall_jj", "CrystalBall_jj", jj_mass, mu_CrystalBall_jj, sigma_CrystalBall_jj, alpha_CrystalBall_jj, n_CrystalBall_jj);
				RooGaussian gauss_jj("gauss_jj", "gauss_jj", jj_mass, mu_CrystalBall_jj, sigma_CrystalBall_jj);
			
				RooRealVar mu_CrystalBall_regjj("mu_CrystalBall_regjj", "#mu (reg)", mean_regjj, mean_regjj-rms_regjj, mean_regjj+rms_regjj, "GeV");
				RooRealVar sigma_CrystalBall_regjj("sigma_CrystalBall_regjj", "#sigma (reg)", rms_regjj, .01*rms_regjj, 5*rms_regjj, "GeV");
				RooRealVar alpha_CrystalBall_regjj("alpha_CrystalBall_regjj", "#alpha (reg)", 1., 0., 30., "GeV");
				RooRealVar n_CrystalBall_regjj("n_CrystalBall_regjj", "n (reg)", 10., 0., 30., "GeV");
				RooCBShape CrystalBall_regjj("CrystalBall_regjj", "CrystalBall_regjj", jj_mass, mu_CrystalBall_regjj, sigma_CrystalBall_regjj, alpha_CrystalBall_regjj, n_CrystalBall_regjj);
				RooGaussian gauss_regjj("gauss_regjj", "gauss_regjj", jj_mass, mu_CrystalBall_regjj, sigma_CrystalBall_regjj);
			
				// Models definitions
				RooAddPdf first_jj("first_jj", "first_jj", pol3_jj, gauss_jj, f0_jj);
				RooAddPdf model_jj("model_jj", "model_jj", pol3_jj, CrystalBall_jj, f0_jj);
			
				RooAddPdf first_regjj("first_regjj", "first_regjj", pol3_regjj, gauss_regjj, f0_regjj);
				RooAddPdf model_regjj("model_regjj", "model_regjj", pol3_regjj, CrystalBall_regjj, f0_regjj);
			
				RooPlot * jj_frame = jj_mass.frame();
				// plot data
				dataset.plotOn(jj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
				// first: gaussian fit
				gauss_jj.fitTo(dataset, Save(), Range(mean_jj-1.5*rms_jj, mean_jj+1.5*rms_jj));
				gauss_regjj.fitTo(dataset, Save(), Range(mean_regjj-1.5*rms_regjj, mean_regjj+1.5*rms_regjj));
				// then, pol3 fit
				pol3_jj.fitTo(dataset, Save());
				pol3_regjj.fitTo(dataset, Save());
				// then pol + gauss fit
				first_jj.fitTo(dataset, Save());
			//	first_jj.plotOn(jj_frame, LineColor(kBlue), LineWidth(2));
				first_regjj.fitTo(datasetreg, Save());
			//	first_regjj.plotOn(jj_frame, LineColor(kBlue), LineWidth(2));
				// finally, CB + pol fit
				RooFitResult *f = model_jj.fitTo(dataset, Save());
				RooFitResult *freg = model_regjj.fitTo(datasetreg, Save());
				model_jj.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
				model_regjj.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList p = f->floatParsFinal();
				RooArgList preg = freg->floatParsFinal();	
				jj_frame->Draw();
				double res = sigma_CrystalBall_jj.getVal()/mu_CrystalBall_jj.getVal()*100.;
				double resreg = sigma_CrystalBall_regjj.getVal()/mu_CrystalBall_regjj.getVal()*100.;
				// printing out plot parameters + creating plots
				plotParameters( &p, c1, 0, jj_frame, true, 1, "CrystalBall+Pol3", 98, 99, 2);
				plotParameters( &preg, c1, 0, jj_frame, true, 3, "test", sigma_CrystalBall_jj.getVal()/mu_CrystalBall_jj.getVal()*100., sigma_CrystalBall_regjj.getVal()/mu_CrystalBall_regjj.getVal()*100., 2);
				c1->Print(Form("pdf/mjj_m%d_CrystalBall_%s.pdf", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("root/mjj_m%d_CrystalBall_%s.root", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("gif/mjj_m%d_CrystalBall_%s.gif", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("png/mjj_m%d_CrystalBall_%s.png", massHypothesis, categoryName[icat].c_str()));
				c1->Clear();
		//		outfile_mjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
				outfile_mjj << setprecision (2) << fixed << categoryName[icat] << "\tCrystalBall"
					<< "\t" << mu_CrystalBall_jj.getVal() << "\t" << sigma_CrystalBall_jj.getVal() << "\t" << sigma_CrystalBall_jj.getVal() / mu_CrystalBall_jj.getVal() * 100.
					<< "\t" << mu_CrystalBall_regjj.getVal() << "\t" << sigma_CrystalBall_regjj.getVal() << "\t" << sigma_CrystalBall_regjj.getVal() / mu_CrystalBall_regjj.getVal() * 100.
					<< "\t" << - (resreg - res) / res * 100.
					<< endl;
			
			// M_jjgg
				// Define pol3
				RooRealVar a0_ggjj("a0_ggjj", "a0_ggjj", 0.5001, 0., 1.);
				RooRealVar a1_ggjj("a1_ggjj", "a1_ggjj", 0.5001, 0., 1.);
				RooRealVar a2_ggjj("a2_ggjj", "a2_ggjj", 0.5001, 0., 1.);
				RooRealVar a3_ggjj("a3_ggjj", "a3_ggjj", 0.5001, 0., 1.);
				RooBernstein pol3_ggjj("pol3_ggjj", "pol3_ggjj", ggjj_mass, RooArgList(a0_ggjj, a1_ggjj, a2_ggjj, a3_ggjj));
				RooRealVar f0_ggjj("f0_ggjj", "f0_ggjj", 0.2, 0.001, .5);
			
				RooRealVar a0_regggjj("a0_regggjj", "a0_regggjj", 0.5001, 0., 1.);
				RooRealVar a1_regggjj("a1_regggjj", "a1_regggjj", 0.5001, 0., 1.);
				RooRealVar a2_regggjj("a2_regggjj", "a2_regggjj", 0.5001, 0., 1.);
				RooRealVar a3_regggjj("a3_regggjj", "a3_regggjj", 0.5001, 0., 1.);
				RooBernstein pol3_regggjj("pol3_regggjj", "pol3_regggjj", ggjj_mass, RooArgList(a0_regggjj, a1_regggjj, a2_regggjj, a3_regggjj));
				RooRealVar f0_regggjj("f0_regggjj", "f0_regggjj", 0.2, 0.001, .5);
			
				// Define gauss and CB (common parameters)
				RooRealVar mu_CrystalBall_ggjj("mu_CrystalBall_ggjj", "#mu", mean_ggjj, mean_ggjj-rms_ggjj, mean_ggjj+rms_ggjj, "GeV");
				RooRealVar sigma_CrystalBall_ggjj("sigma_CrystalBall_ggjj", "#sigma", rms_ggjj, .01*rms_ggjj, 5*rms_ggjj, "GeV");
				RooRealVar alpha_CrystalBall_ggjj("alpha_CrystalBall_ggjj", "#alpha", 1., 0., 30., "GeV");
				RooRealVar n_CrystalBall_ggjj("n_CrystalBall_ggjj", "n", 20., 0., 30., "GeV");
				RooCBShape CrystalBall_ggjj("CrystalBall_ggjj", "CrystalBall_ggjj", ggjj_mass, mu_CrystalBall_ggjj, sigma_CrystalBall_ggjj, alpha_CrystalBall_ggjj, n_CrystalBall_ggjj);
				RooGaussian gauss_ggjj("gauss_ggjj", "gauss_ggjj", ggjj_mass, mu_CrystalBall_ggjj, sigma_CrystalBall_ggjj);
			
				RooRealVar mu_CrystalBall_regggjj("mu_CrystalBall_regggjj", "#mu (reg)", mean_regggjj, mean_regggjj-rms_regggjj, mean_regggjj+rms_regggjj, "GeV");
				RooRealVar sigma_CrystalBall_regggjj("sigma_CrystalBall_regggjj", "#sigma (reg)", rms_regggjj, .01*rms_regggjj, 5*rms_regggjj, "GeV");
				RooRealVar alpha_CrystalBall_regggjj("alpha_CrystalBall_regggjj", "#alpha (reg)", 1., 0., 30., "GeV");
				RooRealVar n_CrystalBall_regggjj("n_CrystalBall_regggjj", "n (reg)", 20., 0., 30., "GeV");
				RooCBShape CrystalBall_regggjj("CrystalBall_regggjj", "CrystalBall_regggjj", ggjj_mass, mu_CrystalBall_regggjj, sigma_CrystalBall_regggjj, alpha_CrystalBall_regggjj, n_CrystalBall_regggjj);
				RooGaussian gauss_regggjj("gauss_regggjj", "gauss_regggjj", ggjj_mass, mu_CrystalBall_regggjj, sigma_CrystalBall_regggjj);
			
				// Models definitions
				RooAddPdf first_ggjj("first_ggjj", "first_ggjj", pol3_ggjj, gauss_ggjj, f0_ggjj);
				RooAddPdf model_ggjj("model_ggjj", "model_ggjj", pol3_ggjj, CrystalBall_ggjj, f0_ggjj);
			
				RooAddPdf first_regggjj("first_regggjj", "first_regggjj", pol3_regggjj, gauss_regggjj, f0_regggjj);
				RooAddPdf model_regggjj("model_regggjj", "model_regggjj", pol3_regggjj, CrystalBall_regggjj, f0_regggjj);
			
				RooPlot * ggjj_frame = ggjj_mass.frame();
				// plot data
				dataset.plotOn(ggjj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(ggjj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
				// first: gaussian fit
				gauss_ggjj.fitTo(dataset, Save(), Range(mean_ggjj-1.5*rms_ggjj, mean_ggjj+1.5*rms_ggjj));
				gauss_regggjj.fitTo(dataset, Save(), Range(mean_regggjj-1.5*rms_regggjj, mean_regggjj+1.5*rms_regggjj));
				// then, pol3 fit
				pol3_ggjj.fitTo(dataset, Save());
				pol3_regggjj.fitTo(dataset, Save());
				// then pol + gauss fit
				first_ggjj.fitTo(dataset, Save());
			//	first_ggjj.plotOn(ggjj_frame, LineColor(kBlue), LineWidth(2));
				first_regggjj.fitTo(datasetreg, Save());
			//	first_regggjj.plotOn(ggjj_frame, LineColor(kBlue), LineWidth(2));
				// finally, CB + pol fit
				RooFitResult *fggjj = model_ggjj.fitTo(dataset, Save());
				RooFitResult *fregggjj = model_regggjj.fitTo(datasetreg, Save());
				model_ggjj.plotOn(ggjj_frame, LineColor(kGreen+3), LineWidth(2));
				model_regggjj.plotOn(ggjj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList pggjj = fggjj->floatParsFinal();
				RooArgList pregggjj = fregggjj->floatParsFinal();	
				ggjj_frame->Draw();
				// printing out plot parameters + creating plots
				plotParameters( &pggjj, c1, 0, ggjj_frame, true, 1, "CrystalBall+Pol3", 98, 99, 2);
				plotParameters( &pregggjj, c1, 0, ggjj_frame, true, 3, "test", sigma_CrystalBall_ggjj.getVal()/mu_CrystalBall_ggjj.getVal()*100., sigma_CrystalBall_regggjj.getVal()/mu_CrystalBall_regggjj.getVal()*100., 2);
				c1->Print(Form("pdf/mggjj_m%d_CrystalBall_%s.pdf", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("root/mggjj_m%d_CrystalBall_%s.root", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("gif/mggjj_m%d_CrystalBall_%s.gif", massHypothesis, categoryName[icat].c_str()));
				c1->Print(Form("png/mggjj_m%d_CrystalBall_%s.png", massHypothesis, categoryName[icat].c_str()));
				c1->Clear();
				double res_ = sigma_CrystalBall_ggjj.getVal()/mu_CrystalBall_ggjj.getVal()*100.;
				double resreg_ = sigma_CrystalBall_regggjj.getVal()/mu_CrystalBall_regggjj.getVal()*100.;
				//		outfile_mjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
				outfile_mggjj << setprecision (2) << fixed << categoryName[icat] << "\tCrystalBall"
					<< "\t" << mu_CrystalBall_ggjj.getVal() << "\t" << sigma_CrystalBall_ggjj.getVal() << "\t" << sigma_CrystalBall_ggjj.getVal() / mu_CrystalBall_ggjj.getVal() * 100.
					<< "\t" << mu_CrystalBall_regggjj.getVal() << "\t" << sigma_CrystalBall_regggjj.getVal() << "\t" << sigma_CrystalBall_regggjj.getVal() / mu_CrystalBall_regggjj.getVal() * 100.
					<< "\t" << - (resreg_ - res_) / res_ * 100.
					<< endl;
		
			} // end of if CRYSTALBALL
	
			if(VOIGT)
			{
				RooRealVar mu_voigt("mu_voigt", "mean", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
				RooRealVar width_voigt("width_voigt", "width", rms_jj, 0.1*rms_jj, 5.*rms_jj, "GeV");
				RooRealVar sigma_voigt("sigma_voigt", "sigma", 5., 0.0001, 20., "GeV");
				RooVoigtian voigt("voigt", "voigt", jj_mass, mu_voigt, width_voigt, sigma_voigt);
			
				RooRealVar mu_voigtreg("mu_voigtreg", "mean (reg)", mean_regjj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
				RooRealVar width_voigtreg("width_voigtreg", "width (reg)", rms_jj, 0.1*rms_jj, 5.*rms_jj, "GeV");
				RooRealVar sigma_voigtreg("sigma_voigtreg", "sigma (reg)", 5., 0.0001, 20., "GeV");
				RooVoigtian voigtreg("voigtreg", "voigtreg", jj_mass, mu_voigtreg, width_voigtreg, sigma_voigtreg);
				RooVoigtian voigtreg2("voigtreg2", "voigtreg2", jj_mass, mu_voigtreg, width_voigtreg, sigma_voigt);
			
				RooPlot * jj_frame = jj_mass.frame();
				dataset.plotOn(jj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
				RooFitResult *f = voigt.fitTo(dataset, Save());
				voigt.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
				RooFitResult * freg = voigtreg.fitTo(datasetreg, Save());
				voigtreg.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList p = f->floatParsFinal();
				RooArgList preg = freg->floatParsFinal();	
				jj_frame->Draw();
			  float fg = 2. * sigma_voigt.getVal() * sqrt(2. * log(2));
			  float fl = 2. * width_voigt.getVal();
			  float fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
			  float fgreg = 2. * sigma_voigtreg.getVal() * sqrt(2. * log(2));
			  float flreg = 2. * width_voigtreg.getVal();
			  float fvreg = 0.5346 * flreg + sqrt(0.2166 * pow(flreg, 2.) + pow(fgreg, 2.));
			  float res=  fv / mu_voigt.getVal() * 100. / (2. * sqrt(2. * log(2.)));
			  float resreg=  fvreg / mu_voigtreg.getVal() * 100. / (2. * sqrt(2. * log(2.)));
				plotParameters( &p, c1, 0, jj_frame, true, 1, "Voigtian", 98, 99, 2);
				plotParameters( &preg, c1, 0, jj_frame, true, 4, "test", res, resreg, 2);
				c1->Print(Form("pdf/mjj_voigt_%s.pdf", categoryName[icat].c_str()));
				c1->Print(Form("root/mjj_voigt_%s.root", categoryName[icat].c_str()));
				c1->Print(Form("gif/mjj_voigt_%s.gif", categoryName[icat].c_str()));
				c1->Print(Form("png/mjj_voigt_%s.png", categoryName[icat].c_str()));
				c1->Clear();
				//		outfile_mjj << "Category\tMethod\tmu\tsigma\tres\tmu_reg\tsigma_reg\tres_reg\timprovement" << endl;
				outfile_mjj << setprecision (2) << fixed << categoryName[icat] << "\tVoigt"
					<< "\t" << mu_voigt.getVal() << "\t" << fv << "\t" << res
					<< "\t" << mu_voigtreg.getVal() << "\t" << fvreg << "\t" << resreg
					<< "\t" << - (resreg - res) / res * 100.
					<< endl;
				
				RooRealVar mu_voigt_("mu_voigt_", "mean", massHypothesis, min_ggjj, max_ggjj, "GeV");
				RooRealVar width_voigt_("width_voigt_", "width", 35., 5., 50., "GeV");
				RooRealVar sigma_voigt_("sigma_voigt_", "sigma", 5., .0, 20., "GeV");
				RooVoigtian voigt_("voigt_", "voigt_", ggjj_mass, mu_voigt_, width_voigt_, sigma_voigt_);
			
				RooRealVar mu_voigt_reg("mu_voigt_reg", "mean (reg)", massHypothesis, min_ggjj, max_ggjj, "GeV");
				RooRealVar width_voigt_reg("width_voigt_reg", "width (reg)", 35., 5., 50., "GeV");
				RooRealVar sigma_voigt_reg("sigma_voigt_reg", "sigma (reg)", 5., 0.0001, 20., "GeV");
				RooVoigtian voigt_reg("voigt_reg", "voigt_reg", ggjj_mass, mu_voigt_reg, width_voigt_reg, sigma_voigt_reg);
			
				RooPlot * ggjj_frame = ggjj_mass.frame();
				dataset.plotOn(ggjj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(ggjj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
				RooFitResult *f_ = voigt_.fitTo(dataset, Save());
				voigt_.plotOn(ggjj_frame, LineColor(kGreen+3), LineWidth(2));
				RooFitResult * freg_ = voigt_reg.fitTo(datasetreg, Save());
				voigt_reg.plotOn(ggjj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList p_ = f_->floatParsFinal();
				RooArgList preg_ = freg_->floatParsFinal();	
				ggjj_frame->Draw();
				float fg_ = 2. * sigma_voigt_.getVal() * sqrt(2. * log(2));
				float fl_ = 2. * width_voigt_.getVal();
				float fv_ = 0.5346 * fl_ + sqrt(0.2166 * pow(fl_, 2.) + pow(fg_, 2.));
				float fgreg_ = 2. * sigma_voigt_reg.getVal() * sqrt(2. * log(2));
				float flreg_ = 2. * width_voigt_reg.getVal();
				float fvreg_ = 0.5346 * flreg_ + sqrt(0.2166 * pow(flreg_, 2.) + pow(fgreg_, 2.));
				float res_=  fv_ / mu_voigt_.getVal() * 100. / (2. * sqrt(2. * log(2.)));
				float resreg_=  fvreg_ / mu_voigt_reg.getVal() * 100. / (2. * sqrt(2. * log(2.)));
				plotParameters( &p_, c1, 0, ggjj_frame, true, 1, "Voigtian", 98, 99, 2);
				plotParameters( &preg_, c1, 0, ggjj_frame, true, 4, "test", res_, resreg_, 2);
				c1->Print(Form("pdf/mggjj_voigt_%s.pdf", categoryName[icat].c_str()));
				c1->Print(Form("root/mggjj_voigt_%s.root", categoryName[icat].c_str()));
				c1->Print(Form("gif/mggjj_voigt_%s.gif", categoryName[icat].c_str()));
				c1->Print(Form("png/mggjj_voigt_%s.png", categoryName[icat].c_str()));
				c1->Clear();
				outfile_mggjj << setprecision (2) << fixed << categoryName[icat] << "\tVoigt"
					<< "\t" << mu_voigt_.getVal() << "\t" << fv_ << "\t" << res_
					<< "\t" << mu_voigt_reg.getVal() << "\t" << fvreg_ << "\t" << resreg_
					<< "\t" << - (resreg_ - res_) / res_ * 100.
					<< endl;
				
			  cout << "fg= " << fg << "\tfl= " << fl << "\tfv= " << fv << endl;
			  cout << "fgreg= " << fgreg << "\tflreg= " << flreg << "\tfvreg= " << fvreg << endl;
			  cout << "res= " << res << "\tresreg= " << resreg << "\timprov= " << fabs(res - resreg)/res * 100.<< endl;
				cout << "fg_= " << fg_ << "\tfl_= " << fl_ << "\tfv_= " << fv_ << endl;
				cout << "fgreg_= " << fgreg_ << "\tflreg_= " << flreg_ << "\tfvreg_= " << fvreg_ << endl;
				cout << "res_= " << res_ << "\tresreg_= " << resreg_ << "\timprov_= " << fabs(res_ - resreg_)/res_ * 100.<< endl;
			}// end of if VOIGT
			
			if(SIMVOIGT)
			{
				RooRealVar mu_voigt("mu_voigt", "mean", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
				RooRealVar width_voigt("width_voigt", "width", rms_jj, 0.1*rms_jj, 5.*rms_jj, "GeV");
				RooRealVar sigma_voigt("sigma_voigt", "sigma", 5., 0.0001, 20., "GeV");
				RooVoigtian voigt("voigt", "voigt", jj_mass, mu_voigt, width_voigt, sigma_voigt);
			
				RooRealVar mu_voigtreg("mu_voigtreg", "mean (reg)", mean_regjj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
				RooRealVar width_voigtreg("width_voigtreg", "width (reg)", rms_jj, 0.1*rms_jj, 5.*rms_jj, "GeV");
				RooRealVar sigma_voigtreg("sigma_voigtreg", "sigma (reg)", 5., 0.0001, 20., "GeV");
				RooVoigtian voigtreg2("voigtreg2", "voigtreg2", jj_mass, mu_voigtreg, width_voigtreg, sigma_voigt);
			
				RooCategory sample("sample", "sample");
				sample.defineType("vanilla");
				sample.defineType("regression");
			
				RooDataSet combData("combdata", "combdata", RooArgList(jj_mass, ggjj_mass), Index(sample), Import("vanilla", dataset), Import("regression", datasetreg));
				RooSimultaneous sim("sim", "sim", sample);
				sim.addPdf(voigt, "vanilla");
				sim.addPdf(voigtreg2, "regression");
			
				RooPlot * jj_frame = jj_mass.frame();
				dataset.plotOn(jj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
				RooFitResult *f = sim.fitTo(combData, Save());
				voigt.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
				voigtreg2.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList p = f->floatParsFinal();
				jj_frame->Draw();
			  float fg = 2. * sigma_voigt.getVal() * sqrt(2. * log(2));
			  float fl = 2. * width_voigt.getVal();
			  float fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
			  float fgreg = 2. * sigma_voigt.getVal() * sqrt(2. * log(2));
			  float flreg = 2. * width_voigtreg.getVal();
			  float fvreg = 0.5346 * flreg + sqrt(0.2166 * pow(flreg, 2.) + pow(fgreg, 2.));
			  float res=  fv / mu_voigt.getVal() * 100. / (2. * sqrt(2. * log(2.)));
			  float resreg=  fvreg / mu_voigtreg.getVal() * 100. / (2. * sqrt(2. * log(2.)));
				plotParameters( &p, c1, 0, jj_frame, false, 1, "Simultaneous Voigtian fit", res, resreg, 2);
				c1->Print(Form("pdf/mjj_simvoigt_%s.pdf", categoryName[icat].c_str()));
				c1->Print(Form("root/mjj_simvoigt_%s.root", categoryName[icat].c_str()));
				c1->Print(Form("gif/mjj_simvoigt_%s.gif", categoryName[icat].c_str()));
				c1->Print(Form("png/mjj_simvoigt_%s.png", categoryName[icat].c_str()));
				c1->Clear();
				outfile_mjj << setprecision (2) << fixed << categoryName[icat] << "\tSimVoigt"
					<< "\t" << mu_voigt.getVal() << "\t" << fv << "\t" << res
					<< "\t" << mu_voigtreg.getVal() << "\t" << fvreg << "\t" << resreg
					<< "\t" << - (resreg - res) / res * 100.
					<< endl;
			
				RooRealVar mu_voigt_("mu_voigt_", "mean", massHypothesis, min_ggjj, max_ggjj, "GeV");
				RooRealVar width_voigt_("width_voigt_", "width", 35., 5., 50., "GeV");
				RooRealVar sigma_voigt_("sigma_voigt_", "sigma", 5., .0, 20., "GeV");
				RooVoigtian voigt_("voigt_", "voigt_", ggjj_mass, mu_voigt_, width_voigt_, sigma_voigt_);
			
				RooRealVar mu_voigt_reg("mu_voigt_reg", "mean (reg)", massHypothesis, min_ggjj, max_ggjj, "GeV");
				RooRealVar width_voigt_reg("width_voigt_reg", "width (reg)", 35., 5., 50., "GeV");
				RooRealVar sigma_voigt_reg("sigma_voigt_reg", "sigma (reg)", 5., 0.0001, 20., "GeV");
				RooVoigtian voigt_reg2("voigt_reg2", "voigt_reg2", ggjj_mass, mu_voigt_reg, width_voigt_reg, sigma_voigt_);
			
				RooSimultaneous sim_("sim_", "sim_", sample);
				sim_.addPdf(voigt_, "vanilla");
				sim_.addPdf(voigt_reg2, "regression");
			
				RooPlot * ggjj_frame = ggjj_mass.frame();
				dataset.plotOn(ggjj_frame, LineColor(kGreen+3), MarkerColor(kGreen+3), XErrorSize(0));
				datasetreg.plotOn(ggjj_frame, LineColor(kRed+2), MarkerColor(kRed+2), MarkerStyle(23), XErrorSize(0));
				RooFitResult *f_ = sim_.fitTo(combData, Save());
				voigt_.plotOn(ggjj_frame, LineColor(kGreen+3), LineWidth(2));
				voigt_reg2.plotOn(ggjj_frame, LineColor(kRed+2), LineWidth(2));
				RooArgList p_ = f_->floatParsFinal();
				ggjj_frame->Draw();
			  float fg_ = 2. * sigma_voigt_.getVal() * sqrt(2. * log(2));
			  float fl_ = 2. * width_voigt_.getVal();
			  float fv_ = 0.5346 * fl_ + sqrt(0.2166 * pow(fl_, 2.) + pow(fg_, 2.));
			  float fg_reg = 2. * sigma_voigt_.getVal() * sqrt(2. * log(2));
			  float fl_reg = 2. * width_voigt_reg.getVal();
			  float fv_reg = 0.5346 * fl_reg + sqrt(0.2166 * pow(fl_reg, 2.) + pow(fg_reg, 2.));
			  float res_=  fv_ / mu_voigt_.getVal() * 100. / (2. * sqrt(2. * log(2.)));
			  float resreg_=  fv_reg / mu_voigt_reg.getVal() * 100. / (2. * sqrt(2. * log(2.)));
				plotParameters( &p_, c1, 0, ggjj_frame, false, 1, "Simultaneous Voigtian fit", res_, resreg_, 2);
				c1->Print(Form("pdf/mggjj_simvoigt_%s.pdf", categoryName[icat].c_str()));
				c1->Print(Form("root/mggjj_simvoigt_%s.root", categoryName[icat].c_str()));
				c1->Print(Form("gif/mggjj_simvoigt_%s.gif", categoryName[icat].c_str()));
				c1->Print(Form("png/mggjj_simvoigt_%s.png", categoryName[icat].c_str()));
				c1->Clear();
				outfile_mggjj << setprecision (2) << fixed << categoryName[icat] << "\tSimVoigt"
					<< "\t" << mu_voigt_.getVal() << "\t" << fv_ << "\t" << res_
					<< "\t" << mu_voigt_reg.getVal() << "\t" << fv_reg << "\t" << resreg_
					<< "\t" << - (resreg_ - res_) / res_ * 100.
					<< endl;
				} // end of if SIMVOIGT
			
		}// end of loop over categories
	
		outfile_mjj.close();
		outfile_mggjj.close();
	
		delete c1;
		c1 = 0;
		infile->Close();
		infilereg->Close();
	} // end of if compareWithRegression
	return 0;
} // end of main

void plotParameters(RooArgList *r2_cat0_param, TCanvas *c, int canvasDivision, RooPlot* frame0, bool isThereSeveralFits, int iclass, string fitfunctionName, float mvaInf, float mvaSup, int precision )
{
//	cout << "entering plotParameters" << endl;
	// Written by H. Brun (November 2011)
	c->cd(canvasDivision);
	TLatex latexLabel;
	latexLabel.SetNDC();
  latexLabel.SetTextSize(0.03);
	std::ostringstream mvaStream;
	mvaStream << setprecision (3) << fixed << mvaInf << " #leq mva < " << mvaSup;
	string mvaString = mvaStream.str();
//  latexLabel.DrawLatex(0.18, 0.96, Form("CMS Private 2011, #sqrt{s} = 7 TeV, Simulation: %s, %s", isThereSeveralFits ? "signal":"background", mvaString.c_str()));
  latexLabel.DrawLatex(0.18, 0.96, "CMS Private 2013, #sqrt{s} = 8 TeV");
//  latexLabel.DrawLatex(0.45, 0.91, "#sqrt{s} = 7 TeV");
//  latexLabel.DrawLatex(0.67, 0.91, isThereSeveralFits ? "Simulation: signal" : "Simulation: background");
  latexLabel.SetTextSize(0.03);
  RooRealVar* obj = new RooRealVar();
  float position = 0.92;
  position -= 0.04 * (iclass);
  if(iclass == 1) latexLabel.DrawLatex(0.55, position, Form("Fit: %s", fitfunctionName.c_str()));
  TIterator *it = (TIterator*) r2_cat0_param->createIterator();

	obj = (RooRealVar*)it->Next();
	while( obj != 0 )
  {
//		cout << "obj->GetName()= " << obj->GetName() << "\tobj->getVal()= " << obj->getVal() << "\tobj->getError()= " << obj->getError() << endl;
   if( strspn("f0_", (char*)obj->GetName())>=3 ) {obj = (RooRealVar*)it->Next(); continue;}
   if( strspn("a0_", (char*)obj->GetName())>=3 ) {obj = (RooRealVar*)it->Next(); continue;}
   if( strspn("a1_", (char*)obj->GetName())>=3 ) {obj = (RooRealVar*)it->Next(); continue;}
   if( strspn("a2_", (char*)obj->GetName())>=3 ) {obj = (RooRealVar*)it->Next(); continue;}
   if( strspn("a3_", (char*)obj->GetName())>=3 ) {obj = (RooRealVar*)it->Next(); continue;}
   if( strspn("n_", (char*)obj->GetName())>=2 ) {obj = (RooRealVar*)it->Next(); continue;}
   if( strspn("alpha_", (char*)obj->GetName())>=6 ) {obj = (RooRealVar*)it->Next(); continue;}
    position -= 0.04;
    std::ostringstream valueStream;
    if( (float)obj->getError() != 0.0 )
    {
      valueStream << setprecision (precision) << fixed << (float)obj->getVal() << " +- " << (float)obj->getError();
//      cout << setprecision (precision) << fixed << (float)obj->getVal() << " +- " << (float)obj->getError() << endl;;
    } else {
       valueStream << setprecision (precision) << fixed << (float)obj->getVal();
    }
    string valueString = valueStream.str();
    latexLabel.DrawLatex(0.60, position, Form("%s = %s %s", obj->GetTitle(), valueString.c_str(), (char*)obj->getUnit()));
  	cout << "obj->GetName()= " << obj->GetName() << "\tobj->getVal()= " << ((RooRealVar*)obj)->getVal() << endl;
		obj = (RooRealVar*)it->Next();
  }
if( (iclass != 1) || (!isThereSeveralFits) )
{
	position -= 0.04;
  std::ostringstream valueStream;
	valueStream << setprecision (precision) << fixed << mvaInf;
	string valueString = valueStream.str();
	latexLabel.DrawLatex(0.60, position, Form("#sigma/#mu= %s %%", valueString.c_str()));
	position -= 0.04;
	valueStream.clear(); valueStream.str(""); valueString = "";
	valueStream << setprecision (precision) << fixed << mvaSup;
	valueString = valueStream.str();
	latexLabel.DrawLatex(0.60, position, Form("#sigma/#mu (reg)= %s %%", valueString.c_str()));
	position -= 0.04;
	valueStream.clear(); valueStream.str(""); valueString = "";
	valueStream << setprecision (precision) << fixed << fabs(mvaInf - mvaSup) / mvaInf * 100.;
	valueString = valueStream.str();
	latexLabel.DrawLatex(0.60, position, Form("improvement= %s %%", valueString.c_str()));
}

	return;
} // end of plotParameters

pair<float,float> sigmaEff(TTree* tree, string var, string cut)
{
	float sigmaEff_ = 0.;
	float mean = 0.;
	int ntot = tree->GetEntries(cut.c_str());
	float variable = 0.;
	tree->SetBranchAddress(var.c_str(), &variable);

  tree->Draw(">>elist", cut.c_str(), "entrylist");
//  tree->Draw(">>elist", "jj_mass < 300", "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");

//	cout << "ntot= " << ntot << endl;
//	cout << "elist->GetN()= " << elist->GetN() << endl;

//	tree->SetEntryList(elist);
	int ievt = elist->GetEntry(0);
	vector<float> values;
	for(int ievtlist = 0 ; ievtlist < ntot ; ievtlist++)
	{
		if(ievt < 0 ) continue;
//		cout << "ievtlist= " << ievtlist << "\tievt= " << ievt << endl;
		tree->GetEntry(ievt);
		mean += variable / (float)ntot;
		values.push_back(variable);
		ievt = elist->Next();
	}// end of event loop

	// sorting values
	sort(values.begin(), values.end());
	// number of entries to be considered
	// 1 sigma = 68.2689492%
	int nInSigma = floor(elist->GetN() * 0.682689492);
//	cout << "nInSigma= " << nInSigma << endl;
	sigmaEff_ = values[elist->GetN()-1] - values[0];
//	cout << "sigmaEff_= " << sigmaEff_ << endl;
	float valueMin = values[0];
	float valueMax = values[elist->GetN()-1];
	float sigma = valueMax - valueMin;
//	cout << "valueMin= " << valueMin << "\tvalueMax= " << valueMax << endl;
	for(int iscan = 0 ; iscan < (elist->GetN() - nInSigma) ; iscan++)
	{
		sigma = values[iscan+nInSigma] - values[iscan];
//		cout << "iscan= " << iscan << "\tvalues[iscan]= " << values[iscan] << "\tvalues[iscan+nInSigma]= " << values[iscan+nInSigma] << "\tsigma= " << sigma << "\tsigmaMax= " << sigmaMax << endl;
		if( sigma < sigmaEff_ )
		{
			valueMin = values[iscan];
			valueMax = values[iscan+nInSigma];
			sigmaEff_ = sigma;
		}
	}
//	sigmaEff_ = sigmaEff_;
	tree->ResetBranchAddresses();

//	cout << "elist->GetN()= " << elist->GetN() << "\tvalues.size()= " << values.size() << endl;
//	cout << var << "\tcut= " << cut << "\tntot= " << ntot << "\tmean= " << mean << endl;
	return make_pair(mean, sigmaEff_);
} // end of sigmaEff

float getFWHM(RooAddPdf *pdf, RooRealVar *var)
{
	TF1 *f1 = (TF1*)pdf->asTF(RooArgList(*var));
//	TCanvas *c = new TCanvas();
//	f1->Draw();
//	c->Print("plaf.pdf");
	float max = f1->GetMaximum();
	float xmax = f1->GetMaximumX();
	float x_hig = f1->GetX(max / 2., xmax, var->getMax());
	float x_low = f1->GetX(max / 2., var->getMin(), xmax);
//	cout << "FWHM" << "\tvar.getMin()= " << var->getMin() << "\tvar.getMax()= " << var->getMax() << endl;
//	cout << "FWHM" << "\tmax= " << max << "\txmax= " << xmax << "\tx_hig= " << x_hig << "\tx_low= " << x_low << endl;
	return x_hig-x_low;
}
