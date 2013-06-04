// C++ headers
#include <string>
#include <sstream>
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
// RooFit headers
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooPlot.h"
#include "RooFitResult.h"
// local files
#include "CMSStyle.C"
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;
using namespace RooFit;

void plotParameters(RooArgList *r2_cat0_param, TCanvas *c, int canvasDivision, RooPlot* frame0, bool isSignal, int iclass, string datasetName, double mvaInf, double mvaSup, int precision);

int main ()
{
	gROOT->Reset();
	gROOT->ProcessLine(".x setTDRStyle.C");
	CMSstyle();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	TFile *infile = TFile::Open("simple.root");
	TTree *intree = (TTree*)infile->Get("Radion_m300_8TeV_nm");
	TFile *infilereg = TFile::Open("simple_reg.root");
	TTree *intreereg = (TTree*)infilereg->Get("Radion_m300_8TeV_nm");
	
	// Observables
	RooRealVar jj_mass("jj_mass", "m_{jj}", 50., 200., "GeV");
	jj_mass.setBins(50);
	RooRealVar ggjj_mass("ggjj_mass", "m_{jj#gamma#gamma}", 200., 400., "GeV");
	ggjj_mass.setBins(50);

	// fit parameters
	RooRealVar mu_voigt("mu_voigt", "mean", 120., 100., 150., "GeV");
	RooRealVar width_voigt("width_voigt", "width", 35., 5., 50., "GeV");
	RooRealVar sigma_voigt("sigma_voigt", "sigma", 5., .0, 20., "GeV");
	RooVoigtian voigt("voigt", "voigt", jj_mass, mu_voigt, width_voigt, sigma_voigt);

	RooRealVar mu_voigtreg("mu_voigtreg", "mean (reg)", 125., 50., 200., "GeV");
	RooRealVar width_voigtreg("width_voigtreg", "width (reg)", 35., 5., 50., "GeV");
	RooRealVar sigma_voigtreg("sigma_voigtreg", "sigma (reg)", 5., 0.0001, 20., "GeV");
	RooVoigtian voigtreg("voigtreg", "voigtreg", jj_mass, mu_voigtreg, width_voigtreg, sigma_voigtreg);

	RooRealVar mu_voigt_("mu_voigt_", "mean", 300., 200., 400., "GeV");
	RooRealVar width_voigt_("width_voigt_", "width", 35., 5., 50., "GeV");
	RooRealVar sigma_voigt_("sigma_voigt_", "sigma", 5., .0, 20., "GeV");
	RooVoigtian voigt_("voigt_", "voigt_", ggjj_mass, mu_voigt_, width_voigt_, sigma_voigt_);

	RooRealVar mu_voigt_reg("mu_voigt_reg", "mean (reg)", 300., 200., 400., "GeV");
	RooRealVar width_voigt_reg("width_voigt_reg", "width (reg)", 35., 5., 50., "GeV");
	RooRealVar sigma_voigt_reg("sigma_voigt_reg", "sigma (reg)", 5., 0.0001, 20., "GeV");
	RooVoigtian voigt_reg("voigt_reg", "voigt_reg", ggjj_mass, mu_voigt_reg, width_voigt_reg, sigma_voigt_reg);

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	RooDataSet dataset("radion", "radion", intree, RooArgList(jj_mass, ggjj_mass));
	RooDataSet datasetreg("radion", "radion", intreereg, RooArgList(jj_mass, ggjj_mass));

	RooPlot * jj_frame = jj_mass.frame();
	dataset.plotOn(jj_frame, LineColor(kGreen+3));
	datasetreg.plotOn(jj_frame, LineColor(kRed+2));
	RooFitResult *f = voigt.fitTo(dataset, Save());
	voigt.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
	RooFitResult * freg = voigtreg.fitTo(datasetreg, Save());
	voigtreg.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
	RooArgList p = f->floatParsFinal();
	RooArgList preg = freg->floatParsFinal();	
	jj_frame->Draw();
	plotParameters( &p, c1, 0, jj_frame, true, 1, "test", 98, 99, 2);
	plotParameters( &preg, c1, 0, jj_frame, true, 4, "test", 98, 99, 4);
	c1->Print("mjj.pdf");
	c1->Clear();

	RooPlot * ggjj_frame = ggjj_mass.frame();
	dataset.plotOn(ggjj_frame, LineColor(kGreen+3));
	datasetreg.plotOn(ggjj_frame, LineColor(kRed+2));
	RooFitResult *f_ = voigt_.fitTo(dataset, Save());
	voigt_.plotOn(ggjj_frame, LineColor(kGreen+3), LineWidth(2));
	RooFitResult * freg_ = voigt_reg.fitTo(datasetreg, Save());
	voigt_reg.plotOn(ggjj_frame, LineColor(kRed+2), LineWidth(2));
	RooArgList p_ = f_->floatParsFinal();
	RooArgList preg_ = freg_->floatParsFinal();	
	ggjj_frame->Draw();
	plotParameters( &p_, c1, 0, ggjj_frame, true, 1, "test", 98, 99, 2);
	plotParameters( &preg_, c1, 0, ggjj_frame, true, 4, "test", 98, 99, 2);
	c1->Print("mggjj.pdf");
	c1->Clear();

	

	delete c1;
	c1 = 0;
//	jj_frame->Delete();
	infile->Close();
	infilereg->Close();

	return 0;
}

void plotParameters(RooArgList *r2_cat0_param, TCanvas *c, int canvasDivision, RooPlot* frame0, bool isSignal, int iclass, string datasetName, double mvaInf, double mvaSup, int precision )
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
//  latexLabel.DrawLatex(0.18, 0.96, Form("CMS Private 2011, #sqrt{s} = 7 TeV, Simulation: %s, %s", isSignal ? "signal":"background", mvaString.c_str()));
  latexLabel.DrawLatex(0.18, 0.96, "CMS Private 2013, #sqrt{s} = 8 TeV");
//  latexLabel.DrawLatex(0.45, 0.91, "#sqrt{s} = 7 TeV");
//  latexLabel.DrawLatex(0.67, 0.91, isSignal ? "Simulation: signal" : "Simulation: background");
  latexLabel.SetTextSize(0.03);
  RooRealVar* obj = new RooRealVar();
  double position = 0.92;
  position -= 0.04 * (iclass);
  if(iclass == 1) latexLabel.DrawLatex(0.55, position,"Fit: Voigtian");
  TIterator *it = (TIterator*) r2_cat0_param->createIterator();

	obj = (RooRealVar*)it->Next();
	while( obj != 0 )
  {
//		cout << "obj->GetName()= " << obj->GetName() << "\tobj->getVal()= " << obj->getVal() << "\tobj->getError()= " << obj->getError() << endl;
//   if( ! strcmp(((char*)obj->GetName()), "mu_voigt") ) continue;
    position -= 0.04;
    std::ostringstream valueStream;
    if( (double)obj->getError() != 0.0 )
    {
      valueStream << setprecision (precision) << fixed << (double)obj->getVal() << " +- " << (double)obj->getError();
//      cout << setprecision (precision) << fixed << (double)obj->getVal() << " +- " << (double)obj->getError() << endl;;
    } else {
       valueStream << setprecision (precision) << fixed << (double)obj->getVal();
    }
    string valueString = valueStream.str();
    latexLabel.DrawLatex(0.60, position, Form("%s = %s %s", obj->GetTitle(), valueString.c_str(), (char*)obj->getUnit()));
		obj = (RooRealVar*)it->Next();
  }
//  cout << "it->Next()->GetName()= " << it->Next()->GetName() << "\tit->Next()->getVal()= " << it->Next()->getVal() << endl;
/*
  position -= 0.04;
  std::ostringstream valueStream;
  valueStream << setprecision (3) << fixed << (double)(frame0->chiSquare(Form("model_%s_class_%i", isSignal?"signal":"background", iclass), datasetName.c_str(), isSignal ? 7 : 5));
  string valueString = valueStream.str();
  latexLabel.DrawLatex(0.60, position, Form("#chi^{2} / ndf = %s", valueString.c_str()));
*/
	return;
}
