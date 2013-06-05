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
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
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

//	TFile *infile = TFile::Open("simple.root");
	TFile *infile = TFile::Open("simple_parton.root");
	TTree *intree = (TTree*)infile->Get("Radion_m300_8TeV_nm");
//	TFile *infilereg = TFile::Open("simple_reg.root");
	TFile *infilereg = TFile::Open("simple_reg_parton.root");
	TTree *intreereg = (TTree*)infilereg->Get("Radion_m300_8TeV_nm");
	
	// Observables
	RooRealVar jj_mass("jj_mass", "m_{jj}", 50., 200., "GeV");
//	RooRealVar jj_mass("jj_mass", "m_{jj}", 80., 200., "GeV");
	jj_mass.setBins(50);
	RooRealVar ggjj_mass("ggjj_mass", "m_{jj#gamma#gamma}", 200., 400., "GeV");
	ggjj_mass.setBins(50);
	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	RooDataSet dataset("radion", "radion", intree, RooArgList(jj_mass, ggjj_mass));
	RooDataSet datasetreg("radion", "radion", intreereg, RooArgList(jj_mass, ggjj_mass));

	double mean_jj = dataset.mean(jj_mass);
	double mean_regjj = datasetreg.mean(jj_mass);
	double rms_jj = dataset.rmsVar(jj_mass)->getVal();
	double rms_regjj = datasetreg.rmsVar(jj_mass)->getVal();

	double mean_ggjj = dataset.mean(ggjj_mass);
	double mean_regggjj = datasetreg.mean(ggjj_mass);
	double rms_ggjj = dataset.rmsVar(ggjj_mass)->getVal();
	double rms_regggjj = datasetreg.rmsVar(ggjj_mass)->getVal();

	bool GAUSS = true;
	bool VOIGT = false;
	bool SIMVOIGT = false;
	bool CRYSTALBALL = false;

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
	dataset.plotOn(jj_frame, LineColor(kGreen+3));
	datasetreg.plotOn(jj_frame, LineColor(kRed+2));
//	RooFitResult *f = gauss.fitTo(dataset, Save(), Range(mean_jj-0.9*rms_jj, mean_jj+1.5*rms_jj));
	RooFitResult *f = gauss.fitTo(dataset, Save(), Range(mean_jj-1.5*rms_jj, mean_jj+1.5*rms_jj));
	gauss.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
	RooFitResult * freg = gaussreg.fitTo(datasetreg, Save(), Range(mean_regjj-1.5*rms_regjj, mean_regjj+1.5*rms_regjj));
	gaussreg.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
	RooArgList p = f->floatParsFinal();
	RooArgList preg = freg->floatParsFinal();	
	jj_frame->Draw();
	double res = sigma_gauss.getVal() / mu_gauss.getVal() * 100.;
	double resreg = sigma_gaussreg.getVal() / mu_gaussreg.getVal() * 100.;
	plotParameters( &p, c1, 0, jj_frame, true, 1, "Gaussian", 98, 99, 2);
	plotParameters( &preg, c1, 0, jj_frame, true, 3, "test", res, resreg, 2);
	c1->Print("pdf/mjj_gauss.pdf");
	c1->Print("root/mjj_gauss.root");
	c1->Print("gif/mjj_gauss.gif");
	c1->Clear();

	RooPlot * ggjj_frame = ggjj_mass.frame();
	dataset.plotOn(ggjj_frame, LineColor(kGreen+3));
	datasetreg.plotOn(ggjj_frame, LineColor(kRed+2));
	RooFitResult *f_ = gauss_.fitTo(dataset, Save(), Range(mean_ggjj-1.5*rms_ggjj, mean_ggjj+1.5*rms_ggjj));
	gauss_.plotOn(ggjj_frame, LineColor(kGreen+3), LineWidth(2));
	RooFitResult * freg_ = gauss_reg.fitTo(datasetreg, Save(), Range(mean_regggjj-1.5*rms_regggjj, mean_regggjj+1.5*rms_regggjj));
	gauss_reg.plotOn(ggjj_frame, LineColor(kRed+2), LineWidth(2));
	RooArgList p_ = f_->floatParsFinal();
	RooArgList preg_ = freg_->floatParsFinal();	
	ggjj_frame->Draw();
	double res_ = sigma_gauss_.getVal() / mu_gauss_.getVal() * 100.;
	double resreg_ = sigma_gauss_reg.getVal() / mu_gauss_reg.getVal() * 100.;
	plotParameters( &p_, c1, 0, ggjj_frame, true, 1, "Gaussian", 98, 99, 2);
	plotParameters( &preg_, c1, 0, ggjj_frame, true, 3, "test", res_, resreg_, 2);
	c1->Print("pdf/mggjj_gauss.pdf");
	c1->Print("root/mggjj_gauss.root");
	c1->Print("gif/mggjj_gauss.gif");
	c1->Clear();

	cout << "res= " << res << "\tresreg= " << resreg << "\timprov= " << fabs(res-resreg)/res*100. << endl;
	cout << "res_= " << res_ << "\tresreg_= " << resreg_ << "\timprov= " << fabs(res_-resreg_)/res_*100. << endl;
}
if( CRYSTALBALL )
{

	RooRealVar mu_CrystalBall("mu_CrystalBall", "mean", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
	RooRealVar sigma_CrystalBall("sigma_CrystalBall", "sigma", rms_jj, .0, 5*rms_jj, "GeV");
	RooRealVar alpha_CrystalBall("alpha_CrystalBall", "alpha", 1., 0., 20., "GeV");
	RooRealVar n_CrystalBall("n_CrystalBall", "n", 20., 0., 30., "GeV");
	RooCBShape CrystalBall("CrystalBall", "CrystalBall", jj_mass, mu_CrystalBall, sigma_CrystalBall, alpha_CrystalBall, n_CrystalBall);
	RooGaussian gauss("gauss", "gauss", jj_mass, mu_CrystalBall, sigma_CrystalBall);

	RooRealVar mu_CrystalBallreg("mu_CrystalBallreg", "mean (reg)", mean_regjj, 50., 200., "GeV");
	RooRealVar sigma_CrystalBallreg("sigma_CrystalBallreg", "sigma (reg)", 5., 0.0001, 50., "GeV");
	RooRealVar alpha_CrystalBallreg("alpha_CrystalBallreg", "alpha (reg)", 35., 5., 50., "GeV");
	RooRealVar n_CrystalBallreg("n_CrystalBallreg", "n (reg)", 20., 5., 90., "GeV");
	RooCBShape CrystalBallreg("CrystalBallreg", "CrystalBallreg", jj_mass, mu_CrystalBallreg, sigma_CrystalBallreg, alpha_CrystalBallreg, n_CrystalBallreg);

	RooPlot * jj_frame = jj_mass.frame();
	dataset.plotOn(jj_frame, LineColor(kGreen+3));
	datasetreg.plotOn(jj_frame, LineColor(kRed+2));
	gauss.fitTo(dataset, Save(), Range(mean_jj-1.5*rms_jj, mean_jj+1.5*rms_jj));
	gauss.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
	RooFitResult *f = CrystalBall.fitTo(dataset, Save());
	CrystalBall.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
//	RooFitResult * freg = CrystalBallreg.fitTo(datasetreg, Save());
//	CrystalBallreg.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
	RooArgList p = f->floatParsFinal();
//	RooArgList preg = freg->floatParsFinal();	
	jj_frame->Draw();
	plotParameters( &p, c1, 0, jj_frame, true, 1, "CrystalBall", 98, 99, 2);
//	plotParameters( &preg, c1, 0, jj_frame, true, 5, "test", 98, 99, 2);
	c1->Print("pdf/mjj_CrystalBall.pdf");
	c1->Print("root/mjj_CrystalBall.root");
	c1->Print("gif/mjj_CrystalBall.gif");
	c1->Clear();

}
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
/*
	RooCategory sample("sample", "sample");
	sample.defineType("vanilla");
	sample.defineType("regression");

	RooDataSet combData("combdata", "combdata", RooArgList(jj_mass, ggjj_mass), Index(sample), Import("vanilla", dataset), Import("regression", datasetreg));
	RooSimultaneous sim("sim", "sim", sample);
	sim.addPdf(voigt, "vanilla");
	sim.addPdf(voigtreg2, "regression");

	RooPlot * jj_frame = jj_mass.frame();
	dataset.plotOn(jj_frame, LineColor(kGreen+3));
	datasetreg.plotOn(jj_frame, LineColor(kRed+2));
	RooFitResult *f = sim.fitTo(combData, Save());
	voigt.plotOn(jj_frame, LineColor(kGreen+3), LineWidth(2));
//	RooFitResult * freg = voigtreg.fitTo(datasetreg, Save());
	voigtreg2.plotOn(jj_frame, LineColor(kRed+2), LineWidth(2));
	RooArgList p = f->floatParsFinal();
//	RooArgList preg = freg->floatParsFinal();	
	jj_frame->Draw();
	plotParameters( &p, c1, 0, jj_frame, true, 1, "test", 98, 99, 2);
//	plotParameters( &preg, c1, 0, jj_frame, true, 4, "test", 98, 99, 4);
	c1->Print("mjj_simvoigt.pdf");
	c1->Clear();
*/


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
  double fg = 2. * sigma_voigt.getVal() * sqrt(2. * log(2));
  double fl = 2. * width_voigt.getVal();
  double fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
  double fgreg = 2. * sigma_voigtreg.getVal() * sqrt(2. * log(2));
  double flreg = 2. * width_voigtreg.getVal();
  double fvreg = 0.5346 * flreg + sqrt(0.2166 * pow(flreg, 2.) + pow(fgreg, 2.));
  double res=  fv / mu_voigt.getVal() * 100. / (2. * sqrt(2. * log(2.)));
  double resreg=  fvreg / mu_voigtreg.getVal() * 100. / (2. * sqrt(2. * log(2.)));
	plotParameters( &p, c1, 0, jj_frame, true, 1, "Voigtian", 98, 99, 2);
	plotParameters( &preg, c1, 0, jj_frame, true, 4, "test", res, resreg, 2);
	c1->Print("pdf/mjj_voigt.pdf");
	c1->Print("root/mjj_voigt.root");
	c1->Print("gif/mjj_voigt.gif");
	c1->Clear();

	RooRealVar mu_voigt_("mu_voigt_", "mean", 300., 200., 400., "GeV");
	RooRealVar width_voigt_("width_voigt_", "width", 35., 5., 50., "GeV");
	RooRealVar sigma_voigt_("sigma_voigt_", "sigma", 5., .0, 20., "GeV");
	RooVoigtian voigt_("voigt_", "voigt_", ggjj_mass, mu_voigt_, width_voigt_, sigma_voigt_);

	RooRealVar mu_voigt_reg("mu_voigt_reg", "mean (reg)", 300., 200., 400., "GeV");
	RooRealVar width_voigt_reg("width_voigt_reg", "width (reg)", 35., 5., 50., "GeV");
	RooRealVar sigma_voigt_reg("sigma_voigt_reg", "sigma (reg)", 5., 0.0001, 20., "GeV");
	RooVoigtian voigt_reg("voigt_reg", "voigt_reg", ggjj_mass, mu_voigt_reg, width_voigt_reg, sigma_voigt_reg);

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
	double fg_ = 2. * sigma_voigt_.getVal() * sqrt(2. * log(2));
	double fl_ = 2. * width_voigt_.getVal();
	double fv_ = 0.5346 * fl_ + sqrt(0.2166 * pow(fl_, 2.) + pow(fg_, 2.));
	double fgreg_ = 2. * sigma_voigt_reg.getVal() * sqrt(2. * log(2));
	double flreg_ = 2. * width_voigt_reg.getVal();
	double fvreg_ = 0.5346 * flreg_ + sqrt(0.2166 * pow(flreg_, 2.) + pow(fgreg_, 2.));
	double res_=  fv_ / mu_voigt_.getVal() * 100. / (2. * sqrt(2. * log(2.)));
	double resreg_=  fvreg_ / mu_voigt_reg.getVal() * 100. / (2. * sqrt(2. * log(2.)));
	plotParameters( &p_, c1, 0, ggjj_frame, true, 1, "Voigtian", 98, 99, 2);
	plotParameters( &preg_, c1, 0, ggjj_frame, true, 4, "test", res_, resreg_, 2);
	c1->Print("pdf/mggjj_voigt.pdf");
	c1->Print("root/mggjj_voigt.root");
	c1->Print("gif/mggjj_voigt.gif");
	c1->Clear();

  cout << "fg= " << fg << "\tfl= " << fl << "\tfv= " << fv << endl;
  cout << "fgreg= " << fgreg << "\tflreg= " << flreg << "\tfvreg= " << fvreg << endl;
  cout << "res= " << res << "\tresreg= " << resreg << "\timprov= " << fabs(res - resreg)/res * 100.<< endl;
/*
  cout << "mean_jj= " << mean_jj << "\trms_jj= " << rms_jj << endl;  cout << "mean_regjj= " << mean_regjj << "\trms_regjj= " << rms_regjj << endl;
	cout << "mean_regjj= " << mean_regjj << "\trms_regjj= " << rms_regjj << endl;
  cout << "## res= " << rms_jj / mean_jj << "\tresreg= " << rms_regjj / mean_regjj << "\timp= " << fabs(rms_jj / mean_jj - rms_regjj / mean_regjj)/(rms_jj / mean_jj)*100. << endl;
*/
	cout << "fg_= " << fg_ << "\tfl_= " << fl_ << "\tfv_= " << fv_ << endl;
	cout << "fgreg_= " << fgreg_ << "\tflreg_= " << flreg_ << "\tfvreg_= " << fvreg_ << endl;
	cout << "res_= " << res_ << "\tresreg_= " << resreg_ << "\timprov_= " << fabs(res_ - resreg_)/res_ * 100.<< endl;
/*
	cout << "mean_ggjj= " << mean_ggjj << "\trms_ggjj= " << rms_ggjj << endl;
	cout << "mean_regggjj= " << mean_regggjj << "\trms_regggjj= " << rms_regggjj << endl;
	cout << "## res= " << rms_ggjj / mean_ggjj << "\tresreg= " << rms_regggjj / mean_regggjj << "\timp= " << fabs(rms_ggjj / mean_ggjj - rms_regggjj / mean_regggjj)/(rms_ggjj / mean_ggjj)*100. << endl;
*/
}
/*
*/


	delete c1;
	c1 = 0;
	infile->Close();
	infilereg->Close();

	return 0;
}

void plotParameters(RooArgList *r2_cat0_param, TCanvas *c, int canvasDivision, RooPlot* frame0, bool isSignal, int iclass, string fitfunctionName, double mvaInf, double mvaSup, int precision )
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
  if(iclass == 1) latexLabel.DrawLatex(0.55, position, Form("Fit: %s", fitfunctionName.c_str()));
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
  	cout << "obj->GetName()= " << obj->GetName() << "\tobj->getVal()= " << ((RooRealVar*)obj)->getVal() << endl;
		obj = (RooRealVar*)it->Next();
  }
if(iclass != 1)
{
	position -= 0.04;
  std::ostringstream valueStream;
	valueStream << setprecision (precision) << fixed << mvaInf;
	string valueString = valueStream.str();
	latexLabel.DrawLatex(0.60, position, Form("res= %s %%", valueString.c_str()));
	position -= 0.04;
	valueStream.clear(); valueStream.str(""); valueString = "";
	valueStream << setprecision (precision) << fixed << mvaSup;
	valueString = valueStream.str();
	latexLabel.DrawLatex(0.60, position, Form("res (reg)= %s %%", valueString.c_str()));
	position -= 0.04;
	valueStream.clear(); valueStream.str(""); valueString = "";
	valueStream << setprecision (precision) << fixed << fabs(mvaInf - mvaSup) / mvaInf * 100.;
	valueString = valueStream.str();
	latexLabel.DrawLatex(0.60, position, Form("improvement= %s %%", valueString.c_str()));
}

/*
  position -= 0.04;
  std::ostringstream valueStream;
  valueStream << setprecision (3) << fixed << (double)(frame0->chiSquare(Form("model_%s_class_%i", isSignal?"signal":"background", iclass), fitfunctionName.c_str(), isSignal ? 7 : 5));
  string valueString = valueStream.str();
  latexLabel.DrawLatex(0.60, position, Form("#chi^{2} / ndf = %s", valueString.c_str()));
*/
	return;
}
