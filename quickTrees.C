// C++ headers
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

int main ()
{
//	TFile *infile = TFile::Open("Radion_m300_8TeV_nm_parton.root");
//	TFile *outfile = new TFile("simple_parton.root", "RECREATE");
	TFile *infile = TFile::Open("Radion_m300_8TeV_nm_genjet.root");
	TFile *outfile = new TFile("simple_genjet.root", "RECREATE");
	TTree *intree = (TTree*)infile->Get("Radion_m300_8TeV_nm");
	TTree *outtree = new TTree("Radion_m300_8TeV_nm", "Radion_m300_8TeV_nm noreg");

	float jj_mass, regjj_mass, ggjj_mass, regggjj_mass;
	int njets_kRadionID, njets_kRadionID_and_CSVM;
	intree->SetBranchAddress("jj_mass", &jj_mass);
	intree->SetBranchAddress("regjj_mass", &regjj_mass);
	intree->SetBranchAddress("ggjj_mass", &ggjj_mass);
	intree->SetBranchAddress("regggjj_mass", &regggjj_mass);
	intree->SetBranchAddress("njets_kRadionID", &njets_kRadionID);
	intree->SetBranchAddress("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM);

	outtree->Branch("jj_mass", &jj_mass, "jj_mass/F");
	outtree->Branch("ggjj_mass", &ggjj_mass, "ggjj_mass/F");
	outtree->Branch("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM, "njets_kRadionID_and_CSVM/I");

	for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
	{
		intree->GetEntry(ievt);
		outtree->Fill();
	}
  outfile->cd();
  outtree->Write();
  outfile->Close();

//	TFile *outfilereg = new TFile("simple_reg_parton.root", "RECREATE");
	TFile *outfilereg = new TFile("simple_reg_genjet.root", "RECREATE");
	TTree *outtreereg = new TTree("Radion_m300_8TeV_nm", "Radion_m300_8TeV_nm reg");
	outtreereg->Branch("jj_mass", &regjj_mass, "jj_mass/F");
	outtreereg->Branch("ggjj_mass", &regggjj_mass, "ggjj_mass/F");
	outtreereg->Branch("njets_kRadionID_and_CSVM", &njets_kRadionID_and_CSVM, "njets_kRadionID_and_CSVM/I");
	for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
	{
		intree->GetEntry(ievt);
		outtreereg->Fill();
	}

	outfilereg->cd();
  outtreereg->Write();
  outfilereg->Close();
  infile->Close();

	return 0;
}
