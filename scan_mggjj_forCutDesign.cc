// get the mass window for 500 GeV hypothesis in case of mgg fit (not enough data stat)
// O. Bondu (Nov. 2013)
// C++ headers
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TH1F.h>
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;

int main()
{
	vector<string> fitStrategy;
	vector<double> eff_1btag;
	vector<double> eff_2btag;

	fitStrategy.push_back("base");
	eff_1btag.push_back(0.492225); eff_2btag.push_back(0.569633);
	fitStrategy.push_back("reg");
	eff_1btag.push_back(0.510863); eff_2btag.push_back(0.60829);
	fitStrategy.push_back("kin");
	eff_1btag.push_back(0.399883); eff_2btag.push_back(0.541969);
	fitStrategy.push_back("regkin");
	eff_1btag.push_back(0.409865); eff_2btag.push_back(0.581347);

for (unsigned int is = 0 ; is < fitStrategy.size() ; is++)
{
	for(int itag = 1 ; itag <= 2 ; itag++)
	{
		vector<double> lowcuts; lowcuts.clear();
		vector<double> highcuts; highcuts.clear();
		if( strcmp("base", fitStrategy[is].c_str()) == 0 )
		{
			if(itag == 1){ lowcuts.push_back(485.); highcuts.push_back(535.); }
			else { lowcuts.push_back(490.); highcuts.push_back(525.); }
		}
		else if( strcmp("reg", fitStrategy[is].c_str()) == 0 )
		{
			if(itag == 1){ lowcuts.push_back(495.); highcuts.push_back(555.); }
			else { lowcuts.push_back(485.); highcuts.push_back(515.); }
		}
		else if( strcmp("kin", fitStrategy[is].c_str()) == 0 )
		{
			if(itag == 1){ lowcuts.push_back(505.); highcuts.push_back(540.); }
			else { lowcuts.push_back(495.); highcuts.push_back(510.); }
		}
		else if( strcmp("regkin", fitStrategy[is].c_str()) == 0 )
		{
			if(itag == 1){ lowcuts.push_back(440.); highcuts.push_back(505.); }
			else { lowcuts.push_back(490.); highcuts.push_back(510.); }
		}
		for(float addeff = 0.00; addeff <= 0.55 ; addeff += 0.05)
		{
//			TFile file(Form("v15bis_%s_mgg_0_massCutVersion0/02013-10-31-Radion_m500_8TeV_nm_m500.root", fitStrategy[is].c_str()));
//			TFile file(Form("v15ter_%s_mgg_0_massCutVersion0/02013-11-01-Radion_m500_8TeV_nm_m500.root", fitStrategy[is].c_str()));
			TFile file(Form("v15quat_%s_mgg_0_massCutVersion0/02013-11-01-Radion_m500_8TeV_nm_m500.root", fitStrategy[is].c_str()));
			TH1F *mggjj = (TH1F*)file.Get(Form("mggjj_%ibtag", itag));
			float effgoal = (itag == 1) ? eff_1btag[is] : eff_2btag[is];
			effgoal += addeff;
			effgoal = min(effgoal, (float)1.);
			int nbins = mggjj->GetNbinsX();
			float ntot= mggjj->Integral(0, nbins+1);
			int lowbin = -1; int highbin = -1;	
			float closesteff = 1.;
			float eff = 0.;
			for(int ibin = 1; ibin < nbins ; ibin++)
			{
				for(int jbin = ibin+1 ; jbin <= nbins ; jbin++)
				{
					eff = (float)mggjj->Integral(ibin, jbin) / (float)ntot;
					if(lowcuts.size() > 0)
					{
//						if(mggjj->GetBinLowEdge(ibin) > lowcuts[lowcuts.size() - 1])
//							continue; // lowcut must not raise
//						if(mggjj->GetBinLowEdge(jbin + 1) < highcuts[highcuts.size() - 1])
//							continue; // highcut must not raise
//						if( (mggjj->GetBinLowEdge(jbin + 1) - mggjj->GetBinLowEdge(ibin)) < (highcuts[highcuts.size() - 1] - lowcuts[lowcuts.size() - 1]) )
//							continue; // window size must not decrease
					}
					if( fabs(eff - effgoal) < closesteff )
					{
						closesteff = fabs(eff - effgoal);
						lowbin = ibin;
						highbin = jbin;
					}
				}
			}
	//		cout << "closesteff= " << closesteff << "\tlowbin= " << lowbin << "\thighbin= " << highbin << "\teff= " << (float)mggjj->Integral(lowbin, highbin) / (float)mggjj->Integral(0, nbins+1) << endl;
			float lowcut = mggjj->GetBinLowEdge(lowbin);
			float highcut = mggjj->GetBinLowEdge(highbin + 1);
			cout << fitStrategy[is] << "\t" << itag << "btag\teffgoal(addeff)= " << effgoal << "\t(" << addeff << ")\teff= " << (float)mggjj->Integral(lowbin, highbin) / (float)mggjj->Integral(0, nbins+1) << "\tlowcut= " << lowcut << "\thighcut= " << highcut << endl; 
			lowcuts.push_back(lowcut);
			highcuts.push_back(highcut);
		}
		cout << endl;
		}
	cout << "###################" << endl;
	}
	
	return 0;
}
