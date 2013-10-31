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


for (int is = 0 ; is < fitStrategy.size() ; is++)
{
	for(float addeff = 0.0; addeff < 0.3 ; addeff += 0.05)
	{
		for(int itag = 1 ; itag <= 2 ; itag++)
		{
		TFile file(Form("v10bis_%s_mgg_0/02013-10-29-Radion_m500_8TeV_nm_m500.root", fitStrategy[is].c_str()));
		TH1F *mggjj = (TH1F*)file.Get(Form("mggjj_%ibtag", itag));
		float effgoal = (itag == 1) ? eff_1btag[is] : eff_2btag[is];
		effgoal += addeff;
//		cout << "effgoal= " << effgoal << endl;
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
		cout << fitStrategy[is] << "\t" << itag << "btag\teffgoal @ 300= " << effgoal << "\teff= " << (float)mggjj->Integral(lowbin, highbin) / (float)mggjj->Integral(0, nbins+1) << "\tlowcut= " << lowcut << "\thighcut= " << highcut << endl; 
	}
	}
cout << "###################" << endl;
}
}
