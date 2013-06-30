// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TEntryList.h"
// C++ headers
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

int main()
{
	TFile *inf = TFile::Open("jetTreeForTraining.root");
	TTree *tr  = (TTree*)inf->Get("jets");
	float jet_pt = 0.;
	tr->SetBranchAddress("jet_pt", &jet_pt);
	TCut preselectionCut("jet_genDR<0.4 && jet_csvBtag > 0.");

	tr->Draw(">>elist", preselectionCut, "entrylist");
	TEntryList *elist = (TEntryList*)gDirectory->Get("elist");

	int ievt = elist->GetEntry(0);
	vector<float> values;
	int ntot = tr->GetEntries(preselectionCut);
	for(int ievtlist = 0 ; ievtlist < ntot ; ievtlist++)
	{
		if(ievt < 0) continue;
		tr->GetEntry(ievt);
		values.push_back(jet_pt);
		ievt = elist->Next();
	} // end of event loop

	sort(values.begin(), values.end());
	cout << "min= " << values[0] << "\tmax= " << values[ntot-1] << endl;

	return 0;
}
