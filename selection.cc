// Radion Selection implementation
// O. Bondu (May 2013)
// TMVA headers
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
// C++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <boost/program_options.hpp>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
// Analysis headers
#include "../h2gglobe/BTagUtils.h"
#include "../KinematicFit/DiJetKinFitter.h"
#include "selection.h"
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	// declare arguments
 	string inputfile;
	string inputtree;
	string outputfile;
	string outputtree;
	string regressionFolder;
	int numberOfRegressionFiles;
	int type; // Same conventions as in h2gglobe: <0 = signal ; =0 = data ; >0 = background
	int SYNC; // mjj and mggjj cuts are different for sync and analysis
	int SYNC_W_PHIL;
	int FULL_DUMP;
	int REMOVE_UNDEFINED_BTAGSF;
	int applyMassCuts;
	int applyPhotonIDControlSample;

	// print out passed arguments
	copy(argv, argv + argc, ostream_iterator<char*>(cout, " ")); cout << endl;
	// argument parsing
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message")
			("inputfile,i", po::value<string>(&inputfile)->default_value("root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_tree_v08/Radion_nm.root"), "input file")
			("inputtree,t", po::value<string>(&inputtree)->default_value("Radion_m300_8TeV_nm"), "input tree")
			("outputtree", po::value<string>(&outputtree)->default_value("Radion_m300_8TeV_nm"), "output tree")
			("outputfile,o", po::value<string>(&outputfile)->default_value("selected.root"), "output file")
			("regressionFolder", po::value<string>(&regressionFolder)->default_value("/afs/cern.ch/user/h/hebda/public/"), "regression folder")
			("numberOfRegressionFiles,r", po::value<int>(&numberOfRegressionFiles)->default_value(2), "number of regression files")
			("type", po::value<int>(&type)->default_value(0), "same conventions as in h2gglobe: <0 = signal ; =0 = data ; >0 = background")
			("sync", po::value<int>(&SYNC)->default_value(0), "mjj and mggjj cuts are overwritten if sync is switched on")
			("removeUndefinedBtagSF", po::value<int>(&REMOVE_UNDEFINED_BTAGSF)->default_value(0), "remove undefined btagSF_M values (should be used only for the limit trees)")
			("applyMassCuts", po::value<int>(&applyMassCuts)->default_value(1), "can switch off mass cuts (e.g. for control plots), prevails other mass cut options if switched off")
			("applyPhotonIDControlSample", po::value<int>(&applyPhotonIDControlSample)->default_value(0), "Invert photon ID CiC cut to populate selection in gjjj instead of ggjj")
			("sync_w_phil", po::value<int>(&SYNC_W_PHIL)->default_value(0), "switch on output for dedicated events")
			("full_dump", po::value<int>(&FULL_DUMP)->default_value(0), "switch on creation of the t.event dump")
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

//	string outputtree = inputtree;

	if(DEBUG) cout << "End of argument parsing" << endl;




	cout << "inputfile= " << inputfile << endl;
	cout << "inputtree= " << inputtree << endl;
	cout << "outputfile= " << outputfile << endl;
	cout << "outputtree= " << outputtree << endl;
	cout << "regressionFolder= " << regressionFolder << endl;

	TFile *infile = TFile::Open(inputfile.c_str());
	TTree *intree = (TTree*)infile->Get(inputtree.c_str());
	TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
	TTree *outtree = new TTree(outputtree.c_str(), Form("%s reduced", outputtree.c_str()));
	ofstream synchrofile, full_dump;
	if(SYNC) synchrofile.open("synchronisation.txt");
	if(FULL_DUMP) full_dump.open(Form("h2gglobe_%s.txt", inputtree.c_str()));

	if(DEBUG) cout << "Setup tree inputs" << endl;
	if(DEBUG) cout << "Setup tree outputs" << endl;
	tree_variables t;
	initialize_variables(&t);

	if(DEBUG) cout << "SetBranchAddresses" << endl;
	setup_intree(intree, &t, type);
	setup_outtree(outtree, &t);

	if(DEBUG) cout << "Setup TRandom3 generator" << endl;
	TRandom3 *r3 = new TRandom3();

	if(DEBUG) cout << "Prepare for regression" << endl;
// prepare for regression
	TMVA::Reader* readerRegres = new TMVA::Reader( "!Color:!Silent" );
	readerRegres->AddVariable( "jet_eta", &t.jet_eta); 
	readerRegres->AddVariable( "jet_emfrac", &t.jet_emfrac);
	readerRegres->AddVariable( "jet_hadfrac", &t.jet_hadfrac);
	readerRegres->AddVariable( "jet_nconstituents", &t.jet_nConstituents_);
	readerRegres->AddVariable( "jet_vtx3dL", &t.jet_secVtx3dL);
	readerRegres->AddVariable( "MET", &t.met_corr_pfmet);
	readerRegres->AddVariable( "jet_dPhiMETJet", &t.jet_dPhiMet_fabs);
//	readerRegres->AddVariable( "rho25", &t.rho); // Added for Phil Oct 28 (later: removed)

// Adding variables
	if(numberOfRegressionFiles != 0 && numberOfRegressionFiles != 2)
	{
		cout << "ERROR: current version must have two regression files or no regression" << endl;
		return 1;
//		readerRegres->BookMVA("BDT", regressionFolder.c_str());
	} else {
		for(int i = 0; i < numberOfRegressionFiles ; i++)
		{
//			readerRegres->BookMVA(Form("BDT_%i", i), Form("%s/TMVARegression_10_Cat%i_BDTG.weights.xml", regressionFolder.c_str(), i)); // Phil, Oct 17
			readerRegres->BookMVA(Form("BDT_%i", i), Form("%s/jetRegressionWeights/TMVARegression_Cat%i_BDTG.weights.xml", regressionFolder.c_str(), i)); // Phil, Oct 28
		}
	}


	int nevents[30] = {0};
	int ilevelmax=0;
	float nevents_w[30] = {0.};
	float nevents_w_btagSF_M[30] = {0.};
	int nevents_sync[30] = {0};
	int nevents1btag[30] = {0};
	int nevents2btag[30] = {0};
	string eventcut[30];
	int njets[30] = {0};
	string jetcut[30];
  int decade = 0;
  int totevents = intree->GetEntries();
  if(DEBUG) totevents = 1;
  cout << "#entries= " << totevents << endl;
  // loop over events
  for(int ievt=0 ; ievt < totevents ; ievt++)
  {
		int ilevel = 0;
		if(DEBUG) cout << "#####\tievt= " << ievt << endl;
    double progress = 10.0*ievt/(1.0*totevents);
    int k = TMath::FloorNint(progress);
    if (k > decade && !DEBUG) cout<<10*k<<" %"<<endl;
    decade = k;

		int njets_kRadionID_ = 0;
		int njets_kRadionID_and_CSVM_ = 0;
    intree->GetEntry(ievt);
		if(DEBUG && SYNC_W_PHIL && !(/*t.event == 6976 ||*/ t.event == 8042 || t.event == 14339 /*|| t.event == 2227 || t.event == 4921 || t.event == 7665 || t.event == 7687 || t.event == 11246 || t.event == 15140 || t.event == 685*/) ) continue;
		if(DEBUG) cout << "#####\tievt= " << ievt << "\trun= " << t.run << "\tlumi= " << t.lumis << "\tevent= " << t.event << endl;
		if( numberOfRegressionFiles != 0 && type < -250 && ((int)t.event % 2 == 0) && !SYNC_W_PHIL) continue; // use regression only on odd events
	
		if(DEBUG) cout << "for MC, get the MC truth hjj system" << endl;
// Compute hjj system
		TLorentzVector gj1, gj2;
		if( t.gr_j1_p4_pt > .01 && t.gr_j2_p4_pt > .01)
		{
			gj1.SetPtEtaPhiM(t.gr_j1_p4_pt, t.gr_j1_p4_eta, t.gr_j1_p4_phi, t.gr_j1_p4_mass);
			gj2.SetPtEtaPhiM(t.gr_j2_p4_pt, t.gr_j2_p4_eta, t.gr_j2_p4_phi, t.gr_j2_p4_mass);
			TLorentzVector hjj = gj1 + gj2;
			t.gr_hjj_p4_pt = hjj.Pt();
			t.gr_hjj_p4_eta = hjj.Eta();
			t.gr_hjj_p4_phi = hjj.Phi();
			t.gr_hjj_p4_mass = hjj.M();
		} else {
			t.gr_hjj_p4_pt = 0.;
			t.gr_hjj_p4_eta = 0.;
			t.gr_hjj_p4_phi = 0.;
			t.gr_hjj_p4_mass = 0.;
		}

		if(DEBUG) cout << "Apply photon ID cuts" << endl;
		// Apply photon ID cuts
		nevents[ilevel]++; eventcut[ilevel] = "Before photon ID";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[0]++;
		float ph1_PFisoA = (t.ph1_pfchargedisogood03 + t.ph1_ecaliso + 2.5 - t.rho * 0.09) * 50. / t.ph1_pt;
		float ph1_PFisoB = (t.ph1_pfchargedisobad04 + t.ph1_ecalisobad + 2.5 - t.rho * 0.23) * 50. / t.ph1_badvtx_Et;
		float ph1_PFisoC = t.ph1_pfchargedisogood03 * 50. / t.ph1_pt;
		float ph2_PFisoA = (t.ph2_pfchargedisogood03 + t.ph2_ecaliso + 2.5 - t.rho * 0.09) * 50. / t.ph2_pt;
		float ph2_PFisoB = (t.ph2_pfchargedisobad04 + t.ph2_ecalisobad + 2.5 - t.rho * 0.23) * 50. / t.ph2_badvtx_Et;
		float ph2_PFisoC = t.ph2_pfchargedisogood03 * 50. / t.ph2_pt;

		if(DEBUG) cout << "t.ph1_pt= " << t.ph1_pt << "\tph1_eta= " << t.ph1_eta << "\tph1_phi= " << t.ph1_phi << "\tph1_r9= " << t.ph1_r9 << "\tph1_SCEta= " << t.ph1_SCEta << endl;
		if(DEBUG) cout << "t.ph2_pt= " << t.ph2_pt << "\tph2_eta= " << t.ph2_eta << "\tph2_phi= " << t.ph2_phi << "\tph2_r9= " << t.ph2_r9 << "\tph2_SCEta= " << t.ph2_SCEta << endl;
		if(DEBUG) cout << "t.ph1_pfchargedisogood03= " << t.ph1_pfchargedisogood03 << "\tph1_ecaliso= " << t.ph1_ecaliso << "\trho= " << t.rho << "\tph1_pt= " << t.ph1_pt << endl;
		if(DEBUG) cout << "t.ph1_pfchargedisobad04= " << t.ph1_pfchargedisobad04 << "\tph1_ecalisobad= " << t.ph1_ecalisobad << "\trho= " << t.rho << "\tph1_badvtx_Et= " << t.ph1_badvtx_Et << endl;
		if(DEBUG) cout << "t.ph1_pfchargedisogood03= " << t.ph1_pfchargedisogood03 << "\tph1_pt= " << t.ph1_pt << endl;
		if(DEBUG) cout << "t.ph2_pfchargedisogood03= " << t.ph2_pfchargedisogood03 << "\tph2_ecaliso= " << t.ph2_ecaliso << "\trho= " << t.rho << "\tph2_pt= " << t.ph2_pt << endl;
		if(DEBUG) cout << "t.ph2_pfchargedisobad04= " << t.ph2_pfchargedisobad04 << "\tph2_ecalisobad= " << t.ph2_ecalisobad << "\trho= " << t.rho << "\tph2_badvtx_Et= " << t.ph2_badvtx_Et << endl;
		if(DEBUG) cout << "t.ph2_pfchargedisogood03= " << t.ph2_pfchargedisogood03 << "\tph2_pt= " << t.ph2_pt << endl;
		if(DEBUG) cout << "ph1_PFisoA= " << ph1_PFisoA << "\tph1_PFisoB= " << ph1_PFisoB << "\tph1_PFisoC= " << ph1_PFisoC << "\tph1_sieie= " << t.ph1_sieie << "\tph1_hoe= " << t.ph1_hoe << "\tph1_isconv= " << t.ph1_isconv << endl;
		if(DEBUG) cout << "ph2_PFisoA= " << ph2_PFisoA << "\tph2_PFisoB= " << ph2_PFisoB << "\tph2_PFisoC= " << ph2_PFisoC << "\tph2_sieie= " << t.ph2_sieie << "\tph2_hoe= " << t.ph2_hoe << "\tph2_isconv= " << t.ph2_isconv << endl;
		if(DEBUG) cout << "t.j1_pt= " << t.j1_pt << "\tj1_csvBtag= " << t.j1_csvBtag << "\tj1_csvMvaBtag= " << t.j1_csvMvaBtag << endl;
		if(DEBUG) cout << "t.j2_pt= " << t.j2_pt << "\tj2_csvBtag= " << t.j2_csvBtag << "\tj2_csvMvaBtag= " << t.j2_csvMvaBtag << endl;
		if(DEBUG) cout << "t.j3_pt= " << t.j3_pt << "\tj3_csvBtag= " << t.j3_csvBtag << "\tj3_csvMvaBtag= " << t.j3_csvMvaBtag << endl;
		if(DEBUG) cout << "t.j4_pt= " << t.j4_pt << "\tj4_csvBtag= " << t.j4_csvBtag << "\tj4_csvMvaBtag= " << t.j4_csvMvaBtag << endl;

		if(DEBUG) cout << "t.ph1_pt= " << t.ph1_pt << "\t(float)(40.*t.PhotonsMass)/(float)120.= " << (float)(40.*t.PhotonsMass)/(float)120. << endl;
		if( t.ph1_pt < (float)(40.*t.PhotonsMass)/(float)120. ) continue;
		nevents[ilevel]++; eventcut[ilevel] = "After floating pt cut for photon 1 (40*mgg/120 GeV)";
		nevents_w[ilevel] += t.evweight; ilevel++;
//		if( t.ph2_pt < 25. ) continue;
		if(DEBUG) cout << "t.ph2_pt= " << t.ph2_pt << "\t(float)(30.*t.PhotonsMass)/(float)120.= " << (float)(30.*t.PhotonsMass)/(float)120. << endl;
		if( t.ph2_pt < (float)(30.*t.PhotonsMass)/(float)120. ) continue; // switching to running pt cut per Hgg recommendations (Nov. 2013)
		nevents[ilevel]++; eventcut[ilevel] = "After fixed pt cut for photon 2 (25 GeV)";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[1]++;
		if(DEBUG) cout << "t.ph1_ciclevel= " << t.ph1_ciclevel << "\tph2_ciclevel= " << t.ph2_ciclevel << endl;
		if( (!applyPhotonIDControlSample) && ((t.ph1_ciclevel < 4) || (t.ph2_ciclevel < 4)) ) continue;
		else if (applyPhotonIDControlSample)
		{
			bool ph1_id = (t.ph1_ciclevel >= 3);
			bool ph2_id = (t.ph2_ciclevel >= 3);
			bool ph1_Lid = (t.ph1_ciclevel >= 0) && (t.ph1_ciclevel < 3);
			bool ph2_Lid = (t.ph2_ciclevel >= 0) && (t.ph2_ciclevel < 3);
			if(ph1_id && ph2_id) continue; // reject gg
			if(ph1_Lid && ph2_Lid) continue; // reject jj
			if( !( (ph1_id && ph2_Lid) || (ph2_id && ph1_Lid)) ) continue; // reject if different from gj or jg
		}
		nevents[ilevel]++; eventcut[ilevel] = "After cic cut on both photons";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[2]++;
		if(DEBUG) cout << "t.PhotonsMass= " << t.PhotonsMass << endl;
		if( (t.PhotonsMass < 100.) || (t.PhotonsMass > 180.) ) continue;
		nevents[ilevel]++; eventcut[ilevel] = "After 100 < mgg < 180";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[3]++;
		t.HT_gg = t.ph1_pt + t.ph2_pt;

		// take only the subset of events where at least two jets remains
		if(DEBUG) cout << "t.njets_passing_kLooseID= " << t.njets_passing_kLooseID << endl;
		if( t.njets_passing_kLooseID < 2 ) continue;
		nevents[ilevel]++; eventcut[ilevel] = "After njet >= 2";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[4]++;
// alternative counting: taking into account only the 4 jets stored !
		int nbjet_tmp = 0;
//		float csv_cut = 0.244; // CSVL
		float csv_cut = 0.679; // CSVM
		for( int ijet = 0 ; ijet < min(t.njets_passing_kLooseID, 15); ijet ++ )
		{
			if( ijet == 0 ){
				if(t.j1_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 1 ){
				if(t.j2_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 2 ){
				if(t.j3_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 3 ){
				if(t.j4_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 4 ){
				if(t.j5_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 5 ){
				if(t.j6_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 6 ){
				if(t.j7_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 7 ){
				if(t.j8_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 8 ){
				if(t.j9_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 9 ){
				if(t.j10_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 10 ){
				if(t.j11_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 11 ){
				if(t.j12_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 12 ){
				if(t.j13_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 13 ){
				if(t.j14_csvBtag > csv_cut)
					nbjet_tmp++; 
			} else if( ijet == 14 ){
				if(t.j15_csvBtag > csv_cut)
					nbjet_tmp++;
			}
		}
		if( DEBUG && t.njets_passing_kLooseID > 4 ) cout << "t.njets_passing_kLooseID= " << t.njets_passing_kLooseID << "\tnbjet_tmp= " << nbjet_tmp << endl;
		if(DEBUG) cout << "nbjet_tmp= " << nbjet_tmp << endl;
		if( nbjet_tmp < 1 ) continue;
		nevents[ilevel]++; eventcut[ilevel] = "After nbjet >= 1";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[5]++;
		if( nbjet_tmp == 1) nevents1btag[10]++; 
		else nevents2btag[10]++;

		jet_variables J;
		initialize_jet_variables( &J );

		TLorentzVector met;
		met.SetPtEtaPhiE(t.met_corr_pfmet, t.met_corr_eta_pfmet, t.met_corr_phi_pfmet, t.met_corr_e_pfmet);

		// loop over jets, store jet info + info on closest genjet / parton (no selection applied)
		if(DEBUG) cout << "t.njets_passing_kLooseID= " << t.njets_passing_kLooseID << endl;
		for( int ijet = 0 ; ijet < min(t.njets_passing_kLooseID, 15); ijet ++ )
		{
			njets[0]++; jetcut[0] = "Before JetID";
			fill_jet_variables( &t, ijet, met);
			t.jet_nConstituents_ = (float) t.jet_nConstituents;
			t.jet_dPhiMet_fabs = fabs(t.jet_dPhiMet);

			if(DEBUG && numberOfRegressionFiles != 0) cout << "input= " << t.jet_pt << "\toutput (BDT_0)= " << readerRegres->EvaluateRegression(Form("BDT_%i", 0))[0] * t.jet_pt << "\toutput (BDT_1)= " << readerRegres->EvaluateRegression(Form("BDT_%i", 1))[0] * t.jet_pt << endl;
			if( REMOVE_UNDEFINED_BTAGSF && t.jet_flavour == 0. ) continue;
//			njets[5]++; jetcut[5] = "After t.jet_csvBtag < 0.";
			if(DEBUG) cout << "now with the regression" << endl;
			if(numberOfRegressionFiles == 0)
				t.jet_regPt = t.jet_pt; // no regression applied
//			else if(numberOfRegressionFiles <= 1)
//				t.jet_regPt = (float)(readerRegres->EvaluateMVA("BDT"));
			else
//				t.jet_regPt = (float)(readerRegres->EvaluateRegression(Form("BDT_%i", t.jet_pt < 80. ? 0 : 1))[0]) * t.jet_pt; // Phil Oct 17
				t.jet_regPt = (float)(readerRegres->EvaluateRegression(Form("BDT_%i", t.jet_pt < 90. ? 0 : 1))[0]) * t.jet_pt; // Phil Oct 28
			t.jet_regkinPt = t.jet_regPt;
			// jet selection
			// ** acceptance + pu id **
			if( t.jet_regPt < 25. ) continue;
			njets[1]++; jetcut[1] = "After jet pt > 25";
			if( fabs(t.jet_eta) > 2.5 ) continue;
			njets[2]++; jetcut[2] = "After jet |eta| < 2.5";
			if( t.jet_betaStarClassic > 0.2 * log( t.nvtx - 0.64) ) continue;
			njets[3]++; jetcut[3] = "After t.jet_betaStarClassic > 0.2 * log( t.nvtx - 0.64)";
			if( t.jet_dR2Mean > 0.06 ) continue;
			njets[4]++; jetcut[4] = "After t.jet_dR2Mean > 0.06";
			if(DEBUG) cout << "Jet is passing selection cuts" << endl;
			// ** call regression to correct the pt **
			// ** store 4-momentum + csv output for combinatorics **
			J.jetPt.push_back(t.jet_pt);
			J.jetbtagSF_M.push_back(t.jet_btagSF_M);
			J.jetflavour.push_back(t.jet_flavour);
			J.jetbtagSFErrorUp_M.push_back(t.jet_btagSFErrorUp_M);
			J.jetbtagSFErrorDown_M.push_back(t.jet_btagSFErrorDown_M);
			J.jetbtagEff_M.push_back(t.jet_btagEff_M);
			J.jetbtagEffError_M.push_back(t.jet_btagEffError_M);
			J.jetdR2Mean.push_back(t.jet_dR2Mean);
			J.jetbetaStarClassic.push_back(t.jet_betaStarClassic);
			J.jetE.push_back(t.jet_e);
			J.jetEta.push_back(t.jet_eta);
			J.jetPhi.push_back(t.jet_phi);
			J.jetCSV.push_back(t.jet_csvBtag);
			J.jetRegPt.push_back(t.jet_regPt);
			J.jetRegKinPt.push_back(t.jet_regkinPt);
			J.jetEmfrac.push_back(t.jet_emfrac);
			J.jetHadfrac.push_back(t.jet_hadfrac);
			J.jetSecVtxPt.push_back(t.jet_secVtxPt);
			J.jetSecVtx3dL.push_back(t.jet_secVtx3dL);
			J.jetDPhiMet.push_back(t.jet_dPhiMet);
			J.jetNConstituents.push_back(t.jet_nConstituents);
// Jet Energy Correction and Jet Energy Resolution
			J.jetJecD_e.push_back(t.jet_jecD_e); J.jetJecD_pt.push_back(t.jet_jecD_pt); J.jetJecD_phi.push_back(t.jet_jecD_phi); J.jetJecD_eta.push_back(t.jet_jecD_eta);
			J.jetJecU_e.push_back(t.jet_jecU_e); J.jetJecU_pt.push_back(t.jet_jecU_pt); J.jetJecU_phi.push_back(t.jet_jecU_phi); J.jetJecU_eta.push_back(t.jet_jecU_eta);
			J.jetJerD_e.push_back(t.jet_jerD_e); J.jetJerD_pt.push_back(t.jet_jerD_pt); J.jetJerD_phi.push_back(t.jet_jerD_phi); J.jetJerD_eta.push_back(t.jet_jerD_eta);
			J.jetJerC_e.push_back(t.jet_jerC_e); J.jetJerC_pt.push_back(t.jet_jerC_pt); J.jetJerC_phi.push_back(t.jet_jerC_phi); J.jetJerC_eta.push_back(t.jet_jerC_eta);
			J.jetJerU_e.push_back(t.jet_jerU_e); J.jetJerU_pt.push_back(t.jet_jerU_pt); J.jetJerU_phi.push_back(t.jet_jerU_phi); J.jetJerU_eta.push_back(t.jet_jerU_eta);


			njets_kRadionID_++;
			if(t.jet_csvBtag > csv_cut) njets_kRadionID_and_CSVM_++;
		} // end of loop over jets
		
		if(DEBUG) cout << "J.jetPt.size()= " << J.jetPt.size() << endl;
		// jet combinatorics
		if( J.jetPt.size() < 2 ) continue;
		nevents[ilevel]++; eventcut[ilevel] = "After njet >=2 passing the jet selection";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[6]++;

		vector<int> btaggedJet;
		for( unsigned int ijet = 0 ; ijet < J.jetPt.size() ; ijet++ )
		{
			if( J.jetCSV[ijet] > csv_cut )
				btaggedJet.push_back(ijet);
		}

		if( btaggedJet.size() < 1 ) continue;
		nevents[ilevel]++; eventcut[ilevel] = "After nbjet >=1 passing the jet selection";
		nevents_w[ilevel] += t.evweight; ilevel++;
		nevents_sync[7]++;
		if( btaggedJet.size() == 1) nevents1btag[12]++; 
		else nevents2btag[12]++;



		int ij1, ij2;
		int ij1Reg, ij2Reg;
		int ij1RegKin, ij2RegKin;
	
		if(DEBUG) cout << "btaggedJet.size()= " << btaggedJet.size() << endl;
		if(DEBUG)
			for(int ijet_=0; ijet_ < (int)btaggedJet.size() ; ijet_++)
				cout << "J.jetPt[btaggedJet[" << ijet_ << "]]= " << J.jetPt[btaggedJet[ijet_]] << endl;
		// if exactly one btag, pick it up, then find the other jet that gives max ptjj
		if( btaggedJet.size() == 1 )
		{
			if(DEBUG) cout << "Entering jet combinatorics: 1btag category" << endl;
			t.category = 1;
			unsigned int ij = btaggedJet[0];
			if(DEBUG) cout << "btaggedJet[0]= " << btaggedJet[0] << endl;
			TLorentzVector j, jreg, jregkin;
			j.SetPtEtaPhiE(J.jetPt[ij], J.jetEta[ij], J.jetPhi[ij], J.jetE[ij]);
			jreg = ((float)J.jetRegPt[ij]/(float)J.jetPt[ij]) * j;
			jregkin = ((float)J.jetRegKinPt[ij]/(float)J.jetPt[ij]) * j;
			int imaxptjj;
			int imaxptjjReg;
			int imaxptjjRegKin;
			float maxptjj = -99.;
			float maxptjjReg = -99.;
			float maxptjjRegKin = -99.;
			for(unsigned int ijet = 0 ; ijet < J.jetPt.size() ; ijet++)
			{
				if( ijet == ij ) continue;
				TLorentzVector tmp_j;
				TLorentzVector tmp_jReg;
				TLorentzVector tmp_jRegKin;
				tmp_j.SetPtEtaPhiE(J.jetPt[ijet], J.jetEta[ijet], J.jetPhi[ijet], J.jetE[ijet]);
				tmp_jReg = ((float)J.jetRegPt[ijet]/(float)J.jetPt[ijet]) * tmp_j;
				tmp_jRegKin = ((float)J.jetRegKinPt[ijet]/(float)J.jetPt[ijet]) * tmp_j;
				TLorentzVector jj = j + tmp_j;
				TLorentzVector jjReg = jreg + tmp_jReg;
				TLorentzVector jjRegKin = jregkin + tmp_jRegKin;
				if( jj.Pt() > maxptjj )
				{
					maxptjj = jj.Pt();
					imaxptjj = ijet;
				}
				if( jjReg.Pt() > maxptjjReg )
				{
					maxptjjReg = jjReg.Pt();
					imaxptjjReg = ijet; 
				}
				if( jjRegKin.Pt() > maxptjjRegKin )
				{
					maxptjjRegKin = jjRegKin.Pt();
					imaxptjjRegKin = ijet; 
				}
			}
			ij1 = ij;
			ij2 = imaxptjj;
			ij1Reg = ij;
			ij2Reg = imaxptjjReg;
			ij1RegKin = ij;
			ij2RegKin = imaxptjjRegKin;
		}
		// if two or more bjets, then loop only over btagged jets
		if( btaggedJet.size() > 1 )
		{
			t.category = 2;
			if(DEBUG) cout << "Entering jet combinatorics: 2btag category" << endl;
			int ij;
			int imaxptjj;
			int imaxptjjReg;
			int imaxptjjRegKin;
			int jmaxptjj;
			int jmaxptjjReg;
			int jmaxptjjRegKin;
			float maxptjj = -99.;
			float maxptjjReg = -99.;
			float maxptjjRegKin = -99.;
			for( unsigned int i = 0 ; i < btaggedJet.size() - 1 ; i++ )
			{
				ij = btaggedJet[i];
				if(DEBUG) cout << "btaggedJet[" << i << "]= " << ij << "\tjetPt[" << ij << "]= " << J.jetPt[ij] << endl;
				TLorentzVector j, jreg, jregkin;
				j.SetPtEtaPhiE(J.jetPt[ij], J.jetEta[ij], J.jetPhi[ij], J.jetE[ij]);
				jreg = ((float)J.jetRegPt[ij]/(float)J.jetPt[ij]) * j;
				jregkin = ((float)J.jetRegKinPt[ij]/(float)J.jetPt[ij]) * j;
				for(unsigned int k = i+1 ; k < btaggedJet.size() ; k++)
				{
					int ijet = btaggedJet[k];
					TLorentzVector tmp_j;
					TLorentzVector tmp_jReg;
					TLorentzVector tmp_jRegKin;
					tmp_j.SetPtEtaPhiE(J.jetPt[ijet], J.jetEta[ijet], J.jetPhi[ijet], J.jetE[ijet]);
					tmp_jReg = ((float)J.jetRegPt[ijet]/(float)J.jetPt[ijet]) * tmp_j;
					tmp_jRegKin = ((float)J.jetRegKinPt[ijet]/(float)J.jetPt[ijet]) * tmp_j;
					TLorentzVector jj = j + tmp_j;
					TLorentzVector jjReg = jreg + tmp_jReg;
					TLorentzVector jjRegKin = jregkin + tmp_jRegKin;
					if(DEBUG) cout << "btaggedJet[" << k << "]= " << btaggedJet[k] << "\tjetPt[" << ijet << "]= " << J.jetPt[ijet] << "\tjj.Pt()= " << jj.Pt() << "\t(maxptjj= " << maxptjj << ")" << endl;
					if( jj.Pt() > maxptjj )
					{
						maxptjj = jj.Pt();
						imaxptjj = ij;
						jmaxptjj = ijet;
					}
					if( jjReg.Pt() > maxptjjReg )
					{
						maxptjjReg = jjReg.Pt();
						imaxptjjReg = ij;
						jmaxptjjReg = ijet; 
					}
					if( jjRegKin.Pt() > maxptjjRegKin )
					{
						maxptjjRegKin = jjRegKin.Pt();
						imaxptjjRegKin = ij;
						jmaxptjjRegKin = ijet; 
					}
				}
			} // end of loop over bjets
			ij1 = imaxptjj;
			ij2 = jmaxptjj;
			ij1Reg = imaxptjjReg;
			ij2Reg = jmaxptjjReg;
			ij1RegKin = imaxptjjRegKin;
			ij2RegKin = jmaxptjjRegKin;
		} // end of if two bjets

		TLorentzVector pho1;
		TLorentzVector pho2;
		TLorentzVector jet1;
		TLorentzVector jet2;
		TLorentzVector regjet1;
		TLorentzVector regjet2;
		TLorentzVector regkinjet1;
		TLorentzVector regkinjet2;
		TLorentzVector kinjet1;
		TLorentzVector kinjet2;
		TLorentzVector pho1_pesD, pho2_pesD, pho1_pesU, pho2_pesU;
		TLorentzVector pho1_perD, pho2_perD, pho1_perU, pho2_perU;
		TLorentzVector jet1_jecD, jet1_jecU, jet1_jerD, jet1_jerC, jet1_jerU;
		TLorentzVector jet2_jecD, jet2_jecU, jet2_jerD, jet2_jerC, jet2_jerU;
		TLorentzVector kinjet1_jecD, kinjet1_jecU, kinjet1_jerD, kinjet1_jerC, kinjet1_jerU;
		TLorentzVector kinjet2_jecD, kinjet2_jecU, kinjet2_jerD, kinjet2_jerC, kinjet2_jerU;
		// FIXME jet systematics for the use of regression still to be implemented
		pho1.SetPtEtaPhiE(t.ph1_pt, t.ph1_eta, t.ph1_phi, t.ph1_e);
		pho2.SetPtEtaPhiE(t.ph2_pt, t.ph2_eta, t.ph2_phi, t.ph2_e);
// the following is the correct thing to call for the day the h2gglobe numbers are to be trusted
/*		pho1_pesD = pho1; pho1_pesD *= t.ph1_pesD_e / t.ph1_e;
		pho2_pesD = pho2; pho2_pesD *= t.ph2_pesD_e / t.ph2_e;
		pho1_pesU = pho1; pho1_pesU *= t.ph1_pesU_e / t.ph1_e;
		pho2_pesU = pho2; pho2_pesU *= t.ph2_pesU_e / t.ph2_e;
		pho1_perD = pho1; pho1_perD *= t.ph1_perD_e / t.ph1_e;
		pho2_perD = pho2; pho2_perD *= t.ph2_perD_e / t.ph2_e;
		pho1_perU = pho1; pho1_perU *= t.ph1_perU_e / t.ph1_e;
		pho2_perU = pho2; pho2_perU *= t.ph2_perU_e / t.ph2_e;*/
// hand-made hard-coded implementation
		pho1_pesD = pho1; pho1_pesD *= (1. - getPESUncertainty(t.ph1_isEB, t.ph1_SCEta, t.ph1_r9));
		pho2_pesD = pho2; pho2_pesD *= (1. - getPESUncertainty(t.ph2_isEB, t.ph2_SCEta, t.ph2_r9));
		pho1_pesU = pho1; pho1_pesU *= (1. + getPESUncertainty(t.ph1_isEB, t.ph1_SCEta, t.ph1_r9));
		pho2_pesU = pho2; pho2_pesU *= (1. + getPESUncertainty(t.ph2_isEB, t.ph2_SCEta, t.ph2_r9));
		pho1_perD = pho1; pho1_perD *= (1. - getPERUncertainty(t.ph1_isEB, t.ph1_SCEta, t.ph1_r9, t.ph1_sigmaEoE, r3));
		pho2_perD = pho2; pho2_perD *= (1. - getPERUncertainty(t.ph2_isEB, t.ph2_SCEta, t.ph2_r9, t.ph2_sigmaEoE, r3));
		pho1_perU = pho1; pho1_perU *= (1. + getPERUncertainty(t.ph1_isEB, t.ph1_SCEta, t.ph1_r9, t.ph1_sigmaEoE, r3));
		pho2_perU = pho2; pho2_perU *= (1. + getPERUncertainty(t.ph2_isEB, t.ph2_SCEta, t.ph2_r9, t.ph2_sigmaEoE, r3));

		jet1.SetPtEtaPhiE(J.jetPt[ij1], J.jetEta[ij1], J.jetPhi[ij1], J.jetE[ij1]);
		jet2.SetPtEtaPhiE(J.jetPt[ij2], J.jetEta[ij2], J.jetPhi[ij2], J.jetE[ij2]);
		jet1_jecD.SetPtEtaPhiE(J.jetJecD_pt[ij1], J.jetJecD_eta[ij1], J.jetJecD_phi[ij1], J.jetJecD_e[ij1]);
		jet2_jecD.SetPtEtaPhiE(J.jetJecD_pt[ij2], J.jetJecD_eta[ij2], J.jetJecD_phi[ij2], J.jetJecD_e[ij2]);
		jet1_jecU.SetPtEtaPhiE(J.jetJecU_pt[ij1], J.jetJecU_eta[ij1], J.jetJecU_phi[ij1], J.jetJecU_e[ij1]);
		jet2_jecU.SetPtEtaPhiE(J.jetJecU_pt[ij2], J.jetJecU_eta[ij2], J.jetJecU_phi[ij2], J.jetJecU_e[ij2]);
		jet1_jerD.SetPtEtaPhiE(J.jetJerD_pt[ij1], J.jetJerD_eta[ij1], J.jetJerD_phi[ij1], J.jetJerD_e[ij1]);
		jet2_jerD.SetPtEtaPhiE(J.jetJerD_pt[ij2], J.jetJerD_eta[ij2], J.jetJerD_phi[ij2], J.jetJerD_e[ij2]);
		jet1_jerC.SetPtEtaPhiE(J.jetJerC_pt[ij1], J.jetJerC_eta[ij1], J.jetJerC_phi[ij1], J.jetJerC_e[ij1]);
		jet2_jerC.SetPtEtaPhiE(J.jetJerC_pt[ij2], J.jetJerC_eta[ij2], J.jetJerC_phi[ij2], J.jetJerC_e[ij2]);
		jet1_jerU.SetPtEtaPhiE(J.jetJerU_pt[ij1], J.jetJerU_eta[ij1], J.jetJerU_phi[ij1], J.jetJerU_e[ij1]);
		jet2_jerU.SetPtEtaPhiE(J.jetJerU_pt[ij2], J.jetJerU_eta[ij2], J.jetJerU_phi[ij2], J.jetJerU_e[ij2]);
		regjet1.SetPtEtaPhiE(J.jetPt[ij1Reg], J.jetEta[ij1Reg], J.jetPhi[ij1Reg], J.jetE[ij1Reg]);
		regjet2.SetPtEtaPhiE(J.jetPt[ij2Reg], J.jetEta[ij2Reg], J.jetPhi[ij2Reg], J.jetE[ij2Reg]);
		regjet1 = ((float)J.jetRegPt[ij1Reg]/(float)J.jetPt[ij1Reg]) * regjet1;
		regjet2 = ((float)J.jetRegPt[ij2Reg]/(float)J.jetPt[ij2Reg]) * regjet2;
		regkinjet1 = regjet1;
		regkinjet2 = regjet2;
		float Hmass = 125.;
		DiJetKinFitter* fitter_jetsH = new DiJetKinFitter( "fitter_jetsH", "fitter_jets", Hmass );
		pair<TLorentzVector,TLorentzVector> jets_kinfitH = fitter_jetsH->fit(regkinjet1, regkinjet2);
		regkinjet1 = jets_kinfitH.first;
		regkinjet2 = jets_kinfitH.second;
		kinjet1 = jet1;
		kinjet2 = jet2;
		kinjet1_jecD = jet1_jecD; kinjet1_jecU = jet1_jecU; kinjet1_jerD = jet1_jerD; kinjet1_jerC = jet1_jerC; kinjet1_jerU = jet1_jerU;
		kinjet2_jecD = jet2_jecD; kinjet2_jecU = jet2_jecU; kinjet2_jerD = jet2_jerD; kinjet2_jerC = jet2_jerC; kinjet2_jerU = jet2_jerU;
		jets_kinfitH = fitter_jetsH->fit(kinjet1, kinjet2);
		kinjet1 = jets_kinfitH.first;
		kinjet2 = jets_kinfitH.second;
		jets_kinfitH = fitter_jetsH->fit(kinjet1_jecD, kinjet2_jecD);
		kinjet1_jecD = jets_kinfitH.first;
		kinjet2_jecD = jets_kinfitH.second;
		jets_kinfitH = fitter_jetsH->fit(kinjet1_jecU, kinjet2_jecU);
		kinjet1_jecU = jets_kinfitH.first;
		kinjet2_jecU = jets_kinfitH.second;
		jets_kinfitH = fitter_jetsH->fit(kinjet1_jerD, kinjet2_jerD);
		kinjet1_jerD = jets_kinfitH.first;
		kinjet2_jerD = jets_kinfitH.second;
		jets_kinfitH = fitter_jetsH->fit(kinjet1_jerC, kinjet2_jerC);
		kinjet1_jerC = jets_kinfitH.first;
		kinjet2_jerC = jets_kinfitH.second;
		jets_kinfitH = fitter_jetsH->fit(kinjet1_jerU, kinjet2_jerU);
		kinjet1_jerU = jets_kinfitH.first;
		kinjet2_jerU = jets_kinfitH.second;

		TLorentzVector jj = jet1 + jet2;
		TLorentzVector regjj = regjet1 + regjet2;
		TLorentzVector regkinjj = regkinjet1 + regkinjet2;
		TLorentzVector kinjj = kinjet1 + kinjet2;
		TLorentzVector gg = pho1 + pho2;
		TLorentzVector ggjj = jj + gg;
		TLorentzVector regggjj = regjj + gg;
		TLorentzVector regkinggjj = regkinjj + gg;
		TLorentzVector kinggjj = kinjj + gg;
		// Photon Energy Scale & Photon Energy Resolution
		TLorentzVector gg_pesD1 = pho1_pesD + pho2;
		TLorentzVector gg_pesD2 = pho1 + pho2_pesD;
		TLorentzVector gg_pesD = pho1_pesD + pho2_pesD;
		TLorentzVector gg_pesU1 = pho1_pesU + pho2;
		TLorentzVector gg_pesU2 = pho1 + pho2_pesU;
		TLorentzVector gg_pesU = pho1_pesU + pho2_pesU;
		TLorentzVector gg_perD1 = pho1_perD + pho2;
		TLorentzVector gg_perD2 = pho1 + pho2_perD;
		TLorentzVector gg_perD = pho1_perD + pho2_perD;
		TLorentzVector gg_perU1 = pho1_perU + pho2;
		TLorentzVector gg_perU2 = pho1 + pho2_perU;
		TLorentzVector gg_perU = pho1_perU + pho2_perU;

		TLorentzVector ggjj_pesD1 = jj + gg_pesD1;
		TLorentzVector ggjj_pesD2 = jj + gg_pesD2;
		TLorentzVector ggjj_pesD = jj + gg_pesD;
		TLorentzVector ggjj_pesU1 = jj + gg_pesU1;
		TLorentzVector ggjj_pesU2 = jj + gg_pesU2;
		TLorentzVector ggjj_pesU = jj + gg_pesU;
		TLorentzVector ggjj_perD1 = jj + gg_perD1;
		TLorentzVector ggjj_perD2 = jj + gg_perD2;
		TLorentzVector ggjj_perD = jj + gg_perD;
		TLorentzVector ggjj_perU1 = jj + gg_perU1;
		TLorentzVector ggjj_perU2 = jj + gg_perU2;
		TLorentzVector ggjj_perU = jj + gg_perU;

		TLorentzVector regggjj_pesD1 = regjj + gg_pesD1;
		TLorentzVector regggjj_pesD2 = regjj + gg_pesD2;
		TLorentzVector regggjj_pesD = regjj + gg_pesD;
		TLorentzVector regggjj_pesU1 = regjj + gg_pesU1;
		TLorentzVector regggjj_pesU2 = regjj + gg_pesU2;
		TLorentzVector regggjj_pesU = regjj + gg_pesU;
		TLorentzVector regggjj_perD1 = regjj + gg_perD1;
		TLorentzVector regggjj_perD2 = regjj + gg_perD2;
		TLorentzVector regggjj_perD = regjj + gg_perD;
		TLorentzVector regggjj_perU1 = regjj + gg_perU1;
		TLorentzVector regggjj_perU2 = regjj + gg_perU2;
		TLorentzVector regggjj_perU = regjj + gg_perU;

		TLorentzVector regkinggjj_pesD1 = regkinjj + gg_pesD1;
		TLorentzVector regkinggjj_pesD2 = regkinjj + gg_pesD2;
		TLorentzVector regkinggjj_pesD = regkinjj + gg_pesD;
		TLorentzVector regkinggjj_pesU1 = regkinjj + gg_pesU1;
		TLorentzVector regkinggjj_pesU2 = regkinjj + gg_pesU2;
		TLorentzVector regkinggjj_pesU = regkinjj + gg_pesU;
		TLorentzVector regkinggjj_perD1 = regkinjj + gg_perD1;
		TLorentzVector regkinggjj_perD2 = regkinjj + gg_perD2;
		TLorentzVector regkinggjj_perD = regkinjj + gg_perD;
		TLorentzVector regkinggjj_perU1 = regkinjj + gg_perU1;
		TLorentzVector regkinggjj_perU2 = regkinjj + gg_perU2;
		TLorentzVector regkinggjj_perU = regkinjj + gg_perU;

		TLorentzVector kinggjj_pesD1 = kinjj + gg_pesD1;
		TLorentzVector kinggjj_pesD2 = kinjj + gg_pesD2;
		TLorentzVector kinggjj_pesD = kinjj + gg_pesD;
		TLorentzVector kinggjj_pesU1 = kinjj + gg_pesU1;
		TLorentzVector kinggjj_pesU2 = kinjj + gg_pesU2;
		TLorentzVector kinggjj_pesU = kinjj + gg_pesU;
		TLorentzVector kinggjj_perD1 = kinjj + gg_perD1;
		TLorentzVector kinggjj_perD2 = kinjj + gg_perD2;
		TLorentzVector kinggjj_perD = kinjj + gg_perD;
		TLorentzVector kinggjj_perU1 = kinjj + gg_perU1;
		TLorentzVector kinggjj_perU2 = kinjj + gg_perU2;
		TLorentzVector kinggjj_perU = kinjj + gg_perU;

		// Jet Energy Correction and Jet Energy Resolution
		TLorentzVector jj_jecD = jet1_jecD + jet2_jecD;
		TLorentzVector kinjj_jecD = kinjet1_jecD + kinjet2_jecD;
		TLorentzVector ggjj_jecD = jj_jecD + gg;
		TLorentzVector kinggjj_jecD = kinjj_jecD + gg;
		TLorentzVector jj_jecU = jet1_jecU + jet2_jecU;
		TLorentzVector kinjj_jecU = kinjet1_jecU + kinjet2_jecU;
		TLorentzVector ggjj_jecU = jj_jecU + gg;
		TLorentzVector kinggjj_jecU = kinjj_jecU + gg;
		TLorentzVector jj_jerD = jet1_jerD + jet2_jerD;
		TLorentzVector kinjj_jerD = kinjet1_jerD + kinjet2_jerD;
		TLorentzVector ggjj_jerD = jj_jerD + gg;
		TLorentzVector kinggjj_jerD = kinjj_jerD + gg;
		TLorentzVector jj_jerC = jet1_jerC + jet2_jerC;
		TLorentzVector kinjj_jerC = kinjet1_jerC + kinjet2_jerC;
		TLorentzVector ggjj_jerC = jj_jerC + gg;
		TLorentzVector kinggjj_jerC = kinjj_jerC + gg;
		TLorentzVector jj_jerU = jet1_jerU + jet2_jerU;
		TLorentzVector kinjj_jerU = kinjet1_jerU + kinjet2_jerU;
		TLorentzVector ggjj_jerU = jj_jerU + gg;
		TLorentzVector kinggjj_jerU = kinjj_jerU + gg;

		t.selection_cut_level = 0;
		t.weight = t.ev_weight;
		t.evweight = t.ev_evweight;
		t.pu_weight = t.ev_pu_weight;
		// adding this for correct yields out of the control plots:
		t.evweight_w_btagSF = t.evweight;
		if( type == -260 ) t.evweight_w_btagSF *= 1.2822;
		if( type  < -250 ) t.evweight_w_btagSF *= eventWeight_2jets("medium", J.jetbtagSF_M[ij1], J.jetbtagSF_M[ij2], J.jetbtagEff_M[ij1], J.jetbtagEff_M[ij2], J.jetCSV[ij1], J.jetCSV[ij2]);
		t.evweight_w_btagSF_reg = t.evweight;
		if( type == -260 ) t.evweight_w_btagSF_reg *= 1.2822;
		if( type  < -250 ) t.evweight_w_btagSF_reg *= eventWeight_2jets("medium", J.jetbtagSF_M[ij1Reg], J.jetbtagSF_M[ij2Reg], J.jetbtagEff_M[ij1Reg], J.jetbtagEff_M[ij2Reg], J.jetCSV[ij1Reg], J.jetCSV[ij2Reg]);
		t.pho1_pt = pho1.Pt();
		t.pho1_e = pho1.E();
		t.pho1_phi = pho1.Phi();
		t.pho1_eta = pho1.Eta();
		t.pho1_mass = pho1.M();
		t.pho1_r9 = t.ph1_r9;
		t.pho1_sieie = t.ph1_sieie;
		t.pho1_hoe = t.ph1_hoe;
		t.pho1_isEB = t.ph1_isEB;
		t.pho1_pfchargedisogood03 = t.ph1_pfchargedisogood03;
		t.pho1_ecaliso = t.ph1_ecaliso;
		t.pho1_pfchargedisobad04 = t.ph1_pfchargedisobad04;
		t.pho1_ecalisobad = t.ph1_ecalisobad;
		t.pho1_badvtx_Et = t.ph1_badvtx_Et;
		t.pho1_PFisoA = (t.ph1_pfchargedisogood03 + t.ph1_ecaliso + 2.5 - t.rho * 0.09) * 50. / t.ph1_pt;
		t.pho1_PFisoB = (t.ph1_pfchargedisobad04 + t.ph1_ecalisobad + 2.5 - t.rho * 0.23) * 50. / t.ph1_badvtx_Et;
		t.pho1_PFisoC = t.ph1_pfchargedisogood03 * 50. / t.ph1_pt;
		t.pho2_pt = pho2.Pt();
		t.pho2_e = pho2.E();
		t.pho2_phi = pho2.Phi();
		t.pho2_eta = pho2.Eta();
		t.pho2_mass = pho2.M();
		t.pho2_r9 = t.ph2_r9;
		t.pho2_sieie = t.ph2_sieie;
		t.pho2_hoe = t.ph2_hoe;
		t.pho2_isEB = t.ph2_isEB;
		t.pho2_pfchargedisogood03 = t.ph2_pfchargedisogood03;
		t.pho2_ecaliso = t.ph2_ecaliso;
		t.pho2_pfchargedisobad04 = t.ph2_pfchargedisobad04;
		t.pho2_ecalisobad = t.ph2_ecalisobad;
		t.pho2_badvtx_Et = t.ph2_badvtx_Et;
		t.pho2_PFisoA = (t.ph2_pfchargedisogood03 + t.ph2_ecaliso + 2.5 - t.rho * 0.09) * 50. / t.ph2_pt;
		t.pho2_PFisoB = (t.ph2_pfchargedisobad04 + t.ph2_ecalisobad + 2.5 - t.rho * 0.23) * 50. / t.ph2_badvtx_Et;
		t.pho2_PFisoC = t.ph2_pfchargedisogood03 * 50. / t.ph2_pt;
		t.jet1_pt = jet1.Pt();
		t.jet1_e = jet1.E();
		t.jet1_phi = jet1.Phi();
		t.jet1_eta = jet1.Eta();
		t.jet1_mass = jet1.M();
		t.jet1_csvBtag = J.jetCSV[ij1];
		t.jet1_btagSF_M = J.jetbtagSF_M[ij1];
		t.jet1_flavour = J.jetflavour[ij1];
		t.jet1_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij1];
		t.jet1_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij1];
		t.jet1_btagEff_M = J.jetbtagEff_M[ij1];
		t.jet1_btagEffError_M = J.jetbtagEffError_M[ij1];
		t.jet1_betaStarClassic = J.jetbetaStarClassic[ij1];
		t.jet1_dR2Mean = J.jetdR2Mean[ij1];
		t.jet2_pt = jet2.Pt();
		t.jet2_e = jet2.E();
		t.jet2_phi = jet2.Phi();
		t.jet2_eta = jet2.Eta();
		t.jet2_mass = jet2.M();
		t.jet2_csvBtag = J.jetCSV[ij2];
		t.jet2_btagSF_M = J.jetbtagSF_M[ij2];
		t.jet2_flavour = J.jetflavour[ij2];
		t.jet2_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij2];
		t.jet2_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij2];
		t.jet2_btagEff_M = J.jetbtagEff_M[ij2];
		t.jet2_btagEffError_M = J.jetbtagEffError_M[ij2];
		t.jet2_betaStarClassic = J.jetbetaStarClassic[ij2];
		t.jet2_dR2Mean = J.jetdR2Mean[ij2];
		t.regjet1_pt = regjet1.Pt();
		t.regjet1_e = regjet1.E();
		t.regjet1_phi = regjet1.Phi();
		t.regjet1_eta = regjet1.Eta();
		t.regjet1_mass = regjet1.M();
		t.regjet1_csvBtag = J.jetCSV[ij1Reg];
		t.regjet1_btagSF_M = J.jetbtagSF_M[ij1Reg];
		t.regjet1_flavour = J.jetflavour[ij1Reg];
		t.regjet1_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij1Reg];
		t.regjet1_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij1Reg];
		t.regjet1_btagEff_M = J.jetbtagEff_M[ij1Reg];
		t.regjet1_btagEffError_M = J.jetbtagEffError_M[ij1Reg];
		t.regjet1_betaStarClassic = J.jetbetaStarClassic[ij1Reg];
		t.regjet1_dR2Mean = J.jetdR2Mean[ij1Reg];
		t.regjet2_pt = regjet2.Pt();
		t.regjet2_e = regjet2.E();
		t.regjet2_phi = regjet2.Phi();
		t.regjet2_eta = regjet2.Eta();
		t.regjet2_mass = regjet2.M();
		t.regjet2_csvBtag = J.jetCSV[ij2Reg];
		t.regjet2_btagSF_M = J.jetbtagSF_M[ij2Reg];
		t.regjet2_flavour = J.jetflavour[ij2Reg];
		t.regjet2_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij2Reg];
		t.regjet2_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij2Reg];
		t.regjet2_btagEff_M = J.jetbtagEff_M[ij2Reg];
		t.regjet2_btagEffError_M = J.jetbtagEffError_M[ij2Reg];
		t.regjet2_betaStarClassic = J.jetbetaStarClassic[ij2Reg];
		t.regjet2_dR2Mean = J.jetdR2Mean[ij2Reg];
		t.regjet1_emfrac = J.jetEmfrac[ij1Reg];
		t.regjet1_hadfrac = J.jetHadfrac[ij1Reg];
		t.regjet1_secVtxPt = J.jetSecVtxPt[ij1Reg];
		t.regjet1_secVtx3dL = J.jetSecVtx3dL[ij1Reg];
		t.regjet1_dPhiMet = J.jetDPhiMet[ij1Reg];
		t.regjet1_nConstituents = J.jetNConstituents[ij1Reg];
		t.regjet2_emfrac = J.jetEmfrac[ij2Reg];
		t.regjet2_hadfrac = J.jetHadfrac[ij2Reg];
		t.regjet2_secVtxPt = J.jetSecVtxPt[ij2Reg];
		t.regjet2_secVtx3dL = J.jetSecVtx3dL[ij2Reg];
		t.regjet2_dPhiMet = J.jetDPhiMet[ij2Reg];
		t.regjet2_nConstituents = J.jetNConstituents[ij2Reg];
		t.regkinjet1_pt = regkinjet1.Pt();
		t.regkinjet1_e = regkinjet1.E();
		t.regkinjet1_phi = regkinjet1.Phi();
		t.regkinjet1_eta = regkinjet1.Eta();
		t.regkinjet1_mass = regkinjet1.M();
		t.regkinjet1_csvBtag = J.jetCSV[ij1RegKin];
		t.regkinjet1_btagSF_M = J.jetbtagSF_M[ij1RegKin];
		t.regkinjet1_flavour = J.jetflavour[ij1RegKin];
		t.regkinjet1_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij1RegKin];
		t.regkinjet1_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij1RegKin];
		t.regkinjet1_btagEff_M = J.jetbtagEff_M[ij1RegKin];
		t.regkinjet1_btagEffError_M = J.jetbtagEffError_M[ij1RegKin];
		t.regkinjet1_betaStarClassic = J.jetbetaStarClassic[ij1RegKin];
		t.regkinjet1_dR2Mean = J.jetdR2Mean[ij1RegKin];
		t.regkinjet2_pt = regkinjet2.Pt();
		t.regkinjet2_e = regkinjet2.E();
		t.regkinjet2_phi = regkinjet2.Phi();
		t.regkinjet2_eta = regkinjet2.Eta();
		t.regkinjet2_mass = regkinjet2.M();
		t.regkinjet2_csvBtag = J.jetCSV[ij2RegKin];
		t.regkinjet2_btagSF_M = J.jetbtagSF_M[ij2RegKin];
		t.regkinjet2_flavour = J.jetflavour[ij2RegKin];
		t.regkinjet2_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij2RegKin];
		t.regkinjet2_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij2RegKin];
		t.regkinjet2_btagEff_M = J.jetbtagEff_M[ij2RegKin];
		t.regkinjet2_btagEffError_M = J.jetbtagEffError_M[ij2RegKin];
		t.regkinjet2_betaStarClassic = J.jetbetaStarClassic[ij2RegKin];
		t.regkinjet2_dR2Mean = J.jetdR2Mean[ij2RegKin];
		t.kinjet1_pt = kinjet1.Pt();
		t.kinjet1_e = kinjet1.E();
		t.kinjet1_phi = kinjet1.Phi();
		t.kinjet1_eta = kinjet1.Eta();
		t.kinjet1_mass = kinjet1.M();
		t.kinjet1_csvBtag = J.jetCSV[ij1];
		t.kinjet1_btagSF_M = J.jetbtagSF_M[ij1];
		t.kinjet1_flavour = J.jetflavour[ij1];
		t.kinjet1_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij1];
		t.kinjet1_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij1];
		t.kinjet1_btagEff_M = J.jetbtagEff_M[ij1];
		t.kinjet1_btagEffError_M = J.jetbtagEffError_M[ij1];
		t.kinjet1_betaStarClassic = J.jetbetaStarClassic[ij1];
		t.kinjet1_dR2Mean = J.jetdR2Mean[ij1];
		t.kinjet2_pt = kinjet2.Pt();
		t.kinjet2_e = kinjet2.E();
		t.kinjet2_phi = kinjet2.Phi();
		t.kinjet2_eta = kinjet2.Eta();
		t.kinjet2_mass = kinjet2.M();
		t.kinjet2_csvBtag = J.jetCSV[ij2];
		t.kinjet2_btagSF_M = J.jetbtagSF_M[ij2];
		t.kinjet2_flavour = J.jetflavour[ij2];
		t.kinjet2_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij2];
		t.kinjet2_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij2];
		t.kinjet2_btagEff_M = J.jetbtagEff_M[ij2];
		t.kinjet2_btagEffError_M = J.jetbtagEffError_M[ij2];
		t.kinjet2_betaStarClassic = J.jetbetaStarClassic[ij2];
		t.kinjet2_dR2Mean = J.jetdR2Mean[ij2];
		t.jj_pt = jj.Pt();
		t.jj_e = jj.E();
		t.jj_phi = jj.Phi();
		t.jj_eta = jj.Eta();
		t.jj_mass = jj.M();
		t.jj_btagSF_M = t.jet1_btagSF_M * t.jet2_btagSF_M;
		t.jj_DR = jet1.DeltaR(jet2);
		t.regjj_pt = regjj.Pt();
		t.regjj_e = regjj.E();
		t.regjj_phi = regjj.Phi();
		t.regjj_eta = regjj.Eta();
		t.regjj_mass = regjj.M();
		t.regjj_btagSF_M = t.regjet1_btagSF_M * t.regjet2_btagSF_M;
		t.regjj_DR = regjet1.DeltaR(regjet2);
		t.regkinjj_pt = regkinjj.Pt();
		t.regkinjj_e = regkinjj.E();
		t.regkinjj_phi = regkinjj.Phi();
		t.regkinjj_eta = regkinjj.Eta();
		t.regkinjj_mass = regkinjj.M();
		t.regkinjj_btagSF_M = t.regkinjet1_btagSF_M * t.regkinjet2_btagSF_M;
		t.regkinjj_DR = regkinjet1.DeltaR(regkinjet2);
		t.kinjj_pt = kinjj.Pt();
		t.kinjj_e = kinjj.E();
		t.kinjj_phi = kinjj.Phi();
		t.kinjj_eta = kinjj.Eta();
		t.kinjj_mass = kinjj.M();
		t.kinjj_btagSF_M = t.kinjet1_btagSF_M * t.kinjet2_btagSF_M;
		t.kinjj_DR = kinjet1.DeltaR(kinjet2);
		t.gg_pt = gg.Pt();
		t.gg_e = gg.E();
		t.gg_phi = gg.Phi();
		t.gg_eta = gg.Eta();
		t.gg_mass = gg.M();
		t.ggjj_pt = ggjj.Pt();
		t.ggjj_e = ggjj.E();
		t.ggjj_phi = ggjj.Phi();
		t.ggjj_eta = ggjj.Eta();
		t.ggjj_mass = ggjj.M();
		t.regggjj_pt = regggjj.Pt();
		t.regggjj_e = regggjj.E();
		t.regggjj_phi = regggjj.Phi();
		t.regggjj_eta = regggjj.Eta();
		t.regggjj_mass = regggjj.M();
		t.regkinggjj_pt = regkinggjj.Pt();
		t.regkinggjj_e = regkinggjj.E();
		t.regkinggjj_phi = regkinggjj.Phi();
		t.regkinggjj_eta = regkinggjj.Eta();
		t.regkinggjj_mass = regkinggjj.M();
		t.kinggjj_pt = kinggjj.Pt();
		t.kinggjj_e = kinggjj.E();
		t.kinggjj_phi = kinggjj.Phi();
		t.kinggjj_eta = kinggjj.Eta();
		t.kinggjj_mass = kinggjj.M();
		t.njets_kLooseID = t.njets_passing_kLooseID;
		t.njets_kLooseID_and_CSVM = t.njets_passing_kLooseID_and_CSVM;
		t.njets_kRadionID = njets_kRadionID_;
		t.njets_kRadionID_and_CSVM = njets_kRadionID_and_CSVM_;
		// Photon Energy Scale & Photon Energy Resolution
		t.gg_mass_pesD1 = gg_pesD1.M();
		t.gg_mass_pesD2 = gg_pesD2.M();
		t.gg_mass_pesD = gg_pesD.M();
		t.gg_mass_pesU1 = gg_pesU1.M();
		t.gg_mass_pesU2 = gg_pesU2.M();
		t.gg_mass_pesU = gg_pesU.M();
		t.gg_mass_perD1 = gg_perD1.M();
		t.gg_mass_perD2 = gg_perD2.M();
		t.gg_mass_perD = gg_perD.M();
		t.gg_mass_perU1 = gg_perU1.M();
		t.gg_mass_perU2 = gg_perU2.M();
		t.gg_mass_perU = gg_perU.M();

		t.ggjj_mass_pesD1 = ggjj_pesD1.M();
		t.ggjj_mass_pesD2 = ggjj_pesD2.M();
		t.ggjj_mass_pesD = ggjj_pesD.M();
		t.ggjj_mass_pesU1 = ggjj_pesU1.M();
		t.ggjj_mass_pesU2 = ggjj_pesU2.M();
		t.ggjj_mass_pesU = ggjj_pesU.M();
		t.ggjj_mass_perD1 = ggjj_perD1.M();
		t.ggjj_mass_perD2 = ggjj_perD2.M();
		t.ggjj_mass_perD = ggjj_perD.M();
		t.ggjj_mass_perU1 = ggjj_perU1.M();
		t.ggjj_mass_perU2 = ggjj_perU2.M();
		t.ggjj_mass_perU = ggjj_perU.M();

		t.regggjj_mass_pesD1 = regggjj_pesD1.M();
		t.regggjj_mass_pesD2 = regggjj_pesD2.M();
		t.regggjj_mass_pesD = regggjj_pesD.M();
		t.regggjj_mass_pesU1 = regggjj_pesU1.M();
		t.regggjj_mass_pesU2 = regggjj_pesU2.M();
		t.regggjj_mass_pesU = regggjj_pesU.M();
		t.regggjj_mass_perD1 = regggjj_perD1.M();
		t.regggjj_mass_perD2 = regggjj_perD2.M();
		t.regggjj_mass_perD = regggjj_perD.M();
		t.regggjj_mass_perU1 = regggjj_perU1.M();
		t.regggjj_mass_perU2 = regggjj_perU2.M();
		t.regggjj_mass_perU = regggjj_perU.M();

		t.regkinggjj_mass_pesD1 = regkinggjj_pesD1.M();
		t.regkinggjj_mass_pesD2 = regkinggjj_pesD2.M();
		t.regkinggjj_mass_pesD = regkinggjj_pesD.M();
		t.regkinggjj_mass_pesU1 = regkinggjj_pesU1.M();
		t.regkinggjj_mass_pesU2 = regkinggjj_pesU2.M();
		t.regkinggjj_mass_pesU = regkinggjj_pesU.M();
		t.regkinggjj_mass_perD1 = regkinggjj_perD1.M();
		t.regkinggjj_mass_perD2 = regkinggjj_perD2.M();
		t.regkinggjj_mass_perD = regkinggjj_perD.M();
		t.regkinggjj_mass_perU1 = regkinggjj_perU1.M();
		t.regkinggjj_mass_perU2 = regkinggjj_perU2.M();
		t.regkinggjj_mass_perU = regkinggjj_perU.M();

		t.kinggjj_mass_pesD1 = kinggjj_pesD1.M();
		t.kinggjj_mass_pesD2 = kinggjj_pesD2.M();
		t.kinggjj_mass_pesD = kinggjj_pesD.M();
		t.kinggjj_mass_pesU1 = kinggjj_pesU1.M();
		t.kinggjj_mass_pesU2 = kinggjj_pesU2.M();
		t.kinggjj_mass_pesU = kinggjj_pesU.M();
		t.kinggjj_mass_perD1 = kinggjj_perD1.M();
		t.kinggjj_mass_perD2 = kinggjj_perD2.M();
		t.kinggjj_mass_perD = kinggjj_perD.M();
		t.kinggjj_mass_perU1 = kinggjj_perU1.M();
		t.kinggjj_mass_perU2 = kinggjj_perU2.M();
		t.kinggjj_mass_perU = kinggjj_perU.M();

		// Jet Energy Correction and Jet Energy Resolution
		t.jj_mass_jecD = jj_jecD.M();
		t.kinjj_mass_jecD = kinjj_jecD.M();
		t.ggjj_mass_jecD = ggjj_jecD.M();
		t.kinggjj_mass_jecD = kinggjj_jecD.M();
		t.jj_mass_jecU = jj_jecU.M();
		t.kinjj_mass_jecU = kinjj_jecU.M();
		t.ggjj_mass_jecU = ggjj_jecU.M();
		t.kinggjj_mass_jecU = kinggjj_jecU.M();
		t.jj_mass_jerD = jj_jerD.M();
		t.kinjj_mass_jerD = kinjj_jerD.M();
		t.ggjj_mass_jerD = ggjj_jerD.M();
		t.kinggjj_mass_jerD = kinggjj_jerD.M();
		t.jj_mass_jerC = jj_jerC.M();
		t.kinjj_mass_jerC = kinjj_jerC.M();
		t.ggjj_mass_jerC = ggjj_jerC.M();
		t.kinggjj_mass_jerC = kinggjj_jerC.M();
		t.jj_mass_jerU = jj_jerU.M();
		t.kinjj_mass_jerU = kinjj_jerU.M();
		t.ggjj_mass_jerU = ggjj_jerU.M();
		t.kinggjj_mass_jerU = kinggjj_jerU.M();
// t.costhetastar
		TLorentzVector Hgg_Rstar(gg);
		TLorentzVector regHgg_Rstar(gg);
		TLorentzVector regkinHgg_Rstar(gg);
		TLorentzVector kinHgg_Rstar(gg);
		Hgg_Rstar.Boost(-ggjj.BoostVector());
		regHgg_Rstar.Boost(-regggjj.BoostVector());
		regkinHgg_Rstar.Boost(-regkinggjj.BoostVector());
		kinHgg_Rstar.Boost(-kinggjj.BoostVector());
		t.costhetastar = Hgg_Rstar.CosTheta();
		t.regcosthetastar = regHgg_Rstar.CosTheta();
		t.regkincosthetastar = regkinHgg_Rstar.CosTheta();
		t.kincosthetastar = kinHgg_Rstar.CosTheta();
// min DR(g, j)
		t.minDRgj = 999999.0;
		t.minDRgregj = 999999.0;
		t.minDRgregkinj = 999999.0;
		t.minDRgkinj = 999999.0;
		t.minDRgj = min(t.minDRgj, (float)pho1.DeltaR(jet1));
		t.minDRgj = min(t.minDRgj, (float)pho1.DeltaR(jet2));
		t.minDRgj = min(t.minDRgj, (float)pho2.DeltaR(jet1));
		t.minDRgj = min(t.minDRgj, (float)pho2.DeltaR(jet2));
		t.minDRgregj = min(t.minDRgregj, (float)pho1.DeltaR(regjet1));
		t.minDRgregj = min(t.minDRgregj, (float)pho1.DeltaR(regjet2));
		t.minDRgregj = min(t.minDRgregj, (float)pho2.DeltaR(regjet1));
		t.minDRgregj = min(t.minDRgregj, (float)pho2.DeltaR(regjet2));
		t.minDRgregkinj = min(t.minDRgregkinj, (float)pho1.DeltaR(regkinjet1));
		t.minDRgregkinj = min(t.minDRgregkinj, (float)pho1.DeltaR(regkinjet2));
		t.minDRgregkinj = min(t.minDRgregkinj, (float)pho2.DeltaR(regkinjet1));
		t.minDRgregkinj = min(t.minDRgregkinj, (float)pho2.DeltaR(regkinjet2));
		t.minDRgkinj = min(t.minDRgkinj, (float)pho1.DeltaR(kinjet1));
		t.minDRgkinj = min(t.minDRgkinj, (float)pho1.DeltaR(kinjet2));
		t.minDRgkinj = min(t.minDRgkinj, (float)pho2.DeltaR(kinjet1));
		t.minDRgkinj = min(t.minDRgkinj, (float)pho2.DeltaR(kinjet2));

		if(REMOVE_UNDEFINED_BTAGSF && (t.regjet1_btagSF_M == -1001 || t.regjet2_btagSF_M == -1001))
		{
			cout << "WARNING: undefined btagSF_M, skipping the t.event:\tevent= " << t.event << "\tregjet1_btagSF_M= " << t.regjet1_btagSF_M << "\tregjet2_btagSF_M= " << t.regjet2_btagSF_M << "\tregjet1_pt= " << t.regjet1_pt << "\tregjet2_pt= " << t.regjet2_pt << endl;
			continue;
			nevents[ilevel]++; eventcut[ilevel] = "After removing undefined btagSF_M";
			nevents_w[ilevel] += t.evweight; nevents_w_btagSF_M[ilevel] += t.evweight * t.regjj_btagSF_M; ilevel++;
		}

// categorisation
		t.selection_cut_level = 0;
		if(SYNC) synchrofile << t.jet1_pt << "\t" << t.jet2_pt << "\t" << t.jj_mass << "\t" << t.ggjj_mass << endl;
		if(t.njets_kRadionID_and_CSVM == 1)
		{
			t.category = 1;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			nevents_sync[8]++;
			eventcut[ilevel] = "1btag category"; ilevel++; ilevel++;
		} else if( t.njets_kRadionID_and_CSVM >=2) {
			ilevel++;
			t.category = 2;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			nevents_sync[9]++;
			eventcut[ilevel] = "2btag category"; ilevel++;
		}



// mjj cut (90/150 for 1btag and 95/140 for 2btags)
		t.selection_cut_level = 2;
		bool pass_mjj = false;
		float min_mjj_1btag, max_mjj_1btag, min_mjj_2btag, max_mjj_2btag;
		min_mjj_1btag = 90.; max_mjj_1btag = 150.; min_mjj_2btag = 95.; max_mjj_2btag = 140.;
		if(SYNC){ min_mjj_1btag = 90.; max_mjj_1btag = 165.; min_mjj_2btag = 95.; max_mjj_2btag = 140.;}
		if(!applyMassCuts){ min_mjj_1btag = 0.; max_mjj_1btag = 14000.; min_mjj_2btag = 0.; max_mjj_2btag = 14000.;}
		pass_mjj = (t.njets_kRadionID_and_CSVM == 1 && (t.regjj_mass < min_mjj_1btag || t.regjj_mass > max_mjj_1btag)) || (t.njets_kRadionID_and_CSVM >= 2 && (t.regjj_mass < min_mjj_2btag || t.regjj_mass > max_mjj_2btag));
		if(pass_mjj)
		{
			outtree->Fill();
			continue;
		}
		if(t.njets_kRadionID_and_CSVM == 1)
		{
			t.category = 1;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			nevents_sync[10]++;
			eventcut[ilevel] = Form("1btag category, after mjj cut (%.1f/%.1f and %.1f/%.1f)",min_mjj_1btag, max_mjj_1btag, min_mjj_2btag, max_mjj_2btag); ilevel++; ilevel++;
		} else if( t.njets_kRadionID_and_CSVM >=2) {
			ilevel++;
			t.category = 2;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			nevents_sync[11]++;
			eventcut[ilevel] = Form("2btag category, after mjj cut (%.1f/%.1f and %.1f/%.1f)",min_mjj_1btag, max_mjj_1btag, min_mjj_2btag, max_mjj_2btag); ilevel++;
		}

// kin fit
// MOVED UPSTREAM, SHOULD BE TRANSPARENT
	
		if(t.njets_kRadionID_and_CSVM == 1)
		{
			t.category = 1;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			eventcut[ilevel] = "1btag category, after kin fit"; ilevel++; ilevel++;
		} else if( t.njets_kRadionID_and_CSVM >=2) {
			ilevel++;
			t.category = 2;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			eventcut[ilevel] = "2btag category, after kin fit"; ilevel++;
		}

// mggjj cut (260/335 and 255/320)
		t.selection_cut_level = 3;
		bool pass_mggjj_cut = false;
		float min_mggjj_1btag, max_mggjj_1btag, min_mggjj_2btag, max_mggjj_2btag;
		min_mggjj_1btag = 260.; max_mggjj_1btag = 335.; min_mggjj_2btag = 255.; max_mggjj_2btag = 320.;
		if(SYNC){ min_mggjj_1btag = 255.; max_mggjj_1btag = 340.; min_mggjj_2btag = 265.; max_mggjj_2btag = 320.;}
		if(!applyMassCuts){ min_mggjj_1btag = 0.; max_mggjj_1btag = 14000.; min_mggjj_2btag = 0.; max_mggjj_2btag = 14000.;}
		pass_mggjj_cut = (t.njets_kRadionID_and_CSVM == 1 && (t.regkinggjj_mass < min_mggjj_1btag || t.regkinggjj_mass > max_mggjj_1btag) ) || (t.njets_kRadionID_and_CSVM >= 2 && (t.regkinggjj_mass < min_mggjj_2btag || t.regkinggjj_mass > max_mggjj_2btag) );
		if( pass_mggjj_cut )
		{
			outtree->Fill();
			continue;
		}
		if(t.njets_kRadionID_and_CSVM == 1)
		{
			t.category = 1;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			nevents_sync[12]++;
			eventcut[ilevel] = Form("1btag category, after mggjj cut (%.1f/%.1f and %.1f/%.1f)", min_mggjj_1btag, max_mggjj_1btag, min_mggjj_2btag, max_mggjj_2btag); ilevel++; ilevel++;
		} else if( t.njets_kRadionID_and_CSVM >=2) {
			ilevel++;
			t.category = 2;
			nevents[ilevel]++;
			nevents_w[ilevel] += t.evweight;
			nevents_sync[13]++;
			eventcut[ilevel] = Form("2btag category, after mggjj cut (%.1f/%.1f and %.1f/%.1f)", min_mggjj_1btag, max_mggjj_1btag, min_mggjj_2btag, max_mggjj_2btag); ilevel++;
		}

		if(SYNC_W_PHIL && (t.event == 1536 || t.event == 1557 || t.event == 1560 || t.event == 7755))
		{
			cout << t.event << endl;
			if(numberOfRegressionFiles == 0)
				cout << t.gg_mass << "\t" <<   t.jj_mass << "\t" << t.ggjj_mass << "\t" <<   t.pho1_pt << "\t" <<   t.pho2_pt << "\t" <<   t.jet1_pt << "\t" <<   t.jet2_pt << "\t" << t.njets_kRadionID_and_CSVM << "\t" <<  t.evweight << "\t" <<  t.pho1_eta << "\t" <<  t.pho1_phi << "\t" <<  t.pho2_eta << "\t" <<  t.pho2_phi << "\t" <<  t.jet1_eta << "\t" <<  t.jet1_phi << "\t" <<  t.jet2_eta << "\t" <<  t.jet2_phi << "\t" << t.vtx_z << "\t" << t.pho1_e << "\t" << t.pho2_e << endl;
			else
				cout << t.gg_mass << "\t" <<   t.regjj_mass << "\t" << t.regggjj_mass << "\t" <<   t.pho1_pt << "\t" <<   t.pho2_pt << "\t" <<   t.regjet1_pt << "\t" <<   t.regjet2_pt << "\t" << t.njets_kRadionID_and_CSVM << "\t" <<  t.evweight << "\t" <<  t.pho1_eta << "\t" <<  t.pho1_phi << "\t" <<  t.pho2_eta << "\t" <<  t.pho2_phi << "\t" <<  t.regjet1_eta << "\t" <<  t.regjet1_phi << "\t" <<  t.regjet2_eta << "\t" <<  t.regjet2_phi << "\t" << t.vtx_z << "\t" << t.pho1_e << "\t" << t.pho2_e << endl;
		}
		if( FULL_DUMP )
			full_dump << "t.run:" << t.run << "\tlumi:" << t.lumis << "\tevent:" << t.event << "\tmGG:" << t.gg_mass << "\tmJJ:" << t.jj_mass << "\tsceta_1:" << 1.0 /*t.ph1_SCEta*/ << "\tsceta_2:" << 1.0 /*t.ph2_SCEta*/ << "\tevcat:" << t.category << "\tdiphoBDT:" << 1.0 << endl;
		ilevelmax=ilevel;
		t.selection_cut_level = 6;
		outtree->Fill();

	} // end of loop over events

	if(DEBUG) cout << "end of loop" << endl;
	if(DEBUG) cout << "ilevelmax= " << ilevelmax << endl;

	for(int i=0 ; i < ilevelmax ; i++)
    cout << "#nevents[" << i << "]= " << nevents[i] << "\t#nevents_w[" << i << "]= " << nevents_w[i] << /*"\t#nevents_w_btagSF_M[" << i << "]= " << nevents_w_btagSF_M[i] << */"\teventcut[" << i << "]= " << eventcut[i] /*<< "\t\t( 1btag= " << nevents1btag[i] << " , 2btag= " << nevents2btag[i] << " ) "*/ << endl;

	if(SYNC)
		for(int i=0 ; i < 14 ; i++)
			cout << nevents_sync[i] << endl;

	if(SYNC) synchrofile.close();
	if(FULL_DUMP) full_dump.close();
  outfile->cd();
  outtree->Write();
  outfile->Close();
  infile->Close();

	if(applyPhotonIDControlSample) cout << "WARNING: you applied the photon ID control sample, please make sure to reweight in (pt, eta) accordingly" << endl;

	return 0;
}

