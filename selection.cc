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
#include <TVector3.h>
#include <TRandom3.h>
#include <TH2F.h>
// Analysis headers
#include "BTagUtils.h"
#include "KinematicFit/DiJetKinFitter.h"
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
    string regressionFilePath;
    int numberOfRegressionFiles;
    int type; // Same conventions as in h2gglobe: <0 = signal ; =0 = data ; >0 = background
    int SYNC; // mjj and mggjj cuts are different for sync and analysis
    int SYNC_W_PHIL;
    int FULL_DUMP;
    int REMOVE_UNDEFINED_BTAGSF;
    int applyMassCuts;
    int applyPhotonIDControlSample;
    int printCutFlow;
    int keep0btag;
    int lambdaReweight;
    int whichPhotonID;
    int iJackknife, nJackknife;
    int FTR14001_style;

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
            ("regressionFilePath", po::value<string>(&regressionFilePath)->default_value("weights/TMVARegression_resonant_BDTG.weights.xml"), "regression file path")
            ("numberOfRegressionFiles,r", po::value<int>(&numberOfRegressionFiles)->default_value(0), "number of regression files")
            ("type", po::value<int>(&type)->default_value(0), "same conventions as in h2gglobe: <0 = signal ; =0 = data ; >0 = background")
            ("sync", po::value<int>(&SYNC)->default_value(0), "mjj and mggjj cuts are overwritten if sync is switched on")
            ("removeUndefinedBtagSF", po::value<int>(&REMOVE_UNDEFINED_BTAGSF)->default_value(0), "remove undefined btagSF_M values (should be used only for the limit trees)")
            ("applyMassCuts", po::value<int>(&applyMassCuts)->default_value(1), "can switch off mass cuts (e.g. for control plots), prevails other mass cut options if switched off")
            ("applyPhotonIDControlSample", po::value<int>(&applyPhotonIDControlSample)->default_value(0), "Invert photon ID CiC cut to populate selection in gjjj instead of ggjj")
            ("sync_w_phil", po::value<int>(&SYNC_W_PHIL)->default_value(0), "switch on output for dedicated events")
            ("full_dump", po::value<int>(&FULL_DUMP)->default_value(0), "switch on creation of the t.event dump")
            ("printCutFlow", po::value<int>(&printCutFlow)->default_value(0), "print cut flow")
            ("keep0btag", po::value<int>(&keep0btag)->default_value(0), "keep 0btag category")
            ("lambdaReweight", po::value<int>(&lambdaReweight)->default_value(-1), "use lambda reweighting (for ggHH sample only)")
            ("whichPhotonID", po::value<int>(&whichPhotonID)->default_value(0), "0= CiC Super Tight, 1= CiC Super Tight with Francois' isolation, 2= photon ID MVA")
            ("nJackknife", po::value<int>(&nJackknife)->default_value(0), "")
            ("iJackknife", po::value<int>(&iJackknife)->default_value(0), "")
            ("FTR14001_style", po::value<int>(&FTR14001_style)->default_value(0), "FTR-14-001 HbbHgg upgrade study cuts")
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

//    string outputtree = inputtree;

    if(DEBUG) cout << "End of argument parsing" << endl;

    if( nJackknife == 1 || (nJackknife >= 2 && iJackknife >= nJackknife) )
    {
        cerr << "Error: you are switching on jacknifing the sample incorrectly:\tnJackknife= " << nJackknife << "\tiJackknife= " << iJackknife << endl;
        return 1;
    }
    if( FTR14001_style == 1 && (whichPhotonID != 0 || numberOfRegressionFiles != 0) )
    {
        cerr << "Error: FTR-14-001 cuts require cut-based photon ID and no regression" << endl;
        return 1;
    }


    cout << "inputfile= " << inputfile << endl;
    cout << "inputtree= " << inputtree << endl;
    cout << "outputfile= " << outputfile << endl;
    cout << "outputtree= " << outputtree << endl;
    cout << "regressionFilePath= " << regressionFilePath << endl;

    TFile *infile = TFile::Open(inputfile.c_str());
    TTree *intree = (TTree*)infile->Get(inputtree.c_str());
    TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
    TTree *outtree = new TTree(outputtree.c_str(), Form("%s reduced", outputtree.c_str()));
    ofstream synchrofile, full_dump, cutFlowFile;
    if(SYNC) synchrofile.open("synchronisation.txt");
    if(FULL_DUMP) full_dump.open(Form("h2gglobe_%s.txt", inputtree.c_str()));
    if(printCutFlow) cutFlowFile.open(Form("cutFlow_%s.dat", outputtree.c_str()));
    // Setup for lambda reweighting for ggHH samples
    TFile *lambdaweightfile = TFile::Open("weights_2D.root");
    TH2F* w_hbb_pt_costhetastar_CS = 0;
    if(lambdaReweight != -1 && type == -2)
    {
        w_hbb_pt_costhetastar_CS = (TH2F*)lambdaweightfile->Get(Form("hFrac_ggHH%i_8TeV_gr_costhetastar_CS_pt_2D", lambdaReweight));
    }

    if(DEBUG) cout << "Setup tree inputs" << endl;
    if(DEBUG) cout << "Setup tree outputs" << endl;
    tree_variables t;
    initialize_variables(&t);

    if(DEBUG) cout << "SetBranchAddresses (inputtree)" << endl;
    setup_intree(intree, &t, type, numberOfRegressionFiles);
    if(DEBUG) cout << "Branch (outtree)" << endl;
    setup_outtree(outtree, &t);

    if(DEBUG) cout << "Setup TRandom3 generator" << endl;
    TRandom3 *r3 = new TRandom3();

    if(DEBUG) cout << "Prepare for regression" << endl;
// prepare for regression
    TMVA::Reader* readerRegres = new TMVA::Reader( "!Color:!Silent" );
    readerRegres->AddVariable( "bjet_pt", &t.jet_pt);
    readerRegres->AddVariable( "bjet_eta", &t.jet_eta);
    readerRegres->AddVariable( "bjet_mt", &t.jet_mt);
    readerRegres->AddVariable( "bjet_phofrac", &t.jet_phofrac);
    readerRegres->AddVariable( "bjet_nhadfrac", &t.jet_nhadfrac);
    readerRegres->AddVariable( "(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptPt) : (-99)", &t.jet_softLeptPt);
    readerRegres->AddVariable( "(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptPtRel) : (-99)", &t.jet_softLeptPtRel);
    readerRegres->AddVariable( "bjet_secVtxM", &t.jet_secVtxM);
    readerRegres->AddVariable( "bjet_secVtx3deL", &t.jet_secVtx3deL);
    readerRegres->AddVariable( "MET", &t.met_corr_pfmet);
    readerRegres->AddVariable( "(abs(bjet_phi-METphi)>3.14159265 ) ? (2*3.14159265-abs(bjet_phi-METphi)) : (abs(bjet_phi-METphi))", &t.jet_dPhiMet);
    readerRegres->AddVariable( "bjet_leadTrackPt", &t.jet_leadTrackPt);
    readerRegres->AddVariable( "bjet_nCharged+bjet_nNeutrals", &t.jet_nConstituents_);
    readerRegres->AddVariable( "rho", &t.rho);
    readerRegres->AddSpectator( "bjet_mufrac", &t.jet_mufrac);
    readerRegres->AddSpectator( "bjet_elefrac", &t.jet_elefrac);
    readerRegres->AddSpectator( "bjet_chadfrac", &t.jet_chadfrac);
    readerRegres->AddSpectator( "(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptDR) : (-99)", &t.jet_softLeptDR);
    readerRegres->AddSpectator( "bjet_secVtxPt", &t.jet_secVtxPt);
    readerRegres->AddSpectator( "bjet_secVtx3dL", &t.jet_secVtx3dL);
    readerRegres->AddSpectator( "bjet_JECUnc", &t.jet_JECUnc);

// Adding variables
    if(numberOfRegressionFiles != 0 && numberOfRegressionFiles != 1)
    {
        cout << "ERROR: current version must have one regression file or no regression" << endl;
        return 1;
    } else {
        for(int i = 0; i < numberOfRegressionFiles ; i++)
        {
            readerRegres->BookMVA("BDTG", Form("%s", regressionFilePath.c_str()));
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

    string flow[30];
    int cutFlow[30] = {0};
    int iflow;

    int totevents = intree->GetEntries();
    int iev0 = 0;
    if(DEBUG) iev0 = 2575;
    if(DEBUG) totevents = iev0+1;
    cout << "#entries= " << totevents << endl;
    // loop over events
    float progress = 0.;
    int k = 0;
    int decade = 0;
    float weight_jackknife = 1.0;
    for(int ievt = iev0 ; ievt < totevents ; ievt++)
    {
        int ilevel = 0;
        if(DEBUG) cout << "#####\tievt= " << ievt << endl;
        progress = 10.0*ievt/(1.0*totevents);
        k = floor(progress);
        if (!DEBUG && k > decade) cout<<10.0*k<<" %"<<endl;
        decade = k;
        // Jackknife stuff :
        if( nJackknife >=2 )
        {
            if( ievt % nJackknife == iJackknife ) continue; // throw away 1 event in nJackknife
            else weight_jackknife = (float)nJackknife/((float)nJackknife - 1.); // to keep consistent weights
        }

        int njets_kRadionID_ = 0;
        int njets_kRadionID_and_CSVM_ = 0;
        intree->GetEntry(ievt);
        if(DEBUG && SYNC_W_PHIL && !(/*t.event == 6976 ||*/ t.event == 8042 || t.event == 14339 /*|| t.event == 2227 || t.event == 4921 || t.event == 7665 || t.event == 7687 || t.event == 11246 || t.event == 15140 || t.event == 685*/) ) continue;
        if(DEBUG) cout << "#####\tievt= " << ievt << "\trun= " << t.run << "\tlumi= " << t.lumis << "\tevent= " << t.event << endl;
    
        if(DEBUG) cout << "for MC, get the MC truth hjj system" << endl;
// Compute hjj system
        TLorentzVector gj1, gj2;
        if( type != 0 && t.gr_j1_p4_pt > .01 && t.gr_j2_p4_pt > .01)
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

        iflow = 0;
        flow[iflow] = "Before photon ID"; cutFlow[iflow]++; iflow++;
        if(DEBUG) cout << "Apply photon ID cuts" << endl;
        // Apply photon ID cuts
        nevents[ilevel]++; eventcut[ilevel] = "Before photon ID";
        nevents_w[ilevel] += t.evweight; ilevel++;
        nevents_sync[0]++;
        t.pho1_PFisoA = (t.ph1_pfchargedisogood03 + t.ph1_ecaliso + 2.5 - t.rho * 0.09) * 50. / t.ph1_pt;
        t.pho1_PFisoB = (t.ph1_pfchargedisobad04 + t.ph1_ecalisobad + 2.5 - t.rho * 0.23) * 50. / t.ph1_badvtx_Et;
        t.pho1_PFisoC = t.ph1_pfchargedisogood03 * 50. / t.ph1_pt;
        t.pho2_PFisoA = (t.ph2_pfchargedisogood03 + t.ph2_ecaliso + 2.5 - t.rho * 0.09) * 50. / t.ph2_pt;
        t.pho2_PFisoB = (t.ph2_pfchargedisobad04 + t.ph2_ecalisobad + 2.5 - t.rho * 0.23) * 50. / t.ph2_badvtx_Et;
        t.pho2_PFisoC = t.ph2_pfchargedisogood03 * 50. / t.ph2_pt;

        if(DEBUG) cout << "t.ph1_pt= " << t.ph1_pt << "\tph1_eta= " << t.ph1_eta << "\tph1_phi= " << t.ph1_phi << "\tph1_r9= " << t.ph1_r9 << "\tph1_SCEta= " << t.ph1_SCEta << "\tph1_isEB= " << t.ph1_isEB << endl;
        if(DEBUG) cout << "t.pho1_PFisoA= " << t.pho1_PFisoA << "\tt.pho1_PFisoB= " << t.pho1_PFisoB << "\tt.pho1_PFisoC= " << t.pho1_PFisoC << "\tph1_sieie= " << t.ph1_sieie << "\tph1_hoe= " << t.ph1_hoe << "\tph1_isconv= " << t.ph1_isconv << "\tt.ph1_r9_cic= " << t.ph1_r9_cic << endl;
        if(DEBUG) cout << "t.ph2_pt= " << t.ph2_pt << "\tph2_eta= " << t.ph2_eta << "\tph2_phi= " << t.ph2_phi << "\tph2_r9= " << t.ph2_r9 << "\tph2_SCEta= " << t.ph2_SCEta << "\tph2_isEB= " << t.ph2_isEB << endl;
        if(DEBUG) cout << "t.pho2_PFisoA= " << t.pho2_PFisoA << "\tt.pho2_PFisoB= " << t.pho2_PFisoB << "\tt.pho2_PFisoC= " << t.pho2_PFisoC << "\tph2_sieie= " << t.ph2_sieie << "\tph2_hoe= " << t.ph2_hoe << "\tph2_isconv= " << t.ph2_isconv << "\tt.ph2_r9_cic= " << t.ph2_r9_cic << endl;

        if(DEBUG) cout << "t.ph1_pt= " << t.ph1_pt << "\t(float)(40.*t.PhotonsMass)/(float)120.= " << (float)(40.*t.PhotonsMass)/(float)120. << endl;
        if( (!FTR14001_style) ? (t.ph1_pt < (float)(40.*t.PhotonsMass)/(float)120.) : (t.ph1_pt < 40.) ) continue;
        nevents[ilevel]++; eventcut[ilevel] = "After floating pt cut for photon 1 (40*mgg/120 GeV)";
        nevents_w[ilevel] += t.evweight; ilevel++;
//        if( t.ph2_pt < 25. ) continue;
        if(DEBUG) cout << "t.ph2_pt= " << t.ph2_pt << "\t(float)(30.*t.PhotonsMass)/(float)120.= " << (float)(30.*t.PhotonsMass)/(float)120. << endl;
        if( (!FTR14001_style) ? (t.ph2_pt < (float)(30.*t.PhotonsMass)/(float)120.) : (t.ph2_pt < 25.) ) continue; // switching to running pt cut per Hgg recommendations (Nov. 2013)
        nevents[ilevel]++; eventcut[ilevel] = "After fixed pt cut for photon 2 (25 GeV)";
        nevents_w[ilevel] += t.evweight; ilevel++;
        nevents_sync[1]++;
        flow[iflow] = "After photon pt cuts"; cutFlow[iflow]++; iflow++;

        if(DEBUG) cout << "t.ph1_ciclevel= " << t.ph1_ciclevel << "\tph2_ciclevel= " << t.ph2_ciclevel << endl;
        if(!applyPhotonIDControlSample)
        {
            if( whichPhotonID == 0 )
            {
                if ((t.ph1_ciclevel < 4) || (t.ph2_ciclevel < 4))  continue;
//                if ((t.ph1_ciclevel < 0) || (t.ph2_ciclevel < 0))  continue;
            }
            else if( whichPhotonID == 1 )
            {
            // switches for CIC 4 vs CIC 1 iso A and B values
            // for hand-made CiC level, by default switching off isolations on trailing photon
                bool noIsoA1 = false; bool noIsoA2 = true;
                bool noIsoB1 = false; bool noIsoB2 = true;
                bool pho1_cic4 = false;
                bool pho2_cic4 = false;
                bool pho1_cic0 = false;
                bool pho2_cic0 = false;
                getHandMadeCiCLevel(&pho1_cic4, &pho2_cic4, &pho1_cic0, &pho2_cic0, &t,  noIsoA1, noIsoA2, noIsoB1, noIsoB2);
                if(DEBUG) cout << "t.ph1_ciclevel= " << t.ph1_ciclevel << "\tt.ph2_ciclevel= " << t.ph2_ciclevel << endl;
                if(DEBUG) cout << "pho1_cic4= " << pho1_cic4 << "\tpho2_cic4= " << pho2_cic4 << "\tpho1_cic0= " << pho1_cic0 << "\tpho2_cic0= " << pho2_cic0 << endl;
                if( !(pho1_cic4) || !(pho2_cic4) ) continue;

            }
            else if( whichPhotonID == 2 )
            {
                bool ph1_id = t.ph1_isEB ? (t.ph1_IDmva > 0.02) : (t.ph1_IDmva > 0.1);
                bool ph2_id = t.ph2_isEB ? (t.ph2_IDmva > 0.02) : (t.ph2_IDmva > 0.1);
                if(!( ph1_id && ph2_id )) continue; 
            }
            else
            {
                cout << "Erm, no valid photon ID was asked for, crashing now for safety reasons" << endl;
                return 3;
            }
        }
        else if (applyPhotonIDControlSample)
        {
            if( whichPhotonID == 0 )
            {
                bool ph1_id = (t.ph1_ciclevel >= 4);
                bool ph2_id = (t.ph2_ciclevel >= 4);
                bool ph1_Lid = (t.ph1_ciclevel >= 0) && (t.ph1_ciclevel < 4);
                bool ph2_Lid = (t.ph2_ciclevel >= 0) && (t.ph2_ciclevel < 4);
                if(ph1_id && ph2_id) continue; // reject gg
                if(ph1_Lid && ph2_Lid) continue; // reject jj
                if( !( (ph1_id && ph2_Lid) || (ph2_id && ph1_Lid)) ) continue; // reject if different from gj or jg
            }
            else if( whichPhotonID == 1 )
            {
                bool noIsoA1 = false; bool noIsoA2 = true;
                bool noIsoB1 = false; bool noIsoB2 = true;
                bool pho1_cic4 = false;
                bool pho2_cic4 = false;
                bool pho1_cic0 = false;
                bool pho2_cic0 = false;
                getHandMadeCiCLevel(&pho1_cic4, &pho2_cic4, &pho1_cic0, &pho2_cic0, &t,  noIsoA1, noIsoA2, noIsoB1, noIsoB2);
                bool ph1_id = pho1_cic4;
                bool ph2_id = pho2_cic4;
                bool ph1_Lid = (pho1_cic0) && !(pho1_cic4);
                bool ph2_Lid = (pho2_cic0) && !(pho2_cic4);
//                assert( pho1_cic4 == (t.ph1_ciclevel >= 4) );
//                assert( pho2_cic4 == (t.ph2_ciclevel >= 4) );
//                assert( pho1_cic0 == (t.ph1_ciclevel >= 0) );
//                assert( pho2_cic0 == (t.ph2_ciclevel >= 0) );
                if(ph1_id && ph2_id) continue; // reject gg
                if(ph1_Lid && ph2_Lid) continue; // reject jj
                if( !( (ph1_id && ph2_Lid) || (ph2_id && ph1_Lid)) ) continue; // reject if different from gj or jg
            }
            else if( whichPhotonID == 2 )
            {
                bool ph1_id = t.ph1_isEB ? (t.ph1_IDmva > 0.02) : (t.ph1_IDmva > 0.1);
                bool ph2_id = t.ph2_isEB ? (t.ph2_IDmva > 0.02) : (t.ph2_IDmva > 0.1);
                bool ph1_Lid = (t.ph1_ciclevel >= 0) && (!ph1_id);
                bool ph2_Lid = (t.ph2_ciclevel >= 0) && (!ph2_id);
                if(ph1_id && ph2_id) continue; // reject gg
                if(ph1_Lid && ph2_Lid) continue; // reject jj
                if( !( (ph1_id && ph2_Lid) || (ph2_id && ph1_Lid)) ) continue; // reject if different from gj or jg
            }
            else
            {
                cout << "Erm, no valid photon ID was asked for, crashing now for safety reasons" << endl;
                return 3;
            }
        }
        flow[iflow] = "After photon cic ID"; cutFlow[iflow]++; iflow++;
        nevents[ilevel]++; eventcut[ilevel] = "After cic cut on both photons";
        nevents_w[ilevel] += t.evweight; ilevel++;
        nevents_sync[2]++;
        if(DEBUG) cout << "t.PhotonsMass= " << t.PhotonsMass << endl;
        if( (!FTR14001_style) ? ((t.PhotonsMass < 100.) || (t.PhotonsMass > 180.)) : ((t.PhotonsMass < 100.) || (t.PhotonsMass > 150.)) ) continue;
        flow[iflow] = "After diphoton mass cut"; cutFlow[iflow]++; iflow++;
        nevents[ilevel]++; eventcut[ilevel] = "After 100 < mgg < 180";
        nevents_w[ilevel] += t.evweight; ilevel++;
        nevents_sync[3]++;
        t.HT_gg = t.ph1_pt + t.ph2_pt;

        // take only the subset of events where at least two jets remains
        if(DEBUG) cout << "t.njets_passing_kLooseID= " << t.njets_passing_kLooseID << endl;
        if( t.njets_passing_kLooseID < 2 ) continue;
        flow[iflow] = "After requirement at least two jets in the tree"; cutFlow[iflow]++; iflow++;
        nevents[ilevel]++; eventcut[ilevel] = "After njet >= 2";
        nevents_w[ilevel] += t.evweight; ilevel++;
        nevents_sync[4]++;
// alternative counting: taking into account only the 4 jets stored !
        int nbjet_tmp = 0;
//        float csv_cut = 0.244; // CSVL
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
        if( (!keep0btag) && nbjet_tmp < 1 ) continue;
        if( FTR14001_style && nbjet_tmp < 2 ) continue; // need at least 2 bjets for FTR14001_style
        flow[iflow] = "After at least one CSVM jet"; cutFlow[iflow]++; iflow++;
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
            if( ! fill_jet_variables( &t, ijet, met, numberOfRegressionFiles) ) break; // don't loose time if this jet is not filled with interesting information
            t.jet_nConstituents_ = (float) t.jet_nConstituents;
            t.jet_dPhiMet_fabs = fabs(t.jet_dPhiMet);

            if(DEBUG && numberOfRegressionFiles != 0) cout << "input= " << t.jet_pt << "\toutput (BDTG)= " << readerRegres->EvaluateRegression("BDTG")[0] << endl;
            if( REMOVE_UNDEFINED_BTAGSF && t.jet_flavour == 0. ) continue;
//            njets[5]++; jetcut[5] = "After t.jet_csvBtag < 0.";
            if(DEBUG) cout << "now with the regression" << endl;
            if(numberOfRegressionFiles == 0)
                t.jet_regPt = t.jet_pt; // no regression applied
            else if(numberOfRegressionFiles == 1)
                t.jet_regPt = (float)(readerRegres->EvaluateRegression("BDTG")[0]);
            t.jet_regkinPt = t.jet_pt; // kinfit is applied on non-regressed jet 4-momenta
            // jet selection
            // ** acceptance + pu id **
            if( (!FTR14001_style) ? (t.jet_regPt < 25.) : (t.jet_regPt < 30.) ) continue;
            njets[1]++; jetcut[1] = "After jet pt > 25";
            if( (!FTR14001_style) ? (fabs(t.jet_eta) > 2.5) : (fabs(t.jet_eta) > 2.4) ) continue;
            njets[2]++; jetcut[2] = "After jet |eta| < 2.5";
//            if( t.jet_betaStarClassic > 0.2 * log( t.nvtx - 0.64) ) continue;
            njets[3]++; jetcut[3] = "After t.jet_betaStarClassic > 0.2 * log( t.nvtx - 0.64)";
//            if( t.jet_dR2Mean > 0.06 ) continue;
            njets[4]++; jetcut[4] = "After t.jet_dR2Mean > 0.06";
            // ** call regression to correct the pt **
            TLorentzVector tmp_jet, pho1, pho2;
            tmp_jet.SetPtEtaPhiE(t.jet_pt, t.jet_eta, t.jet_phi, t.jet_e);
            pho1.SetPtEtaPhiE(t.ph1_pt, t.ph1_eta, t.ph1_phi, t.ph1_e);
            pho2.SetPtEtaPhiE(t.ph2_pt, t.ph2_eta, t.ph2_phi, t.ph2_e);
            if( (!FTR14001_style) ? ((tmp_jet.DeltaR(pho1) < .5) || (tmp_jet.DeltaR(pho2) < .5)) : ((tmp_jet.DeltaR(pho1) < 1.5) || (tmp_jet.DeltaR(pho2) < 1.5)) ) continue;
            njets[5]++; jetcut[5] = "After minDR(g,j) > .5";
            if(DEBUG) cout << "Jet is passing selection cuts" << endl;

            // ** store 4-momentum + csv output for combinatorics **
            J.jetPt.push_back(t.jet_pt);
            J.jetbtagSF_M.push_back(t.jet_btagSF_M);
            J.jetflavour.push_back(t.jet_flavour);
            J.jetbtagSFErrorUp_M.push_back(t.jet_btagSFErrorUp_M);
            J.jetbtagSFErrorDown_M.push_back(t.jet_btagSFErrorDown_M);
            J.jetbtagEff_M.push_back(t.jet_btagEff_M);
            J.jetbtagEffError_M.push_back(t.jet_btagEffError_M);
            J.jetcutbased_wp_level.push_back(t.jet_cutbased_wp_level);
            J.jetsimple_wp_level.push_back(t.jet_simple_wp_level);
            J.jetfull_wp_level.push_back(t.jet_full_wp_level);
            J.jetdR2Mean.push_back(t.jet_dR2Mean);
            J.jetbetaStarClassic.push_back(t.jet_betaStarClassic);
            J.jetE.push_back(t.jet_e);
            J.jetEta.push_back(t.jet_eta);
            J.jetPhi.push_back(t.jet_phi);
            J.jetCSV.push_back(t.jet_csvBtag);
            J.jetRegPt.push_back(t.jet_regPt);
            J.jetRegKinPt.push_back(t.jet_regkinPt);
            J.jetMt.push_back(t.jet_mt);
            J.jetChadfrac.push_back(t.jet_chadfrac);
            J.jetNhadfrac.push_back(t.jet_nhadfrac);
            J.jetPhofrac.push_back(t.jet_phofrac);
            J.jetMufrac.push_back(t.jet_mufrac);
            J.jetElefrac.push_back(t.jet_elefrac);
            J.jetSoftLeptPt.push_back(t.jet_softLeptPt);
            J.jetSoftLeptPtRel.push_back(t.jet_softLeptPtRel);
            J.jetSoftLeptDR.push_back(t.jet_softLeptDR);
            J.jetJECUnc.push_back(t.jet_JECUnc);
            J.jetLeadTrackPt.push_back(t.jet_leadTrackPt);
            J.jetSecVtxPt.push_back(t.jet_secVtxPt);
            J.jetSecVtx3dL.push_back(t.jet_secVtx3dL);
            J.jetSecVtx3deL.push_back(t.jet_secVtx3deL);
            J.jetSecVtxM.push_back(t.jet_secVtxM);
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
        flow[iflow] = "After at least two jets passing the jet selection"; cutFlow[iflow]++; iflow++;
        nevents_w[ilevel] += t.evweight; ilevel++;
        nevents_sync[6]++;

        vector<int> btaggedJet;
        for( unsigned int ijet = 0 ; ijet < J.jetPt.size() ; ijet++ )
        {
            if( J.jetCSV[ijet] > csv_cut )
                btaggedJet.push_back(ijet);
        }

        if( (!keep0btag) && btaggedJet.size() < 1 ) continue;
        if( FTR14001_style && btaggedJet.size() < 2 ) continue; // need at least 2 bjets for FTR14001_style
        if( FTR14001_style && J.jetPt.size() >= 4 ) continue; // need less than 4 jets for FTR14001_style
        // for FTR14001_style there is additional cuts needed : 
        // - lepton veto & DR(electron, photon) > .1 (not possible with the current globe trees)
        // - DR(g,g) < 2.
        // - DR(j,j) < 2.
        // - 70 < mjj < 200
        // Given these depend on the jet combinatorics and that these may change, they will be applied at quicktree level
        nevents[ilevel]++; eventcut[ilevel] = "After nbjet >=1 passing the jet selection";
        flow[iflow] = "After at least one CSVM jet passing the jet selection"; cutFlow[iflow]++; iflow++;
        nevents_w[ilevel] += t.evweight; ilevel++;
        nevents_sync[7]++;
        if( btaggedJet.size() == 1) nevents1btag[12]++; 
        else nevents2btag[12]++;



        int ij1 = 0;
        int ij2 = 0;
        int ij1Reg = 0;
        int ij2Reg = 0;
 
        if(DEBUG) cout << "btaggedJet.size()= " << btaggedJet.size() << endl;
        if(DEBUG)
            for(int ijet_=0; ijet_ < (int)btaggedJet.size() ; ijet_++)
                cout << "J.jetPt[btaggedJet[" << ijet_ << "]]= " << J.jetPt[btaggedJet[ijet_]] << endl;
        // if exactly 0 btag, find the jet pair that gives max ptjj
        if( btaggedJet.size() == 0 )
        {
            if(DEBUG) cout << "Entering jet combinatorics: 0btag category" << endl;
            t.category = 0;
            unsigned int ij = 0;
//            if(DEBUG) cout << "btaggedJet[0]= " << btaggedJet[0] << endl;
            TLorentzVector j, jreg;
            j.SetPtEtaPhiE(J.jetPt[ij], J.jetEta[ij], J.jetPhi[ij], J.jetE[ij]);
            jreg = ((float)J.jetRegPt[ij]/(float)J.jetPt[ij]) * j;
            int imaxptjj;
            int imaxptjjReg;
            float maxptjj = -99.;
            float maxptjjReg = -99.;
            for(unsigned int ijet = 0 ; ijet < J.jetPt.size() ; ijet++)
            {
                if( ijet == ij ) continue;
                TLorentzVector tmp_j;
                TLorentzVector tmp_jReg;
                tmp_j.SetPtEtaPhiE(J.jetPt[ijet], J.jetEta[ijet], J.jetPhi[ijet], J.jetE[ijet]);
                tmp_jReg = ((float)J.jetRegPt[ijet]/(float)J.jetPt[ijet]) * tmp_j;
                TLorentzVector jj = j + tmp_j;
                TLorentzVector jjReg = jreg + tmp_jReg;
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
            }
            ij1 = ij;
            ij2 = imaxptjj;
            ij1Reg = ij;
            ij2Reg = imaxptjjReg;
        }
        if( btaggedJet.size() == 1 )
        {
            if(DEBUG) cout << "Entering jet combinatorics: 1btag category" << endl;
            t.category = 1;
            unsigned int ij = btaggedJet[0];
            if(DEBUG) cout << "btaggedJet[0]= " << btaggedJet[0] << endl;
            TLorentzVector j, jreg;
            j.SetPtEtaPhiE(J.jetPt[ij], J.jetEta[ij], J.jetPhi[ij], J.jetE[ij]);
            jreg = ((float)J.jetRegPt[ij]/(float)J.jetPt[ij]) * j;
            int imaxptjj;
            int imaxptjjReg;
            float maxptjj = -99.;
            float maxptjjReg = -99.;
            for(unsigned int ijet = 0 ; ijet < J.jetPt.size() ; ijet++)
            {
                if( ijet == ij ) continue;
                TLorentzVector tmp_j;
                TLorentzVector tmp_jReg;
                tmp_j.SetPtEtaPhiE(J.jetPt[ijet], J.jetEta[ijet], J.jetPhi[ijet], J.jetE[ijet]);
                tmp_jReg = ((float)J.jetRegPt[ijet]/(float)J.jetPt[ijet]) * tmp_j;
                TLorentzVector jj = j + tmp_j;
                TLorentzVector jjReg = jreg + tmp_jReg;
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
            }
            ij1 = ij;
            ij2 = imaxptjj;
            ij1Reg = ij;
            ij2Reg = imaxptjjReg;
        }
        // if two or more bjets, then loop only over btagged jets
        if( btaggedJet.size() > 1 )
        {
            t.category = 2;
            if(DEBUG) cout << "Entering jet combinatorics: 2btag category" << endl;
            int ij;
            int imaxptjj = 0;
            int imaxptjjReg = 0;
            int jmaxptjj = 0;
            int jmaxptjjReg = 0;
            float maxptjj = -99.;
            float maxptjjReg = -99.;
            for( unsigned int i = 0 ; i < btaggedJet.size() - 1 ; i++ )
            {
                ij = btaggedJet[i];
                if(DEBUG) cout << "btaggedJet[" << i << "]= " << ij << "\tjetPt[" << ij << "]= " << J.jetPt[ij] << endl;
                TLorentzVector j, jreg, jregkin;
                j.SetPtEtaPhiE(J.jetPt[ij], J.jetEta[ij], J.jetPhi[ij], J.jetE[ij]);
                jreg = ((float)J.jetRegPt[ij]/(float)J.jetPt[ij]) * j;
                for(unsigned int k = i+1 ; k < btaggedJet.size() ; k++)
                {
                    int ijet = btaggedJet[k];
                    TLorentzVector tmp_j;
                    TLorentzVector tmp_jReg;
                    tmp_j.SetPtEtaPhiE(J.jetPt[ijet], J.jetEta[ijet], J.jetPhi[ijet], J.jetE[ijet]);
                    tmp_jReg = ((float)J.jetRegPt[ijet]/(float)J.jetPt[ijet]) * tmp_j;
                    TLorentzVector jj = j + tmp_j;
                    TLorentzVector jjReg = jreg + tmp_jReg;
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
                }
            } // end of loop over bjets
            ij1 = imaxptjj;
            ij2 = jmaxptjj;
            ij1Reg = imaxptjjReg;
            ij2Reg = jmaxptjjReg;
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
/*        pho1_pesD = pho1; pho1_pesD *= t.ph1_pesD_e / t.ph1_e;
        pho2_pesD = pho2; pho2_pesD *= t.ph2_pesD_e / t.ph2_e;
        pho1_pesU = pho1; pho1_pesU *= t.ph1_pesU_e / t.ph1_e;
        pho2_pesU = pho2; pho2_pesU *= t.ph2_pesU_e / t.ph2_e;
        pho1_perD = pho1; pho1_perD *= t.ph1_perD_e / t.ph1_e;
        pho2_perD = pho2; pho2_perD *= t.ph2_perD_e / t.ph2_e;
        pho1_perU = pho1; pho1_perU *= t.ph1_perU_e / t.ph1_e;
        pho2_perU = pho2; pho2_perU *= t.ph2_perU_e / t.ph2_e;*/
// hand-made hard-coded implementation
        pho1_pesD = pho1; pho1_pesD *= (1. - getPESUncertainty(t.ph1_isEB, t.ph1_SCEta, t.ph1_r9, t.ph1_pt));
        pho2_pesD = pho2; pho2_pesD *= (1. - getPESUncertainty(t.ph2_isEB, t.ph2_SCEta, t.ph2_r9, t.ph2_pt));
        pho1_pesU = pho1; pho1_pesU *= (1. + getPESUncertainty(t.ph1_isEB, t.ph1_SCEta, t.ph1_r9, t.ph1_pt));
        pho2_pesU = pho2; pho2_pesU *= (1. + getPESUncertainty(t.ph2_isEB, t.ph2_SCEta, t.ph2_r9, t.ph2_pt));
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
        regjet1.SetPtEtaPhiE(J.jetRegPt[ij1Reg], J.jetEta[ij1Reg], J.jetPhi[ij1Reg], J.jetE[ij1Reg]*((float)J.jetRegPt[ij1Reg]/(float)J.jetPt[ij1Reg]));
        regjet2.SetPtEtaPhiE(J.jetRegPt[ij2Reg], J.jetEta[ij2Reg], J.jetPhi[ij2Reg], J.jetE[ij2Reg]*((float)J.jetRegPt[ij2Reg]/(float)J.jetPt[ij2Reg]));
//        regkinjet1 = regjet1;
//        regkinjet2 = regjet2;
//      kinfit is applied on non-regressed jet 4-momentum, but with regression combinatorics
        regkinjet1.SetPtEtaPhiE(J.jetRegKinPt[ij1Reg], J.jetEta[ij1Reg], J.jetPhi[ij1Reg], J.jetE[ij1Reg]);
        regkinjet2.SetPtEtaPhiE(J.jetRegKinPt[ij2Reg], J.jetEta[ij2Reg], J.jetPhi[ij2Reg], J.jetE[ij2Reg]);
        float Hmass_MC = 125.;
        float Hmass_data = 125.03; // From https://twiki.cern.ch/twiki/bin/view/CMSPublic/Hig14009TWiki
//        float Hmass_data = 125.6; // From https://twiki.cern.ch/twiki/bin/view/CMSPublic/Hig13002PubTWiki
        float Hmass = (type==0) ? Hmass_data : Hmass_MC;
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
       
        t.pho1_pt = pho1.Pt();
        t.pho1_e = pho1.E();
        t.pho1_phi = pho1.Phi();
        t.pho1_eta = pho1.Eta();
        t.pho1_mass = pho1.M();
        t.pho1_r9 = t.ph1_r9;
        t.pho1_r9_cic = t.ph1_r9_cic;
        t.pho1_IDmva = t.ph1_IDmva;
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
        t.pho2_r9_cic = t.ph2_r9_cic;
        t.pho2_IDmva = t.ph2_IDmva;
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
        t.jet1_cutbased_wp_level = J.jetcutbased_wp_level[ij1];
        t.jet1_simple_wp_level = J.jetsimple_wp_level[ij1];
        t.jet1_full_wp_level = J.jetfull_wp_level[ij1];
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
        t.jet2_cutbased_wp_level = J.jetcutbased_wp_level[ij2];
        t.jet2_simple_wp_level = J.jetsimple_wp_level[ij2];
        t.jet2_full_wp_level = J.jetfull_wp_level[ij2];
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
        t.regjet1_cutbased_wp_level = J.jetcutbased_wp_level[ij1Reg];
        t.regjet1_simple_wp_level = J.jetsimple_wp_level[ij1Reg];
        t.regjet1_full_wp_level = J.jetfull_wp_level[ij1Reg];
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
        t.regjet2_cutbased_wp_level = J.jetcutbased_wp_level[ij2Reg];
        t.regjet2_simple_wp_level = J.jetsimple_wp_level[ij2Reg];
        t.regjet2_full_wp_level = J.jetfull_wp_level[ij2Reg];
        t.regjet2_betaStarClassic = J.jetbetaStarClassic[ij2Reg];
        t.regjet2_dR2Mean = J.jetdR2Mean[ij2Reg];
        t.regjet1_mt = J.jetMt[ij1Reg];
        t.regjet1_chadfrac = J.jetChadfrac[ij1Reg];
        t.regjet1_nhadfrac = J.jetNhadfrac[ij1Reg];
        t.regjet1_phofrac = J.jetPhofrac[ij1Reg];
        t.regjet1_mufrac = J.jetMufrac[ij1Reg];
        t.regjet1_elefrac = J.jetElefrac[ij1Reg];
        t.regjet1_softLeptPt = J.jetSoftLeptPt[ij1Reg];
        t.regjet1_softLeptPtRel = J.jetSoftLeptPtRel[ij1Reg];
        t.regjet1_softLeptDR = J.jetSoftLeptDR[ij1Reg];
        t.regjet1_leadTrackPt = J.jetLeadTrackPt[ij1Reg];
        t.regjet1_JECUnc = J.jetJECUnc[ij1Reg];
        t.regjet1_secVtxPt = J.jetSecVtxPt[ij1Reg];
        t.regjet1_secVtx3dL = J.jetSecVtx3dL[ij1Reg];
        t.regjet1_secVtx3deL = J.jetSecVtx3deL[ij1Reg];
        t.regjet1_secVtxM = J.jetSecVtxM[ij1Reg];
        t.regjet1_dPhiMet = J.jetDPhiMet[ij1Reg];
        t.regjet1_nConstituents = J.jetNConstituents[ij1Reg];
        t.regjet2_mt = J.jetMt[ij2Reg];
        t.regjet2_chadfrac = J.jetChadfrac[ij2Reg];
        t.regjet2_nhadfrac = J.jetNhadfrac[ij2Reg];
        t.regjet2_phofrac = J.jetPhofrac[ij2Reg];
        t.regjet2_mufrac = J.jetMufrac[ij2Reg];
        t.regjet2_elefrac = J.jetElefrac[ij2Reg];
        t.regjet2_softLeptPt = J.jetSoftLeptPt[ij2Reg];
        t.regjet2_softLeptPtRel = J.jetSoftLeptPtRel[ij2Reg];
        t.regjet2_softLeptDR = J.jetSoftLeptDR[ij2Reg];
        t.regjet2_leadTrackPt = J.jetLeadTrackPt[ij2Reg];
        t.regjet2_JECUnc = J.jetJECUnc[ij2Reg];
        t.regjet2_secVtxPt = J.jetSecVtxPt[ij2Reg];
        t.regjet2_secVtx3dL = J.jetSecVtx3dL[ij2Reg];
        t.regjet2_secVtx3deL = J.jetSecVtx3deL[ij2Reg];
        t.regjet2_secVtxM = J.jetSecVtxM[ij2Reg];
        t.regjet2_dPhiMet = J.jetDPhiMet[ij2Reg];
        t.regjet2_nConstituents = J.jetNConstituents[ij2Reg];
        t.regkinjet1_pt = regkinjet1.Pt();
        t.regkinjet1_e = regkinjet1.E();
        t.regkinjet1_phi = regkinjet1.Phi();
        t.regkinjet1_eta = regkinjet1.Eta();
        t.regkinjet1_mass = regkinjet1.M();
        t.regkinjet1_csvBtag = J.jetCSV[ij1Reg];
        t.regkinjet1_btagSF_M = J.jetbtagSF_M[ij1Reg];
        t.regkinjet1_flavour = J.jetflavour[ij1Reg];
        t.regkinjet1_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij1Reg];
        t.regkinjet1_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij1Reg];
        t.regkinjet1_btagEff_M = J.jetbtagEff_M[ij1Reg];
        t.regkinjet1_btagEffError_M = J.jetbtagEffError_M[ij1Reg];
        t.regkinjet1_cutbased_wp_level = J.jetcutbased_wp_level[ij1Reg];
        t.regkinjet1_simple_wp_level = J.jetsimple_wp_level[ij1Reg];
        t.regkinjet1_full_wp_level = J.jetfull_wp_level[ij1Reg];
        t.regkinjet1_betaStarClassic = J.jetbetaStarClassic[ij1Reg];
        t.regkinjet1_dR2Mean = J.jetdR2Mean[ij1Reg];
        t.regkinjet2_pt = regkinjet2.Pt();
        t.regkinjet2_e = regkinjet2.E();
        t.regkinjet2_phi = regkinjet2.Phi();
        t.regkinjet2_eta = regkinjet2.Eta();
        t.regkinjet2_mass = regkinjet2.M();
        t.regkinjet2_csvBtag = J.jetCSV[ij2Reg];
        t.regkinjet2_btagSF_M = J.jetbtagSF_M[ij2Reg];
        t.regkinjet2_flavour = J.jetflavour[ij2Reg];
        t.regkinjet2_btagSFErrorUp_M = J.jetbtagSFErrorUp_M[ij2Reg];
        t.regkinjet2_btagSFErrorDown_M = J.jetbtagSFErrorDown_M[ij2Reg];
        t.regkinjet2_btagEff_M = J.jetbtagEff_M[ij2Reg];
        t.regkinjet2_btagEffError_M = J.jetbtagEffError_M[ij2Reg];
        t.regkinjet2_cutbased_wp_level = J.jetcutbased_wp_level[ij2Reg];
        t.regkinjet2_simple_wp_level = J.jetsimple_wp_level[ij2Reg];
        t.regkinjet2_full_wp_level = J.jetfull_wp_level[ij2Reg];
        t.regkinjet2_betaStarClassic = J.jetbetaStarClassic[ij2Reg];
        t.regkinjet2_dR2Mean = J.jetdR2Mean[ij2Reg];
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
        t.kinjet1_cutbased_wp_level = J.jetcutbased_wp_level[ij1];
        t.kinjet1_simple_wp_level = J.jetsimple_wp_level[ij1];
        t.kinjet1_full_wp_level = J.jetfull_wp_level[ij1];
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
        t.kinjet2_cutbased_wp_level = J.jetcutbased_wp_level[ij2];
        t.kinjet2_simple_wp_level = J.jetsimple_wp_level[ij2];
        t.kinjet2_full_wp_level = J.jetfull_wp_level[ij2];
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
        t.gg_DR = pho1.DeltaR(pho2);
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
        t.costhetastar_CS = getCosThetaStar_CS(gg, jj);
        t.regcosthetastar_CS = getCosThetaStar_CS(gg, regjj);
        t.regkincosthetastar_CS = getCosThetaStar_CS(gg, regkinjj);
        t.kincosthetastar_CS = getCosThetaStar_CS(gg, kinjj);
        t.dEta_gg_jj = gg.Eta() - jj.Eta();
        t.dEta_gg_regjj = gg.Eta() - regjj.Eta();
        t.dEta_gg_regkinjj = gg.Eta() - regkinjj.Eta();
        t.dEta_gg_kinjj = gg.Eta() - kinjj.Eta();
        t.dPhi_gg_jj = gg.DeltaPhi( jj );
        t.dPhi_gg_regjj = gg.DeltaPhi( regjj );
        t.dPhi_gg_regkinjj = gg.DeltaPhi( regkinjj );
        t.dPhi_gg_kinjj = gg.DeltaPhi( kinjj );
        t.dR_gg_jj = gg.DeltaR( jj );
        t.dR_gg_regjj = gg.DeltaR( regjj );
        t.dR_gg_regkinjj = gg.DeltaR( regkinjj );
        t.dR_gg_kinjj = gg.DeltaR( kinjj );
        // fill gen info that is not simply fast-forwarded
        if( t.gr_hbb_p4_pt > .01 && t.gr_hgg_p4_pt > .01)
        {
            TLorentzVector gr_hbb, gr_hgg;
            gr_hgg.SetPtEtaPhiM(t.gr_hgg_p4_pt, t.gr_hgg_p4_eta, t.gr_hgg_p4_phi, t.gr_hgg_p4_mass);
            gr_hbb.SetPtEtaPhiM(t.gr_hbb_p4_pt, t.gr_hbb_p4_eta, t.gr_hbb_p4_phi, t.gr_hbb_p4_mass);
            t.gr_hbbhgg_costhetastar_CS = getCosThetaStar_CS(gr_hgg, gr_hbb); 
            t.gr_dEta_gg_bb = gr_hgg.Eta() - gr_hbb.Eta();
            t.gr_dPhi_gg_bb = gr_hgg.DeltaPhi( gr_hbb );
            t.gr_dR_gg_bb = gr_hgg.DeltaR( gr_hbb );
        }
        if( t.gr_hjj_p4_pt > .01 && t.gr_hgg_p4_pt > .01)
        {
            TLorentzVector gr_hjj, gr_hgg;
            gr_hgg.SetPtEtaPhiM(t.gr_hgg_p4_pt, t.gr_hgg_p4_eta, t.gr_hgg_p4_phi, t.gr_hgg_p4_mass);
            gr_hjj.SetPtEtaPhiM(t.gr_hjj_p4_pt, t.gr_hjj_p4_eta, t.gr_hjj_p4_phi, t.gr_hjj_p4_mass);
            t.gr_hjjhgg_costhetastar_CS = getCosThetaStar_CS(gr_hgg, gr_hjj); 
            t.gr_dEta_gg_jj = gr_hgg.Eta() - gr_hjj.Eta();
            t.gr_dPhi_gg_jj = gr_hgg.DeltaPhi( gr_hjj );
            t.gr_dR_gg_jj = gr_hgg.DeltaR( gr_hjj );
        }

        t.weight = t.ev_weight;
        t.evweight = t.ev_evweight;
        t.pu_weight = t.ev_pu_weight;
// if the sample is a prompt-fake sample, remove the events with weights > 50, then apply kfactor * spreadfactor (from Francois' studies)
// equivalent to a global increase of the weights of ~14%
        string pf = "_pf";
        if( inputtree.find(pf) != std::string::npos )
        { 
            float kfactor = 1.0;
            float spreadfactor = 1.0;
//            if( t.evweight > 50 ) continue;
            if( t.njets_kRadionID_and_CSVM == 0 )
            {
                kfactor = 1.17;
//                spreadfactor = 3.06;
            } else {
                kfactor = 1.58;
//                spreadfactor = 5.05;
            }
            t.evweight *= kfactor * spreadfactor; 
        }

        // Apply jackknife weight
        t.evweight *= weight_jackknife;
        // adding this for correct yields out of the control plots:
        t.evweight_w_btagSF = t.evweight;
        if( type == -260 ) t.evweight_w_btagSF *= 1.2822;
        if( type  !=   0 ) t.evweight_w_btagSF *= eventWeight_2jets("medium", J.jetbtagSF_M[ij1], J.jetbtagSF_M[ij2], J.jetbtagEff_M[ij1], J.jetbtagEff_M[ij2], J.jetCSV[ij1], J.jetCSV[ij2]);
        if(lambdaReweight != -1 && type == -2) t.evweight_w_btagSF *= w_hbb_pt_costhetastar_CS->GetBinContent(
                                                                            w_hbb_pt_costhetastar_CS->FindBin(fabs(t.gr_hbbhgg_costhetastar_CS), t.gr_hbb_p4_pt) );
        t.evweight_w_btagSF_reg = t.evweight;
        if( type == -260 ) t.evweight_w_btagSF_reg *= 1.2822;
        if( type !=    0 ) t.evweight_w_btagSF_reg *= eventWeight_2jets("medium", J.jetbtagSF_M[ij1Reg], J.jetbtagSF_M[ij2Reg], J.jetbtagEff_M[ij1Reg], J.jetbtagEff_M[ij2Reg], J.jetCSV[ij1Reg], J.jetCSV[ij2Reg]);
        if(lambdaReweight != -1 && type == -2) t.evweight_w_btagSF_reg *= w_hbb_pt_costhetastar_CS->GetBinContent(
                                                                            w_hbb_pt_costhetastar_CS->FindBin(fabs(t.gr_hbbhgg_costhetastar_CS), t.gr_hbb_p4_pt) );

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
        if(t.njets_kRadionID_and_CSVM == 0)
        {
            t.category = 0;
            nevents[ilevel]++;
            nevents_w[ilevel] += t.evweight;
            nevents_sync[8]++;
            eventcut[ilevel] = "0btag category"; ilevel++; ilevel++;
        } else if(t.njets_kRadionID_and_CSVM == 1)
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
        if(t.njets_kRadionID_and_CSVM == 0)
        {
            t.category = 0;
            nevents[ilevel]++;
            nevents_w[ilevel] += t.evweight;
            nevents_sync[10]++;
            eventcut[ilevel] = Form("0btag category, after mjj cut (%.1f/%.1f and %.1f/%.1f)",min_mjj_1btag, max_mjj_1btag, min_mjj_2btag, max_mjj_2btag); ilevel++; ilevel++;
        }    else if(t.njets_kRadionID_and_CSVM == 1)
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
    
        if(t.njets_kRadionID_and_CSVM == 0)
        {
            t.category = 0;
            nevents[ilevel]++;
            nevents_w[ilevel] += t.evweight;
            eventcut[ilevel] = "0btag category, after kin fit"; ilevel++; ilevel++;
        } else if(t.njets_kRadionID_and_CSVM == 1)
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
        if(t.njets_kRadionID_and_CSVM == 0)
        {
            t.category = 0;
            nevents[ilevel]++;
            nevents_w[ilevel] += t.evweight;
            nevents_sync[12]++;
            eventcut[ilevel] = Form("0btag category, after mggjj cut (%.1f/%.1f and %.1f/%.1f)", min_mggjj_1btag, max_mggjj_1btag, min_mggjj_2btag, max_mggjj_2btag); ilevel++; ilevel++;
        } else if(t.njets_kRadionID_and_CSVM == 1)
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

        flow[iflow] = "After all the cuts"; cutFlow[iflow]++; iflow++;

    } // end of loop over events

    if(DEBUG) cout << "end of loop" << endl;
    if(DEBUG) cout << "ilevelmax= " << ilevelmax << endl;

    for(int i=0 ; i < ilevelmax ; i++)
    cout << "#nevents[" << i << "]= " << nevents[i] << "\t#nevents_w[" << i << "]= " << nevents_w[i] << /*"\t#nevents_w_btagSF_M[" << i << "]= " << nevents_w_btagSF_M[i] << */"\teventcut[" << i << "]= " << eventcut[i] /*<< "\t\t( 1btag= " << nevents1btag[i] << " , 2btag= " << nevents2btag[i] << " ) "*/ << endl;

    if(SYNC)
        for(int i=0 ; i < 14 ; i++)
            cout << nevents_sync[i] << endl;

    if(printCutFlow)
        for(int i= 0 ; i < 9 ; i++)
            cutFlowFile << i << "\t" << cutFlow[i] << "\t" << flow[i] << endl;

    if(SYNC) synchrofile.close();
    if(FULL_DUMP) full_dump.close();
    if(printCutFlow) cutFlowFile.close();

  outfile->cd();
  outtree->Write();
  outfile->Close();
  infile->Close();

    if(applyPhotonIDControlSample) cout << "WARNING: you applied the photon ID control sample, please make sure to reweight in (pt, eta) accordingly" << endl;

    return 0;
}

