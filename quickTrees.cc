// C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
// RooFit headers
// local analysis headers
#include "BTagUtils.h"
#include "quickTrees.h"
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
    string whichJet;
    string fitStrategy;
    int cutLevel;
    int mass;
    int removeUndefinedBtagSF;
    int type;
    int massCutVersion;
    int applyPhotonIDControlSample;
    int applyMtotCut;
    int applyMjjCut;
    int applyMggCut;
    string controlSampleWeights;

    // print out passed arguments
    copy(argv, argv + argc, ostream_iterator<char*>(cout, " ")); cout << endl;
    // argument parsing itself
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("inputfile,i", po::value<string>(&inputfile)->default_value("selected.root"), "input file")
            ("inputtree,t", po::value<string>(&inputtree)->default_value("Radion_m300_8TeV_nm"), "input tree")
            ("outputtree", po::value<string>(&outputtree)->default_value("TCVARS"), "output tree")
            ("outputfile,o", po::value<string>(&outputfile)->default_value("minimum.root"), "output file")
            ("type", po::value<int>(&type)->default_value(0), "same conventions as in h2gglobe: <0 = signal ; =0 = data ; >0 = background")
            ("whichJet", po::value<string>(&whichJet)->default_value(""), "which jet to use, base, kin, regkin, reg")
            ("fitStrategy", po::value<string>(&fitStrategy)->default_value("mgg"), "fit strategy to use, mgg or mggjj")
            ("cutLevel", po::value<int>(&cutLevel)->default_value(0), "0= baseline cuts ; -1= (non-res only) only mjj cuts")
            ("mass", po::value<int>(&mass)->default_value(300), "mass hypothesis (for mass cut switches) for resonant. Set to 0 to get the non-resonant cuts")
            ("removeUndefinedBtagSF", po::value<int>(&removeUndefinedBtagSF)->default_value(0), "remove undefined btagSF_M values (should be used only for the limit trees)")
            ("massCutVersion", po::value<int>(&massCutVersion)->default_value(4), "[0-2]= old stuff 3= preapproval values 4= new baseline optimization (Sep. 2014, to be used with CiC photon ID)")
            ("applyPhotonIDControlSample", po::value<int>(&applyPhotonIDControlSample)->default_value(0), "Invert photon ID CiC cut to populate selection in gjjj instead of ggjj")
            ("controlSampleWeights", po::value<string>(&controlSampleWeights)->default_value("scales_2D_pt_data_4GeVbinning.root"), "file containing the weights for creating the control sample")
            ("applyMtotCut", po::value<int>(&applyMtotCut)->default_value(1), "Apply the mtot cut. Should NOT been played with")
            ("applyMjjCut", po::value<int>(&applyMjjCut)->default_value(1), "Apply the mjj cut. Should NOT been played with")
            ("applyMggCut", po::value<int>(&applyMggCut)->default_value(1), "Apply the mgg cut. Should NOT been played with")
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

    if(strcmp(whichJet.c_str(), "base") == 0) whichJet="";

    cout << "inputfile= " << inputfile << endl;
    cout << "inputtree= " << inputtree << endl;
    cout << "outputfile= " << outputfile << endl;
    cout << "outputtree= " << outputtree << endl;
    cout << "whichJet= " << whichJet << endl;
    cout << "cutLevel= " << cutLevel << endl;
    cout << "fitStrategy= " << fitStrategy << endl;
    cout << "mass= " << mass << endl;

    if( ((strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0)) && (massCutVersion == 3) )
    {
        cout << "ERROR: Regression is not supposed to be used for preapproval-like mass cuts" << endl;
        return 1;
    }
    if( applyPhotonIDControlSample != 0 && type != 0 )
    {
        cout << "ERROR: you are trying to apply the control sample weights to something that is not data" << endl;
        return 1;
    }

    TFile *infile = TFile::Open(inputfile.c_str());
    TTree *intree = (TTree*)infile->Get(inputtree.c_str());
    TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
    TTree *outtree = new TTree(outputtree.c_str(), Form("%s minimal", outputtree.c_str()));

    TFile *csWeightFile = TFile::Open(controlSampleWeights.c_str());
//    TH2F *h2D_pt_data = new TH2F("h2D_pt_data", "h2D_pt_data", 35,20.,160.,35,20.,160.);
    TH2F *h2D_pt_data = (TH2F*)csWeightFile->Get("h2D_pt_data");

    tree_variables t;
    setup_intree(intree, &t, whichJet);
    setup_outtree(outtree, &t);    


    int n_0btag = 0;
    int n_1btag = 0;
    int n_2btag = 0;
    int n_1btag_lowMtot = 0;
    int n_2btag_lowMtot = 0;
    float n_w_0btag = 0.;
    float n_w_1btag = 0.;
    float n_w_2btag = 0.;
    float n_w_1btag_lowMtot = 0.;
    float n_w_2btag_lowMtot = 0.;


    bool removeTrainingEvents = 0; //if using regression jets, remove training events from training samples. These are even events in MSSM and ggHH samples.
    float trainingWeight = 1.0; //apply an extra event weight for those samples
    removeTrainingEvents = ((strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0)) && (strncmp(inputtree.c_str(),"MSSM",4)==0 || strncmp(inputtree.c_str(),"ggHH",4)==0);
    if(removeTrainingEvents) trainingWeight = (float)intree->GetEntries()/(float)intree->GetEntries("event%2==1");


    for(int ievt= 0 ; ievt < (int)intree->GetEntries() ; ievt++)
    {
        intree->GetEntry(ievt);

        if(removeTrainingEvents && t.event%2 == 0 ) continue;

        if(removeUndefinedBtagSF)
        {
            if( t.jet1_btagSF_M == -1001 || t.jet2_btagSF_M == -1001) 
            {
                cout << "WARNING: undefined btagSF_M, skipping the t.event:\tevent= " << t.event << "\tjet1_btagSF_M= " << t.jet1_btagSF_M << "\tjet2_btagSF_M= " << t.jet2_btagSF_M << "\tjet1_pt= " << t.jet1_pt << "\tjet2_pt= " << t.jet2_pt << endl;
                continue;
            }
        }

    t.evWeight *= trainingWeight;
    t.evWeight_w_btagSF = t.evWeight;
    t.weight = t.evWeight;
    t.weightBtagSF = -1000;
    t.weightBtagSFerrUp = -1000;
    t.weightBtagSFerrDown = -1000;

    if( type != 0 )
    {
        t.weightBtagSF = eventWeight_2jets("medium", t.jet1_btagSF_M, t.jet2_btagSF_M, t.jet1_btagEff_M, t.jet2_btagEff_M, t.jet1_csvBtag, t.jet2_csvBtag);
        t.weightBtagSFerrUp = eventWeight_error_2jets("medium", t.jet1_btagSF_M, t.jet1_btagSFErrorUp_M, t.jet2_btagSF_M, t.jet2_btagSFErrorUp_M, t.jet1_btagEff_M, t.jet1_btagEffError_M, t.jet2_btagEff_M, t.jet2_btagEffError_M, t.jet1_flavour, t.jet2_flavour, t.jet1_csvBtag, t.jet2_csvBtag);
        t.weightBtagSFerrDown = eventWeight_error_2jets("medium", t.jet1_btagSF_M, t.jet1_btagSFErrorDown_M, t.jet2_btagSF_M, t.jet2_btagSFErrorDown_M, t.jet1_btagEff_M, t.jet1_btagEffError_M, t.jet2_btagEff_M, t.jet2_btagEffError_M, t.jet1_flavour, t.jet2_flavour, t.jet1_csvBtag, t.jet2_csvBtag);
        t.evWeight_w_btagSF *= t.weightBtagSF;
    }

    if( type == -260 ) t.evWeight_w_btagSF *= 1.2822; // m260 is generated with pythia, while the rest is generated with madgraph. This factor is here to compensate the efficiency difference between the two
    if( type <= -250 ) t.evWeight_w_btagSF /= 1000.; // For increased numerical precision for limit settings, not related to actual physics

    if( type == 0 && applyPhotonIDControlSample != 0) t.evWeight_w_btagSF *= h2D_pt_data->GetBinContent(h2D_pt_data->FindBin(t.pho2_pt, t.pho1_pt));


        if( (strcmp("", whichJet.c_str()) == 0) || (strcmp("reg", whichJet.c_str()) == 0) )
            { t.mjj_wokinfit = t.mjj; t.mtot_wokinfit = t.mtot; }
        
// FITTING THE MGG SPECTRUM
        if( strcmp("mgg", fitStrategy.c_str()) == 0 )
        {
            if( massCutVersion == 3 )
            {// FITTING THE MGG SPECTRUM: newer version: preapproval values (18 dec. 2013)
                if( applyMjjCut && (t.mjj_wokinfit < 85. || t.mjj_wokinfit > 155.) ) continue;
                if( !((type == -2 || type == 0) && !applyMtotCut) ){
                    if( strcmp("kin", whichJet.c_str()) != 0 )
                    {
                        if( mass == 260 && (t.mtot_wokinfit < 225. || t.mtot_wokinfit > 280.) ) continue;
                        if( mass == 270 && (t.mtot_wokinfit < 225. || t.mtot_wokinfit > 295.) ) continue;
                        if( mass == 300 && (t.mtot_wokinfit < 255. || t.mtot_wokinfit > 330.) ) continue;
                        if( mass == 350 && (t.mtot_wokinfit < 310. || t.mtot_wokinfit > 395.) ) continue;
                        if( mass == 400 && (t.mtot_wokinfit < 370. || t.mtot_wokinfit > 440.) ) continue;
                        if( mass == 450 && (t.mtot_wokinfit < 410. || t.mtot_wokinfit > 495.) ) continue;
                        if( mass == 500 && (t.mtot_wokinfit < 445. || t.mtot_wokinfit > 535.) ) continue;
                    } else {
                        // First optimization round done by Francois for using kinfit at low mass (see https://github.com/ResonantHbbHgg/Selection/issues/60) 
                        if( mass == 260 && (t.mtot < 225. || t.mtot > 280.) ) continue; // not updated
                        if( mass == 270 && (t.mtot < 260. || t.mtot > 280.) ) continue;
                        if( mass == 300 && (t.mtot < 290. || t.mtot > 310.) ) continue;
                        if( mass == 350 && (t.mtot < 330. || t.mtot > 375.) ) continue;
                        if( mass == 400 && (t.mtot < 380. || t.mtot > 435.) ) continue;
                        if( mass == 450 && (t.mtot < 410. || t.mtot > 495.) ) continue; // not updated
                        if( mass == 500 && (t.mtot < 445. || t.mtot > 535.) ) continue; // not updated
                    }
                } // if the sample is ggHH and applyMtotCut is switch off, then do not apply any mtot cut
            }
            else
            if( massCutVersion == 4 )
            {// FITTING THE MGG SPECTRUM: new baseline optimization (September 2014)
                if( strcmp("", whichJet.c_str()) == 0 || strcmp("kin", whichJet.c_str()) == 0 )
                {
                    // Non-resonant case
                    // From Badder's optimization (Oct. 6) https://indico.cern.ch/event/343170/session/7/contribution/21/material/slides/0.pdf
                    // Nov 6th update: skype conversation with Badder
                    //      2btag = mjj > 90 && mjj < 155 && mtot > 350 && costhetastar_CS > -0.9 && costhetastar_CS < 0.9
                    //      1btag = mjj > 90 && mjj < 155 && mtot > 360 && costhetastar_CS > -0.65 && costhetastar_CS < 0.65
                    //      plots:
                    //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mjj_study_reweight_m0_cat0_ggHH_kin.pdf
                    //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mjj_study_reweight_m0_cat1_ggHH_kin.pdf
                    //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mggjj_kin_study_reweight_m0_cat0_ggHH_kin.pdf
                    //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mggjj_kin_study_reweight_m0_cat1_ggHH_kin.pdf
                    //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_CosTS_CS_study_reweight_m0_cat0_ggHH_kin.pdf
                    //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_CosTS_CS_study_reweight_m0_cat1_ggHH_kin.pdf
                    if( mass == 0 )
                    { 
                        if( t.njets_kRadionID_and_CSVM>=2 )
                        {
                            if( applyMjjCut && (t.mjj_wokinfit < 90. || t.mjj_wokinfit > 155.) ) continue;
                            if( (cutLevel >= 0) && (fabs(t.costhetastar_CS) > .9) ) continue;
                        }
                        else if( t.njets_kRadionID_and_CSVM==1 )
                        {
                            if( applyMjjCut && (t.mjj_wokinfit < 90. || t.mjj_wokinfit > 155.) ) continue;
                            if( (cutLevel >= 0) && (fabs(t.costhetastar_CS) >  .65) ) continue;
                        }
                    }
                    // Resonant case @ low mass
                    // From Phil's optimization (Sept. 4) https://indico.cern.ch/event/333573/session/14/contribution/18/material/slides/0.pdf
                    else if( mass > 0)
                    {
                        if( t.njets_kRadionID_and_CSVM>=2 )
                        {
                            if( applyMjjCut && (t.mjj_wokinfit < 100. || t.mjj_wokinfit > 140.) ) continue;
                        }
                        else if( t.njets_kRadionID_and_CSVM==1 )
                        {
                            if( applyMjjCut && (t.mjj_wokinfit <  90. || t.mjj_wokinfit > 145.) ) continue;
                        }
                        if( applyMtotCut )
                        { // From Francois' optimization (July 24: https://indico.cern.ch/event/327578/session/9/contribution/26/material/slides/0.pdf)
                            if( mass == 260 && (t.mtot < 250. || t.mtot > 270.) ) continue;
                            if( mass == 270 && (t.mtot < 260. || t.mtot > 280.) ) continue;
                            if( mass == 300 && (t.mtot < 290. || t.mtot > 310.) ) continue;
                            if( mass == 350 && (t.mtot < 330. || t.mtot > 375.) ) continue;
                            if( mass == 400 && (t.mtot < 380. || t.mtot > 435.) ) continue;
                            if( mass == 450 && (t.mtot < 430. || t.mtot > 485.) ) continue;
                            if( mass == 500 && (t.mtot < 480. || t.mtot > 535.) ) continue;
                        }
                    }
                }
                if( strcmp("reg", whichJet.c_str()) == 0 || strcmp("regkin", whichJet.c_str()) == 0 )
                {
                    // Regression case: everything comes from Phil
                    // From Phil's optimization (Sept. 4) https://indico.cern.ch/event/333573/session/14/contribution/18/material/slides/0.pdf
                    // Resonant case @ low mass
                    if( mass > 0)
                    {
                        if( t.njets_kRadionID_and_CSVM>=2 )
                        {
                            if( applyMjjCut && (t.mjj_wokinfit < 100. || t.mjj_wokinfit > 140.) ) continue;
                        }
                        else if( t.njets_kRadionID_and_CSVM==1 )
                        {
                            if( applyMjjCut && (t.mjj_wokinfit < 100. || t.mjj_wokinfit > 150.) ) continue;
                        }
                        if( applyMtotCut )
                        { // From Francois' optimization (July 24: https://indico.cern.ch/event/327578/session/9/contribution/26/material/slides/0.pdf)
                            if( mass == 260 && (t.mtot < 250. || t.mtot > 270.) ) continue;
                            if( mass == 270 && (t.mtot < 260. || t.mtot > 280.) ) continue;
                            if( mass == 300 && (t.mtot < 290. || t.mtot > 310.) ) continue;
                            if( mass == 350 && (t.mtot < 330. || t.mtot > 375.) ) continue;
                            if( mass == 400 && (t.mtot < 380. || t.mtot > 435.) ) continue;
                            if( mass == 450 && (t.mtot < 430. || t.mtot > 485.) ) continue;
                            if( mass == 500 && (t.mtot < 480. || t.mtot > 535.) ) continue;
                        }
                    }
                }
            }
        } // END OF FIT MGG SPECTRUM

// FITTING THE MGGJJ SPECTRUM
        if( strcmp("mggjj", fitStrategy.c_str()) == 0 )
        {
            if( massCutVersion == 3)
            { // Cuts for preapproval (18 dec. 2013)
                if( applyMggCut && (t.mgg < 120. || t.mgg > 130.) ) continue;
                if( applyMjjCut && (t.mjj_wokinfit < 90. || t.mjj_wokinfit > 165.) ) continue;
            }
            if( massCutVersion == 4)
            {
             // From Fabricio's optimization (Aug. 15) https://indico.cern.ch/event/335221/contribution/1/material/slides/0.pdf
             // Probable swap of 1btag and 2btag ?
             // Unclear features in the full table (slide 16): would need to redo the study
             // From the same presentation, it notes a window of (122,128) in mgg. We might need to revert to (120,130) if limits degrade.
             // Nov. 6th: given the expected low yield in cat0 decide to go with PAS cuts for cat0 + Fabricio optim for cat 1
                if( t.njets_kRadionID_and_CSVM>=2 )
                {
                    // Fabricio's cuts
                    // if( applyMggCut && (t.mgg < 122. || t.mgg > 128.) ) continue;
                    // if( applyMjjCut && (t.mjj_wokinfit < 95. || t.mjj_wokinfit > 150.) ) continue;
                    // PAS cuts
                    if( applyMggCut && (t.mgg < 120. || t.mgg > 130.) ) continue;
                    if( applyMjjCut && (t.mjj_wokinfit < 90. || t.mjj_wokinfit > 165.) ) continue;
                }
                else if ( t.njets_kRadionID_and_CSVM==1 )
                {
                    // Fabricio's cuts
                    if( applyMggCut && (t.mgg < 122. || t.mgg > 128.) ) continue;
                    if( applyMjjCut && (t.mjj_wokinfit < 85. || t.mjj_wokinfit > 170.) ) continue;
                    // PAS cuts
                    // if( applyMggCut && (t.mgg < 120. || t.mgg > 130.) ) continue;
                    // if( applyMjjCut && (t.mjj_wokinfit < 90. || t.mjj_wokinfit > 165.) ) continue;
                }
            }
        }

// 2D-FITTING OF BOTH MGG AND MJJ
        if( strcmp("2D", fitStrategy.c_str()) == 0 )
      {
            if( massCutVersion == 4)
            { // Cuts from Phil studies, these cuts are not documented yet but see the following for some comparisons:
              // (Sept. 4) https://indico.cern.ch/event/333573/session/14/contribution/18/material/slides/0.pdf)
                if( applyMggCut && (t.mgg < 100. || t.mgg > 180.) ) continue;
                if( applyMjjCut && (t.mjj_wokinfit < 60. || t.mjj_wokinfit > 180.) ) continue;
        //These cuts were taken from Badder's 1D optimization. Should be revisited for 2D extraction?
        // Nov 6th update: skype conversation with Badder
        //      2btag = mjj > 90 && mjj < 155 && mtot > 350 && costhetastar_CS > -0.9 && costhetastar_CS < 0.9
        //      1btag = mjj > 90 && mjj < 155 && mtot > 360 && costhetastar_CS > -0.65 && costhetastar_CS < 0.65
        //      plots:
        //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mjj_study_reweight_m0_cat0_ggHH_kin.pdf
        //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mjj_study_reweight_m0_cat1_ggHH_kin.pdf
        //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mggjj_kin_study_reweight_m0_cat0_ggHH_kin.pdf
        //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_Mggjj_kin_study_reweight_m0_cat1_ggHH_kin.pdf
        //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_CosTS_CS_study_reweight_m0_cat0_ggHH_kin.pdf
        //          https://bmarzocc.web.cern.ch/bmarzocc/ggHH_optimization/cut_CosTS_CS_study_reweight_m0_cat1_ggHH_kin.pdf
        if( mass == 0)
          {
            if( t.njets_kRadionID_and_CSVM>=2 )
              {
            if( (cutLevel >= 0) && (fabs(t.costhetastar_CS) > .9) ) continue;
              }
            else if( t.njets_kRadionID_and_CSVM==1 )
              {
            if( (cutLevel >= 0) && (fabs(t.costhetastar_CS) > .65) ) continue;
              }
          }
        //apply a different mtot cut for resonant search, same mtot cut as 1D analysis
        else if( (mass > 0) && applyMtotCut )
          { // From Francois' optimization (July 24: https://indico.cern.ch/event/327578/session/9/contribution/26/material/slides/0.pdf)
            if( mass == 260 && (t.mtot < 250. || t.mtot > 270.) ) continue;
            if( mass == 270 && (t.mtot < 260. || t.mtot > 280.) ) continue;
            if( mass == 300 && (t.mtot < 290. || t.mtot > 310.) ) continue;
            if( mass == 350 && (t.mtot < 330. || t.mtot > 375.) ) continue;
            if( mass == 400 && (t.mtot < 380. || t.mtot > 435.) ) continue;
            if( mass == 450 && (t.mtot < 430. || t.mtot > 485.) ) continue;
            if( mass == 500 && (t.mtot < 480. || t.mtot > 535.) ) continue;
          }
        }
      }
        if( strcmp("FTR14001", fitStrategy.c_str()) == 0 )
        { // From future studies analysis FTR-13-001 upgraded into FTR-14-001. The non-resonant HbbHgg documentation is in AN 2014/218: 
          // http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2014/218
            if( massCutVersion == 4)
            { // 
                if( applyMggCut && (t.mgg < 100. || t.mgg > 150.) ) continue;
                if( applyMjjCut && (t.mjj_wokinfit < 70. || t.mjj_wokinfit > 200.) ) continue;
                // cut on DR(g,j) > 1.5 already apply at plot-tree level
                if( t.jj_DR > 2. ) continue;
                if( t.gg_DR > 2. ) continue;
            }
        }

// MGG-LIKE SELECTION FOR MAXIME TO PLAY WITH SYSTEMATICS
        if( strcmp("mgg_noCutOnMTot", fitStrategy.c_str()) == 0 )
        {
            if( t.njets_kRadionID_and_CSVM >= 2 )
            {
                if( (strcmp("", whichJet.c_str()) == 0) || (strcmp("kin", whichJet.c_str()) == 0) )
                    if( t.mjj_wokinfit < 95. || t.mjj_wokinfit > 175. ) continue;
                if( (strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
                    if( t.mjj_wokinfit < 90. || t.mjj_wokinfit > 150. ) continue;
            }
            if( t.njets_kRadionID_and_CSVM == 1 )
            {
                if( (strcmp("", whichJet.c_str()) == 0) || (strcmp("kin", whichJet.c_str()) == 0) )
                    if( t.mjj_wokinfit < 100. || t.mjj_wokinfit > 160. ) continue;
                if( (strcmp("reg", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
                    if( t.mjj_wokinfit < 95. || t.mjj_wokinfit > 140. ) continue;
            }
        }


// FILL THE CATEGORY VARIABLE
// FTR has photon-based categories
        if( strcmp("FTR14001", fitStrategy.c_str()) == 0 )
        {
            if( t.pho1_isEB && t.pho2_isEB )
            {
                t.cut_based_ct = 0;
            } else {
                t.cut_based_ct = 1;
            }
        } 
// nonres search has 4 categories
	else if (mass == 0)
	{
	  if( (cutLevel >= 0 && applyMtotCut) )
	  {
	    if( t.njets_kRadionID_and_CSVM >= 2 && t.mtot > 350. ) {t.cut_based_ct = 0; n_2btag++; n_w_2btag += t.evWeight_w_btagSF;}
            if( t.njets_kRadionID_and_CSVM == 1 && t.mtot > 360. ) {t.cut_based_ct = 1; n_1btag++; n_w_1btag += t.evWeight_w_btagSF;}
	    if( t.njets_kRadionID_and_CSVM >= 2 && t.mtot < 350. ) {t.cut_based_ct = 2; n_2btag_lowMtot++; n_w_2btag_lowMtot += t.evWeight_w_btagSF;}
            if( t.njets_kRadionID_and_CSVM == 1 && t.mtot < 360. ) {t.cut_based_ct = 3; n_1btag_lowMtot++; n_w_1btag_lowMtot += t.evWeight_w_btagSF;}
	  }
	  else
	  {
            if( t.njets_kRadionID_and_CSVM >= 2 ) {t.cut_based_ct = 0; n_2btag++; n_w_2btag += t.evWeight_w_btagSF;}
            if( t.njets_kRadionID_and_CSVM == 1 ) {t.cut_based_ct = 1; n_1btag++; n_w_1btag += t.evWeight_w_btagSF;}
	  }
        }
// res search has the usual 2 categories
	else
	{
            if( t.njets_kRadionID_and_CSVM >= 2 ) {t.cut_based_ct = 0; n_2btag++; n_w_2btag += t.evWeight_w_btagSF;}
            if( t.njets_kRadionID_and_CSVM == 1 ) {t.cut_based_ct = 1; n_1btag++; n_w_1btag += t.evWeight_w_btagSF;}
        }

        // to be in sync with Chiara: if kin fit applied store mjj in mjj_wkinfit and no kin fit in mjj
        t.mjj_wkinfit = t.mjj;
        if( (strcmp("kin", whichJet.c_str()) == 0) || (strcmp("regkin", whichJet.c_str()) == 0) )
            t.mjj = t.mjj_wokinfit;

        outtree->Fill();
    }

    if( mass == 0 && cutLevel >= 0 && applyMtotCut )
    {
      cout << "n_1btag= " << n_1btag << "\tn_2btag= " << n_2btag << "\tn_1btag_lowMtot= " << n_1btag_lowMtot << "\tn_2btag_lowMtot= " << n_2btag_lowMtot << endl;
      cout << "n_w_1btag= " << n_w_1btag << "\tn_w_2btag= " << n_w_2btag << "\tn_w_1btag_lowMtot= " << n_w_1btag_lowMtot << "\tn_w_2btag_lowMtot= " << n_w_2btag_lowMtot << endl;
    }
    else
    {
      cout << "n_1btag= " << n_1btag << "\tn_2btag= " << n_2btag << endl;
      cout << "n_w_1btag= " << n_w_1btag << "\tn_w_2btag= " << n_w_2btag << endl;
    }

  outfile->cd();
  outtree->Write();
  outfile->Close();
  infile->Close();

    return 0;
}
