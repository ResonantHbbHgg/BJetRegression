#include "fillPlot2012_radion_commonNtp.h"

#include "../KinematicFit/DiJetKinFitter.h"

#include "../HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"

using namespace TMVA;

using namespace std;

fillPlot2012_radion_commonNtp::fillPlot2012_radion_commonNtp( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer_commonNtp( "Radion", dataset ) {  
  
  // chiara
  usePhilReg    = false;
  useOlivierReg = false;  
  if (usePhilReg && useOlivierReg) cout << "problem! Can not use two regressions at the same time" << endl; 

  // chiara
  fitToGG         = true;
  fitToFourBodies = false;
  if ( (fitToGG && fitToFourBodies) || (!fitToGG && !fitToFourBodies) ) cout << "problem! Can not fit two variables at the same time" << endl; 

  // chiara
  useKinFit = false;

   

  bTaggerType_ = bTaggerType;
  
  setSelectionType(selectionType);

  // signal samples or bkg/data (to apply weights, and to skip or not even events when using the regression) 
  isSignal  = false;
  isSignalG = false;
  TString dataset_tstr(dataset);
  isSignal  = ( dataset_tstr.Contains("Radion_m") || dataset_tstr.Contains("MSSM_H") || dataset_tstr.Contains("Graviton_M") );
  isSignalG = ( dataset_tstr.Contains("Graviton_M") );
  cout << "#### isSignal = " << isSignal << endl;
  cout << "#### isSignalGraviton = " << isSignalG << endl;
  if (isSignal && usePhilReg) { 
    cout << "#### this is dataset " << dataset_tstr << ": since I'm using the regression, I skip even events" << endl; 
    cout << endl;
  }
  
  // number of processed events for signal samples
  nprocessed = 20000.;
  if (isSignal) {
    TString datasetS(dataset);
    if ( datasetS.Contains("Radion_m300") )  { cout << "#### signal, m=300" << endl; nprocessed = 19972.;}
    if ( datasetS.Contains("Radion_m500") )  { cout << "#### signal, m=500" << endl; nprocessed = 19970.;}
    if ( datasetS.Contains("Radion_m700") )  { cout << "#### signal, m=700" << endl; nprocessed = 19969.;}
    if ( datasetS.Contains("Radion_m1000") ) { cout << "#### signal, m=1000"<< endl; nprocessed = 19951.;}
    if ( datasetS.Contains("Radion_m1500") ) { cout << "#### signal, m=1500"<< endl; nprocessed = 19959.;}
    // 
    if ( datasetS.Contains("Radion_m270") )  { cout << "#### signal, m=270"  << endl; nprocessed = 19996.;}
    if ( datasetS.Contains("Radion_m350") )  { cout << "#### signal, m=350"  << endl; nprocessed = 18498.;}
    if ( datasetS.Contains("Radion_m400") )  { cout << "#### signal, m=400"  << endl; nprocessed = 19697.;}
    if ( datasetS.Contains("Radion_m450") )  { cout << "#### signal, m=450"  << endl; nprocessed = 19999.;}
    if ( datasetS.Contains("Radion_m550") )  { cout << "#### signal, m=550"  << endl; nprocessed = 19995.;}
    if ( datasetS.Contains("Radion_m600") )  { cout << "#### signal, m=600"  << endl; nprocessed = 18197.;}
    if ( datasetS.Contains("Radion_m650") )  { cout << "#### signal, m=650"  << endl; nprocessed = 20000.;}
    if ( datasetS.Contains("Radion_m900") )  { cout << "#### signal, m=900"  << endl; nprocessed = 19996.;}
    if ( datasetS.Contains("Radion_m1100") ) { cout << "#### signal, m=1100" << endl; nprocessed = 19400.;}
    if ( datasetS.Contains("Radion_m1200") ) { cout << "#### signal, m=1200" << endl; nprocessed = 20000.;}
    if ( datasetS.Contains("Radion_m1300") ) { cout << "#### signal, m=1300" << endl; nprocessed = 19996.;}
    if ( datasetS.Contains("Radion_m1400") ) { cout << "#### signal, m=1400" << endl; nprocessed = 19999.;}
    //
    if ( datasetS.Contains("MSSM_H_m260_8TeV") ) { cout << "#### signal, m=260, MSSM non RD"  << endl; nprocessed = 300000.;}
    if ( datasetS.Contains("MSSM_H_m350_8TeV") ) { cout << "#### signal, m=350, MSSM non RD"  << endl; nprocessed = 294530.;}
    //
    if ( datasetS.Contains("MSSM_H_m260_RD") ) { cout << "#### signal, m=260, MSSM RD"  << endl; nprocessed = 300000.;}
    if ( datasetS.Contains("MSSM_H_m300_RD") ) { cout << "#### signal, m=300, MSSM RD"  << endl; nprocessed = 299142.;}
    if ( datasetS.Contains("MSSM_H_m350_RD") ) { cout << "#### signal, m=350, MSSM RD"  << endl; nprocessed = 299571.;}
    //
    if ( datasetS.Contains("Graviton_M-300") )  { cout << "#### graviton, m=300" << endl; nprocessed = 49941.;}
    if ( datasetS.Contains("Graviton_M-500") )  { cout << "#### graviton, m=500" << endl; nprocessed = 49905.;}
    if ( datasetS.Contains("Graviton_M-700") )  { cout << "#### graviton, m=700" << endl; nprocessed = 49911.;}
    if ( datasetS.Contains("Graviton_M-1000") ) { cout << "#### graviton, m=1000"<< endl; nprocessed = 49921.;}

    cout << "#### signal sample: processed events = " << nprocessed << endl; 
    cout << endl;
  }

  // jet regression - Phil
  if (usePhilReg) {
    cout << "using Phil's regression" << endl;
    readerRegres = new Reader( "!Color:!Silent" );
    readerRegres->AddVariable( "jet_eta", &fRegr_eta);
    readerRegres->AddVariable( "jet_emfrac", &fRegr_cef);
    readerRegres->AddVariable( "jet_hadfrac",    &fRegr_chf);
    readerRegres->AddVariable( "jet_nconstituents", &fRegr_nconst);
    readerRegres->AddVariable( "jet_vtx3dL", &fRegr_vtx3dl);
    readerRegres->AddVariable( "MET", &fRegr_met);
    readerRegres->AddVariable( "jet_dPhiMETJet", &fRegr_dPhiMet);
    readerRegres->BookMVA("BDTG0","data/regrWeightsPhil/TMVARegression_Cat0_BDTG.weights.xml");
    readerRegres->BookMVA("BDTG1","data/regrWeightsPhil/TMVARegression_Cat1_BDTG.weights.xml");
  }

  // jet regression - Olivier
  if (useOlivierReg) {
    cout << "using Olivier's regression" << endl;
    readerRegres = new Reader( "!Color:!Silent" );
    readerRegres->AddVariable( "jet_pt",             &fRegr_pt);
    readerRegres->AddVariable( "jet_eta",            &fRegr_eta);
    readerRegres->AddVariable( "jet_emfrac",         &fRegr_cef);
    readerRegres->AddVariable( "jet_nConstituents",  &fRegr_nconst);
    readerRegres->AddVariable( "jet_hadfrac",        &fRegr_chf);
    readerRegres->AddVariable( "jet_secVtxPt",       &fRegr_vtxPt);
    readerRegres->AddVariable( "jet_secVtx3dL",      &fRegr_vtx3dl);
    readerRegres->AddVariable( "ev_met_corr_pfmet",  &fRegr_met);
    readerRegres->AddVariable( "jet_dPhiMet",        &fRegr_dPhiMet);
    readerRegres->BookMVA("BDT","data/regrWeightsOlivier_massDep/regression_m300_BDT.weights.xml"); 
  }

  if (!usePhilReg && !useOlivierReg) cout << "no regression is used" << endl;

  // which fit
  cout << endl;
  if (fitToGG)         cout << "#### fitting mgg spectrum"   << endl;
  if (fitToFourBodies) cout << "#### fitting mggjj spectrum" << endl;
  cout << endl;
}

fillPlot2012_radion_commonNtp::~fillPlot2012_radion_commonNtp() {

  outFile_->Close();
  
  if (!tree_) return;
  delete tree_->GetCurrentFile();
}

double fillPlot2012_radion_commonNtp::delta_phi(double phi1, double phi2) {
  
  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}

void fillPlot2012_radion_commonNtp::finalize() {

  this->Init();

  std::string fullFlags = selectionType_ + "_";
  fullFlags+=bTaggerType_;
  this->set_flags(fullFlags);
  this->createOutputFile();

  outFile_->cd();

 
  // ------------------------------------------------------                                                                       
  // histograms for kinematic optimization                                                                                        

  // jets
  TH1D*  h1_njets = new TH1D("njets", "", 11, -0.5, 10.5);
  h1_njets->Sumw2();
  TH1D*  h1_nbjets_loose = new TH1D("nbjets_loose", "", 11, -0.5, 10.5);
  h1_nbjets_loose->Sumw2();
  TH1D*  h1_nbjets_medium = new TH1D("nbjets_medium", "", 11, -0.5, 10.5);
  h1_nbjets_medium->Sumw2();
  TH1D*  h1_nbjets_tight = new TH1D("nbjets_tight", "", 11, -0.5, 10.5);
  h1_nbjets_tight->Sumw2();

  // photons and jets kinematics
  TH1D* h1_ptphot0 = new TH1D("ptphot0", "", 100, 0., 200.);
  h1_ptphot0->Sumw2();
  TH1D* h1_ptphot1 = new TH1D("ptphot1", "", 100, 0., 200.);
  h1_ptphot1->Sumw2();
  TH1D* h1_runptphot0 = new TH1D("runptphot0", "", 100, 0., 200.);
  h1_runptphot0->Sumw2();

  TH1D* h1_ptjet0 = new TH1D("ptjet0", "", 60, 0., 400.);
  h1_ptjet0->Sumw2();
  TH1D* h1_ptjet1 = new TH1D("ptjet1", "", 30, 0., 200.);
  h1_ptjet1->Sumw2();
  TH1D* h1_runptjet0 = new TH1D("runptjet0", "", 60, 0., 400.);
  h1_runptjet0->Sumw2();

  TH1D* h1_etajet0 = new TH1D("etajet0", "", 30, -3., 3.);
  h1_etajet0->Sumw2();
  TH1D* h1_etajet1 = new TH1D("etajet1", "", 30, -3., 3.);
  h1_etajet1->Sumw2();

  // diphoton kine
  TH1D*  h1_ptDiphot = new TH1D("ptDiphot", "", 100, 0., 500.);
  h1_ptDiphot->Sumw2();
  TH1D*  h1_etaDiphot = new TH1D("etaDiphot", "", 40, -10., 10.);
  h1_etaDiphot->Sumw2();
  TH1D* h1_deltaEtaDiphot = new TH1D("deltaEtaDiphot", "", 40, -5., 5.);
  h1_deltaEtaDiphot->Sumw2();

  // dijets kine
  TH1D*  h1_ptDijet = new TH1D("ptDijet", "", 100, 0., 500.);
  h1_ptDijet->Sumw2();
  TH1D*  h1_etaDijet = new TH1D("etaDijet", "", 40, -10., 10.);
  h1_etaDijet->Sumw2();
  TH1D*  h1_ptRatio = new TH1D("ptRatio", "", 100, 0., 3.);
  h1_ptRatio->Sumw2();
  TH1D*  h1_ptDifference = new TH1D("ptDifference", "", 100, -200., 200.);
  h1_ptDifference->Sumw2();
  TH1D* h1_HTjet = new TH1D("HTjet", "", 60, 0., 400.);
  h1_HTjet->Sumw2();
  TH1D* h1_deltaPhiJets = new TH1D("deltaPhiJets", "", 50, -5., 5.);
  h1_deltaPhiJets->Sumw2();
  TH1D* h1_deltaEtaJets = new TH1D("deltaEtaJets", "", 40, -5., 5.);
  h1_deltaEtaJets->Sumw2();
  TH1D* h1_deltaFabsEtaJets = new TH1D("deltaFabsEtaJets", "", 50, -5., 5.);
  h1_deltaFabsEtaJets->Sumw2();

  // photons and jets
  TH1D*  h1_deltaR = new TH1D("deltaR", "", 100, 0., 3.);
  h1_deltaR->Sumw2();
  TH1D*  h1_deltaPhi = new TH1D("deltaPhi", "", 100, 0., 3.1416);
  h1_deltaPhi->Sumw2();
  TH1D*  h1_deltaEta = new TH1D("deltaEta", "", 50, -5., 5.);
  h1_deltaEta->Sumw2();
  TH1D *h1_deltaR_gg    = new TH1D("deltaR_gg",    "", 50, 0., 3.);
  h1_deltaR_gg->Sumw2();
  TH1D *h1_minDeltaR_gb = new TH1D("minDeltaR_gb", "", 50, 0., 3.);
  h1_minDeltaR_gb->Sumw2();

  // masses
  TH1D* h1_mgg_preselG = new TH1D("mgg_preselG", "", 80, 100., 180.);
  h1_mgg_preselG->Sumw2();
  TH1D* h1_mgg_preselJ = new TH1D("mgg_preselJ", "", 80, 100., 180.);
  h1_mgg_preselJ->Sumw2();
  TH1D* h1_mjj_preselJ = new TH1D("mjj_preselJ", "", 200, 0., 500.);
  h1_mjj_preselJ->Sumw2();

  TH1D* h1_mjj_0btag = new TH1D("mjj_0btag", "", 200, 0., 500.);
  h1_mjj_0btag->Sumw2();
  TH1D* h1_mjj_1btag = new TH1D("mjj_1btag", "", 200, 0., 500.);
  h1_mjj_1btag->Sumw2();
  TH1D* h1_mjj_2btag = new TH1D("mjj_2btag", "", 200, 0., 500.);
  h1_mjj_2btag->Sumw2();

  TH1D* h1_mggjj = new TH1D("mggjj", "", 100, 100., 600.);
  h1_mggjj->Sumw2();
  TH1D* h1_mggjj_0btag = new TH1D("mggjj_0btag", "", 100, 100., 600.);
  h1_mggjj_0btag->Sumw2();
  TH1D* h1_mggjj_1btag = new TH1D("mggjj_1btag", "", 100, 100., 600.);
  h1_mggjj_1btag->Sumw2();
  TH1D* h1_mggjj_2btag = new TH1D("mggjj_2btag", "", 100, 100., 600.);
  h1_mggjj_2btag->Sumw2();

  // helicity angles
  TH1D* h1_cosThetaStar = new TH1D("cosThetaStar", "", 50, -1.0001, 1.0001);
  h1_cosThetaStar->Sumw2();

  // kin fit study      
  TH1D* h1_kinfit_chiSquareProbH = new TH1D("kinfit_chiSquareProbH", "", 1000, 0., 1.0001);
  h1_kinfit_chiSquareProbH->Sumw2();
  TH1D* h1_mVstar_kinfit = new TH1D("mVstar_kinfit", "", 100, 100., 600.);
  h1_mVstar_kinfit->Sumw2();
  TH1D* h1_mVstar = new TH1D("mVstar", "", 100, 100., 600.);
  h1_mVstar->Sumw2();
  TH1D* h1_mggjj_kinfit_1btag = new TH1D("mggjj_kinfit_1btag", "", 100, 100., 600.);
  h1_mggjj_kinfit_1btag->Sumw2();
  TH1D* h1_mggjj_kinfit_2btag = new TH1D("mggjj_kinfit_2btag", "", 100, 100., 600.);
  h1_mggjj_kinfit_2btag->Sumw2();
  TH1D* h1_mggjj_nokinfit_1btag = new TH1D("mggjj_nokinfit_1btag", "", 100, 100., 600.);
  h1_mggjj_nokinfit_1btag->Sumw2();
  TH1D* h1_mggjj_nokinfit_2btag = new TH1D("mggjj_nokinfit_2btag", "", 100, 100., 600.);
  h1_mggjj_nokinfit_2btag->Sumw2();

  // 2d plots
  TH2D* h2_mggjj_vs_mjj_1btag = new TH2D("mggjj_vs_mjj_1btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_1btag->Sumw2();
  TH2D* h2_mggjj_vs_mjj_2btag = new TH2D("mggjj_vs_mjj_2btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_2btag->Sumw2();

  TH2D* h2_mggjj_vs_mjj_kinfit_1btag = new TH2D("mggjj_vs_mjj_kinfit_1btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_kinfit_1btag->Sumw2();
  TH2D* h2_mggjj_vs_mjj_kinfit_2btag = new TH2D("mggjj_vs_mjj_kinfit_2btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_kinfit_2btag->Sumw2();


  // -----------------------------------------------------------------------------                                                 
  // variables to be saved in an output tree with selected events 
  float massggnewvtx_t;
  float ptphot1_t, runptphot1_t, ptphot2_t;
  float etaphot1_t, etaphot2_t;
  int   cicphot1_t, cicphot2_t;
  float r9phot1_t, r9phot2_t;
  float ptgg_t, etagg_t, absetagg_t;
  int njets_t;
  float ptcorrjet1_t, ptcorrjet2_t, runptcorrjet1_t;
  float etajet1_t, etajet2_t;
  float btagSFjet1_t, btagSFjet2_t;
  int flavourjet1_t, flavourjet2_t;
  float emfjet1_t, emfjet2_t; 
  float chfjet1_t, chfjet2_t; 
  float nconstjet1_t, nconstjet2_t;
  float vtxPtjet1_t, vtxPtjet2_t;
  float vtx3dljet1_t, vtx3dljet2_t;
  float deltaphijj_t, deltaetajj_t;
  float invmassjet_t, ptjj_t, etajj_t, massggjj_t;
  float massggjjNoKF_t, invmassjetWithKF_t;
  float deltaphijjgg_t, deltaetajjgg_t, deltaRjjgg_t;
  float deltaR_gg_t, minDeltaR_gb_t;
  int btagCategory_t, nbjets_loose_t, nbjets_medium_t, nbjets_tight_t;
  int theGammaCategory_t;
  int theCategory_t;
  float chiSquareProbH_t, absCosThetaStar_t;
  int nvtx_t;
  int theVertex_t;
  int isj1btagged_t, isj2btagged_t;
  float mjj_kin_t, mggjj_kin_t;
  float HT_jet_t, HT_gjet_t;
  float randMgj1_t, randMgj2_t;
  float weight_t;
  float weightNoSF_t;

  /*
  TTree* myTrees = new TTree();
  myTrees->SetName("myTrees");
  myTrees->Branch( "run", &run, "run/F" );
  myTrees->Branch( "event", &event, "event/F" );
  myTrees->Branch( "massggnewvtx", &massggnewvtx_t, "massggnewvtx_t/F" );
  myTrees->Branch( "ptPhot1", &ptphot1_t, "ptphot1_t/F" );
  myTrees->Branch( "runptPhot1", &runptphot1_t, "runptPhot1_t/F" );
  myTrees->Branch( "ptPhot2", &ptphot2_t, "ptphot2_t/F" );
  myTrees->Branch( "etaPhot1", &etaphot1_t, "etaphot1_t/F" );
  myTrees->Branch( "etaPhot2", &etaphot2_t, "etaphot2_t/F" );
  myTrees->Branch( "cicPhot1", &cicphot1_t, "cicphot1_t/I" );
  myTrees->Branch( "cicPhot2", &cicphot2_t, "cicphot2_t/I" );
  myTrees->Branch( "r9Phot1", &r9phot1_t, "r9phot1_t/F" );
  myTrees->Branch( "r9Phot2", &r9phot2_t, "r9phot2_t/F" );
  myTrees->Branch( "ptgg", &ptgg_t, "ptgg_t/F" );
  myTrees->Branch( "etagg", &etagg_t, "etagg_t/F" );
  myTrees->Branch( "absetagg", &absetagg_t, "absetagg_t/F" );
  myTrees->Branch( "njets", &njets_t, "njets_t/I" );
  myTrees->Branch( "runPtCorrJet1", &runptcorrjet1_t, "runptcorrJet1_t/F" );
  myTrees->Branch( "ptCorrJet1", &ptcorrjet1_t, "ptcorrJet1_t/F" );
  myTrees->Branch( "ptCorrJet2", &ptcorrjet2_t, "ptcorrJet2_t/F" );
  myTrees->Branch( "etaJet1", &etajet1_t, "etajet1_t/F" );
  myTrees->Branch( "etaJet2", &etajet2_t, "etajet2_t/F" );
  myTrees->Branch( "btagSFJet1", &btagSFjet1_t, "btagSFjet1_t/F" );
  myTrees->Branch( "btagSFJet2", &btagSFjet2_t, "btagSFjet2_t/F" );
  myTrees->Branch( "flavourJet1", &flavourjet1_t, "flavourjet1_t/I" );
  myTrees->Branch( "flavourJet2", &flavourjet2_t, "flavourjet2_t/I" );
  myTrees->Branch( "emfJet1", &emfjet1_t, "emfjet1_t/F" );
  myTrees->Branch( "emfJet2", &emfjet2_t, "emfjet2_t/F" );
  myTrees->Branch( "chfJet1", &chfjet1_t, "chfjet1_t/F" );
  myTrees->Branch( "chfJet2", &chfjet2_t, "chfjet2_t/F" );
  myTrees->Branch( "nconstJet1", &nconstjet1_t, "nconstjet1_t/F");
  myTrees->Branch( "nconstJet2", &nconstjet2_t, "nconstjet2_t/F");
  myTrees->Branch( "vtxPtJet1", &vtxPtjet1_t, "vtxPtjet1_t/F");
  myTrees->Branch( "vtxPtJet2", &vtxPtjet2_t, "vtxPtjet2_t/F");
  myTrees->Branch( "vtx3dlJet1", &vtx3dljet1_t, "vtx3dljet1_t/F");
  myTrees->Branch( "vtx3dlJet2", &vtx3dljet2_t, "vtx3dljet2_t/F");
  myTrees->Branch( "deltaphijj", &deltaphijj_t, "deltaphijj_t/F");
  myTrees->Branch( "deltaetajj", &deltaetajj_t, "deltaetajj_t/F");
  myTrees->Branch( "mjj", &invmassjet_t, "invmassjet_t/F" );
  myTrees->Branch( "ptjj", &ptjj_t, "ptjj_t/F" );
  myTrees->Branch( "etajj", &etajj_t, "etajj_t/F" );
  myTrees->Branch( "mggjj", &massggjj_t, "massggjj_t/F" );
  myTrees->Branch( "deltaphiggjj", &deltaphijjgg_t, "deltaphijjgg_t/F");
  myTrees->Branch( "deltaetaggjj", &deltaetajjgg_t, "deltaetajjgg_t/F");
  myTrees->Branch( "deltaRggjj", &deltaRjjgg_t, "deltaRjjgg_t/F");
  myTrees->Branch( "btagCategory", &btagCategory_t, "btagCategory_t/I" );
  myTrees->Branch( "theCategory",  &theCategory_t,  "theCategory_t/I" );
  myTrees->Branch( "theGammaCategory",  &theGammaCategory_t,  "theGammaCategory_t/I" );
  myTrees->Branch( "nbjets_loose",  &nbjets_loose_t,  "nbjets_loose_t/I" );
  myTrees->Branch( "nbjets_medium", &nbjets_medium_t, "nbjets_medium_t/I" );
  myTrees->Branch( "nbjets_tight",  &nbjets_tight_t,  "nbjets_tight_t/I" );
  myTrees->Branch( "chiSquareProbH", &chiSquareProbH_t, "chiSquareProbH_t/F" );
  myTrees->Branch( "absCosThetaStar", &absCosThetaStar_t, "absCosThetaStar_t/F" );
  myTrees->Branch( "nvtx", &nvtx_t, "nvtx_t/I" );
  myTrees->Branch( "vertex", &theVertex_t, "theVertex_t/I" );
  myTrees->Branch( "isBtagJet1", &isj1btagged_t, "isj1btagged_t/I" );
  myTrees->Branch( "isBtagJet2", &isj2btagged_t, "isj2btagged_t/I" );
  myTrees->Branch( "mjj_kin", &mjj_kin_t, "mjj_kin_t/F" );
  myTrees->Branch( "mggjj_kin", &mggjj_kin_t, "mggjj_kin_t/F" );
  myTrees->Branch( "HT_jet",  &HT_jet_t,  "HT_jet_t/F" );
  myTrees->Branch( "HT_gjet", &HT_gjet_t, "HT_gjet_t/F" );
  myTrees->Branch( "deltaR_gg",    &deltaR_gg_t,    "deltaR_gg/F" );
  myTrees->Branch( "minDeltaR_gb", &minDeltaR_gb_t, "minDeltaR_gb/F" );
  myTrees->Branch( "MgjRandom1", &randMgj1_t, "MgjRandom1/F" );
  myTrees->Branch( "MgjRandom2", &randMgj2_t, "MgjRandom2/F" );
  myTrees->Branch( "weight", &weight_t, "weight_t/F" );
  myTrees->Branch( "weightNoSF", &weightNoSF_t, "weightNoSF_t/F" );
  */

  // tree to compute final limits in Alexandra's stuff                                                          
  TTree* TCVARS = new TTree("TCVARS","two photon two jet selection");                                                          
  TCVARS->SetName("TCVARS");                                                                                                   
  // TCVARS->Branch( "event",         &event,              "event/F" );
  TCVARS->Branch( "mgg",           &massggnewvtx_t,     "massggnewvtx_t/F" );                                                       
  TCVARS->Branch( "mjj",           &invmassjet_t,       "invmassjet_t/F" );                                                         
  TCVARS->Branch( "mjj_wkinfit",   &invmassjetWithKF_t, "invmassjetWithKF_t/F" );                                                         
  TCVARS->Branch( "mtot",          &massggjj_t,         "massggjj_t/F" );                                                           
  TCVARS->Branch( "mtot_wokinfit", &massggjjNoKF_t,     "massggjjNoKF_t/F" );  
  TCVARS->Branch( "cut_based_ct",  &theCategory_t,      "theCategory_t/I" );                                                        
  TCVARS->Branch( "evWeight",      &weight_t,           "weight_t/F" );   
  // TCVARS->Branch( "evWeightNoSF", &weightNoSF_t,     "weightNoSF_t/F" );   

  // ------------------------------------------------------                                                                       
  // for the kinematic fits, assuming the two jets come from a 125 GeV Higgs                                                      
  float Hmass = 125.;
  DiJetKinFitter* fitter_jetsH = new DiJetKinFitter( "fitter_jetsH", "fitter_jets", Hmass );

  // for the helicity angles study                                                                                                
  HelicityLikelihoodDiscriminant *helicityDiscriminator = new HelicityLikelihoodDiscriminant();
  int seed = 110;
  TRandom3* coin = new TRandom3(seed);
  
  // for random studies
  TRandom* myrand = new TRandom(seed);

  // analysis                                                                                                                     
  if (tree_ == 0) return;

  Long64_t nentries = tree_->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = tree_->GetEntry(jentry);   nbytes += nb;


    // event weight for MC events: xsec, PU and smearings. It correspond to 19.706/fb from pixellumicalc
    weight_t = evweight;  
    weightNoSF_t = evweight;     // this is to keep track of the value before SFs for btagging are applied

    // Phil's regression trained on even events from signal samples : skip them
    if (usePhilReg && isSignal && ((int)event%2 == 0) ) continue;

    // ---------------------------------------------------                                                                        
    // gamma-gamma analysis                                                                                                       

    // photons acceptance                                                                                                         
    if((TMath::Abs(ph1_SCEta)>1.4442&&TMath::Abs(ph1_SCEta)<1.566)||(TMath::Abs(ph2_SCEta)>1.4442&&TMath::Abs(ph2_SCEta)<1.566)
       || TMath::Abs(ph1_SCEta)>2.5 || TMath::Abs(ph2_SCEta)>2.5) continue;  // acceptance                                      

    // photon id  
    bool idphot1(0), idphot2(0);
    idphot1 = (ph1_ciclevel >= cicselection);
    idphot2 = (ph2_ciclevel >= cicselection);
    if(cicselection>0) {
      if(!(idphot1)) continue;
      if(!(idphot2)) continue;
    }

    /*
    // chiara: photonID for control sample
    bool idphot1(0), idphot2(0);
    bool looseidphot1(0), looseidphot2(0);
    idphot1 = (ph1_ciclevel >= cicselection);
    idphot2 = (ph2_ciclevel >= cicselection);
    looseidphot1 = (ph1_ciclevel < cicselection && ph1_ciclevel >= 0);
    looseidphot2 = (ph2_ciclevel < cicselection && ph2_ciclevel >= 0);
    if(cicselection>0) {
      if (idphot1 && idphot2) continue;
      if (!looseidphot1 && !looseidphot2) continue;
      if ( !((idphot1 && looseidphot2) || (idphot2 && looseidphot1)) ) continue;
    }
    */

    // extra cuts on photons: splitting events per photon class if needed                                                         
    int myEB=2;
    if( (fabs(ph1_SCEta)>1.4442 || fabs(ph2_SCEta)>1.4442) ) myEB=0;
    if( (fabs(ph1_SCEta)<1.4442 && fabs(ph2_SCEta)<1.4442) ) myEB=1;

    // the two selected photons for the analysis                                                                                  
    TLorentzVector t4phot1, t4phot2;
    t4phot1.SetPtEtaPhiM(ph1_pt,ph1_eta,ph1_phi,0.);
    t4phot2.SetPtEtaPhiM(ph2_pt,ph2_eta,ph2_phi,0.);
    TVector3 t3phot1, t3phot2;
    t3phot1.SetPtEtaPhi(ph1_pt,ph1_eta,ph1_phi);
    t3phot2.SetPtEtaPhi(ph2_pt,ph2_eta,ph2_phi);
    TLorentzVector t4diPhot;
    t4diPhot.SetPtEtaPhiM( dipho_pt, dipho_eta, dipho_phi, PhotonsMass );
    TVector3 t3diPhot;
    t3diPhot.SetPtEtaPhi( dipho_pt, dipho_eta, dipho_phi );

    // further cuts on photons -----------------------                                                                            

    // invariant mass cut on photons                                                                                              
    if (PhotonsMass<100 || PhotonsMass>180) continue;  

    // control plots to check the preselection                                                                                    
    h1_ptphot0->Fill(ph1_pt, weight_t);
    h1_ptphot1->Fill(ph2_pt, weight_t);
    h1_runptphot0->Fill(ph1_pt*120./PhotonsMass, weight_t);

    // photons pt cuts                                                                                                            
    if(ph1_pt<ptphot1cut * PhotonsMass/120.) continue;         // pt first photon                                                
    if(ph2_pt<ptphot2cut * PhotonsMass/120.) continue;         // pt second photon                                               

    // control plots after preselection on photons (pT, acceptance and ID)                                                        
    h1_mgg_preselG->Fill( PhotonsMass, weight_t );



    // ------------------------------ jets ------------------------------------------
    // jets: preparing vectors with the infos used later on
    ecorrjet[0]     = j1_e;           ecorrjet[1]     = j2_e;           ecorrjet[2]     = j3_e;           ecorrjet[3]     = j4_e;
    ptcorrjet[0]    = j1_pt;          ptcorrjet[1]    = j2_pt;          ptcorrjet[2]    = j3_pt;          ptcorrjet[3]    = j4_pt;
    etajet[0]       = j1_eta;         etajet[1]       = j2_eta;         etajet[2]       = j3_eta;         etajet[3]       = j4_eta;
    phijet[0]       = j1_phi;         phijet[1]       = j2_phi;         phijet[2]       = j3_phi;         phijet[3]       = j4_phi;
    btagjprobjet[0] = j1_jetProbBtag; btagjprobjet[1] = j2_jetProbBtag; btagjprobjet[2] = j3_jetProbBtag; btagjprobjet[3] = j4_jetProbBtag;
    btagcsvjet[0]   = j1_csvBtag;     btagcsvjet[1]   = j2_csvBtag;     btagcsvjet[2]   = j3_csvBtag;     btagcsvjet[3]   = j4_csvBtag;
    cefjet[0]       = j1_emfrac;      cefjet[1]       = j2_emfrac;      cefjet[2]       = j3_emfrac;      cefjet[3]       = j4_emfrac; 
    chfjet[0]       = j1_hadfrac;     chfjet[1]       = j2_hadfrac;     chfjet[2]       = j3_hadfrac;     chfjet[3]       = j4_hadfrac; 
    vtxPtjet[0]     = j1_secVtxPt;    vtxPtjet[1]     = j2_secVtxPt;    vtxPtjet[2]     = j3_secVtxPt;    vtxPtjet[3]     = j4_secVtxPt; 
    vtx3dljet[0]    = j1_secVtx3dL;   vtx3dljet[1]    = j2_secVtx3dL;   vtx3dljet[2]    = j3_secVtx3dL;   vtx3dljet[3]    = j4_secVtx3dL; 
    btagSF[0]       = j1_btagSF;      btagSF[1]       = j2_btagSF;      btagSF[2]       = j3_btagSF;      btagSF[3]       = j4_btagSF;              
    flavour[0]      = j1_flavour;     flavour[1]      = j2_flavour;     flavour[2]      = j3_flavour;     flavour[3]      = j4_flavour; 
    btagEff[0]      = j1_btagEff;     btagEff[1]      = j2_btagEff;     btagEff[2]      = j3_btagEff;     btagEff[3]      = j4_btagEff;     
    nconstjet[0]    = (float)(j1_nNeutrals + j1_nCharged);   
    nconstjet[1]    = (float)(j2_nNeutrals + j2_nCharged);
    nconstjet[2]    = (float)(j3_nNeutrals + j3_nCharged);
    nconstjet[3]    = (float)(j4_nNeutrals + j4_nCharged);

    // applying jet regression
    if (usePhilReg || useOlivierReg) {
      TVector3 tempT3jet, t3met;
      t3met.SetPtEtaPhi(met_corr_pfmet, 0, met_corr_phi_pfmet);   

      for (int ii=0; ii<4; ii++) {
	if (ptcorrjet[ii]<-1) continue;
	fRegr_pt      = ptcorrjet[ii]; 
	fRegr_eta     = etajet[ii];
	fRegr_cef     = cefjet[ii];
	fRegr_nconst  = nconstjet[ii];
	fRegr_chf     = chfjet[ii];
	fRegr_vtxPt   = vtxPtjet[ii];
	fRegr_vtx3dl  = vtx3dljet[ii];
	fRegr_met = met_corr_pfmet;   
	tempT3jet.SetPtEtaPhi(ptcorrjet[ii], etajet[ii], phijet[ii]);
	if (useOlivierReg) fRegr_dPhiMet = tempT3jet.DeltaPhi(t3met);
	if (usePhilReg)    fRegr_dPhiMet = fabs(tempT3jet.DeltaPhi(t3met));

	float thePtCorr = 10000;
	if (usePhilReg) {
	  float corr_factor = ptcorrjet[ii] < 90 ? (float)(readerRegres->EvaluateRegression("BDTG0")[0]) : (float)(readerRegres->EvaluateRegression("BDTG1")[0]);                                                         
	  thePtCorr = ptcorrjet[ii]*corr_factor; 
	}
	if (useOlivierReg) thePtCorr = (float)(readerRegres->EvaluateMVA("BDT"));

	// corrected pT and energy
	float correctionFactor = thePtCorr/ptcorrjet[ii];
	if (correctionFactor<0) cout << "chiara, error: regression < 0!" << endl;

	ptcorrjet[ii] = thePtCorr;
	ecorrjet[ii] *= correctionFactor; 
      }
    }

    // jets, no btagging, passing cut based jetID 
    vector<int> v_puIdJets;
    for (int ij=0; ij<4; ij++) { 
      if ( ptcorrjet[ij]<-1)               continue;
      if ( btagcsvjet[ij]<=0)              continue;   
      if ( ptcorrjet[ij] < ptjetacccut )   continue;
      if ( fabs(etajet[ij])>etajetacccut ) continue;
      if ( !passCutBasedJetId(ij) )        continue;
      v_puIdJets.push_back(ij);
    }

    // jets passing btagging ( + eta/pT/PUid cuts)  
    vector<int> v_looseJP,  v_mediumJP,  v_tightJP;
    vector<int> v_looseCSV, v_mediumCSV, v_tightCSV;
    for (int ij=0; ij<int(v_puIdJets.size());ij++) {
      int index = v_puIdJets[ij];
      if (btagjprobjet[index]>0.275) v_looseJP.push_back(index);
      if (btagjprobjet[index]>0.545) v_mediumJP.push_back(index);
      if (btagjprobjet[index]>0.790) v_tightJP.push_back(index);
      if (btagcsvjet[index]>0.244)   v_looseCSV.push_back(index);
      if (btagcsvjet[index]>0.679)   v_mediumCSV.push_back(index);
      if (btagcsvjet[index]>0.898)   v_tightCSV.push_back(index);
    }

    // control plots (before any cuts on jets)                                                                                    
    h1_njets->Fill( v_puIdJets.size(), weight_t );
    if (bTaggerType_=="JP") {
      h1_nbjets_loose ->Fill( v_looseJP.size(),  weight_t );
      h1_nbjets_medium->Fill( v_mediumJP.size(), weight_t );
      h1_nbjets_tight ->Fill( v_tightJP.size(),  weight_t );
    } else if (bTaggerType_=="CSV") {
      h1_nbjets_loose ->Fill( v_looseCSV.size(),  weight_t );
      h1_nbjets_medium->Fill( v_mediumCSV.size(), weight_t );
      h1_nbjets_tight ->Fill( v_tightCSV.size(),  weight_t );
    } else {
      cout << "this btag algo does not exists" << endl;
    }

    // at least 2 preselected jets                                                                                                
    if (v_puIdJets.size()<2)  continue;


    // choice of analysis jets ---------------------                                                                              
    
    // highest pT btagged jet
    int jet1btag   = -1;
    bool isJustOne = false;
    if (bTaggerType_=="CSV") {
      if (v_mediumCSV.size()>=1) jet1btag = myHighestPtJet(v_mediumCSV);
      if (v_mediumCSV.size()==1) isJustOne = true;
    } else if (bTaggerType_=="JP") {
      if (v_mediumJP.size()>=1) jet1btag = myHighestPtJet(v_mediumJP);
      if (v_mediumJP.size()==1) isJustOne = true;
    } else {
      cout << "this btag algo does not exists" << endl;
    }

    // chiara: da commentare solo x salvare anche 0bjets - init
    // at least 1 btagged jet 
    if ( jet1btag<0 ) continue;
    // chiara: da commentare solo x salvare anche 0bjets - end
   

    // choose the two jets with highest pT(jj), giving preference to btagged ones 
    int jet1 = -1;
    int jet2 = -1;
    int isj1btagged = 0;
    int isj2btagged = 0;
    float maxPtBtag = -999.;



    // chiara: da scommentare solo x salvare anche 0bjets - init 
    /*
    if (jet1btag<0) {
      for (int jetA=0; jetA<(v_puIdJets.size()-1); jetA++) {
	for (int jetB=jetA+1; jetB<v_puIdJets.size(); jetB++) {
	  TLorentzVector t4jetA, t4jetB;
	  int indexA = v_puIdJets[jetA];
	  int indexB = v_puIdJets[jetB];
	  t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	  t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	  TLorentzVector t4AB = t4jetA + t4jetB;
	  if ( t4AB.Pt() > maxPtBtag) {
	    maxPtBtag = t4AB.Pt();
	    jet1 = indexA;
	    jet2 = indexB;
	    isj1btagged = false;
	    isj2btagged = false;
	  }
	}
      }
    }
    */
    // chiara: da scommentare solo x salvare anche 0bjets - end     

    // if (jet1btag>=0) {  // chiara: da scommentare solo x salvare anche 0bjets 


    if( isJustOne ) {  // =1 btagged jets
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
        int index = v_puIdJets[jet];
        if (index==jet1btag) continue;
        TLorentzVector t4jet, t4bjet;
        t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
        t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag],etajet[jet1btag],phijet[jet1btag],ecorrjet[jet1btag]);
	
        TLorentzVector t4_bnotb = t4jet + t4bjet;
        if ( t4_bnotb.Pt() > maxPtBtag) {
          maxPtBtag = t4_bnotb.Pt();
          jet1 = jet1btag;
          jet2 = index;
          isj1btagged = true;
          isj2btagged = false;
        }
      }

    } else {  // >1 btagged jets  
      
      if (bTaggerType_=="CSV") {
        for (int jetA=0; jetA<(v_mediumCSV.size()-1); jetA++) {
          for (int jetB=jetA+1; jetB<v_mediumCSV.size(); jetB++) {
            TLorentzVector t4jetA, t4jetB;
            int indexA = v_mediumCSV[jetA];
            int indexB = v_mediumCSV[jetB];
            t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
            t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
            TLorentzVector t4AB = t4jetA + t4jetB;
	    if ( t4AB.Pt() > maxPtBtag) {
              maxPtBtag = t4AB.Pt();
              jet1 = indexA;
              jet2 = indexB;
              isj1btagged = true;
              isj2btagged = true;
            }
          }
        }
      } else if (bTaggerType_=="JP") {
	for (int jetA=0; jetA<(v_mediumJP.size()-1); jetA++) {
          for (int jetB=jetA+1; jetB<v_mediumJP.size(); jetB++) {
            TLorentzVector t4jetA, t4jetB;
            int indexA = v_mediumJP[jetA];
            int indexB = v_mediumJP[jetB];
            t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
            t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
            TLorentzVector t4AB = t4jetA + t4jetB;
            if ( t4AB.Pt() > maxPtBtag) {
              maxPtBtag = t4AB.Pt();
              jet1 = indexA;
              jet2 = indexB;
              isj1btagged = true;
              isj2btagged = true;
            }
          }
        }
      } else {
        cout << "this btag algo does not exists" << endl;
      }
    }
    // }  // chiara: da scommentare solo x salvare anche 0bjets  
    
    if (jet1<0 || jet2<0) { 
      cout << "problem! at least 1 jet not correctly selected" << endl;
      continue;
    }

    // swap jet1 / jet2 if needed to order in pT 
    if (ptcorrjet[jet1]<ptcorrjet[jet2]) {
      int myjet        = jet1;
      bool ismyjetbtag = isj1btagged;
      jet1 = jet2;
      jet2 = myjet;
      isj1btagged = isj2btagged;
      isj2btagged = ismyjetbtag;
    }

    
    // selected analysis jets                                                                                                     
    TLorentzVector t4jet1, t4jet2;
    t4jet1.SetPtEtaPhiE(ptcorrjet[jet1],etajet[jet1],phijet[jet1],ecorrjet[jet1]);
    t4jet2.SetPtEtaPhiE(ptcorrjet[jet2],etajet[jet2],phijet[jet2],ecorrjet[jet2]);
    TVector3 t3jet1, t3jet2;
    t3jet1.SetPtEtaPhi(ptcorrjet[jet1],etajet[jet1],phijet[jet1]);
    t3jet2.SetPtEtaPhi(ptcorrjet[jet2],etajet[jet2],phijet[jet2]);
    TLorentzVector t4diJet = t4jet1 + t4jet2;
    
    // btag scale factors and flavours for the chosen jets
    float btagSFJet1  = btagSF[jet1];
    float btagSFJet2  = btagSF[jet2];
    float btagEffJet1 = btagEff[jet1];  
    float btagEffJet2 = btagEff[jet2];
    int flavourJet1   = flavour[jet1];
    int flavourJet2   = flavour[jet2];
    float btagOutJet1 = btagcsvjet[jet1]; 
    float btagOutJet2 = btagcsvjet[jet2];

    // regression variables
    float emfJet1    = cefjet[jet1];
    float emfJet2    = cefjet[jet2];  
    float chfJet1    = chfjet[jet1];
    float chfJet2    = chfjet[jet2];  
    float nconstJet1 = nconstjet[jet1];
    float nconstJet2 = nconstjet[jet2];
    float vtxPtJet1  = vtxPtjet[jet1]; 
    float vtxPtJet2  = vtxPtjet[jet2]; 
    float vtx3dlJet1 = vtx3dljet[jet1]; 
    float vtx3dlJet2 = vtx3dljet[jet2]; 

    // invariant mass cut on jets                                                                                                 
    float invMassJJ = t4diJet.M();

    // invariant mass plots after the two preselection on jets and photons                                                        
    h1_mgg_preselJ->Fill( PhotonsMass, weight_t );
    h1_mjj_preselJ->Fill( t4diJet.M(), weight_t );

    
    // -------------------------------------------------------------------------                                                  
    // jet related kinematic variables to further study the jet selection                                                         
    float ptJ1       = t4jet1.Pt();
    float ptJ2       = t4jet2.Pt();
    float etaJ1      = t4jet1.Eta();
    float etaJ2      = t4jet2.Eta();
    float deltaEtaJJ = t4jet1.Eta() - t4jet2.Eta();
    float deltaPhiJJ = t4jet1.DeltaPhi(t4jet2);
    float ptJJ       = t4diJet.Pt();
    float etaJJ      = t4diJet.Eta();

    // jet/photon related kinematic variables                                                                                     
    float deltaR_ggjj   = t4diPhot.DeltaR(t4diJet);
    float deltaPhi_ggjj = t4diPhot.DeltaPhi(t4diJet);
    float deltaEta_ggjj = t4diPhot.Eta() - t4diJet.Eta();

    // gg and gamma,jet
    float deltaR_gg    = t3phot1.DeltaR(t3phot2);
    float deltaR_g1b1  = t3phot1.DeltaR(t3jet1);
    float deltaR_g1b2  = t3phot1.DeltaR(t3jet2);
    float deltaR_g2b1  = t3phot2.DeltaR(t3jet1);
    float deltaR_g2b2  = t3phot2.DeltaR(t3jet2);
    float minDeltaR_gb = deltaR_g1b1;
    if ( deltaR_g1b2 < minDeltaR_gb ) minDeltaR_gb = deltaR_g1b2;
    if ( deltaR_g2b1 < minDeltaR_gb ) minDeltaR_gb = deltaR_g2b1;
    if ( deltaR_g2b2 < minDeltaR_gb ) minDeltaR_gb = deltaR_g2b2;

    // to study bias: random coupling of jets and gammas
    float theNumG  = myrand->Uniform(0,2);
    float theNumJ  = myrand->Uniform(0,2);
    float randMgj1 = 999.;
    float randMgj2 = 999.;
    if (theNumG<1 && theNumJ<1) { 
      randMgj1 = (t4phot1 + t4jet1).M(); 
      randMgj2 = (t4phot2 + t4jet2).M(); 
    };
    if (theNumG>1 && theNumJ<1) {
      randMgj1 = (t4phot2 + t4jet1).M(); 
      randMgj2 = (t4phot1 + t4jet2).M(); 
    };
    if (theNumG<1 && theNumJ>1) {
      randMgj1 = (t4phot1 + t4jet2).M(); 
      randMgj2 = (t4phot2 + t4jet1).M(); 
    };
    if (theNumG>1 && theNumJ>1) {
      randMgj1 = (t4phot2 + t4jet2).M(); 
      randMgj2 = (t4phot1 + t4jet1).M(); 
    }
    
    // further diphoton variables                                                                                                 
    float deltaEtaGG = t4phot1.Eta() - t4phot2.Eta();

    // radion 4vector                                                                                                             
    TLorentzVector t4Radion = t4diJet + t4diPhot;
    float radMass = t4Radion.M();

    // further variables: HT
    float HTgjet = ptJ1 + ptJ2 + ph1_pt + ph2_pt; 
    float HTjet  = ptJ1 + ptJ2;

    // control plots on basic kinematic variables after the two jets and two photons are preselected 
    // jets                                                                                                                       
    h1_ptjet0->Fill( ptJ1, weight_t );
    h1_ptjet1->Fill( ptJ2, weight_t );
    h1_runptjet0->Fill(ptJ1*120./t4diJet.M(), weight);
    h1_HTjet->Fill(HTjet);

    h1_etajet0->Fill( etaJ1, weight_t );
    h1_etajet1->Fill( etaJ2, weight_t );

    h1_ptDijet          -> Fill( ptJJ,  weight_t );
    h1_etaDijet         -> Fill( etaJJ,  weight_t );
    h1_deltaEtaJets     -> Fill( deltaEtaJJ, weight_t );
    h1_deltaFabsEtaJets -> Fill( fabs(deltaEtaJJ), weight_t );
    h1_deltaPhiJets     -> Fill( deltaPhiJJ, weight_t );

    // jets - gammas                                                                                                              
    h1_deltaR       -> Fill( deltaR_ggjj, weight_t );
    h1_deltaPhi     -> Fill( deltaPhi_ggjj, weight_t );
    h1_deltaEta     -> Fill( deltaEta_ggjj, weight_t );
    h1_ptRatio      -> Fill( t4diJet.Pt()/t4diPhot.Pt(), weight_t );
    h1_ptDifference -> Fill( t4diJet.Pt()-t4diPhot.Pt(), weight_t );
    h1_deltaR_gg    -> Fill( deltaR_gg, weight_t );
    h1_minDeltaR_gb -> Fill( minDeltaR_gb, weight_t );

    // gammas                                                                                                                     
    h1_ptDiphot       -> Fill( t4diPhot.Pt(), weight_t );
    h1_etaDiphot      -> Fill( t4diPhot.Eta(), weight_t );
    h1_deltaEtaDiphot -> Fill( deltaEtaGG, weight_t );


    // ------------------------------------------------------------                                                               
    // more refined analyses                                                                                                      

    // perform the kinfit                                                                                                        
    std::pair<TLorentzVector,TLorentzVector> jets_kinfitH = fitter_jetsH->fit(t4jet1, t4jet2);
    float chiSquareProbH = TMath::Prob(fitter_jetsH->getS(), fitter_jetsH->getNDF());

    // helicity angles                                                                                                            
    HelicityLikelihoodDiscriminant::HelicityAngles hangles;
    if( coin->Uniform(1.)<0.5 ) hangles = helicityDiscriminator->computeHelicityAngles(t4phot1, t4phot2, t4jet1, t4jet2);
    else                        hangles = helicityDiscriminator->computeHelicityAngles(t4phot1, t4phot2, t4jet1, t4jet2);
    float cosThetaStar = hangles.helCosThetaStar;
    h1_cosThetaStar -> Fill( hangles.helCosThetaStar, weight_t );

    // refitted jets                                                                                                              
    TLorentzVector jet1_kinfit, jet2_kinfit;
    jet1_kinfit = jets_kinfitH.first;
    jet2_kinfit = jets_kinfitH.second;
    TLorentzVector dijet_kinfit = jet1_kinfit + jet2_kinfit;
    TLorentzVector Vstar_kinfit = dijet_kinfit + t4diPhot;

    // control plots after kin fit                                                                                                
    h1_kinfit_chiSquareProbH->Fill( chiSquareProbH, weight_t );
    h1_mVstar_kinfit   -> Fill( Vstar_kinfit.M(),   weight_t );
    h1_mVstar          -> Fill( t4Radion.M(), weight_t );

    // -------------------------------------------------------------                                                              
    // counting the number of bjets to categorize the events                                                                      
    int btagCategory = -1;
    if (  isj1btagged==1 && isj2btagged==1 ) btagCategory = 2;
    if ( (isj1btagged==1 && isj2btagged==0) || (isj1btagged==0 && isj2btagged==1) ) btagCategory = 1;
    if (  isj1btagged==0 && isj2btagged==0 ) cout << "this should not happen!" << endl;          // chiara: da commentare solo x salvare anche 0bjets
    // if (  isj1btagged==0 && isj2btagged==0 ) btagCategory = 0;                                // chiara: da scommentare solo x salvare anche 0bjets  

    // categorize the events using gammas
    int myR9=-1;
    if(ph1_r9>.94 && ph2_r9>.94) myR9 = 1;
    if(ph1_r9<.94 || ph2_r9<.94) myR9 = 0;
    int theGammaCategory = -1;
    if (myR9==1) theGammaCategory = 0;
    if (myR9==0) theGammaCategory = 1;

    // total category combining photons / jets                                                                                    
    int theCategory = -1;
    if (btagCategory==2) theCategory = 0;
    if (btagCategory==1) theCategory = 1;


    // with kf
    float radMassKF = Vstar_kinfit.M();


    // -------------------------------------------------------------                                                              
    // further cuts for final selection ( = kinfit+regression so far )

    /* versione con ottimizzazione completa
    // mjj optimized window 
    if (fitToGG) {   
      if ( btagCategory==1 && (t4diJet.M()<85. || t4diJet.M()>155.) )  continue;
      if ( btagCategory==2 && (t4diJet.M()<110. || t4diJet.M()>140.) ) continue;
    }
    else if (fitToFourBodies) {   
      if ( btagCategory==1 && (t4diJet.M()<90. || t4diJet.M()>165.) )  continue;
      if ( btagCategory==2 && (t4diJet.M()<115. || t4diJet.M()>140.) ) continue;
    }

    // mggjj optimized window, and mgg (not optimized)
    if (fitToGG) {  
      if (useKinFit) {   
	// chiara, mass hypotesis: 300
	if ( btagCategory==1 && (radMassKF<290. || radMassKF>315.) ) continue;
	if ( btagCategory==2 && (radMassKF<280. || radMassKF>315.) ) continue;
	// chiara, mass hypotesis: 270 
	// if ( btagCategory==1 && (radMassKF<260. || radMassKF>285.) ) continue;
	// if ( btagCategory==2 && (radMassKF<265. || radMassKF>295.) ) continue;
	// chiara, mass hypotesis: 500
	// if ( btagCategory==1 && (radMassKF<475. || radMassKF>540.) ) continue;
	// if ( btagCategory==2 && (radMassKF<475. || radMassKF>540.) ) continue;  // boh. Per questo non ci sono eventi
	                                                                           // metto come 1btag, che sara' suboptimal ma almeno e' larga 
       } else {
	// chiara, mass hypotesis: 300
	if ( btagCategory==1 && (radMass<255. || radMass>330.) ) continue;
	if ( btagCategory==2 && (radMass<270. || radMass>325.) ) continue;
	// chiara, mass hypotesis: 270 
	// if ( btagCategory==1 && (radMass<225. || radMass>295.) ) continue;
	// if ( btagCategory==2 && (radMass<250. || radMass>290.) ) continue;
	// chiara, mass hypotesis: 500
	// if ( btagCategory==1 && (radMass<445. || radMass>535.) ) continue;
	// if ( btagCategory==2 && (radMass<445. || radMass>535.) ) continue;  // boh. Per questo non ci sono eventi
	                                                                           // metto come 1btag, che sara' suboptimal ma almeno e' larga 
       }
    }
    else if (fitToFourBodies) {   
      if( PhotonsMass<120 || PhotonsMass>130) continue;
    }
    */


    // versione per salvare statistica
    // mjj optimized window 
    if (fitToGG) {   
      if ( t4diJet.M()<85. || t4diJet.M()>155. )  continue;
    }
    else if (fitToFourBodies) {   
      if ( t4diJet.M()<90. || t4diJet.M()>165. )  continue;
    }

    // mggjj optimized window, and mgg (not optimized)
    if (fitToGG) {  
      // chiara, mass hypotesis: 260 
      // if ( radMass<225. || radMass>280. ) continue;
      // chiara, mass hypotesis: 270 
      // if ( radMass<225. || radMass>295. ) continue;
      // chiara, mass hypotesis: 300
      // if ( radMass<255. || radMass>330. ) continue;
      // chiara, mass hypotesis: 350
      // if ( radMass<310. || radMass>395. ) continue;
      // chiara, mass hypotesis: 400
      // if ( radMass<370. || radMass>440. ) continue;
      // chiara, mass hypotesis: 450
      // if ( radMass<410. || radMass>495. ) continue;
      // chiara, mass hypotesis: 500
      if ( radMass<445. || radMass>535. ) continue;
    }
    else if (fitToFourBodies) {   
      if( PhotonsMass<120 || PhotonsMass>130) continue;
    }

    // -------------------------------------------------------------   
    // invariant mass of jets by btag category                                                                                    
    if( btagCategory==0 ) {
      h1_mjj_0btag->Fill( invMassJJ, weight_t );
    } else if( btagCategory==1 ) {
      h1_mjj_1btag->Fill( invMassJJ, weight_t );
    } else {
      h1_mjj_2btag->Fill( invMassJJ, weight_t );
    }

    // invariant mass of radion by btag category                                                                                  
    h1_mggjj->Fill( radMass, weight_t );
    if( btagCategory==0 ) {
      h1_mggjj_0btag->Fill( radMass, weight_t );
    } else if( btagCategory==1 ) {
      h1_mggjj_1btag->Fill( radMass, weight_t );
    } else {
      h1_mggjj_2btag->Fill( radMass, weight_t );
    }

    // ---------------------------------------                                                                                    
    // for a comparison of masses with / wo kin fit                                                                               
    if( btagCategory==1 ) {
      h2_mggjj_vs_mjj_1btag        -> Fill( invMassJJ, radMass, weight_t );
      h2_mggjj_vs_mjj_kinfit_1btag -> Fill( invMassJJ, Vstar_kinfit.M(), weight_t );
      h1_mggjj_kinfit_1btag        -> Fill( Vstar_kinfit.M(), weight_t );
      h1_mggjj_nokinfit_1btag      -> Fill( radMass, weight_t );
    } else if (btagCategory==2 ) {
      h2_mggjj_vs_mjj_2btag        -> Fill( invMassJJ, radMass, weight_t );
      h2_mggjj_vs_mjj_kinfit_2btag -> Fill( invMassJJ, Vstar_kinfit.M(), weight_t );
      h1_mggjj_kinfit_2btag        -> Fill( Vstar_kinfit.M(), weight_t );
      h1_mggjj_nokinfit_2btag      -> Fill( radMass, weight_t );
    }


    // chiara
    // use btag SF to weight the events in case of signal.
    // They're dependent on being 1 or 2btag category events
    if (isSignal) {

      if (!isSignalG) {   // chiara: andra' rimosso quando i gravitoni avranno peso btag
	if (flavourJet1==0 || flavourJet2==0) continue;
      }	

      // new implementation
      float theBtagEvWeight;
      if (!isSignalG) theBtagEvWeight =  eventWeight_2jets("medium", btagSFJet1, btagSFJet2, btagEffJet1, btagEffJet2, btagOutJet1, btagOutJet2);    
      else theBtagEvWeight = 1.;   // chiara: andra' rimosso quando i gravitoni avranno peso btag

      // if the regression is used we need a factor x2 since we reject half of the sample 
      // if (usePhilReg)  weight_t = weight_t * theBtagEvWeight * 2. * 19706. / 1000. / nprocessed;
      // if (!usePhilReg) weight_t = weight_t * theBtagEvWeight * 19706. / 1000. / nprocessed;
      if (usePhilReg)  weight_t = weight_t * theBtagEvWeight * 2.;
      if (!usePhilReg) weight_t = weight_t * theBtagEvWeight;
    }


    // filling the tree for selected events                                                                                       
    massggnewvtx_t = PhotonsMass;
    ptphot1_t      = ph1_pt;
    runptphot1_t   = (120./PhotonsMass)*ph1_pt;
    ptphot2_t      = ph2_pt;
    etaphot1_t     = ph1_eta;
    etaphot2_t     = ph2_eta;
    cicphot1_t     = ph1_ciclevel;
    cicphot2_t     = ph2_ciclevel;
    r9phot1_t      = ph1_r9;
    r9phot2_t      = ph2_r9;
    ptgg_t         = dipho_pt;
    etagg_t        = dipho_eta;
    absetagg_t     = fabs(dipho_eta);
    njets_t        = v_puIdJets.size();
    ptcorrjet1_t   = ptJ1;
    ptcorrjet2_t   = ptJ2;
    runptcorrjet1_t = (120./invMassJJ)*ptJ1;
    etajet1_t      = etaJ1;
    etajet2_t      = etaJ2;
    btagSFjet1_t   = btagSFJet1;
    btagSFjet2_t   = btagSFJet2;
    flavourjet1_t  = flavourJet1;
    flavourjet2_t  = flavourJet2;
    emfjet1_t = emfJet1;
    emfjet2_t = emfJet2;
    chfjet1_t = chfJet1;
    chfjet2_t = chfJet2;
    nconstjet1_t = nconstJet1;
    nconstjet2_t = nconstJet2;
    vtxPtjet1_t = vtxPtJet1;
    vtxPtjet2_t = vtxPtJet2;
    vtx3dljet1_t = vtx3dlJet1;
    vtx3dljet2_t = vtx3dlJet2;
    deltaphijj_t   = deltaPhiJJ;
    deltaetajj_t   = deltaEtaJJ;
    invmassjet_t   = invMassJJ;   
    if (useKinFit)  invmassjetWithKF_t = dijet_kinfit.M();  
    if (!useKinFit) invmassjetWithKF_t = invMassJJ;
    ptjj_t         = ptJJ;
    etajj_t        = etaJJ;
    if (!useKinFit) massggjj_t = radMass;      
    if (useKinFit)  massggjj_t = Vstar_kinfit.M();      
    massggjjNoKF_t = radMass;
    deltaphijjgg_t = deltaPhi_ggjj;
    deltaetajjgg_t = deltaEta_ggjj;
    deltaRjjgg_t   = deltaR_ggjj;
    btagCategory_t = btagCategory;
    theCategory_t  = theCategory;
    theGammaCategory_t = theGammaCategory;
    if (bTaggerType_=="JP") {
      nbjets_loose_t  = v_looseJP.size();
      nbjets_medium_t = v_mediumJP.size();
      nbjets_tight_t  = v_tightJP.size();
    } else if (bTaggerType_=="CSV") {
      nbjets_loose_t  = v_looseCSV.size();
      nbjets_medium_t = v_mediumCSV.size();
      nbjets_tight_t  = v_tightCSV.size();
    } else {
      cout << "this btag algo does not exist" << endl;
    }
    chiSquareProbH_t = chiSquareProbH;
    absCosThetaStar_t  = fabs(cosThetaStar);
    nvtx_t      = nvtx;
    theVertex_t = vtx_ind;
    isj1btagged_t = isj1btagged;
    isj2btagged_t = isj2btagged;
    mjj_kin_t = dijet_kinfit.M();
    mggjj_kin_t = Vstar_kinfit.M();
    HT_gjet_t = HTgjet;
    HT_jet_t  = HTjet;
    deltaR_gg_t = deltaR_gg;
    minDeltaR_gb_t = minDeltaR_gb;
    randMgj1_t = randMgj1;
    randMgj2_t = randMgj2;

    // myTrees->Fill();

    TCVARS->Fill();                                                                                                            
    
  } // loop over entries  
  
  outFile_->cd();

  TCVARS->Write();                                                                                                             

  // myTrees->Write();

  h1_njets->Write();
  h1_nbjets_loose->Write();
  h1_nbjets_medium->Write();
  h1_nbjets_tight->Write();

  h1_ptphot0->Write();
  h1_ptphot1->Write();
  h1_runptphot0->Write();

  h1_mgg_preselG->Write();
  h1_mgg_preselJ->Write();
  h1_mjj_preselJ->Write();

  h1_ptjet0->Write();
  h1_ptjet1->Write();
  h1_runptjet0->Write();
  h1_HTjet->Write();
  h1_etajet0->Write();
  h1_etajet1->Write();

  h1_kinfit_chiSquareProbH->Write();

  h1_mjj_0btag->Write();
  h1_mjj_1btag->Write();
  h1_mjj_2btag->Write();

  h1_mggjj->Write();
  h1_mggjj_0btag->Write();
  h1_mggjj_1btag->Write();
  h1_mggjj_2btag->Write();

  h1_ptDiphot->Write();
  h1_etaDiphot->Write();
  h1_deltaEtaDiphot->Write();

  h1_deltaR->Write();
  h1_deltaPhi->Write();
  h1_deltaEta->Write();
  h1_ptDijet->Write();
  h1_etaDijet->Write();
  h1_ptRatio->Write();
  h1_ptDifference->Write();
  h1_deltaR_gg->Write();
  h1_minDeltaR_gb->Write();

  h1_deltaPhiJets->Write();
  h1_deltaEtaJets->Write();
  h1_deltaFabsEtaJets->Write();

  h1_cosThetaStar->Write();

  h1_mVstar->Write();
  h1_mVstar_kinfit->Write();

  h2_mggjj_vs_mjj_1btag        -> Write();
  h2_mggjj_vs_mjj_kinfit_1btag -> Write();
  h1_mggjj_kinfit_1btag        -> Write();
  h1_mggjj_nokinfit_1btag      -> Write();

  h2_mggjj_vs_mjj_2btag        -> Write();
  h2_mggjj_vs_mjj_kinfit_2btag -> Write();
  h1_mggjj_kinfit_2btag        -> Write();
  h1_mggjj_nokinfit_2btag      -> Write();

} // end finalize

Int_t fillPlot2012_radion_commonNtp::GetEntry(Long64_t entry) {

  if (!tree_) return 0;
  return tree_->GetEntry(entry);
}

Long64_t fillPlot2012_radion_commonNtp::LoadTree(Long64_t entry) {

  if (!tree_) return -5;
  Long64_t centry = tree_->LoadTree(entry);
  if (centry < 0) return centry;
  if (tree_->GetTreeNumber() != fCurrent) {
    fCurrent = tree_->GetTreeNumber();
  }
  return centry;
}

void fillPlot2012_radion_commonNtp::Init() {

  // Set branch addresses and branch pointers
  fCurrent = -1;
  tree_->SetMakeClass(1);

  tree_->SetBranchAddress("itype", &itype, &b_itype);
  tree_->SetBranchAddress("run", &run, &b_run);
  tree_->SetBranchAddress("lumis", &lumis, &b_lumis);
  tree_->SetBranchAddress("event", &event, &b_event);
  tree_->SetBranchAddress("weight", &weight, &b_weight);
  tree_->SetBranchAddress("evweight", &evweight, &b_evweight);
  tree_->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  tree_->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
  tree_->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree_->SetBranchAddress("rho", &rho, &b_rho);
  tree_->SetBranchAddress("category", &category, &b_category);
  tree_->SetBranchAddress("met_pfmet", &met_pfmet, &b_met_pfmet);
  tree_->SetBranchAddress("met_phi_pfmet", &met_phi_pfmet, &b_met_phi_pfmet);
  tree_->SetBranchAddress("met_corr_pfmet", &met_corr_pfmet, &b_met_corr_pfmet);
  tree_->SetBranchAddress("met_corr_phi_pfmet", &met_corr_phi_pfmet, &b_met_corr_phi_pfmet);
  tree_->SetBranchAddress("met_corr_eta_pfmet", &met_corr_eta_pfmet, &b_met_corr_eta_pfmet);
  tree_->SetBranchAddress("met_corr_e_pfmet", &met_corr_e_pfmet, &b_met_corr_e_pfmet);
  tree_->SetBranchAddress("ph1_e", &ph1_e, &b_ph1_e);
  tree_->SetBranchAddress("ph2_e", &ph2_e, &b_ph2_e);
  tree_->SetBranchAddress("ph1_pt", &ph1_pt, &b_ph1_pt);
  tree_->SetBranchAddress("ph2_pt", &ph2_pt, &b_ph2_pt);
  tree_->SetBranchAddress("ph1_phi", &ph1_phi, &b_ph1_phi);
  tree_->SetBranchAddress("ph2_phi", &ph2_phi, &b_ph2_phi);
  tree_->SetBranchAddress("ph1_eta", &ph1_eta, &b_ph1_eta);
  tree_->SetBranchAddress("ph2_eta", &ph2_eta, &b_ph2_eta);
  tree_->SetBranchAddress("ph1_r9", &ph1_r9, &b_ph1_r9);
  tree_->SetBranchAddress("ph2_r9", &ph2_r9, &b_ph2_r9);
  tree_->SetBranchAddress("ph1_isPrompt", &ph1_isPrompt, &b_ph1_isPrompt);
  tree_->SetBranchAddress("ph2_isPrompt", &ph2_isPrompt, &b_ph2_isPrompt);
  tree_->SetBranchAddress("ph1_SCEta", &ph1_SCEta, &b_ph1_SCEta);
  tree_->SetBranchAddress("ph2_SCEta", &ph2_SCEta, &b_ph2_SCEta);
  tree_->SetBranchAddress("ph1_SCPhi", &ph1_SCPhi, &b_ph1_SCPhi);
  tree_->SetBranchAddress("ph2_SCPhi", &ph2_SCPhi, &b_ph2_SCPhi);
  tree_->SetBranchAddress("ph1_hoe", &ph1_hoe, &b_ph1_hoe);
  tree_->SetBranchAddress("ph2_hoe", &ph2_hoe, &b_ph2_hoe);
  tree_->SetBranchAddress("ph1_sieie", &ph1_sieie, &b_ph1_sieie);
  tree_->SetBranchAddress("ph2_sieie", &ph2_sieie, &b_ph2_sieie);
  tree_->SetBranchAddress("ph1_pfchargedisogood03", &ph1_pfchargedisogood03, &b_ph1_pfchargedisogood03);
  tree_->SetBranchAddress("ph2_pfchargedisogood03", &ph2_pfchargedisogood03, &b_ph2_pfchargedisogood03);
  tree_->SetBranchAddress("ph1_pfchargedisobad04", &ph1_pfchargedisobad04, &b_ph1_pfchargedisobad04);
  tree_->SetBranchAddress("ph2_pfchargedisobad04", &ph2_pfchargedisobad04, &b_ph2_pfchargedisobad04);
  tree_->SetBranchAddress("ph1_etawidth", &ph1_etawidth, &b_ph1_etawidth);
  tree_->SetBranchAddress("ph2_etawidth", &ph2_etawidth, &b_ph2_etawidth);
  tree_->SetBranchAddress("ph1_phiwidth", &ph1_phiwidth, &b_ph1_phiwidth);
  tree_->SetBranchAddress("ph2_phiwidth", &ph2_phiwidth, &b_ph2_phiwidth);
  tree_->SetBranchAddress("ph1_eseffssqrt", &ph1_eseffssqrt, &b_ph1_eseffssqrt);
  tree_->SetBranchAddress("ph2_eseffssqrt", &ph2_eseffssqrt, &b_ph2_eseffssqrt);
  tree_->SetBranchAddress("ph1_pfchargedisobad03", &ph1_pfchargedisobad03, &b_ph1_pfchargedisobad03);
  tree_->SetBranchAddress("ph2_pfchargedisobad03", &ph2_pfchargedisobad03, &b_ph2_pfchargedisobad03);
  tree_->SetBranchAddress("ph1_sieip", &ph1_sieip, &b_ph1_sieip);
  tree_->SetBranchAddress("ph2_sieip", &ph2_sieip, &b_ph2_sieip);
  tree_->SetBranchAddress("ph1_sipip", &ph1_sipip, &b_ph1_sipip);
  tree_->SetBranchAddress("ph2_sipip", &ph2_sipip, &b_ph2_sipip);
  tree_->SetBranchAddress("ph1_ecaliso", &ph1_ecaliso, &b_ph1_ecaliso);
  tree_->SetBranchAddress("ph2_ecaliso", &ph2_ecaliso, &b_ph2_ecaliso);
  tree_->SetBranchAddress("ph1_ecalisobad", &ph1_ecalisobad, &b_ph1_ecalisobad);
  tree_->SetBranchAddress("ph2_ecalisobad", &ph2_ecalisobad, &b_ph2_ecalisobad);
  tree_->SetBranchAddress("ph1_badvtx_Et", &ph1_badvtx_Et, &b_ph1_badvtx_Et);
  tree_->SetBranchAddress("ph2_badvtx_Et", &ph2_badvtx_Et, &b_ph2_badvtx_Et);
  tree_->SetBranchAddress("ph1_isconv", &ph1_isconv, &b_ph1_isconv);
  tree_->SetBranchAddress("ph2_isconv", &ph2_isconv, &b_ph2_isconv);
  tree_->SetBranchAddress("ph1_ciclevel", &ph1_ciclevel, &b_ph1_ciclevel);
  tree_->SetBranchAddress("ph2_ciclevel", &ph2_ciclevel, &b_ph2_ciclevel);
  tree_->SetBranchAddress("ph1_sigmaEoE", &ph1_sigmaEoE, &b_ph1_sigmaEoE);
  tree_->SetBranchAddress("ph2_sigmaEoE", &ph2_sigmaEoE, &b_ph2_sigmaEoE);
  tree_->SetBranchAddress("ph1_ptoM", &ph1_ptoM, &b_ph1_ptoM);
  tree_->SetBranchAddress("ph2_ptoM", &ph2_ptoM, &b_ph2_ptoM);
  tree_->SetBranchAddress("ph1_isEB", &ph1_isEB, &b_ph1_isEB);
  tree_->SetBranchAddress("ph2_isEB", &ph2_isEB, &b_ph2_isEB);
  tree_->SetBranchAddress("ph1_s4ratio", &ph1_s4ratio, &b_ph1_s4ratio);
  tree_->SetBranchAddress("ph2_s4ratio", &ph2_s4ratio, &b_ph2_s4ratio);
  tree_->SetBranchAddress("ph1_e3x3", &ph1_e3x3, &b_ph1_e3x3);
  tree_->SetBranchAddress("ph2_e3x3", &ph2_e3x3, &b_ph2_e3x3);
  tree_->SetBranchAddress("ph1_e5x5", &ph1_e5x5, &b_ph1_e5x5);
  tree_->SetBranchAddress("ph2_e5x5", &ph2_e5x5, &b_ph2_e5x5);
  tree_->SetBranchAddress("PhotonsMass", &PhotonsMass, &b_PhotonsMass);
  tree_->SetBranchAddress("dipho_E", &dipho_E, &b_dipho_E);
  tree_->SetBranchAddress("dipho_pt", &dipho_pt, &b_dipho_pt);
  tree_->SetBranchAddress("dipho_eta", &dipho_eta, &b_dipho_eta);
  tree_->SetBranchAddress("dipho_phi", &dipho_phi, &b_dipho_phi);
  tree_->SetBranchAddress("dipho_cosThetaStar_CS", &dipho_cosThetaStar_CS, &b_dipho_cosThetaStar_CS);
  tree_->SetBranchAddress("dipho_tanhYStar", &dipho_tanhYStar, &b_dipho_tanhYStar);
  tree_->SetBranchAddress("dipho_Y", &dipho_Y, &b_dipho_Y);
  tree_->SetBranchAddress("vtx_ind", &vtx_ind, &b_vtx_ind);
  tree_->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
  tree_->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
  tree_->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
  tree_->SetBranchAddress("vtx_mva", &vtx_mva, &b_vtx_mva);
  tree_->SetBranchAddress("vtx_mva_2", &vtx_mva_2, &b_vtx_mva_2);
  tree_->SetBranchAddress("vtx_mva_3", &vtx_mva_3, &b_vtx_mva_3);
  tree_->SetBranchAddress("vtx_ptbal", &vtx_ptbal, &b_vtx_ptbal);
  tree_->SetBranchAddress("vtx_ptasym", &vtx_ptasym, &b_vtx_ptasym);
  tree_->SetBranchAddress("vtx_logsumpt2", &vtx_logsumpt2, &b_vtx_logsumpt2);
  tree_->SetBranchAddress("vtx_pulltoconv", &vtx_pulltoconv, &b_vtx_pulltoconv);
  tree_->SetBranchAddress("vtx_prob", &vtx_prob, &b_vtx_prob);
  tree_->SetBranchAddress("njets_passing_kLooseID", &njets_passing_kLooseID, &b_njets_passing_kLooseID);
  tree_->SetBranchAddress("njets_passing_kLooseID_and_CSVL", &njets_passing_kLooseID_and_CSVL, &b_njets_passing_kLooseID_and_CSVL);
  tree_->SetBranchAddress("njets_passing_kLooseID_and_CSVM", &njets_passing_kLooseID_and_CSVM, &b_njets_passing_kLooseID_and_CSVM);
  tree_->SetBranchAddress("njets_passing_kLooseID_and_CSVT", &njets_passing_kLooseID_and_CSVT, &b_njets_passing_kLooseID_and_CSVT);
  tree_->SetBranchAddress("j1_e", &j1_e, &b_j1_e);
  tree_->SetBranchAddress("j1_pt", &j1_pt, &b_j1_pt);
  tree_->SetBranchAddress("j1_phi", &j1_phi, &b_j1_phi);
  tree_->SetBranchAddress("j1_eta", &j1_eta, &b_j1_eta);
  tree_->SetBranchAddress("j1_beta", &j1_beta, &b_j1_beta);
  tree_->SetBranchAddress("j1_betaStar", &j1_betaStar, &b_j1_betaStar);
  tree_->SetBranchAddress("j1_betaStarClassic", &j1_betaStarClassic, &b_j1_betaStarClassic);
  tree_->SetBranchAddress("j1_dR2Mean", &j1_dR2Mean, &b_j1_dR2Mean);
  tree_->SetBranchAddress("j1_csvBtag", &j1_csvBtag, &b_j1_csvBtag);
  tree_->SetBranchAddress("j1_csvMvaBtag", &j1_csvMvaBtag, &b_j1_csvMvaBtag);
  tree_->SetBranchAddress("j1_jetProbBtag", &j1_jetProbBtag, &b_j1_jetProbBtag);
  tree_->SetBranchAddress("j1_tcheBtag", &j1_tcheBtag, &b_j1_tcheBtag);
  tree_->SetBranchAddress("j1_radionMatched", &j1_radionMatched, &b_j1_radionMatched);
  tree_->SetBranchAddress("j1_ptD", &j1_ptD, &b_j1_ptD);
  tree_->SetBranchAddress("j1_nSecondaryVertices", &j1_nSecondaryVertices, &b_j1_nSecondaryVertices);
  tree_->SetBranchAddress("j1_secVtxPt", &j1_secVtxPt, &b_j1_secVtxPt);
  tree_->SetBranchAddress("j1_secVtx3dL", &j1_secVtx3dL, &b_j1_secVtx3dL);
  tree_->SetBranchAddress("j1_secVtx3deL", &j1_secVtx3deL, &b_j1_secVtx3deL);
  tree_->SetBranchAddress("j1_emfrac", &j1_emfrac, &b_j1_emfrac);
  tree_->SetBranchAddress("j1_hadfrac", &j1_hadfrac, &b_j1_hadfrac);
  tree_->SetBranchAddress("j1_ntk", &j1_ntk, &b_j1_ntk);
  tree_->SetBranchAddress("j1_nNeutrals", &j1_nNeutrals, &b_j1_nNeutrals);
  tree_->SetBranchAddress("j1_nCharged", &j1_nCharged, &b_j1_nCharged);
  tree_->SetBranchAddress("j1_axis1", &j1_axis1, &b_j1_axis1);
  tree_->SetBranchAddress("j1_axis2", &j1_axis2, &b_j1_axis2);
  tree_->SetBranchAddress("j1_pull", &j1_pull, &b_j1_pull);
  tree_->SetBranchAddress("j1_flavour", &j1_flavour, &b_j1_flavour);
  tree_->SetBranchAddress("j1_btagSF", &j1_btagSF, &b_j1_btagSF);
  tree_->SetBranchAddress("j1_btagSFErrorUp", &j1_btagSFErrorUp, &b_j1_btagSFErrorUp);
  tree_->SetBranchAddress("j1_btagSFErrorDown", &j1_btagSFErrorDown, &b_j1_btagSFErrorDown);
  tree_->SetBranchAddress("j1_btagEff", &j1_btagEff, &b_j1_btagEff);
  tree_->SetBranchAddress("j1_btagEffError", &j1_btagEffError, &b_j1_btagEffError);
  tree_->SetBranchAddress("j2_e", &j2_e, &b_j2_e);
  tree_->SetBranchAddress("j2_pt", &j2_pt, &b_j2_pt);
  tree_->SetBranchAddress("j2_phi", &j2_phi, &b_j2_phi);
  tree_->SetBranchAddress("j2_eta", &j2_eta, &b_j2_eta);
  tree_->SetBranchAddress("j2_beta", &j2_beta, &b_j2_beta);
  tree_->SetBranchAddress("j2_betaStar", &j2_betaStar, &b_j2_betaStar);
  tree_->SetBranchAddress("j2_betaStarClassic", &j2_betaStarClassic, &b_j2_betaStarClassic);
  tree_->SetBranchAddress("j2_dR2Mean", &j2_dR2Mean, &b_j2_dR2Mean);
  tree_->SetBranchAddress("j2_csvBtag", &j2_csvBtag, &b_j2_csvBtag);
  tree_->SetBranchAddress("j2_csvMvaBtag", &j2_csvMvaBtag, &b_j2_csvMvaBtag);
  tree_->SetBranchAddress("j2_jetProbBtag", &j2_jetProbBtag, &b_j2_jetProbBtag);
  tree_->SetBranchAddress("j2_tcheBtag", &j2_tcheBtag, &b_j2_tcheBtag);
  tree_->SetBranchAddress("j2_radionMatched", &j2_radionMatched, &b_j2_radionMatched);
  tree_->SetBranchAddress("j2_ptD", &j2_ptD, &b_j2_ptD);
  tree_->SetBranchAddress("j2_nSecondaryVertices", &j2_nSecondaryVertices, &b_j2_nSecondaryVertices);
  tree_->SetBranchAddress("j2_secVtxPt", &j2_secVtxPt, &b_j2_secVtxPt);
  tree_->SetBranchAddress("j2_secVtx3dL", &j2_secVtx3dL, &b_j2_secVtx3dL);
  tree_->SetBranchAddress("j2_secVtx3deL", &j2_secVtx3deL, &b_j2_secVtx3deL);
  tree_->SetBranchAddress("j2_emfrac", &j2_emfrac, &b_j2_emfrac);
  tree_->SetBranchAddress("j2_hadfrac", &j2_hadfrac, &b_j2_hadfrac);
  tree_->SetBranchAddress("j2_ntk", &j2_ntk, &b_j2_ntk);
  tree_->SetBranchAddress("j2_nNeutrals", &j2_nNeutrals, &b_j2_nNeutrals);
  tree_->SetBranchAddress("j2_nCharged", &j2_nCharged, &b_j2_nCharged);
  tree_->SetBranchAddress("j2_axis1", &j2_axis1, &b_j2_axis1);
  tree_->SetBranchAddress("j2_axis2", &j2_axis2, &b_j2_axis2);
  tree_->SetBranchAddress("j2_pull", &j2_pull, &b_j2_pull);
  tree_->SetBranchAddress("j2_flavour", &j2_flavour, &b_j2_flavour);
  tree_->SetBranchAddress("j2_btagSF", &j2_btagSF, &b_j2_btagSF);
  tree_->SetBranchAddress("j2_btagSFErrorUp", &j2_btagSFErrorUp, &b_j2_btagSFErrorUp);
  tree_->SetBranchAddress("j2_btagSFErrorDown", &j2_btagSFErrorDown, &b_j2_btagSFErrorDown);
  tree_->SetBranchAddress("j2_btagEff", &j2_btagEff, &b_j2_btagEff);
  tree_->SetBranchAddress("j2_btagEffError", &j2_btagEffError, &b_j2_btagEffError);
  tree_->SetBranchAddress("j3_e", &j3_e, &b_j3_e);
  tree_->SetBranchAddress("j3_pt", &j3_pt, &b_j3_pt);
  tree_->SetBranchAddress("j3_phi", &j3_phi, &b_j3_phi);
  tree_->SetBranchAddress("j3_eta", &j3_eta, &b_j3_eta);
  tree_->SetBranchAddress("j3_beta", &j3_beta, &b_j3_beta);
  tree_->SetBranchAddress("j3_betaStar", &j3_betaStar, &b_j3_betaStar);
  tree_->SetBranchAddress("j3_betaStarClassic", &j3_betaStarClassic, &b_j3_betaStarClassic);
  tree_->SetBranchAddress("j3_dR2Mean", &j3_dR2Mean, &b_j3_dR2Mean);
  tree_->SetBranchAddress("j3_csvBtag", &j3_csvBtag, &b_j3_csvBtag);
  tree_->SetBranchAddress("j3_csvMvaBtag", &j3_csvMvaBtag, &b_j3_csvMvaBtag);
  tree_->SetBranchAddress("j3_jetProbBtag", &j3_jetProbBtag, &b_j3_jetProbBtag);
  tree_->SetBranchAddress("j3_tcheBtag", &j3_tcheBtag, &b_j3_tcheBtag);
  tree_->SetBranchAddress("j3_radionMatched", &j3_radionMatched, &b_j3_radionMatched);
  tree_->SetBranchAddress("j3_ptD", &j3_ptD, &b_j3_ptD);
  tree_->SetBranchAddress("j3_nSecondaryVertices", &j3_nSecondaryVertices, &b_j3_nSecondaryVertices);
  tree_->SetBranchAddress("j3_secVtxPt", &j3_secVtxPt, &b_j3_secVtxPt);
  tree_->SetBranchAddress("j3_secVtx3dL", &j3_secVtx3dL, &b_j3_secVtx3dL);
  tree_->SetBranchAddress("j3_secVtx3deL", &j3_secVtx3deL, &b_j3_secVtx3deL);
  tree_->SetBranchAddress("j3_emfrac", &j3_emfrac, &b_j3_emfrac);
  tree_->SetBranchAddress("j3_hadfrac", &j3_hadfrac, &b_j3_hadfrac);
  tree_->SetBranchAddress("j3_ntk", &j3_ntk, &b_j3_ntk);
  tree_->SetBranchAddress("j3_nNeutrals", &j3_nNeutrals, &b_j3_nNeutrals);
  tree_->SetBranchAddress("j3_nCharged", &j3_nCharged, &b_j3_nCharged);
  tree_->SetBranchAddress("j3_axis1", &j3_axis1, &b_j3_axis1);
  tree_->SetBranchAddress("j3_axis2", &j3_axis2, &b_j3_axis2);
  tree_->SetBranchAddress("j3_pull", &j3_pull, &b_j3_pull);
  tree_->SetBranchAddress("j3_flavour", &j3_flavour, &b_j3_flavour);
  tree_->SetBranchAddress("j3_btagSF", &j3_btagSF, &b_j3_btagSF);
  tree_->SetBranchAddress("j3_btagSFErrorUp", &j3_btagSFErrorUp, &b_j3_btagSFErrorUp);
  tree_->SetBranchAddress("j3_btagSFErrorDown", &j3_btagSFErrorDown, &b_j3_btagSFErrorDown);
  tree_->SetBranchAddress("j3_btagEff", &j3_btagEff, &b_j3_btagEff);
  tree_->SetBranchAddress("j3_btagEffError", &j3_btagEffError, &b_j3_btagEffError);
  tree_->SetBranchAddress("j4_e", &j4_e, &b_j4_e);
  tree_->SetBranchAddress("j4_pt", &j4_pt, &b_j4_pt);
  tree_->SetBranchAddress("j4_phi", &j4_phi, &b_j4_phi);
  tree_->SetBranchAddress("j4_eta", &j4_eta, &b_j4_eta);
  tree_->SetBranchAddress("j4_beta", &j4_beta, &b_j4_beta);
  tree_->SetBranchAddress("j4_betaStar", &j4_betaStar, &b_j4_betaStar);
  tree_->SetBranchAddress("j4_betaStarClassic", &j4_betaStarClassic, &b_j4_betaStarClassic);
  tree_->SetBranchAddress("j4_dR2Mean", &j4_dR2Mean, &b_j4_dR2Mean);
  tree_->SetBranchAddress("j4_csvBtag", &j4_csvBtag, &b_j4_csvBtag);
  tree_->SetBranchAddress("j4_csvMvaBtag", &j4_csvMvaBtag, &b_j4_csvMvaBtag);
  tree_->SetBranchAddress("j4_jetProbBtag", &j4_jetProbBtag, &b_j4_jetProbBtag);
  tree_->SetBranchAddress("j4_tcheBtag", &j4_tcheBtag, &b_j4_tcheBtag);
  tree_->SetBranchAddress("j4_radionMatched", &j4_radionMatched, &b_j4_radionMatched);
  tree_->SetBranchAddress("j4_ptD", &j4_ptD, &b_j4_ptD);
  tree_->SetBranchAddress("j4_nSecondaryVertices", &j4_nSecondaryVertices, &b_j4_nSecondaryVertices);
  tree_->SetBranchAddress("j4_secVtxPt", &j4_secVtxPt, &b_j4_secVtxPt);
  tree_->SetBranchAddress("j4_secVtx3dL", &j4_secVtx3dL, &b_j4_secVtx3dL);
  tree_->SetBranchAddress("j4_secVtx3deL", &j4_secVtx3deL, &b_j4_secVtx3deL);
  tree_->SetBranchAddress("j4_emfrac", &j4_emfrac, &b_j4_emfrac);
  tree_->SetBranchAddress("j4_hadfrac", &j4_hadfrac, &b_j4_hadfrac);
  tree_->SetBranchAddress("j4_ntk", &j4_ntk, &b_j4_ntk);
  tree_->SetBranchAddress("j4_nNeutrals", &j4_nNeutrals, &b_j4_nNeutrals);
  tree_->SetBranchAddress("j4_nCharged", &j4_nCharged, &b_j4_nCharged);
  tree_->SetBranchAddress("j4_axis1", &j4_axis1, &b_j4_axis1);
  tree_->SetBranchAddress("j4_axis2", &j4_axis2, &b_j4_axis2);
  tree_->SetBranchAddress("j4_pull", &j4_pull, &b_j4_pull);
  tree_->SetBranchAddress("j4_flavour", &j4_flavour, &b_j4_flavour);
  tree_->SetBranchAddress("j4_btagSF", &j4_btagSF, &b_j4_btagSF);
  tree_->SetBranchAddress("j4_btagSFErrorUp", &j4_btagSFErrorUp, &b_j4_btagSFErrorUp);
  tree_->SetBranchAddress("j4_btagSFErrorDown", &j4_btagSFErrorDown, &b_j4_btagSFErrorDown);
  tree_->SetBranchAddress("j4_btagEff", &j4_btagEff, &b_j4_btagEff);
  tree_->SetBranchAddress("j4_btagEffError", &j4_btagEffError, &b_j4_btagEffError);
  tree_->SetBranchAddress("JetsMass", &JetsMass, &b_JetsMass);
  tree_->SetBranchAddress("dijet_E", &dijet_E, &b_dijet_E);
  tree_->SetBranchAddress("dijet_Pt", &dijet_Pt, &b_dijet_Pt);
  tree_->SetBranchAddress("dijet_Eta", &dijet_Eta, &b_dijet_Eta);
  tree_->SetBranchAddress("dijet_Phi", &dijet_Phi, &b_dijet_Phi);
  tree_->SetBranchAddress("RadMass", &RadMass, &b_RadMass);
  tree_->SetBranchAddress("radion_E", &radion_E, &b_radion_E);
  tree_->SetBranchAddress("radion_Pt", &radion_Pt, &b_radion_Pt);
  tree_->SetBranchAddress("radion_Eta", &radion_Eta, &b_radion_Eta);
  tree_->SetBranchAddress("radion_Phi", &radion_Phi, &b_radion_Phi);
  tree_->SetBranchAddress("gr_radion_p4_pt", &gr_radion_p4_pt, &b_gr_radion_p4_pt);
  tree_->SetBranchAddress("gr_radion_p4_eta", &gr_radion_p4_eta, &b_gr_radion_p4_eta);
  tree_->SetBranchAddress("gr_radion_p4_phi", &gr_radion_p4_phi, &b_gr_radion_p4_phi);
  tree_->SetBranchAddress("gr_radion_p4_mass", &gr_radion_p4_mass, &b_gr_radion_p4_mass);
  tree_->SetBranchAddress("gr_hgg_p4_pt", &gr_hgg_p4_pt, &b_gr_hgg_p4_pt);
  tree_->SetBranchAddress("gr_hgg_p4_eta", &gr_hgg_p4_eta, &b_gr_hgg_p4_eta);
  tree_->SetBranchAddress("gr_hgg_p4_phi", &gr_hgg_p4_phi, &b_gr_hgg_p4_phi);
  tree_->SetBranchAddress("gr_hgg_p4_mass", &gr_hgg_p4_mass, &b_gr_hgg_p4_mass);
  tree_->SetBranchAddress("gr_hbb_p4_pt", &gr_hbb_p4_pt, &b_gr_hbb_p4_pt);
  tree_->SetBranchAddress("gr_hbb_p4_eta", &gr_hbb_p4_eta, &b_gr_hbb_p4_eta);
  tree_->SetBranchAddress("gr_hbb_p4_phi", &gr_hbb_p4_phi, &b_gr_hbb_p4_phi);
  tree_->SetBranchAddress("gr_hbb_p4_mass", &gr_hbb_p4_mass, &b_gr_hbb_p4_mass);
  tree_->SetBranchAddress("gr_g1_p4_pt", &gr_g1_p4_pt, &b_gr_g1_p4_pt);
  tree_->SetBranchAddress("gr_g1_p4_eta", &gr_g1_p4_eta, &b_gr_g1_p4_eta);
  tree_->SetBranchAddress("gr_g1_p4_phi", &gr_g1_p4_phi, &b_gr_g1_p4_phi);
  tree_->SetBranchAddress("gr_g1_p4_mass", &gr_g1_p4_mass, &b_gr_g1_p4_mass);
  tree_->SetBranchAddress("gr_g2_p4_pt", &gr_g2_p4_pt, &b_gr_g2_p4_pt);
  tree_->SetBranchAddress("gr_g2_p4_eta", &gr_g2_p4_eta, &b_gr_g2_p4_eta);
  tree_->SetBranchAddress("gr_g2_p4_phi", &gr_g2_p4_phi, &b_gr_g2_p4_phi);
  tree_->SetBranchAddress("gr_g2_p4_mass", &gr_g2_p4_mass, &b_gr_g2_p4_mass);
  tree_->SetBranchAddress("gr_b1_p4_pt", &gr_b1_p4_pt, &b_gr_b1_p4_pt);
  tree_->SetBranchAddress("gr_b1_p4_eta", &gr_b1_p4_eta, &b_gr_b1_p4_eta);
  tree_->SetBranchAddress("gr_b1_p4_phi", &gr_b1_p4_phi, &b_gr_b1_p4_phi);
  tree_->SetBranchAddress("gr_b1_p4_mass", &gr_b1_p4_mass, &b_gr_b1_p4_mass);
  tree_->SetBranchAddress("gr_b2_p4_pt", &gr_b2_p4_pt, &b_gr_b2_p4_pt);
  tree_->SetBranchAddress("gr_b2_p4_eta", &gr_b2_p4_eta, &b_gr_b2_p4_eta);
  tree_->SetBranchAddress("gr_b2_p4_phi", &gr_b2_p4_phi, &b_gr_b2_p4_phi);
  tree_->SetBranchAddress("gr_b2_p4_mass", &gr_b2_p4_mass, &b_gr_b2_p4_mass);
  tree_->SetBranchAddress("gr_j1_p4_pt", &gr_j1_p4_pt, &b_gr_j1_p4_pt);
  tree_->SetBranchAddress("gr_j1_p4_eta", &gr_j1_p4_eta, &b_gr_j1_p4_eta);
  tree_->SetBranchAddress("gr_j1_p4_phi", &gr_j1_p4_phi, &b_gr_j1_p4_phi);
  tree_->SetBranchAddress("gr_j1_p4_mass", &gr_j1_p4_mass, &b_gr_j1_p4_mass);
  tree_->SetBranchAddress("gr_j2_p4_pt", &gr_j2_p4_pt, &b_gr_j2_p4_pt);
  tree_->SetBranchAddress("gr_j2_p4_eta", &gr_j2_p4_eta, &b_gr_j2_p4_eta);
  tree_->SetBranchAddress("gr_j2_p4_phi", &gr_j2_p4_phi, &b_gr_j2_p4_phi);
  tree_->SetBranchAddress("gr_j2_p4_mass", &gr_j2_p4_mass, &b_gr_j2_p4_mass);
}

void fillPlot2012_radion_commonNtp::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;
  
  // default values                                                                                                               
  cicselection = 4;
  
  ptphot1cut = 40.;
  ptphot2cut = 30.;
  
  ptjetacccut  = 25.;
  etajetacccut = 2.5;
}

int fillPlot2012_radion_commonNtp::myHighestPtJet(std::vector<int> jets) {

  float firstJpt = -898.;
  int firstJ     = -899;

  for (int ij=0; ij<int(jets.size());ij++) {
    int jetIndex = jets[ij];
    if (ptcorrjet[jetIndex]>firstJpt) {
      firstJpt  = ptcorrjet[jetIndex];
      firstJ    = jetIndex;
    }
  }

  return firstJ;
}

bool fillPlot2012_radion_commonNtp::passCutBasedJetId(int jet) {

  bool isGood = true;

  float thebetastarjet[4], thermsjet[4];  
  thebetastarjet[0] = j1_betaStarClassic;
  thebetastarjet[1] = j2_betaStarClassic;
  thebetastarjet[2] = j3_betaStarClassic;
  thebetastarjet[3] = j4_betaStarClassic;
  thermsjet[0]      = j1_dR2Mean;
  thermsjet[1]      = j2_dR2Mean;
  thermsjet[2]      = j3_dR2Mean;
  thermsjet[3]      = j4_dR2Mean;

  if ( ptcorrjet[jet]<-900) isGood = false;

  if ( fabs(etajet[jet]) < 2.5 ) {
    if ( thebetastarjet[jet] > 0.2 * log( nvtx - 0.64) )  isGood = false;
    if (thermsjet[jet] > 0.06)                            isGood = false;
  } else if (fabs(etajet[jet]) < 3.){
    if ( thermsjet[jet] > 0.05)  isGood =false;
  } else {
    if ( thermsjet[jet] > 0.055) isGood =false;
  }

  return isGood;
}

float fillPlot2012_radion_commonNtp::eventWeight_2jets(std::string WP, const float& j1_SF, const float& j2_SF, const float& j1_eff, const float& j2_eff, const float& j1_csvBtag, const float& j2_csvBtag){

  float csvWP = 0.;
  
  if(WP == "loose")  csvWP = 0.244;
  if(WP == "medium") csvWP = 0.679;
  if(WP == "tight")  csvWP = 0.898;
  
  float weight1 = 1.;
  float weight2 = 1.;

  if(j1_csvBtag > csvWP) weight1 = j1_SF;
  else weight1 = (1-j1_SF*j1_eff)/(1-j1_eff);
  
  if(j2_csvBtag > csvWP) weight2 = j2_SF;
  else weight2 = (1-j2_SF*j2_eff)/(1-j2_eff);

  if(j1_csvBtag <= csvWP && (1-j1_eff) == 0) weight1 = 0.;
  if(j2_csvBtag <= csvWP && (1-j2_eff) == 0) weight2 = 0.;

  return weight1*weight2;
}
