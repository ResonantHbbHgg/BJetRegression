struct tree_variables
{
	// setup tree inputs
	// gen level info
	float gr_radion_p4_pt, gr_radion_p4_eta, gr_radion_p4_phi, gr_radion_p4_mass;
	float gr_hgg_p4_pt, gr_hgg_p4_eta, gr_hgg_p4_phi, gr_hgg_p4_mass;
	float gr_hbb_p4_pt, gr_hbb_p4_eta, gr_hbb_p4_phi, gr_hbb_p4_mass;
	float gr_hjj_p4_pt, gr_hjj_p4_eta, gr_hjj_p4_phi, gr_hjj_p4_mass;
	float gr_g1_p4_pt, gr_g1_p4_eta, gr_g1_p4_phi, gr_g1_p4_mass;
	float gr_g2_p4_pt, gr_g2_p4_eta, gr_g2_p4_phi, gr_g2_p4_mass;
	float gr_b1_p4_pt, gr_b1_p4_eta, gr_b1_p4_phi, gr_b1_p4_mass;
	float gr_b2_p4_pt, gr_b2_p4_eta, gr_b2_p4_phi, gr_b2_p4_mass;
	float gr_j1_p4_pt, gr_j1_p4_eta, gr_j1_p4_phi, gr_j1_p4_mass;
	float gr_j2_p4_pt, gr_j2_p4_eta, gr_j2_p4_phi, gr_j2_p4_mass;
// event variables
	float met_corr_pfmet, met_corr_phi_pfmet, met_corr_eta_pfmet, met_corr_e_pfmet;
	float pu_n, nvtx, rho;
	float weight, evweight, pu_weight;
	float run, lumis, event;
	float ph1_SCEta, ph2_SCEta;
	float ev_weight, ev_evweight, ev_pu_weight;
// object variables
	float ph1_eta, ph2_eta, ph1_pt, ph2_pt, PhotonsMass, ph1_phi, ph2_phi, ph1_e, ph2_e, ph1_r9, ph2_r9, ph1_sieie, ph2_sieie, ph1_hoe, ph2_hoe;
	float ph1_pfchargedisogood03, ph1_ecaliso, ph1_pfchargedisobad04, ph1_ecalisobad, ph1_badvtx_Et, ph1_isconv;
	float ph2_pfchargedisogood03, ph2_ecaliso, ph2_pfchargedisobad04, ph2_ecalisobad, ph2_badvtx_Et, ph2_isconv;
	int ph1_ciclevel, ph2_ciclevel, ph1_isEB, ph2_isEB;

	float j1_e, j1_pt, j1_phi, j1_eta, j1_beta, j1_betaStar, j1_betaStarClassic, j1_dR2Mean, j1_csvBtag, j1_csvMvaBtag, j1_jetProbBtag, j1_tcheBtag, j1_secVtxPt, j1_secVtx3dL, j1_emfrac, j1_hadfrac, j1_axis1, j1_axis2, j1_pull, j1_btagSF_M;
	int j1_ntk, j1_nNeutrals, j1_nCharged;

	float j2_e, j2_pt, j2_phi, j2_eta, j2_beta, j2_betaStar, j2_betaStarClassic, j2_dR2Mean, j2_csvBtag, j2_csvMvaBtag, j2_jetProbBtag, j2_tcheBtag, j2_secVtxPt, j2_secVtx3dL, j2_emfrac, j2_hadfrac, j2_axis1, j2_axis2, j2_pull, j2_btagSF_M;
	int j2_ntk, j2_nNeutrals, j2_nCharged;

	float j3_e, j3_pt, j3_phi, j3_eta, j3_beta, j3_betaStar, j3_betaStarClassic, j3_dR2Mean, j3_csvBtag, j3_csvMvaBtag, j3_jetProbBtag, j3_tcheBtag, j3_secVtxPt, j3_secVtx3dL, j3_emfrac, j3_hadfrac, j3_axis1, j3_axis2, j3_pull, j3_btagSF_M;
	int j3_ntk, j3_nNeutrals, j3_nCharged;

	float j4_e, j4_pt, j4_phi, j4_eta, j4_beta, j4_betaStar, j4_betaStarClassic, j4_dR2Mean, j4_csvBtag, j4_csvMvaBtag, j4_jetProbBtag, j4_tcheBtag, j4_secVtxPt, j4_secVtx3dL, j4_emfrac, j4_hadfrac, j4_axis1, j4_axis2, j4_pull, j4_btagSF_M;
	int j4_ntk, j4_nNeutrals, j4_nCharged;

	float j5_e, j5_pt, j5_phi, j5_eta, j5_beta, j5_betaStar, j5_betaStarClassic, j5_dR2Mean, j5_csvBtag, j5_csvMvaBtag, j5_jetProbBtag, j5_tcheBtag, j5_secVtxPt, j5_secVtx3dL, j5_emfrac, j5_hadfrac, j5_axis1, j5_axis2, j5_pull, j5_btagSF_M;
	int j5_ntk, j5_nNeutrals, j5_nCharged;

	float j6_e, j6_pt, j6_phi, j6_eta, j6_beta, j6_betaStar, j6_betaStarClassic, j6_dR2Mean, j6_csvBtag, j6_csvMvaBtag, j6_jetProbBtag, j6_tcheBtag, j6_secVtxPt, j6_secVtx3dL, j6_emfrac, j6_hadfrac, j6_axis1, j6_axis2, j6_pull, j6_btagSF_M;
	int j6_ntk, j6_nNeutrals, j6_nCharged;

	float j7_e, j7_pt, j7_phi, j7_eta, j7_beta, j7_betaStar, j7_betaStarClassic, j7_dR2Mean, j7_csvBtag, j7_csvMvaBtag, j7_jetProbBtag, j7_tcheBtag, j7_secVtxPt, j7_secVtx3dL, j7_emfrac, j7_hadfrac, j7_axis1, j7_axis2, j7_pull, j7_btagSF_M;
	int j7_ntk, j7_nNeutrals, j7_nCharged;

	float j8_e, j8_pt, j8_phi, j8_eta, j8_beta, j8_betaStar, j8_betaStarClassic, j8_dR2Mean, j8_csvBtag, j8_csvMvaBtag, j8_jetProbBtag, j8_tcheBtag, j8_secVtxPt, j8_secVtx3dL, j8_emfrac, j8_hadfrac, j8_axis1, j8_axis2, j8_pull, j8_btagSF_M;
	int j8_ntk, j8_nNeutrals, j8_nCharged;

	float j9_e, j9_pt, j9_phi, j9_eta, j9_beta, j9_betaStar, j9_betaStarClassic, j9_dR2Mean, j9_csvBtag, j9_csvMvaBtag, j9_jetProbBtag, j9_tcheBtag, j9_secVtxPt, j9_secVtx3dL, j9_emfrac, j9_hadfrac, j9_axis1, j9_axis2, j9_pull, j9_btagSF_M;
	int j9_ntk, j9_nNeutrals, j9_nCharged;

	float j10_e, j10_pt, j10_phi, j10_eta, j10_beta, j10_betaStar, j10_betaStarClassic, j10_dR2Mean, j10_csvBtag, j10_csvMvaBtag, j10_jetProbBtag, j10_tcheBtag, j10_secVtxPt, j10_secVtx3dL, j10_emfrac, j10_hadfrac, j10_axis1, j10_axis2, j10_pull, j10_btagSF_M;
	int j10_ntk, j10_nNeutrals, j10_nCharged;

	float j11_e, j11_pt, j11_phi, j11_eta, j11_beta, j11_betaStar, j11_betaStarClassic, j11_dR2Mean, j11_csvBtag, j11_csvMvaBtag, j11_jetProbBtag, j11_tcheBtag, j11_secVtxPt, j11_secVtx3dL, j11_emfrac, j11_hadfrac, j11_axis1, j11_axis2, j11_pull, j11_btagSF_M;
	int j11_ntk, j11_nNeutrals, j11_nCharged;

	float j12_e, j12_pt, j12_phi, j12_eta, j12_beta, j12_betaStar, j12_betaStarClassic, j12_dR2Mean, j12_csvBtag, j12_csvMvaBtag, j12_jetProbBtag, j12_tcheBtag, j12_secVtxPt, j12_secVtx3dL, j12_emfrac, j12_hadfrac, j12_axis1, j12_axis2, j12_pull, j12_btagSF_M;
	int j12_ntk, j12_nNeutrals, j12_nCharged;

	float j13_e, j13_pt, j13_phi, j13_eta, j13_beta, j13_betaStar, j13_betaStarClassic, j13_dR2Mean, j13_csvBtag, j13_csvMvaBtag, j13_jetProbBtag, j13_tcheBtag, j13_secVtxPt, j13_secVtx3dL, j13_emfrac, j13_hadfrac, j13_axis1, j13_axis2, j13_pull, j13_btagSF_M;
	int j13_ntk, j13_nNeutrals, j13_nCharged;

	float j14_e, j14_pt, j14_phi, j14_eta, j14_beta, j14_betaStar, j14_betaStarClassic, j14_dR2Mean, j14_csvBtag, j14_csvMvaBtag, j14_jetProbBtag, j14_tcheBtag, j14_secVtxPt, j14_secVtx3dL, j14_emfrac, j14_hadfrac, j14_axis1, j14_axis2, j14_pull, j14_btagSF_M;
	int j14_ntk, j14_nNeutrals, j14_nCharged;

	float j15_e, j15_pt, j15_phi, j15_eta, j15_beta, j15_betaStar, j15_betaStarClassic, j15_dR2Mean, j15_csvBtag, j15_csvMvaBtag, j15_jetProbBtag, j15_tcheBtag, j15_secVtxPt, j15_secVtx3dL, j15_emfrac, j15_hadfrac, j15_axis1, j15_axis2, j15_pull, j15_btagSF_M;
	int j15_ntk, j15_nNeutrals, j15_nCharged;

	float jet_e, jet_pt, jet_phi, jet_eta;
	float jet_betaStarClassic, jet_dR2Mean, jet_csvBtag;
	float jet_secVtxPt, jet_secVtx3dL, jet_emfrac, jet_hadfrac, jet_btagSF_M;
	int jet_nNeutrals, jet_nCharged, jet_nConstituents;
	float jet_nConstituents_;
	float jet_dPhiMet_fabs;
	float jet_dPhiMet, jet_regPt, jet_regkinPt;

// setup tree outputs
	float vtx_z;
	float pho1_pt, pho1_e, pho1_phi, pho1_eta, pho1_mass, pho1_r9, pho1_sieie, pho1_hoe;
	float pho2_pt, pho2_e, pho2_phi, pho2_eta, pho2_mass, pho2_r9, pho2_sieie, pho2_hoe;
	int pho1_isEB, pho2_isEB;
	float pho1_pfchargedisogood03, pho1_ecaliso, pho1_pfchargedisobad04, pho1_ecalisobad, pho1_badvtx_Et, pho1_PFisoA, pho1_PFisoB, pho1_PFisoC, pho1_isconv;
	float pho2_pfchargedisogood03, pho2_ecaliso, pho2_pfchargedisobad04, pho2_ecalisobad, pho2_badvtx_Et, pho2_PFisoA, pho2_PFisoB, pho2_PFisoC, pho2_isconv;
	float jet1_pt, jet1_e, jet1_phi, jet1_eta, jet1_mass, jet1_csvBtag, jet1_btagSF_M, jet1_betaStarClassic, jet1_dR2Mean;
	float jet2_pt, jet2_e, jet2_phi, jet2_eta, jet2_mass, jet2_csvBtag, jet2_btagSF_M, jet2_betaStarClassic, jet2_dR2Mean;
	float regjet1_emfrac, regjet1_hadfrac, regjet1_secVtxPt, regjet1_secVtx3dL, regjet1_dPhiMet;
	int regjet1_nConstituents;
	float regjet2_emfrac, regjet2_hadfrac, regjet2_secVtxPt, regjet2_secVtx3dL, regjet2_dPhiMet;
	int regjet2_nConstituents;
	float regjet1_pt, regjet1_e, regjet1_phi, regjet1_eta, regjet1_mass, regjet1_csvBtag, regjet1_btagSF_M, regjet1_betaStarClassic, regjet1_dR2Mean;
	float regjet2_pt, regjet2_e, regjet2_phi, regjet2_eta, regjet2_mass, regjet2_csvBtag, regjet2_btagSF_M, regjet2_betaStarClassic, regjet2_dR2Mean;
	float regkinjet1_pt, regkinjet1_e, regkinjet1_phi, regkinjet1_eta, regkinjet1_mass, regkinjet1_csvBtag, regkinjet1_btagSF_M, regkinjet1_betaStarClassic, regkinjet1_dR2Mean;
	float regkinjet2_pt, regkinjet2_e, regkinjet2_phi, regkinjet2_eta, regkinjet2_mass, regkinjet2_csvBtag, regkinjet2_btagSF_M, regkinjet2_betaStarClassic, regkinjet2_dR2Mean;
	float kinjet1_pt, kinjet1_e, kinjet1_phi, kinjet1_eta, kinjet1_mass, kinjet1_csvBtag, kinjet1_btagSF_M, kinjet1_betaStarClassic, kinjet1_dR2Mean;
	float kinjet2_pt, kinjet2_e, kinjet2_phi, kinjet2_eta, kinjet2_mass, kinjet2_csvBtag, kinjet2_btagSF_M, kinjet2_betaStarClassic, kinjet2_dR2Mean;
	float jj_pt, jj_e, jj_phi, jj_eta, jj_mass, jj_DR, jj_btagSF_M;
	float regjj_pt, regjj_e, regjj_phi, regjj_eta, regjj_mass, regjj_btagSF_M;
	float regkinjj_pt, regkinjj_e, regkinjj_phi, regkinjj_eta, regkinjj_mass, regkinjj_btagSF_M;
	float kinjj_pt, kinjj_e, kinjj_phi, kinjj_eta, kinjj_mass, kinjj_btagSF_M;
	float gg_pt, gg_e, gg_phi, gg_eta, gg_mass;
	float ggjj_pt, ggjj_e, ggjj_phi, ggjj_eta, ggjj_mass, regjj_DR, regkinjj_DR, kinjj_DR;
	float regggjj_pt, regggjj_e, regggjj_phi, regggjj_eta, regggjj_mass;
	float regkinggjj_pt, regkinggjj_e, regkinggjj_phi, regkinggjj_eta, regkinggjj_mass;
	float kinggjj_pt, kinggjj_e, kinggjj_phi, kinggjj_eta, kinggjj_mass;
	int selection_cut_level;
	int category;
	float costhetastar, regcosthetastar, regkincosthetastar, kincosthetastar;
	float minDRgj, minDRgregj, minDRgregkinj, minDRgkinj;
	float HT_gg;

	int njets_passing_kLooseID;
	int njets_passing_kLooseID_and_CSVM;
	int njets_kLooseID, njets_kRadionID;
	int njets_kLooseID_and_CSVM, njets_kRadionID_and_CSVM;
};

void initialize_variables(tree_variables *t)
{
	t->gr_radion_p4_pt = t->gr_radion_p4_eta = t->gr_radion_p4_phi = t->gr_radion_p4_mass = 0.;
	t->gr_hgg_p4_pt = t->gr_hgg_p4_eta = t->gr_hgg_p4_phi = t->gr_hgg_p4_mass = 0.;
	t->gr_hbb_p4_pt = t->gr_hbb_p4_eta = t->gr_hbb_p4_phi = t->gr_hbb_p4_mass = 0.;
	t->gr_hjj_p4_pt = t->gr_hjj_p4_eta = t->gr_hjj_p4_phi = t->gr_hjj_p4_mass = 0.;
	t->gr_g1_p4_pt = t->gr_g1_p4_eta = t->gr_g1_p4_phi = t->gr_g1_p4_mass = 0.;
	t->gr_g2_p4_pt = t->gr_g2_p4_eta = t->gr_g2_p4_phi = t->gr_g2_p4_mass = 0.;
	t->gr_b1_p4_pt = t->gr_b1_p4_eta = t->gr_b1_p4_phi = t->gr_b1_p4_mass = 0.;
	t->gr_b2_p4_pt = t->gr_b2_p4_eta = t->gr_b2_p4_phi = t->gr_b2_p4_mass = 0.;
	t->gr_j1_p4_pt = t->gr_j1_p4_eta = t->gr_j1_p4_phi = t->gr_j1_p4_mass = 0.;
	t->gr_j2_p4_pt = t->gr_j2_p4_eta = t->gr_j2_p4_phi = t->gr_j2_p4_mass = 0.;
	t->selection_cut_level = 0;
	t->category = 0;

	return;
}

void setup_intree(TTree* intree, tree_variables *t, int type)
{
	intree->SetBranchAddress("njets_passing_kLooseID", &t->njets_passing_kLooseID);
	intree->SetBranchAddress("njets_passing_kLooseID_and_CSVM", &t->njets_passing_kLooseID_and_CSVM);
	intree->SetBranchAddress("ph1_eta", &t->ph1_eta);
	intree->SetBranchAddress("ph2_eta", &t->ph2_eta);
	intree->SetBranchAddress("ph1_pt", &t->ph1_pt);
	intree->SetBranchAddress("ph2_pt", &t->ph2_pt);
	intree->SetBranchAddress("ph1_phi", &t->ph1_phi);
	intree->SetBranchAddress("ph2_phi", &t->ph2_phi);
	intree->SetBranchAddress("ph1_e", &t->ph1_e);
	intree->SetBranchAddress("ph2_e", &t->ph2_e);
	intree->SetBranchAddress("ph1_r9", &t->ph1_r9);
	intree->SetBranchAddress("ph2_r9", &t->ph2_r9);
	intree->SetBranchAddress("ph1_sieie", &t->ph1_sieie);
	intree->SetBranchAddress("ph2_sieie", &t->ph2_sieie);
	intree->SetBranchAddress("ph1_hoe", &t->ph1_hoe);
	intree->SetBranchAddress("ph2_hoe", &t->ph2_hoe);
	intree->SetBranchAddress("ph1_isEB", &t->ph1_isEB);
	intree->SetBranchAddress("ph2_isEB", &t->ph2_isEB);
	intree->SetBranchAddress("ph1_pfchargedisogood03", &t->ph1_pfchargedisogood03);
	intree->SetBranchAddress("ph2_pfchargedisogood03", &t->ph2_pfchargedisogood03);
	intree->SetBranchAddress("ph1_ecaliso", &t->ph1_ecaliso);
	intree->SetBranchAddress("ph2_ecaliso", &t->ph2_ecaliso);
	intree->SetBranchAddress("ph1_pfchargedisobad04", &t->ph1_pfchargedisobad04);
	intree->SetBranchAddress("ph2_pfchargedisobad04", &t->ph2_pfchargedisobad04);
	intree->SetBranchAddress("ph1_ecalisobad", &t->ph1_ecalisobad);
	intree->SetBranchAddress("ph2_ecalisobad", &t->ph2_ecalisobad);
	intree->SetBranchAddress("ph1_badvtx_Et", &t->ph1_badvtx_Et);
	intree->SetBranchAddress("ph2_badvtx_Et", &t->ph2_badvtx_Et);
	intree->SetBranchAddress("ph1_isconv", &t->ph1_isconv);
	intree->SetBranchAddress("ph2_isconv", &t->ph2_isconv);
	intree->SetBranchAddress("PhotonsMass", &t->PhotonsMass);
	intree->SetBranchAddress("ph1_ciclevel", &t->ph1_ciclevel);
	intree->SetBranchAddress("ph2_ciclevel", &t->ph2_ciclevel);
	intree->SetBranchAddress("met_corr_pfmet", &t->met_corr_pfmet);
	intree->SetBranchAddress("met_corr_phi_pfmet", &t->met_corr_phi_pfmet);
	intree->SetBranchAddress("met_corr_eta_pfmet", &t->met_corr_eta_pfmet);
	intree->SetBranchAddress("met_corr_e_pfmet", &t->met_corr_e_pfmet);
	intree->SetBranchAddress("pu_n", &t->pu_n);
	intree->SetBranchAddress("nvtx", &t->nvtx);
	intree->SetBranchAddress("rho", &t->rho);
	intree->SetBranchAddress("run", &t->run);
	intree->SetBranchAddress("lumis", &t->lumis);
	intree->SetBranchAddress("event", &t->event);
	intree->SetBranchAddress("ph1_SCEta", &t->ph1_SCEta);
	intree->SetBranchAddress("ph2_SCEta", &t->ph2_SCEta);
	intree->SetBranchAddress("weight", &t->ev_weight);
	intree->SetBranchAddress("evweight", &t->ev_evweight);
	intree->SetBranchAddress("pu_weight", &t->ev_pu_weight);
	intree->SetBranchAddress("vtx_z", &t->vtx_z);

	if( type < -250 )
	{
		intree->SetBranchAddress("gr_radion_p4_pt", &t->gr_radion_p4_pt);
		intree->SetBranchAddress("gr_radion_p4_eta", &t->gr_radion_p4_eta);
		intree->SetBranchAddress("gr_radion_p4_phi", &t->gr_radion_p4_phi);
		intree->SetBranchAddress("gr_radion_p4_mass", &t->gr_radion_p4_mass);
		intree->SetBranchAddress("gr_hgg_p4_pt", &t->gr_hgg_p4_pt);
		intree->SetBranchAddress("gr_hgg_p4_eta", &t->gr_hgg_p4_eta);
		intree->SetBranchAddress("gr_hgg_p4_phi", &t->gr_hgg_p4_phi);
		intree->SetBranchAddress("gr_hgg_p4_mass", &t->gr_hgg_p4_mass);
		intree->SetBranchAddress("gr_hbb_p4_pt", &t->gr_hbb_p4_pt);
		intree->SetBranchAddress("gr_hbb_p4_eta", &t->gr_hbb_p4_eta);
		intree->SetBranchAddress("gr_hbb_p4_phi", &t->gr_hbb_p4_phi);
		intree->SetBranchAddress("gr_hbb_p4_mass", &t->gr_hbb_p4_mass);
		intree->SetBranchAddress("gr_g1_p4_pt", &t->gr_g1_p4_pt);
		intree->SetBranchAddress("gr_g1_p4_eta", &t->gr_g1_p4_eta);
		intree->SetBranchAddress("gr_g1_p4_phi", &t->gr_g1_p4_phi);
		intree->SetBranchAddress("gr_g1_p4_mass", &t->gr_g1_p4_mass);
		intree->SetBranchAddress("gr_g2_p4_pt", &t->gr_g2_p4_pt);
		intree->SetBranchAddress("gr_g2_p4_eta", &t->gr_g2_p4_eta);
		intree->SetBranchAddress("gr_g2_p4_phi", &t->gr_g2_p4_phi);
		intree->SetBranchAddress("gr_g2_p4_mass", &t->gr_g2_p4_mass);
		intree->SetBranchAddress("gr_b1_p4_pt", &t->gr_b1_p4_pt);
		intree->SetBranchAddress("gr_b1_p4_eta", &t->gr_b1_p4_eta);
		intree->SetBranchAddress("gr_b1_p4_phi", &t->gr_b1_p4_phi);
		intree->SetBranchAddress("gr_b1_p4_mass", &t->gr_b1_p4_mass);
		intree->SetBranchAddress("gr_b2_p4_pt", &t->gr_b2_p4_pt);
		intree->SetBranchAddress("gr_b2_p4_eta", &t->gr_b2_p4_eta);
		intree->SetBranchAddress("gr_b2_p4_phi", &t->gr_b2_p4_phi);
		intree->SetBranchAddress("gr_b2_p4_mass", &t->gr_b2_p4_mass);
		intree->SetBranchAddress("gr_j1_p4_pt", &t->gr_j1_p4_pt);
		intree->SetBranchAddress("gr_j1_p4_eta", &t->gr_j1_p4_eta);
		intree->SetBranchAddress("gr_j1_p4_phi", &t->gr_j1_p4_phi);
		intree->SetBranchAddress("gr_j1_p4_mass", &t->gr_j1_p4_mass);
		intree->SetBranchAddress("gr_j2_p4_pt", &t->gr_j2_p4_pt);
		intree->SetBranchAddress("gr_j2_p4_eta", &t->gr_j2_p4_eta);
		intree->SetBranchAddress("gr_j2_p4_phi", &t->gr_j2_p4_phi);
		intree->SetBranchAddress("gr_j2_p4_mass", &t->gr_j2_p4_mass);
	}
	intree->SetBranchAddress("j1_e", &t->j1_e);
	intree->SetBranchAddress("j1_pt", &t->j1_pt);
	intree->SetBranchAddress("j1_phi", &t->j1_phi);
	intree->SetBranchAddress("j1_eta", &t->j1_eta);
	intree->SetBranchAddress("j1_beta", &t->j1_beta);
	intree->SetBranchAddress("j1_betaStar", &t->j1_betaStar);
	intree->SetBranchAddress("j1_betaStarClassic", &t->j1_betaStarClassic);
	intree->SetBranchAddress("j1_dR2Mean", &t->j1_dR2Mean);
	intree->SetBranchAddress("j1_csvBtag", &t->j1_csvBtag);
	intree->SetBranchAddress("j1_csvMvaBtag", &t->j1_csvMvaBtag);
	intree->SetBranchAddress("j1_jetProbBtag", &t->j1_jetProbBtag);
	intree->SetBranchAddress("j1_tcheBtag", &t->j1_tcheBtag);
	intree->SetBranchAddress("j1_secVtxPt", &t->j1_secVtxPt);
	intree->SetBranchAddress("j1_secVtx3dL", &t->j1_secVtx3dL);
	intree->SetBranchAddress("j1_emfrac", &t->j1_emfrac);
	intree->SetBranchAddress("j1_hadfrac", &t->j1_hadfrac);
	intree->SetBranchAddress("j1_btagSF_M", &t->j1_btagSF_M);
	intree->SetBranchAddress("j1_ntk", &t->j1_ntk);
	intree->SetBranchAddress("j1_nNeutrals", &t->j1_nNeutrals);
	intree->SetBranchAddress("j1_nCharged", &t->j1_nCharged);
	intree->SetBranchAddress("j1_axis1", &t->j1_axis1);
	intree->SetBranchAddress("j1_axis2", &t->j1_axis2);
	intree->SetBranchAddress("j1_pull", &t->j1_pull);

	intree->SetBranchAddress("j2_e", &t->j2_e);
	intree->SetBranchAddress("j2_pt", &t->j2_pt);
	intree->SetBranchAddress("j2_phi", &t->j2_phi);
	intree->SetBranchAddress("j2_eta", &t->j2_eta);
	intree->SetBranchAddress("j2_beta", &t->j2_beta);
	intree->SetBranchAddress("j2_betaStar", &t->j2_betaStar);
	intree->SetBranchAddress("j2_betaStarClassic", &t->j2_betaStarClassic);
	intree->SetBranchAddress("j2_dR2Mean", &t->j2_dR2Mean);
	intree->SetBranchAddress("j2_csvBtag", &t->j2_csvBtag);
	intree->SetBranchAddress("j2_csvMvaBtag", &t->j2_csvMvaBtag);
	intree->SetBranchAddress("j2_jetProbBtag", &t->j2_jetProbBtag);
	intree->SetBranchAddress("j2_tcheBtag", &t->j2_tcheBtag);
	intree->SetBranchAddress("j2_secVtxPt", &t->j2_secVtxPt);
	intree->SetBranchAddress("j2_secVtx3dL", &t->j2_secVtx3dL);
	intree->SetBranchAddress("j2_emfrac", &t->j2_emfrac);
	intree->SetBranchAddress("j2_hadfrac", &t->j2_hadfrac);
	intree->SetBranchAddress("j2_btagSF_M", &t->j2_btagSF_M);
	intree->SetBranchAddress("j2_ntk", &t->j2_ntk);
	intree->SetBranchAddress("j2_nNeutrals", &t->j2_nNeutrals);
	intree->SetBranchAddress("j2_nCharged", &t->j2_nCharged);
	intree->SetBranchAddress("j2_axis1", &t->j2_axis1);
	intree->SetBranchAddress("j2_axis2", &t->j2_axis2);
	intree->SetBranchAddress("j2_pull", &t->j2_pull);

	intree->SetBranchAddress("j3_e", &t->j3_e);
	intree->SetBranchAddress("j3_pt", &t->j3_pt);
	intree->SetBranchAddress("j3_phi", &t->j3_phi);
	intree->SetBranchAddress("j3_eta", &t->j3_eta);
	intree->SetBranchAddress("j3_beta", &t->j3_beta);
	intree->SetBranchAddress("j3_betaStar", &t->j3_betaStar);
	intree->SetBranchAddress("j3_betaStarClassic", &t->j3_betaStarClassic);
	intree->SetBranchAddress("j3_dR2Mean", &t->j3_dR2Mean);
	intree->SetBranchAddress("j3_csvBtag", &t->j3_csvBtag);
	intree->SetBranchAddress("j3_csvMvaBtag", &t->j3_csvMvaBtag);
	intree->SetBranchAddress("j3_jetProbBtag", &t->j3_jetProbBtag);
	intree->SetBranchAddress("j3_tcheBtag", &t->j3_tcheBtag);
	intree->SetBranchAddress("j3_secVtxPt", &t->j3_secVtxPt);
	intree->SetBranchAddress("j3_secVtx3dL", &t->j3_secVtx3dL);
	intree->SetBranchAddress("j3_emfrac", &t->j3_emfrac);
	intree->SetBranchAddress("j3_hadfrac", &t->j3_hadfrac);
	intree->SetBranchAddress("j3_btagSF_M", &t->j3_btagSF_M);
	intree->SetBranchAddress("j3_ntk", &t->j3_ntk);
	intree->SetBranchAddress("j3_nNeutrals", &t->j3_nNeutrals);
	intree->SetBranchAddress("j3_nCharged", &t->j3_nCharged);
	intree->SetBranchAddress("j3_axis1", &t->j3_axis1);
	intree->SetBranchAddress("j3_axis2", &t->j3_axis2);
	intree->SetBranchAddress("j3_pull", &t->j3_pull);

	intree->SetBranchAddress("j4_e", &t->j4_e);
	intree->SetBranchAddress("j4_pt", &t->j4_pt);
	intree->SetBranchAddress("j4_phi", &t->j4_phi);
	intree->SetBranchAddress("j4_eta", &t->j4_eta);
	intree->SetBranchAddress("j4_beta", &t->j4_beta);
	intree->SetBranchAddress("j4_betaStar", &t->j4_betaStar);
	intree->SetBranchAddress("j4_betaStarClassic", &t->j4_betaStarClassic);
	intree->SetBranchAddress("j4_dR2Mean", &t->j4_dR2Mean);
	intree->SetBranchAddress("j4_csvBtag", &t->j4_csvBtag);
	intree->SetBranchAddress("j4_csvMvaBtag", &t->j4_csvMvaBtag);
	intree->SetBranchAddress("j4_jetProbBtag", &t->j4_jetProbBtag);
	intree->SetBranchAddress("j4_tcheBtag", &t->j4_tcheBtag);
	intree->SetBranchAddress("j4_secVtxPt", &t->j4_secVtxPt);
	intree->SetBranchAddress("j4_secVtx3dL", &t->j4_secVtx3dL);
	intree->SetBranchAddress("j4_emfrac", &t->j4_emfrac);
	intree->SetBranchAddress("j4_hadfrac", &t->j4_hadfrac);
	intree->SetBranchAddress("j4_btagSF_M", &t->j4_btagSF_M);
	intree->SetBranchAddress("j4_ntk", &t->j4_ntk);
	intree->SetBranchAddress("j4_nNeutrals", &t->j4_nNeutrals);
	intree->SetBranchAddress("j4_nCharged", &t->j4_nCharged);
	intree->SetBranchAddress("j4_axis1", &t->j4_axis1);
	intree->SetBranchAddress("j4_axis2", &t->j4_axis2);
	intree->SetBranchAddress("j4_pull", &t->j4_pull);

	intree->SetBranchAddress("j5_e", &t->j5_e);
	intree->SetBranchAddress("j5_pt", &t->j5_pt);
	intree->SetBranchAddress("j5_phi", &t->j5_phi);
	intree->SetBranchAddress("j5_eta", &t->j5_eta);
	intree->SetBranchAddress("j5_beta", &t->j5_beta);
	intree->SetBranchAddress("j5_betaStar", &t->j5_betaStar);
	intree->SetBranchAddress("j5_betaStarClassic", &t->j5_betaStarClassic);
	intree->SetBranchAddress("j5_dR2Mean", &t->j5_dR2Mean);
	intree->SetBranchAddress("j5_csvBtag", &t->j5_csvBtag);
	intree->SetBranchAddress("j5_csvMvaBtag", &t->j5_csvMvaBtag);
	intree->SetBranchAddress("j5_jetProbBtag", &t->j5_jetProbBtag);
	intree->SetBranchAddress("j5_tcheBtag", &t->j5_tcheBtag);
	intree->SetBranchAddress("j5_secVtxPt", &t->j5_secVtxPt);
	intree->SetBranchAddress("j5_secVtx3dL", &t->j5_secVtx3dL);
	intree->SetBranchAddress("j5_emfrac", &t->j5_emfrac);
	intree->SetBranchAddress("j5_hadfrac", &t->j5_hadfrac);
	intree->SetBranchAddress("j5_btagSF_M", &t->j5_btagSF_M);
	intree->SetBranchAddress("j5_ntk", &t->j5_ntk);
	intree->SetBranchAddress("j5_nNeutrals", &t->j5_nNeutrals);
	intree->SetBranchAddress("j5_nCharged", &t->j5_nCharged);
	intree->SetBranchAddress("j5_axis1", &t->j5_axis1);
	intree->SetBranchAddress("j5_axis2", &t->j5_axis2);
	intree->SetBranchAddress("j5_pull", &t->j5_pull);

	intree->SetBranchAddress("j6_e", &t->j6_e);
	intree->SetBranchAddress("j6_pt", &t->j6_pt);
	intree->SetBranchAddress("j6_phi", &t->j6_phi);
	intree->SetBranchAddress("j6_eta", &t->j6_eta);
	intree->SetBranchAddress("j6_beta", &t->j6_beta);
	intree->SetBranchAddress("j6_betaStar", &t->j6_betaStar);
	intree->SetBranchAddress("j6_betaStarClassic", &t->j6_betaStarClassic);
	intree->SetBranchAddress("j6_dR2Mean", &t->j6_dR2Mean);
	intree->SetBranchAddress("j6_csvBtag", &t->j6_csvBtag);
	intree->SetBranchAddress("j6_csvMvaBtag", &t->j6_csvMvaBtag);
	intree->SetBranchAddress("j6_jetProbBtag", &t->j6_jetProbBtag);
	intree->SetBranchAddress("j6_tcheBtag", &t->j6_tcheBtag);
	intree->SetBranchAddress("j6_secVtxPt", &t->j6_secVtxPt);
	intree->SetBranchAddress("j6_secVtx3dL", &t->j6_secVtx3dL);
	intree->SetBranchAddress("j6_emfrac", &t->j6_emfrac);
	intree->SetBranchAddress("j6_hadfrac", &t->j6_hadfrac);
	intree->SetBranchAddress("j6_btagSF_M", &t->j6_btagSF_M);
	intree->SetBranchAddress("j6_ntk", &t->j6_ntk);
	intree->SetBranchAddress("j6_nNeutrals", &t->j6_nNeutrals);
	intree->SetBranchAddress("j6_nCharged", &t->j6_nCharged);
	intree->SetBranchAddress("j6_axis1", &t->j6_axis1);
	intree->SetBranchAddress("j6_axis2", &t->j6_axis2);
	intree->SetBranchAddress("j6_pull", &t->j6_pull);

	intree->SetBranchAddress("j7_e", &t->j7_e);
	intree->SetBranchAddress("j7_pt", &t->j7_pt);
	intree->SetBranchAddress("j7_phi", &t->j7_phi);
	intree->SetBranchAddress("j7_eta", &t->j7_eta);
	intree->SetBranchAddress("j7_beta", &t->j7_beta);
	intree->SetBranchAddress("j7_betaStar", &t->j7_betaStar);
	intree->SetBranchAddress("j7_betaStarClassic", &t->j7_betaStarClassic);
	intree->SetBranchAddress("j7_dR2Mean", &t->j7_dR2Mean);
	intree->SetBranchAddress("j7_csvBtag", &t->j7_csvBtag);
	intree->SetBranchAddress("j7_csvMvaBtag", &t->j7_csvMvaBtag);
	intree->SetBranchAddress("j7_jetProbBtag", &t->j7_jetProbBtag);
	intree->SetBranchAddress("j7_tcheBtag", &t->j7_tcheBtag);
	intree->SetBranchAddress("j7_secVtxPt", &t->j7_secVtxPt);
	intree->SetBranchAddress("j7_secVtx3dL", &t->j7_secVtx3dL);
	intree->SetBranchAddress("j7_emfrac", &t->j7_emfrac);
	intree->SetBranchAddress("j7_hadfrac", &t->j7_hadfrac);
	intree->SetBranchAddress("j7_btagSF_M", &t->j7_btagSF_M);
	intree->SetBranchAddress("j7_ntk", &t->j7_ntk);
	intree->SetBranchAddress("j7_nNeutrals", &t->j7_nNeutrals);
	intree->SetBranchAddress("j7_nCharged", &t->j7_nCharged);
	intree->SetBranchAddress("j7_axis1", &t->j7_axis1);
	intree->SetBranchAddress("j7_axis2", &t->j7_axis2);
	intree->SetBranchAddress("j7_pull", &t->j7_pull);

	intree->SetBranchAddress("j8_e", &t->j8_e);
	intree->SetBranchAddress("j8_pt", &t->j8_pt);
	intree->SetBranchAddress("j8_phi", &t->j8_phi);
	intree->SetBranchAddress("j8_eta", &t->j8_eta);
	intree->SetBranchAddress("j8_beta", &t->j8_beta);
	intree->SetBranchAddress("j8_betaStar", &t->j8_betaStar);
	intree->SetBranchAddress("j8_betaStarClassic", &t->j8_betaStarClassic);
	intree->SetBranchAddress("j8_dR2Mean", &t->j8_dR2Mean);
	intree->SetBranchAddress("j8_csvBtag", &t->j8_csvBtag);
	intree->SetBranchAddress("j8_csvMvaBtag", &t->j8_csvMvaBtag);
	intree->SetBranchAddress("j8_jetProbBtag", &t->j8_jetProbBtag);
	intree->SetBranchAddress("j8_tcheBtag", &t->j8_tcheBtag);
	intree->SetBranchAddress("j8_secVtxPt", &t->j8_secVtxPt);
	intree->SetBranchAddress("j8_secVtx3dL", &t->j8_secVtx3dL);
	intree->SetBranchAddress("j8_emfrac", &t->j8_emfrac);
	intree->SetBranchAddress("j8_hadfrac", &t->j8_hadfrac);
	intree->SetBranchAddress("j8_btagSF_M", &t->j8_btagSF_M);
	intree->SetBranchAddress("j8_ntk", &t->j8_ntk);
	intree->SetBranchAddress("j8_nNeutrals", &t->j8_nNeutrals);
	intree->SetBranchAddress("j8_nCharged", &t->j8_nCharged);
	intree->SetBranchAddress("j8_axis1", &t->j8_axis1);
	intree->SetBranchAddress("j8_axis2", &t->j8_axis2);
	intree->SetBranchAddress("j8_pull", &t->j8_pull);

	intree->SetBranchAddress("j9_e", &t->j9_e);
	intree->SetBranchAddress("j9_pt", &t->j9_pt);
	intree->SetBranchAddress("j9_phi", &t->j9_phi);
	intree->SetBranchAddress("j9_eta", &t->j9_eta);
	intree->SetBranchAddress("j9_beta", &t->j9_beta);
	intree->SetBranchAddress("j9_betaStar", &t->j9_betaStar);
	intree->SetBranchAddress("j9_betaStarClassic", &t->j9_betaStarClassic);
	intree->SetBranchAddress("j9_dR2Mean", &t->j9_dR2Mean);
	intree->SetBranchAddress("j9_csvBtag", &t->j9_csvBtag);
	intree->SetBranchAddress("j9_csvMvaBtag", &t->j9_csvMvaBtag);
	intree->SetBranchAddress("j9_jetProbBtag", &t->j9_jetProbBtag);
	intree->SetBranchAddress("j9_tcheBtag", &t->j9_tcheBtag);
	intree->SetBranchAddress("j9_secVtxPt", &t->j9_secVtxPt);
	intree->SetBranchAddress("j9_secVtx3dL", &t->j9_secVtx3dL);
	intree->SetBranchAddress("j9_emfrac", &t->j9_emfrac);
	intree->SetBranchAddress("j9_hadfrac", &t->j9_hadfrac);
	intree->SetBranchAddress("j9_btagSF_M", &t->j9_btagSF_M);
	intree->SetBranchAddress("j9_ntk", &t->j9_ntk);
	intree->SetBranchAddress("j9_nNeutrals", &t->j9_nNeutrals);
	intree->SetBranchAddress("j9_nCharged", &t->j9_nCharged);
	intree->SetBranchAddress("j9_axis1", &t->j9_axis1);
	intree->SetBranchAddress("j9_axis2", &t->j9_axis2);
	intree->SetBranchAddress("j9_pull", &t->j9_pull);

	intree->SetBranchAddress("j10_e", &t->j10_e);
	intree->SetBranchAddress("j10_pt", &t->j10_pt);
	intree->SetBranchAddress("j10_phi", &t->j10_phi);
	intree->SetBranchAddress("j10_eta", &t->j10_eta);
	intree->SetBranchAddress("j10_beta", &t->j10_beta);
	intree->SetBranchAddress("j10_betaStar", &t->j10_betaStar);
	intree->SetBranchAddress("j10_betaStarClassic", &t->j10_betaStarClassic);
	intree->SetBranchAddress("j10_dR2Mean", &t->j10_dR2Mean);
	intree->SetBranchAddress("j10_csvBtag", &t->j10_csvBtag);
	intree->SetBranchAddress("j10_csvMvaBtag", &t->j10_csvMvaBtag);
	intree->SetBranchAddress("j10_jetProbBtag", &t->j10_jetProbBtag);
	intree->SetBranchAddress("j10_tcheBtag", &t->j10_tcheBtag);
	intree->SetBranchAddress("j10_secVtxPt", &t->j10_secVtxPt);
	intree->SetBranchAddress("j10_secVtx3dL", &t->j10_secVtx3dL);
	intree->SetBranchAddress("j10_emfrac", &t->j10_emfrac);
	intree->SetBranchAddress("j10_hadfrac", &t->j10_hadfrac);
	intree->SetBranchAddress("j10_btagSF_M", &t->j10_btagSF_M);
	intree->SetBranchAddress("j10_ntk", &t->j10_ntk);
	intree->SetBranchAddress("j10_nNeutrals", &t->j10_nNeutrals);
	intree->SetBranchAddress("j10_nCharged", &t->j10_nCharged);
	intree->SetBranchAddress("j10_axis1", &t->j10_axis1);
	intree->SetBranchAddress("j10_axis2", &t->j10_axis2);
	intree->SetBranchAddress("j10_pull", &t->j10_pull);

	intree->SetBranchAddress("j11_e", &t->j11_e);
	intree->SetBranchAddress("j11_pt", &t->j11_pt);
	intree->SetBranchAddress("j11_phi", &t->j11_phi);
	intree->SetBranchAddress("j11_eta", &t->j11_eta);
	intree->SetBranchAddress("j11_beta", &t->j11_beta);
	intree->SetBranchAddress("j11_betaStar", &t->j11_betaStar);
	intree->SetBranchAddress("j11_betaStarClassic", &t->j11_betaStarClassic);
	intree->SetBranchAddress("j11_dR2Mean", &t->j11_dR2Mean);
	intree->SetBranchAddress("j11_csvBtag", &t->j11_csvBtag);
	intree->SetBranchAddress("j11_csvMvaBtag", &t->j11_csvMvaBtag);
	intree->SetBranchAddress("j11_jetProbBtag", &t->j11_jetProbBtag);
	intree->SetBranchAddress("j11_tcheBtag", &t->j11_tcheBtag);
	intree->SetBranchAddress("j11_secVtxPt", &t->j11_secVtxPt);
	intree->SetBranchAddress("j11_secVtx3dL", &t->j11_secVtx3dL);
	intree->SetBranchAddress("j11_emfrac", &t->j11_emfrac);
	intree->SetBranchAddress("j11_hadfrac", &t->j11_hadfrac);
	intree->SetBranchAddress("j11_btagSF_M", &t->j11_btagSF_M);
	intree->SetBranchAddress("j11_ntk", &t->j11_ntk);
	intree->SetBranchAddress("j11_nNeutrals", &t->j11_nNeutrals);
	intree->SetBranchAddress("j11_nCharged", &t->j11_nCharged);
	intree->SetBranchAddress("j11_axis1", &t->j11_axis1);
	intree->SetBranchAddress("j11_axis2", &t->j11_axis2);
	intree->SetBranchAddress("j11_pull", &t->j11_pull);

	intree->SetBranchAddress("j12_e", &t->j12_e);
	intree->SetBranchAddress("j12_pt", &t->j12_pt);
	intree->SetBranchAddress("j12_phi", &t->j12_phi);
	intree->SetBranchAddress("j12_eta", &t->j12_eta);
	intree->SetBranchAddress("j12_beta", &t->j12_beta);
	intree->SetBranchAddress("j12_betaStar", &t->j12_betaStar);
	intree->SetBranchAddress("j12_betaStarClassic", &t->j12_betaStarClassic);
	intree->SetBranchAddress("j12_dR2Mean", &t->j12_dR2Mean);
	intree->SetBranchAddress("j12_csvBtag", &t->j12_csvBtag);
	intree->SetBranchAddress("j12_csvMvaBtag", &t->j12_csvMvaBtag);
	intree->SetBranchAddress("j12_jetProbBtag", &t->j12_jetProbBtag);
	intree->SetBranchAddress("j12_tcheBtag", &t->j12_tcheBtag);
	intree->SetBranchAddress("j12_secVtxPt", &t->j12_secVtxPt);
	intree->SetBranchAddress("j12_secVtx3dL", &t->j12_secVtx3dL);
	intree->SetBranchAddress("j12_emfrac", &t->j12_emfrac);
	intree->SetBranchAddress("j12_hadfrac", &t->j12_hadfrac);
	intree->SetBranchAddress("j12_btagSF_M", &t->j12_btagSF_M);
	intree->SetBranchAddress("j12_ntk", &t->j12_ntk);
	intree->SetBranchAddress("j12_nNeutrals", &t->j12_nNeutrals);
	intree->SetBranchAddress("j12_nCharged", &t->j12_nCharged);
	intree->SetBranchAddress("j12_axis1", &t->j12_axis1);
	intree->SetBranchAddress("j12_axis2", &t->j12_axis2);
	intree->SetBranchAddress("j12_pull", &t->j12_pull);

	intree->SetBranchAddress("j13_e", &t->j13_e);
	intree->SetBranchAddress("j13_pt", &t->j13_pt);
	intree->SetBranchAddress("j13_phi", &t->j13_phi);
	intree->SetBranchAddress("j13_eta", &t->j13_eta);
	intree->SetBranchAddress("j13_beta", &t->j13_beta);
	intree->SetBranchAddress("j13_betaStar", &t->j13_betaStar);
	intree->SetBranchAddress("j13_betaStarClassic", &t->j13_betaStarClassic);
	intree->SetBranchAddress("j13_dR2Mean", &t->j13_dR2Mean);
	intree->SetBranchAddress("j13_csvBtag", &t->j13_csvBtag);
	intree->SetBranchAddress("j13_csvMvaBtag", &t->j13_csvMvaBtag);
	intree->SetBranchAddress("j13_jetProbBtag", &t->j13_jetProbBtag);
	intree->SetBranchAddress("j13_tcheBtag", &t->j13_tcheBtag);
	intree->SetBranchAddress("j13_secVtxPt", &t->j13_secVtxPt);
	intree->SetBranchAddress("j13_secVtx3dL", &t->j13_secVtx3dL);
	intree->SetBranchAddress("j13_emfrac", &t->j13_emfrac);
	intree->SetBranchAddress("j13_hadfrac", &t->j13_hadfrac);
	intree->SetBranchAddress("j13_btagSF_M", &t->j13_btagSF_M);
	intree->SetBranchAddress("j13_ntk", &t->j13_ntk);
	intree->SetBranchAddress("j13_nNeutrals", &t->j13_nNeutrals);
	intree->SetBranchAddress("j13_nCharged", &t->j13_nCharged);
	intree->SetBranchAddress("j13_axis1", &t->j13_axis1);
	intree->SetBranchAddress("j13_axis2", &t->j13_axis2);
	intree->SetBranchAddress("j13_pull", &t->j13_pull);

	intree->SetBranchAddress("j14_e", &t->j14_e);
	intree->SetBranchAddress("j14_pt", &t->j14_pt);
	intree->SetBranchAddress("j14_phi", &t->j14_phi);
	intree->SetBranchAddress("j14_eta", &t->j14_eta);
	intree->SetBranchAddress("j14_beta", &t->j14_beta);
	intree->SetBranchAddress("j14_betaStar", &t->j14_betaStar);
	intree->SetBranchAddress("j14_betaStarClassic", &t->j14_betaStarClassic);
	intree->SetBranchAddress("j14_dR2Mean", &t->j14_dR2Mean);
	intree->SetBranchAddress("j14_csvBtag", &t->j14_csvBtag);
	intree->SetBranchAddress("j14_csvMvaBtag", &t->j14_csvMvaBtag);
	intree->SetBranchAddress("j14_jetProbBtag", &t->j14_jetProbBtag);
	intree->SetBranchAddress("j14_tcheBtag", &t->j14_tcheBtag);
	intree->SetBranchAddress("j14_secVtxPt", &t->j14_secVtxPt);
	intree->SetBranchAddress("j14_secVtx3dL", &t->j14_secVtx3dL);
	intree->SetBranchAddress("j14_emfrac", &t->j14_emfrac);
	intree->SetBranchAddress("j14_hadfrac", &t->j14_hadfrac);
	intree->SetBranchAddress("j14_btagSF_M", &t->j14_btagSF_M);
	intree->SetBranchAddress("j14_ntk", &t->j14_ntk);
	intree->SetBranchAddress("j14_nNeutrals", &t->j14_nNeutrals);
	intree->SetBranchAddress("j14_nCharged", &t->j14_nCharged);
	intree->SetBranchAddress("j14_axis1", &t->j14_axis1);
	intree->SetBranchAddress("j14_axis2", &t->j14_axis2);
	intree->SetBranchAddress("j14_pull", &t->j14_pull);

	intree->SetBranchAddress("j15_e", &t->j15_e);
	intree->SetBranchAddress("j15_pt", &t->j15_pt);
	intree->SetBranchAddress("j15_phi", &t->j15_phi);
	intree->SetBranchAddress("j15_eta", &t->j15_eta);
	intree->SetBranchAddress("j15_beta", &t->j15_beta);
	intree->SetBranchAddress("j15_betaStar", &t->j15_betaStar);
	intree->SetBranchAddress("j15_betaStarClassic", &t->j15_betaStarClassic);
	intree->SetBranchAddress("j15_dR2Mean", &t->j15_dR2Mean);
	intree->SetBranchAddress("j15_csvBtag", &t->j15_csvBtag);
	intree->SetBranchAddress("j15_csvMvaBtag", &t->j15_csvMvaBtag);
	intree->SetBranchAddress("j15_jetProbBtag", &t->j15_jetProbBtag);
	intree->SetBranchAddress("j15_tcheBtag", &t->j15_tcheBtag);
	intree->SetBranchAddress("j15_secVtxPt", &t->j15_secVtxPt);
	intree->SetBranchAddress("j15_secVtx3dL", &t->j15_secVtx3dL);
	intree->SetBranchAddress("j15_emfrac", &t->j15_emfrac);
	intree->SetBranchAddress("j15_hadfrac", &t->j15_hadfrac);
	intree->SetBranchAddress("j15_btagSF_M", &t->j15_btagSF_M);
	intree->SetBranchAddress("j15_ntk", &t->j15_ntk);
	intree->SetBranchAddress("j15_nNeutrals", &t->j15_nNeutrals);
	intree->SetBranchAddress("j15_nCharged", &t->j15_nCharged);
	intree->SetBranchAddress("j15_axis1", &t->j15_axis1);
	intree->SetBranchAddress("j15_axis2", &t->j15_axis2);
	intree->SetBranchAddress("j15_pull", &t->j15_pull);

	return;
}

void setup_outtree(TTree* outtree, tree_variables *t)
{
	outtree->Branch("category", &t->category, "category/I");
	outtree->Branch("selection_cut_level", &t->selection_cut_level, "selection_cut_level/I");
	outtree->Branch("event", &t->event, "event/F");
	outtree->Branch("vtx_z", &t->vtx_z, "vtx_z/F");
	outtree->Branch("weight", &t->weight, "weight/F");
	outtree->Branch("evweight", &t->evweight, "evweight/F");
	outtree->Branch("pu_weight", &t->pu_weight, "pu_weight/F");
	outtree->Branch("nvtx", &t->nvtx, "nvtx/F");
	outtree->Branch("rho", &t->rho, "rho/F");
	outtree->Branch("met_corr_pfmet", &t->met_corr_pfmet, "met_corr_pfmet/F");
	outtree->Branch("met_corr_phi_pfmet", &t->met_corr_phi_pfmet, "met_corr_phi_pfmet/F");
	outtree->Branch("met_corr_eta_pfmet", &t->met_corr_eta_pfmet, "met_corr_eta_pfmet/F");
	outtree->Branch("met_corr_e_pfmet", &t->met_corr_e_pfmet, "met_corr_e_pfmet/F");
	outtree->Branch("pho1_pt", &t->pho1_pt, "pho1_pt/F");
	outtree->Branch("pho1_e", &t->pho1_e, "pho1_e/F");
	outtree->Branch("pho1_phi", &t->pho1_phi, "pho1_phi/F");
	outtree->Branch("pho1_eta", &t->pho1_eta, "pho1_eta/F");
	outtree->Branch("pho1_mass", &t->pho1_mass, "pho1_mass/F");
	outtree->Branch("pho1_r9", &t->pho1_r9, "pho1_r9/F");
	outtree->Branch("pho1_sieie", &t->pho1_sieie, "pho1_sieie/F");
	outtree->Branch("pho1_hoe", &t->pho1_hoe, "pho1_hoe/F");
	outtree->Branch("pho1_isEB", &t->pho1_isEB, "pho1_isEB/I");
	outtree->Branch("pho1_pfchargedisogood03", &t->pho1_pfchargedisogood03, "pho1_pfchargedisogood03/F");
	outtree->Branch("pho1_ecaliso", &t->pho1_ecaliso, "pho1_ecaliso/F");
	outtree->Branch("pho1_pfchargedisobad04", &t->pho1_pfchargedisobad04, "pho1_pfchargedisobad04/F");
	outtree->Branch("pho1_ecalisobad", &t->pho1_ecalisobad, "pho1_ecalisobad/F");
	outtree->Branch("pho1_badvtx_Et", &t->pho1_badvtx_Et, "pho1_badvtx_Et/F");
	outtree->Branch("pho1_isconv", &t->pho1_isconv, "pho1_isconv/F");
	outtree->Branch("pho1_PFisoA", &t->pho1_PFisoA, "pho1_PFisoA/F");
	outtree->Branch("pho1_PFisoB", &t->pho1_PFisoB, "pho1_PFisoB/F");
	outtree->Branch("pho1_PFisoC", &t->pho1_PFisoC, "pho1_PFisoC/F");
	outtree->Branch("pho2_pt", &t->pho2_pt, "pho2_pt/F");
	outtree->Branch("pho2_e", &t->pho2_e, "pho2_e/F");
	outtree->Branch("pho2_phi", &t->pho2_phi, "pho2_phi/F");
	outtree->Branch("pho2_eta", &t->pho2_eta, "pho2_eta/F");
	outtree->Branch("pho2_mass", &t->pho2_mass, "pho2_mass/F");
	outtree->Branch("pho2_r9", &t->pho2_r9, "pho2_r9/F");
	outtree->Branch("pho2_sieie", &t->pho2_sieie, "pho2_sieie/F");
	outtree->Branch("pho2_hoe", &t->pho2_hoe, "pho2_hoe/F");
	outtree->Branch("pho2_isEB", &t->pho2_isEB, "pho2_isEB/I");
	outtree->Branch("pho2_pfchargedisogood03", &t->pho2_pfchargedisogood03, "pho2_pfchargedisogood03/F");
	outtree->Branch("pho2_ecaliso", &t->pho2_ecaliso, "pho2_ecaliso/F");
	outtree->Branch("pho2_pfchargedisobad04", &t->pho2_pfchargedisobad04, "pho2_pfchargedisobad04/F");
	outtree->Branch("pho2_ecalisobad", &t->pho2_ecalisobad, "pho2_ecalisobad/F");
	outtree->Branch("pho2_badvtx_Et", &t->pho2_badvtx_Et, "pho2_badvtx_Et/F");
	outtree->Branch("pho2_isconv", &t->pho2_isconv, "pho2_isconv/F");
	outtree->Branch("pho2_PFisoA", &t->pho2_PFisoA, "pho2_PFisoA/F");
	outtree->Branch("pho2_PFisoB", &t->pho2_PFisoB, "pho2_PFisoB/F");
	outtree->Branch("pho2_PFisoC", &t->pho2_PFisoC, "pho2_PFisoC/F");
	outtree->Branch("jet1_pt", &t->jet1_pt, "jet1_pt/F");
	outtree->Branch("jet1_pt", &t->jet1_pt, "jet1_pt/F");
	outtree->Branch("jet1_e", &t->jet1_e, "jet1_e/F");
	outtree->Branch("jet1_phi", &t->jet1_phi, "jet1_phi/F");
	outtree->Branch("jet1_eta", &t->jet1_eta, "jet1_eta/F");
	outtree->Branch("jet1_mass", &t->jet1_mass, "jet1_mass/F");
	outtree->Branch("jet1_csvBtag", &t->jet1_csvBtag, "jet1_csvBtag/F");
	outtree->Branch("jet1_btagSF_M", &t->jet1_btagSF_M, "jet1_btagSF_M/F");
	outtree->Branch("jet1_betaStarClassic", &t->jet1_betaStarClassic, "jet1_betaStarClassic/F");
	outtree->Branch("jet1_dR2Mean", &t->jet1_dR2Mean, "jet1_dR2Mean/F");
	outtree->Branch("jet2_pt", &t->jet2_pt, "jet2_pt/F");
	outtree->Branch("jet2_e", &t->jet2_e, "jet2_e/F");
	outtree->Branch("jet2_phi", &t->jet2_phi, "jet2_phi/F");
	outtree->Branch("jet2_eta", &t->jet2_eta, "jet2_eta/F");
	outtree->Branch("jet2_mass", &t->jet2_mass, "jet2_mass/F");
	outtree->Branch("jet2_csvBtag", &t->jet2_csvBtag, "jet2_csvBtag/F");
	outtree->Branch("jet2_btagSF_M", &t->jet2_btagSF_M, "jet2_btagSF_M/F");
	outtree->Branch("jet2_betaStarClassic", &t->jet2_betaStarClassic, "jet2_betaStarClassic/F");
	outtree->Branch("jet2_dR2Mean", &t->jet2_dR2Mean, "jet2_dR2Mean/F");
// storing inputs of the regression for comparison
	outtree->Branch("regjet1_emfrac", &t->regjet1_emfrac, "regjet1_emfrac/F");
	outtree->Branch("regjet1_hadfrac", &t->regjet1_hadfrac, "regjet1_hadfrac/F");
	outtree->Branch("regjet1_secVtxPt", &t->regjet1_secVtxPt, "regjet1_secVtxPt/F");
	outtree->Branch("regjet1_secVtx3dL", &t->regjet1_secVtx3dL, "regjet1_secVtx3dL/F");
	outtree->Branch("regjet1_dPhiMet", &t->regjet1_dPhiMet, "regjet1_dPhiMet/F");
	outtree->Branch("regjet1_nConstituents", &t->regjet1_nConstituents, "regjet1_nConstituents/I");
	outtree->Branch("regjet2_emfrac", &t->regjet2_emfrac, "regjet2_emfrac/F");
	outtree->Branch("regjet2_hadfrac", &t->regjet2_hadfrac, "regjet2_hadfrac/F");
	outtree->Branch("regjet2_secVtxPt", &t->regjet2_secVtxPt, "regjet2_secVtxPt/F");
	outtree->Branch("regjet2_secVtx3dL", &t->regjet2_secVtx3dL, "regjet2_secVtx3dL/F");
	outtree->Branch("regjet2_dPhiMet", &t->regjet2_dPhiMet, "regjet2_dPhiMet/F");
	outtree->Branch("regjet2_nConstituents", &t->regjet2_nConstituents, "regjet2_nConstituents/I");
// regressed / kin fitted jets
	outtree->Branch("regjet1_pt", &t->regjet1_pt, "regjet1_pt/F");
	outtree->Branch("regjet1_e", &t->regjet1_e, "regjet1_e/F");
	outtree->Branch("regjet1_phi", &t->regjet1_phi, "regjet1_phi/F");
	outtree->Branch("regjet1_eta", &t->regjet1_eta, "regjet1_eta/F");
	outtree->Branch("regjet1_mass", &t->regjet1_mass, "regjet1_mass/F");
	outtree->Branch("regjet1_csvBtag", &t->regjet1_csvBtag, "regjet1_csvBtag/F");
	outtree->Branch("regjet1_btagSF_M", &t->regjet1_btagSF_M, "regjet1_btagSF_M/F");
	outtree->Branch("regjet1_betaStarClassic", &t->regjet1_betaStarClassic, "regjet1_betaStarClassic/F");
	outtree->Branch("regjet1_dR2Mean", &t->regjet1_dR2Mean, "regjet1_dR2Mean/F");
	outtree->Branch("regjet2_pt", &t->regjet2_pt, "regjet2_pt/F");
	outtree->Branch("regjet2_e", &t->regjet2_e, "regjet2_e/F");
	outtree->Branch("regjet2_phi", &t->regjet2_phi, "regjet2_phi/F");
	outtree->Branch("regjet2_eta", &t->regjet2_eta, "regjet2_eta/F");
	outtree->Branch("regjet2_mass", &t->regjet2_mass, "regjet2_mass/F");
	outtree->Branch("regjet2_csvBtag", &t->regjet2_csvBtag, "regjet2_csvBtag/F");
	outtree->Branch("regjet2_btagSF_M", &t->regjet2_btagSF_M, "regjet2_btagSF_M/F");
	outtree->Branch("regjet2_betaStarClassic", &t->regjet2_betaStarClassic, "regjet2_betaStarClassic/F");
	outtree->Branch("regjet2_dR2Mean", &t->regjet2_dR2Mean, "regjet2_dR2Mean/F");
	outtree->Branch("regkinjet1_pt", &t->regkinjet1_pt, "regkinjet1_pt/F");
	outtree->Branch("regkinjet1_e", &t->regkinjet1_e, "regkinjet1_e/F");
	outtree->Branch("regkinjet1_phi", &t->regkinjet1_phi, "regkinjet1_phi/F");
	outtree->Branch("regkinjet1_eta", &t->regkinjet1_eta, "regkinjet1_eta/F");
	outtree->Branch("regkinjet1_mass", &t->regkinjet1_mass, "regkinjet1_mass/F");
	outtree->Branch("regkinjet1_csvBtag", &t->regkinjet1_csvBtag, "regkinjet1_csvBtag/F");
	outtree->Branch("regkinjet1_btagSF_M", &t->regkinjet1_btagSF_M, "regkinjet1_btagSF_M/F");
	outtree->Branch("regkinjet1_betaStarClassic", &t->regkinjet1_betaStarClassic, "regkinjet1_betaStarClassic/F");
	outtree->Branch("regkinjet1_dR2Mean", &t->regkinjet1_dR2Mean, "regkinjet1_dR2Mean/F");
	outtree->Branch("regkinjet2_pt", &t->regkinjet2_pt, "regkinjet2_pt/F");
	outtree->Branch("regkinjet2_e", &t->regkinjet2_e, "regkinjet2_e/F");
	outtree->Branch("regkinjet2_phi", &t->regkinjet2_phi, "regkinjet2_phi/F");
	outtree->Branch("regkinjet2_eta", &t->regkinjet2_eta, "regkinjet2_eta/F");
	outtree->Branch("regkinjet2_mass", &t->regkinjet2_mass, "regkinjet2_mass/F");
	outtree->Branch("regkinjet2_csvBtag", &t->regkinjet2_csvBtag, "regkinjet2_csvBtag/F");
	outtree->Branch("regkinjet2_btagSF_M", &t->regkinjet2_btagSF_M, "regkinjet2_btagSF_M/F");
	outtree->Branch("regkinjet2_betaStarClassic", &t->regkinjet2_betaStarClassic, "regkinjet2_betaStarClassic/F");
	outtree->Branch("regkinjet2_dR2Mean", &t->regkinjet2_dR2Mean, "regkinjet2_dR2Mean/F");
	outtree->Branch("kinjet1_pt", &t->kinjet1_pt, "kinjet1_pt/F");
	outtree->Branch("kinjet1_e", &t->kinjet1_e, "kinjet1_e/F");
	outtree->Branch("kinjet1_phi", &t->kinjet1_phi, "kinjet1_phi/F");
	outtree->Branch("kinjet1_eta", &t->kinjet1_eta, "kinjet1_eta/F");
	outtree->Branch("kinjet1_mass", &t->kinjet1_mass, "kinjet1_mass/F");
	outtree->Branch("kinjet1_csvBtag", &t->kinjet1_csvBtag, "kinjet1_csvBtag/F");
	outtree->Branch("kinjet1_btagSF_M", &t->kinjet1_btagSF_M, "kinjet1_btagSF_M/F");
	outtree->Branch("kinjet1_betaStarClassic", &t->kinjet1_betaStarClassic, "kinjet1_betaStarClassic/F");
	outtree->Branch("kinjet1_dR2Mean", &t->kinjet1_dR2Mean, "kinjet1_dR2Mean/F");
	outtree->Branch("kinjet2_pt", &t->kinjet2_pt, "kinjet2_pt/F");
	outtree->Branch("kinjet2_e", &t->kinjet2_e, "kinjet2_e/F");
	outtree->Branch("kinjet2_phi", &t->kinjet2_phi, "kinjet2_phi/F");
	outtree->Branch("kinjet2_eta", &t->kinjet2_eta, "kinjet2_eta/F");
	outtree->Branch("kinjet2_mass", &t->kinjet2_mass, "kinjet2_mass/F");
	outtree->Branch("kinjet2_csvBtag", &t->kinjet2_csvBtag, "kinjet2_csvBtag/F");
	outtree->Branch("kinjet2_btagSF_M", &t->kinjet2_btagSF_M, "kinjet2_btagSF_M/F");
	outtree->Branch("kinjet2_betaStarClassic", &t->kinjet2_betaStarClassic, "kinjet2_betaStarClassic/F");
	outtree->Branch("kinjet2_dR2Mean", &t->kinjet2_dR2Mean, "kinjet2_dR2Mean/F");
	outtree->Branch("jj_pt", &t->jj_pt, "jj_pt/F");
	outtree->Branch("jj_e", &t->jj_e, "jj_e/F");
	outtree->Branch("jj_phi", &t->jj_phi, "jj_phi/F");
	outtree->Branch("jj_eta", &t->jj_eta, "jj_eta/F");
	outtree->Branch("jj_mass", &t->jj_mass, "jj_mass/F");
	outtree->Branch("jj_btagSF_M", &t->jj_btagSF_M, "jj_btagSF_M/F");
	outtree->Branch("jj_DR", &t->jj_DR, "jj_DR/F");
	outtree->Branch("regjj_pt", &t->regjj_pt, "regjj_pt/F");
	outtree->Branch("regjj_e", &t->regjj_e, "regjj_e/F");
	outtree->Branch("regjj_phi", &t->regjj_phi, "regjj_phi/F");
	outtree->Branch("regjj_eta", &t->regjj_eta, "regjj_eta/F");
	outtree->Branch("regjj_mass", &t->regjj_mass, "regjj_mass/F");
	outtree->Branch("regjj_btagSF_M", &t->regjj_btagSF_M, "regjj_btagSF_M/F");
	outtree->Branch("regjj_DR", &t->regjj_DR, "regjj_DR/F");
	outtree->Branch("regkinjj_pt", &t->regkinjj_pt, "regkinjj_pt/F");
	outtree->Branch("regkinjj_e", &t->regkinjj_e, "regkinjj_e/F");
	outtree->Branch("regkinjj_phi", &t->regkinjj_phi, "regkinjj_phi/F");
	outtree->Branch("regkinjj_eta", &t->regkinjj_eta, "regkinjj_eta/F");
	outtree->Branch("regkinjj_mass", &t->regkinjj_mass, "regkinjj_mass/F");
	outtree->Branch("regkinjj_btagSF_M", &t->regkinjj_btagSF_M, "regkinjj_btagSF_M/F");
	outtree->Branch("regkinjj_DR", &t->regkinjj_DR, "regkinjj_DR/F");
	outtree->Branch("kinjj_pt", &t->kinjj_pt, "kinjj_pt/F");
	outtree->Branch("kinjj_e", &t->kinjj_e, "kinjj_e/F");
	outtree->Branch("kinjj_phi", &t->kinjj_phi, "kinjj_phi/F");
	outtree->Branch("kinjj_eta", &t->kinjj_eta, "kinjj_eta/F");
	outtree->Branch("kinjj_mass", &t->kinjj_mass, "kinjj_mass/F");
	outtree->Branch("kinjj_btagSF_M", &t->kinjj_btagSF_M, "kinjj_btagSF_M/F");
	outtree->Branch("kinjj_DR", &t->kinjj_DR, "kinjj_DR/F");
	outtree->Branch("gg_pt", &t->gg_pt, "gg_pt/F");
	outtree->Branch("gg_e", &t->gg_e, "gg_e/F");
	outtree->Branch("gg_phi", &t->gg_phi, "gg_phi/F");
	outtree->Branch("gg_eta", &t->gg_eta, "gg_eta/F");
	outtree->Branch("gg_mass", &t->gg_mass, "gg_mass/F");
	outtree->Branch("ggjj_pt", &t->ggjj_pt, "ggjj_pt/F");
	outtree->Branch("ggjj_e", &t->ggjj_e, "ggjj_e/F");
	outtree->Branch("ggjj_phi", &t->ggjj_phi, "ggjj_phi/F");
	outtree->Branch("ggjj_eta", &t->ggjj_eta, "ggjj_eta/F");
	outtree->Branch("ggjj_mass", &t->ggjj_mass, "ggjj_mass/F");
	outtree->Branch("regggjj_pt", &t->regggjj_pt, "regggjj_pt/F");
	outtree->Branch("regggjj_e", &t->regggjj_e, "regggjj_e/F");
	outtree->Branch("regggjj_phi", &t->regggjj_phi, "regggjj_phi/F");
	outtree->Branch("regggjj_eta", &t->regggjj_eta, "regggjj_eta/F");
	outtree->Branch("regggjj_mass", &t->regggjj_mass, "regggjj_mass/F");
	outtree->Branch("regkinggjj_pt", &t->regkinggjj_pt, "regkinggjj_pt/F");
	outtree->Branch("regkinggjj_e", &t->regkinggjj_e, "regkinggjj_e/F");
	outtree->Branch("regkinggjj_phi", &t->regkinggjj_phi, "regkinggjj_phi/F");
	outtree->Branch("regkinggjj_eta", &t->regkinggjj_eta, "regkinggjj_eta/F");
	outtree->Branch("regkinggjj_mass", &t->regkinggjj_mass, "regkinggjj_mass/F");
	outtree->Branch("kinggjj_pt", &t->kinggjj_pt, "kinggjj_pt/F");
	outtree->Branch("kinggjj_e", &t->kinggjj_e, "kinggjj_e/F");
	outtree->Branch("kinggjj_phi", &t->kinggjj_phi, "kinggjj_phi/F");
	outtree->Branch("kinggjj_eta", &t->kinggjj_eta, "kinggjj_eta/F");
	outtree->Branch("kinggjj_mass", &t->kinggjj_mass, "kinggjj_mass/F");
	outtree->Branch("njets_kLooseID", &t->njets_kLooseID, "njets_kLooseID/I");
	outtree->Branch("njets_kLooseID_and_CSVM", &t->njets_kLooseID_and_CSVM, "njets_kLooseID_and_CSVM/I");
	outtree->Branch("njets_kRadionID", &t->njets_kRadionID, "njets_kRadionID/I");
	outtree->Branch("njets_kRadionID_and_CSVM", &t->njets_kRadionID_and_CSVM, "njets_kRadionID_and_CSVM/I");
	outtree->Branch("costhetastar", &t->costhetastar, "costhetastar/F");
	outtree->Branch("regcosthetastar", &t->regcosthetastar, "regcosthetastar/F");
	outtree->Branch("regkincosthetastar", &t->regkincosthetastar, "regkincosthetastar/F");
	outtree->Branch("kincosthetastar", &t->kincosthetastar, "kincosthetastar/F");
	outtree->Branch("minDRgj", &t->minDRgj, "minDRgj/F");
	outtree->Branch("minDRgregj", &t->minDRgregj, "minDRgregj/F");
	outtree->Branch("minDRgregkinj", &t->minDRgregkinj, "minDRgregkinj/F");
	outtree->Branch("minDRgkinj", &t->minDRgkinj, "minDRgkinj/F");
	outtree->Branch("HT_gg", &t->HT_gg, "HT_gg/F");
// gen level info
	outtree->Branch("gr_radion_p4_pt", &t->gr_radion_p4_pt, "gr_radion_p4_pt/F");
	outtree->Branch("gr_radion_p4_eta", &t->gr_radion_p4_eta, "gr_radion_p4_eta/F");
	outtree->Branch("gr_radion_p4_phi", &t->gr_radion_p4_phi, "gr_radion_p4_phi/F");
	outtree->Branch("gr_radion_p4_mass", &t->gr_radion_p4_mass, "gr_radion_p4_mass/F");
	outtree->Branch("gr_hgg_p4_pt", &t->gr_hgg_p4_pt, "gr_hgg_p4_pt/F");
	outtree->Branch("gr_hgg_p4_eta", &t->gr_hgg_p4_eta, "gr_hgg_p4_eta/F");
	outtree->Branch("gr_hgg_p4_phi", &t->gr_hgg_p4_phi, "gr_hgg_p4_phi/F");
	outtree->Branch("gr_hgg_p4_mass", &t->gr_hgg_p4_mass, "gr_hgg_p4_mass/F");
	outtree->Branch("gr_hbb_p4_pt", &t->gr_hbb_p4_pt, "gr_hbb_p4_pt/F");
	outtree->Branch("gr_hbb_p4_eta", &t->gr_hbb_p4_eta, "gr_hbb_p4_eta/F");
	outtree->Branch("gr_hbb_p4_phi", &t->gr_hbb_p4_phi, "gr_hbb_p4_phi/F");
	outtree->Branch("gr_hbb_p4_mass", &t->gr_hbb_p4_mass, "gr_hbb_p4_mass/F");
	outtree->Branch("gr_hjj_p4_pt", &t->gr_hjj_p4_pt, "gr_hjj_p4_pt/F");
	outtree->Branch("gr_hjj_p4_eta", &t->gr_hjj_p4_eta, "gr_hjj_p4_eta/F");
	outtree->Branch("gr_hjj_p4_phi", &t->gr_hjj_p4_phi, "gr_hjj_p4_phi/F");
	outtree->Branch("gr_hjj_p4_mass", &t->gr_hjj_p4_mass, "gr_hjj_p4_mass/F");
	outtree->Branch("gr_g1_p4_pt", &t->gr_g1_p4_pt, "gr_g1_p4_pt/F");
	outtree->Branch("gr_g1_p4_eta", &t->gr_g1_p4_eta, "gr_g1_p4_eta/F");
	outtree->Branch("gr_g1_p4_phi", &t->gr_g1_p4_phi, "gr_g1_p4_phi/F");
	outtree->Branch("gr_g1_p4_mass", &t->gr_g1_p4_mass, "gr_g1_p4_mass/F");
	outtree->Branch("gr_g2_p4_pt", &t->gr_g2_p4_pt, "gr_g2_p4_pt/F");
	outtree->Branch("gr_g2_p4_eta", &t->gr_g2_p4_eta, "gr_g2_p4_eta/F");
	outtree->Branch("gr_g2_p4_phi", &t->gr_g2_p4_phi, "gr_g2_p4_phi/F");
	outtree->Branch("gr_g2_p4_mass", &t->gr_g2_p4_mass, "gr_g2_p4_mass/F");
	outtree->Branch("gr_b1_p4_pt", &t->gr_b1_p4_pt, "gr_b1_p4_pt/F");
	outtree->Branch("gr_b1_p4_eta", &t->gr_b1_p4_eta, "gr_b1_p4_eta/F");
	outtree->Branch("gr_b1_p4_phi", &t->gr_b1_p4_phi, "gr_b1_p4_phi/F");
	outtree->Branch("gr_b1_p4_mass", &t->gr_b1_p4_mass, "gr_b1_p4_mass/F");
	outtree->Branch("gr_b2_p4_pt", &t->gr_b2_p4_pt, "gr_b2_p4_pt/F");
	outtree->Branch("gr_b2_p4_eta", &t->gr_b2_p4_eta, "gr_b2_p4_eta/F");
	outtree->Branch("gr_b2_p4_phi", &t->gr_b2_p4_phi, "gr_b2_p4_phi/F");
	outtree->Branch("gr_b2_p4_mass", &t->gr_b2_p4_mass, "gr_b2_p4_mass/F");
	outtree->Branch("gr_j1_p4_pt", &t->gr_j1_p4_pt, "gr_j1_p4_pt/F");
	outtree->Branch("gr_j1_p4_eta", &t->gr_j1_p4_eta, "gr_j1_p4_eta/F");
	outtree->Branch("gr_j1_p4_phi", &t->gr_j1_p4_phi, "gr_j1_p4_phi/F");
	outtree->Branch("gr_j1_p4_mass", &t->gr_j1_p4_mass, "gr_j1_p4_mass/F");
	outtree->Branch("gr_j2_p4_pt", &t->gr_j2_p4_pt, "gr_j2_p4_pt/F");
	outtree->Branch("gr_j2_p4_eta", &t->gr_j2_p4_eta, "gr_j2_p4_eta/F");
	outtree->Branch("gr_j2_p4_phi", &t->gr_j2_p4_phi, "gr_j2_p4_phi/F");
	outtree->Branch("gr_j2_p4_mass", &t->gr_j2_p4_mass, "gr_j2_p4_mass/F");

	return;
}

struct jet_variables
{
		std::vector<float> jetPt;
		std::vector<float> jetbtagSF_M;
		std::vector<float> jetbetaStarClassic;
		std::vector<float> jetdR2Mean;
		std::vector<float> jetE;
		std::vector<float> jetEta;
		std::vector<float> jetPhi;
		std::vector<float> jetCSV;
		std::vector<float> jetRegPt;
		std::vector<float> jetRegKinPt;
// regression inputs
		std::vector<float> jetEmfrac;
		std::vector<float> jetHadfrac;
		std::vector<float> jetSecVtxPt;
		std::vector<float> jetSecVtx3dL;
		std::vector<float> jetDPhiMet;
		std::vector<int> jetNConstituents;
};


void initialize_jet_variables( jet_variables * J )
{
		J->jetPt.clear();
		J->jetbtagSF_M.clear();
		J->jetbetaStarClassic.clear();
		J->jetdR2Mean.clear();
		J->jetE.clear();
		J->jetEta.clear();
		J->jetPhi.clear();
		J->jetCSV.clear();
		J->jetRegPt.clear();
		J->jetRegKinPt.clear();
// regression inputs
		J->jetEmfrac.clear();
		J->jetHadfrac.clear();
		J->jetSecVtxPt.clear();
		J->jetSecVtx3dL.clear();
		J->jetDPhiMet.clear();
		J->jetNConstituents.clear();

	return ;
}


void fill_jet_variables( tree_variables * t, int ijet, TLorentzVector met)
{
		TLorentzVector jet;
			if( ijet == 0 )
			{
				t->jet_e = t->j1_e;
				t->jet_pt = t->j1_pt;
				t->jet_phi = t->j1_phi;
				t->jet_eta = t->j1_eta;
				t->jet_betaStarClassic = t->j1_betaStarClassic;
				t->jet_dR2Mean = t->j1_dR2Mean;
				t->jet_csvBtag = t->j1_csvBtag;
				t->jet_secVtxPt = t->j1_secVtxPt;
				t->jet_secVtx3dL = t->j1_secVtx3dL;
				t->jet_emfrac = t->j1_emfrac;
				t->jet_hadfrac = t->j1_hadfrac;
				t->jet_btagSF_M = t->j1_btagSF_M;
				t->jet_betaStarClassic = t->j1_betaStarClassic;
				t->jet_dR2Mean = t->j1_dR2Mean;
				t->jet_nNeutrals = t->j1_nNeutrals;
				t->jet_nCharged = t->j1_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j1_pt, t->j1_eta, t->j1_phi, t->j1_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 1 )
			{
				t->jet_e = t->j2_e;
				t->jet_pt = t->j2_pt;
				t->jet_phi = t->j2_phi;
				t->jet_eta = t->j2_eta;
				t->jet_betaStarClassic = t->j2_betaStarClassic;
				t->jet_dR2Mean = t->j2_dR2Mean;
				t->jet_csvBtag = t->j2_csvBtag;
				t->jet_secVtxPt = t->j2_secVtxPt;
				t->jet_secVtx3dL = t->j2_secVtx3dL;
				t->jet_emfrac = t->j2_emfrac;
				t->jet_hadfrac = t->j2_hadfrac;
				t->jet_btagSF_M = t->j2_btagSF_M;
				t->jet_betaStarClassic = t->j2_betaStarClassic;
				t->jet_dR2Mean = t->j2_dR2Mean;
				t->jet_nNeutrals = t->j2_nNeutrals;
				t->jet_nCharged = t->j2_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j2_pt, t->j2_eta, t->j2_phi, t->j2_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 2 )
			{
				t->jet_e = t->j3_e;
				t->jet_pt = t->j3_pt;
				t->jet_phi = t->j3_phi;
				t->jet_eta = t->j3_eta;
				t->jet_betaStarClassic = t->j3_betaStarClassic;
				t->jet_dR2Mean = t->j3_dR2Mean;
				t->jet_csvBtag = t->j3_csvBtag;
				t->jet_secVtxPt = t->j3_secVtxPt;
				t->jet_secVtx3dL = t->j3_secVtx3dL;
				t->jet_emfrac = t->j3_emfrac;
				t->jet_hadfrac = t->j3_hadfrac;
				t->jet_btagSF_M = t->j3_btagSF_M;
				t->jet_betaStarClassic = t->j3_betaStarClassic;
				t->jet_dR2Mean = t->j3_dR2Mean;
				t->jet_nNeutrals = t->j3_nNeutrals;
				t->jet_nCharged = t->j3_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j3_pt, t->j3_eta, t->j3_phi, t->j3_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 3 )
			{
				t->jet_e = t->j4_e;
				t->jet_pt = t->j4_pt;
				t->jet_phi = t->j4_phi;
				t->jet_eta = t->j4_eta;
				t->jet_betaStarClassic = t->j4_betaStarClassic;
				t->jet_dR2Mean = t->j4_dR2Mean;
				t->jet_csvBtag = t->j4_csvBtag;
				t->jet_secVtxPt = t->j4_secVtxPt;
				t->jet_secVtx3dL = t->j4_secVtx3dL;
				t->jet_emfrac = t->j4_emfrac;
				t->jet_hadfrac = t->j4_hadfrac;
				t->jet_btagSF_M = t->j4_btagSF_M;
				t->jet_betaStarClassic = t->j4_betaStarClassic;
				t->jet_dR2Mean = t->j4_dR2Mean;
				t->jet_nNeutrals = t->j4_nNeutrals;
				t->jet_nCharged = t->j4_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j4_pt, t->j4_eta, t->j4_phi, t->j4_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 4 )
			{
				t->jet_e = t->j5_e;
				t->jet_pt = t->j5_pt;
				t->jet_phi = t->j5_phi;
				t->jet_eta = t->j5_eta;
				t->jet_betaStarClassic = t->j5_betaStarClassic;
				t->jet_dR2Mean = t->j5_dR2Mean;
				t->jet_csvBtag = t->j5_csvBtag;
				t->jet_secVtxPt = t->j5_secVtxPt;
				t->jet_secVtx3dL = t->j5_secVtx3dL;
				t->jet_emfrac = t->j5_emfrac;
				t->jet_hadfrac = t->j5_hadfrac;
				t->jet_btagSF_M = t->j5_btagSF_M;
				t->jet_betaStarClassic = t->j5_betaStarClassic;
				t->jet_dR2Mean = t->j5_dR2Mean;
				t->jet_nNeutrals = t->j5_nNeutrals;
				t->jet_nCharged = t->j5_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j5_pt, t->j5_eta, t->j5_phi, t->j5_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 5 )
			{
				t->jet_e = t->j6_e;
				t->jet_pt = t->j6_pt;
				t->jet_phi = t->j6_phi;
				t->jet_eta = t->j6_eta;
				t->jet_betaStarClassic = t->j6_betaStarClassic;
				t->jet_dR2Mean = t->j6_dR2Mean;
				t->jet_csvBtag = t->j6_csvBtag;
				t->jet_secVtxPt = t->j6_secVtxPt;
				t->jet_secVtx3dL = t->j6_secVtx3dL;
				t->jet_emfrac = t->j6_emfrac;
				t->jet_hadfrac = t->j6_hadfrac;
				t->jet_btagSF_M = t->j6_btagSF_M;
				t->jet_betaStarClassic = t->j6_betaStarClassic;
				t->jet_dR2Mean = t->j6_dR2Mean;
				t->jet_nNeutrals = t->j6_nNeutrals;
				t->jet_nCharged = t->j6_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j6_pt, t->j6_eta, t->j6_phi, t->j6_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 6 )
			{
				t->jet_e = t->j7_e;
				t->jet_pt = t->j7_pt;
				t->jet_phi = t->j7_phi;
				t->jet_eta = t->j7_eta;
				t->jet_betaStarClassic = t->j7_betaStarClassic;
				t->jet_dR2Mean = t->j7_dR2Mean;
				t->jet_csvBtag = t->j7_csvBtag;
				t->jet_secVtxPt = t->j7_secVtxPt;
				t->jet_secVtx3dL = t->j7_secVtx3dL;
				t->jet_emfrac = t->j7_emfrac;
				t->jet_hadfrac = t->j7_hadfrac;
				t->jet_btagSF_M = t->j7_btagSF_M;
				t->jet_betaStarClassic = t->j7_betaStarClassic;
				t->jet_dR2Mean = t->j7_dR2Mean;
				t->jet_nNeutrals = t->j7_nNeutrals;
				t->jet_nCharged = t->j7_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j7_pt, t->j7_eta, t->j7_phi, t->j7_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 7 )
			{
				t->jet_e = t->j8_e;
				t->jet_pt = t->j8_pt;
				t->jet_phi = t->j8_phi;
				t->jet_eta = t->j8_eta;
				t->jet_betaStarClassic = t->j8_betaStarClassic;
				t->jet_dR2Mean = t->j8_dR2Mean;
				t->jet_csvBtag = t->j8_csvBtag;
				t->jet_secVtxPt = t->j8_secVtxPt;
				t->jet_secVtx3dL = t->j8_secVtx3dL;
				t->jet_emfrac = t->j8_emfrac;
				t->jet_hadfrac = t->j8_hadfrac;
				t->jet_btagSF_M = t->j8_btagSF_M;
				t->jet_betaStarClassic = t->j8_betaStarClassic;
				t->jet_dR2Mean = t->j8_dR2Mean;
				t->jet_nNeutrals = t->j8_nNeutrals;
				t->jet_nCharged = t->j8_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j8_pt, t->j8_eta, t->j8_phi, t->j8_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 8 )
			{
				t->jet_e = t->j9_e;
				t->jet_pt = t->j9_pt;
				t->jet_phi = t->j9_phi;
				t->jet_eta = t->j9_eta;
				t->jet_betaStarClassic = t->j9_betaStarClassic;
				t->jet_dR2Mean = t->j9_dR2Mean;
				t->jet_csvBtag = t->j9_csvBtag;
				t->jet_secVtxPt = t->j9_secVtxPt;
				t->jet_secVtx3dL = t->j9_secVtx3dL;
				t->jet_emfrac = t->j9_emfrac;
				t->jet_hadfrac = t->j9_hadfrac;
				t->jet_btagSF_M = t->j9_btagSF_M;
				t->jet_betaStarClassic = t->j9_betaStarClassic;
				t->jet_dR2Mean = t->j9_dR2Mean;
				t->jet_nNeutrals = t->j9_nNeutrals;
				t->jet_nCharged = t->j9_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j9_pt, t->j9_eta, t->j9_phi, t->j9_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 9 )
			{
				t->jet_e = t->j10_e;
				t->jet_pt = t->j10_pt;
				t->jet_phi = t->j10_phi;
				t->jet_eta = t->j10_eta;
				t->jet_betaStarClassic = t->j10_betaStarClassic;
				t->jet_dR2Mean = t->j10_dR2Mean;
				t->jet_csvBtag = t->j10_csvBtag;
				t->jet_secVtxPt = t->j10_secVtxPt;
				t->jet_secVtx3dL = t->j10_secVtx3dL;
				t->jet_emfrac = t->j10_emfrac;
				t->jet_hadfrac = t->j10_hadfrac;
				t->jet_btagSF_M = t->j10_btagSF_M;
				t->jet_betaStarClassic = t->j10_betaStarClassic;
				t->jet_dR2Mean = t->j10_dR2Mean;
				t->jet_nNeutrals = t->j10_nNeutrals;
				t->jet_nCharged = t->j10_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j10_pt, t->j10_eta, t->j10_phi, t->j10_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 10 )
			{
				t->jet_e = t->j11_e;
				t->jet_pt = t->j11_pt;
				t->jet_phi = t->j11_phi;
				t->jet_eta = t->j11_eta;
				t->jet_betaStarClassic = t->j11_betaStarClassic;
				t->jet_dR2Mean = t->j11_dR2Mean;
				t->jet_csvBtag = t->j11_csvBtag;
				t->jet_secVtxPt = t->j11_secVtxPt;
				t->jet_secVtx3dL = t->j11_secVtx3dL;
				t->jet_emfrac = t->j11_emfrac;
				t->jet_hadfrac = t->j11_hadfrac;
				t->jet_btagSF_M = t->j11_btagSF_M;
				t->jet_betaStarClassic = t->j11_betaStarClassic;
				t->jet_dR2Mean = t->j11_dR2Mean;
				t->jet_nNeutrals = t->j11_nNeutrals;
				t->jet_nCharged = t->j11_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j11_pt, t->j11_eta, t->j11_phi, t->j11_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 11 )
			{
				t->jet_e = t->j12_e;
				t->jet_pt = t->j12_pt;
				t->jet_phi = t->j12_phi;
				t->jet_eta = t->j12_eta;
				t->jet_betaStarClassic = t->j12_betaStarClassic;
				t->jet_dR2Mean = t->j12_dR2Mean;
				t->jet_csvBtag = t->j12_csvBtag;
				t->jet_secVtxPt = t->j12_secVtxPt;
				t->jet_secVtx3dL = t->j12_secVtx3dL;
				t->jet_emfrac = t->j12_emfrac;
				t->jet_hadfrac = t->j12_hadfrac;
				t->jet_btagSF_M = t->j12_btagSF_M;
				t->jet_betaStarClassic = t->j12_betaStarClassic;
				t->jet_dR2Mean = t->j12_dR2Mean;
				t->jet_nNeutrals = t->j12_nNeutrals;
				t->jet_nCharged = t->j12_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j12_pt, t->j12_eta, t->j12_phi, t->j12_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 12 )
			{
				t->jet_e = t->j13_e;
				t->jet_pt = t->j13_pt;
				t->jet_phi = t->j13_phi;
				t->jet_eta = t->j13_eta;
				t->jet_betaStarClassic = t->j13_betaStarClassic;
				t->jet_dR2Mean = t->j13_dR2Mean;
				t->jet_csvBtag = t->j13_csvBtag;
				t->jet_secVtxPt = t->j13_secVtxPt;
				t->jet_secVtx3dL = t->j13_secVtx3dL;
				t->jet_emfrac = t->j13_emfrac;
				t->jet_hadfrac = t->j13_hadfrac;
				t->jet_btagSF_M = t->j13_btagSF_M;
				t->jet_betaStarClassic = t->j13_betaStarClassic;
				t->jet_dR2Mean = t->j13_dR2Mean;
				t->jet_nNeutrals = t->j13_nNeutrals;
				t->jet_nCharged = t->j13_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j13_pt, t->j13_eta, t->j13_phi, t->j13_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 13 )
			{
				t->jet_e = t->j14_e;
				t->jet_pt = t->j14_pt;
				t->jet_phi = t->j14_phi;
				t->jet_eta = t->j14_eta;
				t->jet_betaStarClassic = t->j14_betaStarClassic;
				t->jet_dR2Mean = t->j14_dR2Mean;
				t->jet_csvBtag = t->j14_csvBtag;
				t->jet_secVtxPt = t->j14_secVtxPt;
				t->jet_secVtx3dL = t->j14_secVtx3dL;
				t->jet_emfrac = t->j14_emfrac;
				t->jet_hadfrac = t->j14_hadfrac;
				t->jet_btagSF_M = t->j14_btagSF_M;
				t->jet_betaStarClassic = t->j14_betaStarClassic;
				t->jet_dR2Mean = t->j14_dR2Mean;
				t->jet_nNeutrals = t->j14_nNeutrals;
				t->jet_nCharged = t->j14_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j14_pt, t->j14_eta, t->j14_phi, t->j14_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			else
			if( ijet == 14 )
			{
				t->jet_e = t->j15_e;
				t->jet_pt = t->j15_pt;
				t->jet_phi = t->j15_phi;
				t->jet_eta = t->j15_eta;
				t->jet_betaStarClassic = t->j15_betaStarClassic;
				t->jet_dR2Mean = t->j15_dR2Mean;
				t->jet_csvBtag = t->j15_csvBtag;
				t->jet_secVtxPt = t->j15_secVtxPt;
				t->jet_secVtx3dL = t->j15_secVtx3dL;
				t->jet_emfrac = t->j15_emfrac;
				t->jet_hadfrac = t->j15_hadfrac;
				t->jet_btagSF_M = t->j15_btagSF_M;
				t->jet_betaStarClassic = t->j15_betaStarClassic;
				t->jet_dR2Mean = t->j15_dR2Mean;
				t->jet_nNeutrals = t->j15_nNeutrals;
				t->jet_nCharged = t->j15_nCharged;
				t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
				jet.SetPtEtaPhiE(t->j15_pt, t->j15_eta, t->j15_phi, t->j15_e);
				t->jet_dPhiMet = jet.DeltaPhi(met);
			} 
			return;
}
