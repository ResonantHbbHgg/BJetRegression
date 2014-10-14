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
    float gr_hbbhgg_costhetastar_CS, gr_hjjhgg_costhetastar_CS;
    float gr_dEta_gg_bb, gr_dEta_gg_jj;
    float gr_dPhi_gg_bb, gr_dPhi_gg_jj;
    float gr_dR_gg_bb, gr_dR_gg_jj;
// event variables
    float met_corr_pfmet, met_corr_phi_pfmet, met_corr_eta_pfmet, met_corr_e_pfmet;
    float pu_n, nvtx, rho;
    float weight, evweight, pu_weight;
    int run, lumis, event;
    float ph1_SCEta, ph2_SCEta;
    float ev_weight, ev_evweight, ev_pu_weight;
    float evweight_w_btagSF, evweight_w_btagSF_reg;
// object variables
    float ph1_eta, ph2_eta, ph1_pt, ph2_pt, PhotonsMass, ph1_phi, ph2_phi, ph1_e, ph2_e, ph1_r9, ph2_r9, ph1_sieie, ph2_sieie, ph1_hoe, ph2_hoe;
    float ph1_r9_cic, ph2_r9_cic, ph1_IDmva, ph2_IDmva;
    float ph1_pfchargedisogood03, ph1_ecaliso, ph1_pfchargedisobad04, ph1_ecalisobad, ph1_badvtx_Et, ph1_isconv;
    float ph2_pfchargedisogood03, ph2_ecaliso, ph2_pfchargedisobad04, ph2_ecalisobad, ph2_badvtx_Et, ph2_isconv;
    int ph1_ciclevel, ph2_ciclevel, ph1_isEB, ph2_isEB;
    // Photon Energy Scale & Photon Energy Resolution
    float ph1_sigmaEoE, ph2_sigmaEoE;
    float ph1_pesD_e, ph2_pesD_e, ph1_pesU_e, ph2_pesU_e;
    float ph1_perD_e, ph2_perD_e, ph1_perU_e, ph2_perU_e;

    float j1_e, j1_pt, j1_phi, j1_eta, j1_beta, j1_betaStar, j1_betaStarClassic, j1_dR2Mean, j1_csvBtag, j1_csvMvaBtag, j1_jetProbBtag, j1_tcheBtag, j1_secVtxPt, j1_secVtx3dL, j1_secVtx3deL, j1_secVtxM, j1_emfrac, j1_hadfrac, j1_chadfrac, j1_nhadfrac, j1_phofrac, j1_mufrac, j1_elefrac, j1_JECUnc, j1_leadTrackPt, j1_softLeptPt, j1_softLeptPtRel, j1_softLeptDR, j1_axis1, j1_axis2, j1_pull, j1_btagSF_M, j1_btagEff_M, j1_btagSFErrorUp_M, j1_btagSFErrorDown_M, j1_btagEffError_M;
    int j1_ntk, j1_nNeutrals, j1_nCharged, j1_flavour, j1_softLeptIdLooseMu, j1_softLeptIdEle95, j1_cutbased_wp_level, j1_simple_wp_level, j1_full_wp_level;

    float j2_e, j2_pt, j2_phi, j2_eta, j2_beta, j2_betaStar, j2_betaStarClassic, j2_dR2Mean, j2_csvBtag, j2_csvMvaBtag, j2_jetProbBtag, j2_tcheBtag, j2_secVtxPt, j2_secVtx3dL, j2_secVtx3deL, j2_secVtxM, j2_emfrac, j2_hadfrac, j2_chadfrac, j2_nhadfrac, j2_phofrac, j2_mufrac, j2_elefrac, j2_JECUnc, j2_leadTrackPt, j2_softLeptPt, j2_softLeptPtRel, j2_softLeptDR, j2_axis1, j2_axis2, j2_pull, j2_btagSF_M, j2_btagEff_M, j2_btagSFErrorUp_M, j2_btagSFErrorDown_M, j2_btagEffError_M;
    int j2_ntk, j2_nNeutrals, j2_nCharged, j2_flavour, j2_softLeptIdLooseMu, j2_softLeptIdEle95, j2_cutbased_wp_level, j2_simple_wp_level, j2_full_wp_level;

    float j3_e, j3_pt, j3_phi, j3_eta, j3_beta, j3_betaStar, j3_betaStarClassic, j3_dR2Mean, j3_csvBtag, j3_csvMvaBtag, j3_jetProbBtag, j3_tcheBtag, j3_secVtxPt, j3_secVtx3dL, j3_secVtx3deL, j3_secVtxM, j3_emfrac, j3_hadfrac, j3_chadfrac, j3_nhadfrac, j3_phofrac, j3_mufrac, j3_elefrac, j3_JECUnc, j3_leadTrackPt, j3_softLeptPt, j3_softLeptPtRel, j3_softLeptDR, j3_axis1, j3_axis2, j3_pull, j3_btagSF_M, j3_btagEff_M, j3_btagSFErrorUp_M, j3_btagSFErrorDown_M, j3_btagEffError_M;
    int j3_ntk, j3_nNeutrals, j3_nCharged, j3_flavour, j3_softLeptIdLooseMu, j3_softLeptIdEle95, j3_cutbased_wp_level, j3_simple_wp_level, j3_full_wp_level;

    float j4_e, j4_pt, j4_phi, j4_eta, j4_beta, j4_betaStar, j4_betaStarClassic, j4_dR2Mean, j4_csvBtag, j4_csvMvaBtag, j4_jetProbBtag, j4_tcheBtag, j4_secVtxPt, j4_secVtx3dL, j4_secVtx3deL, j4_secVtxM, j4_emfrac, j4_hadfrac, j4_chadfrac, j4_nhadfrac, j4_phofrac, j4_mufrac, j4_elefrac, j4_JECUnc, j4_leadTrackPt, j4_softLeptPt, j4_softLeptPtRel, j4_softLeptDR, j4_axis1, j4_axis2, j4_pull, j4_btagSF_M, j4_btagEff_M, j4_btagSFErrorUp_M, j4_btagSFErrorDown_M, j4_btagEffError_M;
    int j4_ntk, j4_nNeutrals, j4_nCharged, j4_flavour, j4_softLeptIdLooseMu, j4_softLeptIdEle95, j4_cutbased_wp_level, j4_simple_wp_level, j4_full_wp_level;

    float j5_e, j5_pt, j5_phi, j5_eta, j5_beta, j5_betaStar, j5_betaStarClassic, j5_dR2Mean, j5_csvBtag, j5_csvMvaBtag, j5_jetProbBtag, j5_tcheBtag, j5_secVtxPt, j5_secVtx3dL, j5_secVtx3deL, j5_secVtxM, j5_emfrac, j5_hadfrac, j5_chadfrac, j5_nhadfrac, j5_phofrac, j5_mufrac, j5_elefrac, j5_JECUnc, j5_leadTrackPt, j5_softLeptPt, j5_softLeptPtRel, j5_softLeptDR, j5_axis1, j5_axis2, j5_pull, j5_btagSF_M, j5_btagEff_M, j5_btagSFErrorUp_M, j5_btagSFErrorDown_M, j5_btagEffError_M;
    int j5_ntk, j5_nNeutrals, j5_nCharged, j5_flavour, j5_softLeptIdLooseMu, j5_softLeptIdEle95, j5_cutbased_wp_level, j5_simple_wp_level, j5_full_wp_level;

    float j6_e, j6_pt, j6_phi, j6_eta, j6_beta, j6_betaStar, j6_betaStarClassic, j6_dR2Mean, j6_csvBtag, j6_csvMvaBtag, j6_jetProbBtag, j6_tcheBtag, j6_secVtxPt, j6_secVtx3dL, j6_secVtx3deL, j6_secVtxM, j6_emfrac, j6_hadfrac, j6_chadfrac, j6_nhadfrac, j6_phofrac, j6_mufrac, j6_elefrac, j6_JECUnc, j6_leadTrackPt, j6_softLeptPt, j6_softLeptPtRel, j6_softLeptDR, j6_axis1, j6_axis2, j6_pull, j6_btagSF_M, j6_btagEff_M, j6_btagSFErrorUp_M, j6_btagSFErrorDown_M, j6_btagEffError_M;
    int j6_ntk, j6_nNeutrals, j6_nCharged, j6_flavour, j6_softLeptIdLooseMu, j6_softLeptIdEle95, j6_cutbased_wp_level, j6_simple_wp_level, j6_full_wp_level;

    float j7_e, j7_pt, j7_phi, j7_eta, j7_beta, j7_betaStar, j7_betaStarClassic, j7_dR2Mean, j7_csvBtag, j7_csvMvaBtag, j7_jetProbBtag, j7_tcheBtag, j7_secVtxPt, j7_secVtx3dL, j7_secVtx3deL, j7_secVtxM, j7_emfrac, j7_hadfrac, j7_chadfrac, j7_nhadfrac, j7_phofrac, j7_mufrac, j7_elefrac, j7_JECUnc, j7_leadTrackPt, j7_softLeptPt, j7_softLeptPtRel, j7_softLeptDR, j7_axis1, j7_axis2, j7_pull, j7_btagSF_M, j7_btagEff_M, j7_btagSFErrorUp_M, j7_btagSFErrorDown_M, j7_btagEffError_M;
    int j7_ntk, j7_nNeutrals, j7_nCharged, j7_flavour, j7_softLeptIdLooseMu, j7_softLeptIdEle95, j7_cutbased_wp_level, j7_simple_wp_level, j7_full_wp_level;

    float j8_e, j8_pt, j8_phi, j8_eta, j8_beta, j8_betaStar, j8_betaStarClassic, j8_dR2Mean, j8_csvBtag, j8_csvMvaBtag, j8_jetProbBtag, j8_tcheBtag, j8_secVtxPt, j8_secVtx3dL, j8_secVtx3deL, j8_secVtxM, j8_emfrac, j8_hadfrac, j8_chadfrac, j8_nhadfrac, j8_phofrac, j8_mufrac, j8_elefrac, j8_JECUnc, j8_leadTrackPt, j8_softLeptPt, j8_softLeptPtRel, j8_softLeptDR, j8_axis1, j8_axis2, j8_pull, j8_btagSF_M, j8_btagEff_M, j8_btagSFErrorUp_M, j8_btagSFErrorDown_M, j8_btagEffError_M;
    int j8_ntk, j8_nNeutrals, j8_nCharged, j8_flavour, j8_softLeptIdLooseMu, j8_softLeptIdEle95, j8_cutbased_wp_level, j8_simple_wp_level, j8_full_wp_level;

    float j9_e, j9_pt, j9_phi, j9_eta, j9_beta, j9_betaStar, j9_betaStarClassic, j9_dR2Mean, j9_csvBtag, j9_csvMvaBtag, j9_jetProbBtag, j9_tcheBtag, j9_secVtxPt, j9_secVtx3dL, j9_secVtx3deL, j9_secVtxM, j9_emfrac, j9_hadfrac, j9_chadfrac, j9_nhadfrac, j9_phofrac, j9_mufrac, j9_elefrac, j9_JECUnc, j9_leadTrackPt, j9_softLeptPt, j9_softLeptPtRel, j9_softLeptDR, j9_axis1, j9_axis2, j9_pull, j9_btagSF_M, j9_btagEff_M, j9_btagSFErrorUp_M, j9_btagSFErrorDown_M, j9_btagEffError_M;
    int j9_ntk, j9_nNeutrals, j9_nCharged, j9_flavour, j9_softLeptIdLooseMu, j9_softLeptIdEle95, j9_cutbased_wp_level, j9_simple_wp_level, j9_full_wp_level;

    float j10_e, j10_pt, j10_phi, j10_eta, j10_beta, j10_betaStar, j10_betaStarClassic, j10_dR2Mean, j10_csvBtag, j10_csvMvaBtag, j10_jetProbBtag, j10_tcheBtag, j10_secVtxPt, j10_secVtx3dL, j10_secVtx3deL, j10_secVtxM, j10_emfrac, j10_hadfrac, j10_chadfrac, j10_nhadfrac, j10_phofrac, j10_mufrac, j10_elefrac, j10_JECUnc, j10_leadTrackPt, j10_softLeptPt, j10_softLeptPtRel, j10_softLeptDR, j10_axis1, j10_axis2, j10_pull, j10_btagSF_M, j10_btagEff_M, j10_btagSFErrorUp_M, j10_btagSFErrorDown_M, j10_btagEffError_M;
    int j10_ntk, j10_nNeutrals, j10_nCharged, j10_flavour, j10_softLeptIdLooseMu, j10_softLeptIdEle95, j10_cutbased_wp_level, j10_simple_wp_level, j10_full_wp_level;

    float j11_e, j11_pt, j11_phi, j11_eta, j11_beta, j11_betaStar, j11_betaStarClassic, j11_dR2Mean, j11_csvBtag, j11_csvMvaBtag, j11_jetProbBtag, j11_tcheBtag, j11_secVtxPt, j11_secVtx3dL, j11_secVtx3deL, j11_secVtxM, j11_emfrac, j11_hadfrac, j11_chadfrac, j11_nhadfrac, j11_phofrac, j11_mufrac, j11_elefrac, j11_JECUnc, j11_leadTrackPt, j11_softLeptPt, j11_softLeptPtRel, j11_softLeptDR, j11_axis1, j11_axis2, j11_pull, j11_btagSF_M, j11_btagEff_M, j11_btagSFErrorUp_M, j11_btagSFErrorDown_M, j11_btagEffError_M;
    int j11_ntk, j11_nNeutrals, j11_nCharged, j11_flavour, j11_softLeptIdLooseMu, j11_softLeptIdEle95, j11_cutbased_wp_level, j11_simple_wp_level, j11_full_wp_level;

    float j12_e, j12_pt, j12_phi, j12_eta, j12_beta, j12_betaStar, j12_betaStarClassic, j12_dR2Mean, j12_csvBtag, j12_csvMvaBtag, j12_jetProbBtag, j12_tcheBtag, j12_secVtxPt, j12_secVtx3dL, j12_secVtx3deL, j12_secVtxM, j12_emfrac, j12_hadfrac, j12_chadfrac, j12_nhadfrac, j12_phofrac, j12_mufrac, j12_elefrac, j12_JECUnc, j12_leadTrackPt, j12_softLeptPt, j12_softLeptPtRel, j12_softLeptDR, j12_axis1, j12_axis2, j12_pull, j12_btagSF_M, j12_btagEff_M, j12_btagSFErrorUp_M, j12_btagSFErrorDown_M, j12_btagEffError_M;
    int j12_ntk, j12_nNeutrals, j12_nCharged, j12_flavour, j12_softLeptIdLooseMu, j12_softLeptIdEle95, j12_cutbased_wp_level, j12_simple_wp_level, j12_full_wp_level;

    float j13_e, j13_pt, j13_phi, j13_eta, j13_beta, j13_betaStar, j13_betaStarClassic, j13_dR2Mean, j13_csvBtag, j13_csvMvaBtag, j13_jetProbBtag, j13_tcheBtag, j13_secVtxPt, j13_secVtx3dL, j13_secVtx3deL, j13_secVtxM, j13_emfrac, j13_hadfrac, j13_chadfrac, j13_nhadfrac, j13_phofrac, j13_mufrac, j13_elefrac, j13_JECUnc, j13_leadTrackPt, j13_softLeptPt, j13_softLeptPtRel, j13_softLeptDR, j13_axis1, j13_axis2, j13_pull, j13_btagSF_M, j13_btagEff_M, j13_btagSFErrorUp_M, j13_btagSFErrorDown_M, j13_btagEffError_M;
    int j13_ntk, j13_nNeutrals, j13_nCharged, j13_flavour, j13_softLeptIdLooseMu, j13_softLeptIdEle95, j13_cutbased_wp_level, j13_simple_wp_level, j13_full_wp_level;

    float j14_e, j14_pt, j14_phi, j14_eta, j14_beta, j14_betaStar, j14_betaStarClassic, j14_dR2Mean, j14_csvBtag, j14_csvMvaBtag, j14_jetProbBtag, j14_tcheBtag, j14_secVtxPt, j14_secVtx3dL, j14_secVtx3deL, j14_secVtxM, j14_emfrac, j14_hadfrac, j14_chadfrac, j14_nhadfrac, j14_phofrac, j14_mufrac, j14_elefrac, j14_JECUnc, j14_leadTrackPt, j14_softLeptPt, j14_softLeptPtRel, j14_softLeptDR, j14_axis1, j14_axis2, j14_pull, j14_btagSF_M, j14_btagEff_M, j14_btagSFErrorUp_M, j14_btagSFErrorDown_M, j14_btagEffError_M;
    int j14_ntk, j14_nNeutrals, j14_nCharged, j14_flavour, j14_softLeptIdLooseMu, j14_softLeptIdEle95, j14_cutbased_wp_level, j14_simple_wp_level, j14_full_wp_level;

    float j15_e, j15_pt, j15_phi, j15_eta, j15_beta, j15_betaStar, j15_betaStarClassic, j15_dR2Mean, j15_csvBtag, j15_csvMvaBtag, j15_jetProbBtag, j15_tcheBtag, j15_secVtxPt, j15_secVtx3dL, j15_secVtx3deL, j15_secVtxM, j15_emfrac, j15_hadfrac, j15_chadfrac, j15_nhadfrac, j15_phofrac, j15_mufrac, j15_elefrac, j15_JECUnc, j15_leadTrackPt, j15_softLeptPt, j15_softLeptPtRel, j15_softLeptDR, j15_axis1, j15_axis2, j15_pull, j15_btagSF_M, j15_btagEff_M, j15_btagSFErrorUp_M, j15_btagSFErrorDown_M, j15_btagEffError_M;
    int j15_ntk, j15_nNeutrals, j15_nCharged, j15_flavour, j15_softLeptIdLooseMu, j15_softLeptIdEle95, j15_cutbased_wp_level, j15_simple_wp_level, j15_full_wp_level;

    float jet_e, jet_pt, jet_phi, jet_eta;
    float jet_betaStarClassic, jet_dR2Mean, jet_csvBtag;
    float jet_mt, jet_secVtxPt, jet_secVtx3dL, jet_secVtx3deL, jet_secVtxM, jet_emfrac, jet_hadfrac, jet_chadfrac, jet_nhadfrac, jet_phofrac, jet_mufrac, jet_elefrac, jet_JECUnc, jet_leadTrackPt, jet_softLeptPt, jet_softLeptPtRel, jet_softLeptDR, jet_btagSF_M, jet_flavour, jet_btagEff_M, jet_btagSFErrorUp_M, jet_btagSFErrorDown_M, jet_btagEffError_M;
    int jet_nNeutrals, jet_nCharged, jet_nConstituents, jet_cutbased_wp_level, jet_simple_wp_level, jet_full_wp_level;

    float jet_nConstituents_;
    float jet_dPhiMet_fabs;
    float jet_dPhiMet, jet_regPt, jet_regkinPt;
    // Jet Energy Correction and Jet Energy Resolution
    float jet_jecD_e, jet_jecD_pt, jet_jecD_phi, jet_jecD_eta;
    float jet_jecU_e, jet_jecU_pt, jet_jecU_phi, jet_jecU_eta;
    float jet_jerD_e, jet_jerD_pt, jet_jerD_phi, jet_jerD_eta;
    float jet_jerC_e, jet_jerC_pt, jet_jerC_phi, jet_jerC_eta;
    float jet_jerU_e, jet_jerU_pt, jet_jerU_phi, jet_jerU_eta;

    float j1_jecD_e, j1_jecD_pt, j1_jecD_phi, j1_jecD_eta;
    float j1_jecU_e, j1_jecU_pt, j1_jecU_phi, j1_jecU_eta;
    float j1_jerD_e, j1_jerD_pt, j1_jerD_phi, j1_jerD_eta;
    float j1_jerC_e, j1_jerC_pt, j1_jerC_phi, j1_jerC_eta;
    float j1_jerU_e, j1_jerU_pt, j1_jerU_phi, j1_jerU_eta;

    float j2_jecD_e, j2_jecD_pt, j2_jecD_phi, j2_jecD_eta;
    float j2_jecU_e, j2_jecU_pt, j2_jecU_phi, j2_jecU_eta;
    float j2_jerD_e, j2_jerD_pt, j2_jerD_phi, j2_jerD_eta;
    float j2_jerC_e, j2_jerC_pt, j2_jerC_phi, j2_jerC_eta;
    float j2_jerU_e, j2_jerU_pt, j2_jerU_phi, j2_jerU_eta;

    float j3_jecD_e, j3_jecD_pt, j3_jecD_phi, j3_jecD_eta;
    float j3_jecU_e, j3_jecU_pt, j3_jecU_phi, j3_jecU_eta;
    float j3_jerD_e, j3_jerD_pt, j3_jerD_phi, j3_jerD_eta;
    float j3_jerC_e, j3_jerC_pt, j3_jerC_phi, j3_jerC_eta;
    float j3_jerU_e, j3_jerU_pt, j3_jerU_phi, j3_jerU_eta;

    float j4_jecD_e, j4_jecD_pt, j4_jecD_phi, j4_jecD_eta;
    float j4_jecU_e, j4_jecU_pt, j4_jecU_phi, j4_jecU_eta;
    float j4_jerD_e, j4_jerD_pt, j4_jerD_phi, j4_jerD_eta;
    float j4_jerC_e, j4_jerC_pt, j4_jerC_phi, j4_jerC_eta;
    float j4_jerU_e, j4_jerU_pt, j4_jerU_phi, j4_jerU_eta;

    float j5_jecD_e, j5_jecD_pt, j5_jecD_phi, j5_jecD_eta;
    float j5_jecU_e, j5_jecU_pt, j5_jecU_phi, j5_jecU_eta;
    float j5_jerD_e, j5_jerD_pt, j5_jerD_phi, j5_jerD_eta;
    float j5_jerC_e, j5_jerC_pt, j5_jerC_phi, j5_jerC_eta;
    float j5_jerU_e, j5_jerU_pt, j5_jerU_phi, j5_jerU_eta;

    float j6_jecD_e, j6_jecD_pt, j6_jecD_phi, j6_jecD_eta;
    float j6_jecU_e, j6_jecU_pt, j6_jecU_phi, j6_jecU_eta;
    float j6_jerD_e, j6_jerD_pt, j6_jerD_phi, j6_jerD_eta;
    float j6_jerC_e, j6_jerC_pt, j6_jerC_phi, j6_jerC_eta;
    float j6_jerU_e, j6_jerU_pt, j6_jerU_phi, j6_jerU_eta;

    float j7_jecD_e, j7_jecD_pt, j7_jecD_phi, j7_jecD_eta;
    float j7_jecU_e, j7_jecU_pt, j7_jecU_phi, j7_jecU_eta;
    float j7_jerD_e, j7_jerD_pt, j7_jerD_phi, j7_jerD_eta;
    float j7_jerC_e, j7_jerC_pt, j7_jerC_phi, j7_jerC_eta;
    float j7_jerU_e, j7_jerU_pt, j7_jerU_phi, j7_jerU_eta;

    float j8_jecD_e, j8_jecD_pt, j8_jecD_phi, j8_jecD_eta;
    float j8_jecU_e, j8_jecU_pt, j8_jecU_phi, j8_jecU_eta;
    float j8_jerD_e, j8_jerD_pt, j8_jerD_phi, j8_jerD_eta;
    float j8_jerC_e, j8_jerC_pt, j8_jerC_phi, j8_jerC_eta;
    float j8_jerU_e, j8_jerU_pt, j8_jerU_phi, j8_jerU_eta;

    float j9_jecD_e, j9_jecD_pt, j9_jecD_phi, j9_jecD_eta;
    float j9_jecU_e, j9_jecU_pt, j9_jecU_phi, j9_jecU_eta;
    float j9_jerD_e, j9_jerD_pt, j9_jerD_phi, j9_jerD_eta;
    float j9_jerC_e, j9_jerC_pt, j9_jerC_phi, j9_jerC_eta;
    float j9_jerU_e, j9_jerU_pt, j9_jerU_phi, j9_jerU_eta;

    float j10_jecD_e, j10_jecD_pt, j10_jecD_phi, j10_jecD_eta;
    float j10_jecU_e, j10_jecU_pt, j10_jecU_phi, j10_jecU_eta;
    float j10_jerD_e, j10_jerD_pt, j10_jerD_phi, j10_jerD_eta;
    float j10_jerC_e, j10_jerC_pt, j10_jerC_phi, j10_jerC_eta;
    float j10_jerU_e, j10_jerU_pt, j10_jerU_phi, j10_jerU_eta;

    float j11_jecD_e, j11_jecD_pt, j11_jecD_phi, j11_jecD_eta;
    float j11_jecU_e, j11_jecU_pt, j11_jecU_phi, j11_jecU_eta;
    float j11_jerD_e, j11_jerD_pt, j11_jerD_phi, j11_jerD_eta;
    float j11_jerC_e, j11_jerC_pt, j11_jerC_phi, j11_jerC_eta;
    float j11_jerU_e, j11_jerU_pt, j11_jerU_phi, j11_jerU_eta;

    float j12_jecD_e, j12_jecD_pt, j12_jecD_phi, j12_jecD_eta;
    float j12_jecU_e, j12_jecU_pt, j12_jecU_phi, j12_jecU_eta;
    float j12_jerD_e, j12_jerD_pt, j12_jerD_phi, j12_jerD_eta;
    float j12_jerC_e, j12_jerC_pt, j12_jerC_phi, j12_jerC_eta;
    float j12_jerU_e, j12_jerU_pt, j12_jerU_phi, j12_jerU_eta;

    float j13_jecD_e, j13_jecD_pt, j13_jecD_phi, j13_jecD_eta;
    float j13_jecU_e, j13_jecU_pt, j13_jecU_phi, j13_jecU_eta;
    float j13_jerD_e, j13_jerD_pt, j13_jerD_phi, j13_jerD_eta;
    float j13_jerC_e, j13_jerC_pt, j13_jerC_phi, j13_jerC_eta;
    float j13_jerU_e, j13_jerU_pt, j13_jerU_phi, j13_jerU_eta;

    float j14_jecD_e, j14_jecD_pt, j14_jecD_phi, j14_jecD_eta;
    float j14_jecU_e, j14_jecU_pt, j14_jecU_phi, j14_jecU_eta;
    float j14_jerD_e, j14_jerD_pt, j14_jerD_phi, j14_jerD_eta;
    float j14_jerC_e, j14_jerC_pt, j14_jerC_phi, j14_jerC_eta;
    float j14_jerU_e, j14_jerU_pt, j14_jerU_phi, j14_jerU_eta;

    float j15_jecD_e, j15_jecD_pt, j15_jecD_phi, j15_jecD_eta;
    float j15_jecU_e, j15_jecU_pt, j15_jecU_phi, j15_jecU_eta;
    float j15_jerD_e, j15_jerD_pt, j15_jerD_phi, j15_jerD_eta;
    float j15_jerC_e, j15_jerC_pt, j15_jerC_phi, j15_jerC_eta;
    float j15_jerU_e, j15_jerU_pt, j15_jerU_phi, j15_jerU_eta;


// setup tree outputs
    float vtx_z;
    float pho1_pt, pho1_e, pho1_phi, pho1_eta, pho1_mass, pho1_r9, pho1_sieie, pho1_hoe, pho1_r9_cic, pho1_IDmva;
    float pho2_pt, pho2_e, pho2_phi, pho2_eta, pho2_mass, pho2_r9, pho2_sieie, pho2_hoe, pho2_r9_cic, pho2_IDmva;
    int pho1_isEB, pho2_isEB;
    float pho1_pfchargedisogood03, pho1_ecaliso, pho1_pfchargedisobad04, pho1_ecalisobad, pho1_badvtx_Et, pho1_PFisoA, pho1_PFisoB, pho1_PFisoC, pho1_isconv;
    float pho2_pfchargedisogood03, pho2_ecaliso, pho2_pfchargedisobad04, pho2_ecalisobad, pho2_badvtx_Et, pho2_PFisoA, pho2_PFisoB, pho2_PFisoC, pho2_isconv;
    float jet1_pt, jet1_e, jet1_phi, jet1_eta, jet1_mass, jet1_csvBtag, jet1_btagSF_M, jet1_btagEff_M, jet1_btagSFErrorUp_M, jet1_btagSFErrorDown_M, jet1_btagEffError_M, jet1_betaStarClassic, jet1_dR2Mean;
    float jet2_pt, jet2_e, jet2_phi, jet2_eta, jet2_mass, jet2_csvBtag, jet2_btagSF_M, jet2_btagEff_M, jet2_btagSFErrorUp_M, jet2_btagSFErrorDown_M, jet2_btagEffError_M, jet2_betaStarClassic, jet2_dR2Mean;
    int jet1_flavour, jet1_cutbased_wp_level, jet1_simple_wp_level, jet1_full_wp_level;
    int jet2_flavour, jet2_cutbased_wp_level, jet2_simple_wp_level, jet2_full_wp_level;

    float regjet1_mt, regjet1_chadfrac, regjet1_nhadfrac, regjet1_phofrac, regjet1_mufrac, regjet1_elefrac, regjet1_softLeptPt, regjet1_softLeptPtRel, regjet1_softLeptDR, regjet1_leadTrackPt, regjet1_JECUnc, regjet1_secVtxPt, regjet1_secVtx3dL, regjet1_secVtx3deL, regjet1_secVtxM, regjet1_dPhiMet;
    int regjet1_nConstituents;
    float regjet2_mt, regjet2_chadfrac, regjet2_nhadfrac, regjet2_phofrac, regjet2_mufrac, regjet2_elefrac, regjet2_softLeptPt, regjet2_softLeptPtRel, regjet2_softLeptDR, regjet2_leadTrackPt, regjet2_JECUnc, regjet2_secVtxPt, regjet2_secVtx3dL, regjet2_secVtx3deL, regjet2_secVtxM, regjet2_dPhiMet;
    int regjet2_nConstituents;
    float regjet1_pt, regjet1_e, regjet1_phi, regjet1_eta, regjet1_mass, regjet1_csvBtag, regjet1_btagSF_M, regjet1_btagEff_M, regjet1_btagSFErrorUp_M, regjet1_btagSFErrorDown_M, regjet1_btagEffError_M, regjet1_betaStarClassic, regjet1_dR2Mean;
    float regjet2_pt, regjet2_e, regjet2_phi, regjet2_eta, regjet2_mass, regjet2_csvBtag, regjet2_btagSF_M, regjet2_btagEff_M, regjet2_btagSFErrorUp_M, regjet2_btagSFErrorDown_M, regjet2_btagEffError_M, regjet2_betaStarClassic, regjet2_dR2Mean;
    int regjet1_flavour, regjet1_cutbased_wp_level, regjet1_simple_wp_level, regjet1_full_wp_level;
    int regjet2_flavour, regjet2_cutbased_wp_level, regjet2_simple_wp_level, regjet2_full_wp_level;

    float regkinjet1_pt, regkinjet1_e, regkinjet1_phi, regkinjet1_eta, regkinjet1_mass, regkinjet1_csvBtag, regkinjet1_btagSF_M, regkinjet1_btagEff_M, regkinjet1_btagSFErrorUp_M, regkinjet1_btagSFErrorDown_M, regkinjet1_btagEffError_M, regkinjet1_betaStarClassic, regkinjet1_dR2Mean;
    float regkinjet2_pt, regkinjet2_e, regkinjet2_phi, regkinjet2_eta, regkinjet2_mass, regkinjet2_csvBtag, regkinjet2_btagSF_M, regkinjet2_btagEff_M, regkinjet2_btagSFErrorUp_M, regkinjet2_btagSFErrorDown_M, regkinjet2_btagEffError_M, regkinjet2_betaStarClassic, regkinjet2_dR2Mean;
    int regkinjet1_flavour, regkinjet1_cutbased_wp_level, regkinjet1_simple_wp_level, regkinjet1_full_wp_level;
    int regkinjet2_flavour, regkinjet2_cutbased_wp_level, regkinjet2_simple_wp_level, regkinjet2_full_wp_level;

    float kinjet1_pt, kinjet1_e, kinjet1_phi, kinjet1_eta, kinjet1_mass, kinjet1_csvBtag, kinjet1_btagSF_M, kinjet1_btagEff_M, kinjet1_btagSFErrorUp_M, kinjet1_btagSFErrorDown_M, kinjet1_btagEffError_M, kinjet1_betaStarClassic, kinjet1_dR2Mean;
    float kinjet2_pt, kinjet2_e, kinjet2_phi, kinjet2_eta, kinjet2_mass, kinjet2_csvBtag, kinjet2_btagSF_M, kinjet2_btagEff_M, kinjet2_btagSFErrorUp_M, kinjet2_btagSFErrorDown_M, kinjet2_btagEffError_M, kinjet2_betaStarClassic, kinjet2_dR2Mean;
    int kinjet1_flavour, kinjet1_cutbased_wp_level, kinjet1_simple_wp_level, kinjet1_full_wp_level;
    int kinjet2_flavour, kinjet2_cutbased_wp_level, kinjet2_simple_wp_level, kinjet2_full_wp_level;

    float jj_pt, jj_e, jj_phi, jj_eta, jj_mass, jj_DR, jj_btagSF_M;
    float regjj_pt, regjj_e, regjj_phi, regjj_eta, regjj_mass, regjj_btagSF_M;
    float regkinjj_pt, regkinjj_e, regkinjj_phi, regkinjj_eta, regkinjj_mass, regkinjj_btagSF_M;
    float kinjj_pt, kinjj_e, kinjj_phi, kinjj_eta, kinjj_mass, kinjj_btagSF_M;
    float gg_pt, gg_e, gg_phi, gg_eta, gg_mass, gg_DR;
    float ggjj_pt, ggjj_e, ggjj_phi, ggjj_eta, ggjj_mass, regjj_DR, regkinjj_DR, kinjj_DR;
    float regggjj_pt, regggjj_e, regggjj_phi, regggjj_eta, regggjj_mass;
    float regkinggjj_pt, regkinggjj_e, regkinggjj_phi, regkinggjj_eta, regkinggjj_mass;
    float kinggjj_pt, kinggjj_e, kinggjj_phi, kinggjj_eta, kinggjj_mass;
    int selection_cut_level;
    int category;
    float costhetastar, regcosthetastar, regkincosthetastar, kincosthetastar;
    float costhetastar_CS, regcosthetastar_CS, regkincosthetastar_CS, kincosthetastar_CS;
    float minDRgj, minDRgregj, minDRgregkinj, minDRgkinj;
    float HT_gg;
    float dEta_gg_jj, dEta_gg_regjj, dEta_gg_regkinjj, dEta_gg_kinjj;
    float dPhi_gg_jj, dPhi_gg_regjj, dPhi_gg_regkinjj, dPhi_gg_kinjj;
    float dR_gg_jj, dR_gg_regjj, dR_gg_regkinjj, dR_gg_kinjj;
    // Photon Energy Scale & Photon Energy Resolution
    float gg_mass_pesD1, gg_mass_pesU1, gg_mass_perD1, gg_mass_perU1;
    float gg_mass_pesD2, gg_mass_pesU2, gg_mass_perD2, gg_mass_perU2;
    float gg_mass_pesD, gg_mass_pesU, gg_mass_perD, gg_mass_perU;
    float ggjj_mass_pesD1, ggjj_mass_pesU1, ggjj_mass_perD1, ggjj_mass_perU1;
    float ggjj_mass_pesD2, ggjj_mass_pesU2, ggjj_mass_perD2, ggjj_mass_perU2;
    float ggjj_mass_pesD, ggjj_mass_pesU, ggjj_mass_perD, ggjj_mass_perU;
    float kinggjj_mass_pesD1, kinggjj_mass_pesU1, kinggjj_mass_perD1, kinggjj_mass_perU1;
    float kinggjj_mass_pesD2, kinggjj_mass_pesU2, kinggjj_mass_perD2, kinggjj_mass_perU2;
    float kinggjj_mass_pesD, kinggjj_mass_pesU, kinggjj_mass_perD, kinggjj_mass_perU;
    float regggjj_mass_pesD1, regggjj_mass_pesU1, regggjj_mass_perD1, regggjj_mass_perU1;
    float regggjj_mass_pesD2, regggjj_mass_pesU2, regggjj_mass_perD2, regggjj_mass_perU2;
    float regggjj_mass_pesD, regggjj_mass_pesU, regggjj_mass_perD, regggjj_mass_perU;
    float regkinggjj_mass_pesD1, regkinggjj_mass_pesU1, regkinggjj_mass_perD1, regkinggjj_mass_perU1;
    float regkinggjj_mass_pesD2, regkinggjj_mass_pesU2, regkinggjj_mass_perD2, regkinggjj_mass_perU2;
    float regkinggjj_mass_pesD, regkinggjj_mass_pesU, regkinggjj_mass_perD, regkinggjj_mass_perU;
    // Jet Energy Correction and Jet Energy Resolution
    float jj_mass_jecD, jj_mass_jecU, jj_mass_jerD, jj_mass_jerC, jj_mass_jerU;
    float kinjj_mass_jecD, kinjj_mass_jecU, kinjj_mass_jerD, kinjj_mass_jerC, kinjj_mass_jerU;
    float ggjj_mass_jecD, ggjj_mass_jecU, ggjj_mass_jerD, ggjj_mass_jerC, ggjj_mass_jerU;
    float kinggjj_mass_jecD, kinggjj_mass_jecU, kinggjj_mass_jerD, kinggjj_mass_jerC, kinggjj_mass_jerU;


    int njets_passing_kLooseID;
    int njets_passing_kLooseID_and_CSVM;
    int njets_kLooseID, njets_kRadionID;
    int njets_kLooseID_and_CSVM, njets_kRadionID_and_CSVM;
};

void initialize_variables(tree_variables *t)
{
    // to protect gen level variables
    // in case they are not present (ie data or MC backgrounds)
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
    t->gr_hbbhgg_costhetastar_CS = t->gr_hjjhgg_costhetastar_CS = 0.;
    t->gr_dEta_gg_bb = t->gr_dEta_gg_jj = 0.;
    t->gr_dPhi_gg_bb = t->gr_dPhi_gg_jj = 0.;
    t->gr_dR_gg_bb = t->gr_dR_gg_jj = 0.;
    // to protect jet variables related to the regression
    // in case they are not present (ie processing old samples)
    t->j1_secVtxM = t->j1_emfrac = t->j1_hadfrac = t->j1_chadfrac = t->j1_nhadfrac = t->j1_phofrac = t->j1_mufrac = t->j1_elefrac = t->j1_JECUnc = t->j1_leadTrackPt = t->j1_softLeptPt = t->j1_softLeptPtRel = t->j1_softLeptDR = t->j1_softLeptIdLooseMu = t->j1_softLeptIdEle95 = 0.;
    t->j2_secVtxM = t->j2_emfrac = t->j2_hadfrac = t->j2_chadfrac = t->j2_nhadfrac = t->j2_phofrac = t->j2_mufrac = t->j2_elefrac = t->j2_JECUnc = t->j2_leadTrackPt = t->j2_softLeptPt = t->j2_softLeptPtRel = t->j2_softLeptDR = t->j2_softLeptIdLooseMu = t->j2_softLeptIdEle95 = 0.;
    t->j3_secVtxM = t->j3_emfrac = t->j3_hadfrac = t->j3_chadfrac = t->j3_nhadfrac = t->j3_phofrac = t->j3_mufrac = t->j3_elefrac = t->j3_JECUnc = t->j3_leadTrackPt = t->j3_softLeptPt = t->j3_softLeptPtRel = t->j3_softLeptDR = t->j3_softLeptIdLooseMu = t->j3_softLeptIdEle95 = 0.;
    t->j4_secVtxM = t->j4_emfrac = t->j4_hadfrac = t->j4_chadfrac = t->j4_nhadfrac = t->j4_phofrac = t->j4_mufrac = t->j4_elefrac = t->j4_JECUnc = t->j4_leadTrackPt = t->j4_softLeptPt = t->j4_softLeptPtRel = t->j4_softLeptDR = t->j4_softLeptIdLooseMu = t->j4_softLeptIdEle95 = 0.;
    t->j5_secVtxM = t->j5_emfrac = t->j5_hadfrac = t->j5_chadfrac = t->j5_nhadfrac = t->j5_phofrac = t->j5_mufrac = t->j5_elefrac = t->j5_JECUnc = t->j5_leadTrackPt = t->j5_softLeptPt = t->j5_softLeptPtRel = t->j5_softLeptDR = t->j5_softLeptIdLooseMu = t->j5_softLeptIdEle95 = 0.;
    t->j6_secVtxM = t->j6_emfrac = t->j6_hadfrac = t->j6_chadfrac = t->j6_nhadfrac = t->j6_phofrac = t->j6_mufrac = t->j6_elefrac = t->j6_JECUnc = t->j6_leadTrackPt = t->j6_softLeptPt = t->j6_softLeptPtRel = t->j6_softLeptDR = t->j6_softLeptIdLooseMu = t->j6_softLeptIdEle95 = 0.;
    t->j7_secVtxM = t->j7_emfrac = t->j7_hadfrac = t->j7_chadfrac = t->j7_nhadfrac = t->j7_phofrac = t->j7_mufrac = t->j7_elefrac = t->j7_JECUnc = t->j7_leadTrackPt = t->j7_softLeptPt = t->j7_softLeptPtRel = t->j7_softLeptDR = t->j7_softLeptIdLooseMu = t->j7_softLeptIdEle95 = 0.;
    t->j8_secVtxM = t->j8_emfrac = t->j8_hadfrac = t->j8_chadfrac = t->j8_nhadfrac = t->j8_phofrac = t->j8_mufrac = t->j8_elefrac = t->j8_JECUnc = t->j8_leadTrackPt = t->j8_softLeptPt = t->j8_softLeptPtRel = t->j8_softLeptDR = t->j8_softLeptIdLooseMu = t->j8_softLeptIdEle95 = 0.;
    t->j9_secVtxM = t->j9_emfrac = t->j9_hadfrac = t->j9_chadfrac = t->j9_nhadfrac = t->j9_phofrac = t->j9_mufrac = t->j9_elefrac = t->j9_JECUnc = t->j9_leadTrackPt = t->j9_softLeptPt = t->j9_softLeptPtRel = t->j9_softLeptDR = t->j9_softLeptIdLooseMu = t->j9_softLeptIdEle95 = 0.;
    t->j10_secVtxM = t->j10_emfrac = t->j10_hadfrac = t->j10_chadfrac = t->j10_nhadfrac = t->j10_phofrac = t->j10_mufrac = t->j10_elefrac = t->j10_JECUnc = t->j10_leadTrackPt = t->j10_softLeptPt = t->j10_softLeptPtRel = t->j10_softLeptDR = t->j10_softLeptIdLooseMu = t->j10_softLeptIdEle95 = 0.;
    t->j11_secVtxM = t->j11_emfrac = t->j11_hadfrac = t->j11_chadfrac = t->j11_nhadfrac = t->j11_phofrac = t->j11_mufrac = t->j11_elefrac = t->j11_JECUnc = t->j11_leadTrackPt = t->j11_softLeptPt = t->j11_softLeptPtRel = t->j11_softLeptDR = t->j11_softLeptIdLooseMu = t->j11_softLeptIdEle95 = 0.;
    t->j12_secVtxM = t->j12_emfrac = t->j12_hadfrac = t->j12_chadfrac = t->j12_nhadfrac = t->j12_phofrac = t->j12_mufrac = t->j12_elefrac = t->j12_JECUnc = t->j12_leadTrackPt = t->j12_softLeptPt = t->j12_softLeptPtRel = t->j12_softLeptDR = t->j12_softLeptIdLooseMu = t->j12_softLeptIdEle95 = 0.;
    t->j13_secVtxM = t->j13_emfrac = t->j13_hadfrac = t->j13_chadfrac = t->j13_nhadfrac = t->j13_phofrac = t->j13_mufrac = t->j13_elefrac = t->j13_JECUnc = t->j13_leadTrackPt = t->j13_softLeptPt = t->j13_softLeptPtRel = t->j13_softLeptDR = t->j13_softLeptIdLooseMu = t->j13_softLeptIdEle95 = 0.;
    t->j14_secVtxM = t->j14_emfrac = t->j14_hadfrac = t->j14_chadfrac = t->j14_nhadfrac = t->j14_phofrac = t->j14_mufrac = t->j14_elefrac = t->j14_JECUnc = t->j14_leadTrackPt = t->j14_softLeptPt = t->j14_softLeptPtRel = t->j14_softLeptDR = t->j14_softLeptIdLooseMu = t->j14_softLeptIdEle95 = 0.;
    t->j15_secVtxM = t->j15_emfrac = t->j15_hadfrac = t->j15_chadfrac = t->j15_nhadfrac = t->j15_phofrac = t->j15_mufrac = t->j15_elefrac = t->j15_JECUnc = t->j15_leadTrackPt = t->j15_softLeptPt = t->j15_softLeptPtRel = t->j15_softLeptDR = t->j15_softLeptIdLooseMu = t->j15_softLeptIdEle95 = 0.;

    // Other initalizations that we don't want to screw up
    t->selection_cut_level = 0;
    t->category = 0;

    return;
}

void setup_intree(TTree* intree, tree_variables *t, int type, int numberOfRegressionFiles)
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
    intree->SetBranchAddress("ph1_r9_cic", &t->ph1_r9_cic);
    intree->SetBranchAddress("ph2_r9_cic", &t->ph2_r9_cic);
    intree->SetBranchAddress("ph1_IDmva", &t->ph1_IDmva);
    intree->SetBranchAddress("ph2_IDmva", &t->ph2_IDmva);
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
    intree->SetBranchAddress("ph1_sigmaEoE", &t->ph1_sigmaEoE);
    intree->SetBranchAddress("ph2_sigmaEoE", &t->ph2_sigmaEoE);
    intree->SetBranchAddress("ph1_pesD_e", &t->ph1_pesD_e);
    intree->SetBranchAddress("ph2_pesD_e", &t->ph2_pesD_e);
    intree->SetBranchAddress("ph1_pesU_e", &t->ph1_pesU_e);
    intree->SetBranchAddress("ph2_pesU_e", &t->ph2_pesU_e);
    intree->SetBranchAddress("ph1_perD_e", &t->ph1_perD_e);
    intree->SetBranchAddress("ph2_perD_e", &t->ph2_perD_e);
    intree->SetBranchAddress("ph1_perU_e", &t->ph1_perU_e);
    intree->SetBranchAddress("ph2_perU_e", &t->ph2_perU_e);
    intree->SetBranchAddress("ph1_SCEta", &t->ph1_SCEta);
    intree->SetBranchAddress("ph2_SCEta", &t->ph2_SCEta);
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
    intree->SetBranchAddress("weight", &t->ev_weight);
    intree->SetBranchAddress("evweight", &t->ev_evweight);
    intree->SetBranchAddress("pu_weight", &t->ev_pu_weight);
    intree->SetBranchAddress("vtx_z", &t->vtx_z);

    if( type < 0 )
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
    intree->SetBranchAddress("j1_secVtx3deL", &t->j1_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j1_secVtxM", &t->j1_secVtxM);
        intree->SetBranchAddress("j1_emfrac", &t->j1_emfrac);
        intree->SetBranchAddress("j1_hadfrac", &t->j1_hadfrac);
        intree->SetBranchAddress("j1_chadfrac", &t->j1_chadfrac);
        intree->SetBranchAddress("j1_nhadfrac", &t->j1_nhadfrac);
        intree->SetBranchAddress("j1_phofrac", &t->j1_phofrac);
        intree->SetBranchAddress("j1_mufrac", &t->j1_mufrac);
        intree->SetBranchAddress("j1_elefrac", &t->j1_elefrac);
        intree->SetBranchAddress("j1_JECUnc", &t->j1_JECUnc);
        intree->SetBranchAddress("j1_leadTrackPt", &t->j1_leadTrackPt);
        intree->SetBranchAddress("j1_softLeptPt", &t->j1_softLeptPt);
        intree->SetBranchAddress("j1_softLeptPtRel", &t->j1_softLeptPtRel);
        intree->SetBranchAddress("j1_softLeptDR", &t->j1_softLeptDR);
        intree->SetBranchAddress("j1_softLeptIdLooseMu", &t->j1_softLeptIdLooseMu);
        intree->SetBranchAddress("j1_softLeptIdEle95", &t->j1_softLeptIdEle95);
    }
    intree->SetBranchAddress("j1_btagSF_M", &t->j1_btagSF_M);
    intree->SetBranchAddress("j1_flavour", &t->j1_flavour);
    intree->SetBranchAddress("j1_btagSFErrorUp_M", &t->j1_btagSFErrorUp_M);
    intree->SetBranchAddress("j1_btagSFErrorDown_M", &t->j1_btagSFErrorDown_M);
    intree->SetBranchAddress("j1_btagEff_M", &t->j1_btagEff_M);
    intree->SetBranchAddress("j1_btagEffError_M", &t->j1_btagEffError_M);
    intree->SetBranchAddress("j1_cutbased_wp_level", &t->j1_cutbased_wp_level);
    intree->SetBranchAddress("j1_simple_wp_level", &t->j1_simple_wp_level);
    intree->SetBranchAddress("j1_full_wp_level", &t->j1_full_wp_level);
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
    intree->SetBranchAddress("j2_secVtx3deL", &t->j2_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j2_secVtxM", &t->j2_secVtxM);
        intree->SetBranchAddress("j2_emfrac", &t->j2_emfrac);
        intree->SetBranchAddress("j2_hadfrac", &t->j2_hadfrac);
        intree->SetBranchAddress("j2_chadfrac", &t->j2_chadfrac);
        intree->SetBranchAddress("j2_nhadfrac", &t->j2_nhadfrac);
        intree->SetBranchAddress("j2_phofrac", &t->j2_phofrac);
        intree->SetBranchAddress("j2_mufrac", &t->j2_mufrac);
        intree->SetBranchAddress("j2_elefrac", &t->j2_elefrac);
        intree->SetBranchAddress("j2_JECUnc", &t->j2_JECUnc);
        intree->SetBranchAddress("j2_leadTrackPt", &t->j2_leadTrackPt);
        intree->SetBranchAddress("j2_softLeptPt", &t->j2_softLeptPt);
        intree->SetBranchAddress("j2_softLeptPtRel", &t->j2_softLeptPtRel);
        intree->SetBranchAddress("j2_softLeptDR", &t->j2_softLeptDR);
        intree->SetBranchAddress("j2_softLeptIdLooseMu", &t->j2_softLeptIdLooseMu);
        intree->SetBranchAddress("j2_softLeptIdEle95", &t->j2_softLeptIdEle95);
    }
    intree->SetBranchAddress("j2_btagSF_M", &t->j2_btagSF_M);
    intree->SetBranchAddress("j2_flavour", &t->j2_flavour);
    intree->SetBranchAddress("j2_btagSFErrorUp_M", &t->j2_btagSFErrorUp_M);
    intree->SetBranchAddress("j2_btagSFErrorDown_M", &t->j2_btagSFErrorDown_M);
    intree->SetBranchAddress("j2_btagEff_M", &t->j2_btagEff_M);
    intree->SetBranchAddress("j2_btagEffError_M", &t->j2_btagEffError_M);
    intree->SetBranchAddress("j2_cutbased_wp_level", &t->j2_cutbased_wp_level);
    intree->SetBranchAddress("j2_simple_wp_level", &t->j2_simple_wp_level);
    intree->SetBranchAddress("j2_full_wp_level", &t->j2_full_wp_level);
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
    intree->SetBranchAddress("j3_secVtx3deL", &t->j3_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j3_secVtxM", &t->j3_secVtxM);
        intree->SetBranchAddress("j3_emfrac", &t->j3_emfrac);
        intree->SetBranchAddress("j3_hadfrac", &t->j3_hadfrac);
        intree->SetBranchAddress("j3_chadfrac", &t->j3_chadfrac);
        intree->SetBranchAddress("j3_nhadfrac", &t->j3_nhadfrac);
        intree->SetBranchAddress("j3_phofrac", &t->j3_phofrac);
        intree->SetBranchAddress("j3_mufrac", &t->j3_mufrac);
        intree->SetBranchAddress("j3_elefrac", &t->j3_elefrac);
        intree->SetBranchAddress("j3_JECUnc", &t->j3_JECUnc);
        intree->SetBranchAddress("j3_leadTrackPt", &t->j3_leadTrackPt);
        intree->SetBranchAddress("j3_softLeptPt", &t->j3_softLeptPt);
        intree->SetBranchAddress("j3_softLeptPtRel", &t->j3_softLeptPtRel);
        intree->SetBranchAddress("j3_softLeptDR", &t->j3_softLeptDR);
        intree->SetBranchAddress("j3_softLeptIdLooseMu", &t->j3_softLeptIdLooseMu);
        intree->SetBranchAddress("j3_softLeptIdEle95", &t->j3_softLeptIdEle95);
    }
    intree->SetBranchAddress("j3_btagSF_M", &t->j3_btagSF_M);
    intree->SetBranchAddress("j3_flavour", &t->j3_flavour);
    intree->SetBranchAddress("j3_btagSFErrorUp_M", &t->j3_btagSFErrorUp_M);
    intree->SetBranchAddress("j3_btagSFErrorDown_M", &t->j3_btagSFErrorDown_M);
    intree->SetBranchAddress("j3_btagEff_M", &t->j3_btagEff_M);
    intree->SetBranchAddress("j3_btagEffError_M", &t->j3_btagEffError_M);
    intree->SetBranchAddress("j3_cutbased_wp_level", &t->j3_cutbased_wp_level);
    intree->SetBranchAddress("j3_simple_wp_level", &t->j3_simple_wp_level);
    intree->SetBranchAddress("j3_full_wp_level", &t->j3_full_wp_level);
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
    intree->SetBranchAddress("j4_secVtx3deL", &t->j4_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j4_secVtxM", &t->j4_secVtxM);
        intree->SetBranchAddress("j4_emfrac", &t->j4_emfrac);
        intree->SetBranchAddress("j4_hadfrac", &t->j4_hadfrac);
        intree->SetBranchAddress("j4_chadfrac", &t->j4_chadfrac);
        intree->SetBranchAddress("j4_nhadfrac", &t->j4_nhadfrac);
        intree->SetBranchAddress("j4_phofrac", &t->j4_phofrac);
        intree->SetBranchAddress("j4_mufrac", &t->j4_mufrac);
        intree->SetBranchAddress("j4_elefrac", &t->j4_elefrac);
        intree->SetBranchAddress("j4_JECUnc", &t->j4_JECUnc);
        intree->SetBranchAddress("j4_leadTrackPt", &t->j4_leadTrackPt);
        intree->SetBranchAddress("j4_softLeptPt", &t->j4_softLeptPt);
        intree->SetBranchAddress("j4_softLeptPtRel", &t->j4_softLeptPtRel);
        intree->SetBranchAddress("j4_softLeptDR", &t->j4_softLeptDR);
        intree->SetBranchAddress("j4_softLeptIdLooseMu", &t->j4_softLeptIdLooseMu);
        intree->SetBranchAddress("j4_softLeptIdEle95", &t->j4_softLeptIdEle95);
    }
    intree->SetBranchAddress("j4_btagSF_M", &t->j4_btagSF_M);
    intree->SetBranchAddress("j4_flavour", &t->j4_flavour);
    intree->SetBranchAddress("j4_btagSFErrorUp_M", &t->j4_btagSFErrorUp_M);
    intree->SetBranchAddress("j4_btagSFErrorDown_M", &t->j4_btagSFErrorDown_M);
    intree->SetBranchAddress("j4_btagEff_M", &t->j4_btagEff_M);
    intree->SetBranchAddress("j4_btagEffError_M", &t->j4_btagEffError_M);
    intree->SetBranchAddress("j4_cutbased_wp_level", &t->j4_cutbased_wp_level);
    intree->SetBranchAddress("j4_simple_wp_level", &t->j4_simple_wp_level);
    intree->SetBranchAddress("j4_full_wp_level", &t->j4_full_wp_level);
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
    intree->SetBranchAddress("j5_secVtx3deL", &t->j5_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j5_secVtxM", &t->j5_secVtxM);
        intree->SetBranchAddress("j5_emfrac", &t->j5_emfrac);
        intree->SetBranchAddress("j5_hadfrac", &t->j5_hadfrac);
        intree->SetBranchAddress("j5_chadfrac", &t->j5_chadfrac);
        intree->SetBranchAddress("j5_nhadfrac", &t->j5_nhadfrac);
        intree->SetBranchAddress("j5_phofrac", &t->j5_phofrac);
        intree->SetBranchAddress("j5_mufrac", &t->j5_mufrac);
        intree->SetBranchAddress("j5_elefrac", &t->j5_elefrac);
        intree->SetBranchAddress("j5_JECUnc", &t->j5_JECUnc);
        intree->SetBranchAddress("j5_leadTrackPt", &t->j5_leadTrackPt);
        intree->SetBranchAddress("j5_softLeptPt", &t->j5_softLeptPt);
        intree->SetBranchAddress("j5_softLeptPtRel", &t->j5_softLeptPtRel);
        intree->SetBranchAddress("j5_softLeptDR", &t->j5_softLeptDR);
        intree->SetBranchAddress("j5_softLeptIdLooseMu", &t->j5_softLeptIdLooseMu);
        intree->SetBranchAddress("j5_softLeptIdEle95", &t->j5_softLeptIdEle95);
    }
    intree->SetBranchAddress("j5_btagSF_M", &t->j5_btagSF_M);
    intree->SetBranchAddress("j5_flavour", &t->j5_flavour);
    intree->SetBranchAddress("j5_btagSFErrorUp_M", &t->j5_btagSFErrorUp_M);
    intree->SetBranchAddress("j5_btagSFErrorDown_M", &t->j5_btagSFErrorDown_M);
    intree->SetBranchAddress("j5_btagEff_M", &t->j5_btagEff_M);
    intree->SetBranchAddress("j5_btagEffError_M", &t->j5_btagEffError_M);
    intree->SetBranchAddress("j5_cutbased_wp_level", &t->j5_cutbased_wp_level);
    intree->SetBranchAddress("j5_simple_wp_level", &t->j5_simple_wp_level);
    intree->SetBranchAddress("j5_full_wp_level", &t->j5_full_wp_level);
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
    intree->SetBranchAddress("j6_secVtx3deL", &t->j6_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j6_secVtxM", &t->j6_secVtxM);
        intree->SetBranchAddress("j6_emfrac", &t->j6_emfrac);
        intree->SetBranchAddress("j6_hadfrac", &t->j6_hadfrac);
        intree->SetBranchAddress("j6_chadfrac", &t->j6_chadfrac);
        intree->SetBranchAddress("j6_nhadfrac", &t->j6_nhadfrac);
        intree->SetBranchAddress("j6_phofrac", &t->j6_phofrac);
        intree->SetBranchAddress("j6_mufrac", &t->j6_mufrac);
        intree->SetBranchAddress("j6_elefrac", &t->j6_elefrac);
        intree->SetBranchAddress("j6_JECUnc", &t->j6_JECUnc);
        intree->SetBranchAddress("j6_leadTrackPt", &t->j6_leadTrackPt);
        intree->SetBranchAddress("j6_softLeptPt", &t->j6_softLeptPt);
        intree->SetBranchAddress("j6_softLeptPtRel", &t->j6_softLeptPtRel);
        intree->SetBranchAddress("j6_softLeptDR", &t->j6_softLeptDR);
        intree->SetBranchAddress("j6_softLeptIdLooseMu", &t->j6_softLeptIdLooseMu);
        intree->SetBranchAddress("j6_softLeptIdEle95", &t->j6_softLeptIdEle95);
    }
    intree->SetBranchAddress("j6_btagSF_M", &t->j6_btagSF_M);
    intree->SetBranchAddress("j6_flavour", &t->j6_flavour);
    intree->SetBranchAddress("j6_btagSFErrorUp_M", &t->j6_btagSFErrorUp_M);
    intree->SetBranchAddress("j6_btagSFErrorDown_M", &t->j6_btagSFErrorDown_M);
    intree->SetBranchAddress("j6_btagEff_M", &t->j6_btagEff_M);
    intree->SetBranchAddress("j6_btagEffError_M", &t->j6_btagEffError_M);
    intree->SetBranchAddress("j6_cutbased_wp_level", &t->j6_cutbased_wp_level);
    intree->SetBranchAddress("j6_simple_wp_level", &t->j6_simple_wp_level);
    intree->SetBranchAddress("j6_full_wp_level", &t->j6_full_wp_level);
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
    intree->SetBranchAddress("j7_secVtx3deL", &t->j7_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j7_secVtxM", &t->j7_secVtxM);
        intree->SetBranchAddress("j7_emfrac", &t->j7_emfrac);
        intree->SetBranchAddress("j7_hadfrac", &t->j7_hadfrac);
        intree->SetBranchAddress("j7_chadfrac", &t->j7_chadfrac);
        intree->SetBranchAddress("j7_nhadfrac", &t->j7_nhadfrac);
        intree->SetBranchAddress("j7_phofrac", &t->j7_phofrac);
        intree->SetBranchAddress("j7_mufrac", &t->j7_mufrac);
        intree->SetBranchAddress("j7_elefrac", &t->j7_elefrac);
        intree->SetBranchAddress("j7_JECUnc", &t->j7_JECUnc);
        intree->SetBranchAddress("j7_leadTrackPt", &t->j7_leadTrackPt);
        intree->SetBranchAddress("j7_softLeptPt", &t->j7_softLeptPt);
        intree->SetBranchAddress("j7_softLeptPtRel", &t->j7_softLeptPtRel);
        intree->SetBranchAddress("j7_softLeptDR", &t->j7_softLeptDR);
        intree->SetBranchAddress("j7_softLeptIdLooseMu", &t->j7_softLeptIdLooseMu);
        intree->SetBranchAddress("j7_softLeptIdEle95", &t->j7_softLeptIdEle95);
    }
    intree->SetBranchAddress("j7_btagSF_M", &t->j7_btagSF_M);
    intree->SetBranchAddress("j7_flavour", &t->j7_flavour);
    intree->SetBranchAddress("j7_btagSFErrorUp_M", &t->j7_btagSFErrorUp_M);
    intree->SetBranchAddress("j7_btagSFErrorDown_M", &t->j7_btagSFErrorDown_M);
    intree->SetBranchAddress("j7_btagEff_M", &t->j7_btagEff_M);
    intree->SetBranchAddress("j7_btagEffError_M", &t->j7_btagEffError_M);
    intree->SetBranchAddress("j7_cutbased_wp_level", &t->j7_cutbased_wp_level);
    intree->SetBranchAddress("j7_simple_wp_level", &t->j7_simple_wp_level);
    intree->SetBranchAddress("j7_full_wp_level", &t->j7_full_wp_level);
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
    intree->SetBranchAddress("j8_secVtx3deL", &t->j8_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j8_secVtxM", &t->j8_secVtxM);
        intree->SetBranchAddress("j8_emfrac", &t->j8_emfrac);
        intree->SetBranchAddress("j8_hadfrac", &t->j8_hadfrac);
        intree->SetBranchAddress("j8_chadfrac", &t->j8_chadfrac);
        intree->SetBranchAddress("j8_nhadfrac", &t->j8_nhadfrac);
        intree->SetBranchAddress("j8_phofrac", &t->j8_phofrac);
        intree->SetBranchAddress("j8_mufrac", &t->j8_mufrac);
        intree->SetBranchAddress("j8_elefrac", &t->j8_elefrac);
        intree->SetBranchAddress("j8_JECUnc", &t->j8_JECUnc);
        intree->SetBranchAddress("j8_leadTrackPt", &t->j8_leadTrackPt);
        intree->SetBranchAddress("j8_softLeptPt", &t->j8_softLeptPt);
        intree->SetBranchAddress("j8_softLeptPtRel", &t->j8_softLeptPtRel);
        intree->SetBranchAddress("j8_softLeptDR", &t->j8_softLeptDR);
        intree->SetBranchAddress("j8_softLeptIdLooseMu", &t->j8_softLeptIdLooseMu);
        intree->SetBranchAddress("j8_softLeptIdEle95", &t->j8_softLeptIdEle95);
    }
    intree->SetBranchAddress("j8_btagSF_M", &t->j8_btagSF_M);
    intree->SetBranchAddress("j8_flavour", &t->j8_flavour);
    intree->SetBranchAddress("j8_btagSFErrorUp_M", &t->j8_btagSFErrorUp_M);
    intree->SetBranchAddress("j8_btagSFErrorDown_M", &t->j8_btagSFErrorDown_M);
    intree->SetBranchAddress("j8_btagEff_M", &t->j8_btagEff_M);
    intree->SetBranchAddress("j8_btagEffError_M", &t->j8_btagEffError_M);
    intree->SetBranchAddress("j8_cutbased_wp_level", &t->j8_cutbased_wp_level);
    intree->SetBranchAddress("j8_simple_wp_level", &t->j8_simple_wp_level);
    intree->SetBranchAddress("j8_full_wp_level", &t->j8_full_wp_level);
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
    intree->SetBranchAddress("j9_secVtx3deL", &t->j9_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j9_secVtxM", &t->j9_secVtxM);
        intree->SetBranchAddress("j9_emfrac", &t->j9_emfrac);
        intree->SetBranchAddress("j9_hadfrac", &t->j9_hadfrac);
        intree->SetBranchAddress("j9_chadfrac", &t->j9_chadfrac);
        intree->SetBranchAddress("j9_nhadfrac", &t->j9_nhadfrac);
        intree->SetBranchAddress("j9_phofrac", &t->j9_phofrac);
        intree->SetBranchAddress("j9_mufrac", &t->j9_mufrac);
        intree->SetBranchAddress("j9_elefrac", &t->j9_elefrac);
        intree->SetBranchAddress("j9_JECUnc", &t->j9_JECUnc);
        intree->SetBranchAddress("j9_leadTrackPt", &t->j9_leadTrackPt);
        intree->SetBranchAddress("j9_softLeptPt", &t->j9_softLeptPt);
        intree->SetBranchAddress("j9_softLeptPtRel", &t->j9_softLeptPtRel);
        intree->SetBranchAddress("j9_softLeptDR", &t->j9_softLeptDR);
        intree->SetBranchAddress("j9_softLeptIdLooseMu", &t->j9_softLeptIdLooseMu);
        intree->SetBranchAddress("j9_softLeptIdEle95", &t->j9_softLeptIdEle95);
    }
    intree->SetBranchAddress("j9_btagSF_M", &t->j9_btagSF_M);
    intree->SetBranchAddress("j9_flavour", &t->j9_flavour);
    intree->SetBranchAddress("j9_btagSFErrorUp_M", &t->j9_btagSFErrorUp_M);
    intree->SetBranchAddress("j9_btagSFErrorDown_M", &t->j9_btagSFErrorDown_M);
    intree->SetBranchAddress("j9_btagEff_M", &t->j9_btagEff_M);
    intree->SetBranchAddress("j9_btagEffError_M", &t->j9_btagEffError_M);
    intree->SetBranchAddress("j9_cutbased_wp_level", &t->j9_cutbased_wp_level);
    intree->SetBranchAddress("j9_simple_wp_level", &t->j9_simple_wp_level);
    intree->SetBranchAddress("j9_full_wp_level", &t->j9_full_wp_level);
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
    intree->SetBranchAddress("j10_secVtx3deL", &t->j10_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j10_secVtxM", &t->j10_secVtxM);
        intree->SetBranchAddress("j10_emfrac", &t->j10_emfrac);
        intree->SetBranchAddress("j10_hadfrac", &t->j10_hadfrac);
        intree->SetBranchAddress("j10_chadfrac", &t->j10_chadfrac);
        intree->SetBranchAddress("j10_nhadfrac", &t->j10_nhadfrac);
        intree->SetBranchAddress("j10_phofrac", &t->j10_phofrac);
        intree->SetBranchAddress("j10_mufrac", &t->j10_mufrac);
        intree->SetBranchAddress("j10_elefrac", &t->j10_elefrac);
        intree->SetBranchAddress("j10_JECUnc", &t->j10_JECUnc);
        intree->SetBranchAddress("j10_leadTrackPt", &t->j10_leadTrackPt);
        intree->SetBranchAddress("j10_softLeptPt", &t->j10_softLeptPt);
        intree->SetBranchAddress("j10_softLeptPtRel", &t->j10_softLeptPtRel);
        intree->SetBranchAddress("j10_softLeptDR", &t->j10_softLeptDR);
        intree->SetBranchAddress("j10_softLeptIdLooseMu", &t->j10_softLeptIdLooseMu);
        intree->SetBranchAddress("j10_softLeptIdEle95", &t->j10_softLeptIdEle95);
    }
    intree->SetBranchAddress("j10_btagSF_M", &t->j10_btagSF_M);
    intree->SetBranchAddress("j10_flavour", &t->j10_flavour);
    intree->SetBranchAddress("j10_btagSFErrorUp_M", &t->j10_btagSFErrorUp_M);
    intree->SetBranchAddress("j10_btagSFErrorDown_M", &t->j10_btagSFErrorDown_M);
    intree->SetBranchAddress("j10_btagEff_M", &t->j10_btagEff_M);
    intree->SetBranchAddress("j10_btagEffError_M", &t->j10_btagEffError_M);
    intree->SetBranchAddress("j10_cutbased_wp_level", &t->j10_cutbased_wp_level);
    intree->SetBranchAddress("j10_simple_wp_level", &t->j10_simple_wp_level);
    intree->SetBranchAddress("j10_full_wp_level", &t->j10_full_wp_level);
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
    intree->SetBranchAddress("j11_secVtx3deL", &t->j11_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j11_secVtxM", &t->j11_secVtxM);
        intree->SetBranchAddress("j11_emfrac", &t->j11_emfrac);
        intree->SetBranchAddress("j11_hadfrac", &t->j11_hadfrac);
        intree->SetBranchAddress("j11_chadfrac", &t->j11_chadfrac);
        intree->SetBranchAddress("j11_nhadfrac", &t->j11_nhadfrac);
        intree->SetBranchAddress("j11_phofrac", &t->j11_phofrac);
        intree->SetBranchAddress("j11_mufrac", &t->j11_mufrac);
        intree->SetBranchAddress("j11_elefrac", &t->j11_elefrac);
        intree->SetBranchAddress("j11_JECUnc", &t->j11_JECUnc);
        intree->SetBranchAddress("j11_leadTrackPt", &t->j11_leadTrackPt);
        intree->SetBranchAddress("j11_softLeptPt", &t->j11_softLeptPt);
        intree->SetBranchAddress("j11_softLeptPtRel", &t->j11_softLeptPtRel);
        intree->SetBranchAddress("j11_softLeptDR", &t->j11_softLeptDR);
        intree->SetBranchAddress("j11_softLeptIdLooseMu", &t->j11_softLeptIdLooseMu);
        intree->SetBranchAddress("j11_softLeptIdEle95", &t->j11_softLeptIdEle95);
    }
    intree->SetBranchAddress("j11_btagSF_M", &t->j11_btagSF_M);
    intree->SetBranchAddress("j11_flavour", &t->j11_flavour);
    intree->SetBranchAddress("j11_btagSFErrorUp_M", &t->j11_btagSFErrorUp_M);
    intree->SetBranchAddress("j11_btagSFErrorDown_M", &t->j11_btagSFErrorDown_M);
    intree->SetBranchAddress("j11_btagEff_M", &t->j11_btagEff_M);
    intree->SetBranchAddress("j11_btagEffError_M", &t->j11_btagEffError_M);
    intree->SetBranchAddress("j11_cutbased_wp_level", &t->j11_cutbased_wp_level);
    intree->SetBranchAddress("j11_simple_wp_level", &t->j11_simple_wp_level);
    intree->SetBranchAddress("j11_full_wp_level", &t->j11_full_wp_level);
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
    intree->SetBranchAddress("j12_secVtx3deL", &t->j12_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j12_secVtxM", &t->j12_secVtxM);
        intree->SetBranchAddress("j12_emfrac", &t->j12_emfrac);
        intree->SetBranchAddress("j12_hadfrac", &t->j12_hadfrac);
        intree->SetBranchAddress("j12_chadfrac", &t->j12_chadfrac);
        intree->SetBranchAddress("j12_nhadfrac", &t->j12_nhadfrac);
        intree->SetBranchAddress("j12_phofrac", &t->j12_phofrac);
        intree->SetBranchAddress("j12_mufrac", &t->j12_mufrac);
        intree->SetBranchAddress("j12_elefrac", &t->j12_elefrac);
        intree->SetBranchAddress("j12_JECUnc", &t->j12_JECUnc);
        intree->SetBranchAddress("j12_leadTrackPt", &t->j12_leadTrackPt);
        intree->SetBranchAddress("j12_softLeptPt", &t->j12_softLeptPt);
        intree->SetBranchAddress("j12_softLeptPtRel", &t->j12_softLeptPtRel);
        intree->SetBranchAddress("j12_softLeptDR", &t->j12_softLeptDR);
        intree->SetBranchAddress("j12_softLeptIdLooseMu", &t->j12_softLeptIdLooseMu);
        intree->SetBranchAddress("j12_softLeptIdEle95", &t->j12_softLeptIdEle95);
    }
    intree->SetBranchAddress("j12_btagSF_M", &t->j12_btagSF_M);
    intree->SetBranchAddress("j12_flavour", &t->j12_flavour);
    intree->SetBranchAddress("j12_btagSFErrorUp_M", &t->j12_btagSFErrorUp_M);
    intree->SetBranchAddress("j12_btagSFErrorDown_M", &t->j12_btagSFErrorDown_M);
    intree->SetBranchAddress("j12_btagEff_M", &t->j12_btagEff_M);
    intree->SetBranchAddress("j12_btagEffError_M", &t->j12_btagEffError_M);
    intree->SetBranchAddress("j12_cutbased_wp_level", &t->j12_cutbased_wp_level);
    intree->SetBranchAddress("j12_simple_wp_level", &t->j12_simple_wp_level);
    intree->SetBranchAddress("j12_full_wp_level", &t->j12_full_wp_level);
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
    intree->SetBranchAddress("j13_secVtx3deL", &t->j13_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j13_secVtxM", &t->j13_secVtxM);
        intree->SetBranchAddress("j13_emfrac", &t->j13_emfrac);
        intree->SetBranchAddress("j13_hadfrac", &t->j13_hadfrac);
        intree->SetBranchAddress("j13_chadfrac", &t->j13_chadfrac);
        intree->SetBranchAddress("j13_nhadfrac", &t->j13_nhadfrac);
        intree->SetBranchAddress("j13_phofrac", &t->j13_phofrac);
        intree->SetBranchAddress("j13_mufrac", &t->j13_mufrac);
        intree->SetBranchAddress("j13_elefrac", &t->j13_elefrac);
        intree->SetBranchAddress("j13_JECUnc", &t->j13_JECUnc);
        intree->SetBranchAddress("j13_leadTrackPt", &t->j13_leadTrackPt);
        intree->SetBranchAddress("j13_softLeptPt", &t->j13_softLeptPt);
        intree->SetBranchAddress("j13_softLeptPtRel", &t->j13_softLeptPtRel);
        intree->SetBranchAddress("j13_softLeptDR", &t->j13_softLeptDR);
        intree->SetBranchAddress("j13_softLeptIdLooseMu", &t->j13_softLeptIdLooseMu);
        intree->SetBranchAddress("j13_softLeptIdEle95", &t->j13_softLeptIdEle95);
    }
    intree->SetBranchAddress("j13_btagSF_M", &t->j13_btagSF_M);
    intree->SetBranchAddress("j13_flavour", &t->j13_flavour);
    intree->SetBranchAddress("j13_btagSFErrorUp_M", &t->j13_btagSFErrorUp_M);
    intree->SetBranchAddress("j13_btagSFErrorDown_M", &t->j13_btagSFErrorDown_M);
    intree->SetBranchAddress("j13_btagEff_M", &t->j13_btagEff_M);
    intree->SetBranchAddress("j13_btagEffError_M", &t->j13_btagEffError_M);
    intree->SetBranchAddress("j13_cutbased_wp_level", &t->j13_cutbased_wp_level);
    intree->SetBranchAddress("j13_simple_wp_level", &t->j13_simple_wp_level);
    intree->SetBranchAddress("j13_full_wp_level", &t->j13_full_wp_level);
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
    intree->SetBranchAddress("j14_secVtx3deL", &t->j14_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j14_secVtxM", &t->j14_secVtxM);
        intree->SetBranchAddress("j14_emfrac", &t->j14_emfrac);
        intree->SetBranchAddress("j14_hadfrac", &t->j14_hadfrac);
        intree->SetBranchAddress("j14_chadfrac", &t->j14_chadfrac);
        intree->SetBranchAddress("j14_nhadfrac", &t->j14_nhadfrac);
        intree->SetBranchAddress("j14_phofrac", &t->j14_phofrac);
        intree->SetBranchAddress("j14_mufrac", &t->j14_mufrac);
        intree->SetBranchAddress("j14_elefrac", &t->j14_elefrac);
        intree->SetBranchAddress("j14_JECUnc", &t->j14_JECUnc);
        intree->SetBranchAddress("j14_leadTrackPt", &t->j14_leadTrackPt);
        intree->SetBranchAddress("j14_softLeptPt", &t->j14_softLeptPt);
        intree->SetBranchAddress("j14_softLeptPtRel", &t->j14_softLeptPtRel);
        intree->SetBranchAddress("j14_softLeptDR", &t->j14_softLeptDR);
        intree->SetBranchAddress("j14_softLeptIdLooseMu", &t->j14_softLeptIdLooseMu);
        intree->SetBranchAddress("j14_softLeptIdEle95", &t->j14_softLeptIdEle95);
    }
    intree->SetBranchAddress("j14_btagSF_M", &t->j14_btagSF_M);
    intree->SetBranchAddress("j14_flavour", &t->j14_flavour);
    intree->SetBranchAddress("j14_btagSFErrorUp_M", &t->j14_btagSFErrorUp_M);
    intree->SetBranchAddress("j14_btagSFErrorDown_M", &t->j14_btagSFErrorDown_M);
    intree->SetBranchAddress("j14_btagEff_M", &t->j14_btagEff_M);
    intree->SetBranchAddress("j14_btagEffError_M", &t->j14_btagEffError_M);
    intree->SetBranchAddress("j14_cutbased_wp_level", &t->j14_cutbased_wp_level);
    intree->SetBranchAddress("j14_simple_wp_level", &t->j14_simple_wp_level);
    intree->SetBranchAddress("j14_full_wp_level", &t->j14_full_wp_level);
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
    intree->SetBranchAddress("j15_secVtx3deL", &t->j15_secVtx3deL);
    if( numberOfRegressionFiles > 0 )
    {
        intree->SetBranchAddress("j15_secVtxM", &t->j15_secVtxM);
        intree->SetBranchAddress("j15_emfrac", &t->j15_emfrac);
        intree->SetBranchAddress("j15_hadfrac", &t->j15_hadfrac);
        intree->SetBranchAddress("j15_chadfrac", &t->j15_chadfrac);
        intree->SetBranchAddress("j15_nhadfrac", &t->j15_nhadfrac);
        intree->SetBranchAddress("j15_phofrac", &t->j15_phofrac);
        intree->SetBranchAddress("j15_mufrac", &t->j15_mufrac);
        intree->SetBranchAddress("j15_elefrac", &t->j15_elefrac);
        intree->SetBranchAddress("j15_JECUnc", &t->j15_JECUnc);
        intree->SetBranchAddress("j15_leadTrackPt", &t->j15_leadTrackPt);
        intree->SetBranchAddress("j15_softLeptPt", &t->j15_softLeptPt);
        intree->SetBranchAddress("j15_softLeptPtRel", &t->j15_softLeptPtRel);
        intree->SetBranchAddress("j15_softLeptDR", &t->j15_softLeptDR);
        intree->SetBranchAddress("j15_softLeptIdLooseMu", &t->j15_softLeptIdLooseMu);
        intree->SetBranchAddress("j15_softLeptIdEle95", &t->j15_softLeptIdEle95);
    }
    intree->SetBranchAddress("j15_btagSF_M", &t->j15_btagSF_M);
    intree->SetBranchAddress("j15_flavour", &t->j15_flavour);
    intree->SetBranchAddress("j15_btagSFErrorUp_M", &t->j15_btagSFErrorUp_M);
    intree->SetBranchAddress("j15_btagSFErrorDown_M", &t->j15_btagSFErrorDown_M);
    intree->SetBranchAddress("j15_btagEff_M", &t->j15_btagEff_M);
    intree->SetBranchAddress("j15_btagEffError_M", &t->j15_btagEffError_M);
    intree->SetBranchAddress("j15_cutbased_wp_level", &t->j15_cutbased_wp_level);
    intree->SetBranchAddress("j15_simple_wp_level", &t->j15_simple_wp_level);
    intree->SetBranchAddress("j15_full_wp_level", &t->j15_full_wp_level);
    intree->SetBranchAddress("j15_ntk", &t->j15_ntk);
    intree->SetBranchAddress("j15_nNeutrals", &t->j15_nNeutrals);
    intree->SetBranchAddress("j15_nCharged", &t->j15_nCharged);
    intree->SetBranchAddress("j15_axis1", &t->j15_axis1);
    intree->SetBranchAddress("j15_axis2", &t->j15_axis2);
    intree->SetBranchAddress("j15_pull", &t->j15_pull);

    // Jet Energy Correction and Jet Energy Resolution
    intree->SetBranchAddress("j1_jecD_e", &t->j1_jecD_e);
    intree->SetBranchAddress("j1_jecD_pt", &t->j1_jecD_pt);
    intree->SetBranchAddress("j1_jecD_phi", &t->j1_jecD_phi);
    intree->SetBranchAddress("j1_jecD_eta", &t->j1_jecD_eta);
    intree->SetBranchAddress("j1_jecU_e", &t->j1_jecU_e);
    intree->SetBranchAddress("j1_jecU_pt", &t->j1_jecU_pt);
    intree->SetBranchAddress("j1_jecU_phi", &t->j1_jecU_phi);
    intree->SetBranchAddress("j1_jecU_eta", &t->j1_jecU_eta);
    intree->SetBranchAddress("j1_jerD_e", &t->j1_jerD_e);
    intree->SetBranchAddress("j1_jerD_pt", &t->j1_jerD_pt);
    intree->SetBranchAddress("j1_jerD_phi", &t->j1_jerD_phi);
    intree->SetBranchAddress("j1_jerD_eta", &t->j1_jerD_eta);
    intree->SetBranchAddress("j1_jerC_e", &t->j1_jerC_e);
    intree->SetBranchAddress("j1_jerC_pt", &t->j1_jerC_pt);
    intree->SetBranchAddress("j1_jerC_phi", &t->j1_jerC_phi);
    intree->SetBranchAddress("j1_jerC_eta", &t->j1_jerC_eta);
    intree->SetBranchAddress("j1_jerU_e", &t->j1_jerU_e);
    intree->SetBranchAddress("j1_jerU_pt", &t->j1_jerU_pt);
    intree->SetBranchAddress("j1_jerU_phi", &t->j1_jerU_phi);
    intree->SetBranchAddress("j1_jerU_eta", &t->j1_jerU_eta);

    intree->SetBranchAddress("j2_jecD_e", &t->j2_jecD_e);
    intree->SetBranchAddress("j2_jecD_pt", &t->j2_jecD_pt);
    intree->SetBranchAddress("j2_jecD_phi", &t->j2_jecD_phi);
    intree->SetBranchAddress("j2_jecD_eta", &t->j2_jecD_eta);
    intree->SetBranchAddress("j2_jecU_e", &t->j2_jecU_e);
    intree->SetBranchAddress("j2_jecU_pt", &t->j2_jecU_pt);
    intree->SetBranchAddress("j2_jecU_phi", &t->j2_jecU_phi);
    intree->SetBranchAddress("j2_jecU_eta", &t->j2_jecU_eta);
    intree->SetBranchAddress("j2_jerD_e", &t->j2_jerD_e);
    intree->SetBranchAddress("j2_jerD_pt", &t->j2_jerD_pt);
    intree->SetBranchAddress("j2_jerD_phi", &t->j2_jerD_phi);
    intree->SetBranchAddress("j2_jerD_eta", &t->j2_jerD_eta);
    intree->SetBranchAddress("j2_jerC_e", &t->j2_jerC_e);
    intree->SetBranchAddress("j2_jerC_pt", &t->j2_jerC_pt);
    intree->SetBranchAddress("j2_jerC_phi", &t->j2_jerC_phi);
    intree->SetBranchAddress("j2_jerC_eta", &t->j2_jerC_eta);
    intree->SetBranchAddress("j2_jerU_e", &t->j2_jerU_e);
    intree->SetBranchAddress("j2_jerU_pt", &t->j2_jerU_pt);
    intree->SetBranchAddress("j2_jerU_phi", &t->j2_jerU_phi);
    intree->SetBranchAddress("j2_jerU_eta", &t->j2_jerU_eta);

    intree->SetBranchAddress("j3_jecD_e", &t->j3_jecD_e);
    intree->SetBranchAddress("j3_jecD_pt", &t->j3_jecD_pt);
    intree->SetBranchAddress("j3_jecD_phi", &t->j3_jecD_phi);
    intree->SetBranchAddress("j3_jecD_eta", &t->j3_jecD_eta);
    intree->SetBranchAddress("j3_jecU_e", &t->j3_jecU_e);
    intree->SetBranchAddress("j3_jecU_pt", &t->j3_jecU_pt);
    intree->SetBranchAddress("j3_jecU_phi", &t->j3_jecU_phi);
    intree->SetBranchAddress("j3_jecU_eta", &t->j3_jecU_eta);
    intree->SetBranchAddress("j3_jerD_e", &t->j3_jerD_e);
    intree->SetBranchAddress("j3_jerD_pt", &t->j3_jerD_pt);
    intree->SetBranchAddress("j3_jerD_phi", &t->j3_jerD_phi);
    intree->SetBranchAddress("j3_jerD_eta", &t->j3_jerD_eta);
    intree->SetBranchAddress("j3_jerC_e", &t->j3_jerC_e);
    intree->SetBranchAddress("j3_jerC_pt", &t->j3_jerC_pt);
    intree->SetBranchAddress("j3_jerC_phi", &t->j3_jerC_phi);
    intree->SetBranchAddress("j3_jerC_eta", &t->j3_jerC_eta);
    intree->SetBranchAddress("j3_jerU_e", &t->j3_jerU_e);
    intree->SetBranchAddress("j3_jerU_pt", &t->j3_jerU_pt);
    intree->SetBranchAddress("j3_jerU_phi", &t->j3_jerU_phi);
    intree->SetBranchAddress("j3_jerU_eta", &t->j3_jerU_eta);

    intree->SetBranchAddress("j4_jecD_e", &t->j4_jecD_e);
    intree->SetBranchAddress("j4_jecD_pt", &t->j4_jecD_pt);
    intree->SetBranchAddress("j4_jecD_phi", &t->j4_jecD_phi);
    intree->SetBranchAddress("j4_jecD_eta", &t->j4_jecD_eta);
    intree->SetBranchAddress("j4_jecU_e", &t->j4_jecU_e);
    intree->SetBranchAddress("j4_jecU_pt", &t->j4_jecU_pt);
    intree->SetBranchAddress("j4_jecU_phi", &t->j4_jecU_phi);
    intree->SetBranchAddress("j4_jecU_eta", &t->j4_jecU_eta);
    intree->SetBranchAddress("j4_jerD_e", &t->j4_jerD_e);
    intree->SetBranchAddress("j4_jerD_pt", &t->j4_jerD_pt);
    intree->SetBranchAddress("j4_jerD_phi", &t->j4_jerD_phi);
    intree->SetBranchAddress("j4_jerD_eta", &t->j4_jerD_eta);
    intree->SetBranchAddress("j4_jerC_e", &t->j4_jerC_e);
    intree->SetBranchAddress("j4_jerC_pt", &t->j4_jerC_pt);
    intree->SetBranchAddress("j4_jerC_phi", &t->j4_jerC_phi);
    intree->SetBranchAddress("j4_jerC_eta", &t->j4_jerC_eta);
    intree->SetBranchAddress("j4_jerU_e", &t->j4_jerU_e);
    intree->SetBranchAddress("j4_jerU_pt", &t->j4_jerU_pt);
    intree->SetBranchAddress("j4_jerU_phi", &t->j4_jerU_phi);
    intree->SetBranchAddress("j4_jerU_eta", &t->j4_jerU_eta);

    intree->SetBranchAddress("j5_jecD_e", &t->j5_jecD_e);
    intree->SetBranchAddress("j5_jecD_pt", &t->j5_jecD_pt);
    intree->SetBranchAddress("j5_jecD_phi", &t->j5_jecD_phi);
    intree->SetBranchAddress("j5_jecD_eta", &t->j5_jecD_eta);
    intree->SetBranchAddress("j5_jecU_e", &t->j5_jecU_e);
    intree->SetBranchAddress("j5_jecU_pt", &t->j5_jecU_pt);
    intree->SetBranchAddress("j5_jecU_phi", &t->j5_jecU_phi);
    intree->SetBranchAddress("j5_jecU_eta", &t->j5_jecU_eta);
    intree->SetBranchAddress("j5_jerD_e", &t->j5_jerD_e);
    intree->SetBranchAddress("j5_jerD_pt", &t->j5_jerD_pt);
    intree->SetBranchAddress("j5_jerD_phi", &t->j5_jerD_phi);
    intree->SetBranchAddress("j5_jerD_eta", &t->j5_jerD_eta);
    intree->SetBranchAddress("j5_jerC_e", &t->j5_jerC_e);
    intree->SetBranchAddress("j5_jerC_pt", &t->j5_jerC_pt);
    intree->SetBranchAddress("j5_jerC_phi", &t->j5_jerC_phi);
    intree->SetBranchAddress("j5_jerC_eta", &t->j5_jerC_eta);
    intree->SetBranchAddress("j5_jerU_e", &t->j5_jerU_e);
    intree->SetBranchAddress("j5_jerU_pt", &t->j5_jerU_pt);
    intree->SetBranchAddress("j5_jerU_phi", &t->j5_jerU_phi);
    intree->SetBranchAddress("j5_jerU_eta", &t->j5_jerU_eta);

    intree->SetBranchAddress("j6_jecD_e", &t->j6_jecD_e);
    intree->SetBranchAddress("j6_jecD_pt", &t->j6_jecD_pt);
    intree->SetBranchAddress("j6_jecD_phi", &t->j6_jecD_phi);
    intree->SetBranchAddress("j6_jecD_eta", &t->j6_jecD_eta);
    intree->SetBranchAddress("j6_jecU_e", &t->j6_jecU_e);
    intree->SetBranchAddress("j6_jecU_pt", &t->j6_jecU_pt);
    intree->SetBranchAddress("j6_jecU_phi", &t->j6_jecU_phi);
    intree->SetBranchAddress("j6_jecU_eta", &t->j6_jecU_eta);
    intree->SetBranchAddress("j6_jerD_e", &t->j6_jerD_e);
    intree->SetBranchAddress("j6_jerD_pt", &t->j6_jerD_pt);
    intree->SetBranchAddress("j6_jerD_phi", &t->j6_jerD_phi);
    intree->SetBranchAddress("j6_jerD_eta", &t->j6_jerD_eta);
    intree->SetBranchAddress("j6_jerC_e", &t->j6_jerC_e);
    intree->SetBranchAddress("j6_jerC_pt", &t->j6_jerC_pt);
    intree->SetBranchAddress("j6_jerC_phi", &t->j6_jerC_phi);
    intree->SetBranchAddress("j6_jerC_eta", &t->j6_jerC_eta);
    intree->SetBranchAddress("j6_jerU_e", &t->j6_jerU_e);
    intree->SetBranchAddress("j6_jerU_pt", &t->j6_jerU_pt);
    intree->SetBranchAddress("j6_jerU_phi", &t->j6_jerU_phi);
    intree->SetBranchAddress("j6_jerU_eta", &t->j6_jerU_eta);

    intree->SetBranchAddress("j7_jecD_e", &t->j7_jecD_e);
    intree->SetBranchAddress("j7_jecD_pt", &t->j7_jecD_pt);
    intree->SetBranchAddress("j7_jecD_phi", &t->j7_jecD_phi);
    intree->SetBranchAddress("j7_jecD_eta", &t->j7_jecD_eta);
    intree->SetBranchAddress("j7_jecU_e", &t->j7_jecU_e);
    intree->SetBranchAddress("j7_jecU_pt", &t->j7_jecU_pt);
    intree->SetBranchAddress("j7_jecU_phi", &t->j7_jecU_phi);
    intree->SetBranchAddress("j7_jecU_eta", &t->j7_jecU_eta);
    intree->SetBranchAddress("j7_jerD_e", &t->j7_jerD_e);
    intree->SetBranchAddress("j7_jerD_pt", &t->j7_jerD_pt);
    intree->SetBranchAddress("j7_jerD_phi", &t->j7_jerD_phi);
    intree->SetBranchAddress("j7_jerD_eta", &t->j7_jerD_eta);
    intree->SetBranchAddress("j7_jerC_e", &t->j7_jerC_e);
    intree->SetBranchAddress("j7_jerC_pt", &t->j7_jerC_pt);
    intree->SetBranchAddress("j7_jerC_phi", &t->j7_jerC_phi);
    intree->SetBranchAddress("j7_jerC_eta", &t->j7_jerC_eta);
    intree->SetBranchAddress("j7_jerU_e", &t->j7_jerU_e);
    intree->SetBranchAddress("j7_jerU_pt", &t->j7_jerU_pt);
    intree->SetBranchAddress("j7_jerU_phi", &t->j7_jerU_phi);
    intree->SetBranchAddress("j7_jerU_eta", &t->j7_jerU_eta);

    intree->SetBranchAddress("j8_jecD_e", &t->j8_jecD_e);
    intree->SetBranchAddress("j8_jecD_pt", &t->j8_jecD_pt);
    intree->SetBranchAddress("j8_jecD_phi", &t->j8_jecD_phi);
    intree->SetBranchAddress("j8_jecD_eta", &t->j8_jecD_eta);
    intree->SetBranchAddress("j8_jecU_e", &t->j8_jecU_e);
    intree->SetBranchAddress("j8_jecU_pt", &t->j8_jecU_pt);
    intree->SetBranchAddress("j8_jecU_phi", &t->j8_jecU_phi);
    intree->SetBranchAddress("j8_jecU_eta", &t->j8_jecU_eta);
    intree->SetBranchAddress("j8_jerD_e", &t->j8_jerD_e);
    intree->SetBranchAddress("j8_jerD_pt", &t->j8_jerD_pt);
    intree->SetBranchAddress("j8_jerD_phi", &t->j8_jerD_phi);
    intree->SetBranchAddress("j8_jerD_eta", &t->j8_jerD_eta);
    intree->SetBranchAddress("j8_jerC_e", &t->j8_jerC_e);
    intree->SetBranchAddress("j8_jerC_pt", &t->j8_jerC_pt);
    intree->SetBranchAddress("j8_jerC_phi", &t->j8_jerC_phi);
    intree->SetBranchAddress("j8_jerC_eta", &t->j8_jerC_eta);
    intree->SetBranchAddress("j8_jerU_e", &t->j8_jerU_e);
    intree->SetBranchAddress("j8_jerU_pt", &t->j8_jerU_pt);
    intree->SetBranchAddress("j8_jerU_phi", &t->j8_jerU_phi);
    intree->SetBranchAddress("j8_jerU_eta", &t->j8_jerU_eta);

    intree->SetBranchAddress("j9_jecD_e", &t->j9_jecD_e);
    intree->SetBranchAddress("j9_jecD_pt", &t->j9_jecD_pt);
    intree->SetBranchAddress("j9_jecD_phi", &t->j9_jecD_phi);
    intree->SetBranchAddress("j9_jecD_eta", &t->j9_jecD_eta);
    intree->SetBranchAddress("j9_jecU_e", &t->j9_jecU_e);
    intree->SetBranchAddress("j9_jecU_pt", &t->j9_jecU_pt);
    intree->SetBranchAddress("j9_jecU_phi", &t->j9_jecU_phi);
    intree->SetBranchAddress("j9_jecU_eta", &t->j9_jecU_eta);
    intree->SetBranchAddress("j9_jerD_e", &t->j9_jerD_e);
    intree->SetBranchAddress("j9_jerD_pt", &t->j9_jerD_pt);
    intree->SetBranchAddress("j9_jerD_phi", &t->j9_jerD_phi);
    intree->SetBranchAddress("j9_jerD_eta", &t->j9_jerD_eta);
    intree->SetBranchAddress("j9_jerC_e", &t->j9_jerC_e);
    intree->SetBranchAddress("j9_jerC_pt", &t->j9_jerC_pt);
    intree->SetBranchAddress("j9_jerC_phi", &t->j9_jerC_phi);
    intree->SetBranchAddress("j9_jerC_eta", &t->j9_jerC_eta);
    intree->SetBranchAddress("j9_jerU_e", &t->j9_jerU_e);
    intree->SetBranchAddress("j9_jerU_pt", &t->j9_jerU_pt);
    intree->SetBranchAddress("j9_jerU_phi", &t->j9_jerU_phi);
    intree->SetBranchAddress("j9_jerU_eta", &t->j9_jerU_eta);

    intree->SetBranchAddress("j10_jecD_e", &t->j10_jecD_e);
    intree->SetBranchAddress("j10_jecD_pt", &t->j10_jecD_pt);
    intree->SetBranchAddress("j10_jecD_phi", &t->j10_jecD_phi);
    intree->SetBranchAddress("j10_jecD_eta", &t->j10_jecD_eta);
    intree->SetBranchAddress("j10_jecU_e", &t->j10_jecU_e);
    intree->SetBranchAddress("j10_jecU_pt", &t->j10_jecU_pt);
    intree->SetBranchAddress("j10_jecU_phi", &t->j10_jecU_phi);
    intree->SetBranchAddress("j10_jecU_eta", &t->j10_jecU_eta);
    intree->SetBranchAddress("j10_jerD_e", &t->j10_jerD_e);
    intree->SetBranchAddress("j10_jerD_pt", &t->j10_jerD_pt);
    intree->SetBranchAddress("j10_jerD_phi", &t->j10_jerD_phi);
    intree->SetBranchAddress("j10_jerD_eta", &t->j10_jerD_eta);
    intree->SetBranchAddress("j10_jerC_e", &t->j10_jerC_e);
    intree->SetBranchAddress("j10_jerC_pt", &t->j10_jerC_pt);
    intree->SetBranchAddress("j10_jerC_phi", &t->j10_jerC_phi);
    intree->SetBranchAddress("j10_jerC_eta", &t->j10_jerC_eta);
    intree->SetBranchAddress("j10_jerU_e", &t->j10_jerU_e);
    intree->SetBranchAddress("j10_jerU_pt", &t->j10_jerU_pt);
    intree->SetBranchAddress("j10_jerU_phi", &t->j10_jerU_phi);
    intree->SetBranchAddress("j10_jerU_eta", &t->j10_jerU_eta);

    intree->SetBranchAddress("j11_jecD_e", &t->j11_jecD_e);
    intree->SetBranchAddress("j11_jecD_pt", &t->j11_jecD_pt);
    intree->SetBranchAddress("j11_jecD_phi", &t->j11_jecD_phi);
    intree->SetBranchAddress("j11_jecD_eta", &t->j11_jecD_eta);
    intree->SetBranchAddress("j11_jecU_e", &t->j11_jecU_e);
    intree->SetBranchAddress("j11_jecU_pt", &t->j11_jecU_pt);
    intree->SetBranchAddress("j11_jecU_phi", &t->j11_jecU_phi);
    intree->SetBranchAddress("j11_jecU_eta", &t->j11_jecU_eta);
    intree->SetBranchAddress("j11_jerD_e", &t->j11_jerD_e);
    intree->SetBranchAddress("j11_jerD_pt", &t->j11_jerD_pt);
    intree->SetBranchAddress("j11_jerD_phi", &t->j11_jerD_phi);
    intree->SetBranchAddress("j11_jerD_eta", &t->j11_jerD_eta);
    intree->SetBranchAddress("j11_jerC_e", &t->j11_jerC_e);
    intree->SetBranchAddress("j11_jerC_pt", &t->j11_jerC_pt);
    intree->SetBranchAddress("j11_jerC_phi", &t->j11_jerC_phi);
    intree->SetBranchAddress("j11_jerC_eta", &t->j11_jerC_eta);
    intree->SetBranchAddress("j11_jerU_e", &t->j11_jerU_e);
    intree->SetBranchAddress("j11_jerU_pt", &t->j11_jerU_pt);
    intree->SetBranchAddress("j11_jerU_phi", &t->j11_jerU_phi);
    intree->SetBranchAddress("j11_jerU_eta", &t->j11_jerU_eta);

    intree->SetBranchAddress("j12_jecD_e", &t->j12_jecD_e);
    intree->SetBranchAddress("j12_jecD_pt", &t->j12_jecD_pt);
    intree->SetBranchAddress("j12_jecD_phi", &t->j12_jecD_phi);
    intree->SetBranchAddress("j12_jecD_eta", &t->j12_jecD_eta);
    intree->SetBranchAddress("j12_jecU_e", &t->j12_jecU_e);
    intree->SetBranchAddress("j12_jecU_pt", &t->j12_jecU_pt);
    intree->SetBranchAddress("j12_jecU_phi", &t->j12_jecU_phi);
    intree->SetBranchAddress("j12_jecU_eta", &t->j12_jecU_eta);
    intree->SetBranchAddress("j12_jerD_e", &t->j12_jerD_e);
    intree->SetBranchAddress("j12_jerD_pt", &t->j12_jerD_pt);
    intree->SetBranchAddress("j12_jerD_phi", &t->j12_jerD_phi);
    intree->SetBranchAddress("j12_jerD_eta", &t->j12_jerD_eta);
    intree->SetBranchAddress("j12_jerC_e", &t->j12_jerC_e);
    intree->SetBranchAddress("j12_jerC_pt", &t->j12_jerC_pt);
    intree->SetBranchAddress("j12_jerC_phi", &t->j12_jerC_phi);
    intree->SetBranchAddress("j12_jerC_eta", &t->j12_jerC_eta);
    intree->SetBranchAddress("j12_jerU_e", &t->j12_jerU_e);
    intree->SetBranchAddress("j12_jerU_pt", &t->j12_jerU_pt);
    intree->SetBranchAddress("j12_jerU_phi", &t->j12_jerU_phi);
    intree->SetBranchAddress("j12_jerU_eta", &t->j12_jerU_eta);

    intree->SetBranchAddress("j13_jecD_e", &t->j13_jecD_e);
    intree->SetBranchAddress("j13_jecD_pt", &t->j13_jecD_pt);
    intree->SetBranchAddress("j13_jecD_phi", &t->j13_jecD_phi);
    intree->SetBranchAddress("j13_jecD_eta", &t->j13_jecD_eta);
    intree->SetBranchAddress("j13_jecU_e", &t->j13_jecU_e);
    intree->SetBranchAddress("j13_jecU_pt", &t->j13_jecU_pt);
    intree->SetBranchAddress("j13_jecU_phi", &t->j13_jecU_phi);
    intree->SetBranchAddress("j13_jecU_eta", &t->j13_jecU_eta);
    intree->SetBranchAddress("j13_jerD_e", &t->j13_jerD_e);
    intree->SetBranchAddress("j13_jerD_pt", &t->j13_jerD_pt);
    intree->SetBranchAddress("j13_jerD_phi", &t->j13_jerD_phi);
    intree->SetBranchAddress("j13_jerD_eta", &t->j13_jerD_eta);
    intree->SetBranchAddress("j13_jerC_e", &t->j13_jerC_e);
    intree->SetBranchAddress("j13_jerC_pt", &t->j13_jerC_pt);
    intree->SetBranchAddress("j13_jerC_phi", &t->j13_jerC_phi);
    intree->SetBranchAddress("j13_jerC_eta", &t->j13_jerC_eta);
    intree->SetBranchAddress("j13_jerU_e", &t->j13_jerU_e);
    intree->SetBranchAddress("j13_jerU_pt", &t->j13_jerU_pt);
    intree->SetBranchAddress("j13_jerU_phi", &t->j13_jerU_phi);
    intree->SetBranchAddress("j13_jerU_eta", &t->j13_jerU_eta);

    intree->SetBranchAddress("j14_jecD_e", &t->j14_jecD_e);
    intree->SetBranchAddress("j14_jecD_pt", &t->j14_jecD_pt);
    intree->SetBranchAddress("j14_jecD_phi", &t->j14_jecD_phi);
    intree->SetBranchAddress("j14_jecD_eta", &t->j14_jecD_eta);
    intree->SetBranchAddress("j14_jecU_e", &t->j14_jecU_e);
    intree->SetBranchAddress("j14_jecU_pt", &t->j14_jecU_pt);
    intree->SetBranchAddress("j14_jecU_phi", &t->j14_jecU_phi);
    intree->SetBranchAddress("j14_jecU_eta", &t->j14_jecU_eta);
    intree->SetBranchAddress("j14_jerD_e", &t->j14_jerD_e);
    intree->SetBranchAddress("j14_jerD_pt", &t->j14_jerD_pt);
    intree->SetBranchAddress("j14_jerD_phi", &t->j14_jerD_phi);
    intree->SetBranchAddress("j14_jerD_eta", &t->j14_jerD_eta);
    intree->SetBranchAddress("j14_jerC_e", &t->j14_jerC_e);
    intree->SetBranchAddress("j14_jerC_pt", &t->j14_jerC_pt);
    intree->SetBranchAddress("j14_jerC_phi", &t->j14_jerC_phi);
    intree->SetBranchAddress("j14_jerC_eta", &t->j14_jerC_eta);
    intree->SetBranchAddress("j14_jerU_e", &t->j14_jerU_e);
    intree->SetBranchAddress("j14_jerU_pt", &t->j14_jerU_pt);
    intree->SetBranchAddress("j14_jerU_phi", &t->j14_jerU_phi);
    intree->SetBranchAddress("j14_jerU_eta", &t->j14_jerU_eta);

    intree->SetBranchAddress("j15_jecD_e", &t->j15_jecD_e);
    intree->SetBranchAddress("j15_jecD_pt", &t->j15_jecD_pt);
    intree->SetBranchAddress("j15_jecD_phi", &t->j15_jecD_phi);
    intree->SetBranchAddress("j15_jecD_eta", &t->j15_jecD_eta);
    intree->SetBranchAddress("j15_jecU_e", &t->j15_jecU_e);
    intree->SetBranchAddress("j15_jecU_pt", &t->j15_jecU_pt);
    intree->SetBranchAddress("j15_jecU_phi", &t->j15_jecU_phi);
    intree->SetBranchAddress("j15_jecU_eta", &t->j15_jecU_eta);
    intree->SetBranchAddress("j15_jerD_e", &t->j15_jerD_e);
    intree->SetBranchAddress("j15_jerD_pt", &t->j15_jerD_pt);
    intree->SetBranchAddress("j15_jerD_phi", &t->j15_jerD_phi);
    intree->SetBranchAddress("j15_jerD_eta", &t->j15_jerD_eta);
    intree->SetBranchAddress("j15_jerC_e", &t->j15_jerC_e);
    intree->SetBranchAddress("j15_jerC_pt", &t->j15_jerC_pt);
    intree->SetBranchAddress("j15_jerC_phi", &t->j15_jerC_phi);
    intree->SetBranchAddress("j15_jerC_eta", &t->j15_jerC_eta);
    intree->SetBranchAddress("j15_jerU_e", &t->j15_jerU_e);
    intree->SetBranchAddress("j15_jerU_pt", &t->j15_jerU_pt);
    intree->SetBranchAddress("j15_jerU_phi", &t->j15_jerU_phi);
    intree->SetBranchAddress("j15_jerU_eta", &t->j15_jerU_eta);


    return;
}

void setup_outtree(TTree* outtree, tree_variables *t)
{
    outtree->Branch("category", &t->category, "category/I");
    outtree->Branch("selection_cut_level", &t->selection_cut_level, "selection_cut_level/I");
    outtree->Branch("run", &t->run, "run/I");
    outtree->Branch("lumis", &t->lumis, "lumis/I");
    outtree->Branch("event", &t->event, "event/I");
    outtree->Branch("vtx_z", &t->vtx_z, "vtx_z/F");
    outtree->Branch("weight", &t->weight, "weight/F");
    outtree->Branch("evweight", &t->evweight, "evweight/F");
    outtree->Branch("evweight_w_btagSF", &t->evweight_w_btagSF, "evweight_w_btagSF/F");
    outtree->Branch("evweight_w_btagSF_reg", &t->evweight_w_btagSF_reg, "evweight_w_btagSF_reg/F");
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
    outtree->Branch("pho1_r9_cic", &t->pho1_r9_cic, "pho1_r9_cic/F");
    outtree->Branch("pho1_IDmva", &t->pho1_IDmva, "pho1_IDmva/F");
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
    outtree->Branch("pho2_r9_cic", &t->pho2_r9_cic, "pho2_r9_cic/F");
    outtree->Branch("pho2_IDmva", &t->pho2_IDmva, "pho2_IDmva/F");
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
    outtree->Branch("pho1_sigmaEoE", &t->ph1_sigmaEoE, "pho1_sigmaEoE/F");
    outtree->Branch("pho2_sigmaEoE", &t->ph2_sigmaEoE, "pho2_sigmaEoE/F");
    outtree->Branch("pho1_pesD_e", &t->ph1_pesD_e, "pho1_pesD_e/F");
    outtree->Branch("pho2_pesD_e", &t->ph2_pesD_e, "pho2_pesD_e/F");
    outtree->Branch("pho1_pesU_e", &t->ph1_pesU_e, "pho1_pesU_e/F");
    outtree->Branch("pho2_pesU_e", &t->ph2_pesU_e, "pho2_pesU_e/F");
    outtree->Branch("pho1_perD_e", &t->ph1_perD_e, "pho1_perD_e/F");
    outtree->Branch("pho2_perD_e", &t->ph2_perD_e, "pho2_perD_e/F");
    outtree->Branch("pho1_perU_e", &t->ph1_perU_e, "pho1_perU_e/F");
    outtree->Branch("pho2_perU_e", &t->ph2_perU_e, "pho2_perU_e/F");
    outtree->Branch("jet1_pt", &t->jet1_pt, "jet1_pt/F");
    outtree->Branch("jet1_pt", &t->jet1_pt, "jet1_pt/F");
    outtree->Branch("jet1_e", &t->jet1_e, "jet1_e/F");
    outtree->Branch("jet1_phi", &t->jet1_phi, "jet1_phi/F");
    outtree->Branch("jet1_eta", &t->jet1_eta, "jet1_eta/F");
    outtree->Branch("jet1_mass", &t->jet1_mass, "jet1_mass/F");
    outtree->Branch("jet1_csvBtag", &t->jet1_csvBtag, "jet1_csvBtag/F");
    outtree->Branch("jet1_flavour", &t->jet1_flavour, "jet1_flavour/I");
    outtree->Branch("jet1_btagSF_M", &t->jet1_btagSF_M, "jet1_btagSF_M/F");
    outtree->Branch("jet1_btagSFErrorUp_M", &t->jet1_btagSFErrorUp_M, "jet1_btagSFErrorUp_M/F");
    outtree->Branch("jet1_btagSFErrorDown_M", &t->jet1_btagSFErrorDown_M, "jet1_btagSFErrorDown_M/F");
    outtree->Branch("jet1_btagEff_M", &t->jet1_btagEff_M, "jet1_btagEff_M/F");
    outtree->Branch("jet1_btagEffError_M", &t->jet1_btagEffError_M, "jet1_btagEffError_M/F");
    outtree->Branch("jet1_cutbased_wp_level", &t->jet1_cutbased_wp_level, "jet1_cutbased_wp_level/I");
    outtree->Branch("jet1_simple_wp_level", &t->jet1_simple_wp_level, "jet1_simple_wp_level/I");
    outtree->Branch("jet1_full_wp_level", &t->jet1_full_wp_level, "jet1_full_wp_level/I");
    outtree->Branch("jet1_betaStarClassic", &t->jet1_betaStarClassic, "jet1_betaStarClassic/F");
    outtree->Branch("jet1_dR2Mean", &t->jet1_dR2Mean, "jet1_dR2Mean/F");
    outtree->Branch("jet2_pt", &t->jet2_pt, "jet2_pt/F");
    outtree->Branch("jet2_e", &t->jet2_e, "jet2_e/F");
    outtree->Branch("jet2_phi", &t->jet2_phi, "jet2_phi/F");
    outtree->Branch("jet2_eta", &t->jet2_eta, "jet2_eta/F");
    outtree->Branch("jet2_mass", &t->jet2_mass, "jet2_mass/F");
    outtree->Branch("jet2_csvBtag", &t->jet2_csvBtag, "jet2_csvBtag/F");
    outtree->Branch("jet2_flavour", &t->jet2_flavour, "jet2_flavour/I");
    outtree->Branch("jet2_btagSF_M", &t->jet2_btagSF_M, "jet2_btagSF_M/F");
    outtree->Branch("jet2_btagSFErrorUp_M", &t->jet2_btagSFErrorUp_M, "jet2_btagSFErrorUp_M/F");
    outtree->Branch("jet2_btagSFErrorDown_M", &t->jet2_btagSFErrorDown_M, "jet2_btagSFErrorDown_M/F");
    outtree->Branch("jet2_btagEff_M", &t->jet2_btagEff_M, "jet2_btagEff_M/F");
    outtree->Branch("jet2_btagEffError_M", &t->jet2_btagEffError_M, "jet2_btagEffError_M/F");
    outtree->Branch("jet2_cutbased_wp_level", &t->jet2_cutbased_wp_level, "jet2_cutbased_wp_level/I");
    outtree->Branch("jet2_simple_wp_level", &t->jet2_simple_wp_level, "jet2_simple_wp_level/I");
    outtree->Branch("jet2_full_wp_level", &t->jet2_full_wp_level, "jet2_full_wp_level/I");
    outtree->Branch("jet2_betaStarClassic", &t->jet2_betaStarClassic, "jet2_betaStarClassic/F");
    outtree->Branch("jet2_dR2Mean", &t->jet2_dR2Mean, "jet2_dR2Mean/F");
// storing inputs of the regression for comparison
    outtree->Branch("regjet1_mt", &t->regjet1_mt, "regjet1_mt/F");
    outtree->Branch("regjet1_chadfrac", &t->regjet1_chadfrac, "regjet1_chadfrac/F");
    outtree->Branch("regjet1_nhadfrac", &t->regjet1_nhadfrac, "regjet1_nhadfrac/F");
    outtree->Branch("regjet1_phofrac", &t->regjet1_phofrac, "regjet1_phofrac/F");
    outtree->Branch("regjet1_mufrac", &t->regjet1_mufrac, "regjet1_mufrac/F");
    outtree->Branch("regjet1_elefrac", &t->regjet1_elefrac, "regjet1_elefrac/F");
    outtree->Branch("regjet1_softLeptPt", &t->regjet1_softLeptPt, "regjet1_softLeptPt/F");
    outtree->Branch("regjet1_softLeptPtRel", &t->regjet1_softLeptPtRel, "regjet1_softLeptPtRel/F");
    outtree->Branch("regjet1_softLeptDR", &t->regjet1_softLeptDR, "regjet1_softLeptDR/F");
    outtree->Branch("regjet1_leadTrackPt", &t->regjet1_leadTrackPt, "regjet1_leadTrackPt/F");
    outtree->Branch("regjet1_JECUnc", &t->regjet1_JECUnc, "regjet1_JECUnc/F");
    outtree->Branch("regjet1_secVtxPt", &t->regjet1_secVtxPt, "regjet1_secVtxPt/F");
    outtree->Branch("regjet1_secVtx3dL", &t->regjet1_secVtx3dL, "regjet1_secVtx3dL/F");
    outtree->Branch("regjet1_secVtx3deL", &t->regjet1_secVtx3deL, "regjet1_secVtx3deL/F");
    outtree->Branch("regjet1_secVtxM", &t->regjet1_secVtxM, "regjet1_secVtxM/F");
    outtree->Branch("regjet1_dPhiMet", &t->regjet1_dPhiMet, "regjet1_dPhiMet/F");
    outtree->Branch("regjet1_nConstituents", &t->regjet1_nConstituents, "regjet1_nConstituents/I");
    outtree->Branch("regjet2_mt", &t->regjet2_mt, "regjet2_mt/F");
    outtree->Branch("regjet2_chadfrac", &t->regjet2_chadfrac, "regjet2_chadfrac/F");
    outtree->Branch("regjet2_nhadfrac", &t->regjet2_nhadfrac, "regjet2_nhadfrac/F");
    outtree->Branch("regjet2_phofrac", &t->regjet2_phofrac, "regjet2_phofrac/F");
    outtree->Branch("regjet2_mufrac", &t->regjet2_mufrac, "regjet2_mufrac/F");
    outtree->Branch("regjet2_elefrac", &t->regjet2_elefrac, "regjet2_elefrac/F");
    outtree->Branch("regjet2_softLeptPt", &t->regjet2_softLeptPt, "regjet2_softLeptPt/F");
    outtree->Branch("regjet2_softLeptPtRel", &t->regjet2_softLeptPtRel, "regjet2_softLeptPtRel/F");
    outtree->Branch("regjet2_softLeptDR", &t->regjet2_softLeptDR, "regjet2_softLeptDR/F");
    outtree->Branch("regjet2_leadTrackPt", &t->regjet2_leadTrackPt, "regjet2_leadTrackPt/F");
    outtree->Branch("regjet2_JECUnc", &t->regjet2_JECUnc, "regjet2_JECUnc/F");
    outtree->Branch("regjet2_secVtxPt", &t->regjet2_secVtxPt, "regjet2_secVtxPt/F");
    outtree->Branch("regjet2_secVtx3dL", &t->regjet2_secVtx3dL, "regjet2_secVtx3dL/F");
    outtree->Branch("regjet2_secVtx3deL", &t->regjet2_secVtx3deL, "regjet2_secVtx3deL/F");
    outtree->Branch("regjet2_secVtxM", &t->regjet2_secVtxM, "regjet2_secVtxM/F");
    outtree->Branch("regjet2_dPhiMet", &t->regjet2_dPhiMet, "regjet2_dPhiMet/F");
    outtree->Branch("regjet2_nConstituents", &t->regjet2_nConstituents, "regjet2_nConstituents/I");
// regressed / kin fitted jets
    outtree->Branch("regjet1_pt", &t->regjet1_pt, "regjet1_pt/F");
    outtree->Branch("regjet1_e", &t->regjet1_e, "regjet1_e/F");
    outtree->Branch("regjet1_phi", &t->regjet1_phi, "regjet1_phi/F");
    outtree->Branch("regjet1_eta", &t->regjet1_eta, "regjet1_eta/F");
    outtree->Branch("regjet1_mass", &t->regjet1_mass, "regjet1_mass/F");
    outtree->Branch("regjet1_csvBtag", &t->regjet1_csvBtag, "regjet1_csvBtag/F");
    outtree->Branch("regjet1_flavour", &t->regjet1_flavour, "regjet1_flavour/I");
    outtree->Branch("regjet1_btagSF_M", &t->regjet1_btagSF_M, "regjet1_btagSF_M/F");
    outtree->Branch("regjet1_btagSFErrorUp_M", &t->regjet1_btagSFErrorUp_M, "regjet1_btagSFErrorUp_M/F");
    outtree->Branch("regjet1_btagSFErrorDown_M", &t->regjet1_btagSFErrorDown_M, "regjet1_btagSFErrorDown_M/F");
    outtree->Branch("regjet1_btagEff_M", &t->regjet1_btagEff_M, "regjet1_btagEff_M/F");
    outtree->Branch("regjet1_btagEffError_M", &t->regjet1_btagEffError_M, "regjet1_btagEffError_M/F");
    outtree->Branch("regjet1_cutbased_wp_level", &t->regjet1_cutbased_wp_level, "regjet1_cutbased_wp_level/I");
    outtree->Branch("regjet1_simple_wp_level", &t->regjet1_simple_wp_level, "regjet1_simple_wp_level/I");
    outtree->Branch("regjet1_full_wp_level", &t->regjet1_full_wp_level, "regjet1_full_wp_level/I");
    outtree->Branch("regjet1_betaStarClassic", &t->regjet1_betaStarClassic, "regjet1_betaStarClassic/F");
    outtree->Branch("regjet1_dR2Mean", &t->regjet1_dR2Mean, "regjet1_dR2Mean/F");
    outtree->Branch("regjet2_pt", &t->regjet2_pt, "regjet2_pt/F");
    outtree->Branch("regjet2_e", &t->regjet2_e, "regjet2_e/F");
    outtree->Branch("regjet2_phi", &t->regjet2_phi, "regjet2_phi/F");
    outtree->Branch("regjet2_eta", &t->regjet2_eta, "regjet2_eta/F");
    outtree->Branch("regjet2_mass", &t->regjet2_mass, "regjet2_mass/F");
    outtree->Branch("regjet2_csvBtag", &t->regjet2_csvBtag, "regjet2_csvBtag/F");
    outtree->Branch("regjet2_flavour", &t->regjet2_flavour, "regjet2_flavour/I");
    outtree->Branch("regjet2_btagSF_M", &t->regjet2_btagSF_M, "regjet2_btagSF_M/F");
    outtree->Branch("regjet2_btagSFErrorUp_M", &t->regjet2_btagSFErrorUp_M, "regjet2_btagSFErrorUp_M/F");
    outtree->Branch("regjet2_btagSFErrorDown_M", &t->regjet2_btagSFErrorDown_M, "regjet2_btagSFErrorDown_M/F");
    outtree->Branch("regjet2_btagEff_M", &t->regjet2_btagEff_M, "regjet2_btagEff_M/F");
    outtree->Branch("regjet2_btagEffError_M", &t->regjet2_btagEffError_M, "regjet2_btagEffError_M/F");
    outtree->Branch("regjet2_cutbased_wp_level", &t->regjet2_cutbased_wp_level, "regjet2_cutbased_wp_level/I");
    outtree->Branch("regjet2_simple_wp_level", &t->regjet2_simple_wp_level, "regjet2_simple_wp_level/I");
    outtree->Branch("regjet2_full_wp_level", &t->regjet2_full_wp_level, "regjet2_full_wp_level/I");
    outtree->Branch("regjet2_betaStarClassic", &t->regjet2_betaStarClassic, "regjet2_betaStarClassic/F");
    outtree->Branch("regjet2_dR2Mean", &t->regjet2_dR2Mean, "regjet2_dR2Mean/F");
    outtree->Branch("regkinjet1_pt", &t->regkinjet1_pt, "regkinjet1_pt/F");
    outtree->Branch("regkinjet1_e", &t->regkinjet1_e, "regkinjet1_e/F");
    outtree->Branch("regkinjet1_phi", &t->regkinjet1_phi, "regkinjet1_phi/F");
    outtree->Branch("regkinjet1_eta", &t->regkinjet1_eta, "regkinjet1_eta/F");
    outtree->Branch("regkinjet1_mass", &t->regkinjet1_mass, "regkinjet1_mass/F");
    outtree->Branch("regkinjet1_csvBtag", &t->regkinjet1_csvBtag, "regkinjet1_csvBtag/F");
    outtree->Branch("regkinjet1_flavour", &t->regkinjet1_flavour, "regkinjet1_flavour/I");
    outtree->Branch("regkinjet1_btagSF_M", &t->regkinjet1_btagSF_M, "regkinjet1_btagSF_M/F");
    outtree->Branch("regkinjet1_btagSFErrorUp_M", &t->regkinjet1_btagSFErrorUp_M, "regkinjet1_btagSFErrorUp_M/F");
    outtree->Branch("regkinjet1_btagSFErrorDown_M", &t->regkinjet1_btagSFErrorDown_M, "regkinjet1_btagSFErrorDown_M/F");
    outtree->Branch("regkinjet1_btagEff_M", &t->regkinjet1_btagEff_M, "regkinjet1_btagEff_M/F");
    outtree->Branch("regkinjet1_btagEffError_M", &t->regkinjet1_btagEffError_M, "regkinjet1_btagEffError_M/F");
    outtree->Branch("regkinjet1_cutbased_wp_level", &t->regkinjet1_cutbased_wp_level, "regkinjet1_cutbased_wp_level/I");
    outtree->Branch("regkinjet1_simple_wp_level", &t->regkinjet1_simple_wp_level, "regkinjet1_simple_wp_level/I");
    outtree->Branch("regkinjet1_full_wp_level", &t->regkinjet1_full_wp_level, "regkinjet1_full_wp_level/I");
    outtree->Branch("regkinjet1_betaStarClassic", &t->regkinjet1_betaStarClassic, "regkinjet1_betaStarClassic/F");
    outtree->Branch("regkinjet1_dR2Mean", &t->regkinjet1_dR2Mean, "regkinjet1_dR2Mean/F");
    outtree->Branch("regkinjet2_pt", &t->regkinjet2_pt, "regkinjet2_pt/F");
    outtree->Branch("regkinjet2_e", &t->regkinjet2_e, "regkinjet2_e/F");
    outtree->Branch("regkinjet2_phi", &t->regkinjet2_phi, "regkinjet2_phi/F");
    outtree->Branch("regkinjet2_eta", &t->regkinjet2_eta, "regkinjet2_eta/F");
    outtree->Branch("regkinjet2_mass", &t->regkinjet2_mass, "regkinjet2_mass/F");
    outtree->Branch("regkinjet2_csvBtag", &t->regkinjet2_csvBtag, "regkinjet2_csvBtag/F");
    outtree->Branch("regkinjet2_flavour", &t->regkinjet2_flavour, "regkinjet2_flavour/I");
    outtree->Branch("regkinjet2_btagSF_M", &t->regkinjet2_btagSF_M, "regkinjet2_btagSF_M/F");
    outtree->Branch("regkinjet2_btagSFErrorUp_M", &t->regkinjet2_btagSFErrorUp_M, "regkinjet2_btagSFErrorUp_M/F");
    outtree->Branch("regkinjet2_btagSFErrorDown_M", &t->regkinjet2_btagSFErrorDown_M, "regkinjet2_btagSFErrorDown_M/F");
    outtree->Branch("regkinjet2_btagEffError_M", &t->regkinjet2_btagEffError_M, "regkinjet2_btagEffError_M/F");
    outtree->Branch("regkinjet2_cutbased_wp_level", &t->regkinjet2_cutbased_wp_level, "regkinjet2_cutbased_wp_level/I");
    outtree->Branch("regkinjet2_simple_wp_level", &t->regkinjet2_simple_wp_level, "regkinjet2_simple_wp_level/I");
    outtree->Branch("regkinjet2_full_wp_level", &t->regkinjet2_full_wp_level, "regkinjet2_full_wp_level/I");
    outtree->Branch("regkinjet2_btagEff_M", &t->regkinjet2_btagEff_M, "regkinjet2_btagEff_M/F");
    outtree->Branch("regkinjet2_betaStarClassic", &t->regkinjet2_betaStarClassic, "regkinjet2_betaStarClassic/F");
    outtree->Branch("regkinjet2_dR2Mean", &t->regkinjet2_dR2Mean, "regkinjet2_dR2Mean/F");
    outtree->Branch("kinjet1_pt", &t->kinjet1_pt, "kinjet1_pt/F");
    outtree->Branch("kinjet1_e", &t->kinjet1_e, "kinjet1_e/F");
    outtree->Branch("kinjet1_phi", &t->kinjet1_phi, "kinjet1_phi/F");
    outtree->Branch("kinjet1_eta", &t->kinjet1_eta, "kinjet1_eta/F");
    outtree->Branch("kinjet1_mass", &t->kinjet1_mass, "kinjet1_mass/F");
    outtree->Branch("kinjet1_csvBtag", &t->kinjet1_csvBtag, "kinjet1_csvBtag/F");
    outtree->Branch("kinjet1_flavour", &t->kinjet1_flavour, "kinjet1_flavour/I");
    outtree->Branch("kinjet1_btagSF_M", &t->kinjet1_btagSF_M, "kinjet1_btagSF_M/F");
    outtree->Branch("kinjet1_btagSFErrorUp_M", &t->kinjet1_btagSFErrorUp_M, "kinjet1_btagSFErrorUp_M/F");
    outtree->Branch("kinjet1_btagSFErrorDown_M", &t->kinjet1_btagSFErrorDown_M, "kinjet1_btagSFErrorDown_M/F");
    outtree->Branch("kinjet1_btagEffError_M", &t->kinjet1_btagEffError_M, "kinjet1_btagEffError_M/F");
    outtree->Branch("kinjet1_cutbased_wp_level", &t->kinjet1_cutbased_wp_level, "kinjet1_cutbased_wp_level/I");
    outtree->Branch("kinjet1_simple_wp_level", &t->kinjet1_simple_wp_level, "kinjet1_simple_wp_level/I");
    outtree->Branch("kinjet1_full_wp_level", &t->kinjet1_full_wp_level, "kinjet1_full_wp_level/I");
    outtree->Branch("kinjet1_btagEff_M", &t->kinjet1_btagEff_M, "kinjet1_btagEff_M/F");
    outtree->Branch("kinjet1_betaStarClassic", &t->kinjet1_betaStarClassic, "kinjet1_betaStarClassic/F");
    outtree->Branch("kinjet1_dR2Mean", &t->kinjet1_dR2Mean, "kinjet1_dR2Mean/F");
    outtree->Branch("kinjet2_pt", &t->kinjet2_pt, "kinjet2_pt/F");
    outtree->Branch("kinjet2_e", &t->kinjet2_e, "kinjet2_e/F");
    outtree->Branch("kinjet2_phi", &t->kinjet2_phi, "kinjet2_phi/F");
    outtree->Branch("kinjet2_eta", &t->kinjet2_eta, "kinjet2_eta/F");
    outtree->Branch("kinjet2_mass", &t->kinjet2_mass, "kinjet2_mass/F");
    outtree->Branch("kinjet2_csvBtag", &t->kinjet2_csvBtag, "kinjet2_csvBtag/F");
    outtree->Branch("kinjet2_flavour", &t->kinjet2_flavour, "kinjet2_flavour/I");
    outtree->Branch("kinjet2_btagSF_M", &t->kinjet2_btagSF_M, "kinjet2_btagSF_M/F");
    outtree->Branch("kinjet2_btagSFErrorUp_M", &t->kinjet2_btagSFErrorUp_M, "kinjet2_btagSFErrorUp_M/F");
    outtree->Branch("kinjet2_btagSFErrorDown_M", &t->kinjet2_btagSFErrorDown_M, "kinjet2_btagSFErrorDown_M/F");
    outtree->Branch("kinjet2_btagEffError_M", &t->kinjet2_btagEffError_M, "kinjet2_btagEffError_M/F");
    outtree->Branch("kinjet2_cutbased_wp_level", &t->kinjet2_cutbased_wp_level, "kinjet2_cutbased_wp_level/I");
    outtree->Branch("kinjet2_simple_wp_level", &t->kinjet2_simple_wp_level, "kinjet2_simple_wp_level/I");
    outtree->Branch("kinjet2_full_wp_level", &t->kinjet2_full_wp_level, "kinjet2_full_wp_level/I");
    outtree->Branch("kinjet2_btagEff_M", &t->kinjet2_btagEff_M, "kinjet2_btagEff_M/F");
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
    outtree->Branch("gg_DR", &t->gg_DR, "gg_DR/F");
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
    outtree->Branch("costhetastar_CS", &t->costhetastar_CS, "costhetastar_CS/F");
    outtree->Branch("regcosthetastar_CS", &t->regcosthetastar_CS, "regcosthetastar_CS/F");
    outtree->Branch("regkincosthetastar_CS", &t->regkincosthetastar_CS, "regkincosthetastar_CS/F");
    outtree->Branch("kincosthetastar_CS", &t->kincosthetastar_CS, "kincosthetastar_CS/F");
    outtree->Branch("minDRgj", &t->minDRgj, "minDRgj/F");
    outtree->Branch("minDRgregj", &t->minDRgregj, "minDRgregj/F");
    outtree->Branch("minDRgregkinj", &t->minDRgregkinj, "minDRgregkinj/F");
    outtree->Branch("minDRgkinj", &t->minDRgkinj, "minDRgkinj/F");
    outtree->Branch("HT_gg", &t->HT_gg, "HT_gg/F");
    outtree->Branch("dEta_gg_jj", &t->dEta_gg_jj, "dEta_gg_jj/F");
    outtree->Branch("dEta_gg_regjj", &t->dEta_gg_regjj, "dEta_gg_regjj/F");
    outtree->Branch("dEta_gg_regkinjj", &t->dEta_gg_regkinjj, "dEta_gg_regkinjj/F");
    outtree->Branch("dEta_gg_kinjj", &t->dEta_gg_kinjj, "dEta_gg_kinjj/F");
    outtree->Branch("dPhi_gg_jj", &t->dPhi_gg_jj, "dPhi_gg_jj/F");
    outtree->Branch("dPhi_gg_regjj", &t->dPhi_gg_regjj, "dPhi_gg_regjj/F");
    outtree->Branch("dPhi_gg_regkinjj", &t->dPhi_gg_regkinjj, "dPhi_gg_regkinjj/F");
    outtree->Branch("dPhi_gg_kinjj", &t->dPhi_gg_kinjj, "dPhi_gg_kinjj/F");
    outtree->Branch("dR_gg_jj", &t->dR_gg_jj, "dR_gg_jj/F");
    outtree->Branch("dR_gg_regjj", &t->dR_gg_regjj, "dR_gg_regjj/F");
    outtree->Branch("dR_gg_regkinjj", &t->dR_gg_regkinjj, "dR_gg_regkinjj/F");
    outtree->Branch("dR_gg_kinjj", &t->dR_gg_kinjj, "dR_gg_kinjj/F");
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
    outtree->Branch("gr_hbbhgg_costhetastar_CS", &t->gr_hbbhgg_costhetastar_CS, "gr_hbbhgg_costhetastar_CS/F");
    outtree->Branch("gr_hjjhgg_costhetastar_CS", &t->gr_hjjhgg_costhetastar_CS, "gr_hjjhgg_costhetastar_CS/F");
    outtree->Branch("gr_dEta_gg_bb", &t->gr_dEta_gg_bb, "gr_dEta_gg_bb/F");
    outtree->Branch("gr_dEta_gg_jj", &t->gr_dEta_gg_jj, "gr_dEta_gg_jj/F");
    outtree->Branch("gr_dPhi_gg_bb", &t->gr_dPhi_gg_bb, "gr_dPhi_gg_bb/F");
    outtree->Branch("gr_dPhi_gg_jj", &t->gr_dPhi_gg_jj, "gr_dPhi_gg_jj/F");
    outtree->Branch("gr_dR_gg_bb", &t->gr_dR_gg_bb, "gr_dR_gg_bb/F");
    outtree->Branch("gr_dR_gg_jj", &t->gr_dR_gg_jj, "gr_dR_gg_jj/F");
    // Photon Energy Scale & Photon Energy Resolution
    outtree->Branch("gg_mass_pesD1", &t->gg_mass_pesD1, "gg_mass_pesD1/F");
    outtree->Branch("gg_mass_pesU1", &t->gg_mass_pesU1, "gg_mass_pesU1/F");
    outtree->Branch("gg_mass_perD1", &t->gg_mass_perD1, "gg_mass_perD1/F");
    outtree->Branch("gg_mass_perU1", &t->gg_mass_perU1, "gg_mass_perU1/F");
    outtree->Branch("gg_mass_pesD2", &t->gg_mass_pesD2, "gg_mass_pesD2/F");
    outtree->Branch("gg_mass_pesU2", &t->gg_mass_pesU2, "gg_mass_pesU2/F");
    outtree->Branch("gg_mass_perD2", &t->gg_mass_perD2, "gg_mass_perD2/F");
    outtree->Branch("gg_mass_perU2", &t->gg_mass_perU2, "gg_mass_perU2/F");
    outtree->Branch("gg_mass_pesD", &t->gg_mass_pesD, "gg_mass_pesD/F");
    outtree->Branch("gg_mass_pesU", &t->gg_mass_pesU, "gg_mass_pesU/F");
    outtree->Branch("gg_mass_perD", &t->gg_mass_perD, "gg_mass_perD/F");
    outtree->Branch("gg_mass_perU", &t->gg_mass_perU, "gg_mass_perU/F");
    outtree->Branch("regggjj_mass_pesD1", &t->regggjj_mass_pesD1, "regggjj_mass_pesD1/F");
    outtree->Branch("regggjj_mass_pesU1", &t->regggjj_mass_pesU1, "regggjj_mass_pesU1/F");
    outtree->Branch("regggjj_mass_perD1", &t->regggjj_mass_perD1, "regggjj_mass_perD1/F");
    outtree->Branch("regggjj_mass_perU1", &t->regggjj_mass_perU1, "regggjj_mass_perU1/F");
    outtree->Branch("regggjj_mass_pesD2", &t->regggjj_mass_pesD2, "regggjj_mass_pesD2/F");
    outtree->Branch("regggjj_mass_pesU2", &t->regggjj_mass_pesU2, "regggjj_mass_pesU2/F");
    outtree->Branch("regggjj_mass_perD2", &t->regggjj_mass_perD2, "regggjj_mass_perD2/F");
    outtree->Branch("regggjj_mass_perU2", &t->regggjj_mass_perU2, "regggjj_mass_perU2/F");
    outtree->Branch("regggjj_mass_pesD", &t->regggjj_mass_pesD, "regggjj_mass_pesD/F");
    outtree->Branch("regggjj_mass_pesU", &t->regggjj_mass_pesU, "regggjj_mass_pesU/F");
    outtree->Branch("regggjj_mass_perD", &t->regggjj_mass_perD, "regggjj_mass_perD/F");
    outtree->Branch("regggjj_mass_perU", &t->regggjj_mass_perU, "regggjj_mass_perU/F");
    outtree->Branch("regkinggjj_mass_pesD1", &t->regkinggjj_mass_pesD1, "regkinggjj_mass_pesD1/F");
    outtree->Branch("regkinggjj_mass_pesU1", &t->regkinggjj_mass_pesU1, "regkinggjj_mass_pesU1/F");
    outtree->Branch("regkinggjj_mass_perD1", &t->regkinggjj_mass_perD1, "regkinggjj_mass_perD1/F");
    outtree->Branch("regkinggjj_mass_perU1", &t->regkinggjj_mass_perU1, "regkinggjj_mass_perU1/F");
    outtree->Branch("regkinggjj_mass_pesD2", &t->regkinggjj_mass_pesD2, "regkinggjj_mass_pesD2/F");
    outtree->Branch("regkinggjj_mass_pesU2", &t->regkinggjj_mass_pesU2, "regkinggjj_mass_pesU2/F");
    outtree->Branch("regkinggjj_mass_perD2", &t->regkinggjj_mass_perD2, "regkinggjj_mass_perD2/F");
    outtree->Branch("regkinggjj_mass_perU2", &t->regkinggjj_mass_perU2, "regkinggjj_mass_perU2/F");
    outtree->Branch("regkinggjj_mass_pesD", &t->regkinggjj_mass_pesD, "regkinggjj_mass_pesD/F");
    outtree->Branch("regkinggjj_mass_pesU", &t->regkinggjj_mass_pesU, "regkinggjj_mass_pesU/F");
    outtree->Branch("regkinggjj_mass_perD", &t->regkinggjj_mass_perD, "regkinggjj_mass_perD/F");
    outtree->Branch("regkinggjj_mass_perU", &t->regkinggjj_mass_perU, "regkinggjj_mass_perU/F");
    outtree->Branch("ggjj_mass_pesD1", &t->ggjj_mass_pesD1, "ggjj_mass_pesD1/F");
    outtree->Branch("ggjj_mass_pesU1", &t->ggjj_mass_pesU1, "ggjj_mass_pesU1/F");
    outtree->Branch("ggjj_mass_perD1", &t->ggjj_mass_perD1, "ggjj_mass_perD1/F");
    outtree->Branch("ggjj_mass_perU1", &t->ggjj_mass_perU1, "ggjj_mass_perU1/F");
    outtree->Branch("ggjj_mass_pesD2", &t->ggjj_mass_pesD2, "ggjj_mass_pesD2/F");
    outtree->Branch("ggjj_mass_pesU2", &t->ggjj_mass_pesU2, "ggjj_mass_pesU2/F");
    outtree->Branch("ggjj_mass_perD2", &t->ggjj_mass_perD2, "ggjj_mass_perD2/F");
    outtree->Branch("ggjj_mass_perU2", &t->ggjj_mass_perU2, "ggjj_mass_perU2/F");
    outtree->Branch("ggjj_mass_pesD", &t->ggjj_mass_pesD, "ggjj_mass_pesD/F");
    outtree->Branch("ggjj_mass_pesU", &t->ggjj_mass_pesU, "ggjj_mass_pesU/F");
    outtree->Branch("ggjj_mass_perD", &t->ggjj_mass_perD, "ggjj_mass_perD/F");
    outtree->Branch("ggjj_mass_perU", &t->ggjj_mass_perU, "ggjj_mass_perU/F");
    outtree->Branch("kinggjj_mass_pesD1", &t->kinggjj_mass_pesD1, "kinggjj_mass_pesD1/F");
    outtree->Branch("kinggjj_mass_pesU1", &t->kinggjj_mass_pesU1, "kinggjj_mass_pesU1/F");
    outtree->Branch("kinggjj_mass_perD1", &t->kinggjj_mass_perD1, "kinggjj_mass_perD1/F");
    outtree->Branch("kinggjj_mass_perU1", &t->kinggjj_mass_perU1, "kinggjj_mass_perU1/F");
    outtree->Branch("kinggjj_mass_pesD2", &t->kinggjj_mass_pesD2, "kinggjj_mass_pesD2/F");
    outtree->Branch("kinggjj_mass_pesU2", &t->kinggjj_mass_pesU2, "kinggjj_mass_pesU2/F");
    outtree->Branch("kinggjj_mass_perD2", &t->kinggjj_mass_perD2, "kinggjj_mass_perD2/F");
    outtree->Branch("kinggjj_mass_perU2", &t->kinggjj_mass_perU2, "kinggjj_mass_perU2/F");
    outtree->Branch("kinggjj_mass_pesD", &t->kinggjj_mass_pesD, "kinggjj_mass_pesD/F");
    outtree->Branch("kinggjj_mass_pesU", &t->kinggjj_mass_pesU, "kinggjj_mass_pesU/F");
    outtree->Branch("kinggjj_mass_perD", &t->kinggjj_mass_perD, "kinggjj_mass_perD/F");
    outtree->Branch("kinggjj_mass_perU", &t->kinggjj_mass_perU, "kinggjj_mass_perU/F");
    // Jet Energy Correction and Jet Energy Resolution
    outtree->Branch("jj_mass_jecD", &t->jj_mass_jecD, "jj_mass_jecD/F");
    outtree->Branch("jj_mass_jecU", &t->jj_mass_jecU, "jj_mass_jecU/F");
    outtree->Branch("jj_mass_jerD", &t->jj_mass_jerD, "jj_mass_jerD/F");
    outtree->Branch("jj_mass_jerC", &t->jj_mass_jerC, "jj_mass_jerC/F");
    outtree->Branch("jj_mass_jerU", &t->jj_mass_jerU, "jj_mass_jerU/F");
    outtree->Branch("kinjj_mass_jecD", &t->kinjj_mass_jecD, "kinjj_mass_jecD/F");
    outtree->Branch("kinjj_mass_jecU", &t->kinjj_mass_jecU, "kinjj_mass_jecU/F");
    outtree->Branch("kinjj_mass_jerD", &t->kinjj_mass_jerD, "kinjj_mass_jerD/F");
    outtree->Branch("kinjj_mass_jerC", &t->kinjj_mass_jerC, "kinjj_mass_jerC/F");
    outtree->Branch("kinjj_mass_jerU", &t->kinjj_mass_jerU, "kinjj_mass_jerU/F");
    outtree->Branch("ggjj_mass_jecD", &t->ggjj_mass_jecD, "ggjj_mass_jecD/F");
    outtree->Branch("ggjj_mass_jecU", &t->ggjj_mass_jecU, "ggjj_mass_jecU/F");
    outtree->Branch("ggjj_mass_jerD", &t->ggjj_mass_jerD, "ggjj_mass_jerD/F");
    outtree->Branch("ggjj_mass_jerC", &t->ggjj_mass_jerC, "ggjj_mass_jerC/F");
    outtree->Branch("ggjj_mass_jerU", &t->ggjj_mass_jerU, "ggjj_mass_jerU/F");
    outtree->Branch("kinggjj_mass_jecD", &t->kinggjj_mass_jecD, "kinggjj_mass_jecD/F");
    outtree->Branch("kinggjj_mass_jecU", &t->kinggjj_mass_jecU, "kinggjj_mass_jecU/F");
    outtree->Branch("kinggjj_mass_jerD", &t->kinggjj_mass_jerD, "kinggjj_mass_jerD/F");
    outtree->Branch("kinggjj_mass_jerC", &t->kinggjj_mass_jerC, "kinggjj_mass_jerC/F");
    outtree->Branch("kinggjj_mass_jerU", &t->kinggjj_mass_jerU, "kinggjj_mass_jerU/F");
    
    return;
}

struct jet_variables
{
        std::vector<float> jetPt;
        std::vector<float> jetbtagSF_M;
        std::vector<int> jetflavour;
        std::vector<float> jetbtagSFErrorUp_M;
        std::vector<float> jetbtagSFErrorDown_M;
        std::vector<float> jetbtagEff_M;
        std::vector<float> jetbtagEffError_M;
        std::vector<int> jetcutbased_wp_level;
        std::vector<int> jetsimple_wp_level;
        std::vector<int> jetfull_wp_level;
        std::vector<float> jetbetaStarClassic;
        std::vector<float> jetdR2Mean;
        std::vector<float> jetE;
        std::vector<float> jetEta;
        std::vector<float> jetPhi;
        std::vector<float> jetCSV;
        std::vector<float> jetRegPt;
        std::vector<float> jetRegKinPt;
// regression inputs
        std::vector<float> jetMt;
        std::vector<float> jetChadfrac;
        std::vector<float> jetNhadfrac;
        std::vector<float> jetPhofrac;
        std::vector<float> jetMufrac;
        std::vector<float> jetElefrac;
        std::vector<float> jetSoftLeptPt;
        std::vector<float> jetSoftLeptPtRel;
        std::vector<float> jetSoftLeptDR;
        std::vector<float> jetLeadTrackPt;
        std::vector<float> jetJECUnc;
        std::vector<float> jetSecVtxPt;
        std::vector<float> jetSecVtx3dL;
        std::vector<float> jetSecVtx3deL;
        std::vector<float> jetSecVtxM;
        std::vector<float> jetDPhiMet;
        std::vector<int> jetNConstituents;
// Jet Energy Correction and Jet Energy Resolution
        std::vector<float> jetJecD_e, jetJecD_pt, jetJecD_phi, jetJecD_eta;
        std::vector<float> jetJecU_e, jetJecU_pt, jetJecU_phi, jetJecU_eta;
        std::vector<float> jetJerD_e, jetJerD_pt, jetJerD_phi, jetJerD_eta;
        std::vector<float> jetJerC_e, jetJerC_pt, jetJerC_phi, jetJerC_eta;
        std::vector<float> jetJerU_e, jetJerU_pt, jetJerU_phi, jetJerU_eta;
};


void initialize_jet_variables( jet_variables * J )
{
        J->jetPt.clear();
        J->jetbtagSF_M.clear();
        J->jetflavour.clear();
        J->jetbtagSFErrorUp_M.clear();
        J->jetbtagSFErrorDown_M.clear();
        J->jetbtagEff_M.clear();
        J->jetbtagEffError_M.clear();
        J->jetcutbased_wp_level.clear();
        J->jetsimple_wp_level.clear();
        J->jetfull_wp_level.clear();
        J->jetbetaStarClassic.clear();
        J->jetdR2Mean.clear();
        J->jetE.clear();
        J->jetEta.clear();
        J->jetPhi.clear();
        J->jetCSV.clear();
        J->jetRegPt.clear();
        J->jetRegKinPt.clear();
// regression inputs
        J->jetMt.clear();
        J->jetChadfrac.clear();
        J->jetNhadfrac.clear();
        J->jetPhofrac.clear();
        J->jetMufrac.clear();
        J->jetElefrac.clear();
        J->jetSoftLeptPt.clear();
        J->jetSoftLeptPtRel.clear();
        J->jetSoftLeptDR.clear();
        J->jetLeadTrackPt.clear();
        J->jetJECUnc.clear();
        J->jetSecVtxPt.clear();
        J->jetSecVtx3dL.clear();
        J->jetSecVtx3deL.clear();
        J->jetSecVtxM.clear();
        J->jetDPhiMet.clear();
        J->jetNConstituents.clear();
// Jet Energy Correction and Jet Energy Resolution
        J->jetJecD_e.clear(); J->jetJecD_pt.clear(); J->jetJecD_phi.clear(); J->jetJecD_eta.clear();
        J->jetJecU_e.clear(); J->jetJecU_pt.clear(); J->jetJecU_phi.clear(); J->jetJecU_eta.clear();
        J->jetJerD_e.clear(); J->jetJerD_pt.clear(); J->jetJerD_phi.clear(); J->jetJerD_eta.clear();
        J->jetJerC_e.clear(); J->jetJerC_pt.clear(); J->jetJerC_phi.clear(); J->jetJerC_eta.clear();
        J->jetJerU_e.clear(); J->jetJerU_pt.clear(); J->jetJerU_phi.clear(); J->jetJerU_eta.clear();

    return ;
}


bool fill_jet_variables( tree_variables * t, int ijet, TLorentzVector met, int numberOfRegressionFiles)
{
        TLorentzVector jet;
            if( ijet == 0 && t->j1_pt > .01 )
            {
                t->jet_e = t->j1_e;
                t->jet_pt = t->j1_pt;
                t->jet_phi = t->j1_phi;
                t->jet_eta = t->j1_eta;
                t->jet_betaStarClassic = t->j1_betaStarClassic;
                t->jet_dR2Mean = t->j1_dR2Mean;
                t->jet_csvBtag = t->j1_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j1_secVtxPt;
                    t->jet_secVtx3dL = t->j1_secVtx3dL;
                    t->jet_secVtx3deL = t->j1_secVtx3deL;
                    t->jet_secVtxM = t->j1_secVtxM;
            		t->jet_mt = pow(t->j1_e,2)-pow(t->j1_pt*(1+sinh(t->j1_eta)),2) > 0 ? sqrt(pow(t->j1_e,2)-pow(t->j1_pt*(1+sinh(t->j1_eta)),2)) : sqrt(pow(t->j1_pt*(1+sinh(t->j1_eta)),2)-pow(t->j1_e,2)) ;
                    t->jet_emfrac = t->j1_emfrac;
                    t->jet_hadfrac = t->j1_hadfrac;
                    t->jet_chadfrac = t->j1_chadfrac;
                    t->jet_nhadfrac = t->j1_nhadfrac;
                    t->jet_phofrac = t->j1_phofrac;
                    t->jet_mufrac = t->j1_mufrac;
                    t->jet_phofrac = t->j1_phofrac;
                    t->jet_JECUnc = t->j1_JECUnc;
                    t->jet_leadTrackPt = t->j1_leadTrackPt;
                    t->jet_softLeptPt = (t->j1_softLeptIdLooseMu==1 || t->j1_softLeptIdEle95==1) ? (t->j1_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j1_softLeptIdLooseMu==1 || t->j1_softLeptIdEle95==1) ? (t->j1_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j1_softLeptIdLooseMu==1 || t->j1_softLeptIdEle95==1) ? (t->j1_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j1_btagSF_M;
                t->jet_flavour = t->j1_flavour;
                t->jet_btagSFErrorUp_M = t->j1_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j1_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j1_btagEff_M;
                t->jet_btagEffError_M = t->j1_btagEffError_M;
                t->jet_cutbased_wp_level = t->j1_cutbased_wp_level;
                t->jet_simple_wp_level = t->j1_simple_wp_level;
                t->jet_full_wp_level = t->j1_full_wp_level;
                t->jet_betaStarClassic = t->j1_betaStarClassic;
                t->jet_dR2Mean = t->j1_dR2Mean;
                t->jet_nNeutrals = t->j1_nNeutrals;
                t->jet_nCharged = t->j1_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j1_pt, t->j1_eta, t->j1_phi, t->j1_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j1_jecD_e;
                t->jet_jecD_pt = t->j1_jecD_pt;
                t->jet_jecD_phi = t->j1_jecD_phi;
                t->jet_jecD_eta = t->j1_jecD_eta;
                t->jet_jecU_e = t->j1_jecU_e;
                t->jet_jecU_pt = t->j1_jecU_pt;
                t->jet_jecU_phi = t->j1_jecU_phi;
                t->jet_jecU_eta = t->j1_jecU_eta;
                t->jet_jerD_e = t->j1_jerD_e;
                t->jet_jerD_pt = t->j1_jerD_pt;
                t->jet_jerD_phi = t->j1_jerD_phi;
                t->jet_jerD_eta = t->j1_jerD_eta;
                t->jet_jerC_e = t->j1_jerC_e;
                t->jet_jerC_pt = t->j1_jerC_pt;
                t->jet_jerC_phi = t->j1_jerC_phi;
                t->jet_jerC_eta = t->j1_jerC_eta;
                t->jet_jerU_e = t->j1_jerU_e;
                t->jet_jerU_pt = t->j1_jerU_pt;
                t->jet_jerU_phi = t->j1_jerU_phi;
                t->jet_jerU_eta = t->j1_jerU_eta;
            } 
            else
            if( ijet == 1 && t->j2_pt > .01 )
            {
                t->jet_e = t->j2_e;
                t->jet_pt = t->j2_pt;
                t->jet_phi = t->j2_phi;
                t->jet_eta = t->j2_eta;
                t->jet_betaStarClassic = t->j2_betaStarClassic;
                t->jet_dR2Mean = t->j2_dR2Mean;
                t->jet_csvBtag = t->j2_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j2_secVtxPt;
                    t->jet_secVtx3dL = t->j2_secVtx3dL;
                    t->jet_secVtx3deL = t->j2_secVtx3deL;
                    t->jet_secVtxM = t->j2_secVtxM;
             		t->jet_mt = pow(t->j2_e,2)-pow(t->j2_pt*(1+sinh(t->j2_eta)),2) > 0 ? sqrt(pow(t->j2_e,2)-pow(t->j2_pt*(1+sinh(t->j2_eta)),2)) : sqrt(pow(t->j2_pt*(1+sinh(t->j2_eta)),2)-pow(t->j2_e,2)) ;
                    t->jet_emfrac = t->j2_emfrac;
                    t->jet_hadfrac = t->j2_hadfrac;
                    t->jet_chadfrac = t->j2_chadfrac;
                    t->jet_nhadfrac = t->j2_nhadfrac;
                    t->jet_phofrac = t->j2_phofrac;
                    t->jet_mufrac = t->j2_mufrac;
                    t->jet_phofrac = t->j2_phofrac;
                    t->jet_JECUnc = t->j2_JECUnc;
                    t->jet_leadTrackPt = t->j2_leadTrackPt;
                    t->jet_softLeptPt = (t->j2_softLeptIdLooseMu==1 || t->j2_softLeptIdEle95==1) ? (t->j3_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j2_softLeptIdLooseMu==1 || t->j2_softLeptIdEle95==1) ? (t->j3_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j2_softLeptIdLooseMu==1 || t->j2_softLeptIdEle95==1) ? (t->j3_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j2_btagSF_M;
                t->jet_flavour = t->j2_flavour;
                t->jet_btagSFErrorUp_M = t->j2_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j2_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j2_btagEff_M;
                t->jet_btagEffError_M = t->j2_btagEffError_M;
                t->jet_cutbased_wp_level = t->j2_cutbased_wp_level;
                t->jet_simple_wp_level = t->j2_simple_wp_level;
                t->jet_full_wp_level = t->j2_full_wp_level;
                t->jet_betaStarClassic = t->j2_betaStarClassic;
                t->jet_dR2Mean = t->j2_dR2Mean;
                t->jet_nNeutrals = t->j2_nNeutrals;
                t->jet_nCharged = t->j2_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j2_pt, t->j2_eta, t->j2_phi, t->j2_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j2_jecD_e;
                t->jet_jecD_pt = t->j2_jecD_pt;
                t->jet_jecD_phi = t->j2_jecD_phi;
                t->jet_jecD_eta = t->j2_jecD_eta;
                t->jet_jecU_e = t->j2_jecU_e;
                t->jet_jecU_pt = t->j2_jecU_pt;
                t->jet_jecU_phi = t->j2_jecU_phi;
                t->jet_jecU_eta = t->j2_jecU_eta;
                t->jet_jerD_e = t->j2_jerD_e;
                t->jet_jerD_pt = t->j2_jerD_pt;
                t->jet_jerD_phi = t->j2_jerD_phi;
                t->jet_jerD_eta = t->j2_jerD_eta;
                t->jet_jerC_e = t->j2_jerC_e;
                t->jet_jerC_pt = t->j2_jerC_pt;
                t->jet_jerC_phi = t->j2_jerC_phi;
                t->jet_jerC_eta = t->j2_jerC_eta;
                t->jet_jerU_e = t->j2_jerU_e;
                t->jet_jerU_pt = t->j2_jerU_pt;
                t->jet_jerU_phi = t->j2_jerU_phi;
                t->jet_jerU_eta = t->j2_jerU_eta;
            } 
            else
            if( ijet == 2 && t->j3_pt > .01 )
            {
                t->jet_e = t->j3_e;
                t->jet_pt = t->j3_pt;
                t->jet_phi = t->j3_phi;
                t->jet_eta = t->j3_eta;
                t->jet_betaStarClassic = t->j3_betaStarClassic;
                t->jet_dR2Mean = t->j3_dR2Mean;
                t->jet_csvBtag = t->j3_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j3_secVtxPt;
                    t->jet_secVtx3dL = t->j3_secVtx3dL;
                    t->jet_secVtx3deL = t->j3_secVtx3deL;
                    t->jet_secVtxM = t->j3_secVtxM;
             		t->jet_mt = pow(t->j3_e,2)-pow(t->j3_pt*(1+sinh(t->j3_eta)),2) > 0 ? sqrt(pow(t->j3_e,2)-pow(t->j3_pt*(1+sinh(t->j3_eta)),2)) : sqrt(pow(t->j3_pt*(1+sinh(t->j3_eta)),2)-pow(t->j3_e,2)) ;
                    t->jet_emfrac = t->j3_emfrac;
                    t->jet_hadfrac = t->j3_hadfrac;
                    t->jet_chadfrac = t->j3_chadfrac;
                    t->jet_nhadfrac = t->j3_nhadfrac;
                    t->jet_phofrac = t->j3_phofrac;
                    t->jet_mufrac = t->j3_mufrac;
                    t->jet_phofrac = t->j3_phofrac;
                    t->jet_JECUnc = t->j3_JECUnc;
                    t->jet_leadTrackPt = t->j3_leadTrackPt;
                    t->jet_softLeptPt = (t->j3_softLeptIdLooseMu==1 || t->j3_softLeptIdEle95==1) ? (t->j3_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j3_softLeptIdLooseMu==1 || t->j3_softLeptIdEle95==1) ? (t->j3_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j3_softLeptIdLooseMu==1 || t->j3_softLeptIdEle95==1) ? (t->j3_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j3_btagSF_M;
                t->jet_flavour = t->j3_flavour;
                t->jet_btagSFErrorUp_M = t->j3_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j3_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j3_btagEff_M;
                t->jet_btagEffError_M = t->j3_btagEffError_M;
                t->jet_cutbased_wp_level = t->j3_cutbased_wp_level;
                t->jet_simple_wp_level = t->j3_simple_wp_level;
                t->jet_full_wp_level = t->j3_full_wp_level;
                t->jet_betaStarClassic = t->j3_betaStarClassic;
                t->jet_dR2Mean = t->j3_dR2Mean;
                t->jet_nNeutrals = t->j3_nNeutrals;
                t->jet_nCharged = t->j3_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j3_pt, t->j3_eta, t->j3_phi, t->j3_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j3_jecD_e;
                t->jet_jecD_pt = t->j3_jecD_pt;
                t->jet_jecD_phi = t->j3_jecD_phi;
                t->jet_jecD_eta = t->j3_jecD_eta;
                t->jet_jecU_e = t->j3_jecU_e;
                t->jet_jecU_pt = t->j3_jecU_pt;
                t->jet_jecU_phi = t->j3_jecU_phi;
                t->jet_jecU_eta = t->j3_jecU_eta;
                t->jet_jerD_e = t->j3_jerD_e;
                t->jet_jerD_pt = t->j3_jerD_pt;
                t->jet_jerD_phi = t->j3_jerD_phi;
                t->jet_jerD_eta = t->j3_jerD_eta;
                t->jet_jerC_e = t->j3_jerC_e;
                t->jet_jerC_pt = t->j3_jerC_pt;
                t->jet_jerC_phi = t->j3_jerC_phi;
                t->jet_jerC_eta = t->j3_jerC_eta;
                t->jet_jerU_e = t->j3_jerU_e;
                t->jet_jerU_pt = t->j3_jerU_pt;
                t->jet_jerU_phi = t->j3_jerU_phi;
                t->jet_jerU_eta = t->j3_jerU_eta;
            } 
            else
            if( ijet == 3 && t->j4_pt > .01 )
            {
                t->jet_e = t->j4_e;
                t->jet_pt = t->j4_pt;
                t->jet_phi = t->j4_phi;
                t->jet_eta = t->j4_eta;
                t->jet_betaStarClassic = t->j4_betaStarClassic;
                t->jet_dR2Mean = t->j4_dR2Mean;
                t->jet_csvBtag = t->j4_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j4_secVtxPt;
                    t->jet_secVtx3dL = t->j4_secVtx3dL;
                    t->jet_secVtx3deL = t->j4_secVtx3deL;
                    t->jet_secVtxM = t->j4_secVtxM;
             		t->jet_mt = pow(t->j4_e,2)-pow(t->j4_pt*(1+sinh(t->j4_eta)),2) > 0 ? sqrt(pow(t->j4_e,2)-pow(t->j4_pt*(1+sinh(t->j4_eta)),2)) : sqrt(pow(t->j4_pt*(1+sinh(t->j4_eta)),2)-pow(t->j4_e,2)) ;
                    t->jet_emfrac = t->j4_emfrac;
                    t->jet_hadfrac = t->j4_hadfrac;
                    t->jet_chadfrac = t->j4_chadfrac;
                    t->jet_nhadfrac = t->j4_nhadfrac;
                    t->jet_phofrac = t->j4_phofrac;
                    t->jet_mufrac = t->j4_mufrac;
                    t->jet_phofrac = t->j4_phofrac;
                    t->jet_JECUnc = t->j4_JECUnc;
                    t->jet_leadTrackPt = t->j4_leadTrackPt;
                    t->jet_softLeptPt = (t->j4_softLeptIdLooseMu==1 || t->j4_softLeptIdEle95==1) ? (t->j4_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j4_softLeptIdLooseMu==1 || t->j4_softLeptIdEle95==1) ? (t->j4_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j4_softLeptIdLooseMu==1 || t->j4_softLeptIdEle95==1) ? (t->j4_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j4_btagSF_M;
                t->jet_flavour = t->j4_flavour;
                t->jet_btagSFErrorUp_M = t->j4_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j4_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j4_btagEff_M;
                t->jet_btagEffError_M = t->j4_btagEffError_M;
                t->jet_cutbased_wp_level = t->j4_cutbased_wp_level;
                t->jet_simple_wp_level = t->j4_simple_wp_level;
                t->jet_full_wp_level = t->j4_full_wp_level;
                t->jet_betaStarClassic = t->j4_betaStarClassic;
                t->jet_dR2Mean = t->j4_dR2Mean;
                t->jet_nNeutrals = t->j4_nNeutrals;
                t->jet_nCharged = t->j4_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j4_pt, t->j4_eta, t->j4_phi, t->j4_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j4_jecD_e;
                t->jet_jecD_pt = t->j4_jecD_pt;
                t->jet_jecD_phi = t->j4_jecD_phi;
                t->jet_jecD_eta = t->j4_jecD_eta;
                t->jet_jecU_e = t->j4_jecU_e;
                t->jet_jecU_pt = t->j4_jecU_pt;
                t->jet_jecU_phi = t->j4_jecU_phi;
                t->jet_jecU_eta = t->j4_jecU_eta;
                t->jet_jerD_e = t->j4_jerD_e;
                t->jet_jerD_pt = t->j4_jerD_pt;
                t->jet_jerD_phi = t->j4_jerD_phi;
                t->jet_jerD_eta = t->j4_jerD_eta;
                t->jet_jerC_e = t->j4_jerC_e;
                t->jet_jerC_pt = t->j4_jerC_pt;
                t->jet_jerC_phi = t->j4_jerC_phi;
                t->jet_jerC_eta = t->j4_jerC_eta;
                t->jet_jerU_e = t->j4_jerU_e;
                t->jet_jerU_pt = t->j4_jerU_pt;
                t->jet_jerU_phi = t->j4_jerU_phi;
                t->jet_jerU_eta = t->j4_jerU_eta;
            } 
            else
            if( ijet == 4 && t->j5_pt > .01 )
            {
                t->jet_e = t->j5_e;
                t->jet_pt = t->j5_pt;
                t->jet_phi = t->j5_phi;
                t->jet_eta = t->j5_eta;
                t->jet_betaStarClassic = t->j5_betaStarClassic;
                t->jet_dR2Mean = t->j5_dR2Mean;
                t->jet_csvBtag = t->j5_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j5_secVtxPt;
                    t->jet_secVtx3dL = t->j5_secVtx3dL;
                    t->jet_secVtx3deL = t->j5_secVtx3deL;
                    t->jet_secVtxM = t->j5_secVtxM;
             		t->jet_mt = pow(t->j5_e,2)-pow(t->j5_pt*(1+sinh(t->j5_eta)),2) > 0 ? sqrt(pow(t->j5_e,2)-pow(t->j5_pt*(1+sinh(t->j5_eta)),2)) : sqrt(pow(t->j5_pt*(1+sinh(t->j5_eta)),2)-pow(t->j5_e,2)) ;
                    t->jet_emfrac = t->j5_emfrac;
                    t->jet_hadfrac = t->j5_hadfrac;
                    t->jet_chadfrac = t->j5_chadfrac;
                    t->jet_nhadfrac = t->j5_nhadfrac;
                    t->jet_phofrac = t->j5_phofrac;
                    t->jet_mufrac = t->j5_mufrac;
                    t->jet_phofrac = t->j5_phofrac;
                    t->jet_JECUnc = t->j5_JECUnc;
                    t->jet_leadTrackPt = t->j5_leadTrackPt;
                    t->jet_softLeptPt = (t->j5_softLeptIdLooseMu==1 || t->j5_softLeptIdEle95==1) ? (t->j5_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j5_softLeptIdLooseMu==1 || t->j5_softLeptIdEle95==1) ? (t->j5_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j5_softLeptIdLooseMu==1 || t->j5_softLeptIdEle95==1) ? (t->j5_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j5_btagSF_M;
                t->jet_flavour = t->j5_flavour;
                t->jet_btagSFErrorUp_M = t->j5_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j5_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j5_btagEff_M;
                t->jet_btagEffError_M = t->j5_btagEffError_M;
                t->jet_cutbased_wp_level = t->j5_cutbased_wp_level;
                t->jet_simple_wp_level = t->j5_simple_wp_level;
                t->jet_full_wp_level = t->j5_full_wp_level;
                t->jet_betaStarClassic = t->j5_betaStarClassic;
                t->jet_dR2Mean = t->j5_dR2Mean;
                t->jet_nNeutrals = t->j5_nNeutrals;
                t->jet_nCharged = t->j5_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j5_pt, t->j5_eta, t->j5_phi, t->j5_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j5_jecD_e;
                t->jet_jecD_pt = t->j5_jecD_pt;
                t->jet_jecD_phi = t->j5_jecD_phi;
                t->jet_jecD_eta = t->j5_jecD_eta;
                t->jet_jecU_e = t->j5_jecU_e;
                t->jet_jecU_pt = t->j5_jecU_pt;
                t->jet_jecU_phi = t->j5_jecU_phi;
                t->jet_jecU_eta = t->j5_jecU_eta;
                t->jet_jerD_e = t->j5_jerD_e;
                t->jet_jerD_pt = t->j5_jerD_pt;
                t->jet_jerD_phi = t->j5_jerD_phi;
                t->jet_jerD_eta = t->j5_jerD_eta;
                t->jet_jerC_e = t->j5_jerC_e;
                t->jet_jerC_pt = t->j5_jerC_pt;
                t->jet_jerC_phi = t->j5_jerC_phi;
                t->jet_jerC_eta = t->j5_jerC_eta;
                t->jet_jerU_e = t->j5_jerU_e;
                t->jet_jerU_pt = t->j5_jerU_pt;
                t->jet_jerU_phi = t->j5_jerU_phi;
                t->jet_jerU_eta = t->j5_jerU_eta;
            } 
            else
            if( ijet == 5 && t->j6_pt > .01 )
            {
                t->jet_e = t->j6_e;
                t->jet_pt = t->j6_pt;
                t->jet_phi = t->j6_phi;
                t->jet_eta = t->j6_eta;
                t->jet_betaStarClassic = t->j6_betaStarClassic;
                t->jet_dR2Mean = t->j6_dR2Mean;
                t->jet_csvBtag = t->j6_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j6_secVtxPt;
                    t->jet_secVtx3dL = t->j6_secVtx3dL;
                    t->jet_secVtx3deL = t->j6_secVtx3deL;
                    t->jet_secVtxM = t->j6_secVtxM;
             		t->jet_mt = pow(t->j6_e,2)-pow(t->j6_pt*(1+sinh(t->j6_eta)),2) > 0 ? sqrt(pow(t->j6_e,2)-pow(t->j6_pt*(1+sinh(t->j6_eta)),2)) : sqrt(pow(t->j6_pt*(1+sinh(t->j6_eta)),2)-pow(t->j6_e,2)) ;
                    t->jet_emfrac = t->j6_emfrac;
                    t->jet_hadfrac = t->j6_hadfrac;
                    t->jet_chadfrac = t->j6_chadfrac;
                    t->jet_nhadfrac = t->j6_nhadfrac;
                    t->jet_phofrac = t->j6_phofrac;
                    t->jet_mufrac = t->j6_mufrac;
                    t->jet_phofrac = t->j6_phofrac;
                    t->jet_JECUnc = t->j6_JECUnc;
                    t->jet_leadTrackPt = t->j6_leadTrackPt;
                    t->jet_softLeptPt = (t->j6_softLeptIdLooseMu==1 || t->j6_softLeptIdEle95==1) ? (t->j6_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j6_softLeptIdLooseMu==1 || t->j6_softLeptIdEle95==1) ? (t->j6_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j6_softLeptIdLooseMu==1 || t->j6_softLeptIdEle95==1) ? (t->j6_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j6_btagSF_M;
                t->jet_flavour = t->j6_flavour;
                t->jet_btagSFErrorUp_M = t->j6_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j6_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j6_btagEff_M;
                t->jet_btagEffError_M = t->j6_btagEffError_M;
                t->jet_cutbased_wp_level = t->j6_cutbased_wp_level;
                t->jet_simple_wp_level = t->j6_simple_wp_level;
                t->jet_full_wp_level = t->j6_full_wp_level;
                t->jet_betaStarClassic = t->j6_betaStarClassic;
                t->jet_dR2Mean = t->j6_dR2Mean;
                t->jet_nNeutrals = t->j6_nNeutrals;
                t->jet_nCharged = t->j6_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j6_pt, t->j6_eta, t->j6_phi, t->j6_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j6_jecD_e;
                t->jet_jecD_pt = t->j6_jecD_pt;
                t->jet_jecD_phi = t->j6_jecD_phi;
                t->jet_jecD_eta = t->j6_jecD_eta;
                t->jet_jecU_e = t->j6_jecU_e;
                t->jet_jecU_pt = t->j6_jecU_pt;
                t->jet_jecU_phi = t->j6_jecU_phi;
                t->jet_jecU_eta = t->j6_jecU_eta;
                t->jet_jerD_e = t->j6_jerD_e;
                t->jet_jerD_pt = t->j6_jerD_pt;
                t->jet_jerD_phi = t->j6_jerD_phi;
                t->jet_jerD_eta = t->j6_jerD_eta;
                t->jet_jerC_e = t->j6_jerC_e;
                t->jet_jerC_pt = t->j6_jerC_pt;
                t->jet_jerC_phi = t->j6_jerC_phi;
                t->jet_jerC_eta = t->j6_jerC_eta;
                t->jet_jerU_e = t->j6_jerU_e;
                t->jet_jerU_pt = t->j6_jerU_pt;
                t->jet_jerU_phi = t->j6_jerU_phi;
                t->jet_jerU_eta = t->j6_jerU_eta;
            } 
            else
            if( ijet == 6 && t->j7_pt > .01 )
            {
                t->jet_e = t->j7_e;
                t->jet_pt = t->j7_pt;
                t->jet_phi = t->j7_phi;
                t->jet_eta = t->j7_eta;
                t->jet_betaStarClassic = t->j7_betaStarClassic;
                t->jet_dR2Mean = t->j7_dR2Mean;
                t->jet_csvBtag = t->j7_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j7_secVtxPt;
                    t->jet_secVtx3dL = t->j7_secVtx3dL;
                    t->jet_secVtx3deL = t->j7_secVtx3deL;
                    t->jet_secVtxM = t->j7_secVtxM;
             		t->jet_mt = pow(t->j7_e,2)-pow(t->j7_pt*(1+sinh(t->j7_eta)),2) > 0 ? sqrt(pow(t->j7_e,2)-pow(t->j7_pt*(1+sinh(t->j7_eta)),2)) : sqrt(pow(t->j7_pt*(1+sinh(t->j7_eta)),2)-pow(t->j7_e,2)) ;
                    t->jet_emfrac = t->j7_emfrac;
                    t->jet_hadfrac = t->j7_hadfrac;
                    t->jet_chadfrac = t->j7_chadfrac;
                    t->jet_nhadfrac = t->j7_nhadfrac;
                    t->jet_phofrac = t->j7_phofrac;
                    t->jet_mufrac = t->j7_mufrac;
                    t->jet_phofrac = t->j7_phofrac;
                    t->jet_JECUnc = t->j7_JECUnc;
                    t->jet_leadTrackPt = t->j7_leadTrackPt;
                    t->jet_softLeptPt = (t->j7_softLeptIdLooseMu==1 || t->j7_softLeptIdEle95==1) ? (t->j7_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j7_softLeptIdLooseMu==1 || t->j7_softLeptIdEle95==1) ? (t->j7_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j7_softLeptIdLooseMu==1 || t->j7_softLeptIdEle95==1) ? (t->j7_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j7_btagSF_M;
                t->jet_flavour = t->j7_flavour;
                t->jet_btagSFErrorUp_M = t->j7_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j7_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j7_btagEff_M;
                t->jet_btagEffError_M = t->j7_btagEffError_M;
                t->jet_cutbased_wp_level = t->j7_cutbased_wp_level;
                t->jet_simple_wp_level = t->j7_simple_wp_level;
                t->jet_full_wp_level = t->j7_full_wp_level;
                t->jet_betaStarClassic = t->j7_betaStarClassic;
                t->jet_dR2Mean = t->j7_dR2Mean;
                t->jet_nNeutrals = t->j7_nNeutrals;
                t->jet_nCharged = t->j7_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j7_pt, t->j7_eta, t->j7_phi, t->j7_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j7_jecD_e;
                t->jet_jecD_pt = t->j7_jecD_pt;
                t->jet_jecD_phi = t->j7_jecD_phi;
                t->jet_jecD_eta = t->j7_jecD_eta;
                t->jet_jecU_e = t->j7_jecU_e;
                t->jet_jecU_pt = t->j7_jecU_pt;
                t->jet_jecU_phi = t->j7_jecU_phi;
                t->jet_jecU_eta = t->j7_jecU_eta;
                t->jet_jerD_e = t->j7_jerD_e;
                t->jet_jerD_pt = t->j7_jerD_pt;
                t->jet_jerD_phi = t->j7_jerD_phi;
                t->jet_jerD_eta = t->j7_jerD_eta;
                t->jet_jerC_e = t->j7_jerC_e;
                t->jet_jerC_pt = t->j7_jerC_pt;
                t->jet_jerC_phi = t->j7_jerC_phi;
                t->jet_jerC_eta = t->j7_jerC_eta;
                t->jet_jerU_e = t->j7_jerU_e;
                t->jet_jerU_pt = t->j7_jerU_pt;
                t->jet_jerU_phi = t->j7_jerU_phi;
                t->jet_jerU_eta = t->j7_jerU_eta;
            } 
            else
            if( ijet == 7 && t->j8_pt > .01 )
            {
                t->jet_e = t->j8_e;
                t->jet_pt = t->j8_pt;
                t->jet_phi = t->j8_phi;
                t->jet_eta = t->j8_eta;
                t->jet_betaStarClassic = t->j8_betaStarClassic;
                t->jet_dR2Mean = t->j8_dR2Mean;
                t->jet_csvBtag = t->j8_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j8_secVtxPt;
                    t->jet_secVtx3dL = t->j8_secVtx3dL;
                    t->jet_secVtx3deL = t->j8_secVtx3deL;
                    t->jet_secVtxM = t->j8_secVtxM;
             		t->jet_mt = pow(t->j8_e,2)-pow(t->j8_pt*(1+sinh(t->j8_eta)),2) > 0 ? sqrt(pow(t->j8_e,2)-pow(t->j8_pt*(1+sinh(t->j8_eta)),2)) : sqrt(pow(t->j8_pt*(1+sinh(t->j8_eta)),2)-pow(t->j8_e,2)) ;
                    t->jet_emfrac = t->j8_emfrac;
                    t->jet_hadfrac = t->j8_hadfrac;
                    t->jet_chadfrac = t->j8_chadfrac;
                    t->jet_nhadfrac = t->j8_nhadfrac;
                    t->jet_phofrac = t->j8_phofrac;
                    t->jet_mufrac = t->j8_mufrac;
                    t->jet_phofrac = t->j8_phofrac;
                    t->jet_JECUnc = t->j8_JECUnc;
                    t->jet_leadTrackPt = t->j8_leadTrackPt;
                    t->jet_softLeptPt = (t->j8_softLeptIdLooseMu==1 || t->j8_softLeptIdEle95==1) ? (t->j8_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j8_softLeptIdLooseMu==1 || t->j8_softLeptIdEle95==1) ? (t->j8_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j8_softLeptIdLooseMu==1 || t->j8_softLeptIdEle95==1) ? (t->j8_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j8_btagSF_M;
                t->jet_flavour = t->j8_flavour;
                t->jet_btagSFErrorUp_M = t->j8_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j8_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j8_btagEff_M;
                t->jet_btagEffError_M = t->j8_btagEffError_M;
                t->jet_cutbased_wp_level = t->j8_cutbased_wp_level;
                t->jet_simple_wp_level = t->j8_simple_wp_level;
                t->jet_full_wp_level = t->j8_full_wp_level;
                t->jet_betaStarClassic = t->j8_betaStarClassic;
                t->jet_dR2Mean = t->j8_dR2Mean;
                t->jet_nNeutrals = t->j8_nNeutrals;
                t->jet_nCharged = t->j8_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j8_pt, t->j8_eta, t->j8_phi, t->j8_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j8_jecD_e;
                t->jet_jecD_pt = t->j8_jecD_pt;
                t->jet_jecD_phi = t->j8_jecD_phi;
                t->jet_jecD_eta = t->j8_jecD_eta;
                t->jet_jecU_e = t->j8_jecU_e;
                t->jet_jecU_pt = t->j8_jecU_pt;
                t->jet_jecU_phi = t->j8_jecU_phi;
                t->jet_jecU_eta = t->j8_jecU_eta;
                t->jet_jerD_e = t->j8_jerD_e;
                t->jet_jerD_pt = t->j8_jerD_pt;
                t->jet_jerD_phi = t->j8_jerD_phi;
                t->jet_jerD_eta = t->j8_jerD_eta;
                t->jet_jerC_e = t->j8_jerC_e;
                t->jet_jerC_pt = t->j8_jerC_pt;
                t->jet_jerC_phi = t->j8_jerC_phi;
                t->jet_jerC_eta = t->j8_jerC_eta;
                t->jet_jerU_e = t->j8_jerU_e;
                t->jet_jerU_pt = t->j8_jerU_pt;
                t->jet_jerU_phi = t->j8_jerU_phi;
                t->jet_jerU_eta = t->j8_jerU_eta;
            } 
            else
            if( ijet == 8 && t->j9_pt > .01 )
            {
                t->jet_e = t->j9_e;
                t->jet_pt = t->j9_pt;
                t->jet_phi = t->j9_phi;
                t->jet_eta = t->j9_eta;
                t->jet_betaStarClassic = t->j9_betaStarClassic;
                t->jet_dR2Mean = t->j9_dR2Mean;
                t->jet_csvBtag = t->j9_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j9_secVtxPt;
                    t->jet_secVtx3dL = t->j9_secVtx3dL;
                    t->jet_secVtx3deL = t->j9_secVtx3deL;
                    t->jet_secVtxM = t->j9_secVtxM;
             		t->jet_mt = pow(t->j9_e,2)-pow(t->j9_pt*(1+sinh(t->j9_eta)),2) > 0 ? sqrt(pow(t->j9_e,2)-pow(t->j9_pt*(1+sinh(t->j9_eta)),2)) : sqrt(pow(t->j9_pt*(1+sinh(t->j9_eta)),2)-pow(t->j9_e,2)) ;
                    t->jet_emfrac = t->j9_emfrac;
                    t->jet_hadfrac = t->j9_hadfrac;
                    t->jet_chadfrac = t->j9_chadfrac;
                    t->jet_nhadfrac = t->j9_nhadfrac;
                    t->jet_phofrac = t->j9_phofrac;
                    t->jet_mufrac = t->j9_mufrac;
                    t->jet_phofrac = t->j9_phofrac;
                    t->jet_JECUnc = t->j9_JECUnc;
                    t->jet_leadTrackPt = t->j9_leadTrackPt;
                    t->jet_softLeptPt = (t->j9_softLeptIdLooseMu==1 || t->j9_softLeptIdEle95==1) ? (t->j9_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j9_softLeptIdLooseMu==1 || t->j9_softLeptIdEle95==1) ? (t->j9_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j9_softLeptIdLooseMu==1 || t->j9_softLeptIdEle95==1) ? (t->j9_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j9_btagSF_M;
                t->jet_flavour = t->j9_flavour;
                t->jet_btagSFErrorUp_M = t->j9_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j9_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j9_btagEff_M;
                t->jet_btagEffError_M = t->j9_btagEffError_M;
                t->jet_cutbased_wp_level = t->j9_cutbased_wp_level;
                t->jet_simple_wp_level = t->j9_simple_wp_level;
                t->jet_full_wp_level = t->j9_full_wp_level;
                t->jet_betaStarClassic = t->j9_betaStarClassic;
                t->jet_dR2Mean = t->j9_dR2Mean;
                t->jet_nNeutrals = t->j9_nNeutrals;
                t->jet_nCharged = t->j9_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j9_pt, t->j9_eta, t->j9_phi, t->j9_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j9_jecD_e;
                t->jet_jecD_pt = t->j9_jecD_pt;
                t->jet_jecD_phi = t->j9_jecD_phi;
                t->jet_jecD_eta = t->j9_jecD_eta;
                t->jet_jecU_e = t->j9_jecU_e;
                t->jet_jecU_pt = t->j9_jecU_pt;
                t->jet_jecU_phi = t->j9_jecU_phi;
                t->jet_jecU_eta = t->j9_jecU_eta;
                t->jet_jerD_e = t->j9_jerD_e;
                t->jet_jerD_pt = t->j9_jerD_pt;
                t->jet_jerD_phi = t->j9_jerD_phi;
                t->jet_jerD_eta = t->j9_jerD_eta;
                t->jet_jerC_e = t->j9_jerC_e;
                t->jet_jerC_pt = t->j9_jerC_pt;
                t->jet_jerC_phi = t->j9_jerC_phi;
                t->jet_jerC_eta = t->j9_jerC_eta;
                t->jet_jerU_e = t->j9_jerU_e;
                t->jet_jerU_pt = t->j9_jerU_pt;
                t->jet_jerU_phi = t->j9_jerU_phi;
                t->jet_jerU_eta = t->j9_jerU_eta;
            } 
            else
            if( ijet == 9 && t->j10_pt > .01 )
            {
                t->jet_e = t->j10_e;
                t->jet_pt = t->j10_pt;
                t->jet_phi = t->j10_phi;
                t->jet_eta = t->j10_eta;
                t->jet_betaStarClassic = t->j10_betaStarClassic;
                t->jet_dR2Mean = t->j10_dR2Mean;
                t->jet_csvBtag = t->j10_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j10_secVtxPt;
                    t->jet_secVtx3dL = t->j10_secVtx3dL;
                    t->jet_secVtx3deL = t->j10_secVtx3deL;
                    t->jet_secVtxM = t->j10_secVtxM;
             		t->jet_mt = pow(t->j10_e,2)-pow(t->j10_pt*(1+sinh(t->j10_eta)),2) > 0 ? sqrt(pow(t->j10_e,2)-pow(t->j10_pt*(1+sinh(t->j10_eta)),2)) : sqrt(pow(t->j10_pt*(1+sinh(t->j10_eta)),2)-pow(t->j10_e,2)) ;
                    t->jet_emfrac = t->j10_emfrac;
                    t->jet_hadfrac = t->j10_hadfrac;
                    t->jet_chadfrac = t->j10_chadfrac;
                    t->jet_nhadfrac = t->j10_nhadfrac;
                    t->jet_phofrac = t->j10_phofrac;
                    t->jet_mufrac = t->j10_mufrac;
                    t->jet_phofrac = t->j10_phofrac;
                    t->jet_JECUnc = t->j10_JECUnc;
                    t->jet_leadTrackPt = t->j10_leadTrackPt;
                    t->jet_softLeptPt = (t->j10_softLeptIdLooseMu==1 || t->j10_softLeptIdEle95==1) ? (t->j10_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j10_softLeptIdLooseMu==1 || t->j10_softLeptIdEle95==1) ? (t->j10_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j10_softLeptIdLooseMu==1 || t->j10_softLeptIdEle95==1) ? (t->j10_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j10_btagSF_M;
                t->jet_flavour = t->j10_flavour;
                t->jet_btagSFErrorUp_M = t->j10_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j10_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j10_btagEff_M;
                t->jet_btagEffError_M = t->j10_btagEffError_M;
                t->jet_cutbased_wp_level = t->j10_cutbased_wp_level;
                t->jet_simple_wp_level = t->j10_simple_wp_level;
                t->jet_full_wp_level = t->j10_full_wp_level;
                t->jet_betaStarClassic = t->j10_betaStarClassic;
                t->jet_dR2Mean = t->j10_dR2Mean;
                t->jet_nNeutrals = t->j10_nNeutrals;
                t->jet_nCharged = t->j10_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j10_pt, t->j10_eta, t->j10_phi, t->j10_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j10_jecD_e;
                t->jet_jecD_pt = t->j10_jecD_pt;
                t->jet_jecD_phi = t->j10_jecD_phi;
                t->jet_jecD_eta = t->j10_jecD_eta;
                t->jet_jecU_e = t->j10_jecU_e;
                t->jet_jecU_pt = t->j10_jecU_pt;
                t->jet_jecU_phi = t->j10_jecU_phi;
                t->jet_jecU_eta = t->j10_jecU_eta;
                t->jet_jerD_e = t->j10_jerD_e;
                t->jet_jerD_pt = t->j10_jerD_pt;
                t->jet_jerD_phi = t->j10_jerD_phi;
                t->jet_jerD_eta = t->j10_jerD_eta;
                t->jet_jerC_e = t->j10_jerC_e;
                t->jet_jerC_pt = t->j10_jerC_pt;
                t->jet_jerC_phi = t->j10_jerC_phi;
                t->jet_jerC_eta = t->j10_jerC_eta;
                t->jet_jerU_e = t->j10_jerU_e;
                t->jet_jerU_pt = t->j10_jerU_pt;
                t->jet_jerU_phi = t->j10_jerU_phi;
                t->jet_jerU_eta = t->j10_jerU_eta;
            } 
            else
            if( ijet == 10 && t->j11_pt > .01 )
            {
                t->jet_e = t->j11_e;
                t->jet_pt = t->j11_pt;
                t->jet_phi = t->j11_phi;
                t->jet_eta = t->j11_eta;
                t->jet_betaStarClassic = t->j11_betaStarClassic;
                t->jet_dR2Mean = t->j11_dR2Mean;
                t->jet_csvBtag = t->j11_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j11_secVtxPt;
                    t->jet_secVtx3dL = t->j11_secVtx3dL;
                    t->jet_secVtx3deL = t->j11_secVtx3deL;
                    t->jet_secVtxM = t->j11_secVtxM;
             		t->jet_mt = pow(t->j11_e,2)-pow(t->j11_pt*(1+sinh(t->j11_eta)),2) > 0 ? sqrt(pow(t->j11_e,2)-pow(t->j11_pt*(1+sinh(t->j11_eta)),2)) : sqrt(pow(t->j11_pt*(1+sinh(t->j11_eta)),2)-pow(t->j11_e,2)) ;
                    t->jet_emfrac = t->j11_emfrac;
                    t->jet_hadfrac = t->j11_hadfrac;
                    t->jet_chadfrac = t->j11_chadfrac;
                    t->jet_nhadfrac = t->j11_nhadfrac;
                    t->jet_phofrac = t->j11_phofrac;
                    t->jet_mufrac = t->j11_mufrac;
                    t->jet_phofrac = t->j11_phofrac;
                    t->jet_JECUnc = t->j11_JECUnc;
                    t->jet_leadTrackPt = t->j11_leadTrackPt;
                    t->jet_softLeptPt = (t->j11_softLeptIdLooseMu==1 || t->j11_softLeptIdEle95==1) ? (t->j11_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j11_softLeptIdLooseMu==1 || t->j11_softLeptIdEle95==1) ? (t->j11_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j11_softLeptIdLooseMu==1 || t->j11_softLeptIdEle95==1) ? (t->j11_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j11_btagSF_M;
                t->jet_flavour = t->j11_flavour;
                t->jet_btagSFErrorUp_M = t->j11_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j11_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j11_btagEff_M;
                t->jet_btagEffError_M = t->j11_btagEffError_M;
                t->jet_cutbased_wp_level = t->j11_cutbased_wp_level;
                t->jet_simple_wp_level = t->j11_simple_wp_level;
                t->jet_full_wp_level = t->j11_full_wp_level;
                t->jet_betaStarClassic = t->j11_betaStarClassic;
                t->jet_dR2Mean = t->j11_dR2Mean;
                t->jet_nNeutrals = t->j11_nNeutrals;
                t->jet_nCharged = t->j11_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j11_pt, t->j11_eta, t->j11_phi, t->j11_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j11_jecD_e;
                t->jet_jecD_pt = t->j11_jecD_pt;
                t->jet_jecD_phi = t->j11_jecD_phi;
                t->jet_jecD_eta = t->j11_jecD_eta;
                t->jet_jecU_e = t->j11_jecU_e;
                t->jet_jecU_pt = t->j11_jecU_pt;
                t->jet_jecU_phi = t->j11_jecU_phi;
                t->jet_jecU_eta = t->j11_jecU_eta;
                t->jet_jerD_e = t->j11_jerD_e;
                t->jet_jerD_pt = t->j11_jerD_pt;
                t->jet_jerD_phi = t->j11_jerD_phi;
                t->jet_jerD_eta = t->j11_jerD_eta;
                t->jet_jerC_e = t->j11_jerC_e;
                t->jet_jerC_pt = t->j11_jerC_pt;
                t->jet_jerC_phi = t->j11_jerC_phi;
                t->jet_jerC_eta = t->j11_jerC_eta;
                t->jet_jerU_e = t->j11_jerU_e;
                t->jet_jerU_pt = t->j11_jerU_pt;
                t->jet_jerU_phi = t->j11_jerU_phi;
                t->jet_jerU_eta = t->j11_jerU_eta;
            } 
            else
            if( ijet == 11 && t->j12_pt > .01 )
            {
                t->jet_e = t->j12_e;
                t->jet_pt = t->j12_pt;
                t->jet_phi = t->j12_phi;
                t->jet_eta = t->j12_eta;
                t->jet_betaStarClassic = t->j12_betaStarClassic;
                t->jet_dR2Mean = t->j12_dR2Mean;
                t->jet_csvBtag = t->j12_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j12_secVtxPt;
                    t->jet_secVtx3dL = t->j12_secVtx3dL;
                    t->jet_secVtx3deL = t->j12_secVtx3deL;
                    t->jet_secVtxM = t->j12_secVtxM;
             		t->jet_mt = pow(t->j12_e,2)-pow(t->j12_pt*(1+sinh(t->j12_eta)),2) > 0 ? sqrt(pow(t->j12_e,2)-pow(t->j12_pt*(1+sinh(t->j12_eta)),2)) : sqrt(pow(t->j12_pt*(1+sinh(t->j12_eta)),2)-pow(t->j12_e,2)) ;
                    t->jet_emfrac = t->j12_emfrac;
                    t->jet_hadfrac = t->j12_hadfrac;
                    t->jet_chadfrac = t->j12_chadfrac;
                    t->jet_nhadfrac = t->j12_nhadfrac;
                    t->jet_phofrac = t->j12_phofrac;
                    t->jet_mufrac = t->j12_mufrac;
                    t->jet_phofrac = t->j12_phofrac;
                    t->jet_JECUnc = t->j12_JECUnc;
                    t->jet_leadTrackPt = t->j12_leadTrackPt;
                    t->jet_softLeptPt = (t->j12_softLeptIdLooseMu==1 || t->j12_softLeptIdEle95==1) ? (t->j12_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j12_softLeptIdLooseMu==1 || t->j12_softLeptIdEle95==1) ? (t->j12_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j12_softLeptIdLooseMu==1 || t->j12_softLeptIdEle95==1) ? (t->j12_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j12_btagSF_M;
                t->jet_flavour = t->j12_flavour;
                t->jet_btagSFErrorUp_M = t->j12_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j12_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j12_btagEff_M;
                t->jet_btagEffError_M = t->j12_btagEffError_M;
                t->jet_cutbased_wp_level = t->j12_cutbased_wp_level;
                t->jet_simple_wp_level = t->j12_simple_wp_level;
                t->jet_full_wp_level = t->j12_full_wp_level;
                t->jet_betaStarClassic = t->j12_betaStarClassic;
                t->jet_dR2Mean = t->j12_dR2Mean;
                t->jet_nNeutrals = t->j12_nNeutrals;
                t->jet_nCharged = t->j12_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j12_pt, t->j12_eta, t->j12_phi, t->j12_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j12_jecD_e;
                t->jet_jecD_pt = t->j12_jecD_pt;
                t->jet_jecD_phi = t->j12_jecD_phi;
                t->jet_jecD_eta = t->j12_jecD_eta;
                t->jet_jecU_e = t->j12_jecU_e;
                t->jet_jecU_pt = t->j12_jecU_pt;
                t->jet_jecU_phi = t->j12_jecU_phi;
                t->jet_jecU_eta = t->j12_jecU_eta;
                t->jet_jerD_e = t->j12_jerD_e;
                t->jet_jerD_pt = t->j12_jerD_pt;
                t->jet_jerD_phi = t->j12_jerD_phi;
                t->jet_jerD_eta = t->j12_jerD_eta;
                t->jet_jerC_e = t->j12_jerC_e;
                t->jet_jerC_pt = t->j12_jerC_pt;
                t->jet_jerC_phi = t->j12_jerC_phi;
                t->jet_jerC_eta = t->j12_jerC_eta;
                t->jet_jerU_e = t->j12_jerU_e;
                t->jet_jerU_pt = t->j12_jerU_pt;
                t->jet_jerU_phi = t->j12_jerU_phi;
                t->jet_jerU_eta = t->j12_jerU_eta;
            } 
            else
            if( ijet == 12 && t->j13_pt > .01 )
            {
                t->jet_e = t->j13_e;
                t->jet_pt = t->j13_pt;
                t->jet_phi = t->j13_phi;
                t->jet_eta = t->j13_eta;
                t->jet_betaStarClassic = t->j13_betaStarClassic;
                t->jet_dR2Mean = t->j13_dR2Mean;
                t->jet_csvBtag = t->j13_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j13_secVtxPt;
                    t->jet_secVtx3dL = t->j13_secVtx3dL;
                    t->jet_secVtx3deL = t->j13_secVtx3deL;
                    t->jet_secVtxM = t->j13_secVtxM;
             		t->jet_mt = pow(t->j13_e,2)-pow(t->j13_pt*(1+sinh(t->j13_eta)),2) > 0 ? sqrt(pow(t->j13_e,2)-pow(t->j13_pt*(1+sinh(t->j13_eta)),2)) : sqrt(pow(t->j13_pt*(1+sinh(t->j13_eta)),2)-pow(t->j13_e,2)) ;
                    t->jet_emfrac = t->j13_emfrac;
                    t->jet_hadfrac = t->j13_hadfrac;
                    t->jet_chadfrac = t->j13_chadfrac;
                    t->jet_nhadfrac = t->j13_nhadfrac;
                    t->jet_phofrac = t->j13_phofrac;
                    t->jet_mufrac = t->j13_mufrac;
                    t->jet_phofrac = t->j13_phofrac;
                    t->jet_JECUnc = t->j13_JECUnc;
                    t->jet_leadTrackPt = t->j13_leadTrackPt;
                    t->jet_softLeptPt = (t->j13_softLeptIdLooseMu==1 || t->j13_softLeptIdEle95==1) ? (t->j13_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j13_softLeptIdLooseMu==1 || t->j13_softLeptIdEle95==1) ? (t->j13_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j13_softLeptIdLooseMu==1 || t->j13_softLeptIdEle95==1) ? (t->j13_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j13_btagSF_M;
                t->jet_flavour = t->j13_flavour;
                t->jet_btagSFErrorUp_M = t->j13_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j13_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j13_btagEff_M;
                t->jet_btagEffError_M = t->j13_btagEffError_M;
                t->jet_cutbased_wp_level = t->j13_cutbased_wp_level;
                t->jet_simple_wp_level = t->j13_simple_wp_level;
                t->jet_full_wp_level = t->j13_full_wp_level;
                t->jet_betaStarClassic = t->j13_betaStarClassic;
                t->jet_dR2Mean = t->j13_dR2Mean;
                t->jet_nNeutrals = t->j13_nNeutrals;
                t->jet_nCharged = t->j13_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j13_pt, t->j13_eta, t->j13_phi, t->j13_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j13_jecD_e;
                t->jet_jecD_pt = t->j13_jecD_pt;
                t->jet_jecD_phi = t->j13_jecD_phi;
                t->jet_jecD_eta = t->j13_jecD_eta;
                t->jet_jecU_e = t->j13_jecU_e;
                t->jet_jecU_pt = t->j13_jecU_pt;
                t->jet_jecU_phi = t->j13_jecU_phi;
                t->jet_jecU_eta = t->j13_jecU_eta;
                t->jet_jerD_e = t->j13_jerD_e;
                t->jet_jerD_pt = t->j13_jerD_pt;
                t->jet_jerD_phi = t->j13_jerD_phi;
                t->jet_jerD_eta = t->j13_jerD_eta;
                t->jet_jerC_e = t->j13_jerC_e;
                t->jet_jerC_pt = t->j13_jerC_pt;
                t->jet_jerC_phi = t->j13_jerC_phi;
                t->jet_jerC_eta = t->j13_jerC_eta;
                t->jet_jerU_e = t->j13_jerU_e;
                t->jet_jerU_pt = t->j13_jerU_pt;
                t->jet_jerU_phi = t->j13_jerU_phi;
                t->jet_jerU_eta = t->j13_jerU_eta;
            } 
            else
            if( ijet == 13 && t->j14_pt > .01 )
            {
                t->jet_e = t->j14_e;
                t->jet_pt = t->j14_pt;
                t->jet_phi = t->j14_phi;
                t->jet_eta = t->j14_eta;
                t->jet_betaStarClassic = t->j14_betaStarClassic;
                t->jet_dR2Mean = t->j14_dR2Mean;
                t->jet_csvBtag = t->j14_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j14_secVtxPt;
                    t->jet_secVtx3dL = t->j14_secVtx3dL;
                    t->jet_secVtx3deL = t->j14_secVtx3deL;
                    t->jet_secVtxM = t->j14_secVtxM;
            		t->jet_mt = pow(t->j14_e,2)-pow(t->j14_pt*(1+sinh(t->j14_eta)),2) > 0 ? sqrt(pow(t->j14_e,2)-pow(t->j14_pt*(1+sinh(t->j14_eta)),2)) : sqrt(pow(t->j14_pt*(1+sinh(t->j14_eta)),2)-pow(t->j14_e,2)) ;
                    t->jet_emfrac = t->j14_emfrac;
                    t->jet_hadfrac = t->j14_hadfrac;
                    t->jet_chadfrac = t->j14_chadfrac;
                    t->jet_nhadfrac = t->j14_nhadfrac;
                    t->jet_phofrac = t->j14_phofrac;
                    t->jet_mufrac = t->j14_mufrac;
                    t->jet_phofrac = t->j14_phofrac;
                    t->jet_JECUnc = t->j14_JECUnc;
                    t->jet_leadTrackPt = t->j14_leadTrackPt;
                    t->jet_softLeptPt = (t->j14_softLeptIdLooseMu==1 || t->j14_softLeptIdEle95==1) ? (t->j14_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j14_softLeptIdLooseMu==1 || t->j14_softLeptIdEle95==1) ? (t->j14_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j14_softLeptIdLooseMu==1 || t->j14_softLeptIdEle95==1) ? (t->j14_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j14_btagSF_M;
                t->jet_flavour = t->j14_flavour;
                t->jet_btagSFErrorUp_M = t->j14_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j14_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j14_btagEff_M;
                t->jet_btagEffError_M = t->j14_btagEffError_M;
                t->jet_cutbased_wp_level = t->j14_cutbased_wp_level;
                t->jet_simple_wp_level = t->j14_simple_wp_level;
                t->jet_full_wp_level = t->j14_full_wp_level;
                t->jet_betaStarClassic = t->j14_betaStarClassic;
                t->jet_dR2Mean = t->j14_dR2Mean;
                t->jet_nNeutrals = t->j14_nNeutrals;
                t->jet_nCharged = t->j14_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j14_pt, t->j14_eta, t->j14_phi, t->j14_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j14_jecD_e;
                t->jet_jecD_pt = t->j14_jecD_pt;
                t->jet_jecD_phi = t->j14_jecD_phi;
                t->jet_jecD_eta = t->j14_jecD_eta;
                t->jet_jecU_e = t->j14_jecU_e;
                t->jet_jecU_pt = t->j14_jecU_pt;
                t->jet_jecU_phi = t->j14_jecU_phi;
                t->jet_jecU_eta = t->j14_jecU_eta;
                t->jet_jerD_e = t->j14_jerD_e;
                t->jet_jerD_pt = t->j14_jerD_pt;
                t->jet_jerD_phi = t->j14_jerD_phi;
                t->jet_jerD_eta = t->j14_jerD_eta;
                t->jet_jerC_e = t->j14_jerC_e;
                t->jet_jerC_pt = t->j14_jerC_pt;
                t->jet_jerC_phi = t->j14_jerC_phi;
                t->jet_jerC_eta = t->j14_jerC_eta;
                t->jet_jerU_e = t->j14_jerU_e;
                t->jet_jerU_pt = t->j14_jerU_pt;
                t->jet_jerU_phi = t->j14_jerU_phi;
                t->jet_jerU_eta = t->j14_jerU_eta;
            } 
            else
            if( ijet == 14 && t->j15_pt > .01 )
            {
                t->jet_e = t->j15_e;
                t->jet_pt = t->j15_pt;
                t->jet_phi = t->j15_phi;
                t->jet_eta = t->j15_eta;
                t->jet_betaStarClassic = t->j15_betaStarClassic;
                t->jet_dR2Mean = t->j15_dR2Mean;
                t->jet_csvBtag = t->j15_csvBtag;
                if( numberOfRegressionFiles > 0 )
                {
                    t->jet_secVtxPt = t->j15_secVtxPt;
                    t->jet_secVtx3dL = t->j15_secVtx3dL;
                    t->jet_secVtx3deL = t->j15_secVtx3deL;
                    t->jet_secVtxM = t->j15_secVtxM;
             		t->jet_mt = pow(t->j15_e,2)-pow(t->j15_pt*(1+sinh(t->j15_eta)),2) > 0 ? sqrt(pow(t->j15_e,2)-pow(t->j15_pt*(1+sinh(t->j15_eta)),2)) : sqrt(pow(t->j15_pt*(1+sinh(t->j15_eta)),2)-pow(t->j15_e,2)) ;
                    t->jet_emfrac = t->j15_emfrac;
                    t->jet_hadfrac = t->j15_hadfrac;
                    t->jet_chadfrac = t->j15_chadfrac;
                    t->jet_nhadfrac = t->j15_nhadfrac;
                    t->jet_phofrac = t->j15_phofrac;
                    t->jet_mufrac = t->j15_mufrac;
                    t->jet_phofrac = t->j15_phofrac;
                    t->jet_JECUnc = t->j15_JECUnc;
                    t->jet_leadTrackPt = t->j15_leadTrackPt;
                    t->jet_softLeptPt = (t->j15_softLeptIdLooseMu==1 || t->j15_softLeptIdEle95==1) ? (t->j15_softLeptPt) : (-99);
                    t->jet_softLeptPtRel = (t->j15_softLeptIdLooseMu==1 || t->j15_softLeptIdEle95==1) ? (t->j15_softLeptPtRel) : (-99);
                    t->jet_softLeptDR = (t->j15_softLeptIdLooseMu==1 || t->j15_softLeptIdEle95==1) ? (t->j15_softLeptDR) : (-99);
                } else {
                    t->jet_secVtxPt = 0.;
                    t->jet_secVtx3dL = 0.;
                    t->jet_secVtx3deL = 0.;
                    t->jet_secVtxM = 0.;
            		t->jet_mt = 0.;
                    t->jet_emfrac = 0.;
                    t->jet_hadfrac = 0.;
                    t->jet_chadfrac = 0.;
                    t->jet_nhadfrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_mufrac = 0.;
                    t->jet_phofrac = 0.;
                    t->jet_JECUnc = 0.;
                    t->jet_leadTrackPt = 0.;
                    t->jet_softLeptPt = 0.;
                    t->jet_softLeptPtRel = 0.;
                    t->jet_softLeptDR = 0.;
                }
                t->jet_btagSF_M = t->j15_btagSF_M;
                t->jet_flavour = t->j15_flavour;
                t->jet_btagSFErrorUp_M = t->j15_btagSFErrorUp_M;
                t->jet_btagSFErrorDown_M = t->j15_btagSFErrorDown_M;
                t->jet_btagEff_M = t->j15_btagEff_M;
                t->jet_btagEffError_M = t->j15_btagEffError_M;
                t->jet_cutbased_wp_level = t->j15_cutbased_wp_level;
                t->jet_simple_wp_level = t->j15_simple_wp_level;
                t->jet_full_wp_level = t->j15_full_wp_level;
                t->jet_betaStarClassic = t->j15_betaStarClassic;
                t->jet_dR2Mean = t->j15_dR2Mean;
                t->jet_nNeutrals = t->j15_nNeutrals;
                t->jet_nCharged = t->j15_nCharged;
                t->jet_nConstituents = t->jet_nNeutrals + t->jet_nCharged;
                jet.SetPtEtaPhiE(t->j15_pt, t->j15_eta, t->j15_phi, t->j15_e);
                t->jet_dPhiMet = (fabs(jet.Phi()-met.Phi())>3.14159265 ) ? (2*3.14159265-fabs(jet.Phi()-met.Phi())) : (fabs(jet.Phi()-met.Phi()));//to be consistent with training definition
                t->jet_jecD_e = t->j15_jecD_e;
                t->jet_jecD_pt = t->j15_jecD_pt;
                t->jet_jecD_phi = t->j15_jecD_phi;
                t->jet_jecD_eta = t->j15_jecD_eta;
                t->jet_jecU_e = t->j15_jecU_e;
                t->jet_jecU_pt = t->j15_jecU_pt;
                t->jet_jecU_phi = t->j15_jecU_phi;
                t->jet_jecU_eta = t->j15_jecU_eta;
                t->jet_jerD_e = t->j15_jerD_e;
                t->jet_jerD_pt = t->j15_jerD_pt;
                t->jet_jerD_phi = t->j15_jerD_phi;
                t->jet_jerD_eta = t->j15_jerD_eta;
                t->jet_jerC_e = t->j15_jerC_e;
                t->jet_jerC_pt = t->j15_jerC_pt;
                t->jet_jerC_phi = t->j15_jerC_phi;
                t->jet_jerC_eta = t->j15_jerC_eta;
                t->jet_jerU_e = t->j15_jerU_e;
                t->jet_jerU_pt = t->j15_jerU_pt;
                t->jet_jerU_phi = t->j15_jerU_phi;
                t->jet_jerU_eta = t->j15_jerU_eta;
            } 
            else
            {
                t->jet_pt = 0.;
                return false;
            }
            return true;
}

float getPESUncertainty(bool isEB, float sceta, float r9, float pt)
{
    // numbers (in %) are taken from AN 2013/253 v7 (Hgg Moriond Legacy 2014)
//    if( isEB && (fabs(sceta) < 1.) && (r9 < .94) ) return 0.05 * 0.01;
//    else if( isEB && (fabs(sceta) < 1.) && (r9 > .94) ) return 0.05 * 0.01;
//    else if( isEB && (fabs(sceta) > 1.) && (r9 < .94) ) return 0.05 * 0.01;
//    else if( isEB && (fabs(sceta) > 1.) && (r9 > .94) ) return 0.10 * 0.01;
//    else if( !isEB && (fabs(sceta) < 2.) && (r9 < .94) ) return 0.10 * 0.01;
//    else if( !isEB && (fabs(sceta) < 2.) && (r9 > .94) ) return 0.10 * 0.01;
//    else if( !isEB && (fabs(sceta) > 2.) && (r9 < .94) ) return 0.10 * 0.01;
//    else if( !isEB && (fabs(sceta) > 2.) && (r9 > .94) ) return 0.05 * 0.01;
    // numbers (in %) are taken from CMS-PAS-HIG-13-001 (Hgg Moriond 2013: http://cds.cern.ch/record/1530524/files/HIG-13-001-pas.pdf)
    if( pt < 100. )
    {
        if( isEB && (fabs(sceta) < 1.) && (r9 < .94) ) return 0.20 * 0.01;
        else if( isEB && (fabs(sceta) < 1.) && (r9 > .94) ) return 0.20 * 0.01;
        else if( isEB && (fabs(sceta) > 1.) && (r9 < .94) ) return 0.51 * 0.01;
        else if( isEB && (fabs(sceta) > 1.) && (r9 > .94) ) return 0.71 * 0.01;
        else if( !isEB && (fabs(sceta) < 2.) && (r9 < .94) ) return 0.18 * 0.01;
        else if( !isEB && (fabs(sceta) < 2.) && (r9 > .94) ) return 0.88 * 0.01;
        else if( !isEB && (fabs(sceta) > 2.) && (r9 < .94) ) return 0.12 * 0.01;
        else if( !isEB && (fabs(sceta) > 2.) && (r9 > .94) ) return 0.12 * 0.01;
        else return 1.0;
    } else { // Assign 1% on PES if photon pt > 100 GeV
        return 1. * 0.01;
    }
}

float getPERUncertainty(bool isEB, float sceta, float r9, float sigmaEoE, TRandom3 *r)
{
    // numbers (in %) are taken from AN 2013/253 v7 (Hgg Moriond Legacy 2014)
//    if( isEB && (fabs(sceta) < 1.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.05 * 0.01)/sigmaEoE,2) - 1.));
//    else if( isEB && (fabs(sceta) < 1.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.05 * 0.01)/sigmaEoE,2) - 1.));
//    else if( isEB && (fabs(sceta) > 1.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.09 * 0.01)/sigmaEoE,2) - 1.));
//    else if( isEB && (fabs(sceta) > 1.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.10 * 0.01)/sigmaEoE,2) - 1.));
//    else if( !isEB && (fabs(sceta) < 2.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.09 * 0.01)/sigmaEoE,2) - 1.));
//    else if( !isEB && (fabs(sceta) < 2.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.07 * 0.01)/sigmaEoE,2) - 1.));
//    else if( !isEB && (fabs(sceta) > 2.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.06 * 0.01)/sigmaEoE,2) - 1.));
//    else if( !isEB && (fabs(sceta) > 2.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.03 * 0.01)/sigmaEoE,2) - 1.));
    // numbers (in %) are taken from CMS-PAS-HIG-13-001 (Hgg Moriond 2013: http://cds.cern.ch/record/1530524/files/HIG-13-001-pas.pdf)
    if( isEB && (fabs(sceta) < 1.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.25 * 0.01)/sigmaEoE,2) - 1.));
    else if( isEB && (fabs(sceta) < 1.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.23 * 0.01)/sigmaEoE,2) - 1.));
    else if( isEB && (fabs(sceta) > 1.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.60 * 0.01)/sigmaEoE,2) - 1.));
    else if( isEB && (fabs(sceta) > 1.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.72 * 0.01)/sigmaEoE,2) - 1.));
    else if( !isEB && (fabs(sceta) < 2.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.33 * 0.01)/sigmaEoE,2) - 1.));
    else if( !isEB && (fabs(sceta) < 2.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.93 * 0.01)/sigmaEoE,2) - 1.));
    else if( !isEB && (fabs(sceta) > 2.) && (r9 < .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.54 * 0.01)/sigmaEoE,2) - 1.));
    else if( !isEB && (fabs(sceta) > 2.) && (r9 > .94) ) return r->Gaus(0, sigmaEoE * sqrt( pow(1.+(0.36 * 0.01)/sigmaEoE,2) - 1.));
    else return 1.0;
}


float getCosThetaStar_CS(TLorentzVector h1, TLorentzVector h2, float ebeam = 4000)
{// cos theta star angle in the Collins Soper frame
    TLorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
    p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

    TLorentzVector hh;
    hh = h1 + h2;

    TVector3 boost = - hh.BoostVector();
    p1.Boost(boost);
    p2.Boost(boost);
    h1.Boost(boost);

    TVector3 CSaxis = p1.Vect().Unit() - p2.Vect().Unit();
    CSaxis.Unit();

    return cos(   CSaxis.Angle( h1.Vect().Unit() )    );
}

void getHandMadeCiCLevel(bool *ph1_cic4, bool *ph2_cic4, bool *ph1_cic0, bool *ph2_cic0, tree_variables * t,  bool noIsoA1 = false, bool noIsoA2 = false, bool noIsoB1 = false, bool noIsoB2 = false)
{
    *ph1_cic4 = false;
    *ph2_cic4 = false;
    *ph1_cic0 = false;
    *ph2_cic0 = false;
    bool is0NoCuts = true; // if false then cic0 is set to phoLOOSE working point
                           // photon preselection is applied anyway, so no photon ID looser than this should be applied
    // ph1
    if (t->ph1_isEB && t->ph1_r9_cic > .94)
    {
        *ph1_cic4 = 
               (noIsoA1 ? t->pho1_PFisoA < 8.9 : t->pho1_PFisoA < 6. )
            && (noIsoB1 ? t->pho1_PFisoB < 43. : t->pho1_PFisoB < 10.)
            && t->pho1_PFisoC < 3.8
            && t->ph1_sieie < 0.0108
            && t->ph1_hoe < 0.124
            && t->ph1_r9_cic > .94;
        *ph1_cic0 = is0NoCuts ? 1 : (
               t->pho1_PFisoA < 8.9
            && t->pho1_PFisoB < 43.
            && t->pho1_PFisoC < 6.2
            && t->ph1_sieie < 0.0117
            && t->ph1_hoe < 0.137
            && t->ph1_r9_cic > .94);
    }
    else if (t->ph1_isEB && t->ph1_r9_cic < .94)
    {
/*
        std::cout << "cat 1" << std::endl;
        std::cout << "t->ph1_isEB= " << t->ph1_isEB << std::endl;
        std::cout << "t->ph1_r9_cic= " << t->ph1_r9_cic << std::endl;
        std::cout << "t->pho1_PFisoA= " << t->pho1_PFisoA << std::endl;
        std::cout << "t->pho1_PFisoB= " << t->pho1_PFisoB << std::endl;
        std::cout << "t->pho1_PFisoC= " << t->pho1_PFisoC << std::endl;
        std::cout << "t->ph1_sieie= " << t->ph1_sieie << std::endl;
        std::cout << "t->ph1_hoe= " << t->ph1_hoe << std::endl;
        std::cout << "t->ph1_r9_cic= " << t->ph1_r9_cic << std::endl;
*/
        *ph1_cic4 = 
               (noIsoA1 ? t->pho1_PFisoA < 6.3  : t->pho1_PFisoA < 4.7 )
            && (noIsoB1 ? t->pho1_PFisoB < 19.4 : t->pho1_PFisoB < 6.5 )
            && t->pho1_PFisoC < 2.5
            && t->ph1_sieie < 0.0102
            && t->ph1_hoe < 0.092
            && t->ph1_r9_cic > .28;
        *ph1_cic0 = is0NoCuts ? 1 : (
               t->pho1_PFisoA < 6.3
            && t->pho1_PFisoB < 19.4
            && t->pho1_PFisoC < 4.3
            && t->ph1_sieie < 0.0105
            && t->ph1_hoe < 0.14
//            && t->ph1_r9_cic > .25;
            && t->ph1_r9_cic > .28);
    }
    else if (!t->ph1_isEB && t->ph1_r9_cic > .94)
    {
        *ph1_cic4 = 
               (noIsoA1 ? t->pho1_PFisoA < 9.8 : t->pho1_PFisoA < 5.6 )
            && (noIsoB1 ? t->pho1_PFisoB < 24. : t->pho1_PFisoB < 5.6 )
            && t->pho1_PFisoC < 3.1
            && t->ph1_sieie < 0.028
            && t->ph1_hoe < 0.142
            && t->ph1_r9_cic > .94;
        *ph1_cic0 = is0NoCuts ? 1 : (
               t->pho1_PFisoA < 9.8
            && t->pho1_PFisoB < 24.
            && t->pho1_PFisoC < 5.
            && t->ph1_sieie < 0.031
            && t->ph1_hoe < 0.145
//            && t->ph1_r9_cic > .93;
            && t->ph1_r9_cic > .94);
    }
    else if (!t->ph1_isEB && t->ph1_r9_cic < .94)
    {
        *ph1_cic4 = 
               (noIsoA1 ? t->pho1_PFisoA < 6.8 : t->pho1_PFisoA < 3.6 )
            && (noIsoB1 ? t->pho1_PFisoB < 7.9 : t->pho1_PFisoB < 4.4 )
            && t->pho1_PFisoC < 2.2
            && t->ph1_sieie < 0.028
            && t->ph1_hoe < 0.063
            && t->ph1_r9_cic > .24;
        *ph1_cic0 = is0NoCuts ? 1 : (
               t->pho1_PFisoA < 6.8
            && t->pho1_PFisoB < 7.9
            && t->pho1_PFisoC < 4.3
            && t->ph1_sieie < 0.031
            && t->ph1_hoe < 0.143
            && t->ph1_r9_cic > .24);
    }
    // ph2
    if (t->ph2_isEB && t->ph2_r9_cic > .94)
    {
        *ph2_cic4 = 
               (noIsoA2 ? t->pho2_PFisoA < 8.9 : t->pho2_PFisoA < 6. )
            && (noIsoB2 ? t->pho2_PFisoB < 43. : t->pho2_PFisoB < 10.)
            && t->pho2_PFisoC < 3.8
            && t->ph2_sieie < 0.0108
            && t->ph2_hoe < 0.124
            && t->ph2_r9_cic > .94;
        *ph2_cic0 = is0NoCuts ? 1 : (
               t->pho2_PFisoA < 8.9
            && t->pho2_PFisoB < 43.
            && t->pho2_PFisoC < 6.2
            && t->ph2_sieie < 0.0117
            && t->ph2_hoe < 0.137
            && t->ph2_r9_cic > .94);
    }
    else if (t->ph2_isEB && t->ph2_r9_cic < .94)
    {
        *ph2_cic4 = 
               (noIsoA2 ? t->pho2_PFisoA < 6.3  : t->pho2_PFisoA < 4.7 )
            && (noIsoB2 ? t->pho2_PFisoB < 19.4 : t->pho2_PFisoB < 6.5 )
            && t->pho2_PFisoC < 2.5
            && t->ph2_sieie < 0.0102
            && t->ph2_hoe < 0.092
            && t->ph2_r9_cic > .28;
        *ph2_cic0 = is0NoCuts ? 1 : (
               t->pho2_PFisoA < 6.3
            && t->pho2_PFisoB < 19.4
            && t->pho2_PFisoC < 4.3
            && t->ph2_sieie < 0.0105
            && t->ph2_hoe < 0.14
//            && t->ph2_r9_cic > .25;
            && t->ph2_r9_cic > .28);
    }
    else if (!t->ph2_isEB && t->ph2_r9_cic > .94)
    {
        *ph2_cic4 = 
               (noIsoA2 ? t->pho2_PFisoA < 9.8 : t->pho2_PFisoA < 5.6 )
            && (noIsoB2 ? t->pho2_PFisoB < 24. : t->pho2_PFisoB < 5.6 )
            && t->pho2_PFisoC < 3.1
            && t->ph2_sieie < 0.028
            && t->ph2_hoe < 0.142
            && t->ph2_r9_cic > .94;
        *ph2_cic0 = is0NoCuts ? 1 : (
               t->pho2_PFisoA < 9.8
            && t->pho2_PFisoB < 24.
            && t->pho2_PFisoC < 5.
            && t->ph2_sieie < 0.031
            && t->ph2_hoe < 0.145
//            && t->ph2_r9_cic > .93;
            && t->ph2_r9_cic > .94);
    }
    else if (!t->ph2_isEB && t->ph2_r9_cic < .94)
    {
        *ph2_cic4 = 
               (noIsoA2 ? t->pho2_PFisoA < 6.8 : t->pho2_PFisoA < 3.6 )
            && (noIsoB2 ? t->pho2_PFisoB < 7.9 : t->pho2_PFisoB < 4.4 )
            && t->pho2_PFisoC < 2.2
            && t->ph2_sieie < 0.028
            && t->ph2_hoe < 0.063
            && t->ph2_r9_cic > .24;
        *ph2_cic0 = is0NoCuts ? 1 : (
               t->pho2_PFisoA < 6.8
            && t->pho2_PFisoB < 7.9
            && t->pho2_PFisoC < 4.3
            && t->ph2_sieie < 0.031
            && t->ph2_hoe < 0.143
            && t->ph2_r9_cic > .24);
    }
    return;
}
