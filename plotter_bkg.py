#!/usr/bin/env python
#// Small dumb code to play with trees
#// O.Bondu, F. Bojarski (May 2014)
# Various python imports
from os import path
from math import log10, pow
import collections
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1D, TLatex, TLegend, TLorentzVector
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
ROOT.TGaxis.SetMaxDigits(3);

c1 = TCanvas()
CIC1 = "1"
CIC2 = "1"
#afs_plottree = "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/v31_weighted_pf_study_ciclevel_handmade_44_loose_isoAB_pho2_fitToMgg_withKinFit"
afs_plottree = "2014-06-06_selection_noRegression_noMassCut_v12_study_ciclevel_0_0"
eos_tree = "root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_11_tree_08"
Debug = "false"

intL = 19706.
samples = []
# samples.append([ name, typ, dirpath, subdir, file, tree, sample_cut, color, style, label , sigma , N])
#samples.append(["Radion_m270", -270, afs_plottree, "", "Radion_m270_8TeV_noRegression_noMassCut_v*.root", "Radion_m270_8TeV", "evweight_w_btagSF", ROOT.kRed+1, 0, "m_{X} = 270 GeV" , 13.55e-3, 19999])
#samples.append(["Radion_m300", -300, afs_plottree, "", "Radion_m300_8TeV_noRegression_noMassCut_v*.root", "Radion_m300_8TeV", "evweight_w_btagSF", ROOT.kRed+2, 0, "m_{X} = 300 GeV" , 13.55e-3, 19999])
#samples.append(["Radion_m350", -350, afs_plottree, "", "Radion_m350_8TeV_noRegression_noMassCut_v*.root", "Radion_m350_8TeV", "evweight_w_btagSF", ROOT.kRed+3, 0, "m_{X} = 350 GeV" , 13.55e-3, 19999])
#samples.append(["Data", 0, afs_plottree, "", "Data_noRegression_noMassCut_v*.root", "Data", "", ROOT.kBlack, 0, "Data", 1., 1001822])
#samples.append(["qcd_40_8TeV_ff", 30, afs_plottree, "", "qcd_40_8TeV_ff_noRegression_noMassCut_v*.root", "qcd_40_8TeV_ff", "(evweight_w_btagSF) * (evweight_w_btagSF< 10)", ROOT.kCyan+2, 1001, "QCD jets" , 1., 14404429])
#samples.append(["qcd_30_8TeV_ff", 30, afs_plottree, "", "qcd_30_8TeV_ff_noRegression_noMassCut_v*.root", "qcd_30_8TeV_ff", "(evweight_w_btagSF) * (evweight_w_btagSF< 10)", ROOT.kCyan+2, 1001, "QCD jets" , 1., 14404429])
samples.append(["qcd_30_8TeV_pf", 20, afs_plottree, "", "qcd_30_8TeV_pf_noRegression_noMassCut_v*.root", "qcd_30_8TeV_pf", "(1) * (evweight_w_btagSF) * (evweight_w_btagSF< 25)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
samples.append(["qcd_40_8TeV_pf", 20, afs_plottree, "", "qcd_40_8TeV_pf_noRegression_noMassCut_v*.root", "qcd_40_8TeV_pf", "(1) * (evweight_w_btagSF) * (evweight_w_btagSF< 25)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
samples.append(["gjet_20_8TeV_pf", 20, afs_plottree, "", "gjet_20_8TeV_pf_noRegression_noMassCut_v*.root", "gjet_20_8TeV_pf", "(1) *(evweight_w_btagSF)", ROOT.kBlue-4, 2501, "QCD #gamma + jets" , 1., 14404429])
samples.append(["gjet_40_8TeV_pf", 20, afs_plottree, "", "gjet_40_8TeV_pf_noRegression_noMassCut_v*.root", "gjet_40_8TeV_pf", "(1) * (evweight_w_btagSF)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
#samples.append(["diphojet_sherpa_8TeV", 10, afs_plottree, "", "diphojet_sherpa_8TeV_noRegression_noMassCut_v*.root", "diphojet_sherpa_8TeV", "(evweight_w_btagSF)", ROOT.kGreen+2, 1001, "QCD #gamma#gamma + jets" , 1., 14404429])
#samples.append(["DYJetsToLL", 5, afs_plottree, "", "DYJetsToLL_noRegression_noMassCut_v*.root", "DYJetsToLL", "(evweight_w_btagSF)", ROOT.kMagenta+2, 1001, "DY" , 1., 14404429])


#samples.append(["Radion_m270", -270, afs_plottree, "","Radion_m270*.root", "TCVARS", "(reweight)", ROOT.kRed, 0, "M_{X} = 270 GeV", 10e-8, 1])
#
#samples.append(["qcd_40_8TeV_ff", 30, afs_plottree, "", "qcd_40_8TeV_ff_m270.root", "TCVARS", "(reweight)", ROOT.kCyan+2, 1001, "QCD jets" , 1., 14404429])
##samples.append(["qcd_30_8TeV_ff", 30, afs_plottree, "", "qcd_30_8TeV_ff_m270.root", "TCVARS", "(reweight)", ROOT.kCyan+2, 1001, "QCD jets" , 1., 14404429])
#samples.append(["qcd_30_8TeV_pf", 20, afs_plottree, "", "qcd_30_8TeV_pf_m270.root", "TCVARS", "(reweight)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
##samples.append(["qcd_40_8TeV_pf", 20, afs_plottree, "", "qcd_40_8TeV_pf_m270.root", "TCVARS", "(reweight)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
#samples.append(["gjet_20_8TeV_pf", 20, afs_plottree, "", "gjet_20_8TeV_pf_m270.root", "TCVARS", "(reweight)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
#samples.append(["gjet_40_8TeV_pf", 20, afs_plottree, "", "gjet_40_8TeV_pf_m270.root", "TCVARS", "(reweight)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
#samples.append(["diphojet_sherpa_8TeV", 10, afs_plottree, "", "diphojet_sherpa_8TeV_m270.root", "TCVARS", "(reweight)", ROOT.kGreen+2, 1001, "QCD #gamma#gamma + jets" , 1., 14404429])
#samples.append(["DYJetsToLL", 5, afs_plottree, "", "DYJetsToLL_m270.root", "TCVARS", "(reweight)", ROOT.kMagenta+2, 1001, "DY" , 1., 14404429])
#


#####plots.append([ name2, variable, plot_cut, norm, binning, title, additional_info, cutline, cutline2 ])
plots = []
#plots.append(["pho1_pt", "pho1_pt", "", intL, "(100, 0, 500)", "p_{T}^{#gamma1} (GeV)", "", 33.3, ""])
#plots.append(["pho2_pt", "pho2_pt", "", intL, "(100, 0, 500)", "p_{T}^{#gamma2} (GeV)", "", 25., ""])
#plots.append(["jet1_pt", "jet1_pt", "", intL, "(100, 0, 500)", "p_{T}^{jet1} (GeV)", "", 25., ""])
#plots.append(["jet2_pt", "jet2_pt", "", intL, "(100, 0, 500)", "p_{T}^{jet1} (GeV)", "", 25., ""])
#plots.append(["pho1_eta", "pho1_eta", "", intL, "(100, -5, 5)", "eta^{#gamma1}", "", 2.5, -2.5])
#plots.append(["pho2_eta", "pho2_eta", "", intL, "(100, -5, 5)", "eta^{#gamma2}", "", 2.5, -2.5])
#plots.append(["jet1_eta", "jet1_eta", "", intL, "(100, -5, 5)", "eta^{jet1}", "", "", ""])
#plots.append(["jet2_eta", "jet2_eta", "", intL, "(100, -5, 5)", "eta^{jet1}", "", "", ""])
#plots.append(["pho1_phi", "pho1_phi", "", intL, "(100, -5, 5)", "phi^{#gamma1}", "", "", ""])
#plots.append(["pho2_phi", "pho2_phi", "", intL, "(100, -5, 5)", "phi^{#gamma2}", "", "", ""])
#plots.append(["jet1_phi", "jet1_phi", "", intL, "(100, -5, 5)", "phi^{jet1}", "", "", ""])
#plots.append(["jet2_phi", "jet2_phi", "", intL, "(100, -5, 5)", "phi^{jet2}", "", "", ""])
#plots.append(["jj_phi", "jj_phi", "", intL, "(100, -5, 5)", "phi^{jj}", "", "", ""])
#plots.append(["jj_eta", "jj_eta", "", intL, "(100, -5, 5)", "eta^{jj}", "", "", ""])
#plots.append(["gg_phi", "gg_phi", "", intL, "(100, -5, 5)", "phi^{#gamma#gamma}", "", "", ""])
#plots.append(["gg_eta", "gg_eta", "", intL, "(100, -5, 5)", "eta^{#gamma#gamma}", "", "", ""])
#plots.append(["jj_pt", "jj_pt", "", intL, "(500, 0, 500)", "p_{T}^{jj} (GeV)", "", "", ""])
#plots.append(["gg_pt", "gg_pt", "", intL, "(500, 0, 500)", "p_{T}^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_mass", "pho1_mass", "", intL, "(100, -1, 1)", "mass^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_mass", "pho2_mass", "", intL, "(100, -1, 1)", "mass^{#gamma2} (GeV)", "", "", ""])
#plots.append(["jet1_mass", "jet1_mass", "", intL, "(100, 0, 100)", "mass^{jet1} (GeV)", "", "", ""])
#plots.append(["jet2_mass", "jet2_mass", "", intL, "(100, 0, 100)", "mass^{jet2} (GeV)", "", "", ""])
#plots.append(["jj_mass_normData_cut_" + CIC1 + "_" + CIC2, "jj_mass", "", "data", "(50, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
#plots.append(["gg_mass_normData_cut_" + CIC1 + "_" + CIC2, "gg_mass", "", "data", "(40, 100, 180)", "mass^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["jj_mass_norm1", "jj_mass", "", 1, "(50, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
#plots.append(["gg_mass_norm1", "gg_mass", "", 1, "(40, 100, 180)", "mass^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["jj_mass_cut_" + CIC1 + "_" + CIC2, "jj_mass", "", intL, "(50, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
#plots.append(["gg_mass_cut_" + CIC1 + "_" + CIC2, "gg_mass", "", intL, "(40, 100, 180)", "mass^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_e", "pho1_e", "", intL, "(500, 0, 500)", "e^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_e", "pho2_e", "", intL, "(500, 0, 500)", "e^{#gamma2} (GeV)", "", "", ""])
#plots.append(["jet1_e", "jet1_e", "", intL, "(500, 0, 500)", "e^{jet1} (GeV)", "", "", ""])
#plots.append(["jet2_e", "jet2_e", "", intL, "(500, 0, 500)", "e^{jet2} (GeV)", "", "", ""])
#plots.append(["jj_e", "jj_e", "", intL, "(500, 0, 500)", "e^{jj} (GeV)", "", "", ""])
#plots.append(["gg_e", "gg_e", "", intL, "(500, 0, 500)", "e^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_r9", "pho1_r9", "", intL, "(100, 0, 5)", "r9^{#gamma1}", "", "", ""])
#plots.append(["pho2_r9", "pho2_r9", "", intL, "(100, 0, 5)", "r9^{#gamma2}", "", "", ""])
#plots.append(["ggjj_pt", "ggjj_pt", "", intL, "(500, 0, 500)", "p_{T}^{#gamma#gammajj} (GeV)", "", "", ""])
#plots.append(["ggjj_e", "ggjj_e", "", intL, "(1000, 0, 1000)", "e^{#gamma#gammajj} (GeV)", "", "", ""])
#plots.append(["ggjj_mass_normData_cut_" + CIC1 + "_" + CIC2, "ggjj_mass / 1000.", "", "data", "(50, 0, 1.)", "mass^{#gamma#gammajj} (TeV)", "", "", ""])
#plots.append(["ggjj_mass_norm1", "ggjj_mass / 1000.", "", 1, "(50, 0, 1.)", "mass^{#gamma#gammajj} (TeV)", "", "", ""])
#plots.append(["ggjj_mass_cut_" + CIC1 + "_" + CIC2, "ggjj_mass / 1000.", "", intL, "(50, 0, 1.)", "mass^{#gamma#gammajj} (TeV)", "", "", ""])
#plots.append(["ggjj_phi", "ggjj_phi", "", intL, "(100, -5, 5)", "phi^{#gamma#gammajj}", "", "", ""])
#plots.append(["ggjj_eta", "ggjj_eta", "", intL, "(100, -10, 10)", "eta^{#gamma#gammajj}", "", "", ""])
#plots.append(["jet1_betaStarClassic", "jet1_betaStarClassic", "", intL, "(20, 0, 1)", "betaStarClassic^{jet1}", "", "", ""])
#plots.append(["jet2_betaStarClassic", "jet2_betaStarClassic", "", intL, "(20, 0, 1)", "betaStarClassic^{jet2}", "", "", ""])
#plots.append(["jet1_dR2Mean", "jet1_dR2Mean", "", intL, "(80,0, 0.2)", "dR2Mean^{jet1}", "", "", ""])
#plots.append(["jet2_dR2Mean", "jet2_dR2Mean", "", intL, "(80,0, 0.2)", "dR2Mean^{jet2}", "", "", ""])
#plots.append(["pho1_hoe", "pho1_hoe", "", intL, "(100,0, .05)", "hoe^{#gamma1}", "", "", ""])
#plots.append(["pho2_hoe", "pho2_hoe", "", intL, "(100,0, .05)", "hoe^{#gamma2}", "", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_sieie_endcap_allcat", "pho1_sieie", "(pho1_isEB == 0)", intL, "(35, 0.015, 0.05)", "sieie^{#gamma1}", "endcap_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_sieie_endcap_allcat", "pho2_sieie", "(pho2_isEB == 0)", intL, "(35, 0.015, 0.05)", "sieie^{#gamma2}", "endcap_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoA_endcap_allcat", "pho1_PFisoA", "(pho1_isEB == 0)", intL, "(20,0, 10)", "PFisoA^{#gamma1}", "endcap_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoB_endcap_allcat", "pho1_PFisoB", "(pho1_isEB == 0)", intL, "(30,-5, 10)", "PFisoB^{#gamma1}", "endcap_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoA_endcap_allcat", "pho2_PFisoA", "(pho2_isEB == 0)", intL, "(20,0, 10)", "PFisoA^{#gamma2}", "endcap_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoB_endcap_allcat", "pho2_PFisoB", "(pho2_isEB == 0)", intL, "(30,-5, 10)", "PFisoB^{#gamma2}", "endcap_allcat", "", ""])
#
#plots.append([CIC1 + "_" + CIC2 + "pho1_sieie_barrel_allcat", "pho1_sieie", "(pho1_isEB == 1)", intL, "(30, 0, 0.015)", "sieie^{#gamma1}", "barrel_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_sieie_barrel_allcat", "pho2_sieie", "(pho2_isEB == 1)", intL, "(30, 0, 0.015)", "sieie^{#gamma2}", "barrel_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoA_barrel_allcat", "pho1_PFisoA", "(pho1_isEB == 1)", intL, "(20,0, 10)", "PFisoA^{#gamma1}", "barrel_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoB_barrel_allcat", "pho1_PFisoB", "(pho1_isEB == 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma1}", "barrel_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoA_barrel_allcat", "pho2_PFisoA", "(pho2_isEB == 1)", intL, "(20,0, 10)", "PFisoA^{#gamma2}", "barrel_allcat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoB_barrel_allcat", "pho2_PFisoB", "(pho2_isEB == 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma2}", "barrel_allcat", "", ""])
#
#
#plots.append([CIC1 + "_" + CIC2 + "pho1_sieie_endcap_cat0and1", "pho1_sieie", "(pho1_isEB == 0) * (category >= 1)", intL, "(35, 0.015, 0.05)", "sieie^{#gamma1}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_sieie_endcap_cat0and1", "pho2_sieie", "(pho2_isEB == 0) * (category >= 1)", intL, "(35, 0.015, 0.05)", "sieie^{#gamma2}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoA_endcap_cat0and1", "pho1_PFisoA", "(pho1_isEB == 0) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma1}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoB_endcap_cat0and1", "pho1_PFisoB", "(pho1_isEB == 0) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma1}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoA_endcap_cat0and1", "pho2_PFisoA", "(pho2_isEB == 0) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma2}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoB_endcap_cat0and1", "pho2_PFisoB", "(pho2_isEB == 0) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma2}", "endcap_cat0and1", "", ""])
#
#plots.append([CIC1 + "_" + CIC2 + "pho1_sieie_barrel_cat0and1", "pho1_sieie", "(pho1_isEB == 1) * (category >= 1)", intL, "(30, 0, 0.015)", "sieie^{#gamma1}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_sieie_barrel_cat0and1", "pho2_sieie", "(pho2_isEB == 1) * (category >= 1)", intL, "(30, 0, 0.015)", "sieie^{#gamma2}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoA_barrel_cat0and1", "pho1_PFisoA", "(pho1_isEB == 1) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma1}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoB_barrel_cat0and1", "pho1_PFisoB", "(pho1_isEB == 1) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma1}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoA_barrel_cat0and1", "pho2_PFisoA", "(pho2_isEB == 1) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma2}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoB_barrel_cat0and1", "pho2_PFisoB", "(pho2_isEB == 1) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma2}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "_K_factor_isoA_allcat", "pho1_PFisoA", "(((pho2_isEB == 0) && (pho2_sieie > 0.028)) || ((pho2_isEB == 1) && (pho2_sieie > 0.0108)))", intL, "(40, -5, 15)", "isoA^{#gamma1}", ", all cat", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "_K_factor_isoB_allcat", "pho1_PFisoB", "(((pho2_isEB == 0) * (pho2_sieie > 0.028)) + ((pho2_isEB == 1) * (pho2_sieie > 0.0108)))", intL, "(50, -10, 15)", "isoB^{#gamma1}", ", all cat", "", ""])

#plots.append([CIC1 + "_" + CIC2 + "_K_factor_isoA", "pho1_PFisoA", "(((pho2_isEB == 0) && (pho2_sieie > 0.028)) || ((pho2_isEB == 1) && (pho2_sieie > 0.0108))) && (category >= 1)", intL, "(40, -5, 15)", "isoA^{#gamma1}", "", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "_K_factor_isoB", "pho1_PFisoB", "(((pho2_isEB == 0) * (pho2_sieie > 0.028)) + ((pho2_isEB == 1) * (pho2_sieie > 0.0108)))* (category >= 1)", intL, "(50, -10, 15)", "isoB^{#gamma1}", "", "", ""])

plots.append(["K_factor_pho1_isoB_both_ph_cat01_wcut", "pho1_PFisoB", "(((ph1_ciclevel >= 4) && (((pho2_isEB == 0) && (pho2_sieie > 0.028)) || ((pho2_isEB == 1) && (pho2_sieie > 0.0108)))) || ((ph2_ciclevel >= 4) && (((pho1_isEB == 0) && (pho1_sieie > 0.028)) || ((pho1_isEB == 1) && (pho1_sieie > 0.0108))))) && (category >= 1)", intL, "(50, -10, 15)", "isoB^{#gamma1}", "", "", ""])
#plots.append(["0_K_factor_pho1_isoB_both_ph_allcat_nowcut", "pho1_PFisoB", "(((ph1_ciclevel >= 4) && (ph2_ciclevel == 0)) || ((ph2_ciclevel >= 4) && (ph1_ciclevel == 0)))", intL, "(50, -10, 15)", "isoB^{#gamma1}", "", "", ""])
#plots.append(["0_K_factor_pho1_isoB_both_ph_nowcut", "pho1_PFisoB", "((((ph1_ciclevel >= 4) && (ph2_ciclevel == 0)) || ((ph2_ciclevel >= 4) && (ph1_ciclevel == 0))) && (category >= 1))", intL, "(50, -10, 15)", "isoB^{#gamma1}", "", "", ""])
#
#plots.append(["0_K_factor_pho1_isoB_both_ph_allcat_weightcut", "pho1_PFisoB", "((((ph1_ciclevel >= 4) && (ph2_ciclevel == 0)) || ((ph2_ciclevel >= 4) && (ph1_ciclevel == 0))) && (evweight_w_btagSF< 25))", intL, "(50, -10, 15)", "isoB^{#gamma1}", "", "", ""])
#plots.append(["0_K_factor_pho1_isoB_both_ph_weightcut", "pho1_PFisoB", "((((ph1_ciclevel >= 4) && (ph2_ciclevel == 0)) || ((ph2_ciclevel >= 4) && (ph1_ciclevel == 0))) && (category >= 1) && (evweight_w_btagSF< 25)) ", intL, "(50, -10, 15)", "isoB^{#gamma1}", "", "", ""])

#plots.append(["mass_gg_CIC_44", "gg_mass", "(category ==1)", intL, "(20,100, 180)", "mass^{#gamma #gamma}", "", "", ""])

#plots.append(["mass_jj", "jj_mass", "(category ==1)", intL, "(14, 85, 155)", "mass^{jj}", "", "", ""])

#plots.append(["mjj_m270_wkinfit", "mjj_wkinfit", "", 1, "(16,80,160)", "mass^{ jj} with KinFit", "", "", ""])
#plots.append(["mtot_cat1_m400", "mtot", "(cut_based_ct ==1)", intL, "(16,365,445)", "mass^{#gamma #gamma jj} (GeV)", "", "", ""])
#plots.append(["mtot_cat0_m400", "mtot", "(cut_based_ct ==0)", intL, "(16,365,445)", "mass^{#gamma #gamma jj} (GeV)", "", "", ""])

#plots.append(["mass_ggjj", "ggjj_mass", "", 1, "(40,100, 1100)", "mass^{#gamma #gamma jj}", "", "", ""])
#plots.append(["IsoA_pho2_barrel_highr9", "pho2_PFisoA", "(category >=1) * (pho2_isEB == 1) * (pho2_r9 > .94)", intL, "(50, -5, 45)", "IsoA^{#gamma 2}", "", 6, ""])
#plots.append(["IsoA_pho2_barrel_lowr9", "pho2_PFisoA", "(category >=1) * (pho2_isEB == 1) * (pho2_r9 < .94)", intL, "(50, -5, 45)", "IsoA^{#gamma 2}", "", 4.7, ""])
#plots.append(["IsoB_pho2_barrel_highr9", "pho2_PFisoB", "(category >=1) * (pho2_isEB == 1) * (pho2_r9 > .94)", intL, "(50, -5, 45)", "IsoB^{#gamma 2}", "", 10, ""])
#plots.append(["IsoB_pho2_barrel_lowr9", "pho2_PFisoB", "(category >=1) * (pho2_isEB == 1) * (pho2_r9 < .94)", intL, "(50, -5, 45)", "IsoB^{#gamma 2}", "", 6.5, ""])
#plots.append(["IsoA_pho2_endcap_highr9", "pho2_PFisoA", "(category >=1) * (pho2_isEB == 0) * (pho2_r9 > .94)", intL, "(50, -5, 45)", "IsoA^{#gamma 2}", "", 5.6, ""])
#plots.append(["IsoA_pho2_endcap_lowr9", "pho2_PFisoA", "(category >=1) * (pho2_isEB == 0) * (pho2_r9 < .94)", intL, "(50, -5, 45)", "IsoA^{#gamma 2}", "", 3.6, ""])
#plots.append(["IsoB_pho2_endcap_highr9", "pho2_PFisoB", "(category >=1) * (pho2_isEB == 0) * (pho2_r9 > .94)", intL, "(50, -5, 45)", "IsoB^{#gamma 2}", "", 5.6, ""])
#plots.append(["IsoB_pho2_endcap_lowr9", "pho2_PFisoB", "(category >=1) * (pho2_isEB == 0) * (pho2_r9 < .94)", intL, "(50, -5, 45)", "IsoB^{#gamma 2}", "", 4.4, ""])
#plots.append(["pho2_sieie_" + CIC1 +"_" + CIC2 + "_reweight" , "pho2_sieie", "(ph1_ciclevel >= " + CIC1 + ") * (ph2_ciclevel >= " + CIC2 +") * (category >=1)", intL, "(60, 0.004, 0.016)", "sieie^{#gamma 2}", "", 0.0108, ""])


#plots.append([CIC1 + "_" + CIC2 + "pho1_sieie_endcap_cat0and1_CUT_SIEIE", "pho1_sieie", "(pho1_isEB == 0) * (pho2_sieie > 0.028) * (category >= 1)", intL, "(35, 0.015, 0.05)", "sieie^{#gamma1}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_sieie_endcap_cat0and1_CUT_SIEIE", "pho2_sieie", "(pho2_isEB == 0) * (pho2_sieie > 0.028) * (category >= 1)", intL, "(35, 0.015, 0.05)", "sieie^{#gamma2}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoA_endcap_cat0and1_CUT_SIEIE", "pho1_PFisoA", "(pho1_isEB == 0) * (pho2_sieie > 0.028) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma1}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoB_endcap_cat0and1_CUT_SIEIE", "pho1_PFisoB", "(pho1_isEB == 0) * (pho2_sieie > 0.028) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma1}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoA_endcap_cat0and1_CUT_SIEIE", "pho2_PFisoA", "(pho2_isEB == 0) * (pho2_sieie > 0.028) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma2}", "endcap_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoB_endcap_cat0and1_CUT_SIEIE", "pho2_PFisoB", "(pho2_isEB == 0) * (pho2_sieie > 0.028) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma2}", "endcap_cat0and1", "", ""])
#
#plots.append([CIC1 + "_" + CIC2 + "pho1_sieie_barrel_cat0and1_CUT_SIEIE", "pho1_sieie", "(pho1_isEB == 1) * (pho2_sieie > 0.0108) * (category >= 1)", intL, "(30, 0, 0.015)", "sieie^{#gamma1}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_sieie_barrel_cat0and1_CUT_SIEIE", "pho2_sieie", "(pho2_isEB == 1) * (pho2_sieie > 0.0108) * (category >= 1)", intL, "(30, 0, 0.015)", "sieie^{#gamma2}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoA_barrel_cat0and1_CUT_SIEIE", "pho1_PFisoA", "(pho1_isEB == 1) * (pho2_sieie > 0.0108) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma1}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho1_PFisoB_barrel_cat0and1_CUT_SIEIE", "pho1_PFisoB", "(pho1_isEB == 1) * (pho2_sieie > 0.0108) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma1}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoA_barrel_cat0and1_CUT_SIEIE", "pho2_PFisoA", "(pho2_isEB == 1) * (pho2_sieie > 0.0108) * (category >= 1)", intL, "(20,0, 10)", "PFisoA^{#gamma2}", "barrel_cat0and1", "", ""])
#plots.append([CIC1 + "_" + CIC2 + "pho2_PFisoB_barrel_cat0and1_CUT_SIEIE", "pho2_PFisoB", "(pho2_isEB == 1) * (pho2_sieie > 0.0108) * (category >= 1)", intL, "(30,-5, 10)", "PFisoB^{#gamma2}", "barrel_cat0and1", "", ""])
#
#plots.append(["pho1_isconv", "pho1_isconv", "", intL, "(100,0, 1)", "isconv^{#gamma1}", "", "", ""])
#plots.append(["pho2_isconv", "pho2_isconv", "", intL, "(100,0, 1)", "isconv^{#gamma2}", "", "", ""])
#plots.append(["jet1_csvBtag", "jet1_csvBtag", "", intL, "(100, 0, 1)", "csvBtag^{jet1}", "", "", ""])
#plots.append(["jet2_csvBtag", "jet2_csvBtag", "", intL, "(100, 0, 1)", "csvBtag^{jet2}", "", "", ""])


#plots.append(["DeltaEta_gg_jj_eff", "DeltaEta_gg_jj", "", 1, "(40, -5, 5)", "#Delta#eta (#gamma#gamma, jj)", "", "", ""])
#plots.append(["DeltaPhi_gg_jj", "DeltaPhi_gg_jj", "", 1, "(40, -4, 4)", "DeltaPhi_gg_jj", "", "", ""])
#plots.append(["DeltaR_gg_jj", "DeltaR_gg_jj", "", 1, "(50, 0, 8)", "DeltaR_gg_jj", "", "", ""])


for name2, variable, plot_cut, norm, binning, title, additional_info, cutline, cutline2 in plots:
    c1 = TCanvas()
    legend = TLegend(0.45, 0.82, 0.90, 0.93, "")
    legend.SetTextSize(0.025)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    xnbin, xlow, xhigh = map(float, binning.strip().strip("()").split(","))
    ymax = -1
    ymin = 10000000
    firsthistname = ""
    if plot_cut == "": plot_cut = "1"
    hist_signal = {}
    hist_data = {}
    hist_bkg = {}
    label_signal = {}
    label_data = {}
    label_bkg = {}

    for ifile, [ name, typ, dirpath, subdir, file, tree, sample_cut, color, style, label , sigma , N] in enumerate(samples):
#        if Debug:
        print ifile, file, color, style, label
        chain = TChain(tree)
        chain.Add( path.join(dirpath, subdir, file) )
        total_cut = plot_cut
        if sample_cut == "": sample_cut = "1"
        if typ < 0:
            total_cut = "(" + plot_cut + ") * (" + str(sigma) + " * " + str(intL) + ")/" + str(N)
        elif typ == 0:
           total_cut = "(" + plot_cut + ")* (" + sample_cut + ")"
        elif typ > 0:
           total_cut = "(" + plot_cut + ") * (" + sample_cut + ")"
#        if Debug:
#            print "total_cut = ", total_cut
        option = ""
        if ifile != 0:
            option = "same"
        if typ == 0:
            option += "e1"
        chain.Draw(variable + ">>h_tmp" + binning, total_cut, option)
        # Cosmetics
        h = ROOT.gDirectory.Get("h_tmp")
#        if Debug:
#            print type(h)
        h.SetName(name + "_" + name2 + "_" + str(ifile))
        if ifile == 0:
            firsthistname = name + "_" + name2 + "_" + str(ifile)
        h.SetLineWidth(2)
        h.SetLineColor(color)
        h.SetFillColor(color)
        h.SetFillStyle(style)
        h.GetXaxis().SetTitle( title )
        unit = ""
        if title.find("(") != -1:
            unit = title[title.find("(")+1:title.find(")")]
        if norm == 1. or norm == 1:
            h.GetYaxis().SetTitle( "Norm. to unity / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        elif norm == "data" or norm == "Data":
            h.GetYaxis().SetTitle( "Norm. to data / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        else:
            h.GetYaxis().SetTitle( "# events / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        # store histo for redraw in the correct order later
        if typ > 0:
            label_bkg[typ] = label
            if typ not in hist_bkg:
                hist_bkg[typ] = h
            else:
                hist_bkg[typ].Add(h)
        elif typ == 0:
            label_data[typ] = label
            hist_data[typ] = h
        elif typ < 0:
            label_signal[typ] = label
            hist_signal[typ] = h
        del chain, h

#        print hist_bkg
#        print hist_signal
#        print hist_data

    # Sum the backgrounds
    for key in collections.OrderedDict(sorted(hist_bkg.items())):
        for jkey in collections.OrderedDict(sorted(hist_bkg.items())):
            if jkey <= key: continue
#            print key, jkey
            hist_bkg[key].Add(hist_bkg[jkey])

# give the number of events of different category
    for ikey, key in enumerate(collections.OrderedDict(sorted(hist_bkg.items()))):
        if ikey == 0:
            Int_all_bkg =  hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
        elif ikey == 1:
            Int_pp_pf_ff = hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
        elif ikey == 2:
            Int_pf_ff = hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
        elif ikey == 3:
            Int_ff = hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
        
    for ikey, key in enumerate(collections.OrderedDict(sorted(hist_data.items()))):
        if ikey == 0:
            Int_data = hist_data[key].Integral(0, hist_data[key].GetNbinsX() +1)
        else:
            continue

    for ikey, key in enumerate(collections.OrderedDict(sorted(hist_signal.items()))):
        if ikey == 0:
            Int_signal = hist_signal[key].Integral(0, hist_signal[key].GetNbinsX() +1)
        else:
            continue
#    Int_DY = Int_all_bkg - Int_pp_pf_ff
#    Int_pp = Int_pp_pf_ff - Int_pf_ff
#    Int_pf = Int_pf_ff - Int_ff







    # Adjust norm if case happens
    if norm == 1. or norm == 1:
        bkg_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_bkg.items()))):
            if ikey == 0:
                bkg_integral = hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
            hist_bkg[key].Scale( 1. / bkg_integral )
        data_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_data.items()))):
            if ikey == 0:
                data_integral = hist_data[key].Integral(0, hist_data[key].GetNbinsX() +1)
            hist_data[key].Scale( 1. / data_integral )
        signal_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_signal.items()))):
            if ikey == 0:
                signal_integral = hist_signal[key].Integral(0, hist_signal[key].GetNbinsX() +1)
            hist_signal[key].Scale( 1. / signal_integral )
    elif norm == "data" or norm == "Data":
        data_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_data.items()))):
            if ikey == 0:
                data_integral = hist_data[key].Integral(0, hist_data[key].GetNbinsX() +1)
            else:
                continue
        bkg_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_bkg.items()))):
            if ikey == 0:
                bkg_integral = hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
            hist_bkg[key].Scale( data_integral / bkg_integral )
        signal_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_signal.items()))):
            if ikey == 0:
                signal_integral = hist_signal[key].Integral(0, hist_signal[key].GetNbinsX() +1)
            hist_signal[key].Scale( data_integral / signal_integral )
    # redraw in order : background, data, signal, axis
    if len(hist_bkg) + len(hist_data) + len(hist_signal) > 1:
        legend.SetNColumns(2)
    if len(hist_bkg) + len(hist_data) + len(hist_signal) > 6:
        legend.SetNColumns(3)
    for ikey, key in enumerate(collections.OrderedDict(sorted(hist_bkg.items()))):
        if ikey == 0:
            hist_bkg[key].Draw("")
            firsthistname = hist_bkg[key].GetName()
        else:
            hist_bkg[key].Draw("same")
        legend.AddEntry(hist_bkg[key].GetName(), label_bkg[key], "lf")
        ymax = max(ymax, hist_bkg[key].GetMaximum())
        ymin = min(ymin, hist_bkg[key].GetMinimum(0.0))
    for key in hist_data:
        hist_data[key].Draw("e1same")
        legend.AddEntry(hist_data[key].GetName(), label_data[key], "lpe")
        ymax = max(ymax, hist_data[key].GetMaximum())
        ymin = min(ymin, hist_data[key].GetMinimum(0.0))
    for key in collections.OrderedDict(sorted(hist_signal.items())):
        hist_signal[key].Draw("same")
        legend.AddEntry(hist_signal[key].GetName(), label_signal[key], "lf")
        ymax = max(ymax, hist_signal[key].GetMaximum())
        ymin = min(ymin, hist_signal[key].GetMinimum(0.0))
    
#    if Debug:
#        print "ymin = ", ymin, "ymax = ", ymax
    ymin_lin = ymin / 10.
    yrange_lin = ymax - ymin_lin
    ymax_lin = .25 * yrange_lin + ymax
    yrange_log = (log10(ymax) - log10(ymin)) / .77
    ymax_log = pow(10., .25*yrange_log + log10(ymax))
    ymin_log = pow(10., log10(ymin) - .03*yrange_log)

    latexLabel = TLatex()
    latexLabel.SetTextSize(.02)
    latexLabel.SetNDC()
    latexLabel.DrawLatex(.25, .96, "CMS Internal     L = 19.7 fb^{-1}     #sqrt{s} = 8 TeV   " + additional_info)
#    latexLabel.DrawLatex(.15, .9, "#Data = " + str(round(Int_data,1)))
#    latexLabel.DrawLatex(.15, .9, "Eff bkg = " + str(round(Int_all_bkg,1)))
#    latexLabel.DrawLatex(.15, .85, "Eff signal = " + str(round(Int_signal,1)))
#    latexLabel.DrawLatex(.15, .85, "#pp = " + str(round(Int_pp,1)) + ", #pf = " + str(round(Int_pf,1)))
#    latexLabel.DrawLatex(.15, .8, "#ff = " +str (round(Int_ff,1)) + ", #DY = " + str(round(Int_DY,1))) 
#    latexLabel.DrawLatex(.15, .85, additional_info)
#    latexLabel.DrawLatex(.15, .9, "CIC1 = " + CIC1 + " CIC2 =" + CIC2)
#    print "#Data = ", str(round(Int_data,1)), CIC1, CIC2
    ROOT.gPad.RedrawAxis()
    legend.Draw()
    c1.Update()

    line = TLine()
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line2 = TLine()
    line2.SetLineStyle(2)
    line2.SetLineWidth(2)

    h = ROOT.gDirectory.Get(firsthistname)
    h.SetMaximum(ymax_lin)
    h.SetMinimum(ymin_lin)
    if cutline != "":
        line.SetX1(cutline); line.SetY1(ymin_lin); line.SetX2(cutline); line.SetY2(ymax)
        line.Draw("same")
    if cutline2 != "":
        line2.SetX1(cutline2); line2.SetY1(ymin_lin); line2.SetX2(cutline2); line2.SetY2(ymax)
        line2.Draw("same")
    c1.Update()
    c1.Print("pdf/" + name2 + ".pdf")
    c1.Print("gif/" + name2 + ".gif")
    c1.Print("root/" + name2 + ".root")


    c1.SetLogy(1)
    h.SetMaximum(ymax_log)
    h.SetMinimum(ymin_log)
    h.GetYaxis().SetRangeUser(ymin_log, ymax_log)
    if cutline != "":
        line.SetX1(cutline); line.SetY1(ymin_log); line.SetX2(cutline); line.SetY2(ymax)
        line.Draw("same")
    if cutline2 != "":
        line2.SetX1(cutline2); line2.SetY1(ymin_log); line2.SetX2(cutline2); line2.SetY2(ymax)
        line2.Draw("same")
    c1.Update()
    c1.Print("pdf/" + name2 + "_log.pdf")
    c1.Print("gif/" + name2 + "_log.gif")
    c1.Print("root/" + name2 + "_log.root")
    c1.SetLogy(0)

    del c1

