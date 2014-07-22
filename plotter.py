#!/usr/bin/env python
#// Small dumb code to play with trees
#// O.Bondu, F. Bojarski (May 2014)
# Various python imports
from os import path
from math import log10, pow
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1D, TLatex, TLegend, TLorentzVector, TStyle, gStyle
import math as m
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
ROOT.TGaxis.SetMaxDigits(3);

c1 = TCanvas()

afs_plottree = "/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v10/"
eos_tree = "root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_11_tree_08"

intL = 19706.
samples = [] 
# samples.append([ name, dirpath, subdir, file, tree, color, style, label , sigma , N])

samples.append(["Radion_m300", ".", "2014-07-16_selection_noRegression_noMassCut_v12_study_cic_WP_handmade_4_4_loose_isoAB_pho2", "Radion_m300_*.root", "Radion_m300_8TeV", ROOT.kBlack, 0, "Radion M = 300 GeV" , 13.55e-3, 19999])
#samples.append(["Radion_m300", ".", "2014-06-06_selection_noRegression_noMassCut_v12_study_ciclevel_4_4", "Radion_m300_8TeV_noRegression_noMassCut_v12_study_ciclevel_4_4.root", "Radion_m300_8TeV", ROOT.kBlack, 0, "mass = 300 GeV, WP 4" , 13.55e-3, 19999])
#samples.append(["Radion_m300", ".", "2014-06-06_selection_noRegression_noMassCut_v12_study_ciclevel_0_0", "Radion_m300_8TeV_noRegression_noMassCut_v12_study_ciclevel_0_0.root", "Radion_m300_8TeV", ROOT.kBlack, 0, "Radion (M=300 GeV), CIC >= 0" , 13.55e-3, 19999])
#samples.append(["Radion_m300", ".", "2014-06-19_selection_noRegression_noMassCut_v12_study_iso_variables_WP44_handmade", "Radion_m300_8TeV_noRegression_noMassCut_v12_study_iso_variables_WP44_handmade.root", "Radion_m300_8TeV", ROOT.kRed, 0, "hand made WP 4" , 13.55e-3, 19999])
#samples.append(["Radion_m300", ".", "2014-06-18_selection_noRegression_noMassCut_v12_study_iso_variables_WP44", "Radion_m300_8TeV_noRegression_noMassCut_v12_study_iso_variables_WP44.root", "Radion_m300_8TeV", ROOT.kRed, 0, "hand made WP 4" , 13.55e-3, 19999])
#samples.append(["Data", ".", "2014-06-06_selection_noRegression_noMassCut_v12_study_ciclevel_0_1", "Data*.root", "Data", ROOT.kRed, 0, "" , 13.55e-3, 19999])

#CIC1 = "0"
#CIC2 = "1"
#MASSE = "270"
#samples.append(["Data", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "v31_study_ciclevel_" +  CIC1 +"_" + CIC2 +"_fitToMgg_noKinFit/", "Data_m" + MASSE + ".root", "TCVARS", ROOT.kBlack, 20, "Data", 13.55e-3, 19999])
#
#samples.append(["Sherpa", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "v31_study_ciclevel_" +  CIC1 +"_" + CIC2 +"_fitToMgg_noKinFit/", "diphojet_sherpa_8TeV_m" + MASSE + ".root", "TCVARS", ROOT.kRed, 1, "Sherpa", 13.55e-3, 19999])
#
#samples.append(["pouet", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "v31_study_ciclevel_" +  CIC1 +"_" + CIC2 +"_fitToMgg_noKinFit/", "DYJetsToLL_m" + MASSE +".root", "TCVARS", ROOT.kGreen, 1, "DYJet", 13.55e-3, 19999])
#
#samples.append(["pouet", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "v31_study_ciclevel_" +  CIC1 +"_" + CIC2 +"_fitToMgg_noKinFit/", "pf_" + MASSE + ".root", "TCVARS", ROOT.kBlue, 1, "pf", 13.55e-3, 19999])
#samples.append(["pouet", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "v31_study_ciclevel_" +  CIC1 +"_" + CIC2 +"_fitToMgg_noKinFit/", "ff_" + MASSE + ".root", "TCVARS", ROOT.kYellow, 1, "ff", 13.55e-3, 19999])
#
#samples.append(["pouet", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "2014-06-03_selection_noRegression_noMassCut_v10bis_effStudies_RMS_BetaStar", "diphojet_sherpa_8TeV_noRegression_noMassCut_v10bis_effStudies_RMS_BetaStar.root", "diphojet_sherpa_8TeV", ROOT.kBlack, 20, "Sherpa ", 13.55e-3, 19999])

#samples.append(["pouet", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "2014-06-03_selection_noRegression_noMassCut_v10bis_effStudies_RMS_BetaStar_nocut", "diphojet_sherpa_8TeV_noRegression_noMassCut_v10bis_effStudies_RMS_BetaStar_nocut.root", "diphojet_sherpa_8TeV", ROOT.kRed, 20, "Sherpa_nocut ", 13.55e-3, 19999])

#samples.append(["pouet", "", "2014-06-25_selection_noRegression_noMassCut_v12_study_ciclevel_3_0", "Radion_m270_8TeV_noRegression_noMassCut_v12_study_ciclevel_3_0.root", "Radion_m270_8TeV", ROOT.kBlack, 0, "3_0", 13.55e-3, 19999])
#samples.append(["pouet", "", "2014-06-06_selection_noRegression_noMassCut_v12_study_ciclevel_3_1", "Radion_m270_8TeV_noRegression_noMassCut_v12_study_ciclevel_3_1.root", "Radion_m270_8TeV", ROOT.kBlue, 0, "3_1", 13.55e-3, 19999])
#samples.append(["pouet", "", "2014-06-06_selection_noRegression_noMassCut_v12_study_ciclevel_3_2", "Radion_m270_8TeV_noRegression_noMassCut_v12_study_ciclevel_3_2.root", "Radion_m270_8TeV", ROOT.kGreen, 0, "3_2", 13.55e-3, 19999])
#samples.append(["pouet", "", "2014-06-06_selection_noRegression_noMassCut_v12_study_ciclevel_3_3", "Radion_m270_8TeV_noRegression_noMassCut_v12_study_ciclevel_3_3.root", "Radion_m270_8TeV", ROOT.kYellow, 0, "3_3", 13.55e-3, 19999])
#samples.append(["Radion_m300", "", "2014-07-10_selection_noRegression_noMassCut_v12_study_cic_no_cic_cut", "Radion_m300*.root", "Radion_m300_8TeV", ROOT.kRed, 0, "n-1 plot", 13.55e-3, 19999])
#samples.append(["Radion_m300", "", "2014-07-11_selection_noRegression_noMassCut_v12_study_cic_WP_handmade_4_4", "Radion_m300*.root", "Radion_m300_8TeV", ROOT.kBlack, 0, "WP 4_4", 13.55e-3, 19999])

#####plots.append([ name2, variable, cut,critere_eff, norm, Scale Factor, binning, title, additional_info, cutline, cutline2 ])
plots = []

#plots.append(["pho1_cic", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "", "", ""])
#plots.append(["pho1_cic_without_1_isoA", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (-1 < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (-1 < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (-1 < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (-1 < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoA^{#gamma 1}", "", ""])
#plots.append(["pho1_cic_without_1_isoB", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (-1 < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (-1 < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (-1 < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (-1 < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoB^{#gamma 1}", "", ""])
#plots.append(["pho1_cic_without_1_isoC", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (-1 < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (-1 < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (-1 < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (-1 < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoC^{#gamma 1}", "", ""])
#plots.append(["pho1_cic_without_1_sieie", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (-1 < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (-1 < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (-1 < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (-1 < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without sieie^{#gamma 1}", "", ""])
#plots.append(["pho1_cic_without_1_hoe", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (-1 < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (-1 < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (-1 < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (-1 < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without hoe^{#gamma 1}", "", ""])
#plots.append(["pho1_cic_without_1_r9", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (1000 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (1000 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (1000 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (1000 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without r9^{#gamma 1}", "", ""])
#
#plots.append(["pho1_cic_without_2_isoA", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (-1 < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (-1 < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (-1 < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (-1 < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoA^{#gamma 2}", "", ""])
#plots.append(["pho1_cic_without_2_isoB", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (-1 < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (-1 < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (-1 < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (-1 < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoB^{#gamma 2}", "", ""])
#plots.append(["pho1_cic_without_2_isoC", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (-1 < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (-1 < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (-1< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (-1 < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoC^{#gamma 2}", "", ""])
#plots.append(["pho1_cic_without_2_sieie", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (-1 < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (-1 < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (-1 < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (-1 < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without sieie^{#gamma 2}", "", ""])
#plots.append(["pho1_cic_without_2_hoe", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (-1 < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (-1 < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (-1 < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (-1 < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without hoe^{#gamma 2}", "", ""])
#plots.append(["pho1_cic_without_2_r9", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (1000 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (1000 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (1000 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (1000 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without r9^{#gamma 2}", "", ""])
#
#plots.append(["pho1_cic_without_1_isoA", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (-1 < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (-1 < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (-1 < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (-1 < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (pho2_PFisoA < 6) && (pho2_PFisoB < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (pho2_PFisoA < 4.7) && (pho2_PFisoB < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (pho2_PFisoA < 5.6) && (pho2_PFisoB < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (pho2_PFisoA < 3.6) && (pho2_PFisoB < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoA^{#gamma 1}", "", ""])
#
#plots.append(["pho1_cic_without_2_isoA_B", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (pho1_PFisoA < 6) && (pho1_PFisoB < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (pho1_PFisoA < 4.7) && (pho1_PFisoB < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (pho1_PFisoA < 5.6) && (pho1_PFisoB < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (pho1_PFisoA < 3.6) && (pho1_PFisoB < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (-1 < 6) && (-1 < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (-1 < 4.7) && (-1 < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (-1 < 5.6) && (-1 < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (-1 < 3.6) && (-1 < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoA^{#gamma 2}, isoB^{#gamma 2}", "", ""])
#
#plots.append(["pho1_cic_without_1_2_isoA_B", "ph1_ciclevel", "((((pho1_isEB == 1) && (pho1_r9 > 0.94) && (-1 < 6) && (-1 < 10) && (pho1_PFisoC < 3.8) && (pho1_sieie < 0.0108)  && (pho1_hoe < 0.124) && (pho1_r9 > 0.94)  ) || ((pho1_isEB == 1) && (pho1_r9 < 0.94) && (-1 < 4.7) && (-1 < 6.5) && (pho1_PFisoC < 2.5) && (pho1_sieie < 0.0102)  && (pho1_hoe < 0.092) && (pho1_r9 > 0.298)) || ((pho1_isEB == 0) && (pho1_r9 > 0.94) && (-1 < 5.6) && (-1 < 5.6) && (pho1_PFisoC < 3.1) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.142) && (pho1_r9 > 0.94)) || ((pho1_isEB == 0) && (pho1_r9 < 0.94) && (-1 < 3.6) && (-1 < 4.4) && (pho1_PFisoC < 2.2) && (pho1_sieie < 0.028)  && (pho1_hoe < 0.063) && (pho1_r9 > 0.24))) && (((pho2_isEB == 1) && (pho2_r9 > 0.94) && (-1 < 6) && (-1 < 10) && (pho2_PFisoC < 3.8) && (pho2_sieie < 0.0108)  && (pho2_hoe < 0.124) && (pho2_r9 > 0.94)) || ((pho2_isEB == 1) && (pho2_r9< 0.94) && (-1 < 4.7) && (-1 < 6.5) && (pho2_PFisoC < 2.5) && (pho2_sieie < 0.0102)  && (pho2_hoe < 0.092) && (pho2_r9 > 0.298)) || ((pho2_isEB == 0) && (pho2_r9 > 0.94) && (-1 < 5.6) && (-1 < 5.6) && (pho2_PFisoC< 3.1) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.142) && (pho2_r9 > 0.94)) || ((pho2_isEB == 0) && (pho2_r9 < 0.94) && (-1 < 3.6) && (-1 < 4.4) && (pho2_PFisoC < 2.2) && (pho2_sieie < 0.028)  && (pho2_hoe < 0.063) && (pho2_r9 > 0.24)))  && (category >= 1))", "", "", "evweight", "(14, 0, 14)", "ciclevel^{#gamma1}", "without isoA^{#gamma 1}, isoB^{#gamma 1}, isoA^{#gamma 2}, isoB^{#gamma 2}", "", ""])
#
#
#plots.append(["pho1_pt", "pho1_pt", "", "ph1_ciclevel >= 4", 2., "evweight", "(50, 40, 120)", "p_{T}^{#gamma1} (GeV)", "", 33.3, ""])
#plots.append(["pho2_pt", "pho2_pt", "", "(ph2_ciclevel >= 4)*(ph1_ciclevel >= 4)", 2., "evweight", "(25, 20, 120)", "p_{T}^{#gamma2} (GeV)", "", 25., ""])
#plots.append(["jet1_pt", "jet1_pt", "", "", 1., "", "(100, 0, 500)", "p_{T}^{jet1} (GeV)", "", 25., ""])
#plots.append(["jet2_pt", "jet2_pt", "", "", 1., "", "(100, 0, 500)", "p_{T}^{jet1} (GeV)", "", 25., ""])
#plots.append(["pho1_eta", "pho1_eta", "", "", 1., "", "(100, -5, 5)", "eta^{#gamma1}", "", 2.5, -2.5])
#plots.append(["pho2_eta_bidule", "pho2_eta", "", "", 2., "evweight", "(70, 1.3, 2)", "eta^{#gamma2}", "", 2.5, -2.5])
#plots.append(["jet1_eta", "jet1_eta", "", "", 1., "", "(100, -5, 5)", "eta^{jet1}", "", "", ""])
#plots.append(["jet2_eta", "jet2_eta", "", "", 1., "", "(100, -5, 5)", "eta^{jet1}", "", "", ""])
#plots.append(["pho1_phi", "pho1_phi", "", "", 1., "", "(100, -5, 5)", "phi^{#gamma1}", "", "", ""])
#plots.append(["pho2_phi", "pho2_phi", "", "ph2_ciclevel >= 4", 2., "evweight", "(25, -5, 5)", "phi^{#gamma2}", "", "", ""])
#plots.append(["jet1_phi", "jet1_phi", "", "", 1., "", "(100, -5, 5)", "phi^{jet1}", "", "", ""])
#plots.append(["jet2_phi", "jet2_phi", "", "", 1., "", "(100, -5, 5)", "phi^{jet2}", "", "", ""])
#plots.append(["jj_phi", "jj_phi", "", "", 1., "", "(100, -5, 5)", "phi^{jj}", "", "", ""])
#plots.append(["jj_eta", "jj_eta", "", "", 1., "", "(100, -5, 5)", "eta^{jj}", "", "", ""])
#plots.append(["gg_phi", "gg_phi", "", "", 1., "", "(100, -5, 5)", "phi^{#gamma#gamma}", "", "", ""])
#plots.append(["gg_eta", "gg_eta", "", "", 1., "", "(100, -5, 5)", "eta^{#gamma#gamma}", "", "", ""])
#plots.append(["jj_pt", "jj_pt", "", "", 1., "", "(500, 0, 500)", "p_{T}^{jj} (GeV)", "", "", ""])
#plots.append(["gg_pt", "gg_pt", "", "", 1., "", "(500, 0, 500)", "p_{T}^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_mass", "pho1_mass", "", "", 1., "", "(100, -1, 1)", "mass^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_mass", "pho2_mass", "", "ph2_ciclevel >= 4", 2., "evweight", "(16, -.4, .4)", "mass^{#gamma2} (GeV)", "", "", ""])
#plots.append(["jet1_mass", "jet1_mass", "", "", 1., "", "(100, 0, 100)", "mass^{jet1} (GeV)", "", "", ""])
#plots.append(["jet2_mass", "jet2_mass", "", "", 1., "", "(100, 0, 100)", "mass^{jet2} (GeV)", "", "", ""])
#plots.append(["jj_mass", "jj_mass", "", "", 1., "", "(500, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
#plots.append(["gg_mass", "gg_mass", "", "ph2_ciclevel >= 4", 2., "evweight", "(60, 100, 160)", "mass^{#gamma#gamma} (GeV)", "", 100., 180.])
#plots.append(["pho1_e", "pho1_e", "", "", 1., "", "(500, 0, 500)", "e^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_e", "pho2_e", "", "ph2_ciclevel >= 4", 2., "evweight", "(50, 0, 200)", "e^{#gamma2} (GeV)", "", "", ""])
#plots.append(["jet1_e", "jet1_e", "", "", 1., "", "(500, 0, 500)", "e^{jet1} (GeV)", "", "", ""])
#plots.append(["jet2_e", "jet2_e", "", "", 1., "", "(500, 0, 500)", "e^{jet2} (GeV)", "", "", ""])
#plots.append(["jj_e", "jj_e", "", "", 1., "", "(500, 0, 500)", "e^{jj} (GeV)", "", "", ""])
#plots.append(["gg_e", "gg_e", "", "", 1., "", "(500, 0, 500)", "e^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_r9", "pho1_r9", "","", 1., "", "(20, 0, 1)", "r9^{#gamma1}", "", "", ""])
#plots.append(["pho2_r9", "pho2_r9", "", "ph2_ciclevel >= 4",2., "evweight", "(20, 0, 1)", "r9^{#gamma2}", "", "", ""])
#plots.append(["ggjj_pt", "ggjj_pt", "","", 1., "", "(500, 0, 500)", "p_{T}^{#gamma#gammajj} (GeV)", "", "", ""])
#plots.append(["ggjj_e", "ggjj_e", "", "", 1., "", "(1000, 0, 1000)", "e^{#gamma#gammajj} (GeV)", "", "", ""])
#plots.append(["ggjj_mass", "ggjj_mass", "","", 1., "", "(200, 0, 00)", "mass^{#gamma#gammajj} (GeV)", "", "", ""])
#plots.append(["ggjj_phi", "ggjj_phi", "","", 1., "", "(100, -5, 5)", "phi^{#gamma#gammajj}", "", "", ""])
#plots.append(["ggjj_eta", "ggjj_eta", "","", 1., "", "(100, -10, 10)", "eta^{#gamma#gammajj}", "", "", ""])
#plots.append(["pho1_sieie", "pho1_sieie", "","", 1., "", "(100, -0.5, 0.5)", "sieie^{#gamma1}", "", "", ""])
#plots.append(["pho2_sieie", "pho2_sieie", "","ph2_ciclevel >= 4", 2., "", "(20, 0.005, 0.015)", "sieie^{#gamma2}", "", "", ""])
#plots.append(["jet1_betaStarClassic", "jet1_betaStarClassic", "","", 1., "", "(20, 0, 1)", "betaStarClassic^{jet1}", "", "", ""])
#plots.append(["jet2_betaStarClassic", "jet2_betaStarClassic", "","", 1., "", "(20, 0, 1)", "betaStarClassic^{jet2}", "", "", ""])
#plots.append(["jet1_dR2Mean", "jet1_dR2Mean", "njets_kRadionID_and_CSVM == 1","", 1., "", "(80,0, 0.2)", "dR2Mean^{jet1}", "", "", ""])
#plots.append(["jet2_dR2Mean", "jet2_dR2Mean", "", "",1., "", "(80,0, 0.2)", "dR2Mean^{jet2}", "", "", ""])
#plots.append(["pho1_hoe", "pho1_hoe", "","", 1., "", "(100,0, .05)", "hoe^{#gamma1}", "", "", ""])
#plots.append(["pho2_hoe", "pho2_hoe", "","ph2_ciclevel >= 4", 2., "", "(100,0, .05)", "hoe^{#gamma2}", "", "", ""])

#plots.append(["pho2_PFisoA_cut_CIC_1", "pho1_PFisoA", "","ph2_ciclevel >= 1", 2., "", "(20,0, 20)", "PFisoA^{#gamma1}", "", "", ""])
#plots.append(["pho2_PFisoB_cut_CIC_1", "pho1_PFisoB", "","ph2_ciclevel >= 1", 2., "", "(20,-4, 16)", "PFisoB^{#gamma1}", "", "", ""])
#
#plots.append(["pho1_PFisoA_barrel_lowr9", "pho1_PFisoA", "(pho1_isEB == 1) * (pho1_r9 < 0.94) ","", 2., "", "(60,0, 6)", "PFisoA^{#gamma1}", "barrel_lowr9", 4.7, ""])
#plots.append(["pho1_PFisoA_barrel_highr9", "pho1_PFisoA", "(pho1_isEB == 1) * (pho1_r9 > 0.94) ","", 2., "", "(80,0, 8)", "PFisoA^{#gamma1}", "barrel_highr9", 6, ""])
#plots.append(["pho1_PFisoA_endcap_lowr9", "pho1_PFisoA", "(pho1_isEB == 0) * (pho1_r9 < 0.94) ","", 2., "", "(50,0, 5)", "PFisoA^{#gamma1}", "endcap_lowr9", 3.6, ""])
#plots.append(["pho1_PFisoA_endcap_highr9", "pho1_PFisoA", "(pho1_isEB == 0) * (pho1_r9 > 0.94) ","", 2., "", "(70,0, 7)", "PFisoA^{#gamma1}", "endcap_highr9", 5.6, ""])
#
#plots.append(["pho1_PFisoB_barrel_lowr9", "pho1_PFisoB", "(pho1_isEB == 1) * (pho1_r9 < 0.94) ","", 2., "", "(120,-4, 8)", "PFisoB^{#gamma1}", "barrel_lowr9", 6.5, ""])
#plots.append(["pho1_PFisoB_barrel_highr9", "pho1_PFisoB", "(pho1_isEB == 1) * (pho1_r9 > 0.94) ","", 2., "", "(160,-4, 12)", "PFisoB^{#gamma1}", "barrel_highr9", 10, ""])
#plots.append(["pho1_PFisoB_endcap_lowr9", "pho1_PFisoB", "(pho1_isEB == 0) * (pho1_r9 < 0.94) ","", 2., "", "(100,-4, 6)", "PFisoB^{#gamma1}", "endcap_lowr9", 4.4, ""])
#plots.append(["pho1_PFisoB_endcap_highr9", "pho1_PFisoB", "(pho1_isEB == 0) * (pho1_r9 > 0.94) ","", 2., "", "(110,-4, 7)", "PFisoB^{#gamma1}", "endcap_highr9", 5.6, ""])
#
#plots.append(["pho1_PFisoC_barrel_lowr9", "pho1_PFisoC", "(pho1_isEB == 1) * (pho1_r9 < 0.94) ","", 2., "", "(50,-1, 4)", "PFisoC^{#gamma1}", "barrel_lowr9", 2.5, ""])
#plots.append(["pho1_PFisoC_barrel_highr9", "pho1_PFisoC", "(pho1_isEB == 1) * (pho1_r9 > 0.94) ","", 2., "", "(60,-1, 5)", "PFisoC^{#gamma1}", "barrel_highr9", 3.8, ""])
#plots.append(["pho1_PFisoC_endcap_lowr9", "pho1_PFisoC", "(pho1_isEB == 0) * (pho1_r9 < 0.94) ","", 2., "", "(40,-1, 3)", "PFisoC^{#gamma1}", "endcap_lowr9", 2.2, ""])
#plots.append(["pho1_PFisoC_endcap_highr9", "pho1_PFisoC", "(pho1_isEB == 0) * (pho1_r9 > 0.94) ","", 2., "", "(50,-1, 4)", "PFisoC^{#gamma1}", "endcap_highr9", 3.1, ""])
#
#plots.append(["pho1_sieie_barrel_lowr9", "pho1_sieie", "(pho1_isEB == 1) * (pho1_r9 < 0.94) ","", 2., "", "(120,0, 0.0120)", "sieie^{#gamma1}", "barrel_lowr9", 0.0102, ""])
#plots.append(["pho1_sieie_barrel_highr9", "pho1_sieie", "(pho1_isEB == 1) * (pho1_r9 > 0.94) ","", 2., "", "(120,0, 0.0120)", "sieie^{#gamma1}", "barrel_highr9", 0.0108, ""])
#plots.append(["pho1_sieie_endcap_lowr9", "pho1_sieie", "(pho1_isEB == 0) * (pho1_r9 < 0.94) ","", 2., "", "(100,0.02, 0.03)", "sieie^{#gamma1}", "endcap_lowr9", 0.028, ""])
#plots.append(["pho1_sieie_endcap_highr9", "pho1_sieie", "(pho1_isEB == 0) * (pho1_r9 > 0.94) ","", 2., "", "(100,0.02, 0.03)", "sieie^{#gamma1}", "endcap_highr9", 0.028, ""])
#
#plots.append(["pho1_hoe_barrel_lowr9", "pho1_hoe", "(pho1_isEB == 1) * (pho1_r9 < 0.94) ","", 2., "", "(100,0, 0.1)", "hoe^{#gamma1}", "barrel_lowr9", 0.092, ""])
#plots.append(["pho1_hoe_barrel_highr9", "pho1_hoe", "(pho1_isEB == 1) * (pho1_r9 > 0.94) ","", 2., "", "(100, 0.04, 0.14)", "hoe^{#gamma1}", "barrel_highr9", 0.124, ""])
#plots.append(["pho1_hoe_endcap_lowr9", "pho1_hoe", "(pho1_isEB == 0) * (pho1_r9 < 0.94) ","", 2., "", "(80,0, 0.08)", "hoe^{#gamma1}", "endcap_lowr9", 0.063, ""])
#plots.append(["pho1_hoe_endcap_highr9", "pho1_hoe", "(pho1_isEB == 0) * (pho1_r9 > 0.94) ","", 2., "", "(150, 0, 0.150)", "hoe^{#gamma1}", "endcap_highr9", 0.142, ""])
#
#plots.append(["pho1_r9_barrel_lowr9", "pho1_r9", "(pho1_isEB == 1) * (pho1_r9 < 0.94) ","", 2., "", "(90,0.1, 1)", "r9^{#gamma1}", "barrel_lowr9", 0.298, ""])
#plots.append(["pho1_r9_barrel_highr9", "pho1_r9", "(pho1_isEB == 1) * (pho1_r9 > 0.94) ","", 2., "", "(100,0.8, 1)", "r9^{#gamma1}", "barrel_highr9", 0.94, ""])
#plots.append(["pho1_r9_endcap_lowr9", "pho1_r9", "(pho1_isEB == 0) * (pho1_r9 < 0.94) ","", 2., "", "(90,0.1, 1)", "r9^{#gamma1}", "endcap_lowr9", 0.24, ""])
#plots.append(["pho1_r9_endcap_highr9", "pho1_r9", "(pho1_isEB == 0) * (pho1_r9 > 0.94) ","", 2., "", "(100, 0.8, 1)", "r9^{#gamma1}", "endcap_highr9", 0.28, ""])
#
#
#plots.append(["pho2_PFisoA_barrel_lowr9", "pho2_PFisoA", "(pho2_isEB == 1) * (pho2_r9 < 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(80,0, 40)", "PFisoA^{#gamma1}", "barrel_lowr9", 4.7, ""])
plots.append(["pho2_PFisoA_barrel_highr9", "pho2_PFisoA", "(pho2_isEB == 1) * (pho2_r9 > 0.94)","", 2., "", "(15,0, 15)", "PFisoA^{#gamma2}", "barrel_highr9", 6, ""])
#plots.append(["pho2_PFisoA_endcap_lowr9", "pho2_PFisoA", "(pho2_isEB == 0) * (pho2_r9 < 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(80,0, 40)", "PFisoA^{#gamma2}", "endcap_lowr9", 3.6, ""])
#plots.append(["pho2_PFisoA_endcap_highr9", "pho2_PFisoA", "(pho2_isEB == 0) * (pho2_r9 > 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(80,0, 40)", "PFisoA^{#gamma2}", "endcap_highr9", 5.6, ""])

#plots.append(["pho2_PFisoB_barrel_lowr9", "pho2_PFisoB", "(pho2_isEB == 1) * (pho2_r9 < 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(80,0, 40)", "PFisoB^{#gamma2}", "barrel_lowr9", 6.5, ""])
plots.append(["pho2_PFisoB_barrel_highr9", "pho2_PFisoB", "(pho2_isEB == 1) * (pho2_r9 > 0.94)","", 2., "", "(20,0, 20)", "PFisoB^{#gamma2}", "barrel_highr9", 10, ""])
#plots.append(["pho2_PFisoB_endcap_lowr9", "pho2_PFisoB", "(pho2_isEB == 0) * (pho2_r9 < 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(80,0, 40)", "PFisoB^{#gamma2}", "endcap_lowr9", 4.4, ""])
#plots.append(["pho2_PFisoB_endcap_highr9", "pho2_PFisoB", "(pho2_isEB == 0) * (pho2_r9 > 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(80,0, 40)", "PFisoB^{#gamma2}", "endcap_highr9", 5.6, ""])

#plots.append(["pho2_PFisoC_barrel_lowr9", "pho2_PFisoC", "(pho2_isEB == 1) * (pho2_r9 < 0.94) ","", 2., "", "(50,-1, 4)", "PFisoC^{#gamma2}", "barrel_lowr9", 2.5, ""])
#plots.append(["pho2_PFisoC_barrel_highr9", "pho2_PFisoC", "(pho2_isEB == 1) * (pho2_r9 > 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(30,-1, 5)", "PFisoC^{#gamma2}", "barrel_highr9", 3.8, ""])
#plots.append(["pho2_PFisoC_endcap_lowr9", "pho2_PFisoC", "(pho2_isEB == 0) * (pho2_r9 < 0.94) ","", 2., "", "(40,-1, 3)", "PFisoC^{#gamma2}", "endcap_lowr9", 2.2, ""])
#plots.append(["pho2_PFisoC_endcap_highr9", "pho2_PFisoC", "(pho2_isEB == 0) * (pho2_r9 > 0.94) ","", 2., "", "(50,-1, 4)", "PFisoC^{#gamma2}", "endcap_highr9", 3.1, ""])
#
#plots.append(["pho2_sieie_barrel_lowr9", "pho2_sieie", "(pho2_isEB == 1) * (pho2_r9 < 0.94) ","", 2., "", "(120,0, 0.0120)", "sieie^{#gamma2}", "barrel_lowr9", 0.0102, ""])
#plots.append(["pho2_sieie_barrel_highr9", "pho2_sieie", "(pho2_isEB == 1) * (pho2_r9 > 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(25,0.007, 0.0120)", "sieie^{#gamma2}", "barrel_highr9", 0.0108, ""])
#plots.append(["pho2_sieie_endcap_lowr9", "pho2_sieie", "(pho2_isEB == 0) * (pho2_r9 < 0.94) ","", 2., "", "(100,0.02, 0.03)", "sieie^{#gamma2}", "endcap_lowr9", 0.028, ""])
#plots.append(["pho2_sieie_endcap_highr9", "pho2_sieie", "(pho2_isEB == 0) * (pho2_r9 > 0.94) ","", 2., "", "(100,0.02, 0.03)", "sieie^{#gamma2}", "endcap_highr9", 0.028, ""])
#
#plots.append(["pho2_hoe_barrel_lowr9", "pho2_hoe", "(pho2_isEB == 1) * (pho2_r9 < 0.94) ","", 2., "", "(100,0, 0.1)", "hoe^{#gamma2}", "barrel_lowr9", 0.092, ""])
#plots.append(["pho2_hoe_barrel_highr9", "pho2_hoe", "(pho2_isEB == 1) * (pho2_r9 > 0.94) ","(ph1_ciclevel >= 4) * (ph2_ciclevel >= 4)", 2., "", "(28, 0, 0.14)", "hoe^{#gamma2}", "barrel_highr9", 0.124, ""])
#plots.append(["pho2_hoe_endcap_lowr9", "pho2_hoe", "(pho2_isEB == 0) * (pho2_r9 < 0.94) ","", 2., "", "(80,0, 0.08)", "hoe^{#gamma2}", "endcap_lowr9", 0.063, ""])
#plots.append(["pho2_hoe_endcap_highr9", "pho2_hoe", "(pho2_isEB == 0) * (pho2_r9 > 0.94) ","", 2., "", "(100, 0.05, 0.150)", "hoe^{#gamma2}", "endcap_highr9", 0.142, ""])
#
#plots.append(["pho2_r9_barrel_lowr9", "pho2_r9", "(pho2_isEB == 1) * (pho2_r9 < 0.94) ","", 2., "", "(90,0.1, 1)", "r9^{#gamma2}", "barrel_lowr9", 0.298, ""])
#plots.append(["pho2_r9_barrel_highr9", "pho2_r9", "(pho2_isEB == 1) * (pho2_r9 > 0.94) ","", 2., "", "(100,0.8, 1)", "r9^{#gamma2}", "barrel_highr9", 0.94, ""])
#plots.append(["pho2_r9_endcap_lowr9", "pho2_r9", "(pho2_isEB == 0) * (pho2_r9 < 0.94) ","", 2., "", "(90,0.1, 1)", "r9^{#gamma2}", "endcap_lowr9", 0.24, ""])
#plots.append(["pho2_r9_endcap_highr9", "pho2_r9", "(pho2_isEB == 0) * (pho2_r9 > 0.94) ","", 2., "", "(100, 0.8, 1)", "r9^{#gamma2}", "endcap_highr9", 0.28, ""])
#
#plots.append(["toto", "pho2_PFisoB : pho2_PFisoA", "","", 2., "", "(20,-4, 16, 20,-4, 16)", "PFisoA^{#gamma2}", "", "", ""])

#plots.append(["pho2_PFisoC", "pho2_PFisoC", "","ph2_ciclevel >= 4", 2., "", "(24,0, 6)", "PFisoC^{#gamma2}", "", "", ""])
#plots.append(["pho1_isconv", "pho1_isconv", "","", 1., "", "(100,0, 1)", "isconv^{#gamma1}", "", "", ""])
#plots.append(["pho2_isconv", "pho2_isconv", "","ph2_ciclevel >= 4", 2., "evweight", "(50,0, 1)", "isconv^{#gamma2}", "", "", ""])
#plots.append(["jet1_csvBtag", "jet1_csvBtag", "", "", 1., "", "(50, 0.5, 1)", "csvBtag^{jet1}", "", "", ""])
#plots.append(["jet2_csvBtag", "jet2_csvBtag", "", "", 1., "", "(50, 0.5, 1)", "csvBtag^{jet2}", "", "", ""])

#plots.append(["gr_b_DeltaR_min_jet1_cat1", "gr_b_DeltaR_min_jet1", "njets_kRadionID_and_CSVM == 1", "", 2., "evweight", "(50, 0, 1)", "DeltaR(jet1,bquark)", "one btagged jet", 0.4, ""])
#plots.append(["DeltaR_b_jet2_cat1", "gr_b_DeltaR_min_jet2", "njets_kRadionID_and_CSVM == 1", "", 2., "evweight", "(50, 0, 1)", "DeltaR(jet2,bquark)", "one btagged jet", 0.4, ""])
#plots.append(["gr_b_DeltaR_min_jet1_cat0", "gr_b_DeltaR_min_jet1", "njets_kRadionID_and_CSVM >= 2", "", 2., "evweight", "(50, 0, 1)", "DeltaR(jet1,bquark)", "two btagged jet", 0.4, ""])
#plots.append(["gr_b_DeltaR_min_jet2_cat0", "gr_b_DeltaR_min_jet2", "njets_kRadionID_and_CSVM >= 2", "", 2., "evweight", "(50, 0, 1)", "DeltaR(jet2,bquark)", "two btagged jet", 0.4, ""])

#plots.append(["DeltaR_pho1_jet_min", "DeltaR_pho1_jet_min", "", "", 2., "evweight", "(50, 0, 5)", "DeltaR_pho1_jet_min", "", "", ""])
#plots.append(["DeltaR_pho2_jet_min", "DeltaR_pho2_jet_min", "", "ph2_ciclevel >= 4", 2., "evweight", "(40, 0 , 4)", "DeltaR_pho2_jet_min", "", "", ""])
#plots.append(["DeltaR_pho1_pho2", "DeltaR_pho1_pho2", "", "(ph2_ciclevel >= 4) * (ph1_ciclevel >= 4)", 2., "evweight", "(20, 0 , 5)", "DeltaR(pho1,pho2)", "", "", ""])

##plots.append(["aaaaa_test_2D", "pho2_PFisoA : pho2_PFisoB", "", "", 2., "evweight", "(100, -5, 20, 100, -5, 20)", "pho1_PFisoA", "", 33.3, ""])
#plots.append(["mtot_cat0_" + CIC1 +"_"+ CIC2, "mtot", "(cut_based_ct == 0)", "", "", "evWeight", "(25, 350, 450)", "mtot (GeV)", "m= " + MASSE +"GeV", "", ""])
#plots.append(["mtot_cat1_" + CIC1 +"_"+ CIC2, "mtot", "(cut_based_ct == 1)", "", "", "evWeight", "(25, 350, 450)", "mtot (GeV)", "m= " + MASSE +"GeV", "", ""])
#plots.append(["mtot_" + CIC1 +"_"+ CIC2, "mtot", "", "", "", "evWeight", "(25, 350, 450)", "mtot (GeV)", "m= " + MASSE +"GeV", "", ""])

#plots.append(["mjj_cat0_" + CIC1 +"_"+ CIC2, "mjj", "(cut_based_ct == 0)", "", "", "evWeight", "(20, 80, 160)", "mjj (GeV)", "m= " + MASSE +"GeV", "", ""])
#plots.append(["mjj_cat1_" + CIC1 +"_"+ CIC2, "mjj", "(cut_based_ct == 1)", "", "", "evWeight", "(20, 80, 160)", "mjj (GeV)", "m= " + MASSE +"GeV", "", ""])
#plots.append(["mjj_" + CIC1 +"_"+ CIC2, "mjj", "", "", "", "evWeight", "(20, 80, 160)", "mjj (GeV)", "m= " + MASSE +"GeV", "", ""])

#plots.append(["mgg_cat0_" + CIC1 +"_"+ CIC2, "mgg", "(cut_based_ct == 0)", "", "", "evWeight", "(25, 90, 190)", "mgg (GeV)", "m= " + MASSE +"GeV", "", ""])
#plots.append(["mgg_cat1_" + CIC1 +"_"+ CIC2, "mgg", "(cut_based_ct == 1)", "", "", "evWeight", "(25, 90, 190)", "mgg (GeV)", "m= " + MASSE +"GeV", "", ""])
#plots.append(["mgg_" + CIC1 +"_"+ CIC2, "mgg", "", "", "", "evWeight", "(25, 90, 190)", "mgg (GeV)", "m= " + MASSE +"GeV", "", ""])
#plots.append(["ggjj_mass", "ggjj_mass", "", "","",  "evweight", "(50, 100, 700)", "m_ggjj (GeV)", "", "", ""])
Eff = []

for name2, variable, cut, critere_eff,  norm, Scale_Factor, binning, title, additional_info, cutline, cutline2 in plots:
    c1 = TCanvas()
    legend = TLegend(0.45, 0.82, 0.90, 0.93, "")
    legend.SetTextSize(0.025)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    if len(samples) > 1:
        legend.SetNColumns(2)
    xnbin, xlow, xhigh = map(float, binning.strip().strip("()").split(","))

    ynbin = 0
    try :
        xnbin, xlow, xhigh = map(float, binning.strip().strip("()").split(","))
    except :
        xnbin, xlow, xhigh, ynbin, ylow, yhigh = map(float, binning.strip().strip("()").split(","))

    ymax = -1
    ymin = 10000000
    firsthistname = ""
    if cut == "": cut = "1"


    for ifile, [ name, dirpath, subdir, file, tree, color, style, label , sigma , N] in enumerate(samples):
#        print ifile, file, color, style, label
        chain = TChain(tree)
        chain.Add( path.join(dirpath, subdir, file) )
        sample_cut = cut
#        print sample_cut
        print file
#        print chain.GetEntries(cut)
        if chain.GetEntries(cut)==0:
            continue
#        print chain.GetEntries(cut)
        if norm == 1.:
            sample_cut = "(" + sample_cut + ")/" + str( chain.GetEntries() )
        elif norm== 2.:
            sample_cut = "(" + sample_cut + ") * ((" + str(sigma) + " * " + str(intL) + ")/" + str(N) + ")"
        
        option = ""
        if ynbin != 0:
            gStyle.SetPalette(1)
            option = "COLZ"


        if Scale_Factor != "":        
            sample_cut += "* (" + Scale_Factor + ")"   
#        print sample_cut

        latexLabel = TLatex()
        latexLabel.SetTextSize(.03)
        latexLabel.SetNDC()
#        latexLabel.DrawLatex(.20, .85, additional_info)
        ROOT.gPad.RedrawAxis()
#        legend.Draw()
        c1.Update()


        if ifile != 0:
            option = "same"

        chain.Draw(variable + ">>h_tmp_ntot" + binning, sample_cut , "goff")

        # Cosmetics
        h_ntot = ROOT.gDirectory.Get("h_tmp_ntot")
        chain.Draw(variable + ">>h_tmp" + binning, sample_cut + "/" + str(h_ntot.GetMaximum()), option)
        h = ROOT.gDirectory.Get("h_tmp")
        h.SetName(name + "_" + name2 + "_" + str(ifile))
#        print ifile
        if ifile ==0:
            ncut = h.Integral()
        else:
            nwp = h.Integral()
        if ifile == 0:
                firsthistname = name + "_" + name2 + "_" + str(ifile)
        h.SetLineWidth(1)
        h.SetLineColor(color)
        h.GetXaxis().SetTitle( title )
        print type(h)
        unit = ""
        if title.find("(") != -1:
            unit = title[title.find("(")+1:title.find(")")]
        if norm == 1.:
            h.GetYaxis().SetTitle( "Norm. to unity / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        else:
            h.GetYaxis().SetTitle( "# events / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        legend.AddEntry(h.GetName(), label)
        ymax = max(ymax, h.GetMaximum())
#        print "ymax = ", ymax
        ymin = min(ymin, h.GetMinimum(0.0))
#        print "ymin = ", ymin

#calcul ntot, nselect pour calcul efficacite
#        if critere_eff == "":
#            critere_eff = "1"
#        chain.Draw(variable + ">>h_tmp2" + binning, sample_cut + "*(" + critere_eff  + ")" , "goff")
#        h2 = ROOT.gDirectory.Get("h_tmp2")
#        h2.Divide(h_ntot) 
#        h2.SetLineColor(4)
#        h2.Draw("same")   
#        
#        ntot=chain.GetEntries(sample_cut)
#        nselected=chain.GetEntries("(" + sample_cut  + ") * (" + critere_eff + ")")
#        eff = nselected / float(ntot)
#        error = m.sqrt(eff*(1-eff)/ntot)

# calcul correlation entre variable

#        chain.Draw(variable + ">>h_tmp3" + binning, sample_cut, "entrylist")
#        h3 = ROOT.gDirectory.Get("h_tmp3")

#        corr_A_B=0.
#        mean_A=0.
#        mean_A_A=0.
#        mean_B=0.
#        mean_B_B=0.
#        for ievt in xrange(chain.GetEntries()):
#            chain.GetEntry(ievt)
#            if (chain.njets_kRadionID_and_CSVM >= 2):
#               mean_A += chain.jet1_betaStarClassic / ntot
#                mean_A_A += chain.jet1_betaStarClassic * chain.jet1_betaStarClassic / ntot
#                mean_B += chain.jet1_dR2Mean / ntot
#                mean_B_B += chain.jet1_dR2Mean * chain.jet1_dR2Mean / ntot
#                corr_A_B += (chain.jet1_betaStarClassic * chain.jet1_dR2Mean) / ntot
#        sigma_A = m.sqrt(mean_A_A - mean_A * mean_A) 
#        sigma_B = m.sqrt(mean_B_B - mean_B * mean_B)
#        sigma_A_B = (corr_A_B - (mean_A * mean_B)) / ( sigma_A * sigma_B) 
#        print sigma_A_B

#        pho = 0
#        if variable == "DeltaR_gr_pho_pho1_min":
#            pho = 1
#        else:
#            pho = 2
#
#        nocut = "tagada"
#        if color == ROOT.kBlack:
#            nocut = "$\ge$ 0 for pho 1, $\ge$ 0 for pho 2"
#        elif color == ROOT.kRed:
#            nocut = "$\ge$ 1 for pho 1, $\ge$ 0 for pho 2"
#        elif color == ROOT.kBlue:
#            nocut = "$\ge$ 2 for pho 1, $\ge$ 0 for pho 2"
#        elif color == ROOT.kGreen:
#            nocut = "$\ge$ 3 for pho 1, $\ge$ 0 for pho 2"
#        elif color == ROOT.kCyan:
#            nocut = "$\ge$ 4 for pho 1, $\ge$ 0 for pho 2"
#	                   
#	    mass = 0
#    	if name == "Radion_m270":
#		    mass = 270
#        elif name == "Radion_m300":
#            mass = 300
#        elif name == "Radion_m350":
#            mass = 350
#        elif name == "Radion_m400":
#            mass = 400
#
#        Eff.append([mass, pho , nocut , nselected , ntot, eff, error])


        del chain, h



    ymin_lin = ymin / 10.
    yrange_lin = ymax - ymin_lin
    ymax_lin = .25 * yrange_lin + ymax
    yrange_log = (log10(ymax) - log10(ymin)) / .77
    ymax_log = pow(10., .25*yrange_log + log10(ymax))
    ymin_log = pow(10., log10(ymin) - .03*yrange_log)

    latexLabel = TLatex()
    latexLabel.SetTextSize(.025)
    latexLabel.SetNDC()
    latexLabel.DrawLatex(.3, .96, "CMS Internal     L = 19.7 fb^{-1}     #sqrt{s} = 8 TeV   " + additional_info)
#    latexLabel.DrawLatex(.15, .9, "n WP = %s" %round(nwp,1))
#    latexLabel.DrawLatex(.15, .85, "n = " + str(round(ncut, 1)) + " ( + " + str(round(100 * ncut/nwp - 100 ,1)) + "%)")
#    latexLabel.DrawLatex(.2, .85, additional_info)
#    latexLabel.DrawLatex(.2, .85, "n evts = " + str(ntot) + " ,eff = " + str(round(eff,2)))
#    latexLabel.DrawLatex(.45,.75, "blue : proportion of photon with CIC >= 4" )
    ROOT.gPad.RedrawAxis()
    legend.Draw()
    c1.Update()

    line = TLine()
    line.SetLineStyle(2)
    line.SetLineWidth(1)
    line2 = TLine()
    line2.SetLineStyle(2)
    line2.SetLineWidth(1)

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

#pouet = "\\"
#for mass, pho, cut, nselected, ntot, eff, error in Eff:
#    print "      ", cut, "&", round(ntot,1) , pouet + pouet
#    print "      ", "\hline"
#tagadapouet
