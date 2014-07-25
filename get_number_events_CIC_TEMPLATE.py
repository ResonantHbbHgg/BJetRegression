#!/usr/bin/env python
# Various python imports
from os import path
from math import log10, pow
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1D, TLatex, TLegend, TLorentzVector, TStyle, gStyle
import math as m
dirpath_cut = "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/"
subdir_cut = "v31_study_ciclevel_CICLEVEL1_CICLEVEL2_fitToMgg_noKinFit"
file_cut= "TREE_mMASSE.root"
tree_cut = "TCVARS"

dirpath_cut_4_4 = "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/"
subdir_cut_4_4 = "v31_study_ciclevel_4_4_fitToMgg_noKinFit"
file_cut_4_4= "TREE_mMASSE.root"
tree_cut_4_4 = "TCVARS"

n_cut=0.
n_cut_cat0=0.
n_cut_cat1=0.
n_cut_4_4=0.
n_cut_4_4_cat0=0.
n_cut_4_4_cat1=0.

chain_cut = TChain(tree_cut)
chain_cut.Add( path.join(dirpath_cut, subdir_cut, file_cut) )
for ievt in xrange(chain_cut.GetEntries()):
    chain_cut.GetEntry(ievt)
    if (chain_cut.cut_based_ct == 1):
        n_cut_cat1+=chain_cut.evWeight
        n_cut+=chain_cut.evWeight
    elif (chain_cut.cut_based_ct == 0):
        n_cut_cat0+=chain_cut.evWeight
        n_cut+=chain_cut.evWeight

chain_cut_4_4 = TChain(tree_cut_4_4)
chain_cut_4_4.Add( path.join(dirpath_cut_4_4, subdir_cut_4_4, file_cut_4_4) )
for ievt in xrange(chain_cut_4_4.GetEntries()):
    chain_cut_4_4.GetEntry(ievt)
    if (chain_cut_4_4.cut_based_ct == 1):
        n_cut_4_4_cat1+=chain_cut_4_4.evWeight
        n_cut_4_4+=chain_cut_4_4.evWeight
    if (chain_cut_4_4.cut_based_ct == 0):
        n_cut_4_4_cat0+=chain_cut_4_4.evWeight
        n_cut_4_4+=chain_cut_4_4.evWeight



eff_cat0=n_cut_cat0/n_cut_4_4_cat0
eff_cat1=n_cut_cat1/n_cut_4_4_cat1
eff=n_cut/n_cut_4_4

print "  ", "CICLEVEL1", "&", "CICLEVEL2","&", round(n_cut, 1), "&", round(100*(eff-1), 1) , "$\%$", "&", round(n_cut_cat0, 1), "&", round(100*(eff_cat0-1), 1) , "$\%$", "&", round(n_cut_cat1, 1),"&",round(100*(eff_cat1-1), 1), "$\%$", "\\" +"\\"

#print "  ", "\hline"
