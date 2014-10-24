#!/usr/bin/env python
from array import array
import ROOT
from ROOT import gROOT
from ROOT import TFile, TH2F
gROOT.SetBatch()

infile = TFile.Open("2014-10-21_selection_noRegression_noMassCut_v17/Data_noRegression_noMassCut_controlSample_v17.root")
intree = infile.Get("Data")
outfile = TFile("2014-10-21_selection_noRegression_noMassCut_v17/Data_noRegression_noMassCut_controlSampleWeighted_v17.root", "RECREATE")

intree.SetBranchStatus("*", 1)
intree.SetBranchStatus("evweight_w_btagSF", 0)
outtree = intree.CloneTree(0)
weight = array('f', [1.])
outtree.Branch("evweight_w_btagSF", weight, "evweight_w_btagSF/F")

csWeightFile = TFile.Open("scales_2D_pt_data_4GeVbinning.root");
h2D_pt_data = csWeightFile.Get("h2D_pt_data");

for ievt in xrange(intree.GetEntries()):
    intree.GetEntry( ievt )
    weight[0] = h2D_pt_data.GetBinContent(h2D_pt_data.FindBin(intree.pho2_pt, intree.pho1_pt))
    outtree.Fill()

outfile.cd()
outtree.Write()
outfile.Close()
infile.Close()
csWeightFile.Close()
