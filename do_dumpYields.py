#!/usr/bin/env python
from os import path
import ROOT
from ROOT import gROOT, gDirectory
from ROOT import TChain
gROOT.SetBatch()

afs_path = "/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/"

samples = []
#samples.append( [label, samplePath, treeName, weight] )
samples.append( ["Data v17", "v17/2014-10-21_selection_noRegression_noMassCut_v17/Data_noRegression_noMassCut_v17.root", "Data", "1"] )
samples.append( ["Data v18", "v18/2014-10-22_selection_withRegression_noMassCut_v18/Data_withRegression_noMassCut_v18.root", "Data", "1"] )

samples.append( ["ggHH v17", "v17/2014-10-21_selection_noRegression_noMassCut_v17/ggHH_8TeV_noRegression_noMassCut_v17.root", "ggHH_8TeV", "1"] )
samples.append( ["ggHH v18", "v18/2014-10-22_selection_withRegression_noMassCut_v18/ggHH_8TeV_withRegression_noMassCut_v18.root", "ggHH_8TeV", "1"] )

samples.append( ["diphojet v17", "v17/2014-10-21_selection_noRegression_noMassCut_v17/diphojet_sherpa_8TeV_noRegression_noMassCut_v17.root", "diphojet_sherpa_8TeV", "1"] )
samples.append( ["diphojet v18", "v18/2014-10-22_selection_withRegression_noMassCut_v18/diphojet_sherpa_8TeV_withRegression_noMassCut_v18.root", "diphojet_sherpa_8TeV", "1"] )

samples.append( ["radion300 v17", "v17/2014-10-21_selection_noRegression_noMassCut_v17/Radion_m300_8TeV_noRegression_noMassCut_v17.root", "Radion_m300_8TeV", "1"] )
samples.append( ["radion300 v18", "v18/2014-10-22_selection_withRegression_noMassCut_v18/Radion_m300_8TeV_withRegression_noMassCut_v18.root", "Radion_m300_8TeV", "1"] )

cuts = []
#cuts.append( [cutname, cutstring] )
cuts.append( ["cat0", "njets_kRadionID_and_CSVM == 1"] )
cuts.append( ["cat1", "njets_kRadionID_and_CSVM >= 2"] )

for cutname, cut in cuts:
    print "#####", cutname, "#####"
    for label, file, tree, weight in samples:
        chain = TChain(tree)
        chain.Add( path.join(afs_path, file) )
        n = 0. 
        if weight == "" or weight == "1." or weight == "1":
            n = chain.GetEntries(cut)
        else:
            chain.Draw(">>elist", cut, "entrylist")
            elist = gDirectory.Get("elist")
            listEntries = elist.GetN()
            chainEntries = chain.GetEntries()
            treenum = 0 # Be careful with this if there is several trees in the chain
            chain.SetEntryList(elist)
#            print listEntries, chainEntries
            for ientry in xrange(listEntries):
                ievt = elist.GetEntryAndTree(ientry, ROOT.Long(treenum))
                chain.GetEntry( ievt )
                n += chain.evweight_w_btagSF
#        print file, cut, n
        print label, n

