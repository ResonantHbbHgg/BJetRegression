#!/usr/bin/env python
from os import path
from ROOT import gROOT
from ROOT import TChain
gROOT.SetBatch()

eos_prod = "/store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/"
eos_pro2 = "/store/caf/user/bmarzocc/"
#eos_pro2 = "/store/group/phys_higgs/cmshgg/V15_00_11/"
#eos_redu = "/store/group/phys_higgs/cmshgg/reduced/legacy_paper/legacy_paper_reduction_8TeV_v5/mc/"
eos_redu = "/store/cmst3/user/obondu/H2GGLOBE/Radion/reduced/radion_reduction_v11/mc"
eos_tree = "/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_11_tree_08"
samples = []
samples.append(["Radion_m270", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-270_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m270", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-270_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m270", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m270_8TeV"])

samples.append(["Radion_m300", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-300_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m300", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-300_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m300", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m300_8TeV"])

samples.append(["Radion_m350", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-350_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m350", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-350_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m350", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m350_8TeV"])

samples.append(["Radion_m400", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-400_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m400", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-400_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m400", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m400_8TeV"])

samples.append(["Radion_m450", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-450_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m450", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-450_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m450", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m450_8TeV"])

samples.append(["Radion_m500", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-500_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m500", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-500_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m500", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m500_8TeV"])

samples.append(["Radion_m550", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-550_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m550", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-550_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m550", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m550_8TeV"])

samples.append(["Radion_m600", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-600_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m600", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-600_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m600", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m600_8TeV"])

samples.append(["Radion_m650", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-650_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m650", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-650_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m650", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m650_8TeV"])

samples.append(["Radion_m700", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-700_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m700", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-700_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m700", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m700_8TeV"])

samples.append(["Radion_m800", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-800_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v3_AODSIM", "*root", "event"])
samples.append(["Radion_m800", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-800_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v3", "*root", "event"])
samples.append(["Radion_m800", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m800_8TeV"])

samples.append(["Radion_m900", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-900_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m900", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-900_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m900", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m900_8TeV"])

samples.append(["Radion_m1100", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-1100_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m1100", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-1100_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m1100", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m1100_8TeV"])

samples.append(["Radion_m1200", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-1200_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m1200", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-1200_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m1200", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m1200_8TeV"])

samples.append(["Radion_m1300", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-1300_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m1300", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-1300_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m1300", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m1300_8TeV"])

samples.append(["Radion_m1400", "prod", eos_pro2, "RadionToHH_2Gamma_2b_M-1400_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM", "*root", "event"])
samples.append(["Radion_m1400", "redu", eos_redu, "RadionToHH_2Gamma_2b_M-1400_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m1400", "tree", eos_tree, "", "RadionToHH_2Gamma_2b.root", "Radion_m1400_8TeV"])

samples.append(["Radion_m300", "prod", eos_prod, "Summer12_DR53X-PU_RD1_START53_V7N/RadionToHHTo2G2B_M-300_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m300", "redu", eos_redu, "RadionToHHTo2G2B_M-300_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m300", "tree", eos_tree, "", "Radion_nm.root", "Radion_m300_8TeV_nm"])
samples.append(["Radion_m500", "prod", eos_prod, "Summer12_DR53X-PU_RD1_START53_V7N/RadionToHHTo2G2B_M-500_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m500", "redu", eos_redu, "RadionToHHTo2G2B_M-500_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m500", "tree", eos_tree, "", "Radion_nm.root", "Radion_m500_8TeV_nm"])
samples.append(["Radion_m700", "prod", eos_prod, "Summer12_DR53X-PU_RD1_START53_V7N/RadionToHHTo2G2B_M-700_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m700", "redu", eos_redu, "RadionToHHTo2G2B_M-700_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["Radion_m700", "tree", eos_tree, "", "Radion_nm.root", "Radion_m700_8TeV_nm"])

#samples.append(["prod", eos_pro2, "ExtraMaterial/GluGluToHToGG_M-125_StdMaterial", "*root", "event"])
#samples.append(["redu", eos_redu, "Summer12_RD1/GluGluToHToGG_M-125_StdMaterial", "*root", "event"])
samples.append(["ggH", "prod", eos_prod, "Summer12_RD1/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["ggH", "redu", eos_redu, "Summer12_RD1/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["ggH", "tree", eos_tree, "", "SMHiggs.root", "ggh_m125_powheg_8TeV"])
samples.append(["VBF", "prod", eos_prod, "Summer12_RD1/VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["VBF", "redu", eos_redu, "Summer12_RD1/VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["VBF", "tree", eos_tree, "", "SMHiggs.root", "vbf_m125_8TeV"])
samples.append(["WH_", "prod", eos_prod, "Summer12_RD1/WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2", "*root", "event"])
samples.append(["WH_", "redu", eos_redu, "Summer12_RD1/WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2_wh", "*root", "event"])
samples.append(["WH_","tree", eos_tree, "", "SMHiggs.root", "wzh_m125_8TeV_wh"])
samples.append(["ZH_","prod", eos_prod, "Summer12_RD1/WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2", "*root", "event"])
samples.append(["ZH_","redu", eos_redu, "Summer12_RD1/WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2_zh", "*root", "event"])
samples.append(["ZH_","tree", eos_tree, "", "SMHiggs.root", "wzh_m125_8TeV_zh"])
samples.append(["ttH","prod", eos_prod, "Summer12_RD1/TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["ttH","redu", eos_redu, "Summer12_RD1/TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
samples.append(["ttH","tree", eos_tree, "", "SMHiggs.root", "tth_m125_8TeV"])

nprocessed = {}
nreduced = {}
ntrees = {}
treename = {}
for displayName, level, eos_base, sample, file, tree in samples:
    chain = TChain(tree)
    files = "root://eoscms//eos/cms" + path.join(eos_base, sample, file)
    chain.Add(files)
    print "displayName= ", displayName, "level= ", level, "chain.GetEntries()= ", chain.GetEntries()
    if level == "prod":
        nprocessed[displayName] = chain.GetEntries()
    if level == "redu":
        nreduced[displayName] = chain.GetEntries()
    if level == "tree":
        ntrees[displayName] = chain.GetEntries()
        treename[displayName] = tree
#    if level == "tree":
#        nentries_w = 0
#        for ievt in xrange(chain.GetEntries()):
#            chain.GetEntry(ievt)
#            nentries_w += chain.evweight
#        print "displayName= ", displayName, "level= ", "trew", "chain.GetEntries()= ", nentries_w * nprocessed[displayName] / 19706.
    del chain

print "########################################################"
print "treeName, nprocessed, nreduced, ntrees"
print "########################################################"
displayNames = set([ x[0] for x in samples ])
for display in displayNames:
    print treename[display], nprocessed[display], nreduced[display], ntrees[display]
#for display in displayNames:
#    print treename[display], nprocessed[display] / float(nprocessed[display]) * 100, nreduced[display] / float(nprocessed[display]) * 100, ntrees[display] / float(nprocessed[display]) * 100

#samples = []
#samples.append(["root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_11_tree_08/RadionToHH_2Gamma_2b.root", "Radion_m300_8TeV"])
#
##sample = samples[0]
##tree = tree[0]
#for sample, tree in samples:
#    chain = TChain(tree)
#    chain.Add(sample)
#    totyield = 0
#    print "all events: ", chain.GetEntries()
#    for ievt in xrange(chain.GetEntries()):
##    for ievt in xrange(10):
#        chain.GetEntry(ievt)
#        nbjets = 0
#        for ijet in range(1,15):
#            if getattr(chain, "j" + str(ijet) + "_csvBtag") > 0.679:
#                nbjets += 1
##        print "nbjets= ", nbjets
#        if nbjets > 0:
#            totyield += 1
#    del chain
#
#print totyield
#
