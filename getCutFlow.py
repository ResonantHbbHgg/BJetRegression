#!/usr/bin/env python
# Modules imported
import argparse
import numpy
from os import path, listdir
# Argument parser setup
parser = argparse.ArgumentParser()
parser.add_argument("--part", help="1= prod and reduced level, 2= object selection, 3= mass windows, 999= get everything together", nargs='?', const=999, type=int, default=999)
args = parser.parse_args()
# ROOT setup
import ROOT
from ROOT import gROOT, gStyle
from ROOT import TChain, TCanvas, TGraph, TLatex, TLegend, TPad, TGaxis
gROOT.Reset()
gROOT.SetBatch()

if args.part == 1:
    print "part 1: preselection {produced, reduced, optree}-level"
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
    samples.append(["Radion_m1000", "prod", eos_prod, "Summer12_DR53X-PU_RD1_START53_V7N/RadionToHHTo2G2B_M-1000_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
    samples.append(["Radion_m1000", "redu", eos_redu, "RadionToHHTo2G2B_M-1000_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
    samples.append(["Radion_m1000", "tree", eos_tree, "", "Radion_nm.root", "Radion_m1000_8TeV_nm"])
    samples.append(["Radion_m1500", "prod", eos_prod, "Summer12_DR53X-PU_RD1_START53_V7N/RadionToHHTo2G2B_M-1500_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
    samples.append(["Radion_m1500", "redu", eos_redu, "RadionToHHTo2G2B_M-1500_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event"])
    samples.append(["Radion_m1500", "tree", eos_tree, "", "Radion_nm.root", "Radion_m1500_8TeV_nm"])
    
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
#        print "displayName= ", displayName, "level= ", level, "chain.GetEntries()= ", chain.GetEntries()
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
    
#    print "########################################################"
#    print "treeName, nprocessed, nreduced, ntrees"
#    print "########################################################"
    displayNames = set([ x[0] for x in samples ])
    for display in displayNames:
        print treename[display], nprocessed[display], nreduced[display], ntrees[display]
    #for display in displayNames:
    #    print treename[display], nprocessed[display] / float(nprocessed[display]) * 100, nreduced[display] / float(nprocessed[display]) * 100, ntrees[display] / float(nprocessed[display]) * 100
    
if args.part == 2:
    print "part 2: object selection"
    nflow = {}
    cutFlowDir = "."
    cutFlowFiles = [ x for x in listdir(cutFlowDir) if path.isfile(path.join(cutFlowDir,x)) and "cutFlow" in x and ".dat" in x and "part" not in x]
    for file in cutFlowFiles:
        sample = file.replace(".dat", "").replace("cutFlow_", "")
        nflow[sample] = []
        with open(file) as data:
            for line in data:
                nflow[sample].append(line.split("\t")[1])
    for sample in nflow:
        print sample, ", ".join(nflow[sample])
                
        

if args.part == 3:
    print "part 3: mass windows"
    mggjj_cut = {}
    mggjj_cut[260] = [225, 280]
    mggjj_cut[270] = [225, 295]
    mggjj_cut[300] = [255, 330]
    mggjj_cut[350] = [310, 395]
    mggjj_cut[400] = [370, 440]
    cutFlowDir = "."
    selectedFolder = "2014-04-16_selection_noRegression_noMassCut_v10_cutflow/"
    cutFlowFiles = [ x for x in listdir(cutFlowDir) if path.isfile(path.join(cutFlowDir,x)) and "cutFlow" in x and ".dat" in x and "part" not in x]
    for file in cutFlowFiles:
        sample = file.replace(".dat", "").replace("cutFlow_", "")
        if "Radion" in sample and ("m300" in sample or "m500" in sample or "m700" in sample or "m1000" in sample or "m1500" in sample):
            sample = sample.replace("_nm", "")
        mass = int(sample.split("_")[1].strip('m'))
        chain = TChain(sample)
        files = path.join(selectedFolder, sample + "_noRegression_noMassCut_v10_cutflow.root")
        chain.Add(files)
# FIXME
        if not "Radion" in sample:
            masses = [260, 270, 300, 350, 400]
            for mass in masses:
                print sample + "_M" + str(mass), chain.GetEntries(), chain.GetEntries("jj_mass > 85 && jj_mass < 155"), chain.GetEntries("jj_mass > 85 && jj_mass < 155 && ggjj_mass > " + str(mggjj_cut[mass][0]) + " && ggjj_mass < " + str(mggjj_cut[mass][1]))
        else:
            if mass <= 400:
                print sample, chain.GetEntries(), chain.GetEntries("jj_mass > 85 && jj_mass < 155"), chain.GetEntries("jj_mass > 85 && jj_mass < 155 && ggjj_mass > " + str(mggjj_cut[mass][0]) + " && ggjj_mass < " + str(mggjj_cut[mass][1]))
            if mass >= 400:
                print sample, chain.GetEntries(), chain.GetEntries("jj_mass > 90 && jj_mass < 165"), chain.GetEntries("jj_mass > 90 && jj_mass < 165 && gg_mass > 120 && gg_mass < 130")
        
        
 
if args.part == 999:
    print "part 999: get everything together"
# Get all numbers
    flow_low = {}
    flow_high = {}
    flow_hgg = {}
    was400lowDone = False
    for ifile in range(1,4):
        print ifile
        file = "cutFlow_part" + str(ifile) + ".txt"
        with open(file) as data:
            for line in data:
                if "Radion" not in line and "m125" not in line:
                    continue
                sline = line.replace(",", "").split()
                sample = sline[0]
#                print len(sline), sample, sline
                sample = sample.replace("_nm", "")
                mass = int(sample.split("_")[1].split("m")[1])
#                if "m125" in sample:
#                    process = sample.split("_m125")[0] + sample.split("m125")[1].split("8TeV")[1]
#                    base = sample.split("_M")[0]
#                    print sample, process
                if mass > 1100:
                    continue
                if mass == 125 and "Radion" not in sample and sample not in flow_hgg:
                    flow_hgg[sample] = []
                if 200 < mass <= 400 and sample not in flow_low:
                    flow_low[sample] = []
                if mass >= 400 and sample not in flow_high:
                    flow_high[sample] = []
                if mass == 400:
                    print mass, sline, was400lowDone, not was400lowDone
                if mass == 125 and "Radion" not in sample:
                    for iflow_hgg in range(1, len(sline)):
                        if ifile == 2 and (iflow_hgg == 1 or iflow_hgg == 2 or iflow_hgg == 3 or iflow_hgg == 6 or iflow_hgg == 7 or iflow_hgg == 9):
                            continue
                        if ifile == 3 and (iflow_hgg == 1 or iflow_hgg == 2):
                            continue
                        flow_hgg[sample].append( int(sline[iflow_hgg]) )
                if 200 < mass < 400 or (mass == 400 and not was400lowDone and ifile == 3) or (mass == 400 and ifile != 3):
                    for iflow_low in range(1, len(sline)):
                        if ifile == 2 and (iflow_low == 1 or iflow_low == 2 or iflow_low == 3 or iflow_low == 6 or iflow_low == 7 or iflow_low == 9):
                            continue
                        if ifile == 3 and (iflow_low == 1 or iflow_low == 2):
                            continue
                        flow_low[sample].append( int(sline[iflow_low]) )
                if mass == 400 and not was400lowDone and ifile == 3:
                    was400lowDone = True
                    continue
                if mass > 400 or (mass == 400 and was400lowDone and ifile == 3) or (mass == 400 and ifile != 3):
                    for iflow_high in range(1, len(sline)):
                        if ifile == 2 and (iflow_high == 1 or iflow_high == 2 or iflow_high == 3 or iflow_high == 6 or iflow_high == 7 or iflow_high == 9):
                            continue
                        if ifile == 3 and (iflow_high == 1 or iflow_high == 2):
                            continue
                        flow_high[sample].append( int(sline[iflow_high]) )

# organize the numbers
    n_hgg = {}
    x_hgg = {}
    y_hgg = {}
    gr_hgg = {}
    n_low = {}
    x_low = {}
    y_low = {}
    gr_low = {}
    n_high = {}
    x_high = {}
    y_high = {}
    gr_high = {}
    for igraph in range(7):
        n_hgg[igraph] = 0
        x_hgg[igraph] = []
        y_hgg[igraph] = []
        n_low[igraph] = 0
        x_low[igraph] = []
        y_low[igraph] = []
        n_high[igraph] = 0
        x_high[igraph] = []
        y_high[igraph] = []

    base_set = []
    for key in flow_hgg:
        base = key.split("_M")[0]
        base_set.append(base)
        if key not in base:
            flow_hgg[key] = flow_hgg[base] + flow_hgg[key]
#            print base, key, flow_hgg[key]
    base_set = set(base_set)
    for base in base_set:
        del flow_hgg[base]
    offset = {}
    offset["ggh"] = -125
    offset["vbf"] = -124
    offset["_wh"] = -123
    offset["_zh"] = -122
    offset["tth"] = -121
    sigma_wh = 0.7046
    sigma_zh = 0.4153
    w = {}
    w["ggh"] = 1.
    w["vbf"] = 1.
    w["_wh"] = sigma_wh / (sigma_wh + sigma_zh)
    w["_zh"] = sigma_zh / (sigma_wh + sigma_zh)
    w["tth"] = 1.
    for key in flow_hgg:
        if "M300" not in key:
            continue
        off = 0
        wgt = 1.
        for koff in offset:
            if koff in key:
                off = int(offset[koff])
                wgt = w[koff]
        mass = int(key.split("_")[1].split("m")[1])
        print key, off, flow_hgg[key]
        for igraph in range(7):
            x_hgg[igraph].append(mass + off )
            eff = float(flow_hgg[key][igraph]) / float( flow_hgg[key][0] * wgt ) * 100.
            y_hgg[igraph].append( eff )
            n_hgg[igraph] += 1
    for key in flow_low:
        mass = int(key.split("_")[1].split("m")[1])
        print key, flow_low[key]
        for igraph in range(7):
            x_low[igraph].append(mass / 1000.)
            eff = float(flow_low[key][igraph]) / float(flow_low[key][0]) * 100.
            y_low[igraph].append( eff )
            n_low[igraph] += 1
    for key in flow_high:
        mass = int(key.split("_")[1].split("m")[1])
        print key, flow_high[key]
        for igraph in range(7):
            x_high[igraph].append(mass / 1000.)
            eff = float(flow_high[key][igraph]) / float(flow_high[key][0]) * 100.
            y_high[igraph].append( eff )
            n_high[igraph] += 1


    for igraph in range(7):
        y_hgg[igraph] = [b for (a,b) in sorted(zip(x_hgg[igraph], y_hgg[igraph]))]
        x_hgg[igraph] = sorted(x_hgg[igraph])
        x_hgg[igraph] = map(float, x_hgg[igraph])
        y_hgg[igraph] = map(float, y_hgg[igraph])
        x_hgg[igraph] = numpy.asarray(x_hgg[igraph], dtype='float')
        y_hgg[igraph] = numpy.asarray(y_hgg[igraph], dtype='float')
        y_low[igraph] = [b for (a,b) in sorted(zip(x_low[igraph], y_low[igraph]))]
        x_low[igraph] = sorted(x_low[igraph])
        x_low[igraph] = map(float, x_low[igraph])
        y_low[igraph] = map(float, y_low[igraph])
        x_low[igraph] = numpy.asarray(x_low[igraph], dtype='float')
        y_low[igraph] = numpy.asarray(y_low[igraph], dtype='float')
        y_high[igraph] = [b for (a,b) in sorted(zip(x_high[igraph], y_high[igraph]))]
        x_high[igraph] = sorted(x_high[igraph])
        x_high[igraph] = map(float, x_high[igraph])
        y_high[igraph] = map(float, y_high[igraph])
        x_high[igraph] = numpy.asarray(x_high[igraph], dtype='float')
        y_high[igraph] = numpy.asarray(y_high[igraph], dtype='float')

# Now with the plot
    xAxisMin = .250
    xAxisMax = 1.140
    canvasSplit = 30.
    wPad = .82
    gROOT.ProcessLine(".x setTDRStyleMP.C")
    TGaxis.SetMaxDigits(3)
    color = [ROOT.kRed, ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan+3, ROOT.kGreen+2, ROOT.kYellow+3]
    c1 = TCanvas()
    pad1 = TPad("pad1", "pad1", 0., 0., canvasSplit / 100., 1.)
    pad2 = TPad("pad2", "pad2", canvasSplit / 100., 0., 1., 1.)
    pad1.Draw()
    pad2.Draw()
#    gStyle.SetPadBorderMode(0)
#    gStyle.SetFrameBorderMode(0)
    small = .00001
#    c1.Divide(2,1,small,small)
    for igraph in range(1,7):
        pad2.cd()
#        c1.cd(2)
        pad2.SetLeftMargin(small);
        pad2.SetRightMargin(0.05 / (100. * wPad - canvasSplit) * 100.)
        gr_low[igraph] = TGraph(n_low[igraph], x_low[igraph], y_low[igraph])
        gr_low[igraph].SetName("low_" + str(igraph))
        gr_low[igraph].SetMarkerStyle(20)
        gr_low[igraph].SetMarkerColor( color[igraph -1] )
        gr_low[igraph].SetLineColor( color[igraph -1] )
        gr_low[igraph].SetTitle("")
        gr_low[igraph].GetXaxis().SetTitle("m_{X} (TeV)")
        gr_low[igraph].GetYaxis().SetTickLength(0.03 / (100. - canvasSplit) * 100.) # due to assymmetric canvas...
        gr_low[igraph].SetMinimum(0.)
        gr_low[igraph].SetMaximum(100.)
        gr_low[igraph].GetXaxis().SetLimits(xAxisMin, xAxisMax)
        if igraph == 1:
            gr_low[igraph].Draw("ALP")
        else: 
            gr_low[igraph].Draw("LP")
        gr_high[igraph] = TGraph(n_high[igraph], x_high[igraph], y_high[igraph])
        gr_high[igraph].SetName("low_" + str(igraph))
        gr_high[igraph].SetMarkerStyle(21)
        gr_high[igraph].SetMarkerColor( color[igraph -1] )
        gr_high[igraph].SetLineColor( color[igraph -1] )
        gr_high[igraph].SetTitle("")
        gr_high[igraph].GetXaxis().SetTitle("m_{X} (TeV)")
        gr_high[igraph].GetXaxis().SetTickLength(0.03 / (100. - canvasSplit) * 100.) # due to assymmetric canvas...
        gr_high[igraph].GetXaxis().SetNdivisions(510)
        gr_high[igraph].GetYaxis().SetTickLength(0.03 / (100. - canvasSplit) * 100.) # due to assymmetric canvas...
        gr_high[igraph].SetMinimum(0.)
        gr_high[igraph].SetMaximum(100.)
        gr_high[igraph].GetXaxis().SetLimits(xAxisMin, xAxisMax)
        gr_high[igraph].Draw("LP")
#        c1.cd(1)
        pad1.cd()
        pad1.SetRightMargin(small)
        pad1.SetLeftMargin(0.13 / canvasSplit * 100.)
        gr_hgg[igraph] = TGraph(n_hgg[igraph], x_hgg[igraph], y_hgg[igraph])
        gr_hgg[igraph].SetName("low_" + str(igraph))
        gr_hgg[igraph].SetMarkerStyle(20)
        gr_hgg[igraph].SetMarkerColor( color[igraph -1] )
        gr_hgg[igraph].SetLineColor( color[igraph -1] )
        gr_hgg[igraph].SetTitle("")
        gr_hgg[igraph].GetXaxis().SetLabelSize(20)
        gr_hgg[igraph].GetXaxis().SetLimits(-.2, 4.2)
        a = 21
        b = 5
        gr_hgg[igraph].GetXaxis().SetBinLabel(a * (0 + 0) + b, "ggH")
        gr_hgg[igraph].GetXaxis().SetBinLabel(a * (0 + 1) + b, "qqH")
        gr_hgg[igraph].GetXaxis().SetBinLabel(a * (0 + 2) + b, "WH")
        gr_hgg[igraph].GetXaxis().SetBinLabel(a * (0 + 3) + b, "ZH")
        gr_hgg[igraph].GetXaxis().SetBinLabel(a * (0 + 4) + b, "ttH")
        gr_hgg[igraph].GetYaxis().SetTitle("Signal selection efficiency (%)")
        gr_hgg[igraph].GetYaxis().SetTitleOffset(1.04 / canvasSplit * 100.)
        gr_hgg[igraph].GetYaxis().SetTickLength(0.03 / canvasSplit * 100.) # due to assymmetric canvas...
        gr_hgg[igraph].SetMinimum(0.)
        gr_hgg[igraph].SetMaximum(100.)
        if igraph == 1:
            gr_hgg[igraph].Draw("ALP")
        else: 
            gr_hgg[igraph].Draw("LP")
    c1.cd()
    legend = TLegend(0.50, 0.16, 0.80, 0.36, "")
#    legend.SetNColumns(2)
    legend.SetTextSize(0.03)
#    legend.SetTextFont(63) # precision 3
#    legend.SetTextSize(18) # in pixel
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    legend.AddEntry(gr_low[1].GetName(), "# #gamma #geq 2, # jets #geq 2", "lp")
    legend.AddEntry(gr_low[2].GetName(), "#gamma presel. + p_{T} cuts", "lp")
    legend.AddEntry(gr_low[3].GetName(), "#gamma ID", "lp")
    legend.AddEntry(gr_low[4].GetName(), "jet ID", "lp")
    legend.AddEntry(gr_low[5].GetName(), "at least one bjet", "lp")
    legend.AddEntry(gr_low[6].GetName(), "mass cuts", "lp")
    legend.Draw()
    latexLabel = TLatex()
#    latexLabel.SetTextSize(0.03)
    latexLabel.SetTextFont(63)
    latexLabel.SetTextSize(18) # in pixel
    latexLabel.SetNDC()
    latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary   #sqrt{s} = 8 TeV   L = 19.7 fb^{-1}")
    latexLabel.SetTextAngle(45)
    c1.Print("flow.pdf")
    c1.Print("flow.png")
    c1.Print("flow.root")
    del c1
