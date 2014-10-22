#!/usr/bin/env python
from array import array
from os import path
import ROOT
from ROOT import gROOT
from ROOT import TFile
gROOT.SetBatch()

plottreedir = "2014-10-22_selection_withRegression_noMassCut_v18/"
version = "v18"
masses = [270, 300, 350, 400]

for mass in masses:
    infile = TFile.Open( path.join(plottreedir, "Radion_m" + str(mass) + "_8TeV_withRegression_noMassCut_" + version + ".root") )
    intree = infile.Get("Radion_m" + str(mass) + "_8TeV")
    outfile = TFile( path.join(plottreedir, "Radion_m" + str(mass) + "_8TeV_forcedWithRegression_noMassCut_" + version + ".root"), "RECREATE" )
    
    intree.SetBranchStatus("*", 1)
    intree.SetBranchStatus("jj_mass", 0)
    intree.SetBranchStatus("ggjj_mass", 0)
    outtree = intree.CloneTree(0)
    jj_mass = array('f', [1.])
    ggjj_mass = array('f', [1.])
    outtree.Branch("jj_mass", jj_mass, "jj_mass/F")
    outtree.Branch("ggjj_mass", ggjj_mass, "ggjj_mass/F")
    
    for ievt in xrange(intree.GetEntries()):
        intree.GetEntry( ievt )
        jj_mass[0] = intree.regjj_mass
        ggjj_mass[0] = intree.regggjj_mass
        outtree.Fill()
    
    outfile.cd()
    outtree.Write()
    outfile.Close()
    infile.Close()

