#!/bin/bash

./fitMass.exe --inputtree Radion_m500_8TeV_nm --type -500 --inputfile 2013-10-28_selection_noRegression_noMassCut_v01/Radion_m500_8TeV_nm_noRegression_noMassCut_v01.root --inputfile_reg 2013-10-28_selection_PhilRegr1028_noMassCut_v01/Radion_m500_8TeV_nm_PhilRegr1028_noMassCut_v01.root
./fitMass.exe --inputtree Radion_m700_8TeV_nm --type -700 --inputfile 2013-10-28_selection_noRegression_noMassCut_v01/Radion_m700_8TeV_nm_noRegression_noMassCut_v01.root --inputfile_reg 2013-10-28_selection_PhilRegr1028_noMassCut_v01/Radion_m700_8TeV_nm_PhilRegr1028_noMassCut_v01.root
./fitMass.exe --inputtree Radion_m1000_8TeV_nm --type -1000 --inputfile 2013-10-28_selection_noRegression_noMassCut_v01/Radion_m1000_8TeV_nm_noRegression_noMassCut_v01.root --inputfile_reg 2013-10-28_selection_PhilRegr1028_noMassCut_v01/Radion_m1000_8TeV_nm_PhilRegr1028_noMassCut_v01.root
