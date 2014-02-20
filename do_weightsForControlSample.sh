#!/bin/bash

selectionDate="2014-02-17"
selectionVersion="v10"
output="scales_2D_pt_data_4GeVbinning"

#./selection.exe --inputfile root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_tree_v08/Data_Full2012.root --inputtree Data --outputtree Data --outputfile ${selectionDate}_selection_noRegression_noMassCut_${selectionVersion}/Data_noRegression_noMassCut_controlSample_${selectionVersion}.root --numberOfRegressionFiles 0 --type 0 --removeUndefinedBtagSF 0 --applyMassCuts 0 --applyPhotonIDControlSample 1 &> ${selectionDate}_selection_noRegression_noMassCut_${selectionVersion}/Data_noRegression_noMassCut_controlSample_${selectionVersion}.eo

echo "./obtainWeights.exe in ${output}.root" 
./obtainWeights.exe \
--inputfile_data "${selectionDate}_selection_noRegression_noMassCut_${selectionVersion}/Data_noRegression_noMassCut_${selectionVersion}.root" \
--inputfile_CS "${selectionDate}_selection_noRegression_noMassCut_${selectionVersion}/Data_noRegression_noMassCut_controlSample_${selectionVersion}.root" \
--inputtree_data "Data" \
--inputtree_CS "Data" \
--outputfile "${output}.root" &> ${output}.eo
