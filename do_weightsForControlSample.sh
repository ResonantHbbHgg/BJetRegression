#!/bin/bash

selectionDate="2014-10-21"
selectionVersion="v17"
output="scales_2D_pt_data_4GeVbinning"

basedir=""
#basedir="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/${selectionVersion}"

echo "./obtainWeights.exe in ${output}.root" 
./obtainWeights.exe \
--inputfile_data "${basedir}${selectionDate}_selection_noRegression_noMassCut_${selectionVersion}/Data_noRegression_noMassCut_${selectionVersion}.root" \
--inputfile_CS "${basedir}${selectionDate}_selection_noRegression_noMassCut_${selectionVersion}/Data_noRegression_noMassCut_controlSample_${selectionVersion}.root" \
--inputtree_data "Data" \
--inputtree_CS "Data" \
--outputfile "${output}.root" &> ${output}.eo
