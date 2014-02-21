#!/bin/bash

inputversion="v10"
inputfolder="2014-02-17_selection_noRegression_noMassCut_${inputversion}/"

i=-1

for mass in `echo "260 270 300 350 400 450 500 550 600 650 700 800 900 1000 1100"`
do
	inputfile="${inputfolder}/Radion_m${mass}_8TeV_noRegression_noMassCut_${inputversion}.root"
	inputtree="Radion_m${mass}_8TeV"
	if [[ "${mass}" == "260" ]]
	then
		inputfile="${inputfile/Radion/MSSM}"
		inputtree="${inputtree/Radion/MSSM}"
	fi
	./fitMass.exe \
--inputfile ${inputfile} \
--inputtree ${inputtree} \
--massHypothesis ${mass} \
--outputfile resolution_SigmaEff_m${mass}.txt
done

cat resolution_SigmaEff_m*.txt > resolution_SigmaEff_all.txt
rm resolution_SigmaEff_m*.txt 
