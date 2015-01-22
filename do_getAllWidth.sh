#!/bin/bash

inputversion="v22"
inputfolder="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/${inputversion}/2014-12-16_selection_noRegression_noMassCut_${inputversion}"
inputfolderreg="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/${inputversion}/2014-12-16_selection_withRegression_noMassCut_${inputversion}"

i=-1

compareWithRegression=1
for mass in `echo "270 300 350 400"`
do
	inputfile="${inputfolder}/Radion_m${mass}_8TeV_noRegression_noMassCut_${inputversion}.root"
	inputtree="Radion_m${mass}_8TeV"
	inputregfile="${inputfolderreg}/Radion_m${mass}_8TeV_withRegression_noMassCut_${inputversion}.root"
	inputregtree="Radion_m${mass}_8TeV"
	if [[ "${mass}" == "260" ]]
	then
		inputfile="${inputfile/Radion/MSSM}"
		inputtree="${inputtree/Radion/MSSM}"
	fi
	./fitMass.exe \
--inputfile ${inputfile} \
--inputtree ${inputtree} \
--inputregfile ${inputregfile} \
--inputregtree ${inputregtree} \
--massHypothesis ${mass} \
--compareWithRegression ${compareWithRegression} \
--outputfile resolution_SigmaEff_m${mass}.txt
done

cat resolution_SigmaEff_m*.txt > resolution_SigmaEff_all.txt
rm resolution_SigmaEff_m*.txt 
