#!/bin/bash

today=`date +%Y-%m-%d`

inputfile="root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_tree_v06/Radion_Graviton_nm.root"
storingdir="${today}-jetTreeForTraining"
mkdir -p ${storingdir}
for mass in `echo "300 500 700 1000 1500"`
do
	echo "mass= ${mass}"
	inputtree="Radion_m${mass}_8TeV_nm"
	outputfile="${storingdir}/jetTreeForTraining_m${mass}.root"
	./prepareOpTreeInJetTree_forTraining.exe -i ${inputfile} -it ${inputtree} -o ${outputfile} &> ${storingdir}/jetTreeForTraining_m${mass}.log
done


