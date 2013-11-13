#!/bin/bash

version="v18bis"
today=`date +"0%Y-%m-%d"`
set -x

# do mgg trees
#for data in `echo "diphojet_sherpa_8TeV DYJetsToLL gjet_20_8TeV_pf gjet_40_8TeV_pf qcd_30_8TeV_ff qcd_30_8TeV_pf qcd_40_8TeV_ff qcd_40_8TeV_pf"`
#for data in `echo "ggh_m125_minlo_8TeV tth_m125_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh ggh_m125_powheg_8TeV"`
#for data in `echo "ggh_m125_powheg_8TeV"`
#for data in `echo "signal data"`
for data in `echo "signal"`
do
	for mass in `echo "300"`
	do
		for cutLevel in `echo "0"`
		do
			for fitStrategy in `echo "mgg"`
			do
				for whichJet in `echo "base reg kin regkin"`	
				do

					if [[ "${data}" == "data" ]]
					then
						intree="Data"
						removeUndefinedBtagSF=0
						type_=0
					elif [[ "${data}" == "signal" ]]
					then
						intree="Radion_m${mass}_8TeV_nm"
						removeUndefinedBtagSF=1
						type_=-${mass}
					else
						intree="${data}"
						removeUndefinedBtagSF=0
						type_=0 # harmless as it is used only to check btagSF
					fi
					if [[ "${whichJet}"	== "base" ]]
					then
						infile="2013-10-28_selection_noRegression_noMassCut_v01/${intree}_noRegression_noMassCut_v01.root"
					else
						infile="2013-10-28_selection_PhilRegr1028_noMassCut_v01/${intree}_PhilRegr1028_noMassCut_v01.root"
					fi
					outfolder="${version}_${whichJet}_${fitStrategy}_${cutLevel}"
					outfile="${today}-${intree}_m${mass}.root"
					outtree="TCVARS"

for massCutVersion in `echo "0"`
do
					mkdir -p ${outfolder}_massCutVersion${massCutVersion}

./quickTrees.exe \
-i ${infile} \
-it ${intree} \
-o ${outfolder}_massCutVersion${massCutVersion}/${outfile} \
-cutLevel ${cutLevel} \
-m ${mass} \
-fs ${fitStrategy} \
-wj ${whichJet} \
--removeUndefinedBtagSF ${removeUndefinedBtagSF} \
--type ${type_} \
--massCutVersion ${massCutVersion} \
| tee ${outfolder}_massCutVersion${massCutVersion}/${today}-${intree}_m${mass}.log

done # massCutVersion

				done # whichJet
			done # fitStrategy
		done # cutLevel
	done # mass
done # data# do mgg trees



