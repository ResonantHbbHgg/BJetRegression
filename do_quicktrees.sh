#!/bin/bash

version="v12"
today=`date +"0%Y-%m-%d"`

# do mgg trees
for data in `echo "data signal"`
do
	for mass in `echo "300 500"`
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
					elif [[ "${data}" == "signal" ]]
					then
						intree="Radion_m${mass}_8TeV_nm"
						removeUndefinedBtagSF=1
					fi
					if [[ "${whichJet}"	== "base" ]]
					then
						infile="2013-10-28_selection_noRegression_noMassCut_v01/${intree}_noRegression_noMassCut_v01.root"
					else
						infile="2013-10-28_selection_PhilRegr1028_noMassCut_v01/${intree}_PhilRegr1028_noMassCut_v01.root"
					fi
					outfolder="${version}_${whichJet}_${fitStrategy}_${cutLevel}"
					mkdir -p ${outfolder}
					outfile="${today}-${intree}_m${mass}.root"
					outtree="TCVARS"


./quickTrees.exe \
-i ${infile} \
-it ${intree} \
-o ${outfolder}/${outfile} \
-cutLevel ${cutLevel} \
-m ${mass} \
-fs ${fitStrategy} \
-wj ${whichJet} \
--removeUndefinedBtagSF ${removeUndefinedBtagSF} \
| tee ${outfolder}/${today}-${intree}_m${mass}.log

				done # whichJet
			done # fitStrategy
		done # cutLevel
	done # mass
done # data# do mgg trees


for data in `echo "data signal"`
do
	for mass in `echo "300 500 700 1000"`
	do
		for cutLevel in `echo "0"`
		do
			for fitStrategy in `echo "mggjj"`
			do
				for whichJet in `echo "base reg kin regkin"`	
				do

					if [[ "${data}" == "data" ]]
					then
						intree="Data"
						removeUndefinedBtagSF=0
					elif [[ "${data}" == "signal" ]]
					then
						intree="Radion_m${mass}_8TeV_nm"
						removeUndefinedBtagSF=1
					fi
					if [[ "${whichJet}"	== "base" ]]
					then
						infile="2013-10-28_selection_noRegression_noMassCut_v01/${intree}_noRegression_noMassCut_v01.root"
					else
						infile="2013-10-28_selection_PhilRegr1028_noMassCut_v01/${intree}_PhilRegr1028_noMassCut_v01.root"
					fi
					outfolder="${version}_${whichJet}_${fitStrategy}_${cutLevel}"
					mkdir -p ${outfolder}
					outfile="${today}-${intree}_m${mass}.root"
					outtree="TCVARS"


./quickTrees.exe \
-i ${infile} \
-it ${intree} \
-o ${outfolder}/${outfile} \
-cutLevel ${cutLevel} \
-m ${mass} \
-fs ${fitStrategy} \
-wj ${whichJet} \
--removeUndefinedBtagSF ${removeUndefinedBtagSF} \
| tee ${outfolder}/${today}-${intree}_m${mass}.log

				done # whichJet
			done # fitStrategy
		done # cutLevel
	done # mass
done # data
