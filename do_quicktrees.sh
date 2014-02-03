#!/bin/bash

version="v25"
today=`date +"0%Y-%m-%d"`
#set -x

i=-1


##### PREPARE MGGJJ-FIT TREES
kinfitlabel[0]="noKinFit"
kinfitlabel[1]="withKinFit"
kinfitjet[0]="base"
kinfitjet[1]="kin"
for ikin in `seq 0 1`
do
	outfolder="${version}_fitToMggjj_${kinfitlabel[${ikin}]}"
	mkdir -p ${outfolder}
	for sample in `echo "Radion Graviton Data"`
	do
		for mass in `echo "400 450 500 550 600 650 700 800 900 1000 1100"`
		do
			intree=${sample}
			itype="1"
			removeUndefinedBtagSF=0
			if [ "${sample}" == "Radion" ]
			then
				intree="${sample}_m${mass}_8TeV"
				itype="-${mass}"
				removeUndefinedBtagSF=0
			elif [ "${sample}" == "Graviton" ]
			then
				if [ "${mass}" == "500" ] || [ "${mass}" == "700" ] || [ "${mass}" == "1000" ]
				then
					intree="${sample}_m${mass}_8TeV"
					itype="-${mass}"
					removeUndefinedBtagSF=0
				else
					continue
				fi
			elif [ "${sample}" == "Data" ]
			then
				itype="0"
			fi
			i=$((${i} + 1))
			line[${i}]=""
			line[${i}]="${line[${i}]} --inputfile 2014-01-31_selection_noRegression_noMassCut_v06/${intree}_noRegression_noMassCut_v06.root"
			line[${i}]="${line[${i}]} --inputtree ${intree}"
			line[${i}]="${line[${i}]} --outputtree TCVARS"
			line[${i}]="${line[${i}]} --outputfile ${outfolder}/${intree}_m${mass}.root"
			line[${i}]="${line[${i}]} --type ${itype}"
			line[${i}]="${line[${i}]} --whichJet ${kinfitjet[${ikin}]}"
			line[${i}]="${line[${i}]} --fitStrategy mggjj"
			line[${i}]="${line[${i}]} --cutLevel 0"
			line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
			line[${i}]="${line[${i}]} --massCutVersion 3"
			log[${i}]="${outfolder}/${intree}_m${mass}.eo"
		#	echo -e "i= ${i}\tline= ${line[${i}]}"
		done
	done
done
	
##### PREPARE MGG-FIT TREES
outfolder="${version}_fitToMgg_noKinFit"
mkdir -p ${outfolder}
for sample in `echo "Radion Graviton MSSM ggh_m125_powheg_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh tth_m125_8TeV Data"`
do
	for mass in `echo "260 270 300 350 400 450 500"`
	do
		intree=${sample}
		itype="1"
		removeUndefinedBtagSF=0
		if [ "${sample}" == "Radion" ] 
		then
			if [ "${mass}" == "260" ]
			then
				continue
			else
				intree="${sample}_m${mass}_8TeV"
				itype="-${mass}"
				removeUndefinedBtagSF=0
			fi
		elif [ "${sample}" == "Graviton" ]
		then
			if [ "${mass}" == "260" ] || [ "${mass}" == "270" ] || [ "${mass}" == "350" ] || [ "${mass}" == "400" ] || [ "${mass}" == "450" ]
			then
				continue
			else
				intree="${sample}_m${mass}_8TeV"
				itype="-${mass}"
				removeUndefinedBtagSF=0
			fi
		elif [ "${sample}" == "MSSM" ]
		then
			if [ "${mass}" == "270" ] || [ "${mass}" == "400" ] || [ "${mass}" == "450" ] || [ "${mass}" == "500" ]
			then
				continue
			else
				intree="${sample}_m${mass}_8TeV"
				itype="-${mass}"
				removeUndefinedBtagSF=0
			fi
		elif [ "${sample}" == "Data" ]
		then
			itype="0"
		fi
		i=$((${i} + 1))
		line[${i}]=""
		line[${i}]="${line[${i}]} --inputfile 2014-01-31_selection_noRegression_noMassCut_v06/${intree}_noRegression_noMassCut_v06.root"
		line[${i}]="${line[${i}]} --inputtree ${intree}"
		line[${i}]="${line[${i}]} --outputtree TCVARS"
		line[${i}]="${line[${i}]} --outputfile ${outfolder}/${intree}_m${mass}.root"
		line[${i}]="${line[${i}]} --type ${itype}"
		line[${i}]="${line[${i}]} --whichJet base"
		line[${i}]="${line[${i}]} --fitStrategy mgg"
		line[${i}]="${line[${i}]} --cutLevel 0"
		line[${i}]="${line[${i}]} --mass ${mass}"
		line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
		line[${i}]="${line[${i}]} --massCutVersion 3"
		log[${i}]="${outfolder}/${intree}_m${mass}.eo"
	#	echo -e "i= ${i}\tline= ${line[${i}]}"
	done
done

#### PRODUCE EVERYTHING	
itot=${i}
for iline in `seq 0 ${itot}`
do
	echo -e "iline= ${iline} / ${itot}\t\t${line[${iline}]}"
	./quickTrees.exe ${line[${iline}]} &> ${log[${iline}]}
done
