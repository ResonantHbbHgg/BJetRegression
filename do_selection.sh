#!/bin/bash

today=`date +"%Y-%m-%d"`
version="v10_0btag"

eosprefix="root://eoscms//eos/cms"
eospath="/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/"

i=-1

### ##### DATA
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08/Data_Full2012.root"
### tree[${i}]="Data"
### outtree[${i}]="Data"
### typ[${i}]="0"
### 
### ##### SIGNAL
### ### LONG SAMPLE LIST
i=$((${i} + 1))
infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b_0btag.root"
tree[${i}]="Radion_m270_8TeV"
outtree[${i}]="Radion_m270_8TeV"
typ[${i}]="-270"

i=$((${i} + 1))
infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b_0btag.root"
tree[${i}]="Radion_m400_8TeV"
outtree[${i}]="Radion_m400_8TeV"
typ[${i}]="-400"

i=$((${i} + 1))
infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b_0btag.root"
tree[${i}]="Radion_m1100_8TeV_nm"
outtree[${i}]="Radion_m1100_8TeV"
typ[${i}]="-1100"

### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m270_8TeV"
### outtree[${i}]="Radion_m270_8TeV"
### typ[${i}]="-270"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Radion_nm.root"
### tree[${i}]="Radion_m300_8TeV_nm"
### outtree[${i}]="Radion_m300_8TeV"
### typ[${i}]="-300"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m350_8TeV"
### outtree[${i}]="Radion_m350_8TeV"
### typ[${i}]="-350"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m400_8TeV"
### outtree[${i}]="Radion_m400_8TeV"
### typ[${i}]="-400"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m450_8TeV"
### outtree[${i}]="Radion_m450_8TeV"
### typ[${i}]="-450"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Radion_nm.root"
### tree[${i}]="Radion_m500_8TeV_nm"
### outtree[${i}]="Radion_m500_8TeV"
### typ[${i}]="-500"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m550_8TeV"
### outtree[${i}]="Radion_m550_8TeV"
### typ[${i}]="-550"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m600_8TeV"
### outtree[${i}]="Radion_m600_8TeV"
### typ[${i}]="-600"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m650_8TeV"
### outtree[${i}]="Radion_m650_8TeV"
### typ[${i}]="-650"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Radion_nm.root"
### tree[${i}]="Radion_m700_8TeV_nm"
### outtree[${i}]="Radion_m700_8TeV"
### typ[${i}]="-700"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m800_8TeV"
### outtree[${i}]="Radion_m800_8TeV"
### typ[${i}]="-800"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m900_8TeV"
### outtree[${i}]="Radion_m900_8TeV"
### typ[${i}]="-900"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Radion_nm.root"
### tree[${i}]="Radion_m1000_8TeV_nm"
### outtree[${i}]="Radion_m1000_8TeV"
### typ[${i}]="-1000"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m1100_8TeV"
### outtree[${i}]="Radion_m1100_8TeV"
### typ[${i}]="-1100"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m1200_8TeV"
### outtree[${i}]="Radion_m1200_8TeV"
### typ[${i}]="-1200"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m1300_8TeV"
### outtree[${i}]="Radion_m1300_8TeV"
### typ[${i}]="-1300"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/RadionToHH_2Gamma_2b.root"
### tree[${i}]="Radion_m1400_8TeV"
### outtree[${i}]="Radion_m1400_8TeV"
### typ[${i}]="-1400"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Radion_nm.root"
### tree[${i}]="Radion_m1500_8TeV_nm"
### outtree[${i}]="Radion_m1500_8TeV"
### typ[${i}]="-1500"
### 
### ### MSSM samples
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/MSSM_Higgs_RD.root"
### tree[${i}]="MSSM_H_m260_8TeV"
### outtree[${i}]="MSSM_m260_8TeV"
### typ[${i}]="-260"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/MSSM_Higgs_RD.root"
### tree[${i}]="MSSM_H_m300_8TeV"
### outtree[${i}]="MSSM_m300_8TeV"
### typ[${i}]="-300"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/MSSM_Higgs_RD.root"
### tree[${i}]="MSSM_H_m350_8TeV"
### outtree[${i}]="MSSM_m350_8TeV"
### typ[${i}]="-350"
### 
### ###### GRAVITON
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Graviton.root"
### tree[${i}]="Graviton_m300_8TeV"
### outtree[${i}]="Graviton_m300_8TeV"
### typ[${i}]="-300"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Graviton.root"
### tree[${i}]="Graviton_m500_8TeV"
### outtree[${i}]="Graviton_m500_8TeV"
### typ[${i}]="-500"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Graviton.root"
### tree[${i}]="Graviton_m700_8TeV"
### outtree[${i}]="Graviton_m700_8TeV"
### typ[${i}]="-700"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Graviton.root"
### tree[${i}]="Graviton_m1000_8TeV"
### outtree[${i}]="Graviton_m1000_8TeV"
### typ[${i}]="-1000"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Graviton.root"
### tree[${i}]="Graviton_m1500_8TeV"
### outtree[${i}]="Graviton_m1500_8TeV"
### typ[${i}]="-1500"
### 
### ###### SM Higgs
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/SMHiggs.root"
### tree[${i}]="ggh_m125_minlo_8TeV"
### outtree[${i}]="ggh_m125_minlo_8TeV"
### typ[${i}]="-1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/SMHiggs.root"
### tree[${i}]="ggh_m125_powheg_8TeV"
### outtree[${i}]="ggh_m125_powheg_8TeV"
### typ[${i}]="-1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/SMHiggs.root"
### tree[${i}]="vbf_m125_8TeV"
### outtree[${i}]="vbf_m125_8TeV"
### typ[${i}]="-1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/SMHiggs.root"
### tree[${i}]="wzh_m125_8TeV_wh"
### outtree[${i}]="wzh_m125_8TeV_wh"
### typ[${i}]="-1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/SMHiggs.root"
### tree[${i}]="wzh_m125_8TeV_zh"
### outtree[${i}]="wzh_m125_8TeV_zh"
### typ[${i}]="-1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/SMHiggs.root"
### tree[${i}]="tth_m125_8TeV"
### outtree[${i}]="tth_m125_8TeV"
### typ[${i}]="-1"
### 
### ###### Diphoton backgrounds
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="qcd_30_8TeV_ff"
### outtree[${i}]="qcd_30_8TeV_ff"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="qcd_40_8TeV_ff"
### outtree[${i}]="qcd_40_8TeV_ff"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="qcd_30_8TeV_pf"
### outtree[${i}]="qcd_30_8TeV_pf"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="qcd_40_8TeV_pf"
### outtree[${i}]="qcd_40_8TeV_pf"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="gjet_20_8TeV_pf"
### outtree[${i}]="gjet_20_8TeV_pf"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="gjet_40_8TeV_pf"
### outtree[${i}]="gjet_40_8TeV_pf"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="diphojet_sherpa_8TeV"
### outtree[${i}]="diphojet_sherpa_8TeV"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="DYJetsToLL"
### outtree[${i}]="DYJetsToLL"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="LNuGG_FSR_8TeV"
### outtree[${i}]="LNuGG_FSR_8TeV"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="LNuGG_ISR_8TeV"
### outtree[${i}]="LNuGG_ISR_8TeV"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="ttGG_8TeV"
### outtree[${i}]="ttGG_8TeV"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="tGG_8TeV"
### outtree[${i}]="tGG_8TeV"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="TTGJets_8TeV"
### outtree[${i}]="TTGJets_8TeV"
### typ[${i}]="1"
### 
### i=$((${i} + 1))
### infile[${i}]="radion_tree_v08c/Backgrounds.root"
### tree[${i}]="ZGToLLG_8TeV"
### outtree[${i}]="ZGToLLG_8TeV"
### typ[${i}]="1"
### 

itot=${i}
masscutsuffix[0]="noMassCut"
regsuffix[0]="noRegression"
#regsuffix[2]="PhilRegr1017"
#regsuffix[2]="PhilRegr1028"
imasscut=0

for ireg in `echo "0"`
do
	folder="${today}_selection_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_${version}"
	mkdir -p ${folder}
	for isample in `seq 0 ${itot}`
	do
		echo -e "isample= ${isample} / ${itot}\tinfile= ${infile[${isample}]}\ttree= ${tree[${isample}]}\touttree= ${outtree[${isample}]}\ttyp= ${typ[${isample}]}\tireg= ${ireg}"
		file="${outtree[${isample}]}_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_${version}"
		removeZeroFlavour=0
		if [ "${typ[${isample}]}" -lt "-250" ]
		then 
			removeZeroFlavour=1
		fi

./selection.exe \
--inputfile ${eosprefix}${eospath}${infile[${isample}]} \
--inputtree ${tree[${isample}]} \
--outputtree ${outtree[${isample}]} \
--outputfile ${folder}/${file}.root \
--numberOfRegressionFiles 0 \
--type ${typ[${isample}]} \
--removeUndefinedBtagSF ${removeZeroFlavour} \
--applyMassCuts ${imasscut} \
&> ${folder}/${file}.eo

	done # isample
done #ireg



