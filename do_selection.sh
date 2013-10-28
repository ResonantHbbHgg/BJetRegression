#!/bin/bash

today=`date +"%Y-%m-%d"`
version="v01"

eosprefix="root://eoscms//eos/cms"
eospath="/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/"

i=-1

##### DATA
i=$((${i} + 1))
infile[${i}]="radion_tree_v06/Data_Full2012.root"
tree[${i}]="Data"
typ[${i}]="0"

##### SIGNAL
i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Radioon-Graviton_nm.root"
tree[${i}]="Radion_m300_8TeV_nm"
typ[${i}]="-300"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Radioon-Graviton_nm.root"
tree[${i}]="Radion_m500_8TeV_nm"
typ[${i}]="-500"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Radioon-Graviton_nm.root"
tree[${i}]="Radion_m700_8TeV_nm"
typ[${i}]="-700"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Radioon-Graviton_nm.root"
tree[${i}]="Radion_m1000_8TeV_nm"
typ[${i}]="-1000"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Radioon-Graviton_nm.root"
tree[${i}]="Radion_m1500_8TeV_nm"
typ[${i}]="-1500"

##### SM Higgs
i=$((${i} + 1))
infile[${i}]="radion_tree_v07/SMHiggs.root"
tree[${i}]="ggh_m125_minlo_8TeV"
typ[${i}]="-1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/SMHiggs.root"
tree[${i}]="ggh_m125_powheg_8TeV"
typ[${i}]="-1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/SMHiggs.root"
tree[${i}]="vbf_m125_8TeV"
typ[${i}]="-1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/SMHiggs.root"
tree[${i}]="wzh_m125_8TeV_wh"
typ[${i}]="-1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/SMHiggs.root"
tree[${i}]="wzh_m125_8TeV_zh"
typ[${i}]="-1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/SMHiggs.root"
tree[${i}]="tth_m125_8TeV"
typ[${i}]="-1"

##### Diphoton backgrounds
i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="qcd_30_8TeV_ff"
typ[${i}]="1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="qcd_40_8TeV_ff"
typ[${i}]="1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="qcd_30_8TeV_pf"
typ[${i}]="1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="qcd_40_8TeV_pf"
typ[${i}]="1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="gjet_20_8TeV_pf"
typ[${i}]="1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="gjet_40_8TeV_pf"
typ[${i}]="1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="diphojet_sherpa_8TeV"
typ[${i}]="1"

i=$((${i} + 1))
infile[${i}]="radion_tree_v07/Backgrounds.root"
tree[${i}]="DYJetsToLL"
typ[${i}]="1"


itot=${i}
masscutsuffix[0]="noMassCut"
regsuffix[0]="noRegression"
#regsuffix[2]="PhilRegr1017"
regsuffix[2]="PhilRegr1028"
imasscut=0

for ireg in `echo "2 0"`
do
	folder="${today}_selection_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_${version}"
	mkdir -p ${folder}
	for isample in `seq 0 ${itot}`
	do
		echo -e "isample= ${isample}\tinfile= ${infile[${isample}]}\ttree= ${tree[${isample}]}\ttyp= ${typ[${isample}]}\tireg= ${ireg}"
		file="${tree[${isample}]}_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_${version}"

./selection.exe \
-i ${eosprefix}${eospath}${infile[${isample}]} \
-t ${tree[${isample}]} \
-o ${folder}/${file}.root \
-r ${ireg} \
--type ${typ[${isample}]} \
--applyMassCuts ${imasscut} \
&> ${folder}/${file}.eo

	done # isample
done #ireg



