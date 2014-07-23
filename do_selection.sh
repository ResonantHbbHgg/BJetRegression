#!/bin/bash

today=`date +"%Y-%m-%d"`
version="v15_photonID_wo_IsoAB"

eosprefix="root://eoscms//eos/cms"
#eospath="/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/"
eospath="/store/cmst3/group/hbbhgg/H2GGLOBE/Radion/trees/"

keep0btag=1

i=-1

######### DATA
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/Data.root"
####tree[${i}]="Data"
####outtree[${i}]="Data"
####typ[${i}]="0"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/Data.root"
####tree[${i}]="Data"
####outtree[${i}]="Data"
####typ[${i}]="0"
####CS[${i}]="1"
####
##### SIGNAL
### LONG SAMPLE LIST
i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m270_8TeV"
outtree[${i}]="Radion_m270_8TeV"
typ[${i}]="-270"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m300_8TeV_nm"
outtree[${i}]="Radion_m300_8TeV"
typ[${i}]="-300"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m350_8TeV"
outtree[${i}]="Radion_m350_8TeV"
typ[${i}]="-350"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m400_8TeV"
outtree[${i}]="Radion_m400_8TeV"
typ[${i}]="-400"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m450_8TeV"
outtree[${i}]="Radion_m450_8TeV"
typ[${i}]="-450"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m500_8TeV_nm"
outtree[${i}]="Radion_m500_8TeV"
typ[${i}]="-500"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m550_8TeV"
outtree[${i}]="Radion_m550_8TeV"
typ[${i}]="-550"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m600_8TeV"
outtree[${i}]="Radion_m600_8TeV"
typ[${i}]="-600"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m650_8TeV"
outtree[${i}]="Radion_m650_8TeV"
typ[${i}]="-650"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m700_8TeV_nm"
outtree[${i}]="Radion_m700_8TeV"
typ[${i}]="-700"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m800_8TeV"
outtree[${i}]="Radion_m800_8TeV"
typ[${i}]="-800"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m900_8TeV"
outtree[${i}]="Radion_m900_8TeV"
typ[${i}]="-900"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m1000_8TeV_nm"
outtree[${i}]="Radion_m1000_8TeV"
typ[${i}]="-1000"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="Radion_m1100_8TeV"
outtree[${i}]="Radion_m1100_8TeV"
typ[${i}]="-1100"
CS[${i}]="0"

####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Radion_m1200_8TeV"
####outtree[${i}]="Radion_m1200_8TeV"
####typ[${i}]="-1200"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Radion_m1300_8TeV"
####outtree[${i}]="Radion_m1300_8TeV"
####typ[${i}]="-1300"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Radion_m1400_8TeV"
####outtree[${i}]="Radion_m1400_8TeV"
####typ[${i}]="-1400"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Radion_m1500_8TeV_nm"
####outtree[${i}]="Radion_m1500_8TeV"
####typ[${i}]="-1500"
####CS[${i}]="0"
####
### MSSM samples
i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="MSSM_H_m260_8TeV"
outtree[${i}]="MSSM_m260_8TeV"
typ[${i}]="-260"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="MSSM_H_m300_8TeV"
outtree[${i}]="MSSM_m300_8TeV"
typ[${i}]="-300"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/XHH.root"
tree[${i}]="MSSM_H_m350_8TeV"
outtree[${i}]="MSSM_m350_8TeV"
typ[${i}]="-350"
CS[${i}]="0"

######### GRAVITON
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Graviton_m300_8TeV"
####outtree[${i}]="Graviton_m300_8TeV"
####typ[${i}]="-300"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Graviton_m500_8TeV"
####outtree[${i}]="Graviton_m500_8TeV"
####typ[${i}]="-500"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Graviton_m700_8TeV"
####outtree[${i}]="Graviton_m700_8TeV"
####typ[${i}]="-700"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Graviton_m1000_8TeV"
####outtree[${i}]="Graviton_m1000_8TeV"
####typ[${i}]="-1000"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/XHH.root"
####tree[${i}]="Graviton_m1500_8TeV"
####outtree[${i}]="Graviton_m1500_8TeV"
####typ[${i}]="-1500"
####CS[${i}]="0"
####
##### SM Higgs
i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="ggh_m125_minlo_8TeV"
outtree[${i}]="ggh_m125_minlo_8TeV"
typ[${i}]="-1"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="ggh_m125_powheg_8TeV"
outtree[${i}]="ggh_m125_powheg_8TeV"
typ[${i}]="-1"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="vbf_m125_8TeV"
outtree[${i}]="vbf_m125_8TeV"
typ[${i}]="-1"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="wzh_m125_8TeV_wh"
outtree[${i}]="wzh_m125_8TeV_wh"
typ[${i}]="-1"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="wzh_m125_8TeV_zh"
outtree[${i}]="wzh_m125_8TeV_zh"
typ[${i}]="-1"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="tth_m125_8TeV"
outtree[${i}]="tth_m125_8TeV"
typ[${i}]="-1"
CS[${i}]="0"

i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="bbh_m125_8TeV"
outtree[${i}]="bbh_m125_8TeV"
typ[${i}]="-1"
CS[${i}]="0"


##### SM di-Higgs
i=$((${i} + 1))
infile[${i}]="radion_tree_v09/SMHiggs.root"
tree[${i}]="ggHH_8TeV"
outtree[${i}]="ggHH_8TeV"
typ[${i}]="-2"
CS[${i}]="0"


######### Diphoton backgrounds
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/diphoton.root"
####tree[${i}]="qcd_30_8TeV_ff"
####outtree[${i}]="qcd_30_8TeV_ff"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/diphoton.root"
####tree[${i}]="qcd_40_8TeV_ff"
####outtree[${i}]="qcd_40_8TeV_ff"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/diphoton.root"
####tree[${i}]="qcd_30_8TeV_pf"
####outtree[${i}]="qcd_30_8TeV_pf"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/diphoton.root"
####tree[${i}]="qcd_40_8TeV_pf"
####outtree[${i}]="qcd_40_8TeV_pf"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/diphoton.root"
####tree[${i}]="gjet_20_8TeV_pf"
####outtree[${i}]="gjet_20_8TeV_pf"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/diphoton.root"
####tree[${i}]="gjet_40_8TeV_pf"
####outtree[${i}]="gjet_40_8TeV_pf"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/diphoton.root"
####tree[${i}]="diphojet_sherpa_8TeV"
####outtree[${i}]="diphojet_sherpa_8TeV"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/backgrounds.root"
####tree[${i}]="DYJetsToLL"
####outtree[${i}]="DYJetsToLL"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/backgrounds.root"
####tree[${i}]="LNuGG_FSR_8TeV"
####outtree[${i}]="LNuGG_FSR_8TeV"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/backgrounds.root"
####tree[${i}]="LNuGG_ISR_8TeV"
####outtree[${i}]="LNuGG_ISR_8TeV"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/backgrounds.root"
####tree[${i}]="ttGG_8TeV"
####outtree[${i}]="ttGG_8TeV"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/backgrounds.root"
####tree[${i}]="tGG_8TeV"
####outtree[${i}]="tGG_8TeV"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/backgrounds.root"
####tree[${i}]="TTGJets_8TeV"
####outtree[${i}]="TTGJets_8TeV"
####typ[${i}]="1"
####CS[${i}]="0"
####
####i=$((${i} + 1))
####infile[${i}]="radion_tree_v09/backgrounds.root"
####tree[${i}]="ZGToLLG_8TeV"
####outtree[${i}]="ZGToLLG_8TeV"
####typ[${i}]="1"
####CS[${i}]="0"
####

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
        applyCS="${CS[${isample}]}"
        echo -e "isample= ${isample} / ${itot}\tinfile= ${infile[${isample}]}\ttree= ${tree[${isample}]}\touttree= ${outtree[${isample}]}\ttyp= ${typ[${isample}]}\tireg= ${ireg}\tapplyCS= ${applyCS}"
        file="${outtree[${isample}]}_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_${version}"
        if [[ "${applyCS}" == "1" ]]
        then
            file="${outtree[${isample}]}_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_controlSample_${version}"
        fi
        removeZeroFlavour=0
        if [ "${typ[${isample}]}" != "0" ]
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
--applyPhotonIDControlSample ${applyCS} \
--printCutFlow 1 \
--keep0btag ${keep0btag} \
2> ${folder}/${file}.eo | tee ${folder}/${file}.eo | egrep '%|entries'

    mv cutFlow_${tree[${isample}]}.dat ${folder}/

    done # isample
done #ireg



