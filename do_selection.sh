#!/bin/bash

today=`date +"%Y-%m-%d"`
#today="2014-10-22"
version="v20"
eosprefix="root://eoscms//eos/cms"
eospath="/store/cmst3/group/hbbhgg/H2GGLOBE/Radion/trees/"
#eospath="/store/user/hebda/h2gglobe/trees/"

keep0btag=0
whichPhotonID=0 # 0=CiC Super Tight ;; 1=Francois' ;; CiC 2= photonID MVA
numRegFiles=0 #0 for no regression, 1 for regression; if regression is applied, ensure the latest production is used.
regFilePath="weights/TMVARegression_resonant_BDTG.weights.xml" #it is also possible (but not recommended) to use "weights/TMVARegression_SM_BDTG.weights.xml"

## WHAT TO PROCESS
# DATA
doData=0
doDataCS=0
# RESONANT SIGNALS
doRadion=0
doRadionM126=0
doRadionaajj=0
doMSSM=0
doGraviton=0
doGravitonMore=0
# SM Higgs
doSMHiggs=0
doExtraSMHiggs=0
doSMdiHiggs=0
# BACKGROUNDS
doDiphotonBackgrounds=0
doRareBackgrounds=0
# NON-RESONANT SIGNALS
doAnomalousHH=0
doAnomalousHH_C2=1

# going FTR14001_style, default should be 0
FTR14001_style=0
# Keep both the following to 0 if you do not know how to play with this
nJackknife=0
iJackknife=0
if [ ${nJackknife} != 0 ]
then
    version="${version}_jackknifed_${iJackknife}_over_${nJackknife}"
fi

# Initializing the sample list
i=-1

##### DATA
if [ ${doData} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Data_v2.root"
    tree[${i}]="Data"
    outtree[${i}]="Data"
    typ[${i}]="0"
    CS[${i}]="0"
fi

if [ ${doDataCS} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Data_v2.root"
    tree[${i}]="Data"
    outtree[${i}]="Data"
    typ[${i}]="0"
    CS[${i}]="1"
fi

##### SIGNAL
### LONG SAMPLE LIST
if [ ${doRadion} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m270_8TeV"
    outtree[${i}]="Radion_m270_8TeV"
    typ[${i}]="-270"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m300_8TeV"
    outtree[${i}]="Radion_m300_8TeV"
    typ[${i}]="-300"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m350_8TeV"
    outtree[${i}]="Radion_m350_8TeV"
    typ[${i}]="-350"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m400_8TeV"
    outtree[${i}]="Radion_m400_8TeV"
    typ[${i}]="-400"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m450_8TeV"
    outtree[${i}]="Radion_m450_8TeV"
    typ[${i}]="-450"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m500_8TeV"
    outtree[${i}]="Radion_m500_8TeV"
    typ[${i}]="-500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m550_8TeV"
    outtree[${i}]="Radion_m550_8TeV"
    typ[${i}]="-550"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m600_8TeV"
    outtree[${i}]="Radion_m600_8TeV"
    typ[${i}]="-600"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m650_8TeV"
    outtree[${i}]="Radion_m650_8TeV"
    typ[${i}]="-650"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m700_8TeV"
    outtree[${i}]="Radion_m700_8TeV"
    typ[${i}]="-700"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m800_8TeV"
    outtree[${i}]="Radion_m800_8TeV"
    typ[${i}]="-800"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m900_8TeV"
    outtree[${i}]="Radion_m900_8TeV"
    typ[${i}]="-900"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m1000_8TeV"
    outtree[${i}]="Radion_m1000_8TeV"
    typ[${i}]="-1000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m1100_8TeV"
    outtree[${i}]="Radion_m1100_8TeV"
    typ[${i}]="-1100"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m1200_8TeV"
    outtree[${i}]="Radion_m1200_8TeV"
    typ[${i}]="-1200"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m1300_8TeV"
    outtree[${i}]="Radion_m1300_8TeV"
    typ[${i}]="-1300"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m1400_8TeV"
    outtree[${i}]="Radion_m1400_8TeV"
    typ[${i}]="-1400"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Radion_m1500_8TeV"
    outtree[${i}]="Radion_m1500_8TeV"
    typ[${i}]="-1500"
    CS[${i}]="0"
fi

### Radion mH126 samples
if [ ${doRadionM126} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m270_mH126_8TeV"
    outtree[${i}]="Radion_m270_mH126_8TeV"
    typ[${i}]="-270"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m300_mH126_8TeV"
    outtree[${i}]="Radion_m300_mH126_8TeV"
    typ[${i}]="-300"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m500_mH126_8TeV"
    outtree[${i}]="Radion_m500_mH126_8TeV"
    typ[${i}]="-500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m700_mH126_8TeV"
    outtree[${i}]="Radion_m700_mH126_8TeV"
    typ[${i}]="-700"
    CS[${i}]="0"
fi

### Radion aajj samples
if [ ${doRadionaajj} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m270_aajj_8TeV"
    outtree[${i}]="Radion_m270_aajj_8TeV"
    typ[${i}]="-270"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m300_aajj_8TeV"
    outtree[${i}]="Radion_m300_aajj_8TeV"
    typ[${i}]="-300"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m400_aajj_8TeV"
    outtree[${i}]="Radion_m400_aajj_8TeV"
    typ[${i}]="-400"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m450_aajj_8TeV"
    outtree[${i}]="Radion_m450_aajj_8TeV"
    typ[${i}]="-450"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m500_aajj_8TeV"
    outtree[${i}]="Radion_m500_aajj_8TeV"
    typ[${i}]="-500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m550_aajj_8TeV"
    outtree[${i}]="Radion_m550_aajj_8TeV"
    typ[${i}]="-550"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m600_aajj_8TeV"
    outtree[${i}]="Radion_m600_aajj_8TeV"
    typ[${i}]="-600"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m700_aajj_8TeV"
    outtree[${i}]="Radion_m700_aajj_8TeV"
    typ[${i}]="-700"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m800_aajj_8TeV"
    outtree[${i}]="Radion_m800_aajj_8TeV"
    typ[${i}]="-800"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m900_aajj_8TeV"
    outtree[${i}]="Radion_m900_aajj_8TeV"
    typ[${i}]="-900"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Radion_m1000_aajj_8TeV"
    outtree[${i}]="Radion_m1000_aajj_8TeV"
    typ[${i}]="-1000"
    CS[${i}]="0"
fi

### MSSM samples
if [ ${doMSSM} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="MSSM_H_m260_8TeV"
    outtree[${i}]="MSSM_m260_8TeV"
    typ[${i}]="-260"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="MSSM_H_m300_8TeV"
    outtree[${i}]="MSSM_m300_8TeV"
    typ[${i}]="-300"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="MSSM_H_m350_8TeV"
    outtree[${i}]="MSSM_m350_8TeV"
    typ[${i}]="-350"
    CS[${i}]="0"
fi

##### GRAVITON
if [ ${doGraviton} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Graviton_m300_8TeV"
    outtree[${i}]="Graviton_m300_8TeV"
    typ[${i}]="-300"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Graviton_m500_8TeV"
    outtree[${i}]="Graviton_m500_8TeV"
    typ[${i}]="-500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Graviton_m700_8TeV"
    outtree[${i}]="Graviton_m700_8TeV"
    typ[${i}]="-700"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Graviton_m1000_8TeV"
    outtree[${i}]="Graviton_m1000_8TeV"
    typ[${i}]="-1000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_XHH.root"
    tree[${i}]="Graviton_m1500_8TeV"
    outtree[${i}]="Graviton_m1500_8TeV"
    typ[${i}]="-1500"
    CS[${i}]="0"
fi

### More graviton samples (narrow width)
if [ ${doGravitonMore} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Graviton_m270_LR3tev_8TeV"
    outtree[${i}]="Graviton_m270_LR3tev_8TeV"
    typ[${i}]="-270"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Graviton_m300_LR3tev_8TeV"
    outtree[${i}]="Graviton_m300_LR3tev_8TeV"
    typ[${i}]="-300"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Graviton_m350_LR3tev_8TeV"
    outtree[${i}]="Graviton_m350_LR3tev_8TeV"
    typ[${i}]="-350"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Graviton_m400_LR3tev_8TeV"
    outtree[${i}]="Graviton_m400_LR3tev_8TeV"
    typ[${i}]="-400"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Graviton_m450_LR3tev_8TeV"
    outtree[${i}]="Graviton_m450_LR3tev_8TeV"
    typ[${i}]="-450"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/XHH.root"
    tree[${i}]="Graviton_m700_LR3tev_8TeV"
    outtree[${i}]="Graviton_m700_LR3tev_8TeV"
    typ[${i}]="-700"
    CS[${i}]="0"
fi
    
##### SM Higgs
if [ ${doSMHiggs} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_SMHiggs.root"
    tree[${i}]="ggh_m125_powheg_8TeV"
    outtree[${i}]="ggh_m125_powheg_8TeV"
    typ[${i}]="-1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_SMHiggs.root"
    tree[${i}]="vbf_m125_8TeV"
    outtree[${i}]="vbf_m125_8TeV"
    typ[${i}]="-1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_SMHiggs.root"
    tree[${i}]="wzh_m125_8TeV_wh"
    outtree[${i}]="wzh_m125_8TeV_wh"
    typ[${i}]="-1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_SMHiggs.root"
    tree[${i}]="wzh_m125_8TeV_zh"
    outtree[${i}]="wzh_m125_8TeV_zh"
    typ[${i}]="-1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_SMHiggs.root"
    tree[${i}]="tth_m125_8TeV"
    outtree[${i}]="tth_m125_8TeV"
    typ[${i}]="-1"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_SMHiggs.root"
    tree[${i}]="bbh_m125_8TeV"
    outtree[${i}]="bbh_m125_8TeV"
    typ[${i}]="-1"
    CS[${i}]="0"
fi
if [ ${doExtraSMHiggs} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_tree_v10/SMHiggs.root"
    tree[${i}]="ggh_m125_minlo_8TeV"
    outtree[${i}]="ggh_m125_minlo_8TeV"
    typ[${i}]="-1"
    CS[${i}]="0"
fi

##### SM di-Higgs
if [ ${doSMdiHiggs} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_SMHiggs.root"
    tree[${i}]="ggHH_8TeV"
    outtree[${i}]="ggHH_8TeV"
    typ[${i}]="-2"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_1d0_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-500110000"
    CS[${i}]="0"
 fi

##### Diphoton backgrounds
if [ ${doDiphotonBackgrounds} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="qcd_30_8TeV_ff"
    outtree[${i}]="qcd_30_8TeV_ff"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="qcd_40_8TeV_ff"
    outtree[${i}]="qcd_40_8TeV_ff"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="qcd_30_8TeV_pf"
    outtree[${i}]="qcd_30_8TeV_pf"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="qcd_40_8TeV_pf"
    outtree[${i}]="qcd_40_8TeV_pf"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="gjet_20_8TeV_pf"
    outtree[${i}]="gjet_20_8TeV_pf"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="gjet_40_8TeV_pf"
    outtree[${i}]="gjet_40_8TeV_pf"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="diphojet_sherpa_8TeV"
    outtree[${i}]="diphojet_sherpa_8TeV"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="DYJetsToLL"
    outtree[${i}]="DYJetsToLL"
    typ[${i}]="1"
    CS[${i}]="0"
fi
  
if [ ${doRareBackgrounds} == 1 ]
then  
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="LNuGG_FSR_8TeV"
    outtree[${i}]="LNuGG_FSR_8TeV"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="LNuGG_ISR_8TeV"
    outtree[${i}]="LNuGG_ISR_8TeV"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="ttGG_8TeV"
    outtree[${i}]="ttGG_8TeV"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="tGG_8TeV"
    outtree[${i}]="tGG_8TeV"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="TTGJets_8TeV"
    outtree[${i}]="TTGJets_8TeV"
    typ[${i}]="1"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/tree_Bkg_v4.root"
    tree[${i}]="ZGToLLG_8TeV"
    outtree[${i}]="ZGToLLG_8TeV"
    typ[${i}]="1"
    CS[${i}]="0"
fi
    
##### Anomalous non-resonant HH scenarios
##### Type convention using 9 digits, 0 for + and 1 for -
##### -5 (1 sign + XY digits for lambda) (X.YZ digits for Yt) (1 sign + X digits for c2)
if [ ${doAnomalousHH} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_0d0_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-500010000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_10_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-501010000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_15_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-501510000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_20_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-502010000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_2_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_2_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-500210000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-511010000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-511510000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d0_c2_0d0"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d0_c2_0d0_8TeV"
    typ[${i}]="-512010000"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_0d75_c2_0d0"
    outtree[${i}]="ggHH_Lam_0d0_Yt_0d75_c2_0d0_8TeV"
    typ[${i}]="-500007500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d25_c2_0d0"
    outtree[${i}]="ggHH_Lam_0d0_Yt_1d25_c2_0d0_8TeV"
    typ[${i}]="-500012500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_0d75_c2_0d0"
    outtree[${i}]="ggHH_Lam_10_Yt_0d75_c2_0d0_8TeV"
    typ[${i}]="-501007500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d25_c2_0d0"
    outtree[${i}]="ggHH_Lam_10_Yt_1d25_c2_0d0_8TeV"
    typ[${i}]="-501012500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_0d75_c2_0d0"
    outtree[${i}]="ggHH_Lam_15_Yt_0d75_c2_0d0_8TeV"
    typ[${i}]="-501507500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d25_c2_0d0"
    outtree[${i}]="ggHH_Lam_15_Yt_1d25_c2_0d0_8TeV"
    typ[${i}]="-501512500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_0d75_c2_0d0"
    outtree[${i}]="ggHH_Lam_20_Yt_0d75_c2_0d0_8TeV"
    typ[${i}]="-502007500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d25_c2_0d0"
    outtree[${i}]="ggHH_Lam_20_Yt_1d25_c2_0d0_8TeV"
    typ[${i}]="-502012500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_0d75_c2_0d0"
    outtree[${i}]="ggHH_Lam_m10_Yt_0d75_c2_0d0_8TeV"
    typ[${i}]="-511007500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d25_c2_0d0"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d25_c2_0d0_8TeV"
    typ[${i}]="-511012500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_0d75_c2_0d0"
    outtree[${i}]="ggHH_Lam_m15_Yt_0d75_c2_0d0_8TeV"
    typ[${i}]="-511507500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d25_c2_0d0"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d25_c2_0d0_8TeV"
    typ[${i}]="-511512500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_0d75_c2_0d0"
    outtree[${i}]="ggHH_Lam_m20_Yt_0d75_c2_0d0_8TeV"
    typ[${i}]="-512007500"
    CS[${i}]="0"
    
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d25_c2_0d0"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d25_c2_0d0_8TeV"
    typ[${i}]="-512012500"
    CS[${i}]="0"
fi

if [ ${doAnomalousHH_C2} == 1 ]
then
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_0d75_c2_2_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_0d75_c2_2_8TeV"
    typ[${i}]="-500007502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_0d75_c2_3_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_0d75_c2_3_8TeV"
    typ[${i}]="-500007503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_0d75_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_0d75_c2_m2_8TeV"
    typ[${i}]="-500007512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_0d75_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_0d75_c2_m3_8TeV"
    typ[${i}]="-500007513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d0_c2_2_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_1d0_c2_2_8TeV"
    typ[${i}]="-500010002"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d0_c2_3_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_1d0_c2_3_8TeV"
    typ[${i}]="-500010003"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d0_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_1d0_c2_m2_8TeV"
    typ[${i}]="-500010012"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d0_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_0d0_Yt_1d0_c2_m3_8TeV"
    typ[${i}]="-500010013"
    CS[${i}]="0"

#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d25_c2_2_v2"
#    outtree[${i}]="ggHH_Lam_0d0_Yt_1d25_c2_2_8TeV"
#    typ[${i}]="-500012502"
#    CS[${i}]="0"
#
#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d25_c2_3_v2"
#    outtree[${i}]="ggHH_Lam_0d0_Yt_1d25_c2_3_8TeV"
#    typ[${i}]="-500012503"
#    CS[${i}]="0"
#
#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d25_c2_m2_v2"
#    outtree[${i}]="ggHH_Lam_0d0_Yt_1d25_c2_m2_8TeV"
#    typ[${i}]="-500012512"
#    CS[${i}]="0"
#
#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_0d0_Yt_1d25_c2_m3_v2"
#    outtree[${i}]="ggHH_Lam_0d0_Yt_1d25_c2_m3_8TeV"
#    typ[${i}]="-500012513"
#    CS[${i}]="0"
#
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_0d75_c2_2_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_0d75_c2_2_8TeV"
    typ[${i}]="-501007502"
    CS[${i}]="0"

#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_0d75_c2_3_v2"
#    outtree[${i}]="ggHH_Lam_10_Yt_0d75_c2_3_8TeV"
#    typ[${i}]="-501007503"
#    CS[${i}]="0"
#
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_0d75_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_0d75_c2_m2_8TeV"
    typ[${i}]="-501007512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_0d75_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_0d75_c2_m3_8TeV"
    typ[${i}]="-501007513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d0_c2_2_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_1d0_c2_2_8TeV"
    typ[${i}]="-501010002"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d0_c2_3_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_1d0_c2_3_8TeV"
    typ[${i}]="-501010003"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d0_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_1d0_c2_m2_8TeV"
    typ[${i}]="-501010012"
    CS[${i}]="0"

#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d0_c2_m3_v2"
#    outtree[${i}]="ggHH_Lam_10_Yt_1d0_c2_m3_8TeV"
#    typ[${i}]="-501010013"
#    CS[${i}]="0"
#
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d25_c2_2_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_1d25_c2_2_8TeV"
    typ[${i}]="-501012502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d25_c2_3_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_1d25_c2_3_8TeV"
    typ[${i}]="-501012503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d25_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_1d25_c2_m2_8TeV"
    typ[${i}]="-501012512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_10_Yt_1d25_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_10_Yt_1d25_c2_m3_8TeV"
    typ[${i}]="-501012513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_0d75_c2_2_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_0d75_c2_2_8TeV"
    typ[${i}]="-501507502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_0d75_c2_3_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_0d75_c2_3_8TeV"
    typ[${i}]="-501507503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_0d75_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_0d75_c2_m2_8TeV"
    typ[${i}]="-501507512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_0d75_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_0d75_c2_m3_8TeV"
    typ[${i}]="-501507513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d0_c2_2_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d0_c2_2_8TeV"
    typ[${i}]="-501510002"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d0_c2_3_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d0_c2_3_8TeV"
    typ[${i}]="-501510003"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d0_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d0_c2_m2_8TeV"
    typ[${i}]="-501510012"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d0_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d0_c2_m3_8TeV"
    typ[${i}]="-501510013"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d25_c2_2_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d25_c2_2_8TeV"
    typ[${i}]="-501512502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d25_c2_3_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d25_c2_3_8TeV"
    typ[${i}]="-501512503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d25_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d25_c2_m2_8TeV"
    typ[${i}]="-501512512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_15_Yt_1d25_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_15_Yt_1d25_c2_m3_8TeV"
    typ[${i}]="-501512513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_0d75_c2_2_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_0d75_c2_2_8TeV"
    typ[${i}]="-502007502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_0d75_c2_3_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_0d75_c2_3_8TeV"
    typ[${i}]="-502007503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_0d75_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_0d75_c2_m2_8TeV"
    typ[${i}]="-502007512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_0d75_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_0d75_c2_m3_8TeV"
    typ[${i}]="-502007513"
    CS[${i}]="0"

#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d0_c2_2_v2"
#    outtree[${i}]="ggHH_Lam_20_Yt_1d0_c2_2_8TeV"
#    typ[${i}]="-502010002"
#    CS[${i}]="0"
#
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d0_c2_3_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_1d0_c2_3_8TeV"
    typ[${i}]="-502010003"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d0_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_1d0_c2_m2_8TeV"
    typ[${i}]="-502010012"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d0_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_1d0_c2_m3_8TeV"
    typ[${i}]="-502010013"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d25_c2_2_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_1d25_c2_2_8TeV"
    typ[${i}]="-502012502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d25_c2_3_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_1d25_c2_3_8TeV"
    typ[${i}]="-502012503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d25_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_20_Yt_1d25_c2_m2_8TeV"
    typ[${i}]="-502012512"
    CS[${i}]="0"

#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_20_Yt_1d25_c2_m3_v2"
#    outtree[${i}]="ggHH_Lam_20_Yt_1d25_c2_m3_8TeV"
#    typ[${i}]="-502012513"
#    CS[${i}]="0"
#
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_0d75_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_0d75_c2_2_8TeV"
    typ[${i}]="-511007502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_0d75_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_0d75_c2_3_8TeV"
    typ[${i}]="-511007503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_0d75_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_0d75_c2_m2_8TeV"
    typ[${i}]="-511007512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_0d75_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_0d75_c2_m3_8TeV"
    typ[${i}]="-511007513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d0_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d0_c2_2_8TeV"
    typ[${i}]="-511010002"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d0_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d0_c2_3_8TeV"
    typ[${i}]="-511010003"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d0_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d0_c2_m2_8TeV"
    typ[${i}]="-511010012"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d0_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d0_c2_m3_8TeV"
    typ[${i}]="-511010013"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d25_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d25_c2_2_8TeV"
    typ[${i}]="-511012502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d25_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d25_c2_3_8TeV"
    typ[${i}]="-511012503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d25_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d25_c2_m2_8TeV"
    typ[${i}]="-511012512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m10_Yt_1d25_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m10_Yt_1d25_c2_m3_8TeV"
    typ[${i}]="-511012513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_0d75_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_0d75_c2_2_8TeV"
    typ[${i}]="-511507502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_0d75_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_0d75_c2_3_8TeV"
    typ[${i}]="-511507503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_0d75_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_0d75_c2_m2_8TeV"
    typ[${i}]="-511507512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_0d75_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_0d75_c2_m3_8TeV"
    typ[${i}]="-511507513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d0_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d0_c2_2_8TeV"
    typ[${i}]="-511510002"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d0_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d0_c2_3_8TeV"
    typ[${i}]="-511510003"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d0_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d0_c2_m2_8TeV"
    typ[${i}]="-511510012"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d0_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d0_c2_m3_8TeV"
    typ[${i}]="-511510013"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d25_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d25_c2_2_8TeV"
    typ[${i}]="-511512502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d25_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d25_c2_3_8TeV"
    typ[${i}]="-511512503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d25_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d25_c2_m2_8TeV"
    typ[${i}]="-511512512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m15_Yt_1d25_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m15_Yt_1d25_c2_m3_8TeV"
    typ[${i}]="-511512513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_0d75_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_0d75_c2_2_8TeV"
    typ[${i}]="-512007502"
    CS[${i}]="0"

#    i=$((${i} + 1))
#    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
#    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_0d75_c2_3_v2"
#    outtree[${i}]="ggHH_Lam_m20_Yt_0d75_c2_3_8TeV"
#    typ[${i}]="-512007503"
#    CS[${i}]="0"
#
    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_0d75_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_0d75_c2_m2_8TeV"
    typ[${i}]="-512007512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_0d75_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_0d75_c2_m3_8TeV"
    typ[${i}]="-512007513"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d0_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d0_c2_2_8TeV"
    typ[${i}]="-512010002"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d0_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d0_c2_3_8TeV"
    typ[${i}]="-512010003"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d0_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d0_c2_m2_8TeV"
    typ[${i}]="-512010012"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d0_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d0_c2_m3_8TeV"
    typ[${i}]="-512010013"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d25_c2_2_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d25_c2_2_8TeV"
    typ[${i}]="-512012502"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d25_c2_3_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d25_c2_3_8TeV"
    typ[${i}]="-512012503"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d25_c2_m2_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d25_c2_m2_8TeV"
    typ[${i}]="-512012512"
    CS[${i}]="0"

    i=$((${i} + 1))
    infile[${i}]="radion_redu_12_tree_10/ggHH_anomalous_v3.root"
    tree[${i}]="HH_bbaa_8TeV_Lam_m20_Yt_1d25_c2_m3_v2"
    outtree[${i}]="ggHH_Lam_m20_Yt_1d25_c2_m3_8TeV"
    typ[${i}]="-512012513"
    CS[${i}]="0"
fi

itot=${i}
masscutsuffix[0]="noMassCut"
regsuffix[0]="noRegression"
if [[ ${numRegFiles} != 0 ]]
then
    regsuffix[0]="withRegression"
fi 
imasscut=0

for ireg in `echo "0"`
do
    folder="${today}_selection_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_${version}"
    mkdir -p ${folder}
    for isample in `seq 0 ${itot}`
    do
        applyCS="${CS[${isample}]}"
        if [ ${nJackknife} == 0 ]
        then
            printCutFlow="1"
        else
            printCutFlow="0"
        fi
        echo -e "isample= ${isample} / ${itot}\tinfile= ${infile[${isample}]}\ttree= ${tree[${isample}]}\touttree= ${outtree[${isample}]}\ttyp= ${typ[${isample}]}\tireg= ${ireg}\tapplyCS= ${applyCS}"
        file="${outtree[${isample}]}_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_${version}"
        if [[ "${applyCS}" == "1" ]]
        then
            file="${outtree[${isample}]}_${regsuffix[${ireg}]}_${masscutsuffix[${imasscut}]}_controlSample_${version}"
            printCutFlow="0"
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
--regressionFilePath ${regFilePath} \
--numberOfRegressionFiles ${numRegFiles} \
--type ${typ[${isample}]} \
--removeUndefinedBtagSF ${removeZeroFlavour} \
--applyMassCuts ${imasscut} \
--applyPhotonIDControlSample ${applyCS} \
--printCutFlow ${printCutFlow} \
--keep0btag ${keep0btag} \
--whichPhotonID ${whichPhotonID} \
--nJackknife ${nJackknife} \
--iJackknife ${iJackknife} \
--FTR14001_style ${FTR14001_style} \
2> ${folder}/${file}.eo | tee ${folder}/${file}.eo | egrep '%|entries'

    if [[ "${printCutFlow}" == "1" ]]
    then
        mv cutFlow_${outtree[${isample}]}.dat ${folder}/
    fi

    sleep 1 # to avoid eos stress
    done # isample
done #ireg



