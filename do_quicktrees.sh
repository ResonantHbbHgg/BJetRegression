#!/bin/bash

if [ quickTrees.cc -nt quickTrees.exe ] || [ quickTrees.h -nt quickTrees.exe ]
then
    echo "please recompile quickTrees"
    echo "make quickTrees.exe"
    exit 3000 
fi


version="v40"
today=`date +"0%Y-%m-%d"`
#set -x

inputversion="v22"
inputfolder="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v22/2014-12-12_selection_noRegression_noMassCut_v22/"
inputfolderReg="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v22/2014-12-12_selection_withRegression_noMassCut_v22/"
inputversionFTR="v22_FTR"
inputfolderFTR="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v22/2014-11-07_selection_noRegression_noMassCut_v22_FTR/"

# IMPORTANT NOTES:
# FOR NOW THE DEFAULT IS NO REGRESSION
# As of today (Oct. 13) regression will be produced as well by default for low mass when trees are available
 

# WHICH ANALYSIS TO PROCESS
doNonResonant=1
doResonantLowMass=1
doResonantLowMassWithReg=1
doResonantHighMass=1
# WHICH SAMPLES TO PROCESS (default is also running dataCS and diphoton-sherpa, minimum is data + signal)
doTheStrictMinimum=0
doAnomalousHHScenario1=1
doAnomalousHHScenario2=1


# OTHER GLOBAL SETTINGS, IN MOST USE CASE YOU SHOULD NOT TOUCH THIS
cutLevel=0
massCutVersion=4 # From Summer 14 cut update
controlSampleWeights="scales_2D_pt_data_4GeVbinning.root"
applyFTR14001=0
doUnnecessaryDirs=0

# INITIALIZATION OF THE PROCESSING LIST
i=-1

if [ ${doResonantHighMass} == 1 ]
then
    ##### PREPARE MGGJJ-FIT TREES
    kinfitlabel[0]="noKinFit"
    kinfitlabel[1]="withKinFit"
    kinfitjet[0]="base"
    kinfitjet[1]="kin"
    kinfitlist=""
    if [ ${doUnnecessaryDirs} == 0 ]
    then
        kinfitList="1"
    else
        kinfitList="0 1"
    fi
    for ikin in `echo $kinfitList`
    do
        outfolder="${version}_fitToMggjj_${kinfitlabel[${ikin}]}"
        mkdir -p ${outfolder}
        samplelist="Radion Graviton Data"
        if [ ${doTheStrictMinimum} == 0 ]
        then
            samplelist="${samplelist} DataCS diphojet_sherpa_8TeV DYJetsToLL gjet_40_8TeV_pf gjet_20_8TeV_pf qcd_40_8TeV_pf qcd_30_8TeV_pf qcd_40_8TeV_ff qcd_30_8TeV_ff LNuGG_FSR_8TeV LNuGG_ISR_8TeV ttGG_8TeV tGG_8TeV TTGJets_8TeV ZGToLLG_8TeV"
        fi
        for sample in `echo "${samplelist}"`
        do
            for mass in `echo "400 450 500 550 600 650 700 800 900 1000 1100"`
            do
                intree=${sample}
                outtree=${sample}
                itype="1"
                removeUndefinedBtagSF=0
                applyPhotonIDControlSample=0
                suffix=""
                extraline=""
#            extraline="--applyMjjCut 0  --applyMggCut 0 --applyMtotCut 0"
                if [ "${sample}" == "Radion" ]
                then
                    intree="${sample}_m${mass}_8TeV"
                    outtree="${sample}_m${mass}_8TeV"
                    itype="-${mass}"
                    removeUndefinedBtagSF=0
                elif [ "${sample}" == "Graviton" ]
                then
                    if [ "${mass}" == "500" ] || [ "${mass}" == "700" ] || [ "${mass}" == "1000" ]
                    then
                        intree="${sample}_m${mass}_8TeV"
                        outtree="${sample}_m${mass}_8TeV"
                        itype="-${mass}"
                        removeUndefinedBtagSF=0
                    else
                        continue
                    fi
                elif [ "${sample}" == "Data" ]
                then
                    itype="0"
                elif [ "${sample}" == "DataCS" ]
                then
                    itype="0"
                    intree="Data"
                    applyPhotonIDControlSample=1
                    suffix="controlSample_"
                fi
                i=$((${i} + 1))
                line[${i}]=""
                line[${i}]="${line[${i}]} --inputfile ${inputfolder}/${intree}_noRegression_noMassCut_${suffix}${inputversion}.root"
                line[${i}]="${line[${i}]} --inputtree ${intree}"
                line[${i}]="${line[${i}]} --outputtree TCVARS"
                line[${i}]="${line[${i}]} --outputfile ${outfolder}/${outtree}_m${mass}.root"
                line[${i}]="${line[${i}]} --type ${itype}"
                line[${i}]="${line[${i}]} --whichJet ${kinfitjet[${ikin}]}"
                line[${i}]="${line[${i}]} --fitStrategy mggjj"
                line[${i}]="${line[${i}]} --cutLevel ${cutLevel}"
                line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
                line[${i}]="${line[${i}]} --massCutVersion ${massCutVersion}"
                line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
                line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
                line[${i}]="${line[${i}]} ${extraline}"
                log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
            #    echo -e "i= ${i}\tline= ${line[${i}]}"
            done
        done
    done
fi
##### END OF MGGJJ-FIT TREES
    
if [ ${doResonantLowMass} == 1 ] || [ ${doResonantLowMassWithReg} == 1 ]
then
    ##### PREPARE LOW-MASS RESONANCE MGG-FIT TREES
    kinfitlabel[0]="noKinFit"
    kinfitlabel[1]="withKinFit"
    kinfitlabel[2]="withReg"
    kinfitlabel[3]="withRegKinFit"
    kinfitjet[0]="base"
    kinfitjet[1]="kin"
    kinfitjet[2]="reg"
    kinfitjet[3]="regkin"
    kinfitList=""
    if [ ${doResonantLowMass} == 1 ]
    then
	if [ ${doUnnecessaryDirs} == 0 ]
	    then
	    kinfitList="1"
	else
	    kinfitList="0 1"
	fi
    fi
    if [ ${doResonantLowMassWithReg} == 1 ]
    then
	if [ ${doUnnecessaryDirs} == 0 ]
        then
	    kinfitList="$kinfitList 3"
	else
	    kinfitList="$kinfitList 2 3"
	fi
    fi
    for ikin in `echo $kinfitList`
    do
        strategylist="mgg 2D"
        for fitStrategy in `echo "${strategylist}"`
        do
            if [ "$fitStrategy" == "mgg" ]
            then
                outfolder="${version}_fitToMgg_resSearch_${kinfitlabel[${ikin}]}"
            else
                outfolder="${version}_fitTo${fitStrategy}_resSearch_${kinfitlabel[${ikin}]}"
            fi
            mkdir -p ${outfolder}
            samplelist="Radion Graviton MSSM ggh_m125_powheg_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh tth_m125_8TeV bbh_m125_8TeV Data"
            if [ ${doTheStrictMinimum} == 0 ]
            then
                samplelist="${samplelist} DataCS diphojet_sherpa_8TeV DYJetsToLL gjet_40_8TeV_pf gjet_20_8TeV_pf qcd_40_8TeV_pf qcd_30_8TeV_pf qcd_40_8TeV_ff qcd_30_8TeV_ff LNuGG_FSR_8TeV LNuGG_ISR_8TeV ttGG_8TeV tGG_8TeV TTGJets_8TeV ZGToLLG_8TeV"
            fi
            for sample in `echo "${samplelist}"`
            do
                for mass in `echo "260 270 300 350 400 450 500"`
                do
                    intree=${sample}
                    outtree=${sample}
                    itype="1"
                    removeUndefinedBtagSF=0
                    applyPhotonIDControlSample=0
                    suffix=""
                    extraline=""
                    if [ "${sample}" == "Radion" ] 
                    then
                        if [ "${mass}" == "260" ]
                        then
                            continue
                        else
                            intree="${sample}_m${mass}_8TeV"
                            outtree="${sample}_m${mass}_8TeV"
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
                            outtree="${sample}_m${mass}_8TeV"
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
                            outtree="${sample}_m${mass}_8TeV"
                            itype="-${mass}"
                            removeUndefinedBtagSF=0
                        fi
                    elif [ "${sample}" == "Data" ]
                    then
                        itype="0"
                    elif [ "${sample}" == "DataCS" ]
                    then
                        itype="0"
                        applyPhotonIDControlSample=1
                        intree="Data"
                        suffix="controlSample_"
                    fi
                    i=$((${i} + 1))
                    line[${i}]=""
                    if [ "$ikin" == "0" ] || [ "$ikin" == "1" ]
                    then
                        line[${i}]="${line[${i}]} --inputfile ${inputfolder}/${intree}_noRegression_noMassCut_${suffix}${inputversion}.root"
                    else
                        line[${i}]="${line[${i}]} --inputfile ${inputfolderReg}/${intree}_withRegression_noMassCut_${suffix}${inputversion}.root"
                    fi            
                    line[${i}]="${line[${i}]} --inputtree ${intree}"
                    line[${i}]="${line[${i}]} --outputtree TCVARS"
                    line[${i}]="${line[${i}]} --outputfile ${outfolder}/${outtree}_m${mass}.root"
                    line[${i}]="${line[${i}]} --type ${itype}"
                    line[${i}]="${line[${i}]} --whichJet ${kinfitjet[${ikin}]}"
                    line[${i}]="${line[${i}]} --fitStrategy ${fitStrategy}"
                    line[${i}]="${line[${i}]} --cutLevel ${cutLevel}"
                    line[${i}]="${line[${i}]} --mass ${mass}"
                    line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
                    line[${i}]="${line[${i}]} --massCutVersion ${massCutVersion}"
                    line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
                    line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
                    line[${i}]="${line[${i}]} ${extraline}"
                    log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
            #    echo -e "i= ${i}\tline= ${line[${i}]}"
                done # end of loop on mass hyp
            done # end of loop on samples
        done # end of loop on fit strategy
    done # end of loop on kinfit scenarii
fi
##### END OF LOW-MASS RESONANCE TREES

#####
if [ ${doNonResonant} == 1 ]
    then
    ##### PREPARE NON-RESONANT MGG-FIT TREES
    kinfitlabel[0]="noKinFit"
    kinfitlabel[1]="withKinFit"
    kinfitjet[0]="base"
    kinfitjet[1]="kin"
    kinfitlist=""
    if [ ${doUnnecessaryDirs} == 0 ]
    then
        kinfitList="1"
    else
        kinfitList="0 1"
    fi
    for ikin in `echo $kinfitList`
    do
        strategylist="mgg 2D"
        if [ ${applyFTR14001} == 1 ]
        then
            strategylist="FTR14001"
            inputfolder="${inputfolderFTR}"
            inputversion="${inputversionFTR}"
        fi
        for fitStrategy in `echo "${strategylist}"`
        do
            if [ "$fitStrategy" == "mgg" ]
            then
                outfolder="${version}_fitToMgg_nonresSearch_${kinfitlabel[${ikin}]}"
            else
                outfolder="${version}_fitTo${fitStrategy}_nonresSearch_${kinfitlabel[${ikin}]}"
            fi
            mkdir -p ${outfolder}
            samplelist="ggHH_8TeV ggh_m125_powheg_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh tth_m125_8TeV bbh_m125_8TeV Data"
            if [ ${doTheStrictMinimum} == 0 ]
            then
                samplelist="${samplelist} DataCS diphojet_sherpa_8TeV DYJetsToLL gjet_40_8TeV_pf gjet_20_8TeV_pf qcd_40_8TeV_pf qcd_30_8TeV_pf qcd_40_8TeV_ff qcd_30_8TeV_ff LNuGG_FSR_8TeV LNuGG_ISR_8TeV ttGG_8TeV tGG_8TeV TTGJets_8TeV ZGToLLG_8TeV"
            fi
            if [ ${doAnomalousHHScenario1} == 1 ]
            then
                if [ ${doTheStrictMinimum} == 1 ]
                then
                    samplelist=""
                fi
                samplelist="${samplelist} ggHH_Lam_0d0_Yt_0d75_c2_0d0_8TeV ggHH_Lam_0d0_Yt_1d0_c2_0d0_8TeV ggHH_Lam_0d0_Yt_1d25_c2_0d0_8TeV ggHH_Lam_10_Yt_0d75_c2_0d0_8TeV ggHH_Lam_10_Yt_1d0_c2_0d0_8TeV ggHH_Lam_10_Yt_1d25_c2_0d0_8TeV ggHH_Lam_15_Yt_0d75_c2_0d0_8TeV ggHH_Lam_15_Yt_1d0_c2_0d0_8TeV ggHH_Lam_15_Yt_1d25_c2_0d0_8TeV ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV ggHH_Lam_20_Yt_0d75_c2_0d0_8TeV ggHH_Lam_20_Yt_1d0_c2_0d0_8TeV ggHH_Lam_20_Yt_1d25_c2_0d0_8TeV ggHH_Lam_2_Yt_1d0_c2_0d0_8TeV ggHH_Lam_m10_Yt_0d75_c2_0d0_8TeV ggHH_Lam_m10_Yt_1d0_c2_0d0_8TeV ggHH_Lam_m10_Yt_1d25_c2_0d0_8TeV ggHH_Lam_m15_Yt_0d75_c2_0d0_8TeV ggHH_Lam_m15_Yt_1d0_c2_0d0_8TeV ggHH_Lam_m15_Yt_1d25_c2_0d0_8TeV ggHH_Lam_m20_Yt_0d75_c2_0d0_8TeV ggHH_Lam_m20_Yt_1d0_c2_0d0_8TeV ggHH_Lam_m20_Yt_1d25_c2_0d0_8TeV ggHH_Lam_1d0_Yt_0d75_c2_0d0_8TeV ggHH_Lam_1d0_Yt_1d25_c2_0d0_8TeV ggHH_Lam_m2_Yt_1d0_c2_0d0_8TeV ggHH_Lam_3_Yt_1d0_c2_0d0_8TeV ggHH_Lam_5_Yt_1d0_c2_0d0_8TeV"
            fi
            if [ ${doAnomalousHHScenario2} == 1 ]
            then
                if [ ${doTheStrictMinimum} == 1 ]
                then
                    samplelist=""
                fi
                samplelist="${samplelist} ggHH_Lam_0d0_Yt_0d75_c2_2_8TeV  ggHH_Lam_0d0_Yt_0d75_c2_3_8TeV ggHH_Lam_0d0_Yt_0d75_c2_m2_8TeV ggHH_Lam_0d0_Yt_0d75_c2_m3_8TeV ggHH_Lam_0d0_Yt_1d0_c2_2_8TeV ggHH_Lam_0d0_Yt_1d0_c2_3_8TeV ggHH_Lam_0d0_Yt_1d0_c2_m2_8TeV ggHH_Lam_0d0_Yt_1d0_c2_m3_8TeV ggHH_Lam_10_Yt_0d75_c2_2_8TeV ggHH_Lam_10_Yt_0d75_c2_m2_8TeV ggHH_Lam_10_Yt_0d75_c2_m3_8TeV ggHH_Lam_10_Yt_1d0_c2_2_8TeV ggHH_Lam_10_Yt_1d0_c2_3_8TeV ggHH_Lam_10_Yt_1d0_c2_m2_8TeV ggHH_Lam_10_Yt_1d25_c2_2_8TeV ggHH_Lam_10_Yt_1d25_c2_3_8TeV ggHH_Lam_10_Yt_1d25_c2_m2_8TeV ggHH_Lam_10_Yt_1d25_c2_m3_8TeV ggHH_Lam_15_Yt_0d75_c2_2_8TeV ggHH_Lam_15_Yt_0d75_c2_3_8TeV ggHH_Lam_15_Yt_0d75_c2_m2_8TeV ggHH_Lam_15_Yt_0d75_c2_m3_8TeV ggHH_Lam_15_Yt_1d0_c2_2_8TeV ggHH_Lam_15_Yt_1d0_c2_3_8TeV ggHH_Lam_15_Yt_1d0_c2_m2_8TeV ggHH_Lam_15_Yt_1d0_c2_m3_8TeV ggHH_Lam_15_Yt_1d25_c2_2_8TeV ggHH_Lam_15_Yt_1d25_c2_3_8TeV ggHH_Lam_15_Yt_1d25_c2_m2_8TeV ggHH_Lam_15_Yt_1d25_c2_m3_8TeV ggHH_Lam_20_Yt_0d75_c2_2_8TeV ggHH_Lam_20_Yt_0d75_c2_3_8TeV ggHH_Lam_20_Yt_0d75_c2_m2_8TeV ggHH_Lam_20_Yt_0d75_c2_m3_8TeV ggHH_Lam_20_Yt_1d0_c2_3_8TeV ggHH_Lam_20_Yt_1d0_c2_m2_8TeV ggHH_Lam_20_Yt_1d0_c2_m3_8TeV ggHH_Lam_20_Yt_1d25_c2_2_8TeV ggHH_Lam_20_Yt_1d25_c2_3_8TeV ggHH_Lam_20_Yt_1d25_c2_m2_8TeV ggHH_Lam_m10_Yt_0d75_c2_2_8TeV ggHH_Lam_m10_Yt_0d75_c2_3_8TeV ggHH_Lam_m10_Yt_0d75_c2_m2_8TeV ggHH_Lam_m10_Yt_0d75_c2_m3_8TeV ggHH_Lam_m10_Yt_1d0_c2_2_8TeV ggHH_Lam_m10_Yt_1d0_c2_3_8TeV ggHH_Lam_m10_Yt_1d0_c2_m2_8TeV ggHH_Lam_m10_Yt_1d0_c2_m3_8TeV ggHH_Lam_m10_Yt_1d25_c2_2_8TeV ggHH_Lam_m10_Yt_1d25_c2_3_8TeV ggHH_Lam_m10_Yt_1d25_c2_m2_8TeV ggHH_Lam_m10_Yt_1d25_c2_m3_8TeV ggHH_Lam_m15_Yt_0d75_c2_2_8TeV ggHH_Lam_m15_Yt_0d75_c2_3_8TeV ggHH_Lam_m15_Yt_0d75_c2_m2_8TeV ggHH_Lam_m15_Yt_0d75_c2_m3_8TeV ggHH_Lam_m15_Yt_1d0_c2_2_8TeV ggHH_Lam_m15_Yt_1d0_c2_3_8TeV ggHH_Lam_m15_Yt_1d0_c2_m2_8TeV ggHH_Lam_m15_Yt_1d0_c2_m3_8TeV ggHH_Lam_m15_Yt_1d25_c2_2_8TeV ggHH_Lam_m15_Yt_1d25_c2_3_8TeV ggHH_Lam_m15_Yt_1d25_c2_m2_8TeV ggHH_Lam_m15_Yt_1d25_c2_m3_8TeV ggHH_Lam_m20_Yt_0d75_c2_2_8TeV ggHH_Lam_m20_Yt_0d75_c2_m2_8TeV ggHH_Lam_m20_Yt_0d75_c2_m3_8TeV ggHH_Lam_m20_Yt_1d0_c2_2_8TeV ggHH_Lam_m20_Yt_1d0_c2_3_8TeV ggHH_Lam_m20_Yt_1d0_c2_m2_8TeV ggHH_Lam_m20_Yt_1d0_c2_m3_8TeV ggHH_Lam_m20_Yt_1d25_c2_2_8TeV ggHH_Lam_m20_Yt_1d25_c2_3_8TeV ggHH_Lam_m20_Yt_1d25_c2_m2_8TeV ggHH_Lam_m20_Yt_1d25_c2_m3_8TeV ggHH_Lam_1d0_Yt_0d75_c2_m3_8TeV ggHH_Lam_1d0_Yt_1d0_c2_m3_8TeV ggHH_Lam_1d0_Yt_1d25_c2_m3_8TeV ggHH_Lam_1d0_Yt_0d75_c2_m2_8TeV ggHH_Lam_1d0_Yt_1d0_c2_m2_8TeV ggHH_Lam_1d0_Yt_1d25_c2_m2_8TeV ggHH_Lam_1d0_Yt_0d75_c2_2_8TeV ggHH_Lam_1d0_Yt_1d0_c2_2_8TeV ggHH_Lam_1d0_Yt_1d25_c2_2_8TeV ggHH_Lam_1d0_Yt_0d75_c2_3_8TeV ggHH_Lam_1d0_Yt_1d0_c2_3_8TeV ggHH_Lam_1d0_Yt_1d25_c2_3_8TeV ggHH_Lam_0d0_Yt_1d25_c2_m3_8TeV ggHH_Lam_0d0_Yt_1d25_c2_m2_8TeV ggHH_Lam_0d0_Yt_1d25_c2_2_8TeV ggHH_Lam_0d0_Yt_1d25_c2_3_8TeV ggHH_Lam_20_Yt_1d25_c2_m3_8TeV ggHH_Lam_10_Yt_1d0_c2_m3_8TeV ggHH_Lam_20_Yt_1d0_c2_2_8TeV ggHH_Lam_m20_Yt_0d75_c2_3_8TeV ggHH_Lam_10_Yt_0d75_c2_3_8TeV"
            fi

            for sample in `echo "${samplelist}"`
            do
                mass=0
                intree=${sample}
                outtree=${sample}
                itype="1"
                removeUndefinedBtagSF=0
                applyPhotonIDControlSample=0
                suffix=""
                extraline=""
                if [ "${sample}" == "Data" ]
                then
                    itype="0"
                elif [ "${sample}" == "DataCS" ]
                then
                    itype="0"
                    applyPhotonIDControlSample=1
                    intree="Data"
                    suffix="controlSample_"
                elif [[ ${sample} == ggHH_Lam* ]]
                then
                    itype="-500000000"
                fi
                i=$((${i} + 1))
                line[${i}]=""
                line[${i}]="${line[${i}]} --inputfile ${inputfolder}/${intree}_noRegression_noMassCut_${suffix}${inputversion}.root"
                line[${i}]="${line[${i}]} --inputtree ${intree}"
                line[${i}]="${line[${i}]} --outputtree TCVARS"
                line[${i}]="${line[${i}]} --outputfile ${outfolder}/${outtree}_m${mass}.root"
                line[${i}]="${line[${i}]} --type ${itype}"
                line[${i}]="${line[${i}]} --whichJet ${kinfitjet[${ikin}]}"
                line[${i}]="${line[${i}]} --fitStrategy ${fitStrategy}"
                line[${i}]="${line[${i}]} --cutLevel ${cutLevel}"
                line[${i}]="${line[${i}]} --mass ${mass}"
                line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
                line[${i}]="${line[${i}]} --massCutVersion ${massCutVersion}"
                line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
                line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
                line[${i}]="${line[${i}]} ${extraline}"
                log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
            #    echo -e "i= ${i}\tline= ${line[${i}]}"
            done
        done
    done # end of loop on kinfit scenarii
fi


#### PRODUCE EVERYTHING    
itot=${i}
for iline in `seq 0 ${itot}`
do
    echo -e "iline= ${iline} / ${itot}\t\t${line[${iline}]}"
    ./quickTrees.exe ${line[${iline}]} &> ${log[${iline}]}
done
