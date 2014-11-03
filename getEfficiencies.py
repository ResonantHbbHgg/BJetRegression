#!/usr/bin/env python
import numpy
import ROOT
from ROOT import gROOT
from ROOT import TCanvas
from ROOT import TChain
from ROOT import TGraphErrors
from ROOT import TLatex, TLegend
gROOT.SetBatch()
signals = []
#nstrategies = 23#-
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m20_Yt_1d0_c2_0d0_8TeV_m0.root", 18997, 0, 1])#18797
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m15_Yt_1d0_c2_0d0_8TeV_m0.root", 17800, 0, 1])#17700
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m10_Yt_1d0_c2_0d0_8TeV_m0.root", 17798, 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_0d0_Yt_1d0_c2_0d0_8TeV_m0.root", 19800, 0, 1])#19600
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV_m0.root", 20000, 0, 1])
##signals.append(["v34_alpha_fitToMgg_noKinFit/ggHH_Lam_1_8TeV.root", 19706, 1, 0.0000212])
##sanity test ggHH_Lam_1_8TeV_m0.root=ggHH_8TeV_m260.root
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_2_Yt_1d0_c2_0d0_8TeV_m0.root", 18997, 0, 1])#17999
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_10_Yt_1d0_c2_0d0_8TeV_m0.root", 19300, 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_15_Yt_1d0_c2_0d0_8TeV_m0.root", 17599 , 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_20_Yt_1d0_c2_0d0_8TeV_m0.root", 17599,  0, 1])#17399
#
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m20_Yt_0d75_c2_0d0_8TeV_m0.root", 19898, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m15_Yt_0d75_c2_0d0_8TeV_m0.root", 19799, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m10_Yt_0d75_c2_0d0_8TeV_m0.root", 19200, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_0d0_Yt_0d75_c2_0d0_8TeV_m0.root", 19800, 2, 1])
#signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV_m0.root", 20000, 1, 0.0000212])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_10_Yt_0d75_c2_0d0_8TeV_m0.root", 18999, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_15_Yt_0d75_c2_0d0_8TeV_m0.root", 19599, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_20_Yt_0d75_c2_0d0_8TeV_m0.root", 19196, 2, 1])
#
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m20_Yt_1d25_c2_0d0_8TeV_m0.root", 19797, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m15_Yt_1d25_c2_0d0_8TeV_m0.root", 19799, 3, 1])#19699
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m10_Yt_1d25_c2_0d0_8TeV_m0.root", 19099, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_0d0_Yt_1d25_c2_0d0_8TeV_m0.root", 19597, 3, 1])#19397
##signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV.root", 20000, 1, 0.0000212])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_10_Yt_1d25_c2_0d0_8TeV_m0.root", 19197, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_15_Yt_1d25_c2_0d0_8TeV_m0.root", 19299 , 3, 1])#19099
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_20_Yt_1d25_c2_0d0_8TeV_m0.root", 19098, 3, 1])
#
nstrategies = 19
# For the efficiency plot
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_0d75_c2_m3_8TeV.root", 17398, 2, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_0d75_c2_m3_8TeV.root", 16398, 2, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_0d75_c2_m3_8TeV.root", 15699, 2, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_0d75_c2_m3_8TeV.root", 15699, 2, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_0d75_c2_m3_8TeV.root", 20000, 2, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_0d75_c2_m3_8TeV.root", 8199, 2, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_0d75_c2_m3_8TeV.root", 16899, 2, 1])

signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d0_c2_m3_8TeV.root", 15699, 0, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d0_c2_m3_8TeV.root", 16197, 0, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d0_c2_m3_8TeV.root", 13699, 0, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_1d0_c2_m3_8TeV.root", 19900, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d0_c2_m3_8TeV.root", 17197, 0, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d0_c2_m3_8TeV.root", 15496, 0, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d0_c2_m3_8TeV.root", 17899, 0, 1])

signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d25_c2_m3_8TeV.root", 15398, 3, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d25_c2_m3_8TeV.root", 14799, 3, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d25_c2_m3_8TeV.root", 19800, 3, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d25_c2_m3_8TeV.root", 19900, 3, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d25_c2_m3_8TeV.root", 16298, 3, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d25_c2_m3_8TeV.root", 17698, 3, 1])
#nstrategies = 20
# For the efficiency plot
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_0d75_c2_m2_8TeV.root", 9698, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_0d75_c2_m2_8TeV.root", 16399, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_0d75_c2_m2_8TeV.root", 20000, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_0d75_c2_m2_8TeV.root", 16598, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_0d75_c2_m2_8TeV.root", 20000, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_0d75_c2_m2_8TeV.root", 15599, 2, 1])
signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_0d75_c2_m2_8TeV.root", 7197, 2, 1])

#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d0_c2_m2_8TeV.root", 18499, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d0_c2_m2_8TeV.root", 20000, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d0_c2_m2_8TeV.root", 18899, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_1d0_c2_m2_8TeV.root", 18498, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d0_c2_m2_8TeV.root", 17197, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d0_c2_m2_8TeV.root", 17795, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d0_c2_m2_8TeV.root", 16897, 0, 1])

#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d25_c2_m2_8TeV.root", 18499, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d25_c2_m2_8TeV.root", 18499, 3, 1])
##signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d25_c2_m2_8TeV.root", 20000, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d25_c2_m2_8TeV.root", 20000, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d25_c2_m2_8TeV.root", 20000, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d25_c2_m2_8TeV.root", 17698, 3, 1])
#nstrategies = 21
# For the efficiency plot
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_0d75_c2_2_8TeV.root", 20000, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_0d75_c2_2_8TeV.root", 18699, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_0d75_c2_2_8TeV.root", 17499, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_0d75_c2_2_8TeV.root", 17199, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_0d75_c2_2_8TeV.root", 17499, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_0d75_c2_2_8TeV.root", 18699, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_0d75_c2_2_8TeV.root", 20000, 2, 1])
#
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d0_c2_2_8TeV.root", 16598, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d0_c2_2_8TeV.root", 18797, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d0_c2_2_8TeV.root", 8598, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_1d0_c2_2_8TeV.root", 15496, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d0_c2_2_8TeV.root", 7799, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d0_c2_2_8TeV.root", 9698, 0, 1])
##signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d0_c2_2_8TeV.root", 15999, 0, 1])
#
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d25_c2_2_8TeV.root", 16397, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d25_c2_2_8TeV.root", 16295, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d25_c2_2_8TeV.root", 16698, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d25_c2_2_8TeV.root", 7093, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d25_c2_2_8TeV.root", 13699, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d25_c2_2_8TeV.root", 16799, 3, 1])
##nstrategies = 19
# For the efficiency plot
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_0d75_c2_3_8TeV.root", 16498, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_0d75_c2_3_8TeV.root", 19800, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_0d75_c2_3_8TeV.root", 10396, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_0d75_c2_3_8TeV.root", 16398, 2, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_0d75_c2_3_8TeV.root", 20000, 2, 1])
#
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d0_c2_3_8TeV.root", 14198, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d0_c2_3_8TeV.root", 14297, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d0_c2_3_8TeV.root", 14096, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_0d0_Yt_1d0_c2_3_8TeV.root", 11898, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d0_c2_3_8TeV.root", 13897, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d0_c2_3_8TeV.root", 15299, 0, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d0_c2_3_8TeV.root", 15999, 0, 1])
#
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m20_Yt_1d25_c2_3_8TeV.root", 19499, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m15_Yt_1d25_c2_3_8TeV.root", 16499, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_m10_Yt_1d25_c2_3_8TeV.root", 14499, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_10_Yt_1d25_c2_3_8TeV.root", 19900, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_15_Yt_1d25_c2_3_8TeV.root", 7699, 3, 1])
#signals.append(["v35_fitToMgg_noKinFit_0_massCutVersion_4/ggHH_Lam_20_Yt_1d25_c2_3_8TeV.root", 15099, 3, 1])
#
#nstrategies = 23
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m20_Yt_1d0_c2_0d0_8TeV.root", 18997, 0, 1])#18797
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m15_Yt_1d0_c2_0d0_8TeV.root", 17800, 0, 1])#17700
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m10_Yt_1d0_c2_0d0_8TeV.root", 17798, 0, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_0d0_Yt_1d0_c2_0d0_8TeV.root", 19800, 0, 1])#19600
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV.root", 20000, 0, 1])
signals.append(["v34_alpha_fitToMgg_noKinFit/ggHH_Lam_1_8TeV.root", 19706, 1, 0.0000212])
##sanity test ggHH_Lam_1_8TeV.root=ggHH_8TeV_m260.root
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_2_Yt_1d0_c2_0d0_8TeV.root", 18997, 0, 1])#17999
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_10_Yt_1d0_c2_0d0_8TeV.root", 19300, 0, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_15_Yt_1d0_c2_0d0_8TeV.root", 17599 , 0, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_20_Yt_1d0_c2_0d0_8TeV.root", 17599,  0, 1])#17399
#
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m20_Yt_0d75_c2_0d0_8TeV.root", 19898, 2, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m15_Yt_0d75_c2_0d0_8TeV.root", 19799, 2, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m10_Yt_0d75_c2_0d0_8TeV.root", 19200, 2, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_0d0_Yt_0d75_c2_0d0_8TeV.root", 19800, 2, 1])
##signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV.root", 20000, 1, 0.0000212])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_10_Yt_0d75_c2_0d0_8TeV.root", 18999, 2, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_15_Yt_0d75_c2_0d0_8TeV.root", 19599, 2, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_20_Yt_0d75_c2_0d0_8TeV.root", 19196, 2, 1])
#
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m20_Yt_1d25_c2_0d0_8TeV.root", 19797, 3, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m15_Yt_1d25_c2_0d0_8TeV.root", 19799, 3, 1])#19699
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m10_Yt_1d25_c2_0d0_8TeV.root", 19099, 3, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_0d0_Yt_1d25_c2_0d0_8TeV.root", 19597, 3, 1])#19397
##signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV.root", 20000, 1, 0.0000212])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_10_Yt_1d25_c2_0d0_8TeV.root", 19197, 3, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_15_Yt_1d25_c2_0d0_8TeV.root", 19299 , 3, 1])#19099
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_20_Yt_1d25_c2_0d0_8TeV.root", 19098, 3, 1])
#
#signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_m260.root", 96880, 1, 0.0000212])
#signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_0hhh_m260.root", 96880, 0, 0.0000212])
#signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_1hhh_m260.root", 96880, 0, 0.0000212])
#signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_2hhh_m260.root", 96880, 0, 0.0000212])
#signals.append(["v29_fitToMgg_noKinFit/MSSM_m260_8TeV_m260.root", 300000, 0])
#signals.append(["v29_fitToMgg_noKinFit/Radion_m270_8TeV_m270.root", 19996, 0])
#signals.append(["v29_fitToMgg_noKinFit/Radion_m300_8TeV_m300.root", 19972, 0])
#signals.append(["v29_fitToMgg_noKinFit/Radion_m350_8TeV_m350.root", 18498, 0])
#signals.append(["v29_fitToMgg_noKinFit/Radion_m400_8TeV_m400.root", 19697, 0])
#signals.append(["v29_fitToMgg_noKinFit/Radion_m450_8TeV_m450.root", 19999, 0])
#signals.append(["v29_fitToMgg_noKinFit/Radion_m500_8TeV_m500.root", 19970, 0])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m400_8TeV_m400.root", 19697, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m450_8TeV_m450.root", 19999, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m500_8TeV_m500.root", 19970, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m550_8TeV_m550.root", 19995, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m600_8TeV_m600.root", 18197, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m650_8TeV_m650.root", 20000, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m700_8TeV_m700.root", 19969, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m800_8TeV_m800.root", 19999, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m900_8TeV_m900.root", 19996, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m1000_8TeV_m1000.root", 19951, 1])
#signals.append(["v29_fitToMggjj_withKinFit/Radion_m1100_8TeV_m1100.root", 19400, 1])
# radion samples
#signals.append(["v26_fitToMgg_noKinFit/Radion_m270_8TeV_m270.root", 19996])
#signals.append(["v26_fitToMgg_noKinFit/Radion_m300_8TeV_m300.root", 19972])
#signals.append(["v26_fitToMgg_noKinFit/Radion_m350_8TeV_m350.root", 18498])
#signals.append(["v26_fitToMgg_noKinFit/Radion_m400_8TeV_m400.root", 19697])
#signals.append(["v26_fitToMgg_noKinFit/Radion_m450_8TeV_m450.root", 19999])
#signals.append(["v26_fitToMgg_noKinFit/Radion_m500_8TeV_m500.root", 19970])
# graviton samples
#signals.append(["v26_fitToMgg_noKinFit/MSSM_m260_8TeV_m260.root", 300000])
#signals.append(["v26_fitToMgg_noKinFit/MSSM_m300_8TeV_m300.root", 299142])
#signals.append(["v26_fitToMgg_noKinFit/MSSM_m350_8TeV_m350.root", 299571])
# mssm samples
#signals.append(["v26_fitToMgg_noKinFit/Graviton_m300_8TeV_m300.root", 49941])
#signals.append(["v26_fitToMgg_noKinFit/Graviton_m500_8TeV_m500.root", 49905])
#signals.append(["", ])


print signals
lumi = 19706.
n = []
x = []
y = []
ex = []
ey = []




for i in range(nstrategies):
	n.append( 0 )
	x.append( [] )
	y.append( [] )
	ex.append( [] )
	ey.append( [] )

for signal, nprocessed, istrategy, sigma in signals:
	chain = TChain("TCVARS")
	chain.Add(signal)
	for icut in range(2):
		totw = 0
		for ievt in xrange(chain.GetEntries()):
			chain.GetEntry(ievt)
			if icut == 0 and (chain.cut_based_ct == 0 or chain.cut_based_ct == 1):
				totw += chain.evWeight
#				totw += 1
			elif icut == 1 and chain.cut_based_ct == 0:
				totw += chain.evWeight
#				totw += 1
# totw *= 1000. # due to hand-put 1000 factor to weight for limit settings numerical precision
		print totw, nprocessed, lumi
#		eff = float(totw )/ nprocessed / float(sigma)
		eff = totw / lumi  / sigma
#		eff = totw / nprocessed / sigma
		deff = (eff * (1 - eff) / nprocessed)**0.5
		print eff*100, deff*100
#		eff = totw / lumi  / sigma
#		deff = 0.05

# print signal, round(eff * 100, 2), round(deff * 100,2)
# print signal.split("_")[1].split("m")[3], "&", "$", round(eff * 100, 2), "$", "&", "$\pm", round(deff * 100,2), "$"
		lambda_hhh = 1.0
		try:
			lambda_hhh = float(signal.split("ggHH_Lam_")[1].split("_8TeV.root")[0].split("_Yt")[0])
			print "lambda_hhh= ", lambda_hhh, "&", "$", round(eff * 100, 5), "$", "&", "$\pm", round(deff * 100, 5), "$"
		except ValueError:
			try:
				lambda_hhh = float(signal.split("ggHH_Lam_")[1].split("_8TeV.root")[0].split("d0_Yt")[0])
				print "lambda_hhh= ", lambda_hhh, "&", "$", round(eff * 100, 5), "$", "&", "$\pm", round(deff * 100, 5), "$"			
			except ValueError:
				try:
					lambda_hhh = -float(signal.split("ggHH_Lam_m")[1].split("_8TeV.root")[0].split("_Yt")[0])
				except ValueError:
					lambda_hhh = 1.0
				print "lambda_hhh= ", lambda_hhh, "&", "$", round(eff * 100, 5), "$", "&", "$\pm", round(deff * 100, 5), "$"
		n[2*istrategy+icut] = n[2*istrategy+icut]+1
		x[2*istrategy+icut].append( lambda_hhh )
		y[2*istrategy+icut].append(round(eff * 100, 2))
		ex[2*istrategy+icut].append(0.)
		ey[2*istrategy+icut].append(round(deff * 100,2))
	del chain

for istrategy in range(nstrategies):
	x[istrategy] = numpy.asarray(x[istrategy], dtype='f')
	y[istrategy] = numpy.asarray(y[istrategy], dtype='f')
	ex[istrategy] = numpy.asarray(ex[istrategy], dtype='f')
	ey[istrategy] = numpy.asarray(ey[istrategy], dtype='f')
# Removing y errors
	ey[istrategy] = ex[istrategy]
#gROOT.ProcessLine(".x CMSStyle.C")
gROOT.ProcessLine(".x setTDRStyle.C")
c1 = TCanvas()
#legend = TLegend(0.15, 0.64, 0.50, 0.94, "HH #rightarrow #gamma#gammab#bar{b}")
legend = TLegend(0.15, 0.72, 0.50, 0.94, "HH #rightarrow #gamma#gammab#bar{b}")
#lower[0.]{High+Medium Purity}}")
#legend = TLegend(0.15, 0.72, 0.50, 0.92, "#splitline{X#rightarrow HH #rightarrow b#bar{b}#gamma#gamma}{High-purity category}")
#legend = TLegend(0.15, 0.72, 0.50, 0.92, "#splitline{X#rightarrow HH #rightarrow b#bar{b}#gamma#gamma}{Medium-purity category}")
legend.SetTextSize(0.03)
legend.SetFillColor(ROOT.kWhite)
legend.SetLineColor(ROOT.kWhite)
legend.SetShadowColor(ROOT.kWhite)
gr0 = TGraphErrors(n[0],x[0],y[0],ex[0],ey[0]);
gr1 = TGraphErrors(n[1],x[1],y[1],ex[1],ey[1]);
gr2 = TGraphErrors(n[2],x[2],y[2],ex[2],ey[2]);
gr3 = TGraphErrors(n[3],x[3],y[3],ex[3],ey[3]);
gr4 = TGraphErrors(n[4],x[4],y[4],ex[4],ey[4]);
gr5 = TGraphErrors(n[5],x[5],y[5],ex[5],ey[5]);
gr6 = TGraphErrors(n[6],x[6],y[6],ex[6],ey[6]);
gr7 = TGraphErrors(n[7],x[7],y[7],ex[7],ey[7]);

gr0.SetName("gr0")
gr1.SetName("gr1")
gr2.SetName("gr2")
gr3.SetName("gr3")
gr4.SetName("gr4")
gr5.SetName("gr5")
gr6.SetName("gr6")
gr7.SetName("gr7")

gr0.SetMarkerColor(ROOT.kGreen+3)
gr1.SetMarkerColor(ROOT.kGreen+3)
gr2.SetMarkerColor(ROOT.kBlue+2)
gr3.SetMarkerColor(ROOT.kBlue+2)
gr4.SetMarkerColor(ROOT.kRed+3)
gr5.SetMarkerColor(ROOT.kRed+3)
gr6.SetMarkerColor(ROOT.kBlack+2)
gr7.SetMarkerColor(ROOT.kBlack+2)

gr0.SetLineColor(ROOT.kGreen+3)
gr1.SetLineColor(ROOT.kGreen+3)
gr2.SetLineColor(ROOT.kBlue+2)
gr3.SetLineColor(ROOT.kBlue+2)
gr4.SetLineColor(ROOT.kRed+3)
gr5.SetLineColor(ROOT.kRed+3)
gr6.SetLineColor(ROOT.kBlack+2)
gr7.SetLineColor(ROOT.kBlack+2)

#gr0.SetFillColor(ROOT.kGreen-3)
#gr1.SetFillColor(ROOT.kBlue-2)
#gr2.SetFillColor(ROOT.kGreen-3)
#gr3.SetFillColor(ROOT.kBlue-2)
gr0.SetMarkerStyle(20)
gr1.SetMarkerStyle(24)
gr2.SetMarkerStyle(21)
gr3.SetMarkerStyle(25)
gr4.SetMarkerStyle(20)
gr5.SetMarkerStyle(24)
gr6.SetMarkerStyle(20)
gr7.SetMarkerStyle(24)

gr2.SetLineStyle(2)
gr3.SetLineStyle(2)
gr0.SetTitle("")
gr0.GetXaxis().SetTitle("#lambda_{hhh}")
gr0.GetYaxis().SetTitle("Signal selection efficiency (%)")
print min(numpy.amin(y[0]), numpy.amin(y[1])), max(numpy.amax(y[0]), numpy.amax(y[1]))
print min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1]))
gr0.SetMinimum(0.)
gr0.SetMaximum(80.)
gr0.GetXaxis().SetLimits(-25., 25.)
#gr0.SetMaximum(max(numpy.amax(y[0]), numpy.amax(y[1])))
#gr0.SetMinimum(min(numpy.amin(y[0]), numpy.amin(y[1])))
#gr0.GetXaxis().SetLimits(min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1])))
gr0.Draw("ALP");
gr1.Draw("LP");
#gr2.Draw("LP");
#gr3.Draw("LP");
gr4.Draw("LP");
gr5.Draw("LP");
gr6.Draw("LP");
gr7.Draw("LP");

legend.AddEntry(gr0.GetName(), "High+Medium Purity - Yt_1_c2_m3", "lp")
legend.AddEntry(gr1.GetName(), "High Purity - Yt_1", "lp")
#legend.AddEntry(gr2.GetName(), "High+Medium Purity - ggHH", "lp")
#legend.AddEntry(gr3.GetName(), "High Purity - ggHH", "lp")
legend.AddEntry(gr4.GetName(), "High+Medium Purity - Yt_0d75_c2_m3", "lp")
legend.AddEntry(gr5.GetName(), "High Purity -Yt_0d75", "lp")
legend.AddEntry(gr6.GetName(), "High+Medium Purity -Yt_1d25_c2_m3", "lp")
legend.AddEntry(gr7.GetName(), "High Purity -Yt_1d25", "lp")

legend.Draw()
latexLabel = TLatex()
latexLabel.SetTextSize(0.03)
latexLabel.SetNDC()
latexLabel.DrawLatex(0.15, 0.96, "CMS Work in progress")
latexLabel.DrawLatex(0.80, 0.96, "#sqrt{s} = 8 TeV")
c1.Print("2014-11-03_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_m3.pdf")
c1.Print("2014-10-03_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_m3.png")
c1.Print("2014-10-03_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_m3.root")
c1.Print("2014-10-03_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_m3.eps")
#c1.Print("eff_cat0.pdf")
#c1.Print("eff_cat0.png")
#c1.Print("eff_cat0.root")
#c1.Print("eff_cat1.pdf")
#c1.Print("eff_cat1.png")
#c1.Print("eff_cat1.root")
