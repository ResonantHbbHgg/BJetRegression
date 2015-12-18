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
#/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit
#

nstrategies = 29
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_0d75_c2_0d0_8TeV_m0.root", 19898, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_0d75_c2_0d0_8TeV_m0.root", 19799, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_0d75_c2_0d0_8TeV_m0.root", 19200, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_0d75_c2_0d0_8TeV_m0.root", 19800, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_0d75_c2_0d0_8TeV_m0.root", 7599, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_0d75_c2_0d0_8TeV_m0.root", 18999, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_0d75_c2_0d0_8TeV_m0.root", 19599, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_0d75_c2_0d0_8TeV_m0.root", 19196, 2, 1])

signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d0_c2_0d0_8TeV_m0.root", 18997, 0, 1]) #18797
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d0_c2_0d0_8TeV_m0.root", 17800, 0, 1]) #17700
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d0_c2_0d0_8TeV_m0.root", 17798, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m2_Yt_1d0_c2_0d0_8TeV_m0.root", 20000, 0, 1]) #17999
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d0_c2_0d0_8TeV_m0.root", 19800, 0, 1]) #19600
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV_m0.root", 20000, 0, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit/ggHH_Lam_1_8TeV.root", 19706, 1, 0.0000212])
#Ssanity test ggHH_Lam_1_8TeV_m0.root=ggHH_8TeV_m260.root
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_2_Yt_1d0_c2_0d0_8TeV_m0.root", 18997, 0, 1]) #17999
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_3_Yt_1d0_c2_0d0_8TeV_m0.root", 8799, 0, 1]) #17999
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_5_Yt_1d0_c2_0d0_8TeV_m0.root", 17895, 0, 1]) #17999
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d0_c2_0d0_8TeV_m0.root", 19300, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d0_c2_0d0_8TeV_m0.root", 17599 , 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d0_c2_0d0_8TeV_m0.root", 17599,  0, 1]) #17399

signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d25_c2_0d0_8TeV_m0.root", 19797, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d25_c2_0d0_8TeV_m0.root", 19799, 3, 1])#19699
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d25_c2_0d0_8TeV_m0.root", 19099, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d25_c2_0d0_8TeV_m0.root", 19597, 3, 1])#19397
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d25_c2_0d0_8TeV_m0.root", 14299, 3, 1])#19397
#signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV.root", 20000, 1, 0.0000212])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d25_c2_0d0_8TeV_m0.root", 19197, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d25_c2_0d0_8TeV_m0.root", 19299 , 3, 1])#19099
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d25_c2_0d0_8TeV_m0.root", 19098, 3, 1])

#nstrategies = 25
## For the efficiency plot
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_0d75_c2_m3_8TeV_m0.root", 17398, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_0d75_c2_m3_8TeV_m0.root", 16398, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_0d75_c2_m3_8TeV_m0.root", 15699, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_0d75_c2_m3_8TeV_m0.root", 15699, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_0d75_c2_m3_8TeV_m0.root", 19597, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_0d75_c2_m3_8TeV_m0.root", 20000, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_0d75_c2_m3_8TeV_m0.root", 8199, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_0d75_c2_m3_8TeV_m0.root", 16899, 2, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d0_c2_m3_8TeV_m0.root", 15699, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d0_c2_m3_8TeV_m0.root", 16197, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d0_c2_m3_8TeV_m0.root", 13699, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d0_c2_m3_8TeV_m0.root", 19900, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d0_c2_m3_8TeV_m0.root", 16597, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d0_c2_m3_8TeV_m0.root", 17197, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d0_c2_m3_8TeV_m0.root", 15496, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d0_c2_m3_8TeV_m0.root", 17899, 0, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d25_c2_m3_8TeV_m0.root", 15398, 3,1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d25_c2_m3_8TeV_m0.root", 14799, 3,1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d25_c2_m3_8TeV_m0.root", 19800, 3,1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d25_c2_m3_8TeV_m0.root", 18897, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d25_c2_m3_8TeV_m0.root", 16999, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d25_c2_m3_8TeV_m0.root", 19900, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d25_c2_m3_8TeV_m0.root", 16298, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d25_c2_m3_8TeV_m0.root", 39600, 3,1])
#previous 17698
#nstrategies = 19
## For the efficiency plot
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m20_Yt_0d75_c2_m3_8TeV_m0.root", 17398, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m15_Yt_0d75_c2_m3_8TeV_m0.root", 16398, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m10_Yt_0d75_c2_m3_8TeV_m0.root", 15699, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_0d0_Yt_0d75_c2_m3_8TeV_m0.root", 15699, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_10_Yt_0d75_c2_m3_8TeV_m0.root", 20000, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_15_Yt_0d75_c2_m3_8TeV_m0.root", 8199, 2, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_20_Yt_0d75_c2_m3_8TeV_m0.root", 16899, 2, 1])

#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m20_Yt_1d0_c2_m3_8TeV_m0.root", 15699, 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m15_Yt_1d0_c2_m3_8TeV_m0.root", 16197, 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m10_Yt_1d0_c2_m3_8TeV_m0.root", 13699, 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_0d0_Yt_1d0_c2_m3_8TeV_m0.root", 19900, 0, 1])
##signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_10_Yt_1d0_c2_m3_8TeV_m0.root", 17197, 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_15_Yt_1d0_c2_m3_8TeV_m0.root", 15496, 0, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_20_Yt_1d0_c2_m3_8TeV_m0.root", 17899, 0, 1])

#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m20_Yt_1d25_c2_m3_8TeV_m0.root", 15398, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m15_Yt_1d25_c2_m3_8TeV_m0.root", 14799, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_m10_Yt_1d25_c2_m3_8TeV_m0.root", 19800, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_10_Yt_1d25_c2_m3_8TeV_m0.root", 19900, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_15_Yt_1d25_c2_m3_8TeV_m0.root", 16298, 3, 1])
#signals.append(["v35_fitToMgg_nonresSearch_noKinFit/ggHH_Lam_20_Yt_1d25_c2_m3_8TeV_m0.root", 17698, 3, 1])

#nstrategies = 25
## For the efficiency plot
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_0d75_c2_m2_8TeV_m0.root", 9698, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_0d75_c2_m2_8TeV_m0.root", 16399, 2,1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_0d75_c2_m2_8TeV_m0.root", 20000, 2,1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_0d75_c2_m2_8TeV_m0.root", 16598, 2,1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_0d75_c2_m2_8TeV_m0.root", 14198, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_0d75_c2_m2_8TeV_m0.root", 20000, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_0d75_c2_m2_8TeV_m0.root", 15599, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_0d75_c2_m2_8TeV_m0.root", 7197, 2, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d0_c2_m2_8TeV_m0.root", 18499, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d0_c2_m2_8TeV_m0.root", 20000, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d0_c2_m2_8TeV_m0.root", 18899, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d0_c2_m2_8TeV_m0.root", 18498, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d0_c2_m2_8TeV_m0.root", 17896, 0,1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d0_c2_m2_8TeV_m0.root", 17197, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d0_c2_m2_8TeV_m0.root", 17795, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d0_c2_m2_8TeV_m0.root", 16897, 0, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d25_c2_m2_8TeV_m0.root", 18499, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d25_c2_m2_8TeV_m0.root", 18499, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d25_c2_m2_8TeV_m0.root", 20000, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d25_c2_m2_8TeV_m0.root", 19298, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d25_c2_m2_8TeV_m0.root", 19899, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d25_c2_m2_8TeV_m0.root", 20000, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d25_c2_m2_8TeV_m0.root", 20000, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d25_c2_m2_8TeV_m0.root", 17698, 3, 1])
#nstrategies = 25
## For the efficiency plot
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_0d75_c2_2_8TeV_m0.root", 20000, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_0d75_c2_2_8TeV_m0.root", 18699, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_0d75_c2_2_8TeV_m0.root", 17499, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_0d75_c2_2_8TeV_m0.root", 17199, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_0d75_c2_2_8TeV_m0.root", 16698, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_0d75_c2_2_8TeV_m0.root", 17499, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_0d75_c2_2_8TeV_m0.root", 18699, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_0d75_c2_2_8TeV_m0.root", 20000, 2, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d0_c2_2_8TeV_m0.root", 16598, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d0_c2_2_8TeV_m0.root", 18797, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d0_c2_2_8TeV_m0.root", 8598, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d0_c2_2_8TeV_m0.root", 15496, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d0_c2_2_8TeV_m0.root", 11598, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d0_c2_2_8TeV_m0.root", 7799, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d0_c2_2_8TeV_m0.root", 9698, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d0_c2_2_8TeV_m0.root", 15999, 0, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d25_c2_2_8TeV_m0.root", 16397, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d25_c2_2_8TeV_m0.root", 16295, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d25_c2_2_8TeV_m0.root", 16698, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d25_c2_2_8TeV_m0.root", 20000, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d25_c2_2_8TeV_m0.root", 16698, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d25_c2_2_8TeV_m0.root", 7093, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d25_c2_2_8TeV_m0.root", 13699, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d25_c2_2_8TeV_m0.root", 16799, 3, 1])
#nstrategies = 25
# #For the efficiency plot
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_0d75_c2_3_8TeV_m0.root", 19800, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_0d75_c2_3_8TeV_m0.root", 16498, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_0d75_c2_3_8TeV_m0.root", 19800, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_0d75_c2_3_8TeV_m0.root", 10396, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_0d75_c2_3_8TeV_m0.root", 7497, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_0d75_c2_3_8TeV_m0.root", 14398, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_0d75_c2_3_8TeV_m0.root", 16398, 2, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_0d75_c2_3_8TeV_m0.root", 20000, 2, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d0_c2_3_8TeV_m0.root", 14198, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d0_c2_3_8TeV_m0.root", 14297, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d0_c2_3_8TeV_m0.root", 14096, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d0_c2_3_8TeV_m0.root", 11898, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d0_c2_3_8TeV_m0.root", 20000, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d0_c2_3_8TeV_m0.root", 13897, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d0_c2_3_8TeV_m0.root", 15299, 0, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d0_c2_3_8TeV_m0.root", 15999, 0, 1])

#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d25_c2_3_8TeV_m0.root", 19499, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d25_c2_3_8TeV_m0.root", 16499, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d25_c2_3_8TeV_m0.root", 14499, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d25_c2_3_8TeV_m0.root", 18995, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d25_c2_3_8TeV_m0.root", 17498, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d25_c2_3_8TeV_m0.root", 19900, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d25_c2_3_8TeV_m0.root", 7699, 3, 1])
#signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d25_c2_3_8TeV_m0.root", 15099, 3, 1])

#nstrategies = 23
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m20_Yt_1d0_c2_0d0_8TeV.root", 18997, 0, 1])#18797
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m15_Yt_1d0_c2_0d0_8TeV.root", 17800, 0, 1])#17700
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_m10_Yt_1d0_c2_0d0_8TeV.root", 17798, 0, 1])
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_0d0_Yt_1d0_c2_0d0_8TeV.root", 19800, 0, 1])#19600
#signals.append(["v34_alpha_fitToMgg_noKinFit_300_massCutVersion_4/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV.root", 20000, 0, 1])
####not use anymore in the final version of the papaer:
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/ggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV_m0.root", 20000, 1, 1])
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
lumi = 19712.
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
#the following ### commented lines are not in use, only for reference
#	for icut in range(4):
#		for ievt in xrange(chain.GetEntries()):
#			chain.GetEntry(ievt)
#			if icut == 0 and (chain.cut_based_ct == 0 or chain.cut_based_ct == 1):
#				totw += chain.evWeight
##				totw += 1
#			elif icut == 1 and chain.cut_based_ct == 0:
#				totw += chain.evWeight
##				totw += 1
#			if icut == 2 and (chain.cut_based_ct == 2 or chain.cut_based_ct == 3):
#				totw += chain.evWeight
##				totw += 1
#			elif icut == 3 and chain.cut_based_ct == 2:
#				totw += chain.evWeight
##				totw += 1

	for icut in range(2):
		totw = 0
		for ievt in xrange(chain.GetEntries()):
			chain.GetEntry(ievt)
			if icut == 0 and (chain.cut_based_ct == 0 or chain.cut_based_ct ==1):
				totw += chain.evWeight
#				totw += 1
			elif icut == 1 and chain.cut_based_ct == 0:
				totw += chain.evWeight
#				totw += 1
# totw *= 1000. # due to hand-put 1000 factor to weight for limit settings numerical precision
		print totw, nprocessed, lumi
#		eff = float(totw )/ nprocessed / float(sigma)
		eff = totw * 1000./ lumi  / sigma
#		eff = totw / nprocessed / sigma
		deff = (eff * (1 - eff) / nprocessed)**0.5
		print eff*100, deff*100
#		eff = totw / lumi  / sigma
#		deff = 0.05

# print signal, round(eff * 100, 2), round(deff * 100,2)
# print signal.split("_")[1].split("m")[3], "&", "$", round(eff * 100, 2), "$", "&", "$\pm", round(deff * 100,2), "$"
		lambda_hhh = 1.0
		try:
			lambda_hhh = float(signal.split("ggHH_Lam_")[1].split("_8TeV_m0.root")[0].split("_Yt")[0])
			print "lambda_hhh= ", lambda_hhh, "&", "$", round(eff * 100, 5), "$", "&", "$\pm", round(deff * 100, 5), "$"
		except ValueError:
			try:
				lambda_hhh = float(signal.split("ggHH_Lam_")[1].split("_8TeV_m0.root")[0].split("d0_Yt")[0])
				print "lambda_hhh= ", lambda_hhh, "&", "$", round(eff * 100, 5), "$", "&", "$\pm", round(deff * 100, 5), "$"			
			except ValueError:
				try:
					lambda_hhh = -float(signal.split("ggHH_Lam_m")[1].split("_8TeV_m0.root")[0].split("_Yt")[0])
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
#	ey[istrategy] = ex[istrategy]
#gROOT.ProcessLine(".x CMSStyle.C")
gROOT.ProcessLine(".x setTDRStyle.C")
c1 = TCanvas()
legend = TLegend(0.35, 0.75, 0.5, 0.92, " ")
legend.SetTextFont(42)
legend1= TLegend(0.5, 0.75, 0.8, 0.92)
legend1.SetTextFont(42)
legend1.SetHeader("gg #rightarrow HH #rightarrow #gamma#gammab#bar{b},    #it{c_{2}}=0.")

legend.SetTextSize(0.028)
legend.SetFillColor(ROOT.kWhite)
legend.SetLineColor(ROOT.kWhite)
legend.SetShadowColor(ROOT.kWhite)
legend1.SetTextSize(0.028)
legend1.SetFillColor(ROOT.kWhite)
legend1.SetLineColor(ROOT.kWhite)
legend1.SetShadowColor(ROOT.kWhite)

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



gr0.SetMarkerColor(ROOT.kGreen+2)
gr1.SetMarkerColor(ROOT.kGreen+2)
gr2.SetMarkerColor(ROOT.kBlack)
gr3.SetMarkerColor(ROOT.kBlack)
gr4.SetMarkerColor(ROOT.kRed+2)
gr5.SetMarkerColor(ROOT.kRed+2)
gr6.SetMarkerColor(ROOT.kBlue+2)
gr7.SetMarkerColor(ROOT.kBlue+2)

gr0.SetLineColor(ROOT.kGreen+2)
gr1.SetLineColor(ROOT.kGreen+2)
gr2.SetLineColor(ROOT.kBlack)
gr3.SetLineColor(ROOT.kBlack)
gr4.SetLineColor(ROOT.kRed+2)
gr5.SetLineColor(ROOT.kRed+2)
gr6.SetLineColor(ROOT.kBlue+2)
gr7.SetLineColor(ROOT.kBlue+2)

gr0.SetLineStyle(1)
gr1.SetLineStyle(2)
gr2.SetLineStyle(1)
gr3.SetLineStyle(2)
gr4.SetLineStyle(1)
gr5.SetLineStyle(2)
gr6.SetLineStyle(1)
gr7.SetLineStyle(2)

gr0.SetLineWidth(2)
gr1.SetLineWidth(2)
gr2.SetLineWidth(2)
gr3.SetLineWidth(2)
gr4.SetLineWidth(2)
gr5.SetLineWidth(2)
gr6.SetLineWidth(2)
gr7.SetLineWidth(2)

gr0.SetMarkerStyle(23)
gr1.SetMarkerStyle(29)
gr2.SetMarkerStyle(20)
gr3.SetMarkerStyle(29)
gr4.SetMarkerStyle(24)
gr5.SetMarkerStyle(21)
gr6.SetMarkerStyle(25)
gr7.SetMarkerStyle(22)

gr0.SetMarkerSize(0.8)
gr1.SetMarkerSize(1.2)
gr2.SetMarkerSize(1)
gr3.SetMarkerSize(2)
gr4.SetMarkerSize(0.8)
gr5.SetMarkerSize(0.8)
gr6.SetMarkerSize(0.8)
gr7.SetMarkerSize(1.)

gr0.SetTitle("")
#gr0.GetXaxis().SetTitle("#lambda_{hhh}/#lambda_{SM}")
gr0.GetXaxis().SetTitle("#kappa_{#lambda}")
gr0.GetYaxis().SetTitle("Signal selection efficiency (%)")
print min(numpy.amin(y[0]), numpy.amin(y[1])), max(numpy.amax(y[0]), numpy.amax(y[1]))
print min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1]))
gr0.SetMinimum(0.)
gr0.SetMaximum(35.)
gr0.GetXaxis().SetLimits(-25., 25.)
#gr0.SetMaximum(max(numpy.amax(y[0]), numpy.amax(y[1])))
#gr0.SetMinimum(min(numpy.amin(y[0]), numpy.amin(y[1])))
#gr0.GetXaxis().SetLimits(min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1])))
gr0.Draw("ALPE1");
gr1.Draw("LPE1");
gr4.Draw("PE1");
gr5.Draw("PE1");
gr6.Draw("PE1");
gr7.Draw("PE1");
#gr2.Draw("P");
#gr3.Draw("P");
legend.SetEntrySeparation(0.5)
legend1.SetEntrySeparation(0.5)
#legend.SetNColumns(2)
#legend.SetColumnSeparation(2.)
#legend.SetHeader("#splitline{gg #rightarrow HH #rightarrow #gamma#gammab#bar{b}}{c_{2}=0.}")
#legend.SetHeader("#splitline{HH #rightarrow #gamma#gammab#bar{b}}{High-Purity(H):m_{HH}<350.GeV      Medium(M):m_{HH}<360.GeV }")
#legend.AddEntry(gr0.GetName(), "High+Medium Purity m_{hh}>360.GeV y_{t}=1. c_{2}=0.", "lp")
#legend.AddEntry(gr1.GetName(), "High Purity m_{hh}>350.GeV y_{t}=1. c_{2}=0.", "lp")
##legend.AddEntry(gr2.GetName(), "High+Medium Purity - ggHH", "lp")
##legend.AddEntry(gr3.GetName(), "High Purity - ggHH", "lp")
#legend.AddEntry(gr4.GetName(), "High+Medium Purity m_{hh}>360.GeV y_{t}=0.75 c_{2}=0.", "lp")
#legend.AddEntry(gr5.GetName(), "High Purity m_{hh}>350.GeV y_{t}=0.75 c_{2}=0.", "lp")
#legend.AddEntry(gr6.GetName(), "High+Medium Purity m_{hh}>360.GeV y_{t}=1.25 c_{2}=0.", "lp")
#legend.AddEntry(gr7.GetName(), "High Purity m_{hh}>350.GeV y_{t}=1.25 c_{2}=0.", "lp")

legend.AddEntry(gr0.GetName(), "All cat.", "lp")
legend1.AddEntry(gr1.GetName(), "High-Purity cat.     #kappa_{t}=1.", "lp")
#legend.AddEntry(gr2.GetName(), "All cat.  - SM", "p")
#legend.AddEntry(gr3.GetName(), "High-Purity cat. - SM", "p")
legend.AddEntry(gr4.GetName(), "All cat.", "p")
legend1.AddEntry(gr5.GetName(), "High-Purity cat.     #kappa_{t}=0.75", "p")
legend.AddEntry(gr6.GetName(), "All cat.", "p")
legend1.AddEntry(gr7.GetName(), "High-Purity cat.     #kappa_{t}=1.25", "p")


#legend.AddEntry(gr0.GetName(), "High+Medium Purity -Lm-Yt_1_c2_3_v44", "lp")
#legend.AddEntry(gr1.GetName(), "High Purity -Lm- Yt_1", "lp")
##legend.AddEntry(gr2.GetName(), "High+Medium Purity - ggHH", "lp")
##legend.AddEntry(gr3.GetName(), "High Purity - ggHH", "lp")
#legend.AddEntry(gr4.GetName(), "High+Medium Purity -Lm-Yt_0d75_c2_3_v44", "lp")
#legend.AddEntry(gr5.GetName(), "High Purity Lm--Yt_0d75", "lp")
#legend.AddEntry(gr6.GetName(), "High+Medium Purity -Lm-Yt_1d25_c2_3_v44", "lp")
#legend.AddEntry(gr7.GetName(), "High Purity Lm--Yt_1d25", "lp")

legend.Draw()
legend1.Draw()
latexLabel = TLatex()
latexLabel.SetTextSize(0.75 * c1.GetTopMargin())
#latexLabel.SetTextSize(0.03)
latexLabel.SetNDC()
#latexLabel.SetTextFont(42) # helvetica
#latexLabel.DrawLatex(0.45, 0.75, "Low-mass")
#latexLabel.DrawLatex(0.45, 0.70, "c_{2}=0.")

latexLabel.SetTextFont(42) # helvetica
latexLabel.DrawLatex(0.86, 0.96, "8 TeV")
latexLabel.SetTextFont(62) # helvetica bold face
latexLabel.DrawLatex(0.16, 0.88, "CMS")
latexLabel.SetTextFont(52) # helvetica italics
latexLabel.DrawLatex(0.16, 0.84, "Simulation")
latexLabel.SetTextFont(42) # helvetica 
latexLabel.SetTextColor(ROOT.kBlue+2)
latexLabel.DrawLatex(0.75, 0.65,"High-mass")


#latexLabel.DrawLatex(0.15, 0.96, "CMS Simulation")
#latexLabel.DrawLatex(0.80, 0.96, "#sqrt{s} = 8 TeV")
#c1.Print("2015-01-14_lambda_hhh_eff_Summary_Hmass0_massCutVersion_4_c2_3_v44_withKinFit_bis.pdf")
#c1.Print("2015-01-14_lambda_hhh_eff_Summary_Hmass0_massCutVersion_4_c2_3_v44_withKinFit_bis.png")
#c1.Print("2015-01-14_lambda_hhh_eff_Summary_Hmass0_massCutVersion_4_c2_3_v44_withKinFit_bis.root")
#c1.Print("2015-01-14_lambda_hhh_eff_Summary_Hmass0_massCutVersion_4_c2_3_v44_withKinFit_bis.eps")
c1.Print("2015-11-12_lambda_hhh_eff_Summary_Hmass0_massCutVersion_4_c2_0d0_v44_withKinFit_bis.pdf")
c1.Print("2015-11-12_lambda_hhh_eff_Summary_Hmass0_massCutVersion_4_c2_0d0_v44_withKinFit_bis.png")
#c1.Print("2015-11-12_lambda_hhh_eff_Summary_Lmass0_massCutVersion_4_c2_0d0_v44_withKinFit_bis.root")
#c1.Print("2015-11-12_lambda_hhh_eff_Summary_Lmass0_massCutVersion_4_c2_0d0_v44_withKinFit_bis.eps")

#c1.Print("eff_cat0.pdf")
#c1.Print("eff_cat0.png")
#c1.Print("eff_cat0.root")
#c1.Print("eff_cat1.pdf")
#c1.Print("eff_cat1.png")
#c1.Print("eff_cat1.root")
