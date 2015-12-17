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
nstrategies = 4
doSpin2 =  False
if doSpin2:
  nstrategies = 12
# For the efficiency plot
# Note: the number of processed events is used only to compute the error on the efficiencies, the rest of the normalization should have been already included correctly at the time of the processing of the samples in h2gglobe
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/MSSM_m260_8TeV_m260.root", 300000, 0])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/Radion_m270_8TeV_m270.root", 19996, 0])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/Radion_m300_8TeV_m300.root", 19972, 0])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/Radion_m350_8TeV_m350.root", 18498, 0])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/Radion_m400_8TeV_m400.root", 19697, 0])
###signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/Radion_m450_8TeV_m450.root", 19999, 0])
###signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/Radion_m500_8TeV_m500.root", 19970, 0])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m400_8TeV_m400.root", 19697, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m450_8TeV_m450.root", 19999, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m500_8TeV_m500.root", 19970, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m550_8TeV_m550.root", 19995, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m600_8TeV_m600.root", 18197, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m650_8TeV_m650.root", 20000, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m700_8TeV_m700.root", 19969, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m800_8TeV_m800.root", 19999, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m900_8TeV_m900.root", 19996, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m1000_8TeV_m1000.root", 19951, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMggjj_withKinFit/Radion_m1100_8TeV_m1100.root", 19400, 1])
# graviton samples - non narrow, 4 samples
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitTo2D_resSearch_withRegKinFit/Graviton_m300_8TeV_m300.root", 49941, 4])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitTo2D_resSearch_withRegKinFit/Graviton_m500_8TeV_m500.root", 49905, 4])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitToMggjj_withKinFit/Graviton_m500_8TeV_m500.root", 49905, 5])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitToMggjj_withKinFit/Graviton_m700_8TeV_m700.root", 49911, 5])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitToMggjj_withKinFit/Graviton_m1000_8TeV_m1000.root", 49921, 5])
# graviton narrow width - narrow, 5 samples
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitTo2D_resSearch_withRegKinFit/Graviton_m270_LR3tev_8TeV_m270.root", 19798, 8])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitTo2D_resSearch_withRegKinFit/Graviton_m300_LR3tev_8TeV_m300.root", 19996, 8])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitTo2D_resSearch_withRegKinFit/Graviton_m350_LR3tev_8TeV_m350.root", 19999, 8])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitTo2D_resSearch_withRegKinFit/Graviton_m450_LR3tev_8TeV_m450.root", 19999, 8])
signals.append(["/afs/cern.ch/work/o/obondu/Higgs/CMSSW_6_1_1_radion_v2/src/Selection/v44_fitToMggjj_withKinFit/Graviton_m700_LR3tev_8TeV_m700.root", 19999, 9])

# mssm samples
#signals.append(["v26_fitToMgg_noKinFit/Graviton_m300_8TeV_m300.root", 49941])
#signals.append(["v26_fitToMgg_noKinFit/Graviton_m500_8TeV_m500.root", 49905])
#signals.append(["", ])


#print signals
lumi = 19706.

n = []
x = []
y = []
ex = []
ey = []

for i in range(12):
    n.append( 0 )
    x.append( [] )
    y.append( [] )
    ex.append( [] )
    ey.append( [] )


for signal, nprocessed, istrategy in signals:
    if (not doSpin2) and (istrategy >= nstrategies):
      continue
    chain = TChain("TCVARS")
    chain.Add(signal)
    for icut in range(2):
        totw = 0
        for ievt in xrange(chain.GetEntries()):
            chain.GetEntry(ievt)
            if icut == 0:
                totw += chain.evWeight
            elif icut == 1 and chain.cut_based_ct == 0:
                totw += chain.evWeight
        totw *= 1000. # due to hand-put 1000 factor to weight for limit settings numerical precision
        eff = totw / lumi
        deff = (eff * (1 - eff) / nprocessed)**0.5
    #    print signal, round(eff * 100, 2), round(deff * 100,2)
        print signal.split("_")[-1].split(".")[0].split("m")[1], "&", "$", round(eff * 100, 2), "$", "&", "$\pm", round(deff * 100,2), "$"
        n[istrategy+icut*2] = n[istrategy+icut*2]+1
        x[istrategy+icut*2].append(signal.split("_")[-1].split(".")[0].split("m")[1])
        y[istrategy+icut*2].append(round(eff * 100, 2))
        ex[istrategy+icut*2].append(0.)
        if doSpin2 : 
          ey[istrategy+icut*2].append(round(deff * 500,2))
        else :
          ey[istrategy+icut*2].append(round(deff * 100,2))
    del chain

for istrategy in range(nstrategies):
    x[istrategy] = numpy.asarray(x[istrategy], dtype='f')
    y[istrategy] = numpy.asarray(y[istrategy], dtype='f')
    ex[istrategy] = numpy.asarray(ex[istrategy], dtype='f')
    ey[istrategy] = numpy.asarray(ey[istrategy], dtype='f')
    # Removing y errors
#    ey[istrategy] = ex[istrategy] 

if not doSpin2:
# fill with dumb stuff if no sample is present
  for istrategy in range(nstrategies, 12):
    x[istrategy] = x[0]
    y[istrategy] = y[0]
    ex[istrategy] = ex[0]
    ey[istrategy] = ey[0]



#gROOT.ProcessLine(".x CMSStyle.C")
gROOT.ProcessLine(".x setTDRStyle.C")
c1 = TCanvas()
legend = TLegend(0.35, 0.70, 0.90, 0.91)
legend.SetTextFont(42)
legend.SetHeader("gg #rightarrow X #rightarrow HH #rightarrow #gamma#gammab#bar{b}")
if doSpin2:
  legend.SetY1(0.80)
#lower[0.]{High+Medium Purity}}")
#legend = TLegend(0.15, 0.72, 0.50, 0.92, "#splitline{gg#rightarrowX#rightarrow HH #rightarrow b#bar{b}#gamma#gamma}{High-purity category}")
#legend = TLegend(0.15, 0.72, 0.50, 0.92, "#splitline{gg#rightarrowX#rightarrow HH #rightarrow b#bar{b}#gamma#gamma}{Medium-purity category}")
legend.SetTextSize(0.028)
legend.SetFillStyle(0)
legend.SetFillColor(ROOT.kWhite)
legend.SetLineColor(ROOT.kWhite)
legend.SetShadowColor(ROOT.kWhite)
gr0 = TGraphErrors(n[0],x[0],y[0],ex[0],ey[0]);
gr1 = TGraphErrors(n[1],x[1],y[1],ex[1],ey[1]);
gr2 = TGraphErrors(n[2],x[2],y[2],ex[2],ey[2]);
gr3 = TGraphErrors(n[3],x[3],y[3],ex[3],ey[3]);
gr0.SetName("gr0")
gr1.SetName("gr1")
gr2.SetName("gr2")
gr3.SetName("gr3")
gr0.SetMarkerColor(ROOT.kGreen+3)
gr1.SetMarkerColor(ROOT.kBlue+2)
gr2.SetMarkerColor(ROOT.kGreen+3)
gr3.SetMarkerColor(ROOT.kBlue+2)
gr0.SetLineColor(ROOT.kGreen+3)
gr1.SetLineColor(ROOT.kBlue+2)
gr2.SetLineColor(ROOT.kGreen+3)
gr3.SetLineColor(ROOT.kBlue+2)
gr0.SetMarkerStyle(20)
gr1.SetMarkerStyle(21)
gr2.SetMarkerStyle(24)
gr3.SetMarkerStyle(25)
gr2.SetLineStyle(2)
gr3.SetLineStyle(2)
# graviton samples, if the case arise
print n
gr4 = TGraphErrors(n[4],x[4],y[4],ex[4],ey[4]);
gr5 = TGraphErrors(n[5],x[5],y[5],ex[5],ey[5]);
gr6 = TGraphErrors(n[6],x[6],y[6],ex[6],ey[6]);
gr7 = TGraphErrors(n[7],x[7],y[7],ex[7],ey[7]);
gr4.SetName("gr4")
gr5.SetName("gr5")
gr6.SetName("gr6")
gr7.SetName("gr7")
gr4.SetMarkerColor(ROOT.kMagenta+3)
gr5.SetMarkerColor(ROOT.kCyan+2)
gr6.SetMarkerColor(ROOT.kMagenta+3)
gr7.SetMarkerColor(ROOT.kCyan+2)
gr4.SetLineColor(ROOT.kMagenta+3)
gr5.SetLineColor(ROOT.kCyan+2)
gr6.SetLineColor(ROOT.kMagenta+3)
gr7.SetLineColor(ROOT.kCyan+2)
gr4.SetMarkerStyle(33)
gr5.SetMarkerStyle(34)
gr6.SetMarkerStyle(27)
gr7.SetMarkerStyle(28)
gr6.SetLineStyle(2)
gr7.SetLineStyle(2)
gr8 = TGraphErrors(n[8],x[8],y[8],ex[8],ey[8]);
gr9 = TGraphErrors(n[9],x[9],y[9],ex[9],ey[9]);
gr10 = TGraphErrors(n[10],x[10],y[10],ex[10],ey[10]);
gr11 = TGraphErrors(n[11],x[11],y[11],ex[11],ey[11]);
gr8.SetName("gr8")
gr9.SetName("gr9")
gr10.SetName("gr10")
gr11.SetName("gr11")
gr8.SetMarkerColor(ROOT.kRed+3)
gr9.SetMarkerColor(ROOT.kYellow+2)
gr10.SetMarkerColor(ROOT.kRed+3)
gr11.SetMarkerColor(ROOT.kYellow+2)
gr8.SetLineColor(ROOT.kRed+3)
gr9.SetLineColor(ROOT.kYellow+2)
gr10.SetLineColor(ROOT.kRed+3)
gr11.SetLineColor(ROOT.kYellow+2)
gr8.SetMarkerStyle(22)
gr9.SetMarkerStyle(23)
gr10.SetMarkerStyle(26)
gr11.SetMarkerStyle(32)
gr10.SetLineStyle(2)
gr11.SetLineStyle(2)


gr0.SetTitle("")
gr0.GetXaxis().SetLabelSize(0.04)
gr0.GetXaxis().SetTitle("m_{X}^{spin-0} (GeV)")
gr0.GetYaxis().SetTitle("Signal selection efficiency (%)")
#print min(numpy.amin(y[0]), numpy.amin(y[1])), max(numpy.amax(y[0]), numpy.amax(y[1]))
#print min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1]))
gr0.SetMinimum(0.)
gr0.SetMaximum(60.)
gr0.GetXaxis().SetLimits(240., 1109.)
#gr0.SetMaximum(max(numpy.amax(y[0]), numpy.amax(y[1])))
#gr0.SetMinimum(min(numpy.amin(y[0]), numpy.amin(y[1])))
#gr0.GetXaxis().SetLimits(min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1])))
gr0.Draw("ALP");
gr1.Draw("LP");
gr2.Draw("LP");
gr3.Draw("LP");
if doSpin2:
 gr4.Draw("LP"); 
 gr5.Draw("LP"); 
 gr6.Draw("LP"); 
 gr7.Draw("LP"); 
 gr8.Draw("LP"); 
 gr9.Draw("LP"); 
 gr10.Draw("LP"); 
 gr11.Draw("LP"); 
if not doSpin2:
  legend.AddEntry(gr0.GetName(), "All cat. (Low-mass analysis)", "lp")
  legend.AddEntry(gr1.GetName(), "All cat. (High-mass analysis)", "lp")
  legend.AddEntry(gr2.GetName(), "High-Purity cat. (Low-mass analysis)", "lp")
  legend.AddEntry(gr3.GetName(), "High-Purity cat. (High-mass analysis)", "lp")
else:
#  legend.AddEntry(gr0.GetName(), "cat0 + cat1 - low mass radion", "lp")
#  legend.AddEntry(gr1.GetName(), "cat0 + cat1 - high mass radion", "lp")
#  legend.AddEntry(gr2.GetName(), "cat0 - low mass radion", "lp")
#  legend.AddEntry(gr3.GetName(), "cat0 - high mass radion", "lp")
  legend.AddEntry(gr4.GetName(), "cat0 + cat1 - low mass graviton", "lp")
  legend.AddEntry(gr5.GetName(), "cat0 + cat1 - high mass graviton", "lp")
  legend.AddEntry(gr6.GetName(), "cat0 - low mass graviton", "lp")
  legend.AddEntry(gr7.GetName(), "cat0 - high mass graviton", "lp")
  legend.AddEntry(gr8.GetName(), "cat0 + cat1 - low mass graviton narrow", "lp")
  legend.AddEntry(gr9.GetName(), "cat0 + cat1 - high mass graviton narrow", "lp")
  legend.AddEntry(gr10.GetName(), "cat0 - low mass graviton narrow", "lp")
  legend.AddEntry(gr11.GetName(), "cat0 - high mass graviton narrow", "lp")

legend.Draw()
latexLabel = TLatex()
latexLabel.SetTextSize(0.75 * c1.GetTopMargin())
latexLabel.SetNDC()
latexLabel.SetTextFont(42) # helvetica
latexLabel.DrawLatex(0.86, 0.96, "8 TeV")
latexLabel.SetTextFont(61) # helvetica bold face
latexLabel.DrawLatex(0.16, 0.88, "CMS")
latexLabel.SetTextFont(52) # helvetica italics
latexLabel.DrawLatex(0.16, 0.84, "Simulation")
if (doSpin2) :
  latexLabel.DrawLatex(0.37, 0.20, "Limit trees v44 ; Uncertainties #times 5")

c1.Print("eff.pdf")
c1.Print("eff.png")
c1.Print("eff.root")
#c1.Print("eff_cat0.pdf")
#c1.Print("eff_cat0.png")
#c1.Print("eff_cat0.root")
#c1.Print("eff_cat1.pdf")
#c1.Print("eff_cat1.png")
#c1.Print("eff_cat1.root")
