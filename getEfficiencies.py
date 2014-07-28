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
# For the efficiency plot
signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_m260.root", 96880, 1, 0.0000212])
signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_0hhh_m260.root", 96880, 0, 0.0000212])
signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_1hhh_m260.root", 96880, 0, 0.0000212])
signals.append(["v33_ggHH_fitToMgg_noKinFit/ggHH_8TeV_2hhh_m260.root", 96880, 0, 0.0000212])

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


#print signals
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
            if icut == 0:
                totw += chain.evWeight
            elif icut == 1 and chain.cut_based_ct == 0:
                totw += chain.evWeight
#        totw *= 1000. # due to hand-put 1000 factor to weight for limit settings numerical precision
        eff = totw / lumi / sigma
        deff = (eff * (1 - eff) / nprocessed)**0.5
    #    print signal, round(eff * 100, 2), round(deff * 100,2)
    #    print signal.split("_")[1].split("m")[3], "&", "$", round(eff * 100, 2), "$", "&", "$\pm", round(deff * 100,2), "$"
        lambda_hhh = 1.0
        try:
            lambda_hhh = float(signal.split("ggHH_8TeV_")[1].split("_m260.root")[0].split("hhh")[0])
        except ValueError:
            lambda_hhh = 1.0
        print "lambda_hhh= ", lambda_hhh, "&", "$", round(eff * 100, 5), "$", "&", "$\pm", round(deff * 100, 5), "$"

        n[istrategy+icut*2] = n[istrategy+icut*2]+1
        x[istrategy+icut*2].append( lambda_hhh )
        y[istrategy+icut*2].append(round(eff * 100, 2))
        ex[istrategy+icut*2].append(0.)
        ey[istrategy+icut*2].append(round(deff * 100,2))
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
legend = TLegend(0.15, 0.74, 0.50, 0.94, "HH #rightarrow #gamma#gammab#bar{b}")
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
gr0.SetTitle("")
gr0.GetXaxis().SetTitle("#lambda_{hhh}")
gr0.GetYaxis().SetTitle("Signal selection efficiency (%)")
print min(numpy.amin(y[0]), numpy.amin(y[1])), max(numpy.amax(y[0]), numpy.amax(y[1]))
print min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1]))
gr0.SetMinimum(0.)
gr0.SetMaximum(50.)
gr0.GetXaxis().SetLimits(-1., 3.)
#gr0.SetMaximum(max(numpy.amax(y[0]), numpy.amax(y[1])))
#gr0.SetMinimum(min(numpy.amin(y[0]), numpy.amin(y[1])))
#gr0.GetXaxis().SetLimits(min(numpy.amin(x[0]), numpy.amin(x[1])), max(numpy.amax(x[0]), numpy.amax(x[1])))
gr0.Draw("ALP");
gr1.Draw("LP");
gr2.Draw("LP");
gr3.Draw("LP");
legend.AddEntry(gr0.GetName(), "High+Medium Purity - #lambda-reweighted", "lp")
legend.AddEntry(gr1.GetName(), "High+Medium Purity - ggHH", "lp")
legend.AddEntry(gr2.GetName(), "High Purity - #lambda-reweighted", "lp")
legend.AddEntry(gr3.GetName(), "High Purity - ggHH", "lp")
legend.Draw()
latexLabel = TLatex()
latexLabel.SetTextSize(0.03)
latexLabel.SetNDC()
latexLabel.DrawLatex(0.15, 0.96, "CMS Work in progress")
latexLabel.DrawLatex(0.80, 0.96, "#sqrt{s} = 8 TeV")
c1.Print("2014-07-23_lambda_hhh_eff.pdf")
c1.Print("2014-07-23_lambda_hhh_eff.png")
c1.Print("2014-07-23_lambda_hhh_eff.root")
#c1.Print("eff_cat0.pdf")
#c1.Print("eff_cat0.png")
#c1.Print("eff_cat0.root")
#c1.Print("eff_cat1.pdf")
#c1.Print("eff_cat1.png")
#c1.Print("eff_cat1.root")
