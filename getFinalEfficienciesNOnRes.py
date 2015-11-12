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
nstrategies = 19
 ##For the efficiency plot
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_0d75_c2_0d0_8TeV_m0.root", 16498, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_0d75_c2_0d0_8TeV_m0.root", 19800, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_0d75_c2_0d0_8TeV_m0.root", 10396, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_15_Yt_0d75_c2_0d0_8TeV_m0.root", 16398, 2, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_20_Yt_0d75_c2_0d0_8TeV_m0.root", 20000, 2, 1])

signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d0_c2_0d0_8TeV_m0.root", 14198, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d0_c2_0d0_8TeV_m0.root", 14297, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d0_c2_0d0_8TeV_m0.root", 14096, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_0d0_Yt_1d0_c2_0d0_8TeV_m0.root", 11898, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d0_c2_0d0_8TeV_m0.root", 14497, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d0_c2_0d0_8TeV_m0.root", 15299, 0, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d0_c2_0d0_8TeV_m0.root", 15999, 0, 1])

signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m20_Yt_1d25_c2_0d0_8TeV_m0.root", 19499, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m15_Yt_1d25_c2_0d0_8TeV_m0.root", 16499, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_m10_Yt_1d25_c2_0d0_8TeV_m0.root", 14499, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_10_Yt_1d25_c2_0d0_8TeV_m0.root", 19900, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_15_Yt_1d25_c2_0d0_8TeV_m0.root", 7699, 3, 1])
signals.append(["/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitToMgg_nonresSearch_withKinFit/ggHH_Lam_20_Yt_1d25_c2_0d0_8TeV_m0.root", 15099, 3, 1])


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
		eff = totw * 1000./ lumi  / sigma
#		eff = totw / nprocessed / sigma
		deff = (eff * (1 - eff) / nprocessed)**0.5
		print eff*100, deff*100
#		eff = totw / lumi  / sigma
#		deff = 0.05

# print signal, round(eff * 100, 2), round(deff * 100,2)
		print signal
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
	print "in loop for", istrategy
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
gr0.SetMaximum(40.)
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
3
legend.AddEntry(gr0.GetName(), "High+Medium Purity - Yt_1_c2_0d0_v44", "lp")
legend.AddEntry(gr1.GetName(), "High Purity - Yt_1", "lp")
#legend.AddEntry(gr2.GetName(), "High+Medium Purity - ggHH", "lp")
#legend.AddEntry(gr3.GetName(), "High Purity - ggHH", "lp")
legend.AddEntry(gr4.GetName(), "High+Medium Purity - Yt_0d75_c2_0d0_v44", "lp")
legend.AddEntry(gr5.GetName(), "High Purity -Yt_0d75", "lp")
legend.AddEntry(gr6.GetName(), "High+Medium Purity -Yt_1d25_c2_0d0_v44", "lp")
legend.AddEntry(gr7.GetName(), "High Purity -Yt_1d25", "lp")

legend.Draw()
latexLabel = TLatex()
latexLabel.SetTextSize(0.03)
latexLabel.SetNDC()
latexLabel.DrawLatex(0.15, 0.96, "CMS Work in progress")
latexLabel.DrawLatex(0.80, 0.96, "#sqrt{s} = 8 TeV")
c1.Print("2014-11-10_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_0d0_v44_withKinFit.pdf")
c1.Print("2014-11-10_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_0d0_v44_withKinFit.png")
c1.Print("2014-11-10_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_0d0_v44_withKinFit.root")
c1.Print("2014-11-10_lambda_hhh_eff_Summary_mass0_massCutVersion_4_c2_0d0_v44_withKinFit.eps")
#c1.Print("eff_cat0.pdf")
#c1.Print("eff_cat0.png")
#c1.Print("eff_cat0.root")
#c1.Print("eff_cat1.pdf")
#c1.Print("eff_cat1.png")
#c1.Print("eff_cat1.root")
