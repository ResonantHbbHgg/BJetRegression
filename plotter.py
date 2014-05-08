#!/usr/bin/env python
#// Small dumb code to play with trees
#// O.Bondu, F. Bojarski (May 2014)
# Various python imports
from os import path
from math import log10, pow
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TLatex, TLegend, TLorentzVector
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
ROOT.TGaxis.SetMaxDigits(3);

c1 = TCanvas()
afs_plottree = "/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v10/"
eos_tree = "root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_11_tree_08"

samples = []
# samples.append([ name, dirpath, subdir, file, tree, color, style, label ])
samples.append(["Radion_m300", afs_plottree, "2014-02-17_selection_noRegression_noMassCut_v10", "Radion_m300_8TeV_noRegression_noMassCut_v10.root", "Radion_m300_8TeV", ROOT.kBlue, 0, "m_{X} = 260 GeV"])

plots = []
plots.append(["pho1_pt", "pho1_pt", "", 1., "(500, 0, 500)", "p_{T}^{#gamma} (GeV)", ""])


for name, variable, cut, norm, binning, title, additional_info in plots:
    c1 = TCanvas()
    legend = TLegend(0.45, 0.82, 0.90, 0.93, "")
    legend.SetTextSize(0.025)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    if len(samples) > 1:
        legend.SetNColumns(2)
    xnbin, xlow, xhigh = map(float, binning.strip().strip("()").split(","))
    ymax = -1
    ymin = 10000000
    firsthistname = ""
    if cut == "": cut = "1"

    for ifile, [ name, dirpath, subdir, file, tree, color, style, label ] in enumerate(samples):
#        print ifile, file, color, style, label
        chain = TChain(tree)
        chain.Add( path.join(dirpath, subdir, file) )
        sample_cut = cut
        if norm == 1.:
            sample_cut = "(" + sample_cut + ")/" + str( chain.GetEntries() )
        option = ""
        if ifile != 0:
            option = "same"
        chain.Draw(variable + ">>h_tmp" + binning, sample_cut, option)
        # Cosmetics
        h = ROOT.gDirectory.Get("h_tmp")
        h.SetName(name + str(ifile))
        if ifile == 0:
            firsthistname = name + str(ifile)
        h.SetLineWidth(3)
        h.SetLineColor(color)
        h.SetFillColor(color)
        h.SetFillStyle(style)
        h.GetXaxis().SetTitle( title )
        unit = ""
        if norm == 1.:
            if title.find("(") != -1:
                unit = title[title.find("(")+1:title.find(")")]
            h.GetYaxis().SetTitle( "Norm. to unity / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        legend.AddEntry(h.GetName(), label)
        ymax = max(ymax, h.GetMaximum())
        ymin = min(ymin, h.GetMinimum(0.0))
        del chain, h

    ymin_lin = ymin / 10.
    yrange_lin = ymax - ymin_lin
    ymax_lin = .25 * yrange_lin + ymax
    yrange_log = (log10(ymax) - log10(ymin)) / .77
    ymax_log = pow(10., .25*yrange_log + log10(ymax))
    ymin_log = pow(10., log10(ymin) - .03*yrange_log)

    latexLabel = TLatex()
    latexLabel.SetTextSize(.03)
    latexLabel.SetNDC()
    latexLabel.DrawLatex(.25, .96, "CMS Internal     L = 19.7 fb^{-1}     #sqrt{s} = 8 TeV")
    latexLabel.DrawLatex(.20, .85, additional_info)
    ROOT.gPad.RedrawAxis()
    legend.Draw()

    h = ROOT.gDirectory.Get(firsthistname)
    h.SetMaximum(ymax_lin)
    h.SetMinimum(ymin_lin)
    c1.Update()
    c1.Print("pdf/" + name + ".pdf")
    c1.Print("gif/" + name + ".gif")
    c1.Print("root/" + name + ".root")


    c1.SetLogy(1)
    h.SetMaximum(ymax_log)
    h.SetMinimum(ymin_log)
    h.GetYaxis().SetRangeUser(ymin_log, ymax_log)
    c1.Update()
    c1.Print("pdf/" + name + "_log.pdf")
    c1.Print("gif/" + name + "_log.gif")
    c1.Print("root/" + name + "_log.root")
    c1.SetLogy(0)

    del c1

