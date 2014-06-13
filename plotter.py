#!/usr/bin/env python
#// Small dumb code to play with trees
#// O.Bondu, F. Bojarski (May 2014)
# Various python imports
from os import path
from math import log10, pow
import collections
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1D, TLatex, TLegend, TLorentzVector
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
ROOT.TGaxis.SetMaxDigits(3);

c1 = TCanvas()
afs_plottree = "/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v14/2014-06-12_selection_noRegression_noMassCut_v14"
eos_tree = "root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_11_tree_08"

intL = 19706.
samples = []
# samples.append([ name, typ, dirpath, subdir, file, tree, sample_cut, color, style, label , sigma , N])
#samples.append(["Radion_m300", -300, afs_plottree, "", "Radion_m300_8TeV_noRegression_noMassCut_v*.root", "Radion_m300_8TeV", "evweight_w_btagSF", ROOT.kRed, 0, "m_{X} = 300 GeV" , 13.55e-3, 19999])
samples.append(["Data", 0, afs_plottree, "", "Data_noRegression_noMassCut_v*.root", "Data", "", ROOT.kBlack, 0, "Data", 1., 1001822])
samples.append(["qcd_40_8TeV_ff", 30, afs_plottree, "", "qcd_40_8TeV_ff_noRegression_noMassCut_v*.root", "qcd_40_8TeV_ff", "(evweight_w_btagSF) * (evweight_w_btagSF < 50)", ROOT.kCyan+2, 1001, "QCD jets" , 1., 14404429])
samples.append(["qcd_30_8TeV_ff", 30, afs_plottree, "", "qcd_30_8TeV_ff_noRegression_noMassCut_v*.root", "qcd_30_8TeV_ff", "(evweight_w_btagSF) * (evweight_w_btagSF < 50)", ROOT.kCyan+2, 1001, "QCD jets" , 1., 14404429])
samples.append(["qcd_30_8TeV_pf", 20, afs_plottree, "", "qcd_30_8TeV_pf_noRegression_noMassCut_v*.root", "qcd_30_8TeV_pf", "(evweight_w_btagSF) * (evweight_w_btagSF < 50)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
samples.append(["qcd_40_8TeV_pf", 20, afs_plottree, "", "qcd_40_8TeV_pf_noRegression_noMassCut_v*.root", "qcd_40_8TeV_pf", "(evweight_w_btagSF) * (evweight_w_btagSF < 50)", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
samples.append(["gjet_20_8TeV_pf", 20, afs_plottree, "", "gjet_20_8TeV_pf_noRegression_noMassCut_v*.root", "gjet_20_8TeV_pf", "evweight_w_btagSF", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
samples.append(["gjet_40_8TeV_pf", 20, afs_plottree, "", "gjet_40_8TeV_pf_noRegression_noMassCut_v*.root", "gjet_40_8TeV_pf", "evweight_w_btagSF", ROOT.kBlue-4, 1001, "QCD #gamma + jets" , 1., 14404429])
samples.append(["diphojet_sherpa_8TeV", 10, afs_plottree, "", "diphojet_sherpa_8TeV_noRegression_noMassCut_v*.root", "diphojet_sherpa_8TeV", "evweight_w_btagSF", ROOT.kGreen+2, 1001, "QCD #gamma#gamma + jets" , 1., 14404429])
samples.append(["DYJetsToLL", 5, afs_plottree, "", "DYJetsToLL_noRegression_noMassCut_v*.root", "DYJetsToLL", "evweight_w_btagSF", ROOT.kMagenta+2, 1001, "DY" , 1., 14404429])


#####plots.append([ name2, variable, plot_cut, norm, binning, title, additional_info, cutline, cutline2 ])
plots = []
#plots.append(["pho1_pt", "pho1_pt", "", intL, "(100, 0, 500)", "p_{T}^{#gamma1} (GeV)", "", 33.3, ""])
#plots.append(["pho2_pt", "pho2_pt", "", intL, "(100, 0, 500)", "p_{T}^{#gamma2} (GeV)", "", 25., ""])
#plots.append(["jet1_pt", "jet1_pt", "", intL, "(100, 0, 500)", "p_{T}^{jet1} (GeV)", "", 25., ""])
#plots.append(["jet2_pt", "jet2_pt", "", intL, "(100, 0, 500)", "p_{T}^{jet1} (GeV)", "", 25., ""])
#plots.append(["pho1_eta", "pho1_eta", "", intL, "(100, -5, 5)", "eta^{#gamma1}", "", 2.5, -2.5])
#plots.append(["pho2_eta", "pho2_eta", "", intL, "(100, -5, 5)", "eta^{#gamma2}", "", 2.5, -2.5])
#plots.append(["jet1_eta", "jet1_eta", "", intL, "(100, -5, 5)", "eta^{jet1}", "", "", ""])
#plots.append(["jet2_eta", "jet2_eta", "", intL, "(100, -5, 5)", "eta^{jet1}", "", "", ""])
#plots.append(["pho1_phi", "pho1_phi", "", intL, "(100, -5, 5)", "phi^{#gamma1}", "", "", ""])
#plots.append(["pho2_phi", "pho2_phi", "", intL, "(100, -5, 5)", "phi^{#gamma2}", "", "", ""])
#plots.append(["jet1_phi", "jet1_phi", "", intL, "(100, -5, 5)", "phi^{jet1}", "", "", ""])
#plots.append(["jet2_phi", "jet2_phi", "", intL, "(100, -5, 5)", "phi^{jet2}", "", "", ""])
#plots.append(["jj_phi", "jj_phi", "", intL, "(100, -5, 5)", "phi^{jj}", "", "", ""])
#plots.append(["jj_eta", "jj_eta", "", intL, "(100, -5, 5)", "eta^{jj}", "", "", ""])
#plots.append(["gg_phi", "gg_phi", "", intL, "(100, -5, 5)", "phi^{#gamma#gamma}", "", "", ""])
#plots.append(["gg_eta", "gg_eta", "", intL, "(100, -5, 5)", "eta^{#gamma#gamma}", "", "", ""])
#plots.append(["jj_pt", "jj_pt", "", intL, "(500, 0, 500)", "p_{T}^{jj} (GeV)", "", "", ""])
#plots.append(["gg_pt", "gg_pt", "", intL, "(500, 0, 500)", "p_{T}^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_mass", "pho1_mass", "", intL, "(100, -1, 1)", "mass^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_mass", "pho2_mass", "", intL, "(100, -1, 1)", "mass^{#gamma2} (GeV)", "", "", ""])
#plots.append(["jet1_mass", "jet1_mass", "", intL, "(100, 0, 100)", "mass^{jet1} (GeV)", "", "", ""])
#plots.append(["jet2_mass", "jet2_mass", "", intL, "(100, 0, 100)", "mass^{jet2} (GeV)", "", "", ""])
plots.append(["jj_mass_normData", "jj_mass", "", "data", "(50, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
plots.append(["gg_mass_normData", "gg_mass", "", "data", "(40, 100, 180)", "mass^{#gamma#gamma} (GeV)", "", "", ""])
plots.append(["jj_mass_norm1", "jj_mass", "", 1, "(50, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
plots.append(["gg_mass_norm1", "gg_mass", "", 1, "(40, 100, 180)", "mass^{#gamma#gamma} (GeV)", "", "", ""])
plots.append(["jj_mass", "jj_mass", "", intL, "(50, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
plots.append(["gg_mass", "gg_mass", "", intL, "(40, 100, 180)", "mass^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_e", "pho1_e", "", intL, "(500, 0, 500)", "e^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_e", "pho2_e", "", intL, "(500, 0, 500)", "e^{#gamma2} (GeV)", "", "", ""])
#plots.append(["jet1_e", "jet1_e", "", intL, "(500, 0, 500)", "e^{jet1} (GeV)", "", "", ""])
#plots.append(["jet2_e", "jet2_e", "", intL, "(500, 0, 500)", "e^{jet2} (GeV)", "", "", ""])
#plots.append(["jj_e", "jj_e", "", intL, "(500, 0, 500)", "e^{jj} (GeV)", "", "", ""])
#plots.append(["gg_e", "gg_e", "", intL, "(500, 0, 500)", "e^{#gamma#gamma} (GeV)", "", "", ""])
#plots.append(["pho1_r9", "pho1_r9", "", intL, "(100, 0, 5)", "r9^{#gamma1}", "", "", ""])
#plots.append(["pho2_r9", "pho2_r9", "", intL, "(100, 0, 5)", "r9^{#gamma2}", "", "", ""])
#plots.append(["ggjj_pt", "ggjj_pt", "", intL, "(500, 0, 500)", "p_{T}^{#gamma#gammajj} (GeV)", "", "", ""])
#plots.append(["ggjj_e", "ggjj_e", "", intL, "(1000, 0, 1000)", "e^{#gamma#gammajj} (GeV)", "", "", ""])
plots.append(["ggjj_mass_normData", "ggjj_mass / 1000.", "", "data", "(50, 0, 1.)", "mass^{#gamma#gammajj} (TeV)", "", "", ""])
plots.append(["ggjj_mass_norm1", "ggjj_mass / 1000.", "", 1, "(50, 0, 1.)", "mass^{#gamma#gammajj} (TeV)", "", "", ""])
plots.append(["ggjj_mass", "ggjj_mass / 1000.", "", intL, "(50, 0, 1.)", "mass^{#gamma#gammajj} (TeV)", "", "", ""])
#plots.append(["ggjj_phi", "ggjj_phi", "", intL, "(100, -5, 5)", "phi^{#gamma#gammajj}", "", "", ""])
#plots.append(["ggjj_eta", "ggjj_eta", "", intL, "(100, -10, 10)", "eta^{#gamma#gammajj}", "", "", ""])
#plots.append(["pho1_sieie", "pho1_sieie", "", intL, "(100, -0.5, 0.5)", "sieie^{#gamma1}", "", "", ""])
#plots.append(["pho2_sieie", "pho2_sieie", "", intL, "(100, -0.5, 0.5)", "sieie^{#gamma2}", "", "", ""])
#plots.append(["jet1_betaStarClassic", "jet1_betaStarClassic", "", intL, "(20, 0, 1)", "betaStarClassic^{jet1}", "", "", ""])
#plots.append(["jet2_betaStarClassic", "jet2_betaStarClassic", "", intL, "(20, 0, 1)", "betaStarClassic^{jet2}", "", "", ""])
#plots.append(["jet1_dR2Mean", "jet1_dR2Mean", "", intL, "(80,0, 0.2)", "dR2Mean^{jet1}", "", "", ""])
#plots.append(["jet2_dR2Mean", "jet2_dR2Mean", "", intL, "(80,0, 0.2)", "dR2Mean^{jet2}", "", "", ""])
#plots.append(["pho1_hoe", "pho1_hoe", "", intL, "(100,0, .05)", "hoe^{#gamma1}", "", "", ""])
#plots.append(["pho2_hoe", "pho2_hoe", "", intL, "(100,0, .05)", "hoe^{#gamma2}", "", "", ""])
#plots.append(["pho1_PFisoA", "pho1_PFisoA", "", intL, "(100,0, 10)", "PFisoA^{#gamma1}", "", "", ""])
#plots.append(["pho1_PFisoB", "pho1_PFisoB", "", intL, "(150,-5, 10)", "PFisoB^{#gamma1}", "", "", ""])
#plots.append(["pho2_PFisoA", "pho2_PFisoA", "", intL, "(100,0, 10)", "PFisoA^{#gamma2}", "", "", ""])
#plots.append(["pho2_PFisoB", "pho2_PFisoB", "", intL, "(150,-5, 10)", "PFisoB^{#gamma2}", "", "", ""])
#plots.append(["pho1_isconv", "pho1_isconv", "", intL, "(100,0, 1)", "isconv^{#gamma1}", "", "", ""])
#plots.append(["pho2_isconv", "pho2_isconv", "", intL, "(100,0, 1)", "isconv^{#gamma2}", "", "", ""])
#plots.append(["jet1_csvBtag", "jet1_csvBtag", "", intL, "(100, 0, 1)", "csvBtag^{jet1}", "", "", ""])
#plots.append(["jet2_csvBtag", "jet2_csvBtag", "", intL, "(100, 0, 1)", "csvBtag^{jet2}", "", "", ""])




for name2, variable, plot_cut, norm, binning, title, additional_info, cutline, cutline2 in plots:
    c1 = TCanvas()
    legend = TLegend(0.45, 0.82, 0.90, 0.93, "")
    legend.SetTextSize(0.025)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    xnbin, xlow, xhigh = map(float, binning.strip().strip("()").split(","))
    ymax = -1
    ymin = 10000000
    firsthistname = ""
    if plot_cut == "": plot_cut = "1"
    hist_signal = {}
    hist_data = {}
    hist_bkg = {}
    label_signal = {}
    label_data = {}
    label_bkg = {}

    for ifile, [ name, typ, dirpath, subdir, file, tree, sample_cut, color, style, label , sigma , N] in enumerate(samples):
#        print ifile, file, color, style, label
        chain = TChain(tree)
        chain.Add( path.join(dirpath, subdir, file) )
        total_cut = plot_cut
        if sample_cut == "": sample_cut = "1"
        if typ < 0:
            total_cut = "(" + plot_cut + ") * (" + str(sigma) + " * " + str(intL) + ")/" + str(N)
        elif typ == 0:
           total_cut = "(" + plot_cut + ")"
        elif typ > 0:
           total_cut = "(" + plot_cut + ") * (" + sample_cut + ")"
        option = ""
        if ifile != 0:
            option = "same"
        if typ == 0:
            option += "e1"
        chain.Draw(variable + ">>h_tmp" + binning, total_cut, option)
        # Clsosmetics
        h = ROOT.gDirectory.Get("h_tmp")
        h.SetName(name + "_" + name2 + "_" + str(ifile))
        if ifile == 0:
            firsthistname = name + "_" + name2 + "_" + str(ifile)
        h.SetLineWidth(2)
        h.SetLineColor(color)
        h.SetFillColor(color)
        h.SetFillStyle(style)
        h.GetXaxis().SetTitle( title )
        unit = ""
        if title.find("(") != -1:
            unit = title[title.find("(")+1:title.find(")")]
        if norm == 1. or norm == 1:
            h.GetYaxis().SetTitle( "Norm. to unity / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        elif norm == "data" or norm == "Data":
            h.GetYaxis().SetTitle( "Norm. to data / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        else:
            h.GetYaxis().SetTitle( "# events / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        # store histo for redraw in the correct order later
        if typ > 0:
            label_bkg[typ] = label
            if typ not in hist_bkg:
                hist_bkg[typ] = h
            else:
                hist_bkg[typ].Add(h)
        elif typ == 0:
            label_data[typ] = label
            hist_data[typ] = h
        elif typ < 0:
            label_signal[typ] = label
            hist_signal[typ] = h
        del chain, h

#        print hist_bkg
#        print hist_signal
#        print hist_data

    # Sum the backgrounds
    for key in collections.OrderedDict(sorted(hist_bkg.items())):
        for jkey in collections.OrderedDict(sorted(hist_bkg.items())):
            if jkey <= key: continue
            hist_bkg[key].Add(hist_bkg[jkey])
    # Adjust norm if case happens
    if norm == 1. or norm == 1:
        bkg_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_bkg.items()))):
            if ikey == 0:
                bkg_integral = hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
            hist_bkg[key].Scale( 1. / bkg_integral )
        data_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_data.items()))):
            if ikey == 0:
                data_integral = hist_data[key].Integral(0, hist_data[key].GetNbinsX() +1)
            hist_data[key].Scale( 1. / data_integral )
        signal_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_signal.items()))):
            if ikey == 0:
                signal_integral = hist_signal[key].Integral(0, hist_signal[key].GetNbinsX() +1)
            hist_signal[key].Scale( 1. / signal_integral )
    elif norm == "data" or norm == "Data":
        data_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_data.items()))):
            if ikey == 0:
                data_integral = hist_data[key].Integral(0, hist_data[key].GetNbinsX() +1)
            else:
                continue
        bkg_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_bkg.items()))):
            if ikey == 0:
                bkg_integral = hist_bkg[key].Integral(0, hist_bkg[key].GetNbinsX() + 1)
            hist_bkg[key].Scale( data_integral / bkg_integral )
        signal_integral = 1.
        for ikey, key in enumerate(collections.OrderedDict(sorted(hist_signal.items()))):
            if ikey == 0:
                signal_integral = hist_signal[key].Integral(0, hist_signal[key].GetNbinsX() +1)
            hist_signal[key].Scale( data_integral / signal_integral )
    # redraw in order : background, data, signal, axis
    if len(hist_bkg) + len(hist_data) + len(hist_signal) > 1:
        legend.SetNColumns(2)
    if len(hist_bkg) + len(hist_data) + len(hist_signal) > 6:
        legend.SetNColumns(3)
    for ikey, key in enumerate(collections.OrderedDict(sorted(hist_bkg.items()))):
        if ikey == 0:
            hist_bkg[key].Draw("")
            firsthistname = hist_bkg[key].GetName()
        else:
            hist_bkg[key].Draw("same")
        legend.AddEntry(hist_bkg[key].GetName(), label_bkg[key], "lf")
        ymax = max(ymax, hist_bkg[key].GetMaximum())
        ymin = min(ymin, hist_bkg[key].GetMinimum(0.0))
    for key in hist_data:
        hist_data[key].Draw("e1same")
        legend.AddEntry(hist_data[key].GetName(), label_data[key], "lpe")
        ymax = max(ymax, hist_data[key].GetMaximum())
        ymin = min(ymin, hist_data[key].GetMinimum(0.0))
    for key in collections.OrderedDict(sorted(hist_signal.items())):
        hist_signal[key].Draw("same")
        legend.AddEntry(hist_signal[key].GetName(), label_signal[key], "lf")
        ymax = max(ymax, hist_signal[key].GetMaximum())
        ymin = min(ymin, hist_signal[key].GetMinimum(0.0))
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
    c1.Update()

    line = TLine()
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line2 = TLine()
    line2.SetLineStyle(2)
    line2.SetLineWidth(2)

    h = ROOT.gDirectory.Get(firsthistname)
    h.SetMaximum(ymax_lin)
    h.SetMinimum(ymin_lin)
    if cutline != "":
        line.SetX1(cutline); line.SetY1(ymin_lin); line.SetX2(cutline); line.SetY2(ymax)
        line.Draw("same")
    if cutline2 != "":
        line2.SetX1(cutline2); line2.SetY1(ymin_lin); line2.SetX2(cutline2); line2.SetY2(ymax)
        line2.Draw("same")
    c1.Update()
    c1.Print("pdf/" + name2 + ".pdf")
    c1.Print("gif/" + name2 + ".gif")
    c1.Print("root/" + name2 + ".root")


    c1.SetLogy(1)
    h.SetMaximum(ymax_log)
    h.SetMinimum(ymin_log)
    h.GetYaxis().SetRangeUser(ymin_log, ymax_log)
    if cutline != "":
        line.SetX1(cutline); line.SetY1(ymin_log); line.SetX2(cutline); line.SetY2(ymax)
        line.Draw("same")
    if cutline2 != "":
        line2.SetX1(cutline2); line2.SetY1(ymin_log); line2.SetX2(cutline2); line2.SetY2(ymax)
        line2.Draw("same")
    c1.Update()
    c1.Print("pdf/" + name2 + "_log.pdf")
    c1.Print("gif/" + name2 + "_log.gif")
    c1.Print("root/" + name2 + "_log.root")
    c1.SetLogy(0)

    del c1

