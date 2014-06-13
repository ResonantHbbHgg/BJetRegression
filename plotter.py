#!/usr/bin/env python
#// Small dumb code to play with trees
#// O.Bondu, F. Bojarski (May 2014)
# Various python imports
from os import path
from math import log10, pow
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
# samples.append([ name, typ, dirpath, subdir, file, tree, color, style, label , sigma , N])
samples.append(["Radion_m300", -300, afs_plottree, "", "Radion_m300_8TeV_noRegression_noMassCut_v*.root", "Radion_m300_8TeV", ROOT.kBlue, 0, "m_{X} = 300 GeV" , 13.55e-3, 19999])
samples.append(["Data", 0, afs_plottree, "", "Data_noRegression_noMassCut_v*.root", "Data", ROOT.kBlack, 0, "Data" , 13.55e-3, 19999])

#####plots.append([ name2, variable, cut, norm, binning, title, additional_info, cutline, cutline2 ])
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
#plots.append(["jj_mass", "jj_mass", "", intL, "(500, 0, 500)", "mass^{jj} (GeV)", "", "", ""])
#plots.append(["gg_mass", "gg_mass", "", intL, "(100, 0, 200)", "mass^{#gamma#gamma} (GeV)", "", 100., 180.])
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
plots.append(["ggjj_mass", "ggjj_mass", "", intL, "(50, 0, 1000)", "mass^{#gamma#gammajj} (GeV)", "", "", ""])
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




for name2, variable, cut, norm, binning, title, additional_info, cutline, cutline2 in plots:
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


    for ifile, [ name, typ, dirpath, subdir, file, tree, color, style, label , sigma , N] in enumerate(samples):
#        print ifile, file, color, style, label
        chain = TChain(tree)
        chain.Add( path.join(dirpath, subdir, file) )
        sample_cut = cut
        if norm == 1.:
            sample_cut = "(" + sample_cut + ")/" + str( chain.GetEntries() )
        else:
            sample_cut = "(" + sample_cut + ") * (" + str(sigma) + " * " + str(intL) + ")/" + str(N)
        option = ""
        if ifile != 0:
            option = "same"
        chain.Draw(variable + ">>h_tmp" + binning, sample_cut, option)
        # Clsosmetics
        h = ROOT.gDirectory.Get("h_tmp")
        h.SetName(name + "_" + name2 + "_" + str(ifile))
        if ifile == 0:
            firsthistname = name + "_" + name2 + "_" + str(ifile)
        h.SetLineWidth(3)
        h.SetLineColor(color)
        h.SetFillColor(color)
        h.SetFillStyle(style)
        h.GetXaxis().SetTitle( title )
        unit = ""
        if title.find("(") != -1:
            unit = title[title.find("(")+1:title.find(")")]
        if norm == 1.:
            h.GetYaxis().SetTitle( "Norm. to unity / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        else:
            h.GetYaxis().SetTitle( "# events / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
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

