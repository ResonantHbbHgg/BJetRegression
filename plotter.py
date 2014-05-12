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
afs_plottree = "/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v10/"
eos_tree = "root://eoscms//eos/cms/store/cmst3/user/obondu/H2GGLOBE/Radion/trees/radion_redu_11_tree_08"

intL = 19706.
samples = []
# samples.append([ name, dirpath, subdir, file, tree, color, style, label , sigma , N])
samples.append(["Radion_m300", afs_plottree, "2014-02-17_selection_noRegression_noMassCut_v10", "Radion_m300_8TeV_noRegression_noMassCut_v10.root", "Radion_m300_8TeV", ROOT.kBlue, 0, "m_{X} = 300 GeV" , 13.55e-3, 19999])
samples.append(["Radion_m300", "/afs/cern.ch/work/f/fbojarsk/CMSSW_6_1_1/src/BJetRegression/", "2014-05-12_selection_noRegression_noMassCut_v10bis_effStudies", "Radion_m300_8TeV_noRegression_noMassCut_v10bis_effStudies.root", "Radion_m300_8TeV", ROOT.kRed, 0, "m_{X} = 300 GeV sans cut", 13.55e-3, 19999])

plots = []
plots.append(["pho1_pt", "pho1_pt", "", 2., "(100, 0, 500)", "p_{T}^{#gamma1} (GeV)", ""])
plots.append(["pho2_pt", "pho2_pt", "", 2., "(100, 0, 500)", "p_{T}^{#gamma2} (GeV)", ""])
plots.append(["jet1_pt", "jet1_pt", "", 2., "(100, 0, 500)", "p_{T}^{jet1} (GeV)", ""])
plots.append(["jet2_pt", "jet2_pt", "", 2., "(100, 0, 500)", "p_{T}^{jet1} (GeV)", ""])
plots.append(["pho1_eta", "pho1_eta", "", 2., "(100, -5, 5)", "eta^{#gamma1}", ""])
plots.append(["pho2_eta", "pho2_eta", "", 2., "(100, -5, 5)", "eta^{#gamma2}", ""])
plots.append(["jet1_eta", "jet1_eta", "", 2., "(100, -5, 5)", "eta^{jet1}", ""])
plots.append(["jet2_eta", "jet2_eta", "", 2., "(100, -5, 5)", "eta^{jet1}", ""])
plots.append(["pho1_phi", "pho1_phi", "", 2., "(100, -5, 5)", "phi^{#gamma1}", ""])
plots.append(["pho2_phi", "pho2_phi", "", 2., "(100, -5, 5)", "phi^{#gamma2}", ""])
plots.append(["jet1_phi", "jet1_phi", "", 2., "(100, -5, 5)", "phi^{jet1}", ""])
plots.append(["jet2_phi", "jet2_phi", "", 2., "(100, -5, 5)", "phi^{jet2}", ""])
plots.append(["jj_phi", "jj_phi", "", 2., "(100, -5, 5)", "phi^{jj}", ""])
plots.append(["jj_eta", "jj_eta", "", 2., "(100, -5, 5)", "eta^{jj}", ""])
plots.append(["gg_phi", "gg_phi", "", 2., "(100, -5, 5)", "phi^{#gamma#gamma}", ""])
plots.append(["gg_eta", "gg_eta", "", 2., "(100, -5, 5)", "eta^{#gamma#gamma}", ""])
plots.append(["jj_pt", "jj_pt", "", 2., "(500, 0, 500)", "p_{T}^{jj} (GeV)", ""])
plots.append(["gg_pt", "gg_pt", "", 2., "(500, 0, 500)", "p_{T}^{#gamma#gamma} (GeV)", ""])
plots.append(["pho1_mass", "pho1_mass", "", 2., "(100, -1, 1)", "mass^{#gamma1} (GeV)", ""])
plots.append(["pho2_mass", "pho2_mass", "", 2., "(100, -1, 1)", "mass^{#gamma2} (GeV)", ""])
plots.append(["jet1_mass", "jet1_mass", "", 2., "(100, 0, 100)", "mass^{jet1} (GeV)", ""])
plots.append(["jet2_mass", "jet2_mass", "", 2., "(100, 0, 100)", "mass^{jet2} (GeV)", ""])
plots.append(["jj_mass", "jj_mass", "", 2., "(500, 0, 500)", "mass^{jj} (GeV)", ""])
plots.append(["gg_mass", "gg_mass", "", 2., "(100, 0, 200)", "mass^{#gamma#gamma} (GeV)", ""])
plots.append(["pho1_e", "pho1_e", "", 2., "(500, 0, 500)", "e^{#gamma1} (GeV)", ""])
plots.append(["pho2_e", "pho2_e", "", 2., "(500, 0, 500)", "e^{#gamma2} (GeV)", ""])
plots.append(["jet1_e", "jet1_e", "", 2., "(500, 0, 500)", "e^{jet1} (GeV)", ""])
plots.append(["jet2_e", "jet2_e", "", 2., "(500, 0, 500)", "e^{jet2} (GeV)", ""])
plots.append(["jj_e", "jj_e", "", 2., "(500, 0, 500)", "e^{jj} (GeV)", ""])
plots.append(["gg_e", "gg_e", "", 2., "(500, 0, 500)", "e^{#gamma#gamma} (GeV)", ""])
plots.append(["pho1_r9", "pho1_r9", "", 2., "(100, 0, 5)", "r9^{#gamma1}", ""])
plots.append(["pho2_r9", "pho2_r9", "", 2., "(100, 0, 5)", "r9^{#gamma2}", ""])
plots.append(["ggjj_pt", "ggjj_pt", "", 2., "(500, 0, 500)", "p_{T}^{#gamma#gammajj} (GeV)", ""])
plots.append(["ggjj_e", "ggjj_e", "", 2., "(1000, 0, 1000)", "e^{#gamma#gammajj} (GeV)", ""])
plots.append(["ggjj_mass", "ggjj_mass", "", 2., "(200, 0, 00)", "mass^{#gamma#gammajj} (GeV)", ""])
plots.append(["ggjj_phi", "ggjj_phi", "", 2., "(100, -5, 5)", "phi^{#gamma#gammajj}", ""])
plots.append(["ggjj_eta", "ggjj_eta", "", 2., "(100, -10, 10)", "eta^{#gamma#gammajj}", ""])
plots.append(["pho1_sieie", "pho1_sieie", "", 2., "(100, -0.5, 0.5)", "sieie^{#gamma1}", ""])
plots.append(["pho2_sieie", "pho2_sieie", "", 2., "(100, -0.5, 0.5)", "sieie^{#gamma2}", ""])
plots.append(["jet1_betaStarClassic", "jet1_betaStarClassic", "", 2., "(20, 0, 1)", "betaStarClassic^{jet1}", ""])
plots.append(["jet2_betaStarClassic", "jet2_betaStarClassic", "", 2., "(20, 0, 1)", "betaStarClassic^{jet2}", ""])
plots.append(["jet1_dR2Mean", "jet1_dR2Mean", "", 2., "(25,0, 0.1)", "dR2Mean^{jet1}", ""])
plots.append(["jet2_dR2Mean", "jet2_dR2Mean", "", 2., "(25,0, 0.1)", "dR2Mean^{jet2}", ""])
plots.append(["pho1_hoe", "pho1_hoe", "", 2., "(100,0, .05)", "hoe^{#gamma1}", ""])
plots.append(["pho2_hoe", "pho2_hoe", "", 2., "(100,0, .05)", "hoe^{#gamma2}", ""])
plots.append(["pho1_PFisoA", "pho1_PFisoA", "", 2., "(100,0, 10)", "PFisoA^{#gamma1}", ""])
plots.append(["pho1_PFisoB", "pho1_PFisoB", "", 2., "(150,-5, 10)", "PFisoB^{#gamma1}", ""])
plots.append(["pho2_PFisoA", "pho2_PFisoA", "", 2., "(100,0, 10)", "PFisoA^{#gamma2}", ""])
plots.append(["pho2_PFisoB", "pho2_PFisoB", "", 2., "(150,-5, 10)", "PFisoB^{#gamma2}", ""])
plots.append(["pho1_isconv", "pho1_isconv", "", 2., "(100,0, 1)", "isconv^{#gamma1}", ""])
plots.append(["pho2_isconv", "pho2_isconv", "", 2., "(100,0, 1)", "isconv^{#gamma2}", ""])
plots.append(["jet1_csvBtag", "jet1_csvBtag", "", 2., "(100, 0, 1)", "csvBtag^{jet1}", ""])
plots.append(["jet2_csvBtag", "jet2_csvBtag", "", 2., "(100, 0, 1)", "csvBtag^{jet2}", ""])




for name2, variable, cut, norm, binning, title, additional_info in plots:
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


    for ifile, [ name, dirpath, subdir, file, tree, color, style, label , sigma , N] in enumerate(samples):
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
        h.SetName(name + str(ifile))
        if ifile == 0:
            firsthistname = name + str(ifile)
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

    if name2 == "jet1_pt":
        line = TLine(25, ymin, 25, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")

    if name2 == "jet2_pt":
        line = TLine(25, ymin, 25, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME") 

    if name2 == "pho1_pt":
        line = TLine(33.3, ymin, 33.3, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")

    if name2 == "pho2_pt":
        line = TLine(25, ymin, 25, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")

    if name2 == "jet1_eta":
        line = TLine(2.5, ymin, 2.5, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")
        line2 = TLine(-2.5, ymin, -2.5, ymax)
        line2.SetLineStyle(2)
        line2.Draw("same")


    if name2 == "jet2_eta":
        line = TLine(2.5, ymin, 2.5, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")
        line2 = TLine(-2.5, ymin, -2.5, ymax)
        line2.SetLineStyle(2)
        line2.Draw("same")

    if name2 == "pho1_eta":
        line = TLine(2.5, ymin, 2.5, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")
        line2 = TLine(-2.5, ymin, -2.5, ymax)
        line2.SetLineStyle(2)
        line2.Draw("same")
 
    if name2 == "pho2_eta":
        line = TLine(2.5, ymin, 2.5, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")
        line2 = TLine(-2.5, ymin, -2.5, ymax)
        line2.SetLineStyle(2)
        line2.Draw("same")

    if name2 == "gg_mass":
        line = TLine(100, ymin,100, ymax)
        line.SetLineStyle(2)
        line.Draw("SAME")
        line2 = TLine(180, ymin, 180, ymax)
        line2.SetLineStyle(2)
        line2.Draw("same")

    h = ROOT.gDirectory.Get(firsthistname)
    h.SetMaximum(ymax_lin)
    h.SetMinimum(ymin_lin)
    c1.Update()
    c1.Print("pdf/" + name2 + ".pdf")
    c1.Print("gif/" + name2 + ".gif")
    c1.Print("root/" + name2 + ".root")


    c1.SetLogy(1)
    h.SetMaximum(ymax_log)
    h.SetMinimum(ymin_log)
    h.GetYaxis().SetRangeUser(ymin_log, ymax_log)
    c1.Update()
    c1.Print("pdf/" + name2 + "_log.pdf")
    c1.Print("gif/" + name2 + "_log.gif")
    c1.Print("root/" + name2 + "_log.root")
    c1.SetLogy(0)

    del c1

