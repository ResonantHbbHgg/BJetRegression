#!/usr/bin/env python
from math import log10, pow
import sys
from os import path
import ROOT
from ROOT import gROOT, TGaxis, TLegend, TLatex
from ROOT import TChain, TCanvas, TH1D
gROOT.SetBatch()
here = sys.modules[__name__]


eos_prod = "/store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/"
eos_redu = "/store/cmst3/user/obondu/H2GGLOBE/Radion/reduced/radion_reduction_v11/mc"
samples = []
#samples.append(["Radion_m300", "prod", eos_prod, "Summer12_DR53X-PU_RD1_START53_V7N/RadionToHHTo2G2B_M-300_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kBlue, 3002, "m_{X} = 300 GeV"])
#samples.append(["ggH_m125", "prod", eos_prod, "Summer12_RD1/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kRed, 3345, "ggH"])
#samples.append(["VBF_m125", "prod", eos_prod, "Summer12_RD1/VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kGreen, 3315, "VBF"])
#samples.append(["VH_m125", "prod", eos_prod, "Summer12_RD1/WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2", "*root", "event", ROOT.kMagenta, 3325, "VH"])
#samples.append(["ttH_m125", "prod", eos_prod, "Summer12_RD1/TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kCyan, 3354, "ttH"])

samples.append(["Radion_m300", "redu", eos_redu, "RadionToHHTo2G2B_M-300_TuneZ2star_8TeV-nm-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kBlue, 3002, "m_{X} = 300 GeV"])
samples.append(["ggH_m125", "redu", eos_redu, "Summer12_RD1/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kRed, 3345, "ggH"])
#samples.append(["VBF_m125", "redu", eos_redu, "Summer12_RD1/VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kGreen, 3315, "VBF"])
#samples.append(["VH_m125", "redu", eos_redu, "Summer12_RD1/WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2", "*root", "event", ROOT.kMagenta, 3325, "VH"])
samples.append(["ttH_m125", "redu", eos_redu, "Summer12_RD1/TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1", "*root", "event", ROOT.kCyan, 3354, "ttH"])

#gROOT.ProcessLine(".x CMSStyle.C")
gROOT.ProcessLine(".x setTDRStyle.C")
TGaxis.SetMaxDigits(3)
nprocessed = {}
nreduced = {}

plots = []
#plots.append(["nphotons", "pho_n", "", 1., "(10, 0, 10)", "# photons", ""])
#plots.append(["njets", "jet_algoPF1_n", "", 1., "(30, 0, 150)", "# jets", ""])
#plots.append(["njets_pt10", "jet_algoPF1_n", "jet_algoPF1_p4.Pt() > 10", 1., "(50, 0, 200)", "# jets", "p_{T}^{jet} > 10 GeV"])
#plots.append(["njets_pt25", "jet_algoPF1_n", "jet_algoPF1_p4.Pt() > 25", 1., "(50, 0, 200)", "# jets", "p_{T}^{jet} > 25 GeV"])
#plots.append(["njets_pt25_1csvm", "jet_algoPF1_n", "jet_algoPF1_p4.Pt() > 25 && jet_algoPF1_csvBtag > 0", 1., "(50, 0, 200)", "# jets", "p_{T}^{jet} > 25 GeV, #geq 1 btag"])
#plots.append(["pho_ecalsumetconedr03", "pho_ecalsumetconedr03[0]", "", 1., "(25, 0, 10)", "ECAL iso (GeV)", ""])

plots.append(["pho_eta", "pho_p4[0].Eta()", "", 1., "(60, -3., 3.)", "#eta^{#gamma}", ""])
plots.append(["pho_pt", "pho_p4[0].Pt()", "", 1., "(50, 0, 200)", "p_{T}^{#gamma} (GeV)", ""])
plots.append(["pho_sieie", "pho_sieie[0]", "", 1., "(25, 0, .05)", "#sigma_{i#etai#eta}", ""])
plots.append(["pho_hoe", "pho_hoe[0]", "", 1., "(25, 0, 1)", "photon H/E", ""])
plots.append(["pho_ecaliso", "pho_ecalsumetconedr03[0] - 0.012*pho_p4[0].Et()", "", 1., "(25, 0, 10)", "I_{ecal} (GeV)", ""])
plots.append(["pho_hcaliso", "pho_hcalsumetconedr03[0] - 0.005*pho_p4[0].Et()", "", 1., "(25, 0, 10)", "I_{hcal} (GeV)", ""])
plots.append(["pho_trkiso", "pho_trksumpthollowconedr03[0] - 0.002*pho_p4[0].Et()", "", 1., "(25, 0, 10)", "I_{trk} (GeV)", ""])

plots.append(["phos_eta", "pho_p4.Eta()", "", 1., "(60, -3., 3.)", "#eta^{#gamma}", ""])
plots.append(["phos_pt", "pho_p4.Pt()", "", 1., "(50, 0, 200)", "p_{T}^{#gamma} (GeV)", ""])
plots.append(["phos_sieie", "pho_sieie", "", 1., "(25, 0, .05)", "#sigma_{i#etai#eta}", ""])
plots.append(["phos_hoe", "pho_hoe", "", 1., "(25, 0, 1)", "photon H/E", ""])
plots.append(["phos_ecaliso", "pho_ecalsumetconedr03 - 0.012*pho_p4.Et()", "", 1., "(25, 0, 10)", "I_{ecal} (GeV)", ""])
plots.append(["phos_hcaliso", "pho_hcalsumetconedr03 - 0.005*pho_p4.Et()", "", 1., "(25, 0, 10)", "I_{hcal} (GeV)", ""])
plots.append(["phos_trkiso", "pho_trksumpthollowconedr03 - 0.002*pho_p4.Et()", "", 1., "(25, 0, 10)", "I_{trk} (GeV)", ""])


for ivariable, [name, variable, cut, norm, binning, title, additional_info] in enumerate(plots):
    c1 = TCanvas()
    legend = TLegend(0.45, 0.82, 0.90, 0.93, "")
    legend.SetTextSize(0.025)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    if len(samples) > 1:
        legend.SetNColumns(2)
    xnbins, xlow, xhigh = map(float, binning.strip().strip("()").split(","))
    ymax = -1
    ymin = 10000000
    firsthistname = ""
    if cut == "": cut = "1."

    for ifile, [displayName, level, eos_base, sample, file, tree, color, style, label] in enumerate(samples):
        chain = TChain(tree)
        files = "root://eoscms//eos/cms" + path.join(eos_base, sample, file)
        chain.Add(files)
        print "ifile= ", ifile, "displayName= ", displayName, "level= ", level, "chain.GetEntries()= ", chain.GetEntries()
        if ivariable == 0:
            print "\tfiles= ", files
        sample_cut = cut
        if norm == 1.:
            sample_cut = "(" + sample_cut + ")/" + str( float(chain.GetEntries()) )
        option = ""
        if ifile != 0:
            option = "same"
        # Cosmetics
        unit = ""
        histname = "h_" + variable + str(ifile)
        if level == "prod":
            nprocessed[displayName] = chain.GetEntries()
        if level == "redu":
            nreduced[displayName] = chain.GetEntries()
    #    if level == "tree":
    #        nentries_w = 0
    #        for ievt in xrange(chain.GetEntries()):
    #            chain.GetEntry(ievt)
    #            nentries_w += chain.evweight
    #        print "displayName= ", displayName, "level= ", "trew", "chain.GetEntries()= ", nentries_w * nprocessed[displayName] / 19706.
        chain.Draw(variable + ">>h_tmp" + binning, sample_cut, option)
        h = ROOT.gDirectory.Get("h_tmp")
        h.SetName(histname)
        if ifile == 0:
            firsthistname = histname
        h.SetLineWidth(3)
        h.SetLineColor(color)
        h.SetFillColor(color)
        h.SetFillStyle(style)
        h.GetXaxis().SetTitle( title )
        if norm == 1.:
            if title.find("(") != -1:
                unit = title[title.find("(")+1:title.find(")")]
            h.GetYaxis().SetTitle( "Norm. to unity / ( " + str(((xhigh - xlow) / xnbins)) + " " + unit + " )")
        legend.AddEntry(h.GetName(), label, "lf")
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
    latexLabel.DrawLatex(.25, .96, "CMS Private   #sqrt{s} = 8TeV")
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
    del c1, legend, latexLabel


