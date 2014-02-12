from ROOT import gROOT
from ROOT import TChain
gROOT.SetBatch()

samples = []
#samples.append(["2014-01-31_selection_noRegression_noMassCut_v06/DYJetsToLL_noRegression_noMassCut_v06.root", "DYJetsToLL"])
#samples.append(["2014-01-31_selection_noRegression_noMassCut_v06/LNuGG_FSR_8TeV_noRegression_noMassCut_v06.root", "LNuGG_FSR_8TeV"])
#samples.append(["2014-01-31_selection_noRegression_noMassCut_v06/LNuGG_ISR_8TeV_noRegression_noMassCut_v06.root", "LNuGG_ISR_8TeV"])
#samples.append(["2014-01-31_selection_noRegression_noMassCut_v06/tGG_8TeV_noRegression_noMassCut_v06.root", "tGG_8TeV"])
#samples.append(["2014-01-31_selection_noRegression_noMassCut_v06/ttGG_8TeV_noRegression_noMassCut_v06.root", "ttGG_8TeV"])
#samples.append(["2014-01-31_selection_noRegression_noMassCut_v06/TTGJets_8TeV_noRegression_noMassCut_v06.root", "TTGJets_8TeV"])
samples.append(["2014-02-12_selection_noRegression_noMassCut_v09/Radion_m300_8TeV_noRegression_noMassCut_v09.root", "Radion_m300_8TeV"])
#samples.append(["v28_fitToMgg_noKinFit/Radion_m300_8TeV_m300.root", "TCVARS"])

for sample, tree in samples:
	chain = TChain(tree)
	chain.Add(sample)
	totyield = 0
	for ievt in xrange(chain.GetEntries()):
		chain.GetEntry(ievt)
#		totyield += chain.evweight
#		if chain.cut_based_ct == 0:
#			totyield += chain.evWeight * 1000.
		if chain.evweight_w_btagSF < 10. and 85. < chain.jj_mass < 155. and 255. < chain.ggjj_mass < 330. and chain.category == 2:
			totyield += chain.evweight_w_btagSF
	del chain
	print sample.split("/")[1].split("_noRegression")[0], round(totyield,2)

