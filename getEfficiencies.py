from ROOT import gROOT
from ROOT import TChain
gROOT.SetBatch()

signals = []
# radion samples
signals.append(["v26_fitToMgg_noKinFit/Radion_m270_8TeV_m270.root", 19996])
signals.append(["v26_fitToMgg_noKinFit/Radion_m300_8TeV_m300.root", 19972])
signals.append(["v26_fitToMgg_noKinFit/Radion_m350_8TeV_m350.root", 18498])
signals.append(["v26_fitToMgg_noKinFit/Radion_m400_8TeV_m400.root", 19697])
signals.append(["v26_fitToMgg_noKinFit/Radion_m450_8TeV_m450.root", 19999])
signals.append(["v26_fitToMgg_noKinFit/Radion_m500_8TeV_m500.root", 19970])
# graviton samples
signals.append(["v26_fitToMgg_noKinFit/MSSM_m260_8TeV_m260.root", 300000])
signals.append(["v26_fitToMgg_noKinFit/MSSM_m300_8TeV_m300.root", 299142])
signals.append(["v26_fitToMgg_noKinFit/MSSM_m350_8TeV_m350.root", 299571])
# mssm samples
signals.append(["v26_fitToMgg_noKinFit/Graviton_m300_8TeV_m300.root", 49941])
signals.append(["v26_fitToMgg_noKinFit/Graviton_m500_8TeV_m500.root", 49905])
#signals.append(["", ])


#print signals
lumi = 19706.

for signal, nprocessed in signals:
	chain = TChain("TCVARS")
	chain.Add(signal)
	totw = 0
	for ievt in xrange(chain.GetEntries()):
		chain.GetEntry(ievt)
		totw += chain.evWeight
	totw *= 1000. # due to hand-put 1000 factor to weight for limit settings numerical precision
	eff = totw / lumi
	deff = (eff * (1 - eff) / nprocessed)**0.5
#	print signal, round(eff * 100, 2), round(deff * 100,2)
	print signal.split("_")[3].split("m")[1], "&", "$", round(eff * 100, 2), "$", "&", "$\pm", round(deff * 100,2), "$"
	del chain
