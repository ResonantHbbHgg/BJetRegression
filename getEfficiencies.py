import numpy
from ROOT import gROOT
from ROOT import TCanvas
from ROOT import TChain
from ROOT import TGraphErrors
from ROOT import TLatex

gROOT.SetBatch()

signals = []
# For the efficiency plot
signals.append(["v29_fitToMgg_noKinFit/MSSM_m260_8TeV_m260.root", 300000])
signals.append(["v29_fitToMgg_noKinFit/Radion_m270_8TeV_m270.root", 19996])
signals.append(["v29_fitToMgg_noKinFit/Radion_m300_8TeV_m300.root", 19972])
signals.append(["v29_fitToMgg_noKinFit/Radion_m350_8TeV_m350.root", 18498])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m400_8TeV_m400.root", 19697])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m450_8TeV_m450.root", 19999])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m500_8TeV_m500.root", 19970])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m550_8TeV_m550.root", 19995])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m600_8TeV_m600.root", 18197])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m650_8TeV_m650.root", 20000])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m700_8TeV_m700.root", 19969])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m800_8TeV_m800.root", 19999])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m900_8TeV_m900.root", 19996])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m1000_8TeV_m1000.root", 19951])
signals.append(["v29_fitToMggjj_withKinFit/Radion_m1100_8TeV_m1100.root", 19400])

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
n = 0
x = []
y = []
ex = []
ey = []

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
    n = n+1
    x.append(signal.split("_")[3].split("m")[1])
    y.append(round(eff * 100, 2))
    ex.append(0.)
    ey.append(round(deff * 100,2))
    del chain

x = numpy.asarray(x, dtype='f')
y = numpy.asarray(y, dtype='f')
ex = numpy.asarray(ex, dtype='f')
ey = numpy.asarray(ey, dtype='f')
# Removing y errors
ey = ex 


#gROOT.ProcessLine(".x CMSStyle.C")
gROOT.ProcessLine(".x setTDRStyle.C")
c1 = TCanvas()
gr = TGraphErrors(n,x,y,ex,ey);
gr.SetMarkerStyle(20)
gr.SetTitle("")
gr.GetXaxis().SetTitle("M_{X} (GeV)")
gr.GetYaxis().SetTitle("Signal selection efficiency (%)")
gr.Draw("ALP");
latexLabel = TLatex()
latexLabel.SetTextSize(0.03)
latexLabel.SetNDC()
latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary   #sqrt{s} = 8 TeV   L = 19.7 fb^{-1}")
c1.Print("eff.pdf")
c1.Print("eff.root")
