import os, re
import ROOT
ROOT.gROOT.SetBatch()

class result(object):
	def __init__(self, method, mass, var, cat, fwhm, mean):
		self.method = method
		self.mass = mass
		self.var = var
		self.cat = cat
		self.fwhm = fwhm
		self.mean = mean

	def __str__(self):
		return "method= " + str(self.method) + "\tmass= " + str(self.mass) + "\tvar= " + str(self.var) + "\tcat= " + str(self.cat) + "\twidth= " + str(self.fwhm) + "\tmean= " + str(self.mean)

def getFWHM(workspace, var_name, pdf_name):
	var = workspace.var(var_name)
	pdf = workspace.pdf(pdf_name)
	f1 = pdf.asTF(ROOT.RooArgList(var))
	max = f1.GetMaximum()
	xmax = f1.GetMaximumX()
	x_hig = f1.GetX(max / 2., xmax, var.getMax())
	x_low = f1.GetX(max / 2., var.getMin(), xmax)
	return x_hig-x_low


mypath = "/afs/cern.ch/work/o/obondu/Higgs/HIG-13-032/low_mass/"
lowmass_dir = [ d for d in os.listdir(mypath) if os.path.isdir(os.path.join(mypath,d)) and "radlim" in d]
lowmass_fil = [ os.path.join(mypath, d, "hgg.mH125.6_8TeV.inputsig.root") for d in lowmass_dir]

mypath = "/afs/cern.ch/work/o/obondu/Higgs/HIG-13-032/high_mass/"
higmass_fil = [os.path.join(mypath, f) for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f)) and "inputsig" in f]

allresults = []

for filename in lowmass_fil:
	file = ROOT.TFile(filename)
	ws = file.Get("w_all")
	dir = filename.split('/')[-2].split('_')[0]
	mass = re.findall('\d+', dir)[0]
	mass = int(mass)
	res = result("fwhm", mass, "mgg", "cat0", getFWHM(ws, "mgg", "CMS_hgg_sig_cat0"), 125.)
	allresults.append(res)
	res = result("fwhm", mass, "mgg", "cat1", getFWHM(ws, "mgg", "CMS_hgg_sig_cat1"), 125.)
	allresults.append(res)
	

for filename in higmass_fil:
	file = ROOT.TFile(filename)
	ws = file.Get("w_all")
	mass = os.path.basename(filename).split('.')[1]
	mass = re.findall('\d+', mass)[0]
	mass = int(mass)
	res = result("fwhm", mass, "mtot", "cat0", getFWHM(ws, "mtot", "CMS_hgg_sig_cat0"), mass)
	allresults.append(res)
	res = result("fwhm", mass, "mtot", "cat1", getFWHM(ws, "mtot", "CMS_hgg_sig_cat1"), mass)
	allresults.append(res)

for res in allresults:
	print res
