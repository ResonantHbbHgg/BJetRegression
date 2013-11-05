#!/usr/bin/env python

from math import floor, ceil

def myround(x, base=5):
    return int(base * round(float(x)/base))
def upround(x, base=5):
    return int(base * ceil(float(x)/base))
def downround(x, base=5):
    return int(base * floor(float(x)/base))



file1='Radion_m300_performanceSummary_mggjj_nokin.txt'
file2='Radion_m500_performanceSummary_mggjj_nokin.txt'
file3='Radion_m300_performanceSummary_mggjj_kin.txt'
file4='Radion_m500_performanceSummary_mggjj_kin.txt'

# [i][0] is 1btag, [i][1] is 2btag
lowcut = [[0 for x in xrange(2)] for x in xrange(4)]
higcut = [[0 for x in xrange(2)] for x in xrange(4)]
# 0 is reg
lowcut [0] = [250., 265.]
higcut [0] = [330., 330.]
# 1 is base
lowcut [1] = [255., 250.]
higcut [1] = [330., 325.]
# 2 is regkin
lowcut [2] = [290., 285.]
higcut [2] = [315., 315.]
# 3 is kin
lowcut [3] = [290., 285.]
higcut [3] = [315., 315.]

jets = ["base", "reg", "kin", "regkin"]
#jets = ["base", "reg"]

for jet in jets:
	if jet == "base":   ireg = 0; file01 = file1; file02 = file2; ijet = 1
	if jet == "reg":    ireg = 1; file01 = file1; file02 = file2; ijet = 0
	if jet == "kin":    ireg = 0; file01 = file3; file02 = file4; ijet = 3
	if jet == "regkin": ireg = 1; file01 = file3; file02 = file4; ijet = 2
#	print jet, ireg, file01, file02
	#continue
	fileA=open(file01, 'r')
	fileB=open(file02, 'r')
	while True:
		line1 = fileA.readline()
		line1 = line1.rstrip()
		if not line1: break
#		print line1
		line2 = fileB.readline()
		line2 = line2.rstrip()
		if not line2: break
#		print line2
#		continue
		method   = line1.split("\t")[1]
		category = line1.split("\t")[0]
		if method != "CrystalBall": continue
		if category == "allcat": continue
		if category == "1btag": icat = 0
		else: icat = 1
#		if category == "1btag": continue
		cat1     = line1.split("\t")[0]
		cat2     = line2.split("\t")[0]
		mu1      = float(line1.split("\t")[2 + ireg*3])
		mu2      = float(line2.split("\t")[2 + ireg*3])
		sig1     = float(line1.split("\t")[3 + ireg*3])
		sig2     = float(line2.split("\t")[3 + ireg*3])
#		print line1
#		print line2
#		print cat1,mu1,sig1
#		print cat2,mu2,sig2
		L = higcut[ijet][icat] - lowcut[ijet][icat]
		x = lowcut[ijet][icat] + L/2
		mu = mu2-mu1
		k = sig2 / sig1
		print "for {} and category {} shift by {} and enlarge by a factor {:.3f}, [{:.1f}, {:.1f}] ---> [{:.3f}, {:.3f}] ---> [{:.1f}, {:.1f}]".format(jet, cat1, x, k, x-L/2, x+L/2, x+mu-L*k/2, x+mu+L*k/2, downround(x+mu-L*k/2), upround(x+mu+L*k/2))
	fileA.close()
	fileB.close()
	print "######################\n"
	
	

