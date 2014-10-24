import operator

class result(object):
	def __init__(self, method, mass, var, cat, width, mean):
		self.method = method
		self.mass = mass
		self.var = var
		self.cat = cat
		self.width = width
		self.mean = mean
	def __str__(self):
		return "method= " + str(self.method) + "\tmass= " + str(self.mass) + "\tvar= " + str(self.var) + "\tcat= " + str(self.cat) + "\twidth= " + str(self.width) + "\tmean= " + str(self.mean)
	def __getitem__(self, i):
		if i == 0: return self.method
		if i == 1: return self.mass
		if i == 2: return self.var
		if i == 3: return self.cat
		if i == 4: return self.width
		if i == 5: return self.mean

def getAllResults(filename):
	allresults = []
	with open(filename) as file:
		for line in file:
			l = line.split()
			imethod = l.index('method=') + 1
			imass = l.index('mass=') + 1
			ivar = l.index('var=') + 1
			icat = l.index('cat=') + 1
			iwidth = l.index('width=') + 1
			imean = l.index('mean=') + 1
			res = result(l[imethod], int(l[imass]), l[ivar], l[icat], float(l[iwidth]), float(l[imean]))
			allresults.append(res)
	return allresults
	
def filterResults(allResults):
	allResults = [res for res in allResults if res.var != "mgg" and res.method != "CrystalBall" ]
#	allResults.sort(cmp=lambda a,b: cmp(a.mass, b.mass))
	return allResults
			
def getResultsForMass(allResults, mass):
	return [res for res in allResults if res.mass == mass]	

def sortNicely(allResults):
# var, cat, method
	allResults.sort(key=operator.itemgetter(2,3,0))
	return allResults

filename = "resolution_SigmaEff_all.txt"
allResults = getAllResults(filename)
filteredResults = filterResults(allResults)
for mass in map(int, "260 270 300 350 400 450 500 550 600 650 700 800 900 1000 1100".split()):
	line = str(mass)
	for res in sortNicely(getResultsForMass(filteredResults, mass)):
		line += " & "
		line += str(round(res.width, 2))
#		print res
	line += "\\\\"
	print line

