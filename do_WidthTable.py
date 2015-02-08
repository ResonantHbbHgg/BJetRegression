import operator

class result(object):
    def __init__(self, method, mass, var, cat, width, mean, res, regression=False):
        self.method = method
        self.mass = mass
        self.var = var
        self.cat = cat
        self.width = width
        self.mean = mean
        self.res = res
        self.regression = regression

    def __str__(self):
        return "method= " + str(self.method) + "\tmass= " + str(self.mass) + "\tvar= " + str(self.var) + "\tcat= " + str(self.cat) + "\twidth= " + str(self.width) + "\tmean= " + str(self.mean) + "\tres= " + str(self.res) + "\tregression= ", str(self.regression)
    def __getitem__(self, i):
        if i == 0: return self.method
        if i == 1: return self.mass
        if i == 2: return self.var
        if i == 3: return self.cat
        if i == 4: return self.width
        if i == 5: return self.mean
        if i == 6: return self.res
        if i == 7: return self.regression

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
            ires = l.index('res=') + 1
            res = result(l[imethod], int(l[imass]), l[ivar], l[icat], float(l[iwidth]), float(l[imean]), float(l[ires]))
            allresults.append(res)
    return allresults

def getFromPerformanceSummary(mass_, var_):
    filename = "performanceSummary_m" + str(mass_) + "_" + var_ + ".txt"
    allresults = []
    with open(filename) as file:
        for line in file:
            l = line.split()
            # Skip the first header line
            if "Category" in l[0]:
                continue
            # Fill no regression first
            cat_ = l[0]
            method_ = l[1]
            mean_ = float(l[2])
            width_ = float(l[3])
            res_ = float(l[4])
            res = result( method_, mass_, var_, cat_, width_, mean_, res_, False)
            allresults.append(res)
            # With regression
            mean_ = float(l[5])
            width_ = float(l[6])
            res_ = float(l[7])
            res = result( method_, mass_, var_, cat_, width_, mean_, res_, True)
            allresults.append(res)
    return allresults
    
def filterResults(allResults, isRegression=False):
    allResults = [res for res in allResults if res.cat != "allcat" and res.regression == isRegression ]
#    allResults = [res for res in allResults if res.var != "mgg" and res.method != "CrystalBall" ]
#    allResults.sort(cmp=lambda a,b: cmp(a.mass, b.mass))
    return allResults
            
def getResultsForMass(allResults, mass):
    return [res for res in allResults if res.mass == mass]    

def sortNicely(allResults):
    # var, cat, method
    allResults.sort(key=operator.itemgetter(2,3,0))
    return allResults

allResults = []
for mass in map(int, "270 300 350 400".split()):
    for spectrum in ["mjj", "mggjj"]:
        for res in getFromPerformanceSummary(mass, spectrum) :
            allResults.append( res )

print "without regression"
filteredResults = filterResults(allResults, False)
for mass in map(int, "270 300 350 400".split()):
    line = str(mass)
    for res in sortNicely(getResultsForMass(filteredResults, mass)):
        line += " & "
        line += str(round(res.res, 2))
#        print res
    line += "\\\\"
    print line

print "with regression"
filteredResults = filterResults(allResults, True)
for mass in map(int, "270 300 350 400".split()):
    line = str(mass)
    for res in sortNicely(getResultsForMass(filteredResults, mass)):
        line += " & "
        line += str(round(res.res, 2))
#        print res
    line += "\\\\"
    print line

#filename = "resolution_SigmaEff_all.txt"
#allResults = getAllResults(filename)
#filteredResults = filterResults(allResults)
##for mass in map(int, "260 270 300 350 400 450 500 550 600 650 700 800 900 1000 1100".split()):
#for mass in map(int, "270 300 350 400".split()):
#    line = str(mass)
#    for res in sortNicely(getResultsForMass(filteredResults, mass)):
#        line += " & "
#        line += str(round(res.width, 2))
##        print res
#    line += "\\\\"
#    print line
#
