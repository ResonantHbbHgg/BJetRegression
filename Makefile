CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
TMVA = -L${ROOTSYS}lib -lTMVA

all: prepareOpTreeInJetTree_forTraining.exe trainRegression.exe selection.exe fitMass.exe quickTrees.exe

prepareOpTreeInJetTree_forTraining.exe: prepareOpTreeInJetTree_forTraining.C
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) prepareOpTreeInJetTree_forTraining.C -o prepareOpTreeInJetTree_forTraining.exe

trainRegression.exe: trainRegression.C
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) trainRegression.C -o trainRegression.exe

selection.exe: selection.C
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) selection.C -o selection.exe

quickTrees.exe: quickTrees.C
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) quickTrees.C -o quickTrees.exe

fitMass.exe: fitMass.C
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) $(ROOFITLIBS) fitMass.C -o fitMass.exe

clean:
	rm *.exe
