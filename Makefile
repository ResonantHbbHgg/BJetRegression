CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
TMVA = -L${ROOTSYS}lib -lTMVA

all: prepareOpTreeInJetTree_forTraining.exe trainRegression.exe selection.exe fitMass.exe quickTrees.exe undoVarTRansformNorm.exe

prepareOpTreeInJetTree_forTraining.exe: prepareOpTreeInJetTree_forTraining.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) prepareOpTreeInJetTree_forTraining.cc -o prepareOpTreeInJetTree_forTraining.exe

trainRegression.exe: trainRegression.cc
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) trainRegression.cc -o trainRegression.exe

selection.exe: selection.cc
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) selection.cc -o selection.exe

quickTrees.exe: quickTrees.cc
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) quickTrees.cc -o quickTrees.exe

fitMass.exe: fitMass.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) $(ROOFITLIBS) fitMass.cc -o fitMass.exe

undoVarTRansformNorm.exe: undoVarTRansformNorm.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) undoVarTRansformNorm.cc -o undoVarTRansformNorm.exe

clean:
	rm *.exe
