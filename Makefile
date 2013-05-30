CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
TMVA = -L${ROOTSYS}lib -lTMVA

all: prepareOpTree.exe trainRegression.exe

prepareOpTree.exe: prepareOpTree.C
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) prepareOpTree.C -o prepareOpTree.exe

trainRegression.exe: trainRegression.C
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) trainRegression.C -o trainRegression.exe

