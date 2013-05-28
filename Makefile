CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats

all: prepareOpTree.exe

prepareOpTree.exe: prepareOpTree.C
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) prepareOpTree.C -o prepareOpTree.exe



