CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
# boost
BOOSTFLAGS = -I${BOOST_ROOT}include/boost-1_48
BOOSTLIBS = -L${BOOST_ROOT}lib -lboost_program_options-gcc43-mt-1_48

TMVA = -L${ROOTSYS}lib -lTMVA

all: prepareOpTreeInJetTree_forTraining.exe trainRegression.exe selection.exe quickTrees.exe undoVarTRansformNorm.exe prepareCheckOverTraining.exe

prepareOpTreeInJetTree_forTraining.exe: prepareOpTreeInJetTree_forTraining.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) prepareOpTreeInJetTree_forTraining.cc -o prepareOpTreeInJetTree_forTraining.exe

prepareCheckOverTraining.exe: prepareCheckOverTraining.cc
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) prepareCheckOverTraining.cc -o prepareCheckOverTraining.exe

trainRegression.exe: trainRegression.cc
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) trainRegression.cc -o trainRegression.exe

TFitParticleEtEtaPhi.o: ../KinematicFit/TFitParticleEtEtaPhi.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) -c ../KinematicFit/TFitParticleEtEtaPhi.cc -o TFitParticleEtEtaPhi.o

TKinFitter.o: ../KinematicFit/TKinFitter.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) -c ../KinematicFit/TKinFitter.cc -o TKinFitter.o

TFitConstraintM.o: ../KinematicFit/TFitConstraintM.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) -c ../KinematicFit/TFitConstraintM.cc -o TFitConstraintM.o

TAbsFitParticle.o: ../KinematicFit/TAbsFitParticle.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) -c ../KinematicFit/TAbsFitParticle.cc -o TAbsFitParticle.o

DiJetKinFitter.o: ../KinematicFit/DiJetKinFitter.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) -c ../KinematicFit/DiJetKinFitter.cc -o DiJetKinFitter.o

TAbsFitConstraint.o: ../KinematicFit/TAbsFitConstraint.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) -c ../KinematicFit/TAbsFitConstraint.cc -o TAbsFitConstraint.o

selection.o: selection.cc selection.h
	$(CC) $(TMVA) $(CCFLAGS) $(ROOTFLAGS) $(BOOSTFLAGS) -c selection.cc -o selection.o

selection.exe: selection.o DiJetKinFitter.o TKinFitter.o TFitParticleEtEtaPhi.o TAbsFitParticle.o TFitConstraintM.o TAbsFitConstraint.o
	$(CC) $(TMVA)  $(ROOTLIBS) $(BOOSTLIBS) selection.o DiJetKinFitter.o TKinFitter.o TFitParticleEtEtaPhi.o TAbsFitParticle.o TFitConstraintM.o TAbsFitConstraint.o -o selection.exe

BTagUtils.o: ../h2gglobe/BTagUtils.cc ../h2gglobe/BTagUtils.h
	$(CC) $(CCFLAGS) $(ROOTFLAGS) -c ../h2gglobe/BTagUtils.cc -o BTagUtils.o

quickTrees.o: quickTrees.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(BOOSTFLAGS) -c quickTrees.cc -o quickTrees.o

quickTrees.exe: quickTrees.o BTagUtils.o
	$(CC) $(ROOTLIBS) $(BOOSTLIBS) quickTrees.o BTagUtils.o -o quickTrees.exe

fitMass.exe: fitMass.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) $(ROOFITLIBS) fitMass.cc -o fitMass.exe

undoVarTRansformNorm.exe: undoVarTRansformNorm.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) undoVarTRansformNorm.cc -o undoVarTRansformNorm.exe

clean:
	rm *.exe; rm *.o
