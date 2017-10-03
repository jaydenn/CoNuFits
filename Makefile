CC = mpif90
FLAGS = -O3 -DMPI -DOMPI_SKIP_MPICXX
LIBS = -L../MultiNest_v3.9 -lnest3 -lmpi -llapack -lgfortran -L/usr/local/lib -lgsl -lgslcblas -lstdc++ -lm
INCLUDE = -I./src/inc
OBJECTS = src/ConuFits.o src/fileIO.o src/calcRates.o src/nuFlux.o src/detectorFunctions.o src/formfactorSI.o src/monteCarlo.o src/likelihood.o src/nuRate.o src/SMrate.o src/BSMrate.o src/sterileOscillation.o src/sterileRate.o src/NSIrate.o src/NSIglobalFit.o

default: ConuFits

ConuFits: $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

src/%.o: src/%.cpp
	$(CC) $< $(FLAGS) ${INCLUDE} -c -o $@

clean:
	-rm src/*.o
	-rm -f ./ConuFits
