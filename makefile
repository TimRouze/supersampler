CXX=g++
CC=gcc



all: sub_sampler comparator sortCSV



#Other compilation flags
CFLAGS_BLIGHT+= -Ofast -flto -march=native -mtune=native -std=c++17 -pipe -lz -fopenmp -msse4

#Needed object files
BLO=utils.o


sub_sampler: sub_sampler.o decycling.o utils.o
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

sub_sampler.o: SubSampler.cpp SubSampler.h
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

comparator: comparator.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

comparator.o: Comparator.cpp Comparator.h
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)


decycling: decycling.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

decycling.o: Decycling.cpp Decycling.h
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

utils.o: utils.cpp $(INC) utils.h
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

sortCSV: sortCSV.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

sortCSV.o: sort_csv.cpp
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

clean:
	rm -rf *.o



rebuild: clean all
