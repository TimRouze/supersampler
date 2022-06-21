CXX=g++
CC=gcc



#Here snippet should be your program
all: sub_sampler comparator
#bench_blight is only here for benchmark purpose and can be removed



#Your compilation flags
CFLAGS_CUSTOM= -O9000  -std=c++17

#Other compilation flags
CFLAGS_LZ4+= -w -Wall -std=gnu99 -DUSE_THREADS  -fstrict-aliasing -Iext $(DEFS)
CFLAGS_BLIGHT+= -DNDEBUG -Ofast -flto -march=native -mtune=native -g -std=c++17 -pipe -lz -fopenmp -msse4 -Ilz4 -pg

#Needed object files
LZ4O=lz4/lz4frame.o lz4/lz4.o lz4/xxhash.o lz4/lz4hc.o
BLO=blight.o utils.o


sub_sampler: sub_sampler.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

sub_sampler.o: SubSampler.cpp
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

comparator: comparator.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

comparator.o: Comparator.cpp
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

#Blight object files (BLO)
utils.o: utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

blight.o: blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)
	
	
	
#Lz4 object files (LZ4O)
lz4/lz4frame.o: lz4/lz4frame.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)

lz4/lz4.o: lz4/lz4.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)

lz4/lz4hc.o: lz4/lz4hc.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)
	
lz4/xxhash.o: lz4/xxhash.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4) 




clean:
	rm -rf *.o
	


rebuild: clean all
