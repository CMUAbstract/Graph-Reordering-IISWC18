ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

ifdef PD
PD = -DPD
endif

ifdef BYTE
CODE = -DBYTE
else ifdef NIBBLE
CODE = -DNIBBLE
else
CODE = -DBYTERLE
endif

ifdef LOWMEM
MEM = -DLOWMEM
endif

ifdef ALIGNED
ALIGN = -DALIGNED
endif

OPENMP=1

#compilers
ifdef CILK
PCC = g++
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -O2 -DCILK $(INTT) $(INTE) $(CODE) $(PD) $(MEM)
PLFLAGS = -fcilkplus -lcilkrts

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -std=c++14 -O3 -DCILKP $(INTT) $(INTE) $(CODE) $(PD) $(MEM)

else ifdef OPENMP
PCC = g++
PCFLAGS = -std=c++14 -fopenmp -march=native -O3 -DOPENMP $(INTT) $(INTE) $(CODE) $(PD) $(MEM) $(ALIGN)
#PCFLAGS = -std=c++14 -fopenmp -march=native -O0 -g -DOPENMP $(INTT) $(INTE) $(CODE) $(PD) $(MEM)

else
PCC = g++
PCFLAGS = -std=c++14 -O2 $(INTT) $(INTE) $(CODE) $(PD) $(MEM)
endif

COMMON= ligra.h graph.h compressedVertex.h vertex.h utils.h IO.h parallel.h gettime.h timer.h index_map.h maybe.h sequence.h edgeMap_utils.h binary_search.h quickSort.h blockRadixSort.h transpose.h parseCommandLine.h byte.h byteRLE.h nibble.h byte-pd.h byteRLE-pd.h nibble-pd.h vertexSubset.h pvector.h encoder.C decoder.C sliding_queue.h platform_atomics.h


#ALL= encoder decoder BFS BC BellmanFord Components Radii PageRank PageRankDelta BFSCC BFS-Bitvector KCore MIS Triangle CF
ALL= BFS BC BC-4iter BellmanFord BellmanFord-push Radii PageRank PageRankDelta preprocessor preprocessor-weighted preprocessor-undirected KCore CF Components MIS preprocessor-densewrite preprocessor-bipartite

all: $(ALL)

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

$(COMMON):
	ln -s ../ligra/$@ .

.PHONY : clean

clean :
	rm -f *.o $(ALL)

cleansrc :
	rm -f *.o $(ALL)
	rm $(COMMON)
