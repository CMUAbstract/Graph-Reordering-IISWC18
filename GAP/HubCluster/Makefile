# See LICENSE.txt for license details.

ifneq ($(DEBUG),)
    CXX_FLAGS += -std=c++11 -O0 -g -Wall
else
    CXX_FLAGS += -std=c++11 -O3 -g -Wall
endif

ifneq ($(PERFCTR),)
    CXX_FLAGS += -pthread -DLIKWID_PERFMON -I/usr/local/include/ -L/usr/local/lib -llikwid -lm
endif

ifneq ($(SIMULATE),) 
    CXX_FLAGS += -DSIMULATE -I/home/vignesh/Graph_Specialization/Simulator/include
endif

PAR_FLAG = -fopenmp -march=native

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -openmp 
endif

ifneq (,$(findstring sunCC,$(CXX)))
	CXX_FLAGS = -std=c++11 -xO3 -m64 -xtarget=native
	PAR_FLAG = -xopenmp
endif

ifneq ($(SERIAL), 1)
	CXX_FLAGS += $(PAR_FLAG)
endif

#KERNELS = bc bfs-push pr sssp sssp-bf bfs preprocessor preprocessor-weighted cc
KERNELS = bc pr sssp bfs cc preprocessor-pr preprocessor-sssp preprocessor-bfs preprocessor-sssp 
SUITE = $(KERNELS) graph-reorderer #converter

.PHONY: all
all: $(SUITE)

% : %.cc *.h
	$(CXX) $(CXX_FLAGS) $< -o $@

# Testing
#include test/test.mk

# Benchmark Automation
#include benchmark/bench.mk
#include benchmark/bench_medium.mk


.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
