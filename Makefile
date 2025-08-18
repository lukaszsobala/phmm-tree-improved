# Makefile for PHMM-Tree with OpenMP parallelization
# Automatically detects available compilers and enables OpenMP

# Compiler detection
CC  := $(shell command -v gcc 2>/dev/null  || command -v clang 2>/dev/null  || echo cc)
CXX := $(shell command -v g++ 2>/dev/null  || command -v clang++ 2>/dev/null || echo c++)

# Base flags
CFLAGS_BASE   = -O3 -Wall -std=c99
CXXFLAGS_BASE = -O3 -Wall -std=gnu++17
LDFLAGS_BASE  = -lm

# OpenMP detection and flags (check with C++ compiler because we link C++)
OPENMP_FLAG := $(shell $(CXX) -fopenmp -E - </dev/null >/dev/null 2>&1 && echo "-fopenmp" || echo "")

ifeq ($(OPENMP_FLAG),-fopenmp)
CFLAGS   = $(CFLAGS_BASE) $(OPENMP_FLAG) -DOPENMP_ENABLED
CXXFLAGS = $(CXXFLAGS_BASE) $(OPENMP_FLAG) -DOPENMP_ENABLED
LDFLAGS  = $(LDFLAGS_BASE) $(OPENMP_FLAG)
$(info OpenMP support detected - enabling parallel compilation)
else
CFLAGS   = $(CFLAGS_BASE)
CXXFLAGS = $(CXXFLAGS_BASE)
LDFLAGS  = $(LDFLAGS_BASE)
$(info OpenMP not available - compiling sequential version)
endif

# Thread detection at compile time
CFLAGS   += -DAUTO_THREAD_DETECTION
CXXFLAGS += -DAUTO_THREAD_DETECTION

# Source files

SOURCES = fitch.c phylip.c dist.c neighbor.c upgma.c kitsch.c \
		  class_functions.cpp process_matrices.cpp public_functions.cpp \
		  process_alignments.cpp process_hmms.cpp process_usearch.cpp process_prc.cpp \
		  HMMTree.cpp phylip_draw_tree.cpp hhsuite.cpp

OBJECTS = $(SOURCES:.c=.o)
OBJECTS := $(OBJECTS:.cpp=.o)

# Target executable
TARGET = phmm-tree

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# Compile C files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Compile C++ files  
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJECTS) $(TARGET)

# Test compilation with different thread counts
test-parallel: $(TARGET)
	@echo "Testing with 1 thread:"
	OMP_NUM_THREADS=1 ./$(TARGET) || true
	@echo "Testing with 2 threads:"
	OMP_NUM_THREADS=2 ./$(TARGET) || true
	@echo "Testing with 4 threads:"
	OMP_NUM_THREADS=4 ./$(TARGET) || true

# Show compilation info
info:
	@echo "Compiler: $(CC)"
	@echo "OpenMP flag: $(OPENMP_FLAG)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"

# Force rebuild
rebuild: clean all

.PHONY: all clean test-parallel info rebuild
