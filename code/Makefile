# Compiler settings
CXX = g++
NVCC = nvcc
CXXFLAGS = -std=c++17 -Wall -Wextra -Iinclude -fopenmp 
LDFLAGS = -lgomp  
# Source files
MAIN_SRC = main.cpp
SEQ_SRC = src/sequential.cpp
CUDA_SRC = src/graph_cuda.cu
OMP_SRC = src/graph_omp.cpp


# Object files
MAIN_OBJ = $(MAIN_SRC:.cpp=.o)
SEQ_OBJ = $(SEQ_SRC:.cpp=.o)
CUDA_OBJ = $(CUDA_SRC:.cu=.o)
OMP_OBJ = $(OMP_SRC:.cpp=.o)

# Output binary
BIN = maxflow

# Default target
all: $(BIN)

# Link all object files into one binary
$(BIN): $(MAIN_OBJ) $(SEQ_OBJ) $(CUDA_OBJ) $(OMP_OBJ)
	$(NVCC) -std=c++17 -o $@ $^ $(LDFLAGS)

# Compile C++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA files
%.o: %.cu
	$(NVCC) -Xcompiler -fPIC -std=c++17 -Iinclude -c $< -o $@

# Clean up
clean:
	rm -f *.o $(BIN)

.PHONY: all clean
