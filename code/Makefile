# Compiler settings
CXX = g++
NVCC = nvcc
CXXFLAGS = -std=c++17 -Wall -Wextra -Iinclude
NVCCFLAGS = -O3 --cudart shared -rdc=true

# Source files
MAIN_SRC = main.cpp
SEQ_SRC = src/sequential.cpp
CUDA_SRC = src/graph_cuda.cu

# Object files
MAIN_OBJ = $(MAIN_SRC:.cpp=.o)
SEQ_OBJ = $(SEQ_SRC:.cpp=.o)
CUDA_OBJ = $(CUDA_SRC:.cu=.o)

# Output binary
BIN = maxflow

# Default target
all: $(BIN)

# Link all object files into one binary
$(BIN): $(MAIN_OBJ) $(SEQ_OBJ) $(CUDA_OBJ)
	$(NVCC) -std=c++17 -rdc=true -o $@ $^

# Compile C++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA files
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -std=c++17 -Iinclude -c $< -o $@

# Clean up
clean:
	rm -f *.o $(BIN)

.PHONY: all clean
