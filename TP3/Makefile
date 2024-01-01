CC = gcc
CFLAGS = -O3 -ftree-vectorize -msse4 -mavx -march=native -funroll-loops -finline-functions 

SHELL = /bin/sh
BIN_NAME = cuda

CXX = nvcc
LD  = nvcc

CXXFLAGS = -O3 -g -std=c++11 -arch=sm_35 -Wno-deprecated-gpu-targets

SRC_DIR = src/
BIN_DIR = bin
BUILD_DIR = build
SRC = $(wildcard $(SRC_DIR)*.cu)
OBJ = $(patsubst src/%.cu,build/%.o,$(SRC))
BIN = $(BIN_NAME)

vpath %.cu $(SRC_DIR)

.DEFAULT_GOAL = all

$(BUILD_DIR)/%.o: %.cu
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

$(BIN_DIR)/$(BIN_NAME): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJ) $(LIBS)

checkdirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

all: checkdirs $(BIN_DIR)/$(BIN_NAME) MDseq.exe MDpar.exe

MDseq.exe: $(SRC_DIR)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS)-pg  $(SRC_DIR)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC_DIR)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS)-pg $(SRC_DIR)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm -rf slurm* MDpar.exe MDseq.exe gmon.out cp_average.txt cp_output.txt cp_traj.xyz lixo relatorio.txt

load:
	module load cuda/11.3.1
	module load gcc/7.2.0

run:
	sbatch run.sh

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	sbatch test.sh

runseqWperf: MDseq.exe
	srun --partition=cpar perf stat ./MDseq.exe < inputdata.txt

runparWperf: MDpar.exe
	@export OMP_NUM_THREADS=40; srun --partition=cpar perf stat ./MDpar.exe < inputdata.txt

test: 
	@diff -s cp_average.txt Correct_res/cp_average.txt
	@diff -s cp_traj.xyz Correct_res/cp_traj.xyz
	@diff -s --suppress-common-lines cp_output.txt Correct_res/cp_output.txt 