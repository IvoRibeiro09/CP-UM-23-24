CC = gcc
SRC = src/
CFLAGS = -O3 -ftree-vectorize -msse4 -mavx -march=native -funroll-loops -finline-functions 

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS)-pg  $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS)-pg $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm -rf slurm* MDpar.exe MDseq.exe gmon.out cp_average.txt cp_output.txt cp_traj.xyz lixo relatorio.txt

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=40;\
	./MDpar.exe < inputdata.txt

runseqWperf: MDseq.exe
	srun --partition=cpar perf stat ./MDseq.exe < inputdata.txt

runparWperf: MDpar.exe
	@export OMP_NUM_THREADS=40; srun --partition=cpar perf stat ./MDpar.exe < inputdata.txt

test: 
	@diff -s cp_average.txt Correct_res/cp_average.txt
	@diff -s cp_traj.xyz Correct_res/cp_traj.xyz
	@diff -s --suppress-common-lines cp_output.txt Correct_res/cp_output.txt 