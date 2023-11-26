CC = gcc
SRC = src/
CFLAGS = -O3 -fopenmp -ftree-vectorize -msse4 -mavx -pg 


.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm -rf slurm* MDpar.exe MDseq.exe gmon.out cp_average.txt cp_output.txt cp_traj.xyz lixo relatorio.txt

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=16;\
	./MDpar.exe < inputdata.txt

runseqperf: MDseq.exe
	 perf stat ./MDseq.exe < inputdata.txt

runparperf: MDpar.exe
	export OMP_NUM_THREADS=16;\
	 perf stat ./MDpar.exe < inputdata.txt

test: 
	@diff -s cp_average.txt Correct_res/cp_average.txt
	@diff -s cp_traj.xyz Correct_res/cp_traj.xyz
	@diff -s --suppress-common-lines cp_output.txt Correct_res/cp_output.txt 