CC = gcc
SRC = src/
CFLAGS = -O3 -fno-omit-frame-pointer -ffast-math -funroll-loops -ftree-vectorize -msse4 -mavx -pg
LIBS = -lm

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp $(LIBS) -o MD.exe


clean:
	rm -rf MD.exe gmon.out cp_average.txt cp_output.txt cp_traj.xyz

run:
	srun --partition=cpar perf stat ./MD.exe < inputdata.txt

run1:
	perf stat ./MD.exe < inputdata.txt

test: 
	@diff -s cp_average.txt Correct_res/cp_average.txt
	@diff -s cp_traj.xyz Correct_res/cp_traj.xyz
	@diff -s --suppress-common-lines cp_output.txt Correct_res/cp_output.txt 