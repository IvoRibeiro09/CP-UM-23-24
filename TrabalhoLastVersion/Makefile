CC = gcc
SRC = src/
CFLAGS = # none

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -O0 -pg -lm -o buggfixer.exe
	$(CC) $(CFLAGS) $(SRC)MD.cpp -funroll-loops -O2 -pg -lm -o current.exe
	$(CC) $(CFLAGS) $(SRC)MD.cpp -funroll-loops -O2 -DDONT_USE_VOL -lm -o extraflag.exe
	$(CC) $(CFLAGS) $(SRC)MD.cpp -funroll-loops -Ofast -DDONT_USE_VOL -lm -o bernardo.exe
	$(CC) $(CFLAGS) $(SRC)MD.cpp -funroll-loops -Ofast -lm -o bernardo2.exe

clean:
	rm ./buggfixer.exe ./Dtvol.exe ./current.exe

run:
	./MD.exe < inputdata.txt
