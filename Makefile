FLAGS = -std=c99 -Wall

all: complex.c fft.c test.c
	gcc complex.c $(FLAGS) -c
	gcc fft.c $(FLAGS) -c 
	gcc test.c complex.o fft.o $(FLAGS) -o test -lm
	gcc benchmark.c complex.o fft.o $(FLAGS) -o benchmark -lm

clean:
	rm complex.o fft.o test benchmark
