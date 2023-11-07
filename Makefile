CC=gcc
MPICC=mpicc

FLAG=-g -O3

gbm2d: alloc.o atopkge.o complex.o errpkge.o fft.o franuni.o getpars.o maingray.o SegyRead.o sinc.o MpiFrame.o syscalls.o subcalls.o split.o
	$(MPICC) $(FLAG) -o gbm2d MpiFrame.o alloc.o atopkge.o complex.o errpkge.o fft.o franuni.o getpars.o maingray.o SegyRead.o sinc.o syscalls.o subcalls.o split.o -lm


MpiFrame.o: MpiFrame.c
	$(MPICC) $(FLAG) -c MpiFrame.c

alloc.o: alloc.c
	$(CC) $(FLAG) -c alloc.c

atopkge.o: atopkge.c
	$(CC) $(FLAG) -c atopkge.c

complex.o: complex.c
	$(CC) $(FLAG) -c complex.c

errpkge.o: errpkge.c
	$(CC) $(FLAG) -c errpkge.c

fft.o: fft.c
	$(CC) $(FLAG) -c fft.c

franuni.o: franuni.c
	$(CC) $(FLAG) -c franuni.c

getpars.o: getpars.c
	$(CC) $(FLAG) -c getpars.c

maingray.o: maingray.c
	$(CC) $(FLAG) -c maingray.c

SegyRead.o: SegyRead.c
	$(CC) $(FLAG) -c SegyRead.c

sinc.o: sinc.c
	$(CC) $(FLAG) -c sinc.c

syscalls.o: syscalls.c
	$(CC) $(FLAG) -c syscalls.c

subcalls.o: subcalls.c
	$(CC) $(FLAG) -c subcalls.c

split.o: split.c
	$(CC) $(FLAG) -c split.c

accimage.o: accimage.c
	$(CC) $(FLAG) -c accimage.c

clean:
	rm -f *.o gbm2d







