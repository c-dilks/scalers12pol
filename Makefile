
all: sca_read_bin.o

sca_read_bin.o :	sca_read_bin.c
	gcc read_scalers.c sca_read_bin.c -o sca_read_bin.o

clean:	
	/bin/rm -f sca_read_bin.o
