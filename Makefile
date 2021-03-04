#修改优化选项会影响测试结果
CC = gcc
CFLAGS = -lm -O2

FF = gfortran
FFLAGS = -O2

all:stream2 stream2_mu

stream2:stream2.c
	$(CC) $(CFLAGS) -o $@ $^

stream2_mu:stream2.c
	$(CC) $(CFLAGS) -fopenmp -o $@ $^

stream2_f:stream2.f mysecond.c
	$(CC) $(CFLAGS) -c mysecond.c
	$(FF) $(FFLAGS) stream2.f mysecond.o -o $@

clean:
	rm -f stream2_mu stream2 *.o
	rm -f stream2_f fort.3