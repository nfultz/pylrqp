CFLAGS ?= -Wall -Wextra -std=c99 -fPIC
IFLAGS ?= -I/usr/include/python2.7 -I/usr/local/lib/python2.7/site-packages/numpy/core/include
LFLAGS ?= -lblas -llapack -lm -lpython

# ---- Link ---------------------------
_LowRankQP.so:  LowRankQP.o
	$(CC) $(CFLAGS) LowRankQP.o LowRankQPmodule.o $(LFLAGS) -shared -o _LowRankQP.so

# ---- gcc C compile ------------------
LowRankQP.o:  LowRankQP.c LowRankQPmodule.c Makefile
	$(CC) -c LowRankQP.c LowRankQPmodule.c $(IFLAGS) $(CFLAGS)


clean:
	rm *.o *.so
