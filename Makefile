
# ---- Link ---------------------------
_LowRankQP.so:  LowRankQP.o
	gcc -shared -fPIC LowRankQP.o LowRankQPmodule.o -lblas -llapack -o _LowRankQP.so

# ---- gcc C compile ------------------
LowRankQP.o:  LowRankQP.c LowRankQPmodule.c Makefile
	gcc -shared -fPIC -c LowRankQP.c LowRankQPmodule.c -I/usr/include/python2.7 -I/usr/lib/python2.7/dist-packages/numpy/core/include/numpy
