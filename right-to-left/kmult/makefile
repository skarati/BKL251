INCDRS = -I include/

SRCFLS = source/gfe1xAdd.S 		\
	 source/gfe1xSq.S 		\
	 source/gfe1xMult.S  		\
	 source/gfe1xMultConst.S 	\
	 source/conditionalSwap.S 	\
	 source/ladderStep.S 	\
	 binKLmain.c
         
OBJFLS = source/gfe1xAdd.o 		\
	 source/gfe1xSq.o 		\
	 source/gfe1xMult.o  		\
	 source/gfe1xMultConst.o 	\
	 source/conditionalSwap.o 	\
	 source/ladderStep.o		\
	 binKLmain.o

EXE    = binKLmain

#CFLAGS = -march=haswell -mtune=haswell -mavx2 -m64 -O3 -funroll-loops -fomit-frame-pointer
CFLAGS = -march=skylake -mtune=skylake -m64 -mavx2 -msse4.2 -mpclmul -no-pie -O3 -fomit-frame-pointer

CC     = gcc
LL     = gcc

$(EXE): $(OBJFLS)
	$(LL) -o $@ $(OBJFLS) -lm -no-pie

.c.o:
	$(CC) $(INCDRS) $(CFLAGS) -o $@ -c $<

clean:
	-rm $(EXE)
	-rm $(OBJFLS)
