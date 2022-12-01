# quick makefile for testing coloring algorithm
#
# make sure to have opencilk's clang compiler 
# and if you do, check the directory path
# 
# run any executable with an mtx file as argument
# to imitate results

CC = gcc
CLANG = /opt/opencilk/bin/clang
GFLAGS = -O3 -w -g
CFLAGS = -fopencilk $(GFLAGS)
MPFLAGS = -fopenmp $(GFLAGS)
PFLAGS = -pthread $(GFLAGS)
DEPS = mmio.h
OBJ = mmio.o
OBJSEQ = ColoringSCCSequential.o
OBJC = ColoringSCCopenCilk.o
OBJMP = ColoringSCCopenMP.o
OBJP = ColoringSCCpthreads.o

all: ColoringSCCSequential ColoringSCCopenMP ColoringSCCpthreads ColoringSCCopenCilk

$(OBJ) : mmio.c $(DEPS)
	$(CC) -c -o $@ $< $(GFLAGS)

$(OBJSEQ) : ColoringSCCSequential.c $(DEPS)
	$(CC) -c -o $@ $< $(GFLAGS)

$(OBJMP) : ColoringSCCopenMP.c $(DEPS)
	$(CC) -c -o $@ $< $(MPFLAGS)

$(OBJP) : ColoringSCCpthreads.c $(DEPS)
	$(CC) -c -o $@ $< $(PFLAGS)

$(OBJC) : ColoringSCCopenCilk.c $(DEPS)
	$(CLANG) -c -o $@ $< $(CFLAGS)

ColoringSCCSequential : $(OBJ) $(OBJSEQ) 
	$(CC) -o $@ $^ $(GFLAGS)

ColoringSCCopenMP: $(OBJ) $(OBJMP)
	$(CC) -o $@ $^ $(MPFLAGS)

ColoringSCCpthreads: $(OBJ) $(OBJP)
	$(CC) -o $@ $^ $(PFLAGS)

ColoringSCCopenCilk: $(OBJ) $(OBJC)
	$(CLANG) -o $@ $^ $(CFLAGS)

run:
	@./ColoringSCCSequential ./mtxfiles/celegansneural.mtx
	@echo ""
	@./ColoringSCCopenCilk ./mtxfiles/celegansneural.mtx
	@echo ""
	@./ColoringSCCopenMP ./mtxfiles/celegansneural.mtx
	@echo ""
	@./ColoringSCCpthreads ./mtxfiles/celegansneural.mtx

.PHONY : clean
clean : 
	@rm ColoringSCCSequential ColoringSCCopenCilk ColoringSCCopenMP ColoringSCCpthreads $(OBJ) $(OBJSEQ) $(OBJC) $(OBJMP) $(OBJP)