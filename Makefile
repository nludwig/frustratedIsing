CC = icc
CCFLAGS = -O3 -xHost -std=c99 -qopenmp
#CC = gcc
#CCFLAGS = -O3 -march=native -mtune=native -std=c99 -fopenmp
OBJ = main.o MCIcore.o io.o pcg_basic.o utility.o shapes.o energy.o forRandOrder.o
LFLAGS = -lfftw3 -lm

MCIsingFI : $(OBJ)
	$(CC) $(CCFLAGS) $(OBJ) $(LFLAGS) -o MCIsingFIv27
main.o : main.c headerBundle.h
	$(CC) $(CCFLAGS) -c main.c $(LFLAGS)
MCIcore.o : MCIcore.c MCIcore.h
	$(CC) $(CCFLAGS) -c MCIcore.c $(LFLAGS)
io.o : io.c dataStructs.h io.h utility.h utility.o
	$(CC) $(CCFLAGS) -c io.c $(LFLAGS)
pcg_basic.o : pcg_basic.c pcg_basic.h
	$(CC) $(CCFLAGS) -c pcg_basic.c $(LFLAGS)
utility.o : utility.c utility.h
	$(CC) $(CCFLAGS) -c utility.c $(LFLAGS)
shapes.o : shapes.c shapes.h
	$(CC) $(CCFLAGS) -c shapes.c $(LFLAGS)
energy.o : energy.c energy.h
	$(CC) $(CCFLAGS) -c energy.c $(LFLAGS)
forRandOrder.o : forRandOrder.c forRandOrder.h
	$(CC) $(CCFLAGS) -c forRandOrder.c $(LFLAGS)
clean :
	rm MCIsingFIv27 $(OBJ)
