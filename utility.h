#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include "dataStructs.h"

/* general */
char** split(char* delimStr, const char* delimiter, int* numWords);
void swapD(double* a, double* b);
void swapMultD(double** v, int narr, int L, int a, int b);
void qsortMultD(double** v, int narr, int L, int left, int right, int (*comp)(double*,double*));
void freeAll(void** ptrs, int nptrs);
uint16_t ndigits(long int number);
int iCmp(int* a, int* b);
int lattCmp(lattice* a, lattice* b);
int lattCmpRev(lattice* a, lattice* b);
int uint32_tCmp(uint32_t* a, uint32_t* b);
int doubleMatrixCmp(double** a, double**b, int d); //0 if diff, 1 if same
int intArrMin(int* arr, int larr);
int intArrMax(int* arr, int larr);

/* lattice related */
void getPos(int n, int* pos, int* L_in);
void getPosLat(lattice config, lattice site, int* pos);
int getSite(int* pos, int* L_in);
int getSitePBC(int* pos, int* L_in);
lattice getSiteLat(lattice config, int* pos);
int getNN(int part, int e, int* L_in);
int getNNPBC(int part, int e, int* L_in);
lattice getNNLat(lattice config, int part, int e);
void setCharge(lattice c);
double* checkCharge(lattice c, int lc, double* ret);
//void testRng(int numRNs, pcg32_random_t* rng);
void initPara(para p);
para copyPara(para p);

/* bc related */
double pbc(double dist, int vectDir);
void wrapIntoL(int* pos, int* L_in);
void wrapIntoL_double(double* pos, int* L_in);
double dist_intvect_nopbc(int* vector0, int* vector1);
double dist_vect(double* vector0, double* vector1);
void diff_lat(double* dr, int part_i, int part_j);
double dist_lat(int part_i, int part_j);
