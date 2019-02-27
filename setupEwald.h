#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <assert.h>

//compile w/ icc for ~50x speedup

struct parameters {
        double sigma;
};
typedef struct parameters* para;

void getPos(int n, int* pos, int* L_in);
double pbc(double dst, int e, int* L);
void diff_lat(double* dr, unsigned long i, unsigned long j, int* L);
void setupEwald(double*** ewald, int* L, para p);
void setupEwald_omp(double*** ewald, int* L, para p);
void dumpEwald(double** ewald, int* L, FILE* f);
void loadEwald(double*** ewald, int N, char** lines, int nlines);
