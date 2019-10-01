#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <assert.h>

double* setupEwald(int* L, double sigma);
void dumpEwald(double* ewald, int* L, FILE* f);
int** setupPos();
void getPos(int n, int* pos);
void minImageVect(int* pos, int* L_in);
