#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <fftw3.h>
#include <assert.h>
#include <stdbool.h>
#include "dataStructs.h"
#include "utility.h"

int** buildCube(int l);
int** buildRectangle(int* l);
int* getShapeInterface(int** shape, int nsites, int* ninterface, int** Lcube, int* border);
int getSA(int** shape, int nsites, int** Lcube, int* border);
int** getShapeInterfacePos(int** shape, int nsites, int* ninterface, int** Lcube, int* border);
int uniqueint(int** arr, int larr);
int uniqueuint32_t(uint32_t** arr, int larr);
int uniquelattice(lattice** arr, int larr);
int* getShapeBounds(int** shape, int nsites);
int* embedShapeInBorder(int** shape, int nsites, int* border);
int** buildNeighCube(double Rc);
int diff(int** arr1, int larr1, int** arr2, int larr2);
//int diff(int* arr1, int larr1, int* arr2, int larr2);
double* matrixOnVector(double** M, double* v);
void matrixOnVectorInPlace(double** M, double* vout);
double** matrixOnMatrix(double** M1, double** M2);
void matrixOnMatrixInPlace(double** M, double** Mout);
void copyMatrixInPlace(double** M, double** Mout, int d);
double** copyMatrix(double** M, int d);
double*** buildPiO2RotationMatricies(int d);
