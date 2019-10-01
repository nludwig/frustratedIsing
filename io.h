#pragma once

#include "cLibsBundle.h"
#include "dataStructs.h"
#include "switches.h"
#include "shapes.h"
#include "utility.h"

void dumpParameters(para p, su s, umb umbr, int coreNum, FILE* f, FILE* myerr);
int* readSetParameters(para p, char**** linesByType, int coreNum, FILE* f, FILE* myerr);
void dumpData(lattice config, su solutes, FILE* out);
int readData(char*** lines, FILE* f);
void loadSu(char** lines, int nlines, su* s, su* ts, char** soluteNameHolder, double** solComHolder, \
                 double** unwrappedSolComHolder, double** orientationHH, double*** orientationH,      \
                 double** solRelPos, double*** solRelPosPtrs, double** shlRelPos,                      \
                 double*** shlRelPosPtrs, lattice** solCurPos, lattice** shlCurPos, bool** hydrH, FILE* myerr);
void dumpSuCom(su s, int i, FILE* f);
void dumpSu(su s, FILE* f);
void dumpAllSuComs(su s, FILE* f);
void loadUmbrellas(char** lines, int nlines, umb umbrellaSprings, su solutes, para parameters, FILE* myumbout, FILE* myerr);
void dumpUmbrellas(umb u, su s, para parameters, FILE* f, FILE* myerr);
int readLines(char* header, char*** soluteLines, char** inLines, int ninLines, FILE* inF, FILE* myerr);
void dumpLammps(lattice config, su solutes, int ts, int coreNum, FILE* out, FILE* myerr);
void dumpFI(lattice config, su solutes, int ts, FILE* out, FILE* myerr);
