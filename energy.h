#pragma once

#include "cLibsBundle.h"
#include "MCIcore.h"
#include "pcg_basic.h"   //RNG
#include "utility.h"
#include "dataStructs.h"

/* energy */
double energy_ising_full(lattice config, su s, double* E, para parameters);
double energy_ising_su(lattice config, su s, para parameters);
double energy_ising_su_local(su s, para parameters, lattice* subset, int lsubset);
double energy_full(lattice config, su s, double* E, double* ewald, para parameters);
double energy_local(su s, double* E, para parameters, lattice* subset, int lsubset);
double energy_lr(lattice c, double* ewald);
int getEwaldIndexForRij(int i, int j);
double energy_lr_nflip(lattice c, double* ewald, uint32_t* k, int nk);
double energy_umbr(umb umbrellaSprings);
