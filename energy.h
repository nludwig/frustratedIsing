#pragma once

#include "cLibsBundle.h"
#include "MCIcore.h"
#include "pcg_basic.h"   //RNG
#include "utility.h"
#include "dataStructs.h"

/* energy */
double energy_coul_full(lattice config, double** ewald, para parameters);
double energy_ising_full(lattice config, su s, double* E, para parameters);
double energy_ising_su(lattice config, su s, para parameters);
double energy_ising_su_local(su s, para parameters, lattice* subset, int lsubset);
double energy_full(lattice config, su s, double* E, double** ewald, para parameters);
double energy_local(su s, double* E, para parameters, lattice* subset, int lsubset);
//void test_energy_coul(lattice config, fftw_complex* fftout, double* full, fftw_plan* fftplan, int opt, pcg32_random_t* rng, para parameters);
void test_energy_coul_simpCubic(lattice config, double** ewald, pcg32_random_t* rng, para parameters, FILE* myerr);
void test_energy_coul_rand(lattice config, double** ewald, pcg32_random_t* rng, para parameters, FILE* myerr);
double delta_energy_lr_singleflip(lattice c, double** ewald, int k);
double delta_energy_lr_doubleflip(lattice c, double** ewald, int k, int l);
double delta_energy_lr_nflip(lattice c, double** ewald, uint32_t* k, int nk);
double energy_lr(lattice c, double** ewald);
double energy_lr_singleflip(lattice c, double** ewald, int k);
double energy_lr_doubleflip(lattice c, double** ewald, int k, int l);
double energy_lr_nflip(lattice c, double** ewald, uint32_t* k, int nk);
double energy_umbr(umb umbrellaSprings);
