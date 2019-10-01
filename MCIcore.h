#pragma once

#include "headerBundle.h"

/* global var.s 
 * double pi=acos(-1);
 * int N;
 * int Nsu;
 * int Ntotsusites;
 * int Nc;
 * int d;
 * int shlL;
 * int nPlusSites
 * int nMinusSites
 * int nZeroSites
 * double totSvPlus
 * double totSvMinus
 * int types;
 * int inshlopt;
 * int* L;
 * int* bc;      0: no PBC; 1: PBC
 * FILE* myout;
 * FILE* myerr;
 * FILE* mvstats=NULL;
 * FILE* mySuComOut=NULL;
 * FILE* myumbout=NULL;
 */

/* solute (su) related */
int setSuRelPos(lattice c, su solutes, int innerShellOption, int neutralOverall, int coreNum, pcg32_random_t* rng, FILE* myerr);
void updSuCurPos(lattice c, su solutes, int ind, double* com);
void updLatSuStatus(lattice config, su solutes);
void updLatHydrophobicStatus(lattice config, su solutes);
void updSuCurPos_latSuStatus(lattice config, su solutes, int ind); //wraps updSuCurPos&updLatSuStatus
bool setHydrophobicExist(su solutes);
int checkSuOverlap(su solutes, int ind, int innerShellOption);
//int findVacated(lattice config, lattice trialConfig, su solutes, su trialSolutes, int ind, lattice** v, int lv);

/* core */
//int setLatPos(lattice config);
int maintainChargeNeut(lattice config, uint32_t* ii, int lii, int chgSwitch, int* dPlus, int* dMinus, int* dZero, int coreNum, pcg32_random_t* rng, FILE* myerr); //lc len c or -1 for N, dZero 1-rho or -1 if use nZeroSites
void neutralizeStep(lattice c, int rand, int* del, double* charges, int ntypes);
void setLatInitSu(lattice config); //init lat to no solutes
int willParaLatBeChargeNeut(int coreNum); //0 no; 1 yes
double* setLatVal_rand(lattice config, su s, int coreNum, pcg32_random_t* rng, FILE* myerr);
double* setLatVal_data(lattice config, char** lines, int nlines, int coreNum, FILE* myerr);
void** setLatNeigh_coul_twoway(lattice config, para parameters, FILE* myerr);
double* setupEwald(int* L, double sigma, int nCoresP);
int localSubset(lattice* core, int lcore, lattice** subset);
lattice* numToLat(lattice config, uint32_t* nums, uint32_t numnums);
uint32_t* latToNum(lattice config, lattice* sites, int nsites);
int spinSwapMoveNN(lattice config, lattice trialconfig, su s, double* en_oldLR, int site, double* ewald, para parameters, pcg32_random_t* rng, FILE* myerr);
int spinSwapMove(lattice config, lattice trialconfig, su s, double* en_oldLR, int site, double* ewald, para parameters, pcg32_random_t* rng, FILE* myerr);
int clusterMove(lattice config, lattice trialconfig, su s, double* E, const uint32_t site, double* ewald, para parameters, int coreNum, pcg32_random_t* rng, FILE* myerr);
uint32_t* buildCluster(lattice config, const uint32_t site0, const uint32_t site1, const uint32_t rotnAxis, uint32_t* nInCluster, double* suBlocks, para parameters, int coreNum, pcg32_random_t* rng);
double* comOrientationNudge(lattice c, su s, int ind, double* trialCom, double** trialOrientation, pcg32_random_t* rng, FILE* myerr);
int suMove(lattice config, lattice trialconfig, su solutes, su trialsolutes, umb umbrella, umb trialumbrella, double* E, int solInd, double* ewald, int innerShellOption, bool rotationOpt, double sweepsMult, double kTeff, para parameters, int coreNum, pcg32_random_t* rng, FILE* myerr);
uint32_t buildSuMoveCore(lattice c, lattice tc, su s, su ts, int ind, uint32_t** core); //union
uint32_t addNNToCore(lattice c, su s, uint32_t** core, uint32_t lcore);
uint32_t* addNNToCore_toNewArr(lattice c, su s, uint32_t* core, uint32_t* lcore);
int findVacated(lattice config, lattice trialConfig, su solutes, su trialSolutes, int ind, uint32_t** v); //intersection
int relaxSvNearSu(lattice tconfig, su ts, uint32_t* core, uint32_t lcore, double* ewald, double sweepsMult, double kTeff, para p, int coreNum, pcg32_random_t* rng, FILE* myerr);
#if EQUI_INNER_DUMP
int* MCstep_flipswap(lattice config, lattice trialconfig, double* Ein, su solutes, su trialsolutes, umb umbrella, umb trialumbrella, double* ewald, int innerShellOption, bool rotationOpt, double sweepsMult, double kTeff, para parameters, int coreNum, pcg32_random_t* rng, int enFreq, int dumpFreq, int suComFreq, int umbrFreq, int dataFreq, FILE* dumpout, char* prefix, int outerN, bool outputStyleFI, FILE* myout, FILE* myumbout, FILE* mySuComOut, FILE* myerr);
#endif
#if !EQUI_INNER_DUMP
int* MCstep_flipswap(lattice config, lattice trialconfig, double* Ein, su solutes, su trialsolutes, umb umbrella, umb trialumbrella, double* ewald, int innerShellOption, bool rotationOpt, double sweepsMult, double kTeff, para parameters, int coreNum, pcg32_random_t* rng, FILE* myout, FILE* myumbout, FILE* mySuComOut, FILE* myerr);
#endif
int* MCstep_cluster(lattice config, lattice trialconfig, double* Ein, su solutes, su trialsolutes, umb umbrella, umb trialumbrella, double* ewald, int innerShellOption, bool rotationOpt, double sweepsMult, double kTeff, para parameters, int coreNum, pcg32_random_t* rng, FILE* mvstats, FILE* myerr);

/* spoof */
void runSpoof(lattice config, su s, double* ewald, para parameters, umb umbrella, FILE* f, FILE* myout);
int prepForReadTrj(int* nlines, int* nchars, FILE* f); //returns niter
int readTrjStep(char*** lines, char* linesHolder, FILE* f); //returns nchars
int findSpace(char* s, int spaceNum);
void lammpstrjStepToLat(lattice c, char** lines);
void fitrjStepToLat(lattice c, char** lines);
