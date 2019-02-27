#pragma once

/* data structures */
struct parameters {
        double J;
        double Q;
        double kT;
        double beta;
        double beta_eff;
        double sigma;
        double coulCutoff;
        double e_cut;
        double e_self; //const: shows up in output only
        int fftOn;
};
typedef struct parameters* para;
struct latticeSite {
        int nCoulNeigh;
        int su;        //-1: no solute; !(-1): solute #
        bool hydrophobic;
        struct latticeSite** nearNeigh;
        struct latticeSite** coulNeigh;
        double* val;
        double* dInvErfcNeigh;
};
typedef struct latticeSite* lattice;
struct solute {
        double suStep;
        double surfType;
        double isiMult;
        double mvProb;
        int nsites;
        int ninshl;
        int hlw;
        int nMotion;
        int* linDim;
        double* com;
        double* unwrappedCom;
        double** orientation;
        double** relPos;
        double** inShlRelPos;
        struct latticeSite** currPos;
        struct latticeSite** inShlCurrPos;
        char* shape;
        bool* hydrophobic;
};
typedef struct solute* su;
struct umbrellas {
        int npairs;
        su** pair;
        double* k;
        double* r0;
        char** energyIOUnits;
};
typedef struct umbrellas* umb;
