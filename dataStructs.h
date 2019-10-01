#pragma once

// data structures
struct parameters {
        double J; //Ising strength
        double Q; //electrostatic strength
        double kT; //thermal energy
        double beta; //inverse thermal energy
        double betaEff_clst; //effective beta for formation of clusters
        double sigma; //Ewald short- long-length scale cutoff/Gaussian size
        double coulCutoff; //short-range (real space) Coulomb cutoff distance
        double e_cut; //const: energy of Coulomb interaction @ cutoff
        double e_self; //const: shows up in output only
        int fftOn; //0 -> no long-range Ewald sum
};
typedef struct parameters* para;
struct latticeSite {
        int nCoulNeigh; //number of neighbors for short-range part of Coulomb
        int su; //-1: no solute; !(-1): solute #
        bool hydrophobic; //is site hydrophobic?
        struct latticeSite** nearNeigh; //set of 2*d nearest neighbors of site
        struct latticeSite** coulNeigh; //set of nCoulNeigh Coulomb neighbors of site
        double* val; //valence of site
        double* dInvErfcNeigh; //1/r*erfc(r/sigma) of each Coulomb neighbor of site
};
typedef struct latticeSite* lattice;
struct solute {
        double suStep; //range of su step size (uniform over [0,1) * suStep)
        double surfType; //solute site valence
        double isiMult; //solute Ising valence multiplier
        double mvProb; //probability to attempt su move given su is activated
        int nsites; //number site in su
        int ninshl; //number of sites in optional su shell
        int hlw; //is su hollow? (not implemented)
        int nMotion; //in how many dimensions is su allowed to move? (0->stationary, 1->z only, ...)
        int* linDim; //size of su (ie. if rect, length in x,y,z)
        double* com; //current su COM
        double* unwrappedCom; //current su COM (not wrapped into box; useful for diffusion)
        double** orientation; //rotation matrix
        double** relPos; //position of each site relative to su COM position
        double** inShlRelPos; //position of each shell site """"
        struct latticeSite** currPos; //lattice site occupied by each site based on current COM
        struct latticeSite** inShlCurrPos; //lattice shell site occupied by each site based on current COM
        char* shape; //cube, rect, sphere, ...
        bool* hydrophobic; //is site [i] of su hydrophobic?
};
typedef struct solute* su;
struct umbrellas {
        int npairs; //number of umbrellas
        su** pair; //pair i: su j,k
        double* k; //spring constant of pair i
        double* r0; //equilibrium spring length of pair i
        char** energyIOUnits; //kT?
};
typedef struct umbrellas* umb;
