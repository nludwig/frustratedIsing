#include "energy.h"

//functions for computing the energy given a lattice
//configuration and set of parameters. minor variations
//in functions to allow for different interaction energies
//controlled by pre-processor switches found in switches.h.
//
//for standard, FI usage, the main functions that will be
//called are
//-energy_full() which computes the energy
//of the entire lattice;
//-energy_lr() which computes the long-range/Fourier space
//part of the electrostatic energy of the entire lattice;
//-energy_local() which computes the short-range/real space
//part of the energy of a local region of the lattice
//(quicker than computing entire lattice energy when only
//local changes have been made, such as a swap move, cluster
//move or solute move);
//-energy_lr_nflip() which computes the long-range/Fourier
//space part of the electrostatic energy for the n input
//charges, and is also quicker than recomputing for the
//entire lattice

////////////////////////////////////////
//ENERGY FUNCTIONS WITH ISING ENERGIES//
////////////////////////////////////////

////////////////////////////////////////////
//USE RAW Q+, Q- TO COMPUTE ISING ENERGIES//
////////////////////////////////////////////
#if !ISI_USE_SIGN
//compute Ising "fluid-fluid" energy of entire lattice
double energy_ising_full(lattice c, su s, double* E, para p) {
        extern int N,d;
        const int NNN=2*d;
        double e_i=0.0,e_isu=0.0;
        for(int i=0;i<N;++i) {
                double tmpi=0.0,tmpisu=0.0;
                if(c[i].su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) tmpi+=*neigh->val;
                                else              tmpisu+=*neigh->val * s[neigh->su].isiMult;
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) tmpisu+=*neigh->val * s[c[i].su].isiMult;
                                else              continue;
                        }
                }
                e_i+=*c[i].val * tmpi;
                e_isu+=*c[i].val * tmpisu;
        }
        e_i*=0.5; //correct for double-counting
        e_isu*=0.5;
        if(E!=NULL) {
                E[0]=-(p->J)*e_i;
                E[1]=-(p->J)*e_isu;
                return E[0]+E[1];
        }
        return -(p->J)*(e_i+e_isu);
}

//compute Ising "fluid-solute" energy of entire lattice
double energy_ising_su(lattice c, su s, para p) {
        extern int N,d;
        const int NNN=2*d;
        double e_isu=0.0;
        for(int i=0;i<N;++i) {
                double tmpisu=0.0;
                if(c[i].su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) continue;
                                else              tmpisu+=*neigh->val * s[neigh->su].isiMult;
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) tmpisu+=*neigh->val * s[c[i].su].isiMult;
                                else              continue;
                        }
                }
                e_isu+=*c[i].val * tmpisu;
        }
        return -(p->J)*e_isu*0.5;
}

//compute Ising "fluid-solute" energy of local region
double energy_ising_su_local(su s, para p, lattice* subset, int lsubset) {
        extern int d;
        const int NNN=2*d;
        double e_isu=0.0;
        for(int i=0;i<lsubset;++i) {
                double tmpisu=0.0;
                if(subset[i]->su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) continue;
                                else              tmpisu+=*neigh->val * s[neigh->su].isiMult;
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) tmpisu+=*neigh->val * s[subset[i]->su].isiMult;
                                else              continue;
                        }
                }
                e_isu+=*subset[i]->val * tmpisu;
        }
        return -(p->J)*e_isu*0.5;
}

//compute Ising and electrostatic energy of entire lattice
double energy_full(lattice c, su s, double* E, double* ewald, para p) {
        extern double pi;
        #if SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma)) - p->e_cut;
        #endif
        #if !SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma));
        #endif
        extern int N,d;
        const int NNN=2*d;
        double e_isvsv=0.0,e_isvsu=0.0,e_isusu=0.0,e_csrsvsv=0.0,e_csrsvsu=0.0,e_csrsusu=0.0;
        double e_ihh=0.0,e_isvh=0.0,e_isuh=0.0;
        for(int i=0;i<N;++i) {
                double tmpisvsv=0.0,tmpisvsu=0.0,tmpisusu=0.0,tmpcsrsvsv=0.0,tmpcsrsvsu=0.0,tmpcsrsusu=0.0;
                double tmpihh=0.0,tmpisvh=0.0,tmpisuh=0.0;
                if(c[i].su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) {
                                        tmpisvsv+=*neigh->val;
                                        #if COUL_ON
                                        tmpcsrsvsv+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                        #endif
                                }
                                else {
                                        if(neigh->hydrophobic==false) {
                                                tmpisvsu+=*neigh->val * s[neigh->su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                tmpisvh+=-fabs(*neigh->val * *c[i].val);
                                        }
                                }
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(c[i].hydrophobic==false) {
                                                tmpisvsu+=*neigh->val * s[c[i].su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                tmpisvh+=-fabs(*neigh->val * *c[i].val);
                                        }
                                }
                                else {
                                        if(c[i].hydrophobic==false) {
                                                if(neigh->hydrophobic==false) {
                                                        tmpisusu+=*neigh->val * s[neigh->su].isiMult*s[c[i].su].isiMult;
                                                        #if COUL_ON
                                                        tmpcsrsusu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                        #endif
                                                }
                                                else { //one false
                                                        tmpisuh+=-fabs(*neigh->val * *c[i].val);
                                                }
                                        }
                                        else {
                                                if(neigh->hydrophobic==false) { //one false
                                                        tmpisuh+=-fabs(*neigh->val * *c[i].val);
                                                }
                                                else { //both hydrophobic
                                                        tmpihh+=fabs(*neigh->val * *c[i].val);
                                                }
                                        }
                                }
                        }
                }
                e_isvsv+=*c[i].val * tmpisvsv;
                e_isvsu+=*c[i].val * tmpisvsu;
                e_isusu+=*c[i].val * tmpisusu;
                e_ihh+=tmpihh;
                e_isvh+=tmpisvh;
                e_isuh+=tmpisuh;
                #if COUL_ON     //Frustrated Ising
                if(c[i].hydrophobic==false) {
                        if(c[i].su==-1) {
                                for(int j=0;j<c[i].nCoulNeigh;++j) {
                                        if(c[i].coulNeigh[j]->hydrophobic==true) continue;
                                        const lattice neigh=c[i].coulNeigh[j];
                                        const double rInvErfc=c[i].dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsv+=*neigh->val * rInvErfc;
                                        else              tmpcsrsvsu+=*neigh->val * rInvErfc;
                                }
                        }
                        else {
                                for(int j=0;j<c[i].nCoulNeigh;++j) {
                                        if(c[i].coulNeigh[j]->hydrophobic==true) continue;
                                        const lattice neigh=c[i].coulNeigh[j];
                                        const double rInvErfc=c[i].dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsu+=*neigh->val * rInvErfc;
                                        else              tmpcsrsusu+=*neigh->val * rInvErfc;
                                }
                        }
                }
                e_csrsvsv+=*c[i].val * tmpcsrsvsv;
                e_csrsvsu+=*c[i].val * tmpcsrsvsu;
                e_csrsusu+=*c[i].val * tmpcsrsusu;
                #endif
        }
        double e_clr=0.0;
        #if FFT_ON
        e_clr=(p->Q)*energy_lr(c,ewald);
        #endif
        e_isvsv*=-(p->J)*0.5; //correct for double-counting
        e_isvsu*=-(p->J)*0.5;
        e_isusu*=-(p->J)*0.5;
        e_csrsvsv*=(p->Q)*0.5;
        e_csrsvsu*=(p->Q)*0.5;
        e_csrsusu*=(p->Q)*0.5;
        e_ihh*=-(p->J)*0.5;
        e_isvh*=-(p->J)*0.5;
        e_isuh*=-(p->J)*0.5;
        if(E!=NULL) {
                E[0]=e_isvsv;
                E[1]=e_isvsu;
                E[2]=e_isusu;
                E[3]=e_csrsvsv;
                E[4]=e_csrsvsu;
                E[5]=e_csrsusu;
                E[6]=e_clr;
                E[7]=e_ihh;
                E[8]=e_isvh;
                E[9]=e_isuh;
                return E[0]+E[1]+E[2]+E[3]+E[4]+E[5]+E[6]+E[7]+E[8]+E[9];
        }
        return e_isvsv+e_isvsu+e_isusu+e_csrsvsv+e_csrsvsu+e_csrsusu+e_clr+e_ihh+e_isvh+e_isuh;
}

//compute Ising and electrostatic energy of local region
double energy_local(su s, double* E, para p, lattice* subset, int lsubset) {
        extern double pi;
        #if SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma)) - p->e_cut;
        #endif
        #if !SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma));
        #endif
        extern int d;
        const int NNN=2*d;
        double e_isvsv=0.0,e_isvsu=0.0,e_isusu=0.0,e_csrsvsv=0.0,e_csrsvsu=0.0,e_csrsusu=0.0;
        double e_ihh=0.0,e_isvh=0.0,e_isuh=0.0;
        for(int i=0;i<lsubset;++i) {
                double tmpisvsv=0.0,tmpisvsu=0.0,tmpisusu=0.0,tmpcsrsvsv=0.0,tmpcsrsvsu=0.0,tmpcsrsusu=0.0;
                double tmpihh=0.0,tmpisvh=0.0,tmpisuh=0.0;
                if(subset[i]->su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) {
                                        tmpisvsv+=*neigh->val;
                                        #if COUL_ON
                                        tmpcsrsvsv+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                        #endif
                                }
                                else {
                                        if(neigh->hydrophobic==false) {
                                                tmpisvsu+=*neigh->val * s[neigh->su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                tmpisvh+=-fabs(*neigh->val * *subset[i]->val);
                                        }
                                }
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(subset[i]->hydrophobic==false) {
                                                tmpisvsu+=*neigh->val * s[subset[i]->su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                tmpisvh+=-fabs(*neigh->val * *subset[i]->val);
                                        }
                                }
                                else {
                                        if(subset[i]->hydrophobic==false) {
                                                if(neigh->hydrophobic==false) {
                                                        tmpisusu+=*neigh->val * s[neigh->su].isiMult * s[subset[i]->su].isiMult;
                                                        #if COUL_ON
                                                        tmpcsrsusu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                        #endif
                                                }
                                                else { //one false
                                                        tmpisuh+=-fabs(*neigh->val * *subset[i]->val);
                                                }
                                        }
                                        else {
                                                if(neigh->hydrophobic==false) { //one false
                                                        tmpisuh+=-fabs(*neigh->val * *subset[i]->val);
                                                }
                                                else { //both hydrophobic
                                                        tmpihh+=fabs(*neigh->val * *subset[i]->val);
                                                }
                                        }
                                }
                        }
                }
                e_isvsv+=*subset[i]->val * tmpisvsv;
                e_isvsu+=*subset[i]->val * tmpisvsu;
                e_isusu+=*subset[i]->val * tmpisusu;
                e_ihh+=tmpihh;
                e_isvh+=tmpisvh;
                e_isuh+=tmpisuh;
                #if COUL_ON     //Frustrated Ising
                if(subset[i]->hydrophobic==false) {
                        if(subset[i]->su==-1) {
                                for(int j=0;j<subset[i]->nCoulNeigh;++j) {
                                        if(subset[i]->coulNeigh[j]->hydrophobic==true) continue;
                                        const lattice neigh=subset[i]->coulNeigh[j];
                                        const double rInvErfc=subset[i]->dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsv+=*neigh->val * rInvErfc;
                                        else              tmpcsrsvsu+=*neigh->val * rInvErfc;
                                }
                        }
                        else {
                                for(int j=0;j<subset[i]->nCoulNeigh;++j) {
                                        if(subset[i]->coulNeigh[j]->hydrophobic==true) continue;
                                        const lattice neigh=subset[i]->coulNeigh[j];
                                        const double rInvErfc=subset[i]->dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsu+=*neigh->val * rInvErfc;
                                        else              tmpcsrsusu+=*neigh->val * rInvErfc;
                                }
                        }
                }
                e_csrsvsv+=*subset[i]->val * tmpcsrsvsv;
                e_csrsvsu+=*subset[i]->val * tmpcsrsvsu;
                e_csrsusu+=*subset[i]->val * tmpcsrsusu;
                #endif
        }
        e_isvsv*=-(p->J)*0.5; //correct for double-counting
        e_isvsu*=-(p->J)*0.5;
        e_isusu*=-(p->J)*0.5;
        e_csrsvsv*=(p->Q)*0.5;
        e_csrsvsu*=(p->Q)*0.5;
        e_csrsusu*=(p->Q)*0.5;
        e_ihh*=-(p->J)*0.5;
        e_isvh*=-(p->J)*0.5;
        e_isuh*=-(p->J)*0.5;
        if(E!=NULL) {
                E[0]=e_isvsv;
                E[1]=e_isvsu;
                E[2]=e_isusu;
                E[3]=e_csrsvsv;
                E[4]=e_csrsvsu;
                E[5]=e_csrsusu;
                E[7]=e_ihh;
                E[8]=e_isvh;
                E[9]=e_isuh;
                return E[0]+E[1]+E[2]+E[3]+E[4]+E[5]+E[7]+E[8]+E[9];
        }
        return e_isvsv+e_isvsu+e_isusu+e_csrsvsv+e_csrsvsu+e_csrsusu+e_ihh+e_isvh+e_isuh;
}
#endif

//////////////////////////////////////////////
//USE SIGN OF SPIN TO COMPUTE ISING ENERGIES//
//(generalized to replace q+,q- with macros //
// defined in header:                       //
// q+ -> ISI_SPIN_UP                        //
// q- -> ISI_SPIN_DN                        //
//////////////////////////////////////////////
#if ISI_USE_SIGN
//compute Ising "fluid-fluid" energy of entire lattice
double energy_ising_full(lattice c, su s, double* E, para p) {
        extern int N,d;
        const int NNN=2*d;
        double e_i=0.0,e_isu=0.0;
        for(int i=0;i<N;++i) {
                double tmpi=0.0,tmpisu=0.0;
                if(c[i].su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(*neigh->val>0.0)      tmpi+=ISI_SPIN_UP;
                                        else if(*neigh->val<0.0) tmpi+=ISI_SPIN_DN;
                                }
                                else {
                                        if(*neigh->val>0.0)      tmpisu+=ISI_SPIN_UP*s[neigh->su].isiMult;
                                        else if(*neigh->val<0.0) tmpisu+=ISI_SPIN_DN*s[neigh->su].isiMult;
                                }
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(*neigh->val>0.0)      tmpisu+=ISI_SPIN_UP*s[c[i].su].isiMult;
                                        else if(*neigh->val<0.0) tmpisu+=ISI_SPIN_DN*s[c[i].su].isiMult;
                                }
                                else continue; //both solute sites
                        }
                }
                if(*c[i].val>0.0) {
                        e_i+=ISI_SPIN_UP*tmpi;
                        e_isu+=ISI_SPIN_UP*tmpisu;
                }
                else if(*c[i].val<0.0) {
                        e_i+=ISI_SPIN_DN*tmpi;
                        e_isu+=ISI_SPIN_DN*tmpisu;
                }
        }
        e_i*=0.5; //correct for double-counting
        e_isu*=0.5;
        if(E!=NULL) {
                E[0]=-(p->J)*e_i;
                E[1]=-(p->J)*e_isu;
                return E[0]+E[1];
        }
        return -(p->J)*(e_i+e_isu);
}

//compute Ising "fluid-solute" energy of entire lattice
double energy_ising_su(lattice c, su s, para p) {
        extern int N,d;
        const int NNN=2*d;
        double e_isu=0.0;
        for(int i=0;i<N;++i) {
                double tmpisu=0.0;
                if(c[i].su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) continue;
                                else {
                                        if(*neigh->val>0.0)      tmpisu+=ISI_SPIN_UP*s[neigh->su].isiMult;
                                        else if(*neigh->val<0.0) tmpisu+=ISI_SPIN_DN*s[neigh->su].isiMult;
                                }
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(*neigh->val>0.0)      tmpisu+=ISI_SPIN_UP*s[c[i].su].isiMult;
                                        else if(*neigh->val<0.0) tmpisu+=ISI_SPIN_DN*s[c[i].su].isiMult;
                                }
                                else continue; //both solute sites
                        }
                }
                if(*c[i].val>0.0) {
                        e_isu+=ISI_SPIN_UP*tmpisu;
                }
                else if(*c[i].val<0.0) {
                        e_isu+=ISI_SPIN_DN*tmpisu;
                }
        }
        return -(p->J)*e_isu*0.5;
}

//compute Ising "fluid-solute" energy of local region
double energy_ising_su_local(su s, para p, lattice* subset, int lsubset) {
        extern int d;
        const int NNN=2*d;
        double e_isu=0.0;
        for(int i=0;i<lsubset;++i) {
                double tmpisu=0.0;
                if(subset[i]->su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) continue;
                                else {
                                        if(*neigh->val>0.0)      tmpisu+=ISI_SPIN_UP*s[neigh->su].isiMult;
                                        else if(*neigh->val<0.0) tmpisu+=ISI_SPIN_DN*s[neigh->su].isiMult;
                                }
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(*neigh->val>0.0)      tmpisu+=ISI_SPIN_UP*s[subset[i]->su].isiMult;
                                        else if(*neigh->val<0.0) tmpisu+=ISI_SPIN_DN*s[subset[i]->su].isiMult;
                                }
                                else continue; //both solute sites
                        }
                }
                if(*subset[i]->val>0.0) {
                        e_isu+=ISI_SPIN_UP*tmpisu;
                }
                else if(*subset[i]->val<0.0) {
                        e_isu+=ISI_SPIN_DN*tmpisu;
                }
        }
        return -(p->J)*e_isu*0.5;
}

//compute Ising and electrostatic energy of entire lattice
double energy_full(lattice c, su s, double* E, double* ewald, para p) {
        extern double pi;
        #if SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma)) - p->e_cut;
        #endif
        #if !SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma));
        #endif
        extern int N,d;
        const int NNN=2*d;
        double e_isvsv=0.0,e_isvsu=0.0,e_isusu=0.0,e_csrsvsv=0.0,e_csrsvsu=0.0,e_csrsusu=0.0;
        double e_ihh=0.0,e_isvh=0.0,e_isuh=0.0;
        for(int i=0;i<N;++i) {
                double tmpisvsv=0.0,tmpisvsu=0.0,tmpisusu=0.0,tmpcsrsvsv=0.0,tmpcsrsvsu=0.0,tmpcsrsusu=0.0;
                double tmpihh=0.0,tmpisvh=0.0,tmpisuh=0.0;
                if(c[i].su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(*neigh->val>0.0)      tmpisvsv+=ISI_SPIN_UP;
                                        else if(*neigh->val<0.0) tmpisvsv+=ISI_SPIN_DN;
                                        #if COUL_ON
                                        tmpcsrsvsv+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                        #endif
                                }
                                else {
                                        if(neigh->hydrophobic==false) {
                                                if(*neigh->val>0.0)      tmpisvsu+=ISI_SPIN_UP*s[c[i].su].isiMult;
                                                else if(*neigh->val<0.0) tmpisvsu+=ISI_SPIN_DN*s[c[i].su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                double tmpenergy=0.0;
                                                if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                else                    tmpenergy+=ISI_SPIN_DN;
                                                if(*c[i].val>0.0)       tmpenergy*=ISI_SPIN_UP;
                                                else                    tmpenergy*=ISI_SPIN_DN;
                                                tmpisvh+=-fabs(tmpenergy);
                                        }
                                }
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=c[i].nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(c[i].hydrophobic==false) {
                                                if(*neigh->val>0.0)      tmpisvsu+=ISI_SPIN_UP*s[neigh->su].isiMult;
                                                else if(*neigh->val<0.0) tmpisvsu+=ISI_SPIN_DN*s[neigh->su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                double tmpenergy=0.0;
                                                if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                else                    tmpenergy+=ISI_SPIN_DN;
                                                if(*c[i].val>0.0)       tmpenergy*=ISI_SPIN_UP;
                                                else                    tmpenergy*=ISI_SPIN_DN;
                                                tmpisvh+=-fabs(tmpenergy);
                                        }
                                }
                                else {
                                        if(c[i].hydrophobic==false) {
                                                if(neigh->hydrophobic==false) {
                                                        if(*neigh->val>0.0)      tmpisusu+=ISI_SPIN_UP*s[neigh->su].isiMult*s[c[i].su].isiMult;
                                                        else if(*neigh->val<0.0) tmpisusu+=ISI_SPIN_DN*s[neigh->su].isiMult*s[c[i].su].isiMult;
                                                        #if COUL_ON
                                                        tmpcsrsusu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                        #endif
                                                }
                                                else { //one false
                                                        double tmpenergy=0.0;
                                                        if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                        else                    tmpenergy+=ISI_SPIN_DN;
                                                        if(*c[i].val>0.0)       tmpenergy*=ISI_SPIN_UP;
                                                        else                    tmpenergy*=ISI_SPIN_DN;
                                                        tmpisuh+=-fabs(tmpenergy);
                                                }
                                        }
                                        else {
                                                double tmpenergy=0.0;
                                                if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                else                    tmpenergy+=ISI_SPIN_DN;
                                                if(*c[i].val>0.0)       tmpenergy*=ISI_SPIN_UP;
                                                else                    tmpenergy*=ISI_SPIN_DN;
                                                if(neigh->hydrophobic==false) { //one false
                                                        tmpisuh+=-fabs(tmpenergy);
                                                }
                                                else { //both hydrophobic
                                                        tmpihh+=fabs(tmpenergy);
                                                }
                                        }
                                }
                        }
                }
                if(*c[i].val>0.0) {
                        e_isvsv+=ISI_SPIN_UP * tmpisvsv;
                        e_isvsu+=ISI_SPIN_UP * tmpisvsu;
                        e_isusu+=ISI_SPIN_UP * tmpisusu;
                }
                else if(*c[i].val<0.0) {
                        e_isvsv+=ISI_SPIN_DN * tmpisvsv;
                        e_isvsu+=ISI_SPIN_DN * tmpisvsu;
                        e_isusu+=ISI_SPIN_DN * tmpisusu;
                }
                e_ihh+=tmpihh;
                e_isvh+=tmpisvh;
                e_isuh+=tmpisuh;
                #if COUL_ON     //Frustrated Ising
                if(c[i].hydrophobic==false) {
                        if(c[i].su==-1) {
                                for(int j=0;j<c[i].nCoulNeigh;++j) {
                                        const lattice neigh=c[i].coulNeigh[j];
                                        const double rInvErfc=c[i].dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsv+=*neigh->val * rInvErfc;
                                        else              tmpcsrsvsu+=*neigh->val * rInvErfc;
                                }
                        }
                        else {
                                for(int j=0;j<c[i].nCoulNeigh;++j) {
                                        const lattice neigh=c[i].coulNeigh[j];
                                        const double rInvErfc=c[i].dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsu+=*neigh->val * rInvErfc;
                                        else              tmpcsrsusu+=*neigh->val * rInvErfc;
                                }
                        }
                }
                e_csrsvsv+=*c[i].val * tmpcsrsvsv;
                e_csrsvsu+=*c[i].val * tmpcsrsvsu;
                e_csrsusu+=*c[i].val * tmpcsrsusu;
                #endif
        }
        double e_clr=0.0;
        #if FFT_ON
        e_clr=(p->Q)*energy_lr(c,ewald);
        #endif
        e_isvsv*=-(p->J)*0.5; //correct for double-counting
        e_isvsu*=-(p->J)*0.5;
        e_isusu*=-(p->J)*0.5;
        e_csrsvsv*=(p->Q)*0.5;
        e_csrsvsu*=(p->Q)*0.5;
        e_csrsusu*=(p->Q)*0.5;
        e_ihh*=-(p->J)*0.5;
        e_isvh*=-(p->J)*0.5;
        e_isuh*=-(p->J)*0.5;
        if(E!=NULL) {
                E[0]=e_isvsv;
                E[1]=e_isvsu;
                E[2]=e_isusu;
                E[3]=e_csrsvsv;
                E[4]=e_csrsvsu;
                E[5]=e_csrsusu;
                E[6]=e_clr;
                E[7]=e_ihh;
                E[8]=e_isvh;
                E[9]=e_isuh;
                return E[0]+E[1]+E[2]+E[3]+E[4]+E[5]+E[6]+E[7]+E[8]+E[9];
        }
        return e_isvsv+e_isvsu+e_isusu+e_csrsvsv+e_csrsvsu+e_csrsusu+e_clr+e_ihh+e_isvh+e_isuh;
}

//compute Ising and electrostatic energy of local region
double energy_local(su s, double* E, para p, lattice* subset, int lsubset) {
        extern double pi;
        #if SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma)) - p->e_cut;
        #endif
        #if !SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma));
        #endif
        extern int d;
        const int NNN=2*d;
        double e_isvsv=0.0,e_isvsu=0.0,e_isusu=0.0,e_csrsvsv=0.0,e_csrsvsu=0.0,e_csrsusu=0.0;
        double e_ihh=0.0,e_isvh=0.0,e_isuh=0.0;
        for(int i=0;i<lsubset;++i) {
                double tmpisvsv=0.0,tmpisvsu=0.0,tmpisusu=0.0,tmpcsrsvsv=0.0,tmpcsrsvsu=0.0,tmpcsrsusu=0.0;
                double tmpihh=0.0,tmpisvh=0.0,tmpisuh=0.0;
                if(subset[i]->su==-1) {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(*neigh->val>0.0)      tmpisvsv+=ISI_SPIN_UP;
                                        else if(*neigh->val<0.0) tmpisvsv+=ISI_SPIN_DN;
                                        #if COUL_ON
                                        tmpcsrsvsv+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                        #endif
                                }
                                else {
                                        if(neigh->hydrophobic==false) {
                                                if(*neigh->val>0.0)      tmpisvsu+=ISI_SPIN_UP*s[neigh->su].isiMult;
                                                else if(*neigh->val<0.0) tmpisvsu+=ISI_SPIN_DN*s[neigh->su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                double tmpenergy=0.0;
                                                if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                else                    tmpenergy+=ISI_SPIN_DN;
                                                if(*subset[i]->val>0.0) tmpenergy*=ISI_SPIN_UP;
                                                else                    tmpenergy*=ISI_SPIN_DN;
                                                tmpisvh+=-fabs(tmpenergy);
                                        }
                                }
                        }
                }
                else {
                        for(int j=0;j<NNN;++j) {
                                const lattice neigh=subset[i]->nearNeigh[j];
                                if(neigh->su==-1) {
                                        if(subset[i]->hydrophobic==false) {
                                                if(*neigh->val>0.0)      tmpisvsu+=ISI_SPIN_UP*s[subset[i]->su].isiMult;
                                                else if(*neigh->val<0.0) tmpisvsu+=ISI_SPIN_DN*s[subset[i]->su].isiMult;
                                                #if COUL_ON
                                                tmpcsrsvsu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                #endif
                                        }
                                        else {
                                                double tmpenergy=0.0;
                                                if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                else                    tmpenergy+=ISI_SPIN_DN;
                                                if(*subset[i]->val>0.0) tmpenergy*=ISI_SPIN_UP;
                                                else                    tmpenergy*=ISI_SPIN_DN;
                                                tmpisvh+=-fabs(tmpenergy);
                                        }
                                }
                                else {
                                        if(subset[i]->hydrophobic==false) {
                                                if(neigh->hydrophobic==false) {
                                                        if(*neigh->val>0.0)      tmpisusu+=ISI_SPIN_UP*s[neigh->su].isiMult*s[subset[i]->su].isiMult;
                                                        else if(*neigh->val<0.0) tmpisusu+=ISI_SPIN_DN*s[neigh->su].isiMult*s[subset[i]->su].isiMult;
                                                        #if COUL_ON
                                                        tmpcsrsusu+=*neigh->val * rInvErfcNN; //NNs not stored in coulNeighs
                                                        #endif
                                                }
                                                else { //one false
                                                        double tmpenergy=0.0;
                                                        if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                        else                    tmpenergy+=ISI_SPIN_DN;
                                                        if(*subset[i]->val>0.0) tmpenergy*=ISI_SPIN_UP;
                                                        else                    tmpenergy*=ISI_SPIN_DN;
                                                        tmpisuh+=-fabs(tmpenergy);
                                                }
                                        }
                                        else {
                                                double tmpenergy=0.0;
                                                if(*neigh->val>0.0)     tmpenergy+=ISI_SPIN_UP;
                                                else                    tmpenergy+=ISI_SPIN_DN;
                                                if(*subset[i]->val>0.0) tmpenergy*=ISI_SPIN_UP;
                                                else                    tmpenergy*=ISI_SPIN_DN;
                                                if(neigh->hydrophobic==false) { //one false
                                                        tmpisuh+=-fabs(tmpenergy);
                                                }
                                                else {
                                                        tmpihh+=fabs(tmpenergy);
                                                }
                                        }
                                }
                        }
                }
                if(*subset[i]->val>0.0) {
                        e_isvsv+=ISI_SPIN_UP * tmpisvsv;
                        e_isvsu+=ISI_SPIN_UP * tmpisvsu;
                        e_isusu+=ISI_SPIN_UP * tmpisusu;
                }
                else if (*subset[i]->val<0.0) {
                        e_isvsv+=ISI_SPIN_DN * tmpisvsv;
                        e_isvsu+=ISI_SPIN_DN * tmpisvsu;
                        e_isusu+=ISI_SPIN_DN * tmpisusu;
                }
                e_ihh+=tmpihh;
                e_isvh+=tmpisvh;
                e_isuh+=tmpisuh;
                #if COUL_ON     //Frustrated Ising
                if(subset[i]->hydrophobic==false) {
                        if(subset[i]->su==-1) {
                                for(int j=0;j<subset[i]->nCoulNeigh;++j) {
                                        const lattice neigh=subset[i]->coulNeigh[j];
                                        const double rInvErfc=subset[i]->dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsv+=*neigh->val * rInvErfc;
                                        else              tmpcsrsvsu+=*neigh->val * rInvErfc;
                                }
                        }
                        else {
                                for(int j=0;j<subset[i]->nCoulNeigh;++j) {
                                        const lattice neigh=subset[i]->coulNeigh[j];
                                        const double rInvErfc=subset[i]->dInvErfcNeigh[j];
                                        if(neigh->su==-1) tmpcsrsvsu+=*neigh->val * rInvErfc;
                                        else              tmpcsrsusu+=*neigh->val * rInvErfc;
                                }
                        }
                }
                e_csrsvsv+=*subset[i]->val * tmpcsrsvsv;
                e_csrsvsu+=*subset[i]->val * tmpcsrsvsu;
                e_csrsusu+=*subset[i]->val * tmpcsrsusu;
                #endif
        }
        e_isvsv*=-(p->J)*0.5; //correct for double-counting
        e_isvsu*=-(p->J)*0.5;
        e_isusu*=-(p->J)*0.5;
        e_csrsvsv*=(p->Q)*0.5;
        e_csrsvsu*=(p->Q)*0.5;
        e_csrsusu*=(p->Q)*0.5;
        e_ihh*=-(p->J)*0.5;
        e_isvh*=-(p->J)*0.5;
        e_isuh*=-(p->J)*0.5;
        if(E!=NULL) {
                E[0]=e_isvsv;
                E[1]=e_isvsu;
                E[2]=e_isusu;
                E[3]=e_csrsvsv;
                E[4]=e_csrsvsu;
                E[5]=e_csrsusu;
                E[7]=e_ihh;
                E[8]=e_isvh;
                E[9]=e_isuh;
                return E[0]+E[1]+E[2]+E[3]+E[4]+E[5]+E[7]+E[8]+E[9];
        }
        return e_isvsv+e_isvsu+e_isusu+e_csrsvsv+e_csrsvsu+e_csrsusu+e_ihh+e_isvh+e_isuh;
}

#endif

//////////////////////////////
//NON-ISING ENERGY FUNCTIONS//
//////////////////////////////

//compute long-range/Fourier space part of electrostatic
//interaction over entire lattice
double energy_lr(lattice c, double* ewald) {
        #if !FFT_ON
        return 0.0;
        #endif
        extern int N;
        double Elr=0.0;
        for(int j=0;j<N-1;++j) {
                if(c[j].hydrophobic==true) continue;
                double Ej=0.0;
                for(int l=j+1;l<N;++l) {
                        if(c[l].hydrophobic==true) continue;
                        const int i=getEwaldIndexForRij(j,l);
                        Ej+=*c[l].val*ewald[i];
                }
                Elr+=*c[j].val*Ej;
        }
        Elr+=ewald[0]*N/2.0;
        extern double pi;
        return 4.0*pi*Elr/N;
}

//i is center
int getEwaldIndexForRij(int i, int j) {
        extern int d;
        int ri[d],rj[d],rij[d];
        getPos(i,ri);
        getPos(j,rj);
        for(int e=0;e<d;++e) {
                rij[e]=rj[e]-ri[e];
        }
        const int n=getWrappedSite(rij,NULL);
        return n;
}

//compute long-range/Fourier space part of electrostatic
//interaction over spins k
double energy_lr_nflip(lattice c, double* ewald, uint32_t* k, int nk) {
        #if !FFT_ON
        return 0.0;
        #endif
        for(int i=1;i<nk;++i) assert(k[i-1]<k[i]);
        extern int N;
        double Elr[nk];
        for(int i=0;i<nk;++i) Elr[i]=0.0;
        for(uint32_t j=0;j<k[0];++j) {
                if(c[j].hydrophobic==true) continue;
                const double cjv=*c[j].val;
                for(int l=0;l<nk;++l) {
                        const uint32_t kl=k[l];
                        const int i=getEwaldIndexForRij(kl,j);
                        Elr[l]+=cjv*ewald[i];
                }
        }
        for(int i=1;i<nk;++i) {
                for(uint32_t j=k[i-1]+1;j<k[i];++j) {
                        if(c[j].hydrophobic==true) continue;
                        const double cjv=*c[j].val;
                        for(int l=0;l<i;++l) {
                                const uint32_t kl=k[l];
                                const int i=getEwaldIndexForRij(kl,j);
                                Elr[l]+=cjv*ewald[i];
                        }
                        for(int l=i;l<nk;++l) {
                                const uint32_t kl=k[l];
                                const int i=getEwaldIndexForRij(kl,j);
                                Elr[l]+=cjv*ewald[i];
                        }
                }
        }
        for(uint32_t j=k[nk-1]+1;j<N;++j) {
                if(c[j].hydrophobic==true) continue;
                const double cjv=*c[j].val;
                for(int l=0;l<nk;++l) {
                        const uint32_t kl=k[l];
                        const int i=getEwaldIndexForRij(kl,j);
                        Elr[l]+=cjv*ewald[i];
                }
        }

        double ElrSum=0.0;
        for(int i=0;i<nk;++i) {
                const uint32_t ki=k[i];
                if(c[ki].hydrophobic==true) continue;
                ElrSum+=*c[ki].val*(Elr[i]+*c[ki].val*ewald[0]);
                if(i<nk-1) {
                        for(int j=i+1;j<nk;++j) {
                                const uint32_t kj=k[j];
                                if(c[kj].hydrophobic==true) continue;
                                const int l=getEwaldIndexForRij(ki,kj);
                                const double Elrkikj=(*c[ki].val)*ewald[l]*(*c[kj].val);
                                ElrSum+=Elrkikj;
                        }
                }
        }
        extern double pi;
        return 4.0*pi*ElrSum/N;
}

//compute energy of umbrellas u
double energy_umbr(umb u) {
        if(u==NULL) return 0.0;
        double e_u=0.0;
        for(int i=0;i<u->npairs;++i) {
                double* com0=(u->pair[i][0])->com;
                double* com1=(u->pair[i][1])->com;
                const double r0=u->r0[i];
                const double sep=dist_vect(com0,com1);
                const double dr=sep-r0;
                e_u+=dr*dr*u->k[i];
        }
        return e_u/2.0;
}
