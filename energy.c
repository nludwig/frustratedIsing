#include "energy.h"

////////////////////////////////////////
//ENERGY FUNCTIONS WITH ISING ENERGIES//
////////////////////////////////////////

////////////////////////////////////////////
//USE RAW Q+, Q- TO COMPUTE ISING ENERGIES//
////////////////////////////////////////////
#if !ISI_USE_SIGN
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

double energy_full(lattice c, su s, double* E, double** ewald, para p) {
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

double energy_full(lattice c, su s, double* E, double** ewald, para p) {
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
double energy_coul_full(lattice c, double** ewald, para p) {
        extern double pi;
        #if SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma)) - p->e_cut;
        #endif
        #if !SHIFT_COUL
        const double rInvErfcNN=erfc(1.0/(p->sigma));
        #endif
        extern int N,d;
        const int NNN=2*d;
        double e_c=0.0;
        for(int i=0;i<N;++i) {
                double tmp=0.0;
                for(int j=0;j<NNN;++j) {
                        const lattice neigh=c[i].nearNeigh[j];
                        if(c[i].su!=-1 && neigh->su!=-1) continue;
                        tmp+=*neigh->val * rInvErfcNN;
                }
                for(int j=0;j<c[i].nCoulNeigh;++j) {
                        const lattice neigh=c[i].coulNeigh[j];
                        if(c[i].su!=-1 && neigh->su!=-1) continue;
                        const double rInvErfc=c[i].dInvErfcNeigh[j];
                        tmp+=*neigh->val * rInvErfc;
                }
                e_c+=*c[i].val * tmp;
        }
        double e_cfft=0.0;
        #if FFT_ON
        e_cfft=energy_lr(c,ewald);
        #endif
        return (p->Q)*(0.5*e_c+e_cfft);
}

void test_energy_coul_simpCubic(lattice c, double** ewald, pcg32_random_t* rng, para p, FILE* myerr) {
        fprintf(myerr,"test_energy_coul_simpCubic...\n");
        const double Msc=1.748;
        double* chg;
        const int opt=0;
        if(opt==0)      chg=setLatVal_simpCubic(c);
        //else if(opt==1) chg=setLatVal_spacedSimpCubic(c);
        //else if(opt==2) chg=setLatVal_spacedShiftedSimpCubic(c);
        //else if(opt==3) chg=setLatVal_simpCubicWDefect(c,rng);
                
        assert(chg[0]-chg[1]==0 /*neutral lattice*/);
        free(chg); chg=NULL;
        const double EcRpFFT=energy_coul_full(c,ewald,p);
        const double EcFFT=(p->Q)*energy_lr(c,ewald);
        const double EcR=EcRpFFT-EcFFT;
        double EcTot;
        if(opt==0)      EcTot=EcR+EcFFT-(p->e_self);
        //else if(opt==1) EcTot=EcR+EcFFT-(p->e_self)/8.0;
        //else if(opt==2) EcTot=EcR+EcFFT-(p->e_self)/8.0;
        //else if(opt==3) EcTot=EcR+EcFFT-(p->e_self);

        extern int N;
        const double r0=1.0;
        double EcLat=-N*(p->Q)/(r0)*Msc/2.0; // /2.0 because otherwise double-counting (think of 1/2 \sum q_i*V_i)
        if(opt==1 || opt==2) EcLat/=8.0;

        fprintf(myerr,"Ec_expect\t%f\n"
                      "Ec_predic\t%f\n"
                      "Ec_ratio\t%f\n",
                      EcLat,EcTot,EcTot/EcLat);
        fprintf(myerr,"EcR\t%f\n"
                      "EcFFT\t%f\n"
                      "EcSelf\t%f\n",
                      EcR,EcFFT,(opt==0 || opt==3) ? (p->e_self) : (p->e_self)/8.0);
        fprintf(myerr,"M_expect\t%f\n"
                      "M_predic\t%f\n",
                      Msc,EcTot/(-N*(p->Q)/(r0)));
                      
        fprintf(myerr,"test_energy_coul: ending... ");
        double relErr=(EcTot-EcLat)/EcLat;
        relErr*=(relErr<0.0 ? -1.0 : 1.0);
        if(relErr<=0.01) fprintf(myerr,"success: relErr=%f <= 0.01\n",relErr);
        else {
                fprintf(myerr,"FAIL: relErr=%f <= 0.01. USE CAUTION WHEN PROCEEDING:"
                              " MAY WANT TO RECONSIDER PARAMETER SETTINGS (sigma in particular)\n",relErr);
        }
        fflush(myerr);
}

void test_energy_coul_rand(lattice c, double** ewald, pcg32_random_t* rng, para p, FILE* myerr) {
        fprintf(myerr,"test_energy_coul_rand...\n");
        const double Msc=1.748;
        double* chg;
        const int opt=0;
        if(opt==0)      chg=setLatVal_rand(c,rng,myerr);
        //else if(opt==1) chg=setLatVal_spacedSimpCubic(c);
        //else if(opt==2) chg=setLatVal_spacedShiftedSimpCubic(c);
        //else if(opt==3) chg=setLatVal_simpCubicWDefect(c,rng);
                
        assert(chg[0]-chg[1]==0 /*neutral lattice*/);
        free(chg); chg=NULL;
        const double EcRpFFT=energy_coul_full(c,ewald,p);
        const double EcFFT=(p->Q)*energy_lr(c,ewald);
        const double EcR=EcRpFFT-EcFFT;
        double EcTot;
        if(opt==0)      EcTot=EcR+EcFFT-(p->e_self);
        //else if(opt==1) EcTot=EcR+EcFFT-(p->e_self)/8.0;
        //else if(opt==2) EcTot=EcR+EcFFT-(p->e_self)/8.0;
        //else if(opt==3) EcTot=EcR+EcFFT-(p->e_self);

        extern int N;
        const double r0=1.0;
        double EcLat=-N*(p->Q)/(r0)*Msc/2.0; // /2.0 because otherwise double-counting (think of 1/2 \sum q_i*V_i)
        if(opt==1 || opt==2) EcLat/=8.0;

        fprintf(myerr,"Ec_expect\t%f\n"
                      "Ec_predic\t%f\n"
                      "Ec_ratio\t%f\n",
                      EcLat,EcTot,EcTot/EcLat);
        fprintf(myerr,"EcR\t%f\n"
                      "EcFFT\t%f\n"
                      "EcSelf\t%f\n",
                      EcR,EcFFT,(opt==0 || opt==3) ? (p->e_self) : (p->e_self)/8.0);
        fprintf(myerr,"M_expect\t%f\n"
                      "M_predic\t%f\n",
                      Msc,EcTot/(-N*(p->Q)/(r0)));
        fprintf(myerr,"test_energy_coul: ending... ");
        double relErr=(EcTot-EcLat)/EcLat;
        relErr*=(relErr<0.0 ? -1.0 : 1.0);
        if(relErr<=0.01) fprintf(myerr,"success: relErr=%f <= 0.01\n",relErr);
        else {
                fprintf(myerr,"FAIL: relErr=%f <= 0.01. USE CAUTION WHEN PROCEEDING:"
                              "MAY WANT TO RECONSIDER PARAMETER SETTINGS (sigma in particular)\n",relErr);
        }
        fflush(myerr);
}

//assumes spin has NOT already been flipped!
//(ie. new spin = -1*ck)
//to change this convention, just multiply
//through deltaCk by -1
double delta_energy_lr_singleflip(lattice c, double** ewald, int k) {
        extern int N;
        double deltaElr=0.0;
        for(int j=0;j<k;++j) {
                const int kmj=k-j-1;
                if(c[j].hydrophobic==true) continue;
                deltaElr+=*c[j].val*ewald[j][kmj];
        }
        for(int j=k+1;j<N;++j) {
                const int jmk=j-k-1;
                if(c[j].hydrophobic==true) continue;
                deltaElr+=*c[j].val*ewald[k][jmk];
        }
        //assumes spin has NOT already been flipped!
        const double deltaCk=(-(*c[k].val)) - (*c[k].val);
        deltaElr*=deltaCk;
        extern double pi;
        return 4.0*pi*deltaElr/N;
}

//assumes spins have NOT already been flipped!
//(ie. new spins = -1*ck , -1*cl)
//to change this convention, just multiply
//through deltaCk,deltaCl by -1
double delta_energy_lr_doubleflip(lattice c, double** ewald, int k, int l) {
        assert(k<l);
        extern int N;
        double deltaElrk=0.0, deltaElrl=0.0;
        for(int j=0;j<k;++j) {
                if(c[j].hydrophobic==true) continue;
                const int kmj=k-j-1;
                deltaElrk+=*c[j].val*ewald[j][kmj];
                const int lmj=l-j-1;
                deltaElrl+=*c[j].val*ewald[j][lmj];
        }
        //middle portion has k flipped already, but l same
        for(int j=k+1;j<l;++j) {
                if(c[j].hydrophobic==true) continue;
                const int jmk=j-k-1;
                deltaElrk+=*c[j].val*ewald[k][jmk];
                const int lmj=l-j-1;
                deltaElrl+=*c[j].val*ewald[j][lmj];
        }
        //last portion has both k and l flipped
        for(int j=l+1;j<N;++j) {
                if(c[j].hydrophobic==true) continue;
                const int jmk=j-k-1;
                deltaElrk+=*c[j].val*ewald[k][jmk];
                const int jml=j-l-1;
                deltaElrl+=*c[j].val*ewald[l][jml];
        }
        //assumes spins have NOT already been flipped!
        const double deltaCk=(-(*c[k].val)) - (*c[k].val);
        const double deltaCl=(-(*c[l].val)) - (*c[l].val);
        deltaElrk*=deltaCk;
        deltaElrl*=deltaCl;
        extern double pi;
        return 4.0*pi*(deltaElrk+deltaElrl)/N;
}

//generalizes delta_energy_lr_doubleflip
//assumes spins have NOT already been flipped!
//(ie. new spins = -1*ck , -1*cl)
//to change this convention, just multiply
//through deltaCk,deltaCl by -1
double delta_energy_lr_nflip(lattice c, double** ewald, uint32_t* k, int nk) {
        for(int i=1;i<nk;++i) assert(k[i-1]<k[i]);
        extern int N;
        double deltaElr[nk];
        for(int i=0;i<nk;++i) deltaElr[i]=0.0;
        for(uint32_t j=0;j<k[0];++j) {
                if(c[j].hydrophobic==true) continue;
                const double cjv=*c[j].val;
                for(int l=0;l<nk;++l) {
                        const uint32_t kl=k[l];
                        const uint32_t klmj=kl-j-1u;
                        deltaElr[l]+=cjv*ewald[j][klmj];
                }
        }
        for(int i=1;i<nk;++i) {
                for(uint32_t j=k[i-1]+1;j<k[i];++j) {
                        if(c[j].hydrophobic==true) continue;
                        const double cjv=*c[j].val;
                        for(int l=0;l<i;++l) {
                                const uint32_t kl=k[l];
                                const uint32_t jmkl=j-kl-1u;
                                deltaElr[l]+=cjv*ewald[kl][jmkl];
                        }
                        for(int l=i;l<nk;++l) {
                                const uint32_t kl=k[l];
                                const uint32_t klmj=kl-j-1u;
                                deltaElr[l]+=cjv*ewald[j][klmj];
                        }
                }
        }
        for(uint32_t j=k[nk-1]+1;j<N;++j) {
                if(c[j].hydrophobic==true) continue;
                const double cjv=*c[j].val;
                for(int l=0;l<nk;++l) {
                        const uint32_t kl=k[l];
                        const uint32_t jmkl=j-kl-1u;
                        deltaElr[l]+=cjv*ewald[kl][jmkl];
                }
        }
        
        double deltaElrSum=0.0;
        for(int i=0;i<nk;++i) { //assumes spins have NOT already been flipped!
                const uint32_t ki=k[i];
                const double deltaCki=(-(*c[ki].val)) - (*c[ki].val);
                deltaElrSum+=deltaCki*deltaElr[i];
        }
        extern double pi;
        return 4.0*pi*deltaElrSum/N;
}

double energy_lr(lattice c, double** ewald) {
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
                        const int lmj=l-j-1;
                        Ej+=*c[l].val*ewald[j][lmj];
                }
                Elr+=*c[j].val*Ej;
        }
        Elr+=ewald[N-1][0]*N/2.0;
        extern double pi;
        return 4.0*pi*Elr/N;
}

double energy_lr_singleflip(lattice c, double** ewald, int k) {
        #if !FFT_ON
        return 0.0;
        #endif
        if(c[k].hydrophobic==true) return 0.0;
        extern int N;
        double Elr=0.0;
        for(int j=0;j<k;++j) {
                if(c[j].hydrophobic==true) continue;
                const int kmj=k-j-1;
                Elr+=*c[j].val*ewald[j][kmj];
        }
        for(int j=k+1;j<N;++j) {
                if(c[j].hydrophobic==true) continue;
                const int jmk=j-k-1;
                Elr+=*c[j].val*ewald[k][jmk];
        }
        Elr*=*c[k].val;
        extern double pi;
        return 4.0*pi*Elr/N;
}

double energy_lr_doubleflip(lattice c, double** ewald, int k, int l) {
        #if !FFT_ON
        return 0.0;
        #endif
        assert(k<l);
        extern int N;
        double Elrk=0.0, Elrl=0.0;
        for(int j=0;j<k;++j) {
                if(c[j].hydrophobic==true) continue;
                const int kmj=k-j-1;
                Elrk+=*c[j].val*ewald[j][kmj];
                const int lmj=l-j-1;
                Elrl+=*c[j].val*ewald[j][lmj];
        }
        //middle portion has k flipped already, but l same
        for(int j=k+1;j<l;++j) {
                if(c[j].hydrophobic==true) continue;
                const int jmk=j-k-1;
                Elrk+=*c[j].val*ewald[k][jmk];
                const int lmj=l-j-1;
                Elrl+=*c[j].val*ewald[j][lmj];
        }
        //last portion has both k and l flipped
        for(int j=l+1;j<N;++j) {
                if(c[j].hydrophobic==true) continue;
                const int jmk=j-k-1;
                Elrk+=*c[j].val*ewald[k][jmk];
                const int jml=j-l-1;
                Elrl+=*c[j].val*ewald[l][jml];
        }
        Elrk+=*c[k].val*ewald[N-1][0];
        Elrl+=*c[l].val*ewald[N-1][0];
        if(c[k].hydrophobic==false) Elrk*=*c[k].val;
        if(c[l].hydrophobic==false) Elrl*=*c[l].val;
        double Elrkl=0.0;
        if(c[k].hydrophobic==false && c[k].hydrophobic==false) {
                Elrkl=(*c[k].val)*ewald[k][k-l-1]*(*c[l].val);
        }
        extern double pi;
        return 4.0*pi*(Elrk+Elrl+Elrkl)/N;
        //return 4.0*pi*(Elrk+Elrl)/N;
}

double energy_lr_nflip(lattice c, double** ewald, uint32_t* k, int nk) {
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
                        const uint32_t klmj=kl-j-1u;
                        Elr[l]+=cjv*ewald[j][klmj];
                }
        }
        for(int i=1;i<nk;++i) {
                for(uint32_t j=k[i-1]+1;j<k[i];++j) {
                        if(c[j].hydrophobic==true) continue;
                        const double cjv=*c[j].val;
                        for(int l=0;l<i;++l) {
                                const uint32_t kl=k[l];
                                const uint32_t jmkl=j-kl-1u;
                                Elr[l]+=cjv*ewald[kl][jmkl];
                        }
                        for(int l=i;l<nk;++l) {
                                const uint32_t kl=k[l];
                                const uint32_t klmj=kl-j-1u;
                                Elr[l]+=cjv*ewald[j][klmj];
                        }
                }
        }
        for(uint32_t j=k[nk-1]+1;j<N;++j) {
                if(c[j].hydrophobic==true) continue;
                const double cjv=*c[j].val;
                for(int l=0;l<nk;++l) {
                        const uint32_t kl=k[l];
                        const uint32_t jmkl=j-kl-1u;
                        Elr[l]+=cjv*ewald[kl][jmkl];
                }
        }

        double ElrSum=0.0;
        for(int i=0;i<nk;++i) {
                const uint32_t ki=k[i];
                if(c[ki].hydrophobic==true) continue;
                ElrSum+=*c[ki].val*(Elr[i]+*c[ki].val*ewald[N-1][0]);
                //ElrSum+=*c[ki].val*Elr[i];
                if(i<nk-1) {
                        for(int j=i+1;j<nk;++j) {
                                const uint32_t kj=k[j];
                                if(c[kj].hydrophobic==true) continue;
                                const uint32_t kjmki=kj-ki-1u;
                                const double Elrkikj=(*c[ki].val)*ewald[ki][kjmki]*(*c[kj].val);
                                ElrSum+=Elrkikj;
                        }
                }
        }
        extern double pi;
        return 4.0*pi*ElrSum/N;
}

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
