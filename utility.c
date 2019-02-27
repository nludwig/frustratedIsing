#include "utility.h"

/*******************************/
/********    utility    ********/
/*******************************/

// returns *nwords words which were separated in
// str by elements from sep
//
// nwords: either one-element int* or NULL:
//      NULL -> don't return nwords
char** split(char* str, const char* sep, int* nwords) {
        int len=strlen(str);
        char* mem=malloc(len+1);
        assert(mem!=NULL /*malloc err*/);
        memset(mem,'\0',len+1);
        char* temp[len];
        for(int i=0;i<len;++i)  temp[i]=NULL;
        char* ptr=str;
        int i=0;
        while(1) {
                const int n=strcspn(ptr,sep);
                temp[i]=mem+(ptr-str);
                strncpy(temp[i],ptr,n);
                ++i;
                ptr+=n+1;
                if(*(ptr-1)=='\0') break;
        }
        char** words=malloc(i*sizeof(*words));
        assert(words!=NULL /*malloc err*/);
        for(int j=0;j<i;++j) words[j]=temp[j];
        if(nwords!=NULL)     *nwords=i;
        return words;
}

void swapD(double* a, double* b) {
        double t=*a;
        *a=*b;
        *b=t;
}

// supporting function for qsortMult
// adapted from K&R page 120
void swapMultD(double** v, int narr, int L, int a, int b) {
        assert(a<L && b<L);
        double tmp;
        for(int i=0;i<narr;++i) {
                tmp=v[i][a];
                v[i][a]=v[i][b];
                v[i][b]=tmp;
                return;
        }
}

// qsort multiple arrays of same length based on the first array
// adapted from K&R page 120
void qsortMult(double** v, int narr, int L, int l, int r, int (*comp)(double*, double*)) {
        int i, last;
        if(l>=r) return;
        swapMultD(v,narr,L,l,(l+r)/2);
        last=l;
        for(i=l+1;i<=r;++i) {
                if((*comp)(v[i],v[l]) < 0) swapMultD(v,narr,L,++last,i);
        }
        swapMultD(v,narr,L,l,last);
        qsortMult(v,narr,L,l,last-1,comp);
        qsortMult(v,narr,L,last+1,r,comp);
}

void freeAll(void** ptrs, int l) {
        for(int i=0;i<l;++i) {
                free(ptrs[i]);
                ptrs[i]=NULL;
        }
}

int iabs(int a) {
        if(a>=0) return a;
        else     return -a;
}

uint16_t ndigits(long int num) {
        if(num==0l) return 1u;
        if(num<0l)  num*=-1l;
        const uint32_t n=(uint32_t)num;
        uint32_t pow=1u,ndigits=1u;
        while(n>9u*pow) {
                pow*=10u;
                ++ndigits;
        }
        return ndigits;
}

int iCmp(int* a, int* b) {
        return ( *(int*)a > *(int*)b ) - ( *(int*)a < *(int*)b );
}

int lattCmp(lattice* a, lattice* b) {
        return ( *(lattice*)a > *(lattice*)b ) - ( *(lattice*)a < *(lattice*)b );
}

int lattCmpRev(lattice* a, lattice* b) {
        return  ( *(lattice*)a < *(lattice*)b ) - ( *(lattice*)a > *(lattice*)b );
}

int uint32_tCmp(uint32_t* a, uint32_t* b) {
        return ( *(uint32_t*)a > *(uint32_t*)b ) - ( *(uint32_t*)a < *(uint32_t*)b );
}

int doubleMatrixCmp(double** a, double** b, int d) {
        for(int i=0;i<d;++i) {
                for(int j=0;j<d;++j) {
                        if(a[i][j]!=b[i][j]) return 0;
                }
        }
        return 1;
}

//set all parameters to 0 by default
void initPara(para p) {
        p->J=0.0;
        p->Q=0.0;
        p->kT=0.0;
        p->beta=0.0;
        p->beta_eff=0.0;
        p->sigma=0.0;
        p->coulCutoff=0.0;
        p->e_cut=0.0;
        p->e_self=0.0;
        p->fftOn=0;
}

int intArrMin(int* a, int l) {
        int m=a[0];
        for(int i=1;i<l;++i) if(a[i]<m) m=a[i];
        return m;
}

int intArrMax(int* a, int l) {
        int m=a[0];
        for(int i=1;i<l;++i) if(a[i]>m) m=a[i];
        return m;
}

/*******************************/
/****** lattice related   ******/
/*******************************/

void getPos(int n, int* pos, int* L_in) {
        if(L_in==NULL) {
                extern int* L;
                L_in=L;
        }
        extern int d;
        int prod=1;
        n+=1; //req.d to work
        for(int e=0;e<d;++e) prod*=L_in[e];
        for(int e=0;e<d;++e) {
                prod/=L_in[e];
                int i=0;
                while(i*prod < n) ++i;
                --i;    //go one too far
                pos[e]=i;
                n-=i*prod;
        }
}

void getPosLat(lattice c, lattice cn, int* pos) {
        const int n=(int)(cn-c);
        getPos(n,pos,NULL);
}

int getSite(int* pos, int* L_in) {
        extern int d;
        if(L_in==NULL) {
                extern int* L;
                L_in=L;
        }
        int n=0,prod=1;
        for(int e=d-1;e>=0;--e) { //ROW MAJOR
                n+=pos[e]*prod;
                prod*=L_in[e];
        }
        //for(int e=0;e<d;++e) { COL MAJOR
        //        int prod=1;
        //        for(int j=0;j<e;++j) prod*=L[j];
        //        n+=pos[e]*prod;
        //}
        return n;
}

int getSitePBC(int* pos, int* L_in) {
        wrapIntoL(pos,L_in);
        return getSite(pos,L_in);
}

lattice getSiteLat(lattice c, int* pos) {
        return c+getSitePBC(pos,NULL);
}

int getNN(int ind, int e, int* L_in) {
        extern int d;
        int nnpos[d];
        getPos(ind,nnpos,L_in);
        int sign=1;
        if(e<0) {
                e=-1*e-1;       //[-d,-1] are negative directions of [0,d): -1:0; -2:1; ... -d:(d-1)
                sign=-1;
        }
        nnpos[e]+=sign;
        int nn=getSite(nnpos,L_in);
        return nn;
}

int getNNPBC(int ind, int e, int* L_in) {
        extern int d;
        int nnpos[d];
        getPos(ind,nnpos,L_in);
        int sign=1;
        if(e<0) {
                e=-1*e-1;       //[-d,-1] are negative directions of [0,d): -1:0; -2:1; ... -d:(d-1)
                sign=-1;
        }
        nnpos[e]+=sign;
        int nn=getSitePBC(nnpos,L_in);
        return nn;
}

lattice getNNLat(lattice c, int ind, int e) {
        return c+getNNPBC(ind,e,NULL);
}

void setCharge(lattice c) {
        extern int N,nPlusSvSites,nMinusSvSites,nZeroSvSites;
        extern int nPlusSuSites,nMinusSuSites,nHydrSuSites;
        extern double totSvPlus,totSvMinus;
        extern double totSuPlus,totSuMinus;
        nPlusSvSites=nMinusSvSites=nZeroSvSites=0;
        nPlusSuSites=nMinusSuSites=nHydrSuSites=0;
        totSvPlus=totSvMinus=0.0;
        totSuPlus=totSuMinus=0.0;
        for(int i=0;i<N;++i) {
                if(c[i].su==-1) {
                        if(*c[i].val>0.0) {
                                ++nPlusSvSites;
                                totSvPlus+=*c[i].val;
                        }
                        else if(*c[i].val<0.0) {
                                ++nMinusSvSites;
                                totSvMinus+=*c[i].val;
                        }
                        else ++nZeroSvSites;
                }
                else {
                        if(c[i].hydrophobic==false) {
                                if(*c[i].val>0.0) {
                                        ++nPlusSuSites;
                                        totSuPlus+=*c[i].val;
                                }
                                else if(*c[i].val<0.0) {
                                        ++nMinusSuSites;
                                        totSuMinus+=*c[i].val;
                                }
                        }
                        else ++nHydrSuSites;
                }
        }
}

double* checkCharge(lattice c, int lc, double* ret) {
        extern int N;
        if(lc==-1) lc=N;
        int currPlusSvSites=0,currMinusSvSites=0,currZeroSvSites=0,currHydrSvSites=0;
        int currPlusSuSites=0,currMinusSuSites=0,currZeroSuSites=0,currHydrSuSites=0;
        double currSvPlus=0.0,currSvMinus=0.0;
        double currSuPlus=0.0,currSuMinus=0.0;
        for(int i=0;i<lc;++i) {
                if(c[i].su==-1) {
                        if(c[i].hydrophobic==false) {
                                if(*c[i].val>0.0) {
                                        ++currPlusSvSites;
                                        currSvPlus+=*c[i].val;
                                }
                                else if(*c[i].val<0.0) {
                                        ++currMinusSvSites;
                                        currSvMinus+=*c[i].val;
                                }
                                else ++currZeroSvSites;
                        }
                        else    ++currHydrSvSites;
                }
                else {
                        if(c[i].hydrophobic==false) {
                                if(*c[i].val>0.0) {
                                        ++currPlusSuSites;
                                        currSuPlus+=*c[i].val;
                                }
                                else if(*c[i].val<0.0) {
                                        ++currMinusSuSites;
                                        currSuMinus+=*c[i].val;
                                }
                                else ++currZeroSuSites;
                        }
                        else ++currHydrSuSites;
                }
        }
        if(ret==NULL) {
                ret=malloc(12*sizeof(*ret));
                assert(ret!=NULL /*malloc*/);
        }
        ret[0]=currPlusSvSites;
        ret[1]=currMinusSvSites;
        ret[2]=currZeroSvSites;
        ret[3]=currHydrSvSites;
        ret[4]=currSvPlus;
        ret[5]=currSvMinus;
        ret[6]=currPlusSuSites;
        ret[7]=currMinusSuSites;
        ret[8]=currZeroSuSites;
        ret[9]=currHydrSuSites;
        ret[10]=currSuPlus;
        ret[11]=currSuMinus;
        return ret;
}

para copyPara(para p) {
        para pNew=malloc(1*sizeof(*pNew));
        assert(pNew!=NULL /*malloc*/);
//as of 181010:
//struct parameters {
//        double J;
//        double Q;
//        double kT;
//        double beta;
//        double beta_eff;
//        double sigma;
//        double coulCutoff;
//        double e_cut;
//        double e_self; //const: shows up in output only
//        int fftOn;
//};
        pNew->J=p->J;
        pNew->Q=p->Q;
        pNew->kT=p->kT;
        pNew->beta=p->beta;
        pNew->beta_eff=p->beta_eff;
        pNew->sigma=p->sigma;
        pNew->coulCutoff=p->coulCutoff;
        pNew->e_cut=p->e_cut;
        pNew->e_self=p->e_self;
        pNew->fftOn=p->fftOn;
        return pNew;
}

/*******************************/
/********  bc related   ********/
/*******************************/

//double pbc(double dst, int e)
//{
//        extern int* L;
//        dst -= L[e] * round(dst / L[e]);
//        return dst;
//}

double pbc(double dst, int e) {
        extern int* L;
        while(dst<-L[e]/2) dst+=L[e];
        while(dst>=L[e]/2) dst-=L[e];
        return dst;
}

void wrapIntoL(int* pos, int* L_in) {
        if(L_in==NULL) {
                extern int* L;
                L_in=L;
        }
        extern int d;
        extern int* bc;
        for(int e=0;e<d;++e) {
                if(pos[e]>=0 && pos[e]<L_in[e]) continue;
                assert(bc[2*e]==1 || bc[2*e+1]==1 /*wrapIntoL: bc.s must be set*/);
                while(pos[e]<0)        pos[e]+=L_in[e];
                while(pos[e]>=L_in[e]) pos[e]-=L_in[e];
                assert(pos[e]>=0 /*wrapIntoL: pos still out of (below) box*/);
                assert(pos[e]<L_in[e] /*wrapIntoL: pos still out of (above) box*/);
        }
}

void wrapIntoL_double(double* pos, int* L_in) {
        if(L_in==NULL) {
                extern int* L;
                L_in=L;
        }
        extern int d;
        extern int* bc;
        for(int e=0;e<d;++e) {
                if(pos[e]>=0.0 && pos[e]<L_in[e])    continue;
                assert(bc[2*e]==1 || bc[2*e+1]==1 /*bc.s not set*/);
                while(pos[e]<0.0 && bc[2*e]==1)           pos[e]+=L_in[e];
                while(pos[e]>=L_in[e] && bc[2*e+1]==1)  pos[e]-=L_in[e];
                assert(pos[e]>=0.0 /*pos still out of (below) box*/);
                assert(pos[e]<L_in[e] /*pos still out of (above) box*/);
        }
}

double dist_intvect_nopbc(int* u, int* v) {
        int r=0.0;
        extern int d;
        for(int e=0;e<d;++e) {
                const int dif=u[e]-v[e];
                r+=dif*dif;
        }
        return sqrt((double)r);
}

double dist_vect(double* u, double* v) {
        double r=0.0;
        extern int d;
        extern int* bc;
        for(int e=0;e<d;++e) {
                double dif=u[e]-v[e];
                if(bc[2*e]==1) dif=pbc(round(dif), e);
                r+=dif*dif;
        }
        return sqrt(r);
}

//dr must have d entries
void diff_lat(double* dr, int i, int j) {
        extern int d;
        int u[d];
        int v[d];
        getPos(i,u,NULL);
        getPos(j,v,NULL);
        for(int e=0;e<d;++e) {
                dr[e]=(double)(u[e]-v[e]);
                dr[e]=pbc(dr[e],e);
        }
}

double dist_lat(int i, int j) {
        extern int d;
        int u[d];
        int v[d];
        getPos(i,u,NULL);
        getPos(j,v,NULL);
        double ud[d];
        double vd[d];
        for(int e=0;e<d;++e) {
                ud[e]=(double)u[e];
                vd[e]=(double)v[e];
        }
        const double r=dist_vect(ud,vd);
        return r;
}