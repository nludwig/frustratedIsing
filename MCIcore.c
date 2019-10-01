#include "MCIcore.h"

/*******************************/
/******* solute (su) related ***/
/*******************************/

int setSuRelPos(lattice c, su s, int inshlopt, int neutralOverall, int coreNum, pcg32_random_t* rng, FILE* myerr) {
        extern int Nsu,shlL,d;
        int comPH=0;
        int dct;
        //build relPos for each solute, update curPos
        fprintf(myerr,"setSuRelPos: starting routine\n");
        fflush(myerr);
        for(int i=0;i<Nsu;++i) {
                if(strcmp(s[i].shape,"cube")==0) {
                        int** cube=buildCube(s[i].linDim[0]);
                        int* Lcube=NULL;
                        int border[d];
                        int ninterface;
                        for(int e=0;e<d;++e) border[e]=shlL;
                        int** interfacePos=getShapeInterfacePos(cube,s[i].nsites,&ninterface,&Lcube,border);
                        if(s[i].ninshl!=ninterface) {
                                fprintf(myerr,"setSuRelPos: s[%d].ninshl:%d  !=  ninterface:%d. setting val to ninterface\n",i,s[i].ninshl,ninterface);
                                fflush(myerr);
                                s[i].ninshl=ninterface;
                        }
                        for(int j=0;j<s[i].nsites;++j) {
                                for(int e=0;e<d;++e) {
                                        s[i].relPos[j][e]=cube[j][e]-(s[i].linDim[0]-1)/2.0;
                                }
                        }
                        for(int j=0;j<s[i].ninshl;++j) {
                                for(int e=0;e<d;++e) {
                                        s[i].inShlRelPos[j][e]=interfacePos[j][e]-(s[i].linDim[0]-1)/2.0;
                                }
                        }
                        free(Lcube);
                        free(*cube); free(cube);
                        free(*interfacePos); free(interfacePos);
                }
                else if(strcmp(s[i].shape,"rectangle")==0 || strcmp(s[i].shape,"rect")==0) {
                        int** rect=buildRectangle(s[i].linDim);
                        int* Lcube=NULL;
                        int border[d];
                        int ninterface=0;
                        for(int e=0;e<d;++e) border[e]=shlL;
                        int** interfacePos=getShapeInterfacePos(rect,s[i].nsites,&ninterface,&Lcube,border);
                        if(s[i].ninshl!=ninterface) {
                                fprintf(myerr,"setSuRelPos: s[%d].ninshl:%d  !=  ninterface:%d. setting val to ninterface\n",i,s[i].ninshl,ninterface);
                                s[i].ninshl=ninterface;
                        }
                        for(int j=0;j<s[i].nsites;++j) {
                                for(int e=0;e<d;++e) {
                                        s[i].relPos[j][e]=rect[j][e]-(s[i].linDim[e]-1)/2.0-border[e];
                                }
                        }
                        for(int j=0;j<s[i].ninshl;++j) {
                                for(int e=0;e<d;++e) {
                                        s[i].inShlRelPos[j][e]=interfacePos[j][e]-(s[i].linDim[e]-1)/2.0-border[e];
                                }
                        }
                        #if TEST_BUILD
                        fprintf(myerr,"linDim: %d\t%d\t%d\t\tborder: %d\t%d\t%d\n",s[i].linDim[0],s[i].linDim[1],s[i].linDim[2],border[0],border[1],border[2]);
                        fprintf(myerr,"rect x\ty\tz\t\t\n");
                        for(int j=0;j<s[i].nsites;++j) {
                                fprintf(myerr,"%d\t%d\t%d\n",rect[j][0],rect[j][1],rect[j][2]);
                        }
                        fprintf(myerr,"interf x\ty\tz\t\t\n");
                        for(int j=0;j<s[i].ninshl;++j) {
                                fprintf(myerr,"%d\t%d\t%d\t\t%.1f\t%.1f\t%.1f\n",
                                interfacePos[j][0],interfacePos[j][1],interfacePos[j][2],
                                interfacePos[j][0]-(s[i].linDim[0]-1)/2.0-border[0],
                                interfacePos[j][1]-(s[i].linDim[1]-1)/2.0-border[1],
                                interfacePos[j][2]-(s[i].linDim[2]-1)/2.0-border[2]);
                        }
                        #endif
                        free(Lcube);
                        free(*rect); free(rect);
                        if(ninterface>0) {
                                free(*interfacePos);
                                free(interfacePos);
                        }
                }
                else if(strcmp(s[i].shape,"sphere")==0 || strcmp(s[i].shape,"sph")==0) {
                        //find number of cells in sphere
                        int nSph=0;
                        const int lCube=2*s[i].linDim[0];
                        int** cube=buildCube(lCube);
                        const int nCube=lCube*lCube*lCube;
                        const double R2=s[i].linDim[0]*s[i].linDim[0];
                        for(int j=0;j<nCube;++j) {
                                double dist2=0.0;
                                for(int e=0;e<d;++e) {
                                        const double re=cube[j][e]-(lCube-1)/2.0;
                                        dist2+=re*re;
                                }
                                if(dist2<=R2) ++nSph;
                        }
                        fprintf(myerr,"setSuRelPos nSph: %d\n",nSph);
                        //build relPos's of sphere
                        int* sphereHold=malloc(nSph*d*sizeof(*sphereHold));
                        int** sphere=malloc(nSph*sizeof(*sphere));
                        assert(sphereHold!=NULL && sphere!=NULL /*malloc*/);
                        for(int j=0;j<nSph;++j) sphere[j]=sphereHold+j*d;
                        int k=0;
                        for(int j=0;j<nCube;++j) {
                                double dist2=0.0;
                                for(int e=0;e<d;++e) {
                                        const double re=cube[j][e]-(lCube-1)/2.0;
                                        dist2+=re*re;
                                }
                                if(dist2<=R2) {
                                        for(int e=0;e<d;++e) {
                                                sphere[k][e]=cube[j][e];
                                                s[i].relPos[k][e]=cube[j][e]-(lCube-1)/2.0;
                                        }
                                        ++k;
                                }
                                if(k>=nSph) break;
                        }
                        //build interface of sphere
                        int* Lcube=NULL;
                        int border[d];
                        int ninterface;
                        for(int e=0;e<d;++e) border[e]=shlL;
                        int** interfacePos=getShapeInterfacePos(sphere,s[i].nsites,&ninterface,&Lcube,border);
                        if(s[i].ninshl!=ninterface) {
                                fprintf(myerr,"setSuRelPos: s[%d].ninshl:%d  !=  ninterface:%d. setting val to ninterface\n",i,s[i].ninshl,ninterface);
                                s[i].ninshl=ninterface;
                        }
                        for(int j=0;j<s[i].ninshl;++j) {
                                for(int e=0;e<d;++e) {
                                        s[i].inShlRelPos[j][e]=interfacePos[j][e]-(s[i].linDim[0]-1)/2.0;
                                }
                        }
                        free(Lcube);
                        free(*cube); free(cube);
                        free(*interfacePos); free(interfacePos);
                }
                else {
                        fprintf(myerr,"shape '%s' not yet implemented. try cube, rect, or sphere\n",s[i].shape);
                        exit(1);
                }
                updSuCurPos(c,s,i,(double*)NULL);
        }
        //update the lattice
        updLatSuStatus(c,s);
        updLatHydrophobicStatus(c,s);

        //assume all su sites of same sign have same valence
        for(int i=0;i<Nsu;++i) {
                for(int j=i+1;j<Nsu;++j) {
                        if((s[i].surfType>0.0 && s[j].surfType>0.0) ||
                           (s[i].surfType<0.0 && s[j].surfType<0.0)) {
                                assert(s[i].surfType==s[j].surfType);
                        }
                }
        }
        double suPlusCharge=0.0,suMinusCharge=0.0;
        for(int i=0;i<Nsu;++i) {
                if(s[i].surfType>0.0) {
                        suPlusCharge=s[i].surfType;
                }
                else if(s[i].surfType<0.0) {
                        suMinusCharge=s[i].surfType;
                }
                if(suPlusCharge!=0.0 && suMinusCharge!=0.0) break;
        }
        //update values of sites in su's
        extern double plusCharge,minusCharge;
        int delPlus=0,delMinus=0,delZero=0,delSuPlus=0,delSuMinus=0;
        for(int i=0;i<Nsu;++i) {
                if(checkSuOverlap(s,i,inshlopt)==1) {
                        fprintf(myerr,"setSuRelPos: solutes with given COM, size overlap. Exiting.\n");
                        return 1;
                }
                for(int j=0;j<s[i].nsites;++j) {
                        const double origVal=*s[i].currPos[j]->val;
                        *s[i].currPos[j]->val=s[i].surfType;
                        if(s[i].surfType==suPlusCharge) delSuPlus+=1; //for su
                        else                            delSuMinus+=1; //for su
                        if(origVal==plusCharge)         delPlus-=1;
                        else if(origVal==minusCharge)   delMinus-=1;
                        else if(origVal==0.0)           delZero-=1;
                        else if(origVal==suPlusCharge)  delSuPlus-=1;
                        else if(origVal==suMinusCharge) delSuMinus-=1;
                }
        }
        const int ntypes=5;
        int del[]={delPlus,delMinus,delZero,delSuPlus,delSuMinus};
        double charges[]={plusCharge,minusCharge,0.0,suPlusCharge,suMinusCharge};
        fprintf(myerr,"setSuRelPot: placing su.s lead to delPlus,delMinus,delZero,delSuPlus,delSuMinus: ");
        for(int i=0;i<ntypes;++i) fprintf(myerr,"%d,",del[i]);
        fprintf(myerr,"\n");

        //neutralize the lattice if necessary
        extern int N;
        int n[]={0,0,0,0,0};
        for(int i=0;i<N;++i) {
                for(int j=0;j<ntypes;++j) {
                        if(*c[i].val==charges[j]) {
                                n[j]+=1;
                                break;
                        }
                }
        }
        fprintf(myerr,"setSuRelPot: nPlus,nMinus,nZero,nSuPlus,nSuMinus: ");
        for(int i=0;i<ntypes;++i) fprintf(myerr,"%d,",n[i]);
        fprintf(myerr,"\n");
        const uint32_t maxIter=10u*(uint32_t)ceil(N*(1.0+log(N)));
        for(uint32_t i=0;i<maxIter;++i) {
                bool endloop=true;
                for(int j=0;j<ntypes;++j) {
                        if(del[j]!=0) {
                                endloop=false;
                                break;
                        }
                }
                if(endloop==true) {
                        fprintf(myerr,"setSuRelPos neut loop received endloop==true. i: %u\n",i);
                        break;
                }
                //use RNG index to remove pos. bias for refill
                const uint32_t rand=pcg32_boundedrand_r(rng,N);
                if(c[rand].su==-1 && c[rand].hydrophobic==false) {
                        neutralizeStep(c,rand,del,charges,ntypes);
                }
                if(i==maxIter-1u) fprintf(myerr,"setSuRelPos reached last step in neut loop");
        }
        fprintf(myerr,"setSuRelPot: post neutralization delPlus,delMinus,delZero,delSuPlus,delSuMinus: ");
        for(int i=0;i<ntypes;++i) fprintf(myerr,"%d,",del[i]);
        fprintf(myerr,"\n");
        for(int j=0;j<ntypes;++j) {
                n[j]=0;
        }
        for(int i=0;i<N;++i) {
                for(int j=0;j<ntypes;++j) {
                        if(*c[i].val==charges[j]) {
                                n[j]+=1;
                                break;
                        }
                }
        }
        fprintf(myerr,"setSuRelPot: nPlus,nMinus,nZero,nSuPlus,nSuMinus: ");
        for(int i=0;i<ntypes;++i) fprintf(myerr,"%d,",n[i]);
        fprintf(myerr,"\n");
        setCharge(c,coreNum);
        double ret=0.0;
        for(int i=0;i<ntypes;++i) ret+=n[i]*charges[i];
        //for(int i=0;i<ntypes;++i) ret+=del[i]*charges[i];
        return (int)round(ret);
}

//each type goes to next in line
//(ie. + -> -, - -> 0, ...)
//rather than set order
//(ie. +,-,0,..)
//to further randomize IC
//del = [delPlus,delMinus,delZero,delSuPlus,delSuMinus,..]
//charges = [plusCharge,minusCharge,0.0,suPlusCharge,suMinusCharge]
void neutralizeStep(lattice c, int rand, int* del, double* charges, int ntypes) {
        const double v=*c[rand].val;
        int type=-1;
        for(int i=0;i<ntypes;++i) {
                if(v==charges[i]) {
                        type=i;
                        break;
                }
        }
        assert(type!=-1);
        if(del[type]>0) {
                for(int i=1;i<ntypes;++i) {
                        const int j=(type+i)%ntypes;
                        if(del[j]<0) {
                                *c[rand].val=charges[j];
                                del[j]+=1;
                                del[type]-=1;
                                break;
                        }
                }
        }
}

/* take COM & rel. pos.s, and find
 * positions on lattice
 */
/* com either double* of length d or
 * NULL; NULL -> use s[ind].com as com
 */
void updSuCurPos(lattice c, su s, int ind, double* com) {
        extern int d;
        int ipos[d];
        if(com==NULL) com=s[ind].com;
        for(int i=0;i<s[ind].nsites;++i) { //update solute pos.s
                double* rotatedRelPos=matrixOnVector(s[ind].orientation,s[ind].relPos[i]);
                for(int e=0;e<d;++e) ipos[e]=lround(com[e]+rotatedRelPos[e]);
                s[ind].currPos[i]=getSiteLat(c,ipos);
                free(rotatedRelPos);
        }
        for(int i=0;i<s[ind].ninshl;++i) { //update shell pos.s
                double* rotatedRelPos=matrixOnVector(s[ind].orientation,s[ind].inShlRelPos[i]);
                for(int e=0;e<d;++e) ipos[e]=lround(com[e]+rotatedRelPos[e]);
                s[ind].inShlCurrPos[i]=getSiteLat(c,ipos);
                free(rotatedRelPos);
        }
}

void updLatSuStatus(lattice c, su s) {
        extern int N,Nsu;
        for(int i=0;i<N;++i) c[i].su=-1;
        for(int i=0;i<Nsu;++i) {
                for(int j=0;j<s[i].nsites;++j) {
                        s[i].currPos[j]->su=i;
                }
        }
}

void updLatHydrophobicStatus(lattice c, su s) {
        extern int N,Nsu;
        for(int i=0;i<N;++i) c[i].hydrophobic=false;
        for(int i=0;i<Nsu;++i) {
                for(int j=0;j<s[i].nsites;++j) {
                        const bool hydr=s[i].hydrophobic[j];
                        if(hydr==true) {
                                s[i].currPos[j]->hydrophobic=hydr;
                        }
                }
        }
}

/* wrapper for updSoluCurrPos & updLattSolStatus */
void updSuCurPos_latSuStatus(lattice c, su s, int ind) {
        updSuCurPos(c,s,ind,(double*)NULL);
        updLatSuStatus(c,s);
}

bool setHydrophobicExist(su s) {
        extern int Nsu;
        bool ret=false;
        for(int i=0;i<Nsu;++i) {
                for(int j=0;j<s[i].nsites;++j) {
                        const bool hij=s[i].hydrophobic[j];
                        if(hij==false) {
                                ;
                        }
                        else {
                                ret=true;
                                break;
                        }
                }
                if(ret==true) {
                        break;
                }
        }
        return ret;
}

/* check solutes are non-overlapping with s[ind]; 
 * returns 0 if no overlap, 1 if overlap 
 */
/* inshlopt: how to deal with overlap
 *      of solute inner shells. Current options:
 ** 0: overlap of inner shells allowed
 ** 1: overlap of inner shells prohibited
 */
/* two-stage comparison:
 ** heuristic: quickly skip solutes that
 **     have large separation (defined by
 **     closest approach of d-cubes along
 **     diagonal)
 ** algorithmic: sort-search to confirm
 **     solutes do not share the same site
 */
int checkSuOverlap(su s, int ind, int inshlopt) {
        extern int Nsu,d,shlL;
        if(Nsu<=1) return 0;
        //rigorous, general comparison performed via sort-search: setting up sorted array
        lattice* sindsites=malloc((s[ind].nsites+inshlopt*s[ind].ninshl)*sizeof(*sindsites));
        assert(sindsites!=NULL /*malloc*/);
        for(int i=0;i<s[ind].nsites;++i)          sindsites[i]=s[ind].currPos[i];
        for(int i=0;i<inshlopt*s[ind].ninshl;++i) sindsites[i+s[ind].nsites]=s[ind].inShlCurrPos[i];
        qsort(sindsites,s[ind].nsites+inshlopt*s[ind].ninshl,sizeof(lattice),(int (*)(const void*,const void*))lattCmp);

        double r2,maxsize,minsize;
        double rvec[d];
        for(int i=0;i<Nsu;++i) {
                if(i==ind) continue;
                //heuristic goes here
                //if still here, need to do a more rigorous check: sort-search
                for(int j=0;j<s[i].nsites;++j) {
                        void* bsearchRet=bsearch((s[i].currPos+j),sindsites,(s[ind].nsites+inshlopt*s[ind].ninshl),
                                                 sizeof(lattice),(int (*)(const void*,const void*))lattCmp);
                        if(bsearchRet!=NULL) {
                                free(sindsites);
                                return 1;
                        }
                }
                for(int j=0;j<inshlopt*s[i].ninshl;++j) {
                        void* bsearchRet=bsearch((s[i].inShlCurrPos+j),sindsites,(s[ind].nsites+inshlopt*s[ind].ninshl),
                                                 sizeof(lattice),(int (*)(const void*,const void*))lattCmp);
                        if(bsearchRet!=NULL) {
                                free(sindsites);
                                return 1;
                        }
                }
        }
        free(sindsites);
        return 0;
}
                /*
                //heuristic to bail out early if solutes are far apart or very close
                //related to dist. of closest approach for cubes along a diagonal
                //see pg.14 NBL lab notes 1608-
                r2=0.0;
                for(int e=0;e<d;++e)
                {
                        rvec[e]=pbc((s[ind].com[e]-s[i].com[e]),e);
                        r2+=rvec[e]*rvec[e];
                }
                if(strcmp(s[i].shape,"cube")==0)
                {
                        minsize=(double)(s[ind].linDim[0]+s[i].linDim[0]+4*inshlopt*shlL);
                        maxsize=minsize;
                }
                if(strcmp(s[i].shape,"rectangle")==0 || strcmp(s[i].shape,"rect")==0)
                {
                        const int minLinDimi=intArrMin(s[i].linDim,d);
                        const int maxLinDimi=intArrMax(s[i].linDim,d);
                        const int minLinDimind=intArrMin(s[ind].linDim,d);
                        const int maxLinDimind=intArrMax(s[ind].linDim,d);
                        minsize=(double)(minLinDimi+minLinDimind+4*inshlopt*shlL);
                        maxsize=(double)(maxLinDimi+maxLinDimind+4*inshlopt*shlL);
                }
                minsize=minsize*minsize*1.0/4.0;
                maxsize=maxsize*maxsize;
                if(4.0*r2 > d*maxsize)  continue;
                else if(r2 < minsize)
                {
                        free(sindsites);
                        return 1;
                }
                */

/*******************************/
/********     core      ********/
/*******************************/

//lc length of c or -1 -> use N as lc
//dZero desired number of zeros in system
//      or -1 if use nZeroSites
//dPlus, dMinus, dZero desired quantities or NULL
//      if use default, ie. extern int nPlusSites
//sending lc=0 effectively skips checkCharge()
//the chgOrig is a hack to allow for use of this
//function for maintaining spins on solute moves
int maintainChargeNeut(lattice c, uint32_t* ii, int lii, int chgSwitch, int* dPlus, int* dMinus, int* dZero, int coreNum, pcg32_random_t* rng, FILE* myerr) {
        extern int N;
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites;
        assert(lii<=N);
        extern double plusCharge,minusCharge;
        const uint32_t maxIter=10u*(uint32_t)ceil(lii*(1.0+log(lii)));
        double* chg;
        double* chgOrig=NULL;
        if(chgSwitch==1)      chg=checkCharge(c,N,NULL);
        else if(chgSwitch==0) chg=checkCharge(c,0,NULL); //returns all 0's
        int delPlus,delMinus,delZero;
        if(dPlus==NULL)  delPlus=chg[0]-nPlusSvSites[coreNum];
        else             delPlus=chg[0]-(*dPlus);
        if(dMinus==NULL) delMinus=chg[1]-nMinusSvSites[coreNum];
        else             delMinus=chg[1]-(*dMinus);
        if(dZero==NULL)  delZero=chg[2]-nZeroSvSites[coreNum];
        else             delZero=chg[2]-(*dZero);

        if(delZero<0) {
                for(uint32_t i=0u;i<maxIter;++i) {
                        if(delZero==0) break;
                        const uint32_t rand=pcg32_boundedrand_r(rng,lii);
                        const int ind=ii[rand];
                        if(c[ind].su!=-1 || c[ind].hydrophobic==true) continue;
                        const double tempv=*c[ind].val;
                        if(tempv==0.0)     continue;
                        else if(tempv>0.0) delPlus-=1;
                        else if(tempv<0.0) delMinus-=1;
                        *c[ind].val=0.0;
                        delZero+=1;
                }
        }
        if(delZero>0) {
                for(uint32_t i=0u;i<maxIter;++i) {
                        if(delZero==0) break;
                        const uint32_t rand=pcg32_boundedrand_r(rng,lii);
                        const int ind=ii[rand];
                        if(c[ind].su!=-1 || c[ind].hydrophobic==true) continue;
                        if(*c[ind].val!=0.0) continue;
                        else if(delMinus<0) {
                                *c[ind].val=minusCharge;
                                delMinus+=1;
                        }
                        else if(delPlus<0) {
                                *c[ind].val=plusCharge;
                                delPlus+=1;
                        }
                        else { //default to delMinus
                                *c[ind].val=minusCharge;
                                delMinus+=1;
                        }
                        delZero-=1;
                }
        }
        if(delZero!=0) {
                fprintf(myerr,"maintainChargeNeut: couldn't get delZero (%d) == 0\n",delZero);
                free(chg);
                return -1;
        }

        int delc=delPlus-delMinus;
        if(delc!=0) {
                if(delc%2 != 0) {
                        fprintf(myerr,"maintainChargeNeut: odd delc:%d (tnplus-nplus)\n",delc);
                        free(chg);
                        return -1;
                }
                const int fromInt = (delc>0) ? 1 : -1;
                const int toInt = -1*fromInt;
                const double from = (delc>0) ? plusCharge : minusCharge;
                const double to = (delc>0) ? minusCharge : plusCharge;
                int maxSwitches=0;
                for(int k=0;k<lii;++k) {
                        const int ind=ii[k];
                        if(c[ind].su!=-1 || c[ind].hydrophobic==true) continue;
                        if(*c[ind].val==from) maxSwitches+=2*fromInt;
                }
                if(fromInt*maxSwitches > fromInt*delc) {
                        for(uint32_t i=0;i<maxIter;++i) {
                                if(delc==0) break;
                                const uint32_t rand=pcg32_boundedrand_r(rng,lii);   //roll a number [0,lc)
                                const int ind=ii[rand];
                                if(c[ind].su!=-1 || c[ind].hydrophobic==true) continue;
                                if(*c[ind].val==from) {    //use RNG index to remove pos. bias for refill
                                        *c[ind].val=to;
                                        delc+=2*toInt;
                                }
                        }
                }
                else if(fromInt*maxSwitches < fromInt*delc) { //fromInt* will make sign +
                        fprintf(myerr,"maintainChargeNeut: best case delc:%d (tnplus-nplus)\n",delc-maxSwitches);
                        free(chg);
                        return -1;
                }
                else if(fromInt*maxSwitches == fromInt*delc) {
                        for(int k=0;k<lii;++k) {
                                const int ind=ii[k];
                                if(c[ind].su!=-1 || c[ind].hydrophobic==true) continue;
                                *c[ind].val=to;
                        }
                }
        }

        chg=checkCharge(c,N,chg);
        delPlus=chg[0]-nPlusSvSites[coreNum];
        delMinus=chg[1]-nMinusSvSites[coreNum];
        delZero=chg[2]-nZeroSvSites[coreNum];
        if(delPlus!=0 || delMinus!=0 || delZero!=0) {
                fprintf(myerr,"maintainChargeNeut: delPlus:%d delMinus:%d delZero:%d\n",delPlus,delMinus,delZero);
                free(chg);
                return -1;
        }
        free(chg);
        return 0;
}

//set lattice to init value: no solutes
void setLatInitSu(lattice c) {
        extern int N;
        for(int i=0;i<N;++i) c[i].su=-1;
}

int willParaLatBeChargeNeut(int coreNum) {
        extern int *nPlusSvSites,*nMinusSvSites,*nPlusSuSites,*nMinusSuSites;
        extern double plusCharge,minusCharge;
        const double totalPlusCharge=(double)(nPlusSvSites[coreNum]+nPlusSuSites[coreNum])*plusCharge;
        const double totalMinusCharge=(double)(nMinusSvSites[coreNum]+nMinusSuSites[coreNum])*minusCharge;
        if(totalPlusCharge==-totalMinusCharge) return 1;
        else                                   return 0;
}

double* setLatVal_rand(lattice c, su s, int coreNum, pcg32_random_t* rng, FILE* myerr) {
        extern int N,Nsu;
        extern int* bc;
        extern double plusCharge,minusCharge;
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites,*nPlusSuSites,*nMinusSuSites,*nHydrSuSites;
        int currPlusSvSites=0,currMinusSvSites=0,currZeroSvSites=0;
        int currPlusSuSites=0,currMinusSuSites=0,currHydrSuSites=0;
        int plusSuSiteNum=0;
        int minusSuSiteNum=0;
        for(int i=0;i<N;++i) {
                const uint32_t rand=pcg32_boundedrand_r(rng,3);  // pick a number [0,3): 0->-1; 1->0; 2->1
                if(rand==0u) {
                        if(currMinusSvSites<nMinusSvSites[coreNum]) {
                                *c[i].val=minusCharge;
                                ++currMinusSvSites;
                        }
                        else if(currZeroSvSites<nZeroSvSites[coreNum]) {
                                *c[i].val=0.0;
                                ++currZeroSvSites;
                        }
                        else if(currPlusSuSites<nPlusSuSites[coreNum]) {
                                int sumPrevSuSites=0;
                                bool gotSite=false;
                                for(int j=0;j<Nsu;++j) {
                                        if(s[j].surfType>0.0) {
                                                if(plusSuSiteNum<s[j].nsites+sumPrevSuSites) {
                                                        *c[i].val=s[j].surfType;
                                                        //c[i].su=j;
                                                        plusSuSiteNum+=1;
                                                        gotSite=true;
                                                        ++currPlusSuSites;
                                                        break;
                                                }
                                                else {
                                                        sumPrevSuSites+=s[j].nsites;
                                                }
                                        }
                                }
                                if(gotSite==false) {
                                        fprintf(stderr,"setLatVal_rand: couldn't get su surfType. exiting\n");
                                        exit(1);
                                }
                        }
                        else if(currMinusSuSites<nMinusSuSites[coreNum]) {
                                int sumPrevSuSites=0;
                                bool gotSite=false;
                                for(int j=0;j<Nsu;++j) {
                                        if(s[j].surfType>0.0) {
                                                if(minusSuSiteNum<s[j].nsites+sumPrevSuSites) {
                                                        *c[i].val=s[j].surfType;
                                                        //c[i].su=j;
                                                        minusSuSiteNum+=1;
                                                        gotSite=true;
                                                        ++currMinusSuSites;
                                                        break;
                                                }
                                                else {
                                                        sumPrevSuSites+=s[j].nsites;
                                                }
                                        }
                                }
                                if(gotSite==false) {
                                        fprintf(stderr,"setLatVal_rand: couldn't get su surfType. exiting\n");
                                        exit(1);
                                }
                        }
                        else if(currHydrSuSites<nHydrSuSites[coreNum]) {
                                *c[i].val=plusCharge;
                                c[i].su=0;
                                c[i].hydrophobic=true;
                                ++currHydrSuSites;
                        }
                        else if(currPlusSvSites<nPlusSvSites[coreNum]) {
                                *c[i].val=plusCharge;
                                ++currPlusSvSites;
                        }
                        else {
                                fprintf(myerr,"setLatVal_rand: ran out of sites to place\n"
                                              "i,c+sv,c-sv,c0sv,c+su,c-su,chsu:%d,%d,%d,%d,%d,%d,%d\n"
                                              "c+sv,c-sv,c0sv,c+su,c-su,chsu:%d,%d,%d,%d,%d,%d\n",
                                        i,currPlusSvSites,currMinusSvSites,currZeroSvSites,currPlusSuSites,
                                        currMinusSuSites,currHydrSuSites,nPlusSvSites[coreNum],nMinusSvSites[coreNum],
                                        nZeroSvSites[coreNum],nPlusSuSites[coreNum],nMinusSuSites[coreNum],nHydrSuSites[coreNum]);
                                exit(1);
                        }
                }
                else if(rand==1u) {
                        if(currZeroSvSites<nZeroSvSites[coreNum]) {
                                *c[i].val=0.0;
                                ++currZeroSvSites;
                        }
                        else if(currPlusSuSites<nPlusSuSites[coreNum]) {
                                int sumPrevSuSites=0;
                                bool gotSite=false;
                                for(int j=0;j<Nsu;++j) {
                                        if(s[j].surfType>0.0) {
                                                if(plusSuSiteNum<s[j].nsites+sumPrevSuSites) {
                                                        *c[i].val=s[j].surfType;
                                                        //c[i].su=j;
                                                        plusSuSiteNum+=1;
                                                        gotSite=true;
                                                        ++currPlusSuSites;
                                                        break;
                                                }
                                                else {
                                                        sumPrevSuSites+=s[j].nsites;
                                                }
                                        }
                                }
                                if(gotSite==false) {
                                        fprintf(stderr,"setLatVal_rand: couldn't get su surfType. exiting\n");
                                        exit(1);
                                }
                        }
                        else if(currMinusSuSites<nMinusSuSites[coreNum]) {
                                int sumPrevSuSites=0;
                                bool gotSite=false;
                                for(int j=0;j<Nsu;++j) {
                                        if(s[j].surfType>0.0) {
                                                if(minusSuSiteNum<s[j].nsites+sumPrevSuSites) {
                                                        *c[i].val=s[j].surfType;
                                                        //c[i].su=j;
                                                        minusSuSiteNum+=1;
                                                        gotSite=true;
                                                        ++currMinusSuSites;
                                                        break;
                                                }
                                                else {
                                                        sumPrevSuSites+=s[j].nsites;
                                                }
                                        }
                                }
                                if(gotSite==false) {
                                        fprintf(stderr,"setLatVal_rand: couldn't get su surfType. exiting\n");
                                        exit(1);
                                }
                        }
                        else if(currHydrSuSites<nHydrSuSites[coreNum]) {
                                *c[i].val=plusCharge;
                                c[i].su=0;
                                c[i].hydrophobic=true;
                                ++currHydrSuSites;
                        }
                        else if(currPlusSvSites<nPlusSvSites[coreNum]) {
                                *c[i].val=plusCharge;
                                ++currPlusSvSites;
                        }
                        else if(currMinusSvSites<nMinusSvSites[coreNum]) {
                                *c[i].val=minusCharge;
                                ++currMinusSvSites;
                        }
                        else {
                                fprintf(myerr,"setLatVal_rand: ran out of sites to place\n"
                                              "i,c+sv,c-sv,c0sv,c+su,c-su,chsu:%d,%d,%d,%d,%d,%d,%d\n"
                                              "c+sv,c-sv,c0sv,c+su,c-su,chsu:%d,%d,%d,%d,%d,%d\n",
                                        i,currPlusSvSites,currMinusSvSites,currZeroSvSites,currPlusSuSites,
                                        currMinusSuSites,currHydrSuSites,nPlusSvSites[coreNum],nMinusSvSites[coreNum],
                                        nZeroSvSites[coreNum],nPlusSuSites[coreNum],nMinusSuSites[coreNum],nHydrSuSites[coreNum]);
                                exit(1);
                        }
                }
                else if(rand==2u) {
                        if(currPlusSvSites<nPlusSvSites[coreNum]) {
                                *c[i].val=plusCharge;
                                ++currPlusSvSites;
                        }
                        else if(currMinusSvSites<nMinusSvSites[coreNum]) {
                                *c[i].val=minusCharge;
                                ++currMinusSvSites;
                        }
                        else if(currZeroSvSites<nZeroSvSites[coreNum]) {
                                *c[i].val=0.0;
                                ++currZeroSvSites;
                        }
                        else if(currPlusSuSites<nPlusSuSites[coreNum]) {
                                int sumPrevSuSites=0;
                                bool gotSite=false;
                                for(int j=0;j<Nsu;++j) {
                                        if(s[j].surfType>0.0) {
                                                if(plusSuSiteNum<s[j].nsites+sumPrevSuSites) {
                                                        *c[i].val=s[j].surfType;
                                                        //c[i].su=j;
                                                        plusSuSiteNum+=1;
                                                        gotSite=true;
                                                        ++currPlusSuSites;
                                                        break;
                                                }
                                                else {
                                                        sumPrevSuSites+=s[j].nsites;
                                                }
                                        }
                                }
                                if(gotSite==false) {
                                        fprintf(stderr,"setLatVal_rand: couldn't get su surfType. exiting\n");
                                        exit(1);
                                }
                        }
                        else if(currMinusSuSites<nMinusSuSites[coreNum]) {
                                int sumPrevSuSites=0;
                                bool gotSite=false;
                                for(int j=0;j<Nsu;++j) {
                                        if(s[j].surfType>0.0) {
                                                if(minusSuSiteNum<s[j].nsites+sumPrevSuSites) {
                                                        *c[i].val=s[j].surfType;
                                                        //c[i].su=j;
                                                        minusSuSiteNum+=1;
                                                        gotSite=true;
                                                        ++currMinusSuSites;
                                                        break;
                                                }
                                                else {
                                                        sumPrevSuSites+=s[j].nsites;
                                                }
                                        }
                                }
                                if(gotSite==false) {
                                        fprintf(stderr,"setLatVal_rand: couldn't get su surfType. exiting\n");
                                        exit(1);
                                }
                        }
                        else if(currHydrSuSites<nHydrSuSites[coreNum]) {
                                *c[i].val=plusCharge;
                                c[i].su=0;
                                c[i].hydrophobic=true;
                                ++currHydrSuSites;
                        }
                        else {
                                fprintf(myerr,"setLatVal_rand: ran out of sites to place\n"
                                              "i,c+sv,c-sv,c0sv,c+su,c-su,chsu:%d,%d,%d,%d,%d,%d,%d\n"
                                              "c+sv,c-sv,c0sv,c+su,c-su,chsu:%d,%d,%d,%d,%d,%d\n",
                                        i,currPlusSvSites,currMinusSvSites,currZeroSvSites,currPlusSuSites,
                                        currMinusSuSites,currHydrSuSites,nPlusSvSites[coreNum],nMinusSvSites[coreNum],
                                        nZeroSvSites[coreNum],nPlusSuSites[coreNum],nMinusSuSites[coreNum],nHydrSuSites[coreNum]);
                                exit(1);
                        }
                }
                else exit(1);
        }
        fprintf(myerr,"setLatVal_rand: rng gave nSvPlus,nSvMinus,nSvZero,nSuPlus,nSuMinus,nSuHydr: %d,%d,%d,%d,%d,%d "
                      "N:%d nSvPlus/N,nSvMinus/N,nSvZero/N,nSuPlus/N,nSuMinus/N,nSuHydr/N:%f,%f,%f,%f,%f,%f\n"
                      "input from .para file nSvPlus,nSvMinus,nSvZero,nSuPlus,nSuMinus,nSuHydr: %d,%d,%d,%d,%d,%d, with +,-: %f,%f.\n",
                      currPlusSvSites,currMinusSvSites,currZeroSvSites,currPlusSuSites,currMinusSuSites,currHydrSuSites,
                      N,(double)(currPlusSvSites)/N,(double)(currMinusSvSites)/N,(double)(currZeroSvSites)/N,(double)(currPlusSuSites)/N,
                      (double)(currMinusSuSites)/N,(double)(currHydrSuSites)/N,nPlusSvSites[coreNum],nMinusSvSites[coreNum],
                      nZeroSvSites[coreNum],nPlusSuSites[coreNum],nMinusSuSites[coreNum],nHydrSuSites[coreNum],plusCharge,minusCharge);

        return NULL;
}

double* setLatVal_data(lattice c, char** lines, int nlines, int coreNum, FILE* myerr) {
        extern int N;
        int i=0;
        while( strcmp(lines[i],"Cells") != 0 ) { 
                ++i;
                assert(i<nlines /* couldnt find Cells header*/);
        }
        ++i;
        char** words=NULL;
        int nwords[1]={0};
        for(int j=0; j<N; ++j) {
                assert((i+j)<nlines);
                words=split(lines[i+j], "\t", nwords);
                if(*nwords==1) {
                        fprintf(myerr,"lines[i+j-1]:%s\nlines[i+j]:%s\nwords[0]:%s\ti:%d\tj:%d\tnlines:%d\n",lines[i+j-1],lines[i+j],words[0],i,j,nlines);
                }
                *c[j].val=atof(words[1]);
                free(*words);   free(words);    words=NULL;
        }
        setCharge(c,coreNum);
        return checkCharge(c,N,NULL);
}

//equivalent to rho_xyz = cos(pi x)cos(pi y)cos(pi z)
//FT has only one non-zero term: k=(L/2,L/2,L/2)
double* setLatVal_simpCubic(lattice c, int coreNum) {
        extern int N,d;
        extern int* L;
        assert(d==3 /*setLatVal_simpCubic*/);
        assert(N==L[0]*L[1]*L[2] /*setLatVal_simpCubic*/);
        int ct=0;
        double spin=1.0;
        for(int i=0;i<L[0];++i) {
                spin*=-1.0;
                for(int j=0;j<L[1];++j) {
                        spin*=-1.0;
                        for(int k=0;k<L[2];++k) {
                                spin*=-1.0;
                                assert(spin==1.0 || spin==-1.0 /*setLatVal_simpCubic*/);
                                *c[ct].val=spin;
                                ++ct;
                        }
                }
        }
        setCharge(c,coreNum);
        return checkCharge(c,N,NULL);
}

//equivalent to rho_xyz = cos(pi x)cos(pi y)cos(pi z)
//FT has only one non-zero term: k=(L/2,L/2,L/2)
double* setLatVal_simpCubicWDefect(lattice c, int coreNum, pcg32_random_t* rng) {
        extern int N,d;
        extern int* L;
        assert(d==3 /*setLatVal_simpCubic*/);
        assert(N==L[0]*L[1]*L[2] /*setLatVal_simpCubic*/);
        int ct=0;
        double spin=1.0;
        for(int i=0;i<L[0];++i) {
                spin*=-1.0;
                for(int j=0;j<L[1];++j) {
                        spin*=-1.0;
                        for(int k=0;k<L[2];++k) {
                                spin*=-1.0;
                                assert(spin==1.0 || spin==-1.0 /*setLatVal_simpCubic*/);
                                *c[ct].val=spin;
                                ++ct;
                        }
                }
        }

        for(int i=0;i<2;++i) {
                //introduce one defect
                const uint32_t ind0=pcg32_boundedrand_r(rng,N);
                const int ind1=getWrappedNN((int)ind0,1,NULL);
                const int tmp=*c[ind1].val;
                *c[ind1].val=*c[ind0].val;
                *c[ind0].val=tmp;
        }

        setCharge(c,coreNum);
        return checkCharge(c,N,NULL);
}

void** setLatNeigh_coul_twoway(lattice c, para p, FILE* myerr) {
        extern int N,d;
        int totCoulNeigh=0;
        #if COUL_ON
        const double Rc=(p->sigma)*(p->coulCutoff);
        #endif
        for(int i=0;i<N;++i) { //get totNN
                #if COUL_ON
                if(Rc>=sqrt(2.0)) {
                        for(int j=0;j<N;++j) {
                                if(i==j) continue;
                                const double r=dist_lat(i,j);
                                if((r <= Rc) && (r > 1.0)) ++totCoulNeigh;
                        }
                }
                #endif
                c[i].nearNeigh=NULL;
                c[i].coulNeigh=NULL;
                c[i].nCoulNeigh=0;
                c[i].dInvErfcNeigh=NULL;
        }
        fprintf(myerr,"setLatNeigh_coul: NNs: %d\tcoul neighbors:%d\ttotal neighbors (two-way): %d\n",2*d*N,totCoulNeigh,2*d*N+totCoulNeigh);

        //init
        const int NNN=2*d;       //number nearest neighbors is 2 * d
        lattice* neighHolder=malloc((NNN*N+totCoulNeigh)*sizeof(*neighHolder));
        assert(neighHolder!=NULL /*malloc*/);
        double* distInvErfcHolder=NULL;
        if(totCoulNeigh!=0) {
                distInvErfcHolder=malloc(totCoulNeigh*sizeof(*distInvErfcHolder));
                assert(distInvErfcHolder!=NULL /*malloc*/);
        }
        #if COUL_ON
        int** cube=NULL;
        int* pos=NULL;
        int* neighPos=NULL;
        int* origin=NULL;
        int Ncube=1;
        if(totCoulNeigh!=0) {
                for(int e=0;e<d;++e) Ncube*=(2*(int)ceil(Rc)+1);
                cube=buildNeighCube(Rc);
                pos=malloc(d*sizeof(*pos));
                neighPos=malloc(d*sizeof(*pos));
                origin=malloc(d*sizeof(*pos));
                assert(pos!=NULL && neighPos!=NULL && origin!=NULL /*malloc*/);
                for(int e=0;e<d;++e) origin[e]=0;
        }
        #endif
        extern double pi;
        int coulNeighCt=0;
        for(int i=0;i<N;++i) {
                c[i].nearNeigh=(neighHolder+i*NNN+coulNeighCt);  //first NNN of each neigh block are NN
                for(int e=-d;e<d;++e) { // [-d,-1] neg. dir., [0,d) pos. w/ -1:0,...,-d:(d-1)
                        neighHolder[(e+d)+i*NNN+coulNeighCt]=getNNLat(c,i,e);
                }
                #if COUL_ON
                if(totCoulNeigh!=0) {
                        getPos(i,pos);
                        for(int j=0;j<Ncube;++j) {
                                const double r=dist_intvect_nopbc(cube[j],origin);
                                if((r <= Rc) && (r > 1.0)) {
                                        for(int e=0;e<d;++e) neighPos[e]=pos[e]+cube[j][e];
                                        const int neigh=getWrappedSite(neighPos,NULL);


                                        if(c[i].nCoulNeigh==0) {
                                                c[i].coulNeigh=(neighHolder+(i+1)*NNN+coulNeighCt);  //first NNN of each neigh block are NN
                                                c[i].dInvErfcNeigh=(distInvErfcHolder+coulNeighCt);
                                        }
                                        neighHolder[(i+1)*NNN+coulNeighCt]=(c+neigh);    //first NNN of each neigh block are NN
                                        #if SHIFT_COUL
                                        distInvErfcHolder[coulNeighCt]=erfc(r/(p->sigma))/r - p->e_cut;  //truncated&shifted
                                        #endif
                                        #if !SHIFT_COUL
                                        distInvErfcHolder[coulNeighCt]=erfc(r/(p->sigma))/r;             //truncated
                                        #endif
                                        ++coulNeighCt;
                                        ++(c[i].nCoulNeigh);
                                }
                        }
                }
                #endif
        }
        #if COUL_ON
        if(totCoulNeigh!=0) {
                free(*cube);
                free(cube);
                free(pos);
                free(neighPos);
                free(origin);
        }
        if(totCoulNeigh!=coulNeighCt) {
                fprintf(myerr,"setLatNeigh_coul_twoway(): totCoulNeigh!=coulNeighCt: %d!=%d\n",totCoulNeigh,coulNeighCt);
        }
        #endif
        void** ptrs=malloc(2*sizeof(*ptrs));
        assert(ptrs!=NULL /*malloc*/);
        ptrs[0]=(void*)neighHolder;
        ptrs[1]=(void*)distInvErfcHolder;
        return ptrs;
}

double* setupEwald(int* L, double sigma, int nCoresP) {
        fprintf(stderr,"setupEwald: fn called\n");
        fflush(stderr);
        if(L==NULL) {
                extern int* L;
        }
        const int d=3;
        assert(d==3 /*setupEwald assumes...*/);
        const int N=L[0]*L[1]*L[2];
        //const double err=10e-10; //hardcode
        const double pi=acos(-1);
        const double s2=sigma*sigma;
        double* ewald=malloc(N*sizeof(*ewald));
        assert(ewald!=NULL);
        double complex* cexpValsHolderX=malloc(L[0]*L[0]*sizeof(*cexpValsHolderX));
        double complex* cexpValsHolderY=malloc(L[1]*L[1]*sizeof(*cexpValsHolderY));
        double complex* cexpValsHolderZ=malloc(L[2]*L[2]*sizeof(*cexpValsHolderZ));
        assert(cexpValsHolderX!=NULL && cexpValsHolderY!=NULL && cexpValsHolderZ!=NULL);
        double complex** cexpValsX=malloc(L[0]*sizeof(*cexpValsX)); //add precompute of vectors
        double complex** cexpValsY=malloc(L[1]*sizeof(*cexpValsY));
        double complex** cexpValsZ=malloc(L[2]*sizeof(*cexpValsZ));
        assert(cexpValsX!=NULL && cexpValsY!=NULL && cexpValsZ!=NULL);
        int i=0;
        for(int k=0;k<L[0];++k) {
                cexpValsX[k]=cexpValsHolderX+i;
                for(int x=-L[0]/2;x<L[0]/2;++x) {
                        cexpValsX[k][x+L[0]/2]=cexp(2.0*pi*I/L[0]*k*x);
                        ++i;
                }
        }
        i=0;
        for(int k=0;k<L[1];++k) {
                cexpValsY[k]=cexpValsHolderY+i;
                for(int y=-L[1]/2;y<L[1]/2;++y) {
                        cexpValsY[k][y+L[1]/2]=cexp(2.0*pi*I/L[1]*k*y);
                        ++i;
                }
        }
        i=0;
        for(int k=0;k<L[2];++k) {
                cexpValsZ[k]=cexpValsHolderZ+i;
                for(int z=-L[2]/2;z<L[2]/2;++z) {
                        cexpValsZ[k][z+L[2]/2]=cexp(2.0*pi*I/L[2]*k*z);
                        ++i;
                }
        }
        fprintf(stderr,"setupEwald: pre main loop\n");
        fflush(stderr);
        const int j=0;
        #pragma omp parallel for num_threads(nCoresP)
        for(int l=j+1;l<N;++l) {
                int rjl[d];
                getPos(l,rjl);
                minImageVect(rjl,L);
                const int x=rjl[0]+L[0]/2;
                const int y=rjl[1]+L[1]/2;
                const int z=rjl[2]+L[2]/2;
                ewald[l]=0.0;
                double cx,cy,cz;
                for(int kx=0;kx<L[0];++kx) {
                        const double kx2=(double)(kx*kx)/(L[0]*L[0]);
                        if(kx!=0) cx=1.0;
                        else      cx=0.5;
                        for(int ky=0;ky<L[1];++ky) {
                                const double ky2=(double)(ky*ky)/(L[1]*L[1]);
                                if(ky!=0) cy=1.0;
                                else      cy=0.5;
                                for(int kz=0;kz<L[2];++kz) {
                                        if(kx!=0 || ky!=0 || kz!=0) {
                                                const double kz2=(double)(kz*kz)/(L[2]*L[2]);
                                                if(kz!=0) cz=1.0;
                                                else      cz=0.5;
                                                const double k2=4.0*pi*pi*(kx2+ky2+kz2);
                                                //use binary numbers to repr permutations
                                                //0-> +1
                                                //1-> -1
                                                const double complex ekr000=cexpValsX[kx][x]*
                                                                            cexpValsY[ky][y]*
                                                                            cexpValsZ[kz][z];
                                                const double complex ekr001=cexpValsX[kx][x]*
                                                                            cexpValsY[ky][y]*
                                                                            conj(cexpValsZ[kz][z]);
                                                const double complex ekr010=cexpValsX[kx][x]*
                                                                            conj(cexpValsY[ky][y])*
                                                                            cexpValsZ[kz][z];
                                                const double complex ekr011=cexpValsX[kx][x]*
                                                                            conj(cexpValsY[ky][y])*
                                                                            conj(cexpValsZ[kz][z]);
                                                const double Ek=exp(-k2*s2/4.0)/k2*creal(ekr000+ekr001+ekr010+ekr011);
                                                ewald[l]+=2.0*cx*cy*cz*Ek;
                                        }
                                }
                        }
                }
        }
        fprintf(stderr,"setupEwald: post main loop\n");
        fflush(stderr);
        //diag element: 'self E'
        const double x=0.0;
        const double y=0.0;
        const double z=0.0;
        ewald[0]=0.0;
        double cx,cy,cz;
        for(int kx=0;kx<L[0];++kx) {
                const double kx2=(double)(kx*kx)/(L[0]*L[0]);
                const double xkx=x*kx/L[0];
                if(kx!=0) cx=1.0;
                else      cx=0.5;
                for(int ky=0;ky<L[1];++ky) {
                        const double ky2=(double)(ky*ky)/(L[1]*L[1]);
                        const double yky=y*ky/L[1];
                        if(ky!=0) cy=1.0;
                        else      cy=0.5;
                        for(int kz=0;kz<L[2];++kz) {
                                if(kx!=0 || ky!=0 || kz!=0) {
                                        const double kz2=(double)(kz*kz)/(L[2]*L[2]);
                                        const double zkz=z*kz/L[2];
                                        if(kz!=0) cz=1.0;
                                        else      cz=0.5;
                                        const double k2=4.0*pi*pi*(kx2+ky2+kz2);
                                        //use binary numbers to repr permutations
                                        //0-> +1
                                        //1-> -1
                                        const double kr000=2.0*pi*(xkx+yky+zkz);
                                        const double kr001=2.0*pi*(xkx+yky-zkz);
                                        const double kr010=2.0*pi*(xkx-yky+zkz);
                                        const double kr011=2.0*pi*(xkx-yky-zkz);
                                        const double Ek=exp(-k2*s2/4.0)/k2*
                                                creal(cexp(I*kr000)+cexp(I*kr001)+cexp(I*kr010)+cexp(I*kr011));
                                        ewald[0]+=2.0*cx*cy*cz*Ek;
                                }
                        }
                }
        }
        fprintf(stderr,"setupEwald: work done, pre-cleanup\n");
        fflush(stderr);
        free(cexpValsHolderX); free(cexpValsX);
        free(cexpValsHolderY); free(cexpValsY);
        free(cexpValsHolderZ); free(cexpValsZ);
        fprintf(stderr,"setupEwald: post-cleanup. returning sucessfully\n");
        fflush(stderr);
        return ewald;
}

/* create a local subset of the lattice formed from
 * the neighbors of core; pass core, length of core,
 * and a NULL ptr. subset. subset will be allocated
 * memory and built in the function; caller must
 * free. Returns length of subset or -1 on error.
 */
int localSubset(lattice* core, int lcore, lattice** subset) {
        //compile local subset of core: neighs of core
        extern int d;
        int lsubset=(2*d+1)*lcore; //2*d: NNs; +1: the core elements themselves
        for(int i=0;i<lcore;++i) lsubset+=core[i]->nCoulNeigh;
        *subset=malloc(lsubset*sizeof(**subset));
        assert(*subset!=NULL /*malloc*/);
        int k=0;
        for(int i=0;i<lcore;++i) {
                (*subset)[k++]=core[i];
                for(int e=0;e<2*d;++e)                 (*subset)[k++]=core[i]->nearNeigh[e];
                for(int j=0;j<core[i]->nCoulNeigh;++j) (*subset)[k++]=core[i]->coulNeigh[j];
        }
        while(k<lsubset) (*subset)[k++]=NULL;   //in case there are leftovers
        //sort & trim subset: only one copy of each site allowed
        qsort(*subset,lsubset,sizeof(lattice),(int (*)(const void*,const void*))lattCmp);
        lsubset=uniquelattice(subset,lsubset);
        *subset=realloc(*subset,lsubset*sizeof(**subset));
        assert(*subset!=NULL /*realloc*/);
        return lsubset;
}

lattice* numToLat(lattice c, uint32_t* nums, uint32_t nnums) {
        lattice* lat=malloc(nnums*sizeof(*lat));
        assert(lat!=NULL /*malloc*/);
        for(uint32_t i=0u;i<nnums;++i) lat[i]=c+nums[i];
        return lat;
}

uint32_t* latToNum(lattice c, lattice* sites, int nsites) {
        uint32_t* nums=malloc(nsites*sizeof(*nums));
        assert(nums!=NULL);
        for(int i=0;i<nsites;++i) nums[i]=sites[i]-c;
        return nums;
}

int spinSwapMoveNN(lattice c, lattice tc, su s, double* E, int ind, double* ewald, para p, pcg32_random_t* rng, FILE* myerr) {
        extern int d,N;
        uint32_t ind2=pcg32_boundedrand_r(rng,2*d);
        lattice nn=getNNLat(c,ind,(int)ind2);
        if(nn->su!=-1 || *c[ind].val==*nn->val) return 0;

        //check for swap
        const int lcore=2;
        lattice core[lcore];
        core[0]=c+ind;
        core[1]=nn;
        lattice* subset=NULL;   //subset alloc'd in localSubset
        int lsubset=localSubset(core,lcore,&subset);
        assert(lsubset>0 /*localSubset output*/);

        uint32_t flipped[2];
        flipped[0]=(ind<ind2) ? ind : ind2; //min
        flipped[1]=(ind>ind2) ? ind : ind2; //max
        double Eoldloc[ENERGY_LENGTH];
        energy_local(s,Eoldloc,p,subset,lsubset);
        Eoldloc[6]=(p->Q)*energy_lr_nflip(c,ewald,flipped,2);
        swapD(c[ind].val, nn->val);
        double Enewloc[ENERGY_LENGTH];
        energy_local(s,Enewloc,p,subset,lsubset);
        Enewloc[6]=(p->Q)*energy_lr_nflip(c,ewald,flipped,2);
        double deltaEloc[ENERGY_LENGTH];
        double deltaen=0.0;
        for(int i=0;i<ENERGY_LENGTH;++i) {
                deltaEloc[i]=Enewloc[i]-Eoldloc[i];
                deltaen+=deltaEloc[i];
        }

        #if TEST_BUILD
        //test if local & full give same deltaen:
        double EnewFull[ENERGY_LENGTH];
        const double fen_new=energy_full(c,s,EnewFull,ewald,p);
        swapD(c[ind].val, nn->val);
        double EoldFull[ENERGY_LENGTH];
        const double fen_old=energy_full(c,s,EoldFull,ewald,p);
        swapD(c[ind].val, nn->val);
        if(fabs(fen_new-fen_old-deltaen) >= 0.00001) {
                fprintf(myerr,"spinSwapMove: full deltaen: %f\t\tlocal deltaen:%f\tdiff:%f\n",fen_new-fen_old,deltaen,fen_new-fen_old-deltaen);
                fprintf(myerr,"\tfull\tlocal\n0:\t%f\t%f\n1:\t%f\t%f\n2:\t%f\t%f\n3:\t%f\t%f\n4:\t%f\t%f\n" \
                              "5:\t%f\t%f\n6:\t%f\t%f\n7:\t%f\t%f\n8:\t%f\t%f\n9:\t%f\t%f\n",  \
                        EnewFull[0]-EoldFull[0],deltaEloc[0],EnewFull[1]-EoldFull[1],deltaEloc[1],EnewFull[2]-EoldFull[2],deltaEloc[2], \
                        EnewFull[3]-EoldFull[3],deltaEloc[3],EnewFull[4]-EoldFull[4],deltaEloc[4],EnewFull[5]-EoldFull[5],deltaEloc[5], \
                        EnewFull[6]-EoldFull[6],deltaEloc[6],EnewFull[7]-EoldFull[7],deltaEloc[7],EnewFull[8]-EoldFull[8],deltaEloc[8], \
                        EnewFull[9]-EoldFull[9],deltaEloc[9]       );
                fprintf(myerr,"core spins: %f\t%f\n",*core[0]->val,*core[1]->val);
                //fprintf(myerr,"i\tx,y,z\tval\taddr\tlsubset %d\n",lsubset);
                //for(int i=0;i<lsubset;++i) {
                //        const int num=subset[i]-c;
                //        int tmppos[d];
                //        getPos(num,tmppos);
                //        fprintf(myerr,"%d\t%d,%d,%d\t%f\t%p\n",num,tmppos[0],tmppos[1],tmppos[2],*subset[i]->val,subset[i]);
                //}
                swapD(c[ind].val, nn->val);
                return 0;
        }
        #endif

        free(subset);
        if(deltaen<=0.0) {
                if(E!=NULL) {
                        for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaEloc[i];
                }
                if(tc!=NULL) swapD(tc[ind].val, tc[ind2].val);
                return 1;
        }

        //check Boltz. prob
        const double pr=(double)pcg32_random_r(rng)/UINT32_MAX; //pr in [0,1)
        const double boltz=exp(-(p->beta) * deltaen);
        if(pr<=boltz) {
                if(E!=NULL) {
                        for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaEloc[i];
                }
                if(tc!=NULL) swapD(tc[ind].val, tc[ind2].val);
                return 1;
        }
        else {
                swapD(c[ind].val, nn->val);
                return 0;
        }
}

/* NOTE: contrary to usual convention,
 * 1 => move sucessful; 0 => move rejected 
 *
 * Assumes ind is not in solute (this assumption
 * holds when sv vs su move is picked based
 * on whether the cell is in solute or not: see
 * MCstep)
 */
int spinSwapMove(lattice c, lattice tc, su s, double* E, int ind, double* ewald, para p, pcg32_random_t* rng, FILE* myerr) {
        extern int d,N;
        uint32_t ind2=pcg32_boundedrand_r(rng, N);
        if(c[ind2].su!=-1 || *c[ind].val==*c[ind2].val) return 0;

        //check for swap
        const int lcore=2;
        lattice core[lcore];
        core[0]=c+ind;
        core[1]=c+ind2;
        lattice* subset=NULL;   //subset alloc'd in localSubset
        int lsubset=localSubset(core,lcore,&subset);
        assert(lsubset>0 /*localSubset output*/);

        uint32_t flipped[2]; //=[min(ind,ind2),max(")]
        if(ind<ind2) {
                flipped[0]=ind;
                flipped[1]=ind2;
        }
        else {
                flipped[0]=ind2;
                flipped[1]=ind;
        }
        double Eoldloc[ENERGY_LENGTH];
        energy_local(s,Eoldloc,p,subset,lsubset);
        Eoldloc[6]=(p->Q)*energy_lr_nflip(c,ewald,flipped,2);
        swapD(c[ind].val, c[ind2].val);
        double Enewloc[ENERGY_LENGTH];
        energy_local(s,Enewloc,p,subset,lsubset);
        Enewloc[6]=(p->Q)*energy_lr_nflip(c,ewald,flipped,2);
        double deltaEloc[ENERGY_LENGTH];
        double deltaen=0.0;
        for(int i=0;i<ENERGY_LENGTH;++i) {
                deltaEloc[i]=Enewloc[i]-Eoldloc[i];
                deltaen+=deltaEloc[i];
        }

        #if TEST_BUILD
        //test if local & full give same deltaen:
        double EnewFull[ENERGY_LENGTH];
        const double fen_new=energy_full(c,s,EnewFull,ewald,p);
        swapD(c[ind].val, c[ind2].val);
        double EoldFull[ENERGY_LENGTH];
        const double fen_old=energy_full(c,s,EoldFull,ewald,p);
        swapD(c[ind].val, c[ind2].val);
        if(fabs(fen_new-fen_old-deltaen) >= 0.00001) {
                fprintf(myerr,"spinSwapMove: full deltaen: %f\t\tlocal deltaen:%f\tdiff:%f\n",fen_new-fen_old,deltaen,fen_new-fen_old-deltaen);
                fprintf(myerr,"\tfull\tlocal\n0:\t%f\t%f\n1:\t%f\t%f\n2:\t%f\t%f\n3:\t%f\t%f\n" \
                              "4:\t%f\t%f\n5:\t%f\t%f\n6:\t%f\t%f\n7:\t%f\t%f\n8:\t%f\t%f\n9:\t%f\t%f\n",  \
                        EnewFull[0]-EoldFull[0],deltaEloc[0],EnewFull[1]-EoldFull[1],deltaEloc[1],EnewFull[2]-EoldFull[2],deltaEloc[2], \
                        EnewFull[3]-EoldFull[3],deltaEloc[3],EnewFull[4]-EoldFull[4],deltaEloc[4],EnewFull[5]-EoldFull[5],deltaEloc[5], \
                        EnewFull[6]-EoldFull[6],deltaEloc[6],EnewFull[7]-EoldFull[7],deltaEloc[7],EnewFull[8]-EoldFull[8],deltaEloc[8], \
                        EnewFull[9]-EoldFull[9],deltaEloc[9]);
                fprintf(myerr,"core spins: %f\t%f\n",*core[0]->val,*core[1]->val);
                //fprintf(myerr,"i\tx,y,z\tval\taddr\tlsubset %d\n",lsubset);
                //for(int i=0;i<lsubset;++i) {
                //        const int num=subset[i]-c;
                //        int tmppos[d];
                //        getPos(num,tmppos);
                //        fprintf(myerr,"%d\t%d,%d,%d\t%f\t%p\n",num,tmppos[0],tmppos[1],tmppos[2],*subset[i]->val,subset[i]);
                //}
                swapD(c[ind].val, c[ind2].val);
                return 0;
        }
        #endif

        free(subset);
        if(deltaen<=0.0) {
                if(E!=NULL) {
                        for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaEloc[i];
                }
                if(tc!=NULL) swapD(tc[ind].val, tc[ind2].val);
                return 1;
        }

        //check Boltz. prob
        const double pr=(double)pcg32_random_r(rng)/UINT32_MAX; //pr in [0,1)
        const double boltz=exp(-(p->beta) * deltaen);
        if(pr<=boltz) {
                if(E!=NULL) {
                        for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaEloc[i];
                }
                if(tc!=NULL) swapD(tc[ind].val, tc[ind2].val);
                return 1;
        }
        else {
                swapD(c[ind].val, c[ind2].val);
                return 0;
        }
}

/* make a cluster move as discussed in Grousson, Viot 2001.
 * ind is seed0; find seed1 by randomly selecting a center/
 * pivor and rotation axis, then rotating seed1 \pi about
 * that pivot/axis. If seed0, seed1 have opposite sign,
 * start building clusters.
 * Clusters are built by considering NNs of each current
 * seed. We wish to construct clusters of same spins, and
 * each cluster with the same number of spins, so that a
 * cluster move will not change the net charge of the sys.
 * To avoid excessive numbers of operations to build clusters,
 * we build by adding NNs, rather than all neighbors. Walk
 * through NNs of seeds in random order; if roll p in [0,1)
 * less than 1-exp(-4*beta_eff*J), add the considered NNs to
 * each cluster & to the stack. Once all possible additions
 * are attempted for a given seed, move to the next seed
 * pair in the stack. Once clusters are built, we need to
 * compare energies as usual to determine whether to accept
 * the move. Because of the exotic nature of this cluster
 * move, the probability of acceptance is:
 *
 * 1            if      (\deltaE_1^{ji}+(1-\beta_eff/\beta)*\deltaE_0^{ji})<0
 * exp(-\beta(\deltaE_1^{ji}+(1-\beta_eff/\beta)*\deltaE_0^{ji}))        else
 *
 * where E_0 is the NN part of the Ising Hamiltonian, and
 * E_1 is the other part (here, the screened Coulomb part).
 */

//built for fully occupied lattice (ie. each
//site it +1 or -1). need to rewrite if want
//to expand to case with defects (0's) or
//more general solvent types (+2, -1, or di-
//fferent shapes, etc etc)
int clusterMove(lattice c, lattice tc, su s, double* E, const uint32_t ind0, double* ewald, para p, int coreNum, pcg32_random_t* rng, FILE* myerr) {
        extern int d,N,Nsu;
        extern int* L;
        extern int* bc;

        //boring, low-success, correct way to find pair
        const uint32_t rot=pcg32_boundedrand_r(rng, d);     //roll rotation axis
        uint32_t cent=pcg32_boundedrand_r(rng, N);          //roll pivot
        while(cent==ind0) cent=pcg32_boundedrand_r(rng, N); //get valid center
        int seed0[d];
        int seed1[d];
        int centpos[d];
        getPos(ind0,seed0);
        getPos(cent,centpos);
        for(uint32_t e=0;e<d;++e) {
                seed1[e]=seed0[e];
                if(e!=rot) seed1[e]-=2*(seed0[e]-centpos[e]);
        }
        const uint32_t ind1=getWrappedSite(seed1,NULL);
        if(c[ind1].su!=-1 || *c[ind1].val==*c[ind0].val) return 0;
        uint32_t seedPair[3];
        seedPair[0]=ind0;
        seedPair[1]=ind1;
        seedPair[2]=rot;

        double suBlocks=0.0;
        uint32_t nInCluster=0u;
        uint32_t* inCluster=buildCluster(tc,seedPair[0],seedPair[1],seedPair[2],&nInCluster,&suBlocks,p,coreNum,rng);
        const double E_IsuClustMvConstraint=2.0*(p->J)*suBlocks; //see note @ buildCluster()
        if(nInCluster==0u) return 0;
        double clusterE=0.0;
        double deltaE[ENERGY_LENGTH];
        for(int i=0;i<ENERGY_LENGTH;++i) deltaE[i]=0.0;
        #if !COUL_ON //pure Ising: accept cluster
        double E_IsuOld=0.0;
        double E_IsuNew=0.0;
        if(Nsu>0) {
                lattice* latInCluster=numToLat(c,inCluster,nInCluster);
                lattice* subset=NULL;
                int lsubset=localSubset(latInCluster,nInCluster,&subset);
                E_IsuOld=energy_ising_su_local(s,p,subset,lsubset);
                free(latInCluster); free(subset);
                latInCluster=numToLat(tc,inCluster,nInCluster);
                lsubset=localSubset(latInCluster,nInCluster,&subset);
                E_IsuNew=energy_ising_su_local(s,p,subset,lsubset); //just use s for isimult
                free(latInCluster); free(subset);
                deltaE[3]=E_IsuNew-E_IsuOld+E_IsuClustMvConstraint; //su
                clusterE=deltaE[3];

                ////test if local & full give same deltaen:
                ////double fen_new=energy_full(tc,s,NULL,ffto,fftp,p);
                ////double fen_old=energy_full(c,s,NULL,ffto,fftp,p);
                //double fen_new_isu=energy_ising_su(tc,s,p);
                //double fen_old_isu=energy_ising_su(c,s,p);
                //fprintf(myerr,"clusterMove: local cluster en: %f\tfull isu:%f\n",clusterE,fen_new_isu-fen_old_isu+E_IsuClustMvConstraint);
                ////fprintf(myerr,"clusterMove: local cluster en: %f\tfull clusteren:%f\tfull isu:%f\n",clusterE,fen_new-fen_old,fen_new_isu-fen_old_isu);
                //fprintf(myerr,"(local new,old)(full new,old)(constraint): (%f,%f)(%f,%f)(%f)\tsuBlocks: %f\n",E_IsuNew,E_IsuOld,fen_new_isu,fen_old_isu,E_IsuClustMvConstraint,suBlocks);
        }
        #endif
        #if COUL_ON     //frustrated Ising
        lattice* latInCluster=numToLat(c,inCluster,nInCluster);
        lattice* subset=NULL;
        int lsubset=localSubset(latInCluster,nInCluster,&subset);
        double Eold[ENERGY_LENGTH];
        energy_local(s,Eold,p,subset,lsubset);
        free(latInCluster); free(subset);
        latInCluster=numToLat(tc,inCluster,nInCluster);
        lsubset=localSubset(latInCluster,nInCluster,&subset);
        double Enew[ENERGY_LENGTH];
        energy_local(s,Enew,p,subset,lsubset); //just use s for isimult
        free(latInCluster); free(subset);
        Eold[6]=Enew[6]=0.0;
        #if FFT_ON
        Eold[6]=(p->Q)*energy_lr_nflip(c,ewald,inCluster,nInCluster);
        Enew[6]=(p->Q)*energy_lr_nflip(tc,ewald,inCluster,nInCluster);
        #endif
        //E: 0:ising svsv  1:ising svsu  2:ising susu  3:csr svsv  4:csr svsu  5: csr susu  6: clr ..
        for(int i=0;i<ENERGY_LENGTH;++i) deltaE[i]=Enew[i]-Eold[i];
        deltaE[1]+=E_IsuClustMvConstraint;
        clusterE=(1.0-(p->kT)*(p->beta_eff))*deltaE[0];
        for(int i=1;i<ENERGY_LENGTH;++i) clusterE+=deltaE[i];

        #if TEST_BUILD
        //full E for cmp w local scheme
        double fEold[ENERGY_LENGTH];
        energy_full(c,s,fEold,ewald,p);
        double fEnew[ENERGY_LENGTH];
        energy_full(tc,s,fEnew,ewald,p);
        double fdeltaE[ENERGY_LENGTH];
        for(int i=0;i<ENERGY_LENGTH;++i) fdeltaE[i]=fEnew[i]-fEold[i];
        fdeltaE[1]+=E_IsuClustMvConstraint;
        double fclusterE=(1.0-(p->kT)*(p->beta_eff))*fdeltaE[0];
        for(int i=1;i<ENERGY_LENGTH;++i) fclusterE+=fdeltaE[i];
        if(fabs(fclusterE-clusterE) >= 0.00001) {
                fprintf(myerr,"full clusterE,local,diff: %f,%f,%f\n",fclusterE,clusterE,fclusterE-clusterE);
                fprintf(myerr,"\tfull\tlocal\n0:\t%f\t%f\n1:\t%f\t%f\n2:\t%f\t%f\n3:\t%f\t%f\n"           \
                              "4:\t%f\t%f\n5:\t%f\t%f\n6:\t%f\t%f\n7:\t%f\t%f\n8:\t%f\t%f\n9:\t%f\t%f\n",  \
                        fdeltaE[0],deltaE[0],fdeltaE[1],deltaE[1],fdeltaE[2],deltaE[2],fdeltaE[3],deltaE[3],\
                        fdeltaE[4],deltaE[4],fdeltaE[5],deltaE[5],fdeltaE[6],deltaE[6],fdeltaE[7],deltaE[7], \
                        fdeltaE[8],deltaE[8],fdeltaE[9],deltaE[9]                                             );
        }
        fprintf(myerr,"i\tx,y,z\tval\tlsubset %d\n",lsubset);
        for(int i=0;i<nInCluster;++i) {
                const int num=subset[i]-c;
                int tmppos[d];
                getPos(inCluster[i],tmppos);
                fprintf(myerr,"%d\t%d,%d,%d\t%f\n",inCluster[i],tmppos[0],tmppos[1],tmppos[2],*c[inCluster[i]].val);
        }
        #endif
        #endif
        free(inCluster);

        if(clusterE<=0.0) { //accepted
                for(int i=0;i<N;++i) *c[i].val=*tc[i].val;
                if(E!=NULL) {
                        for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaE[i];
                }
                return (int)nInCluster;
        }

        // check Boltz. prob
        const double pr=(double)pcg32_random_r(rng)/UINT32_MAX; // pr in [0,1)
        const double boltz=exp(-(p->beta) * clusterE);
        //fprintf(myerr,"nInCluster: %d\t(beta-beta_eff)*deltaE_I: %f\tbeta*deltaE_C:%f\tbeta*clusterE: %f\tboltz: %f\tp: %f\n",(nInCluster),(p->beta-p->beta_eff)*deltaE_I,(p->beta)*deltaE_C,(p->beta)*clusterE,boltz,pr);
        if(pr<=boltz) { // accepted
                for(int i=0;i<N;++i) *c[i].val=*tc[i].val;
                if(E!=NULL) {
                        for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaE[i];
                }
                return (int)nInCluster;
        }
        else {
                for(int i=0;i<N;++i) *tc[i].val=*c[i].val;
                return 0;
        }
}

//suBlocks tracks number of submoves blocked because ONLY one
//of the NNs is a su. When this occurs, the submove is auto.
//rejected, but there is a cost incurred. Half of the cost
//is tracked in the solute part of the Ising interaction.
//The other half is the spin that might have been flipped
//had its partner not been in a cluster. The number of these
//spins which are the same as their cluster are tracked by
//adding 1 to suBlocks; those which are opposite their clu-
//ster are tracked by subtracting 1 from suBlocks. The energy
//cost is computed in clusterMove and is:
//+2*J*suBlocks

//built for fully occupied lattice (ie. each
//site it +1 or -1). need to rewrite if want
//to expand to case with defects (0's) or
//more general solvent types (+2, -1, or di-
//fferent shapes, etc etc)
uint32_t* buildCluster(lattice c, const uint32_t ind0, const uint32_t ind1, const uint32_t axis, uint32_t* nInCluster, double* suBlocks, para p, int coreNum, pcg32_random_t* rng) {
        double spin0=*c[ind0].val, spin1=*c[ind1].val;
        extern int d,N;
        extern int *nPlusSvSites,*nMinusSvSites;
        const int stackL = (nPlusSvSites[coreNum]>=nMinusSvSites[coreNum]) ? nPlusSvSites[coreNum] : nMinusSvSites[coreNum];
        uint32_t* inCluster=malloc(N*sizeof(*inCluster));
        uint32_t* s0=malloc(stackL*sizeof(*s0));
        uint32_t* s1=malloc(stackL*sizeof(*s1));
        assert(inCluster!=NULL && s0!=NULL && s1!=NULL /*malloc*/);
        for(int i=0;i<N;++i) inCluster[i]=UINT32_MAX; //init for safety
                
        #if COUL_ON     //for FI
        const double boltz=1.0-exp(-4.0*(p->beta_eff)*(p->J));
        #endif
        #if !COUL_ON    //for pure I
        const double boltz=1.0-exp(-4.0*(p->beta)*(p->J));
        #endif
        *nInCluster=0u;
        s0[0]=ind0;
        s1[0]=ind1;
        inCluster[(*nInCluster)++]=s0[0];
        inCluster[(*nInCluster)++]=s1[0];
        uint32_t i=0u, end=1u;
        while(end!=i) {
                *c[s0[i]].val=spin1;
                *c[s1[i]].val=spin0;
                uint32_t nn0[2*d], nn1[2*d];
                uint32_t lnn=0u;
                for(int e=0;e<2*d;++e) { //d=3 => nn->[-z,-y,-x,x,y,z]
                        nn0[lnn]=(uint32_t)(c[s0[i]].nearNeigh[e]-c);
                        const int erot=(e-d!=axis && e-d!=-axis-1) ? 2*d-e-1 : e; //axis->e 0->2 1->1 2->0; consider rotn by pi abt axis
                        nn1[lnn]=(uint32_t)(c[s1[i]].nearNeigh[erot]-c);
                        if(c[nn0[lnn]].su == -1 && c[nn1[lnn]].su == -1) {
                                if( *c[nn0[lnn]].val == spin0 && *c[nn1[lnn]].val == spin1         \
                                    && bsearch((nn0+lnn),inCluster,*nInCluster,sizeof(uint32_t),    \
                                                (int (*)(const void*,const void*))uint32_tCmp)==NULL \
                                    && bsearch((nn1+lnn),inCluster,*nInCluster,sizeof(uint32_t),      \
                                                 (int (*)(const void*,const void*))uint32_tCmp)==NULL  )   ++lnn;
                        }
                        else { //if only one in su, update suBlocks ctr
                                if(c[nn0[lnn]].su != -1 && c[nn1[lnn]].su != -1) {
                                        continue; //both su: drop out
                                }
                                uint32_t* nnNotSu;
                                double* nonSuClusterSpin;
                                if(c[nn0[lnn]].su != -1) {
                                        nnNotSu=nn1;
                                        nonSuClusterSpin=&spin1;
                                }
                                else {
                                        nnNotSu=nn0;
                                        nonSuClusterSpin=&spin0;
                                }
                                if(*c[nnNotSu[lnn]].val == *nonSuClusterSpin) *suBlocks+=1.0;
                                else                                          *suBlocks-=1.0;
                        }
                }
                for(uint32_t j=0;j<lnn;++j) {
                        const double pr=(double)pcg32_random_r(rng)/UINT32_MAX; // pr in [0,1)
                        if(pr <= boltz) { //add to cluster
                                s0[end]=nn0[j];
                                s1[end]=nn1[j];
                                inCluster[(*nInCluster)++]=s0[end];
                                inCluster[(*nInCluster)++]=s1[end];
                                ++end;
                        }
                }
                ++i;
                qsort(inCluster,*nInCluster,sizeof(uint32_t),(int (*)(const void*,const void*))uint32_tCmp);
        }
        free(s0); free(s1);
        inCluster=realloc(inCluster,*nInCluster*sizeof(*inCluster));
        assert(inCluster!=NULL /*realloc*/);
        return inCluster;
}

//if su is asymmetric and rotates, may need nudge back onto lattice
double* comOrientationNudge(lattice c, su s, int ind, double* tCom, double** tOrientation, pcg32_random_t* rng, FILE* myerr) {
        extern int d;
        const int nsites=s[ind].nsites, ninshl=s[ind].ninshl;
        double** absPos=malloc((nsites+ninshl)*sizeof(*absPos));
        double* absPosH=malloc((nsites+ninshl)*d*sizeof(*absPosH));
        assert(absPos!=NULL && absPosH!=NULL);
        for(int i=0;i<nsites;++i) {
                absPos[i]=absPosH+i*d;
                for(int e=0;e<d;++e) absPos[i][e]=s[ind].relPos[i][e];
                matrixOnVectorInPlace(tOrientation,absPos[i]);
                for(int e=0;e<d;++e) absPos[i][e]+=tCom[e];
        }
        for(int i=0;i<ninshl;++i) {
                absPos[i+nsites]=absPosH+(i+nsites)*d;
                for(int e=0;e<d;++e) absPos[i+nsites][e]=s[ind].inShlRelPos[i][e];
                matrixOnVectorInPlace(tOrientation,absPos[i+nsites]);
                for(int e=0;e<d;++e) absPos[i+nsites][e]+=tCom[e];
        }
        double* delta=malloc(d*sizeof(*delta));
        assert(delta!=NULL);
        for(int e=0;e<d;++e) delta[e]=0.0;
        bool needsNudge=false;
        for(int i=0;i<nsites;++i) {
                double r[d];
                for(int e=0;e<d;++e) {
                        r[e]=fabs(fmod(absPos[i][e],1.0));
                        if(r[e]!=0.0) {
                                needsNudge=true;
                                if(delta[e]!=0.0) {
                                        if(fabs(r[e]-delta[e])>0.00001) {
                                                fprintf(myerr,"comOrientationReqNudge: unexpected case, delta[e]!=r[e]\ne,delta,r %d,%.2f,%.2f\n",e,delta[e],r[e]);
                                                exit(1);
                                        }
                                }
                                else delta[e]=r[e];
                        }
                }
        }
        free(absPosH); free(absPos);
        if(needsNudge==true) {
                for(int e=0;e<d;++e) {
                        uint32_t rn=pcg32_boundedrand_r(rng,2); //0 -> negative, 1 -> positive
                        if(rn==0) delta[e]*=-1.0;
                        else      continue;
                }
                return delta;
        }
        else {
                free(delta);
                return NULL;
        }
}

//NOTE: contrary to usual convention,
//1 => move sucessful; 0 => move rejected 
//
//suMove (soluteMove) tuned for use in Frustrated Ising
//system. Observation: charged solutes will have
//tendency to be surrounded by, first, a nearly
//non-fluctuating few layers of solvent, and next,
//by a further few fluctuating layers of solvent.
//Thus, need to make certain that trial moves have
//config.s satisfying these observations or energy
//of trial moves will be high & acceptance ratio
//low. The move will consist of moving the solute
//and an adjacent few layers of solvent ('inner
//shell')
//
//Then the energy shall be calculated, etc.
//
//return 0 -> move rejected
//return 1 -> move accepted
//return -1 -> some function call failed
//Designed to exit gracefully on return 0 or 1,
//merely to spit out error message and return -1
//on a fail (no cleanup, etc.)
//
//sweepsMult and kTeff, if both non-zero, determine
//the behavior of solvent relaxation built into
//the solute move.
//
//the solvent relaxation is supposed to help increase
//the probability of successful solute moves. this
//only happens when kTeff > kT. At low T, solutes have
//trouble moving because they disrupt ordered oscill-
//ations of solvent, but by "locally heating" the
//solvent, we can, for the right parameters, increase
//the probability of a successful move.

int suMove(lattice c, lattice tc, su s, su ts, umb u, umb tu, double* E, int ind, double* ewald, int inshlopt, bool rotation, double sweepsMult, double kTeff, para p, int coreNum, pcg32_random_t* rng, FILE* myerr) {
        extern int d,N;
        void* ptrs[10]; //hardcode
        int nptrs=0;
        //roll dist., orientation, dir
        double trialCom[d];
        double oldUnwrappedCom[d];
        double unwrappedCom[d];
        for(int e=0;e<d;++e) {
                trialCom[e]=s[ind].com[e];
                unwrappedCom[e]=s[ind].unwrappedCom[e];
                oldUnwrappedCom[e]=s[ind].unwrappedCom[e];
        }
        const int changeD=s[ind].nMotion;
        const int strtD=d-changeD;
        int j=0;
        for(int e=strtD;e<d;++e) { //dist in suStep*[-0.5,0.5)
                const double rn=(double)pcg32_random_r(rng)/UINT32_MAX-0.5;
                unwrappedCom[e]+=round(s[ind].suStep*rn);
                trialCom[e]=unwrappedCom[e];
                wrapIntoL_double(trialCom, NULL); //PBC
                if(fabs(trialCom[e]-s[ind].com[e])<=DBL_EPSILON) ++j;
        }

        double** oldOrientation=copyMatrix(s[ind].orientation,d);
        double** trialOrientation=copyMatrix(s[ind].orientation,d);
        ptrs[nptrs++]=(void*)*oldOrientation;
        ptrs[nptrs++]=(void*)oldOrientation;
        ptrs[nptrs++]=(void*)*trialOrientation;
        ptrs[nptrs++]=(void*)trialOrientation;
        double oldCom[d];
        for(int e=0;e<d;++e) oldCom[e]=s[ind].com[e];
        double*** R=NULL;
        if(rotation==true) {
                R=buildPiO2RotationMatricies(d);
                int* eOrd=forRandOrder(d,rng);
                for(int e=0;e<d;++e) {
                        const int ee=eOrd[e];
                        const uint32_t rn=pcg32_boundedrand_r(rng,4); //hit each dir. w/ rotation [0,3pi/2]
                        for(uint32_t i=0;i<rn;++i) matrixOnMatrixInPlace(R[ee],trialOrientation);
                }
                free(eOrd);
                copyMatrixInPlace(trialOrientation,ts[ind].orientation,d);
                ptrs[nptrs++]=(void*)**R;
                ptrs[nptrs++]=(void*)*R;
                ptrs[nptrs++]=(void*)R;
                if(j==changeD && doubleMatrixCmp(trialOrientation,oldOrientation,d)==1) {
                        freeAll(ptrs,nptrs); //proposed move is: no move
                        return 0;
                }
        }
        else if(j==changeD) { //proposed move is: no move
                freeAll(ptrs,nptrs);
                return 0;
        }
        //if su is asymmetric and rotates, may need nudge back onto lattice
        double* nudge=comOrientationNudge(c,s,ind,trialCom,trialOrientation,rng,myerr);
        if(nudge!=NULL) {
                #if TEST_BUILD
                fprintf(myerr,"nudged tCom %.2f,%.2f,%.2f by %.2f,%.2f,%.2f\n",trialCom[0],trialCom[1],trialCom[2],nudge[0],nudge[1],nudge[2]);
                #endif
                for(int e=0;e<d;++e) {
                        unwrappedCom[e]+=nudge[e];
                        trialCom[e]=unwrappedCom[e];
                        wrapIntoL_double(trialCom, NULL); //PBC
                }
                free(nudge);
        }
        for(int e=0;e<d;++e) {
                ts[ind].unwrappedCom[e]=unwrappedCom[e];
                ts[ind].com[e]=trialCom[e];
        }

        //setup ts & tc for proposed move
        updSuCurPos(tc,ts,ind,(double*)NULL);
        extern bool hydrophobicExist;
        if(checkSuOverlap(ts,ind,inshlopt)==1) {
                for(int k=0;k<N;++k) *tc[k].val=*c[k].val;
                for(int e=0;e<d;++e) {
                        ts[ind].unwrappedCom[e]=oldUnwrappedCom[e];
                        ts[ind].com[e]=oldCom[e];
                }
                copyMatrixInPlace(oldOrientation,ts[ind].orientation,d);
                updSuCurPos_latSuStatus(tc,ts,ind); //waste?
                if(hydrophobicExist==true) updLatHydrophobicStatus(tc,ts);
                freeAll(ptrs,nptrs);
                #if TEST_BUILD
                fprintf(myerr,"suMove: overlap w another solute\n");
                #endif
                return 0;
        }

        //build local region for suMove. used for local E calc as
        //well as for sv relaxtion pre&post su move if desired
        uint32_t* core=NULL;
        uint32_t lcore=buildSuMoveCore(c,tc,s,ts,ind,&core);
        uint32_t lcoreE=lcore;
        lattice* coreE=NULL;
        double EsvRelaxi0[ENERGY_LENGTH];
        double EsvRelaxi1[ENERGY_LENGTH];
        #if TEST_BUILD
        double EsvRelaxi0Full[ENERGY_LENGTH];
        double EsvRelaxi1Full[ENERGY_LENGTH];
        #endif
        for(int i=0;i<ENERGY_LENGTH;++i) {
                EsvRelaxi0[i]=0.0;
                EsvRelaxi1[i]=0.0;
                #if TEST_BUILD
                EsvRelaxi0Full[i]=0.0;
                EsvRelaxi1Full[i]=0.0;
                #endif
        }
        if(sweepsMult!=0.0 && kTeff!=0.0) { //relax solvent local to solute
                lattice* coreLat=numToLat(tc,core,(int)lcore);
                lattice* suLocalRegionLat=NULL;
                lcore=(uint32_t)localSubset(coreLat,(int)lcore,&suLocalRegionLat);
                free(core);
                core=latToNum(tc,suLocalRegionLat,(int)lcore);
                //setup for E computation
                uint32_t* corePlusNNNum=addNNToCore_toNewArr(tc,ts,core,&lcoreE);
                lattice* corePlusNNLat=numToLat(tc,core,lcoreE);
                lcoreE=localSubset(corePlusNNLat,lcoreE,&coreE);
                energy_local(s,EsvRelaxi0,p,coreE,lcoreE); //E_i
                EsvRelaxi0[6]=(p->Q)*energy_lr_nflip(tc,ewald,core,lcore);
                #if TEST_BUILD
                energy_full(tc,ts,EsvRelaxi0Full,ewald,p);
                #endif
                const int nSucMoves=relaxSvNearSu(tc,ts,core,lcore,ewald,sweepsMult,kTeff,p,coreNum,rng,myerr);
                energy_local(s,EsvRelaxi1,p,coreE,lcoreE); //E_f
                #if TEST_BUILD
                energy_full(tc,ts,EsvRelaxi1Full,ewald,p);
                #endif
                EsvRelaxi1[6]=(p->Q)*energy_lr_nflip(tc,ewald,core,lcore);
                free(coreLat);
                free(suLocalRegionLat);
                free(corePlusNNNum);
                free(corePlusNNLat);
        }
        else {
                lattice* coreLat=numToLat(tc,core,(int)lcore);
                lcoreE=localSubset(coreLat,(int)lcore,&coreE);
                free(coreLat);
        }
        ptrs[nptrs++]=(void*)core;
        ptrs[nptrs++]=(void*)coreE;
        double Esui[ENERGY_LENGTH];
        double EsuiSum=energy_local(s,Esui,p,coreE,lcoreE);
        Esui[6]=(p->Q)*energy_lr_nflip(tc,ewald,core,lcore);
        EsuiSum+=Esui[6];
        #if TEST_BUILD
        double EsuiFull[ENERGY_LENGTH];
        double EsuiSumFull=energy_full(tc,ts,EsuiFull,ewald,p);
        #endif

        //trial lattice: move the su
        updLatSuStatus(tc,ts);
        if(hydrophobicExist==true) updLatHydrophobicStatus(tc,ts);

        //find proposed vacacted sites using ts on c
        uint32_t* vacated=NULL;
        const int lv=findVacated(c,tc,s,ts,ind,&vacated);
        ptrs[nptrs++]=(void*)vacated;
        if(lv==0) { //happens when rotates in place
                #if TEST_BUILD
                fprintf(myerr,"suMove: findVacated returns 0 sites vacated\n");
                #endif
                for(int k=0;k<N;++k) *tc[k].val=*c[k].val;
                for(int e=0;e<d;++e) {
                        ts[ind].unwrappedCom[e]=oldUnwrappedCom[e];
                        ts[ind].com[e]=oldCom[e];
                }
                copyMatrixInPlace(oldOrientation,ts[ind].orientation,d);
                updSuCurPos_latSuStatus(tc,ts,ind);
                if(hydrophobicExist==true) updLatHydrophobicStatus(tc,ts);
                freeAll(ptrs,nptrs);
                return 0;
        }

        const int nsites=s[ind].nsites;
        const int ninshl=s[ind].ninshl;
        for(int i=0;i<nsites;++i) *ts[ind].currPos[i]->val=ts[ind].surfType;
        for(int i=0;i<ninshl;++i) {
                if(ts[ind].inShlCurrPos[i]->su==-1 && s[ind].inShlCurrPos[i]->su==-1) {
                        *ts[ind].inShlCurrPos[i]->val=*s[ind].inShlCurrPos[i]->val;
                }
        }

        double* tcchg=checkCharge(tc,N,NULL);
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites;
        int desiredPlus=nPlusSvSites[coreNum]-tcchg[0],desiredMinus=nMinusSvSites[coreNum]-tcchg[1],desiredZero=nZeroSvSites[coreNum]-tcchg[2];
        free(tcchg);

        int neutRet=maintainChargeNeut(tc,vacated,lv,0,&desiredPlus,&desiredMinus,&desiredZero,coreNum,rng,myerr); //giving lc=0 -> skips internal checkCharge()
        if(neutRet!=0) { //neut failed
                for(int k=0;k<N;++k) *tc[k].val=*c[k].val;
                for(int e=0;e<d;++e) {
                        ts[ind].unwrappedCom[e]=oldUnwrappedCom[e];
                        ts[ind].com[e]=oldCom[e];
                }
                copyMatrixInPlace(oldOrientation,ts[ind].orientation,d);
                updSuCurPos_latSuStatus(tc,ts,ind);
                if(hydrophobicExist==true) updLatHydrophobicStatus(tc,ts);
                freeAll(ptrs,nptrs);
                #if TEST_BUILD
                fprintf(myerr,"suMove: neut failed; move refused\n");
                #endif
                return 0;
        }

        //get energy post su move
        double Esuf[ENERGY_LENGTH];
        double EsufSum=energy_local(s,Esuf,p,coreE,lcoreE);
        Esuf[6]=(p->Q)*energy_lr_nflip(tc,ewald,core,lcore);
        EsufSum+=Esuf[6];
        #if TEST_BUILD
        double EsufFull[ENERGY_LENGTH];
        double EsufSumFull=energy_full(tc,ts,EsufFull,ewald,p);
        #endif

        //post-su relaxation moves if desired. already built core
        double EsvRelaxf0[ENERGY_LENGTH];
        double EsvRelaxf1[ENERGY_LENGTH];
        #if TEST_BUILD
        double EsvRelaxf0Full[ENERGY_LENGTH];
        double EsvRelaxf1Full[ENERGY_LENGTH];
        #endif
        for(int i=0;i<ENERGY_LENGTH;++i) {
                EsvRelaxf0[i]=0.0;
                EsvRelaxf1[i]=0.0;
                #if TEST_BUILD
                EsvRelaxf0Full[i]=0.0;
                EsvRelaxf1Full[i]=0.0;
                #endif
        }
        if(sweepsMult!=0.0 && kTeff!=0.0) { //relax solvent local to solute
                energy_local(s,EsvRelaxf0,p,coreE,lcoreE); //E_i
                EsvRelaxf0[6]=(p->Q)*energy_lr_nflip(tc,ewald,core,lcore);
                #if TEST_BUILD
                energy_full(tc,ts,EsvRelaxf0Full,ewald,p);
                #endif
                const int nSucMoves=relaxSvNearSu(tc,ts,core,lcore,ewald,sweepsMult,kTeff,p,coreNum,rng,myerr);
                energy_local(s,EsvRelaxf1,p,coreE,lcoreE); //E_f
                EsvRelaxf1[6]=(p->Q)*energy_lr_nflip(tc,ewald,core,lcore);
                #if TEST_BUILD
                energy_full(tc,ts,EsvRelaxf1Full,ewald,p);
                #endif
        }

        const double deltaEumb=energy_umbr(tu)-energy_umbr(u);
        const double b=p->beta;
        const double be=1.0/(p->kT+kTeff);
        double betadeltaen=b*deltaEumb;
        double betadeltaE[ENERGY_LENGTH];
        double deltaE[ENERGY_LENGTH];
        double deltaESum=b*deltaEumb;
        for(int i=0;i<ENERGY_LENGTH;++i) {
                if(sweepsMult!=0.0 && kTeff!=0.0) {
                        deltaE[i]=EsvRelaxf1[i]-EsvRelaxi0[i];
                }
                else {
                        deltaE[i]=Esuf[i]-Esui[i];
                }
                deltaESum+=deltaE[i];
                betadeltaE[i]=b*(Esuf[i]-Esui[i])+ \
                              (b-be)*(EsvRelaxf1[i]-EsvRelaxf0[i]+EsvRelaxi1[i]-EsvRelaxi0[i]);
                betadeltaen+=betadeltaE[i];
        }

        //cmp to full energy:
        #if TEST_BUILD
        double betadeltaenFull=b*deltaEumb;
        double betadeltaEFull[ENERGY_LENGTH];
        double deltaEFull[ENERGY_LENGTH];
        double deltaEFullSum=b*deltaEumb;
        for(int i=0;i<ENERGY_LENGTH;++i) {
                if(sweepsMult!=0.0 && kTeff!=0.0) {
                        deltaEFull[i]=EsvRelaxf1Full[i]-EsvRelaxi0Full[i];
                }
                else {
                        deltaEFull[i]=EsufFull[i]-EsuiFull[i];
                }
                deltaEFullSum+=deltaEFull[i];
                betadeltaEFull[i]=b*(EsufFull[i]-EsuiFull[i])+ \
                              (b-be)*(EsvRelaxf1Full[i]-EsvRelaxf0Full[i]+EsvRelaxi1Full[i]-EsvRelaxi0Full[i]);
                betadeltaenFull+=betadeltaEFull[i];
        }
        if(fabs(betadeltaenFull-betadeltaen) >= 0.00001 || fabs(deltaEFullSum-deltaESum) >= 0.00001) {
                fprintf(myerr,"suMove(): betadeltaenFull,betadeltaen,diff: %f,%f,%f\n",betadeltaenFull,betadeltaen,betadeltaenFull-betadeltaen);
                fprintf(myerr,"deltaEFullSum,deltaESum,diff: %f,%f,%f\n",deltaEFullSum,deltaESum,deltaEFullSum-deltaESum);
                fprintf(myerr,"\tfull\tlocal\tf-l\n0:\t%f\t%f\t%f\n1:\t%f\t%f\t%f\n2:\t%f\t%f\t%f\n3:\t%f\t%f\t%f\n4:\t%f\t%f\t%f\n" \
                              "5:\t%f\t%f\t%f\n6:\t%f\t%f\t%f\n7:\t%f\t%f\t%f\n8:\t%f\t%f\t%f\n9:\t%f\t%f\t%f\n", \
                        betadeltaEFull[0],betadeltaE[0],betadeltaEFull[0]-betadeltaE[0],betadeltaEFull[1],betadeltaE[1],betadeltaEFull[1]-betadeltaE[1], \
                        betadeltaEFull[2],betadeltaE[2],betadeltaEFull[2]-betadeltaE[2],betadeltaEFull[3],betadeltaE[3],betadeltaEFull[3]-betadeltaE[3], \
                        betadeltaEFull[4],betadeltaE[4],betadeltaEFull[4]-betadeltaE[4],betadeltaEFull[5],betadeltaE[5],betadeltaEFull[5]-betadeltaE[5], \
                        betadeltaEFull[6],betadeltaE[6],betadeltaEFull[6]-betadeltaE[6],betadeltaEFull[7],betadeltaE[7],betadeltaEFull[7]-betadeltaE[7], \
                        betadeltaEFull[8],betadeltaE[8],betadeltaEFull[8]-betadeltaE[8],betadeltaEFull[9],betadeltaE[9],betadeltaEFull[9]-betadeltaE[9]);
                fprintf(myerr,"corei\tval,su\ttval,tsu\ttval-val\n");
                for(int i=0;i<lcore;++i) {
                        int ci=(int)core[i];
                        int tmppos[d];
                        getPos(ci,tmppos);
                        fprintf(myerr,"%d\t%d,%d,%d\t%.1f,%d\t%.1f,%d\t%.1f\n",ci,tmppos[0],tmppos[1],tmppos[2],*c[ci].val,c[ci].su,*tc[ci].val,tc[ci].su,*tc[ci].val-*c[ci].val);
                }
                //assert(lsubset_old==lsubset_new);
                //fprintf(myerr,"subsi\tval,su\ttval,tsu\ttval-val\n");
                //for(int i=0;i<lsubset_old;++i) {
                //        if(fabs(*subset_new[i]->val-*subset_old[i]->val)>0.00001 || subset_new[i]->su!=subset_old[i]->su) {
                //                fprintf(myerr,"%d\t%.1f,%d\t%.1f,%d\t%.1f\n",i,*subset_old[i]->val,subset_old[i]->su,*subset_new[i]->val,subset_new[i]->su,*subset_new[i]->val-*subset_old[i]->val);
                //        }
                //}
                fprintf(myerr,"latti\tval,su\ttval,tsu\ttval-val\n");
                for(int i=0;i<N;++i) {
                        if(fabs(*tc[i].val-*c[i].val)>0.00001 || tc[i].su!=c[i].su) {
                                int tmppos[d];
                                getPos(i,tmppos);
                                fprintf(myerr,"%d\t%d,%d,%d\t%.1f,%d\t%.1f,%d\t%.1f\n",i,tmppos[0],tmppos[1],tmppos[2],*c[i].val,c[i].su,*tc[i].val,tc[i].su,*tc[i].val-*c[i].val);
                        }
                }
                for(int k=0;k<N;++k) *tc[k].val=*c[k].val;
                for(int e=0;e<d;++e) {
                        ts[ind].unwrappedCom[e]=oldUnwrappedCom[e];
                        ts[ind].com[e]=oldCom[e];
                }
                copyMatrixInPlace(oldOrientation,ts[ind].orientation,d);
                updSuCurPos_latSuStatus(tc,ts,ind);
                if(hydrophobicExist==true) updLatHydrophobicStatus(tc,ts);
                freeAll(ptrs,nptrs);
                return 0;
        }
        #endif

        //accept/reject config?
        if(betadeltaen<=0.0) { // accepted
                for(int k=0;k<N;++k) *c[k].val=*tc[k].val;
                for(int e=0;e<d;++e) {
                        s[ind].unwrappedCom[e]=unwrappedCom[e];
                        s[ind].com[e]=trialCom[e];
                }
                copyMatrixInPlace(trialOrientation,s[ind].orientation,d);
                updSuCurPos_latSuStatus(c,s,ind);
                if(hydrophobicExist==true) updLatHydrophobicStatus(c,s);
                for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaE[i];
                freeAll(ptrs,nptrs);
                return 1;
        }

        //check Boltz. prob;
        const double pr=(double)pcg32_random_r(rng)/UINT32_MAX; //pr in [0,1)
        const double boltz=exp(-betadeltaen);
        if(pr<=boltz) { //accepted
                for(int k=0;k<N;++k) *c[k].val=*tc[k].val;
                for(int e=0;e<d;++e) {
                        s[ind].unwrappedCom[e]=unwrappedCom[e];
                        s[ind].com[e]=trialCom[e];
                }
                copyMatrixInPlace(trialOrientation,s[ind].orientation,d);
                updSuCurPos_latSuStatus(c,s,ind);
                if(hydrophobicExist==true) updLatHydrophobicStatus(c,s);
                for(int i=0;i<ENERGY_LENGTH;++i) E[i]+=deltaE[i];
                freeAll(ptrs,nptrs);
                return 1;
        }
        else { //rejected
                for(int k=0;k<N;++k) *tc[k].val=*c[k].val;
                for(int e=0;e<d;++e) {
                        ts[ind].unwrappedCom[e]=oldUnwrappedCom[e];
                        ts[ind].com[e]=oldCom[e];
                }
                copyMatrixInPlace(oldOrientation,ts[ind].orientation,d);
                updSuCurPos_latSuStatus(tc,ts,ind);
                if(hydrophobicExist==true) updLatHydrophobicStatus(tc,ts);
                freeAll(ptrs,nptrs);
                return 0;
        }
}

//union
uint32_t buildSuMoveCore(lattice c, lattice tc, su s, su ts, int ind, uint32_t** core) {
        const int nsites=s[ind].nsites, ninshl=s[ind].ninshl;
        int lcore=2*(nsites+ninshl);
        (*core)=malloc(lcore*sizeof(**core));
        assert((*core)!=NULL /*malloc*/);
        for(int i=0;i<nsites;++i) {
                (*core)[i]=s[ind].currPos[i]-c;
                (*core)[i+lcore/2]=ts[ind].currPos[i]-tc;
        }
        for(int i=0;i<ninshl;++i) {
                (*core)[i+nsites]=s[ind].inShlCurrPos[i]-c;
                (*core)[i+nsites+lcore/2]=ts[ind].inShlCurrPos[i]-tc;
        }
        qsort((*core),lcore,sizeof(uint32_t),(int (*)(const void*,const void*))uint32_tCmp);
        lcore=uniqueuint32_t(core,lcore);
        (*core)=realloc((*core),lcore*sizeof(**core));
        assert((*core)!=NULL /*realloc*/);
        return lcore;
}

uint32_t addNNToCore(lattice c, su s, uint32_t** core, uint32_t lcore) {
        extern int d;
        (*core)=realloc((*core),((2*d+1)*lcore)*sizeof(**core));
        assert((*core)!=NULL /*realloc*/);
        int j=lcore;
        for(int i=0;i<lcore;++i) {
                for(int e=0;e<2*d;++e) {
                        (*core)[j]=(uint32_t)(c[(*core)[i]].nearNeigh[e]-c);
                        ++j;
                }
        }
        qsort((*core),(2*d+1)*lcore,sizeof(uint32_t),(int (*)(const void*,const void*))uint32_tCmp);
        lcore=uniqueuint32_t(core,(2*d+1)*lcore);
        (*core)=realloc((*core),lcore*sizeof(**core));
        assert((*core)!=NULL /*realloc*/);
        return lcore;
}

uint32_t* addNNToCore_toNewArr(lattice c, su s, uint32_t* core, uint32_t* lcore) {
        uint32_t* coreNew=malloc((*lcore)*sizeof(*coreNew));
        for(int i=0;i<*lcore;++i) coreNew[i]=core[i];
        *lcore=addNNToCore(c,s,&coreNew,*lcore);
        return coreNew;
}

//intersection
int findVacated(lattice c, lattice tc, su s, su ts, int ind, uint32_t** v) {
        const int nsites=s[ind].nsites, ninshl=s[ind].ninshl;
        const int lmax=nsites+ninshl;
        (*v)=malloc(lmax*sizeof(**v));
        uint32_t* tspos=malloc(lmax*sizeof(*tspos));
        assert((*v)!=NULL && tspos!=NULL /*malloc*/);
        for(int i=0;i<nsites;++i) tspos[i]=ts[ind].currPos[i]-tc;
        for(int i=0;i<ninshl;++i) tspos[i+nsites]=ts[ind].inShlCurrPos[i]-tc;
                
        qsort(tspos,lmax,sizeof(uint32_t),(int (*)(const void*,const void*))uint32_tCmp);

        int lv=0;
        for(int i=0;i<nsites;++i) {
                const uint32_t site=s[ind].currPos[i]-c;
                void* bsearchRet=bsearch(&site,tspos,lmax,sizeof(uint32_t),
                                         (int (*)(const void*,const void*))uint32_tCmp);
                if(bsearchRet != NULL) continue; //if found, not vacated
                (*v)[lv++]=site;
        }
        for(int i=0;i<ninshl;++i) {
                if(s[ind].inShlCurrPos[i]->su!=-1) continue; //another su occupies inShl site
                const uint32_t site=s[ind].inShlCurrPos[i]-c;
                void* bsearchRet=bsearch(&site,tspos,lmax,sizeof(uint32_t),
                                         (int (*)(const void*,const void*))uint32_tCmp);
                if(bsearchRet != NULL) continue; //if found, not vacated
                (*v)[lv++]=site;
        }
        free(tspos);
        (*v)=realloc((*v),lv*sizeof(**v));
        if(lv>0) assert(*v!=NULL /*realloc*/);
        return lv;
}

int relaxSvNearSu(lattice tc, su ts, uint32_t* core, uint32_t lcore, double* ewald, double sweepsMult, double kTeff, para p, int coreNum, pcg32_random_t* rng, FILE* myerr) {
        const uint32_t itera=(uint32_t)round(lcore*sweepsMult);
        para pEff=copyPara(p);
        pEff->kT+=kTeff;
        pEff->beta=1.0/pEff->kT;
        int nSucMoves=0;
        for(uint32_t i=0;i<itera;++i) {
                const uint32_t icore=pcg32_boundedrand_r(rng,lcore);
                const int ind=(int)core[icore];
                if(tc[ind].su==-1) { //not in solute: standard swap/flip move
                        const int tmp=spinSwapMoveNN(tc,NULL,ts,NULL,ind,ewald,pEff,rng,myerr);
                        nSucMoves+=tmp;
                }
                //do nothing if in solute
        }
        free(pEff);
        extern int N;
        double* chg=checkCharge(tc,N,NULL);
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites;
        const int delPlus=chg[0]-nPlusSvSites[coreNum];
        const int delMinus=chg[1]-nMinusSvSites[coreNum];
        const int delZero=chg[2]-nZeroSvSites[coreNum];
        if(delPlus!=0 || delMinus!=0 || delZero!=0) {
                fprintf(myerr,"relaxSvNearSu: delPlus:%d delMinus:%d delZero:%d ; move refused\n",delPlus,delMinus,delZero);
        }
        free(chg);
        return nSucMoves;
}

#if EQUI_INNER_DUMP
int* MCstep_flipswap(lattice c, lattice tc, double* Ein, su s, su ts, umb u, umb tu, double* ewald, int inshlopt, bool rotation, double sweepsMult, double kTeff, para p, int coreNum, pcg32_random_t* rng, int enFreq, int dumpFreq, int suComFreq, int umbrFreq, int dataFreq, FILE* dumpout, char* prefix, int outerN, bool outputStyleFI, FILE* myout, FILE* myumbout, FILE* mySuComOut, FILE* myerr) {
#endif
#if !EQUI_INNER_DUMP
int* MCstep_flipswap(lattice c, lattice tc, double* Ein, su s, su ts, umb u, umb tu, double* ewald, int inshlopt, bool rotation, double sweepsMult, double kTeff, para p, int coreNum, pcg32_random_t* rng, FILE* myout, FILE* myumbout, FILE* mySuComOut, FILE* myerr) {
#endif
        #if EQUI_INNER_DUMP
        extern int d;
        #endif
        extern int N,Nsu;
        const int mvtypes=2*3;  //hardcode
        int* nmoves=malloc(mvtypes*sizeof(*nmoves)); //={sucSvMvs,svMvs,sucSuMvs,suMvs,sucClustMvs,clustMvs}
        assert(nmoves!=NULL /*malloc*/);
        for(int i=0;i<mvtypes;++i)      nmoves[i]=0;    // hardcode

        double* E=NULL;
        #if COUL_ON
        if(Ein==NULL) {
                E=malloc(ENERGY_LENGTH*sizeof(*E));
                assert(E!=NULL /*malloc*/);
                energy_full(c,s,E,ewald,p);
        }
        else    E=Ein;
        #endif

        for(int i=0;i<N;++i) { //N proposed moves 'one step'; ~1MD timestep
                const uint32_t ind=pcg32_boundedrand_r(rng,N);
                if(c[ind].su==-1) { //not in solute: standard swap/flip move
                        const int mvSuccess=spinSwapMove(c,tc,s,E,ind,ewald,p,rng,myerr);
                        if(mvSuccess>0) ++nmoves[0];
                        ++nmoves[1];
                }
                else { //in solute: solute move
                        const double roll=(double)pcg32_random_r(rng)/UINT32_MAX;
                        if(roll<=s[c[ind].su].mvProb) {
                                #if TEST_BUILD
                                double Etmp[ENERGY_LENGTH];
                                energy_full(c,s,Etmp,ewald,p);
                                for(int ii=0;ii<ENERGY_LENGTH;++ii) {
                                        if(fabs(Etmp[ii]-E[ii])>=0.00001) {
                                                fprintf(myerr,"i,Etmp!=E %d,%f!=%f\n",ii,Etmp[ii],E[ii]);
                                        }
                                }
                                #endif
                                const int mvSuccess=suMove(c,tc,s,ts,u,tu,E,c[ind].su,ewald,inshlopt,rotation,sweepsMult,kTeff,p,coreNum,rng,myerr);
                                if(mvSuccess>0) ++nmoves[2];
                                ++nmoves[3];
                        }
                }

                #if EQUI_INNER_DUMP
                const int j=outerN*N+i;
                if(j%enFreq==0) {
                        #if !COUL_ON
                        double Eholder[ENERGY_LENGTH];
                        E=Eholder;
                        #endif
                        #if TEST_BUILD
                        double Etmp[ENERGY_LENGTH];
                        energy_full(c,s,Etmp,ewald,p);
                        for(int ii=0;ii<ENERGY_LENGTH;++ii) {
                                if(fabs(Etmp[ii]-E[ii])>=0.00001) {
                                        fprintf(myerr,"i,Etmp!=E %d,%f!=%f\n",ii,Etmp[ii],E[ii]);
                                }
                        }
                        #endif
                        const double eUmb=energy_umbr(u);
                        const double eSelf=-p->e_self;
                        fprintf(myout,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",j,E[0]+E[1]+E[2]+E[3]+E[4]+E[5]+E[6]+E[7]+E[8]+E[9]+eUmb+eSelf,E[0],E[1],E[2],E[3],E[4],E[5],E[6],E[7],E[8],E[9],eUmb,eSelf);
                        #if !COUL_ON
                        E=NULL;
                        #endif
                }
                if(dumpFreq!=0 && j%dumpFreq==0) {
                        if(outputStyleFI==true) dumpFI(c,s,j,dumpout,myerr);
                        else                    dumpLammps(c,s,j,coreNum,dumpout,myerr);
                }
                if(s!=NULL && suComFreq!=0 && j%suComFreq==0) dumpAllSuComs(s,mySuComOut);
                if(u!=NULL && umbrFreq!=0 && j%umbrFreq==0) {
                        fprintf(myumbout,"%d\t",j);
                        for(int k=0;k<(u->npairs);++k) {
                                const double r_scal=dist_vect((u->pair[k][0])->com,(u->pair[k][1])->com);
                                fprintf(myumbout,"%f\t",r_scal);
                        }
                        fprintf(myumbout,"\n");
                }
                if(dataFreq!=0 && j%dataFreq==0) {
                        const uint16_t lprefix=(uint16_t)strlen(prefix); //strlen ret pos of '\0'
                        const uint16_t ldotdata=(uint16_t)strlen(".data");
                        const uint16_t ndig=ndigits((long int)i);
                        char* datafstr=(char*)malloc(lprefix+1u+ndig+ldotdata);
                        assert(datafstr!=NULL /*malloc*/);
                        strncpy(datafstr,prefix,lprefix);
                        const int snpfret=snprintf(datafstr+lprefix,ndig+2u,"_%u",i);
                        assert(snpfret==ndig+1 /*snprintf*/);
                        strncat(datafstr,".data",ldotdata);
                        FILE* dataout=fopen(datafstr,"w");
                        assert(dataout!=NULL /*fopen*/);
                        dumpData(c,s,dataout);
                        assert(fclose(dataout)==0);
                        free(datafstr); datafstr=NULL;
                }
                #endif
        }
        double* chg=checkCharge(c,N,NULL);
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites;
        const int delPlus=chg[0]-nPlusSvSites[coreNum];
        const int delMinus=chg[1]-nMinusSvSites[coreNum];
        const int delZero=chg[2]-nZeroSvSites[coreNum];
        if(delPlus!=0 || delMinus!=0 || delZero!=0) {
                fprintf(myerr,"MCstep_flipswap: delPlus:%d delMinus:%d delZero:%d ; move refused\n",delPlus,delMinus,delZero);
        }
        free(chg);
        return nmoves;
}

int* MCstep_cluster(lattice c, lattice tc, double* Ein, su s, su ts, umb u, umb tu, double* ewald, int inshlopt, bool rotation, double sweepsMult, double kTeff, para p, int coreNum, pcg32_random_t* rng, FILE* mvstats, FILE* myerr) {
        extern int N,Nsu;
        const int mvtypes=2*3;  //hardcode
        const int cStepsPerSweep=100;     //hardcode
        int* nmoves=malloc(mvtypes*sizeof(*nmoves)); //={sucSolvMvs,svMvs,sucSuMvs,suMvs,sucClustMvs,clustMvs}
        assert(nmoves!=NULL /*malloc*/);
        for(int i=0;i<mvtypes;++i) nmoves[i]=0;    //hardcode
        double* E=NULL;
        #if COUL_ON
        if(Ein==NULL) {
                E=malloc(ENERGY_LENGTH*sizeof(*E));
                assert(E!=NULL /*malloc*/);
                energy_full(c,s,E,ewald,p);
        }
        else    E=Ein;
        #endif
        for(int i=0;i<cStepsPerSweep;++i) {
                const uint32_t ind=pcg32_boundedrand_r(rng,N);
                if(c[ind].su==-1) {      //not in solute: cluster move
                        const int nInCluster=clusterMove(c,tc,s,E,ind,ewald,p,coreNum,rng,myerr);
                        if(nInCluster>0) ++nmoves[4];
                        ++nmoves[5];
                        if(i%cStepsPerSweep==0) fprintf(mvstats,"%d\n",nInCluster);
                }
                else {                   //in solute: solute move
                        const double roll=(double)pcg32_random_r(rng)/UINT32_MAX;
                        if(roll<=s[c[ind].su].mvProb) {
                                const int mvSuccess=suMove(c,tc,s,ts,u,tu,E,c[ind].su,ewald,inshlopt,rotation,sweepsMult,kTeff,p,coreNum,rng,myerr);
                                if(mvSuccess>0) ++nmoves[2];
                                ++nmoves[3];
                        }
                }
        }
        double* chg=checkCharge(c,N,NULL);
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites;
        const int delPlus=chg[0]-nPlusSvSites[coreNum];
        const int delMinus=chg[1]-nMinusSvSites[coreNum];
        const int delZero=chg[2]-nZeroSvSites[coreNum];
        if(delPlus!=0 || delMinus!=0 || delZero!=0) {
                fprintf(myerr,"MCstep_cluster: delPlus:%d delMinus:%d delZero:%d ; move refused\n",delPlus,delMinus,delZero);
        }
        free(chg);
        return nmoves;
}

/*******************************/
/********    SPOOF      ********/
/*******************************/

//compare trj header with input N, L
//loop:
//-read step
//-put config onto lattice
//-calc & output E
//cleanup
void runSpoof(lattice c, su s, double* ewald, para p, umb u, FILE* f, FILE* myout) {
        int nlines,nchars;
        const int niter=prepForReadTrj(&nlines,&nchars,f);
        extern int N;
        char** lines=malloc(N*sizeof(*lines));
        char* linesHolder=malloc(2*nchars/niter); //double avg amt of chars per iter
        assert(lines!=NULL && linesHolder!=NULL /*malloc*/);
        *lines=linesHolder;
        
        fprintf(myout,"#step\tE_tot\tE_isvsv\tE_isvsu\tE_isusu\tE_csrsvsv\tE_csrsvsu\tE_csrsusu\tE_clr\tE_ihh\tE_isvh\tE_isuh\tE_umb\tE_self\n");
        for(int i=0;i<niter;++i) {
                const int ncharsPerIter=readTrjStep(&lines,linesHolder,f);
                fitrjStepToLat(c,lines);
                //lammpstrjStepToLat(c,lines);

                double E[ENERGY_LENGTH];
                energy_full(c,s,E,ewald,p);
                const double eUmb=energy_umbr(u);
                const double eSelf=-p->e_self;
                fprintf(myout,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,E[0]+E[1]+E[2]+E[3]+E[4]+E[5]+E[6]+E[7]+E[8]+E[9]+eUmb+eSelf,E[0],E[1],E[2],E[3],E[4],E[5],E[6],E[7],E[8],E[9],eUmb,eSelf);
        }
        free(linesHolder);
        free(lines);
}

int prepForReadTrj(int* nlines, int* nchars, FILE* f) {
        const int lbuf=1024;
        char* buf=malloc(lbuf);
        assert(buf!=NULL /*malloc*/);
        memset(buf,'\0',lbuf);

        *nlines=0;
        *nchars=0;
        assert(fseek(f,0,SEEK_SET) == 0);
        while(fgets(buf,lbuf,f) != NULL) {
                const int l=strlen(buf);
                const char end=buf[l-1];
                if(end!='\n' && end!='\0') {
                        fprintf(stderr,"prepForTrjRead: buffer overflow\n");
                        exit(1);
                }
                *(nchars)+=l;
                ++(*nlines);
        }
        assert(fseek(f,0,SEEK_SET) == 0);

        extern int N,d;
        extern int* L;
        int Ltrj[d];
        int i=0;
        while(fgets(buf,lbuf,f) != NULL) {
                if(i==3) assert(N==atoi(buf) /*diff. inp N and trj N*/);
                else if(i>=5 && i<5+d) Ltrj[i-5]=atoi(buf+2);
                else if(i==5+d) {
                        int ct=0;
                        for(int e=0;e<d;++e) if(L[e]==Ltrj[e]) ++ct;
                        if(ct!=d) {
                                fprintf(stderr,"prepForTrjRead: L & trj L do not match\nL\tL_trj\n");
                                for(int e=0;e<d;++e) fprintf(stderr,"%d\t%d\n",L[e],Ltrj[e]);
                                exit(1);
                        }
                        break;
                }
                ++i;
        }
        free(buf);      buf=NULL;
        const int header=9;     //hardcode
        return *nlines/(N+header);
}

int readTrjStep(char*** lines, char* linesHolder, FILE* f) {
        const int lbuf=1024;
        char* buf=malloc(lbuf);
        assert(buf!=NULL /*malloc*/);
        memset(buf,'\0',lbuf);
        extern int N;
        int i=0,l=0,tmp=0;
        while(fgets(buf,lbuf,f) != NULL) {
                tmp=strlen(buf);
                assert(buf[tmp-1]=='\n' /*buffer overflow*/);
                buf[tmp-1]='\0';
                (*lines)[i]=linesHolder+l;
                *(*lines)[i]='\0'; //so strncat will know where to cat
                strncat((*lines)[i],buf,tmp);
                l+=tmp;
                if(++i>=N) break;
        }
        const int header=9;     //hardcode
        i=0;
        while(fgets(buf,lbuf,f) != NULL) if(++i>=header) break;
        free(buf);
        return l;
}

int findSpace(char* s, int spaceNum) {
        const int l=strlen(s);
        int spaceCt=0;
        for(int i=0;i<l;++i) {
                if(spaceCt==spaceNum) return i;
                if(s[i]==' ') ++spaceCt;
        }
        return -1;
}

//assumes .lammpstrj is in row-major order;
//dumpLammps() outputs in this order
void lammpstrjStepToLat(lattice c, char** lines) {
        extern int N;
        for(int i=0;i<N;++i) {
                const int j=findSpace(lines[i],1);
                const double v=atof(lines[i]+j);
                *c[i].val = (v==1.0) ? v : -1.0;
        }
}

//assumes .fitrj is in row-major order;
//dumpLammps() outputs in this order
void fitrjStepToLat(lattice c, char** lines) {
        extern int N;
        for(int i=0;i<N;++i) {
                const int j=findSpace(lines[i],1);
                const double v=atoi(lines[i]+j);
                if(v==0) {
                        *c[i].val=-1.0;
                }
                else if (v==1) {
                        *c[i].val=1.0;
                }
                else if (v==2) {
                        *c[i].val=0.0;
                }
        }
}
