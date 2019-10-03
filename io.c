#include "io.h"

//input and output functions
//such as for loading .data & .para files or for
//dumping .data, .para, .lammpstrj, or .fitrj files

//******************************
//*******       IO      ********
//******************************

//dumps .para files out. pairs with readSetParameters() output
void dumpParameters(para p, su s, umb u, int coreNum, FILE* f, FILE* myerr) {
        extern int N,Nsu,d,shlL,types,inshlopt;
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites,*nPlusSuSites,*nMinusSuSites,*nHydrSuSites;
        extern double plusCharge,minusCharge;
        extern int* L;
        extern int* bc;
        fprintf(f,"force no\n"); //hardcode: calc. calc'able param.s rather than reading
        fprintf(f,"N %d\n"                                                                                 \
                  "Nsu %d\n"                                                                                \
                  "d %d\n"                                                                                   \
                  "shlL %d\n"                                                                                 \
                  "nPlusSvSites %d\n"                                                                          \
                  "nMinusSvSites %d\n"                                                                          \
                  "nZeroSvSites %d\n"                                                                            \
                  "nPlusSuSites %d\n"                                                                             \
                  "nMinusSuSites %d\n"                                                                             \
                  "nHydrSuSites %d\n"                                                                               \
                  "plusCharge %f\n"                                                                                  \
                  "minusCharge %f\n"                                                                                  \
                  "types %d\n"                                                                                         \
                  "inshlopt %d\n",                                                                                      \
                  N,Nsu,d,shlL,nPlusSvSites[coreNum],nMinusSvSites[coreNum],nZeroSvSites[coreNum],                       \
                  nPlusSuSites[coreNum],nMinusSuSites[coreNum],nHydrSuSites[coreNum],plusCharge,minusCharge,types,inshlopt);
        fprintf(f,"L ");
        for(int e=0;e<d;++e) fprintf(f,"%d ",L[e]);
        fprintf(f,"\n");
        fprintf(f,"bc ");
        for(int e=0;e<d;++e) fprintf(f,"%d %d ",bc[2*e],bc[2*e+1]);
        fprintf(f,"\n");
        fprintf(f,"J %f\n"                                                                \
                  "Q %f\n"                                                                 \
                  "kT %f\n"                                                                 \
                  "beta %f\n"                                                                \
                  "betaEff_clst %f\n"                                                         \
                  "sigma %f\n"                                                                 \
                  "coulCutoff %f\n"                                                             \
                  "fftOn %d\n"                                                                   \
                  "e_cut %f\n"                                                                    \
                  ,p->J,p->Q,p->kT,p->beta,p->betaEff_clst,p->sigma,p->coulCutoff,p->fftOn,p->e_cut);
        if(Nsu>0)   dumpSu(s,f);
        if(u!=NULL) dumpUmbrellas(u,s,p,f,myerr);
        fflush(f);
}

//reads .para files in. pairs with dumpParameters() output
//
//!!REQUIRES d COMES BEFORE L,bc!!
int* readSetParameters(para p, char**** linesByType, int coreNum, FILE* f, FILE* myerr) {
        extern int N,Nsu,d,shlL,types,inshlopt;
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites,*nPlusSuSites,*nMinusSuSites,*nHydrSuSites;
        extern double plusCharge,minusCharge;
        extern int* L;
        extern int* bc;

        const int posTypes=2; //hardcode: solutes & umbrellas
        //first nLines.. is length; first of lines.. is NULL
        int* nLinesByType=malloc((1+posTypes)*sizeof(*nLinesByType));
        *linesByType=malloc((1+posTypes)*sizeof(**linesByType));
        assert(nLinesByType!=NULL || *linesByType!=NULL /*malloc*/);
        for(int i=0;i<1+posTypes;++i) {
                (*linesByType)[i]=NULL;
                nLinesByType[i]=0;
        }

        int i=0;
        const int lbuf=32768;
        char buf[lbuf];
        memset(buf,'\0',lbuf);
        char** words=NULL;
        int nwords[1]={0};

        long int ph=0l;
        int force=0;    //switch: 0 -> calc. certain param.s; 1 -> force use read input
        while(fgets(buf, lbuf, f) != NULL) {
                words=split(buf, " ", nwords);
                if( strcmp(words[0],"force")==0  && strcmp(words[1],"yes")==0 ) {
                        force=1;
                }
                else if( strcmp(words[0],"N")==0 )   N=atoi(words[1]);
                else if( strcmp(words[0],"Nsu")==0 ) Nsu=atoi(words[1]);
                else if( strcmp(words[0],"d")==0 ) {
                        d=atoi(words[1]);
                }
                else if( strcmp(words[0],"shlL")==0 )          shlL=atoi(words[1]);
                else if( strcmp(words[0],"nPlusSvSites")==0 )  nPlusSvSites[coreNum]=atoi(words[1]);
                else if( strcmp(words[0],"nMinusSvSites")==0 ) nMinusSvSites[coreNum]=atoi(words[1]);
                else if( strcmp(words[0],"nZeroSvSites")==0 )  nZeroSvSites[coreNum]=atoi(words[1]);
                else if( strcmp(words[0],"nPlusSuSites")==0 )  nPlusSuSites[coreNum]=atoi(words[1]);
                else if( strcmp(words[0],"nMinusSuSites")==0 ) nMinusSuSites[coreNum]=atoi(words[1]);
                else if( strcmp(words[0],"nHydrSuSites")==0 )  nHydrSuSites[coreNum]=atoi(words[1]);
                else if( strcmp(words[0],"plusCharge")==0 )    plusCharge=atof(words[1]);
                else if( strcmp(words[0],"minusCharge")==0 )   minusCharge=atof(words[1]);
                else if( strcmp(words[0],"types")==0 )         types=atoi(words[1]);
                else if( strcmp(words[0],"inshlopt")==0 )      inshlopt=atoi(words[1]);
                else if( strcmp(words[0],"L")==0 ) {
                        for(int e=0;e<d;++e) L[e]=atoi(words[1+e]);
                }
                else if( strcmp(words[0],"bc")==0 ) {
                        for(int e=0;e<2*d;++e) bc[e]=atoi(words[1+e]);
                }
                else if( strcmp(words[0],"J")==0 )                 p->J=atof(words[1]);
                else if( strcmp(words[0],"Q")==0 )                 p->Q=atof(words[1]);
                else if( strcmp(words[0],"kT")==0 )                p->kT=atof(words[1]);
                else if( strcmp(words[0],"beta")==0 && force==1 )  p->beta=atof(words[1]);
                else if( strcmp(words[0],"betaEff_clst")==0 )      p->betaEff_clst=atof(words[1]);
                else if( strcmp(words[0],"sigma")==0 )             p->sigma=atof(words[1]);
                else if( strcmp(words[0],"coulCutoff")==0 )        p->coulCutoff=atof(words[1]);
                else if( strcmp(words[0],"fftOn")==0 )             p->fftOn=atoi(words[1]);
                else if( strcmp(words[0],"e_cut")==0 && force==1 ) p->e_cut=atof(words[1]);
                else if( strcmp(words[0],"Solutes\n")==0 ) {
                        ++nLinesByType[0];
                        ph=ftell(f);
                        assert(fseek(f,0l,SEEK_SET)==0);
                        nLinesByType[1]=readLines("Solutes",*(linesByType)+1,(char**)NULL,0,f,myerr);
                        assert(fseek(f,ph,SEEK_SET)==0);
                }
                else if( strcmp(words[0],"Umbrellas\n")==0 ) {
                        ++nLinesByType[0];
                        ph=ftell(f);
                        assert(fseek(f,0l,SEEK_SET)==0);
                        nLinesByType[2]=readLines("Umbrellas",*(linesByType)+2,(char**)NULL,0,f,myerr);
                        assert(fseek(f,ph,SEEK_SET)==0);
                }
                i+=1;
                free(*words);   free(words);    words=NULL;
        }

        //calculable parameters (0.0 input -> calculate based on other given):
        if(p->beta==0.0 || force==0)    p->beta=1.0/(p->kT);
        if(p->betaEff_clst==0.0 || force==0) {
                if(p->kT<=5.0) p->betaEff_clst=1.0/5.0;     //hardcode from G&V 2001
                else           p->betaEff_clst=1.0/(p->kT); //T_eff ~5kT
        }

        extern double pi;
        #if SHIFT_COUL
        if(p->e_cut==0.0 || force==0) p->e_cut=erfc(p->coulCutoff)/(p->coulCutoff*p->sigma);
        #endif
        #if !SHIFT_COUL
        p->e_cut=0.0;
        #endif
        #if COUL_ON
        p->e_self=N*(p->Q)*(p->e_cut/2.0+1.0/(sqrt(pi)*p->sigma));
        #endif
        #if !COUL_ON
        p->e_self=0.0;
        #endif
        fprintf(myerr,"readSetParameters: done. e_self: %f\n",p->e_self);

        return nLinesByType;
}


//Output data file which can be read by readData()
//to allow a simulation to be restarted with an
//equivalent configuration.
//Data file of format:
//
//Frustrated Ising data file v#.##
//
//N cells
//N_t types
//N_s solutes
//
//d dimensions
//e_0lo e_0hi
//e_1lo e_1hi
//...
//e_(d-1)lo e_(d-1)hi
//
//Cells
//i type_i e_0i e_1i ... e_(d-1)i solute_i
//...
//
//Solutes
//i surfType_i com_0i ... com_(d-1)i step_i nsites_i ninshl_i linDim_i hlw_i nMotion_i shape_i
//...
void dumpData(lattice c, su s, FILE* f) {
        extern int N,types,Nsu,d;
        extern int* L;
        fprintf(f,"Frustrated Ising data file v1.00\n\n" \
                  "%d cells\n"                            \
                  "%d types\n"                             \
                  "%d solutes\n\n", N, types, Nsu           );
        fprintf(f,"%d dimensions\n",d);
        for(int e=0;e<d;++e) fprintf(f,"%d\t%d\n",0,L[e]); //lower bound hardcoded 0
                
        fprintf(f,"\nCells\n");
        for(int i=0;i<N;++i) {
                fprintf(f,"%d\t%d\t",i,(int)lround(*c[i].val));
                int pos[d];
                getPos(i,pos);
                for(int e=0;e<d;++e) fprintf(f,"%d\t",pos[e]);
                fprintf(f,"%d\n",c[i].su);
        }
        if(Nsu>0) dumpSu(s,f);
        fflush(f);
}

//Read data from .data file output by dumpData() into memory.
int readData(char*** lines, FILE* f) {
        const int lbuf=32768;
        char* buf=malloc(lbuf);
        assert(buf!=NULL /*malloc*/);
        memset(buf,'\0',lbuf);
        int nlines=0, totchars=0;
        int fseekret=fseek(f,0,SEEK_SET);
        assert(fseekret==0);
        while(fgets(buf,lbuf,f) != NULL) {
                totchars+=strlen(buf);
                ++nlines;
        }
        fseekret=fseek(f,0,SEEK_SET);
        assert(fseekret==0);
        *lines=malloc(nlines*sizeof(**lines));
        char* linesStore=malloc(totchars*sizeof(*linesStore));
        assert(*lines!=NULL && linesStore!=NULL /*malloc*/);
        int i=0,l=0,tmp=0;
        *(*lines)=linesStore;
        **(*lines)='\0';
        while(fgets(buf,lbuf,f) != NULL) {
                (*lines)[i]=linesStore+l;
                tmp=strlen(buf);
                assert(buf[tmp-1]=='\n' /*overflow buffer*/);
                buf[tmp-1]='\0'; //replace '\n'
                strncat((*lines)[i],buf,tmp);
                l+=tmp;
                i+=1;
        }
        free(buf);
        return nlines;
}

//dump unwrapped solute i center of mass
void dumpSuCom(su s, int i, FILE* f) {
        extern int d;
        for(int e=0;e<d;++e) fprintf(f,"%f ",s[i].unwrappedCom[e]);
}

//load solutes from .data or .para file
void loadSu(char** lines, int nlines, su* s, su* ts, char** sNH, double** comH, double** unwrappedComH,          \
                 double** orientationHH, double*** orientationH, double** suRelPos, double*** suRelPosPtrs,       \
                 double** shlRelPos, double*** shlRelPosPtrs, lattice** suCurPos, lattice** shlCurPos, bool** hydrH, FILE* myerr) {
        const int maxShapeLett=10; //"rectangle" + "\0"
        extern int Nsu,d;
        assert(nlines>0 && lines!=NULL /*passed empty lines*/);
        assert(strcmp(lines[0],"Solutes")==0 /*incorrect first line*/);
        assert(nlines == Nsu+1 /*check Nsu*/);
        (*s)=malloc(Nsu*sizeof(**s));
        (*ts)=malloc(Nsu*sizeof(**ts));
        (*comH)=malloc(2*(Nsu*d)*sizeof(**comH));
        (*unwrappedComH)=malloc(2*(Nsu*d)*sizeof(**unwrappedComH));
        (*orientationH)=malloc(2*(Nsu*d)*sizeof(**orientationH));
        (*orientationHH)=malloc(2*(Nsu*d*d)*sizeof(**orientationHH));
        (*sNH)=malloc(Nsu*maxShapeLett);
        assert((*s)!=NULL && (*ts)!=NULL && (*comH)!=NULL && (*unwrappedComH)!=NULL &&   \
               (*orientationH)!=NULL && (*orientationHH)!=NULL && (*sNH)!=NULL /*malloc*/ );
        char* sNHPH=(*sNH);
        memset((*sNH),'\0',Nsu*maxShapeLett);
        fprintf(myerr,"loadSu: got past solutes malloc\n");
        //line-by-line load
        char** words;
        int nwords[1]={0};
        extern int shlL;
        const int nSuPara=9+d+d*d;  //hardcode
        char** hydrStrings=malloc(Nsu*sizeof(*hydrStrings));
        assert(hydrStrings!=NULL);
        for(int i=0;i<Nsu;++i) {
                words=split(lines[i+1], " ", nwords);
                assert(*nwords==nSuPara /*incorrect # lines*/);
                uint8_t j=0u;

                int numTypeAndSize;
                char** typeAndSize=split(words[++j], ",", &numTypeAndSize);
                strcpy((*sNH),typeAndSize[0]);
                (*ts)[i].shape = (*s)[i].shape = (*sNH);
                *sNH+=strlen(*sNH)+1;
                if(strcmp((*s)[i].shape,"cube")==0) {
                        int* linDimHolder=malloc(1*sizeof(*linDimHolder));
                        assert(linDimHolder!=NULL /*malloc*/);
                        *linDimHolder=atoi(typeAndSize[1]);
                        (*ts)[i].linDim = (*s)[i].linDim = linDimHolder;
                        (*s)[i].nsites=1;
                        (*s)[i].ninshl=1;
                        for(int e=0;e<d;++e) {
                                (*s)[i].nsites *= (*s)[i].linDim[0];
                                (*s)[i].ninshl *= ((*s)[i].linDim[0] + 2*shlL);
                        }
                        (*s)[i].ninshl -= (*s)[i].nsites;
                        (*ts)[i].nsites=(*s)[i].nsites;
                        (*ts)[i].ninshl=(*s)[i].ninshl;
                }
                else if(strcmp((*s)[i].shape,"rectangle")==0 || strcmp((*s)[i].shape,"rect")==0) {
                        int* linDimHolder=malloc(d*sizeof(*linDimHolder));
                        assert(linDimHolder!=NULL /*malloc*/);
                        for(int e=0;e<d;++e) linDimHolder[e]=atoi(typeAndSize[1+e]);
                        (*ts)[i].linDim = (*s)[i].linDim = linDimHolder;
                        (*s)[i].nsites=1;
                        (*s)[i].ninshl=1;
                        for(int e=0;e<d;++e) {
                                (*s)[i].nsites *= (*s)[i].linDim[e];
                                (*s)[i].ninshl *= ((*s)[i].linDim[e] + 2*shlL);
                        }
                        (*s)[i].ninshl -= (*s)[i].nsites;
                        (*ts)[i].nsites=(*s)[i].nsites;
                        (*ts)[i].ninshl=(*s)[i].ninshl;

                        int** rect=buildRectangle((*s)[i].linDim);
                        int* Lcube=NULL;
                        int border[d];
                        int ninterface=0;
                        //fprintf(myerr,"shlL: %d\n",shlL);
                        for(int e=0;e<d;++e) border[e]=shlL;
                        int** interfacePos=getShapeInterfacePos(rect,(*s)[i].nsites,&ninterface,&Lcube,border);
                        if((*s)[i].ninshl!=ninterface) {
                                (*s)[i].ninshl=ninterface;
                                (*ts)[i].ninshl=ninterface;
                        }
                        free(Lcube);
                        free(*rect); free(rect);
                        if(ninterface>0) {
                                free(*interfacePos);
                                free(interfacePos);
                        }
                }
                else if(strcmp((*s)[i].shape,"sphere")==0 || strcmp((*s)[i].shape,"sph")==0) {
                        int* linDimHolder=malloc(1*sizeof(*linDimHolder));
                        assert(linDimHolder!=NULL /*malloc*/);
                        *linDimHolder=atoi(typeAndSize[1]);
                        (*ts)[i].linDim = (*s)[i].linDim = linDimHolder;
                        //find number of cells in sphere
                        int nSph=0;
                        const int lCube=2*(*s)[i].linDim[0];
                        int** cube=buildCube(lCube);
                        const int nCube=lCube*lCube*lCube;
                        const double R2=(*s)[i].linDim[0]*(*s)[i].linDim[0];
                        for(int j=0;j<nCube;++j) {
                                double dist2=0.0;
                                for(int e=0;e<d;++e) {
                                        const double re=cube[j][e]-(lCube-1)/2.0;
                                        dist2+=re*re;
                                }
                                if(dist2<=R2) ++nSph;
                        }
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
                                        }
                                        ++k;
                                }
                        }
                        //build interface of sphere
                        int* Lcube=NULL;
                        int border[d];
                        int ninterface;
                        for(int e=0;e<d;++e) border[e]=shlL;
                        int** interfacePos=getShapeInterfacePos(sphere,nSph,&ninterface,&Lcube,border);

                        (*s)[i].nsites=nSph;
                        (*s)[i].ninshl=ninterface;
                        (*ts)[i].nsites=(*s)[i].nsites;
                        (*ts)[i].ninshl=(*s)[i].ninshl;
                        free(Lcube);
                        free(*cube); free(cube);
                        free(*interfacePos); free(interfacePos);
                }
                else {
                        fprintf(myerr,"shape '%s' not yet implemented. try cube, rect, or sphere. i:%d\n",(*s)[i].shape,i);
                        //fprintf(myerr,"%d\n",strcmp((*s)[i].shape,"rect")==0);
                        //fprintf(myerr,"%d\n",strcmp((*s)[i].shape,"rectangle")==0 || strcmp((*s)[i].shape,"rect")==0);
                        //fprintf(myerr,"%s\n",typeAndSize[0]);
                        exit(1);
                }

                (*s)[i].orientation=(*orientationH)+i*d;
                (*ts)[i].orientation=(*orientationH)+i*d+Nsu*d;
                for(int e=0;e<d;++e) {
                        (*orientationH)[e+i*d]=(*orientationHH)+(e+i*d)*d;
                        (*orientationH)[e+i*d+Nsu*d]=(*orientationHH)+(e+i*d+Nsu*d)*d;
                        for(int k=0;k<d;++k) {
                                (*ts)[i].orientation[e][k]=(*s)[i].orientation[e][k]=atof(words[++j]);
                        }
                }

                (*ts)[i].surfType = (*s)[i].surfType = atof(words[++j]);

                (*s)[i].com=(*comH)+i*d;
                (*s)[i].unwrappedCom=(*unwrappedComH)+i*d;
                (*ts)[i].com=(*comH)+i*d+Nsu*d;
                (*ts)[i].unwrappedCom=(*unwrappedComH)+i*d+Nsu*d;
                for(int e=0;e<d;++e) {
                        (*comH)[e+i*d+Nsu*d]=(*comH)[e+i*d]=(*unwrappedComH)[e+i*d+Nsu*d]=(*unwrappedComH)[e+i*d]=atof(words[++j]);
                }
                wrapIntoL_double((*comH)+i*d,NULL); //PBC
                wrapIntoL_double((*comH)+i*d+Nsu*d,NULL); //PBC
                
                (*ts)[i].suStep   = (*s)[i].suStep   = atoi(words[++j]);
                (*ts)[i].hlw      = (*s)[i].hlw      = atoi(words[++j]);
                (*ts)[i].nMotion  = (*s)[i].nMotion  = atoi(words[++j]);
                (*ts)[i].isiMult  = (*s)[i].isiMult  = atof(words[++j]);
                (*ts)[i].mvProb   = (*s)[i].mvProb   = atof(words[++j]);

                char* hydrString=malloc(((*s)[i].nsites+1)*sizeof(*hydrString));
                assert(hydrString!=NULL);
                ++j;
                //capital F -> none hydrophobic/organic
                char pfix[2];
                pfix[0]=words[j][0];
                pfix[1]='\0';
                if(strcmp(pfix,"F")==0) {
                        for(int ii=0;ii<(*s)[i].nsites;++ii) hydrString[ii]='f';
                        hydrString[(*s)[i].nsites]='\0';
                        hydrStrings[i]=hydrString;
                }
                else {
                        for(int ii=0;ii<(*s)[i].nsites;++ii) hydrString[ii]=words[j][ii];
                        hydrString[(*s)[i].nsites]='\0';
                        hydrStrings[i]=hydrString;
                }
                
                free(*words);  free(words);   words=NULL;
                free(*typeAndSize);     free(typeAndSize);
        }
        (*sNH)=sNHPH;
        #if TEST_BUILD
        fprintf(myerr,"su #\tnsites\tninshl\n");
        for(int i=0;i<Nsu;++i) {
                fprintf(myerr,"%d\t%d\t%d\n",i,(*s)[i].nsites,(*s)[i].ninshl);
        }
        #endif

        extern int Ntotsusites; //extern int
        Ntotsusites=0;
        int Ntotinshl=0;
        for(int i=0;i<Nsu;++i) {
                Ntotsusites+=(*s)[i].nsites;
                Ntotinshl+=(*s)[i].ninshl;
        }

        *suRelPos=malloc(1*Ntotsusites*d*sizeof(**suRelPos));
        *suRelPosPtrs=malloc(1*Ntotsusites*sizeof(**suRelPosPtrs));
        *shlRelPos=malloc(1*Ntotinshl*d*sizeof(**shlRelPos));
        *shlRelPosPtrs=malloc(1*Ntotinshl*sizeof(**shlRelPosPtrs));
        *suCurPos=malloc(2*Ntotsusites*sizeof(**suCurPos));
        *shlCurPos=malloc(2*Ntotinshl*sizeof(**shlCurPos));
        *hydrH=malloc(1*Ntotsusites*sizeof(**hydrH));
        assert((*suRelPos!=NULL) && (*suRelPosPtrs!=NULL)  &&         \
               (*shlRelPos!=NULL) && (*shlRelPosPtrs!=NULL) &&         \
               (*shlCurPos!=NULL) && (*suCurPos!=NULL) && (*hydrH!=NULL));
        fprintf(myerr,"loadSu: got past su posBlock mallocs\n");
        int ct=0;
        int ctShl=0;
        for(int i=0;i<Nsu;++i) {
                (*s)[i].currPos=(*suCurPos)+ct;
                (*ts)[i].currPos=(*suCurPos)+ct+Ntotsusites;
                for(int j=0;j<(*s)[i].nsites;++j) {
                        (*suRelPosPtrs)[ct+j]=(*suRelPos)+(ct+j)*d;
                }
                (*ts)[i].relPos=(*s)[i].relPos=(*suRelPosPtrs)+ct;
                (*ts)[i].hydrophobic=(*s)[i].hydrophobic=(*hydrH)+ct;
                ct+=(*s)[i].nsites;

                (*s)[i].inShlCurrPos=(*shlCurPos)+ctShl;
                (*ts)[i].inShlCurrPos=(*shlCurPos)+ctShl+Ntotinshl;
                for(int j=0;j<(*s)[i].ninshl;++j) {
                        (*shlRelPosPtrs)[ctShl+j]=(*shlRelPos)+(ctShl+j)*d;
                }
                (*ts)[i].inShlRelPos=(*s)[i].inShlRelPos=(*shlRelPosPtrs)+ctShl;
                ctShl+=(*s)[i].ninshl;
        }
        extern bool hydrophobicExist;
        for(int i=0;i<Nsu;++i) {
                for(int j=0;j<(*s)[i].nsites;++j) {
                        if(hydrStrings[i][j]=='t' || hydrStrings[i][j]=='T' || hydrStrings[i][j]=='1') {
                                (*s)[i].hydrophobic[j]=true;
                                hydrophobicExist=true;
                        }
                        else if(hydrStrings[i][j]=='f' || hydrStrings[i][j]=='F' || hydrStrings[i][j]=='0') {
                                (*s)[i].hydrophobic[j]=false;
                        }
                        else if(hydrStrings[i][j]=='\0') {
                                fprintf(myerr,"loadSu: error, reached end of hydrStrings[%d][%d] w nsites %d. hydrStrings[i]: %s. exiting\n",i,j,(*s)[i].nsites,hydrStrings[i]);
                                exit(1);
                        }
                        else {
                                fprintf(myerr,"loadSu: error, hydrStrings[%d][%d] incorrect input %c exiting\n",i,j,hydrStrings[i][j]);
                                exit(1);
                        }
                }
        }
        #if TEST_BUILD
        for(int i=0;i<Nsu;++i) {
                for(int j=0;j<(*s)[i].nsites;++j) {
                        fprintf(myerr,"%c\t",hydrStrings[i][j]);
                        if((*s)[i].hydrophobic[j]==false) fprintf(myerr,"false\n");
                        else fprintf(myerr,"true\n");
                }
        }
        #endif
        for(int i=0;i<Nsu;++i) free(hydrStrings[i]);
        free(hydrStrings);
}

//dump solutes from memory into .para or .data file
void dumpSu(su s, FILE* f) {
        extern int Nsu,d;
        fprintf(f,"\nSolutes\n");
        for(int i=0;i<Nsu;++i) {
                fprintf(f,"%d %s",i,s[i].shape);
                //assert(strcmp(s[i].shape,"sphere")!=0 && strcmp(s[i].shape,"custom")!=0 /*tb implemented*/);
                if(strcmp(s[i].shape,"cube")==0) {
                        fprintf(f,",%d ",s[i].linDim[0]);
                }
                else if(strcmp(s[i].shape,"rectangle")==0 || strcmp(s[i].shape,"rect")==0) {
                        for(int e=0;e<d;++e) fprintf(f,",%d",s[i].linDim[e]);
                        fprintf(f," ");
                        for(int e=0;e<d;++e) {
                                for(int j=0;j<d;++j) {
                                        fprintf(f,"%d ",(int)s[i].orientation[e][j]);
                                }
                        }
                }
                else if(strcmp(s[i].shape,"sphere")==0 || strcmp(s[i].shape,"sph")==0) {
                        fprintf(f,",%d ",s[i].linDim[0]);
                        for(int e=0;e<d;++e) {
                                for(int j=0;j<d;++j) {
                                        fprintf(f,"%d ",(int)s[i].orientation[e][j]);
                                }
                        }
                }
                else {
                        fprintf(stderr,"dumpSu: didn't recognize su type '%s'. exiting\n",s[i].shape);
                        exit(1);
                }

                fprintf(f,"%f ",s[i].surfType);
                dumpSuCom(s,i,f);
                fprintf(f,"%f %d %d %f %f ", \
                        s[i].suStep,           \
                        s[i].hlw,               \
                        s[i].nMotion,            \
                        s[i].isiMult,             \
                        s[i].mvProb                );
                for(int j=0;j<s[i].nsites;++j) {
                        if(s[i].hydrophobic[j]==false) fprintf(f,"f");
                        else                           fprintf(f,"t");
                }
                fprintf(f,"\n");
        }
        fprintf(f,"End Solutes\n");
        fflush(f);
}

//dump all unwrapped solute centers of mass
void dumpAllSuComs(su s, FILE*f) {
        extern int Nsu;
        for(int i=0;i<Nsu;++i) {
                dumpSuCom(s,i,f);
                fprintf(f,"\n");
        }
        fprintf(f,"\n");
}

//load umbrellas from .para file
void loadUmbrellas(char** lines, int nlines, umb u, su s, para p, FILE* myumbout, FILE* myerr) {
        assert(strcmp(lines[0],"Umbrellas") == 0 /*wrong starting line*/);
        const int maxUnitsLett=3+1;
        u->npairs=nlines-1;
        u->pair=malloc((u->npairs)*sizeof(*u->pair));
        su* pairHolder=malloc(2*(u->npairs)*sizeof(*pairHolder));
        u->k=malloc((u->npairs)*sizeof(*u->k));
        u->energyIOUnits=malloc((u->npairs)*sizeof(*u->energyIOUnits));
        char* energyIOUnitsHolder=malloc((u->npairs * maxUnitsLett)*sizeof(energyIOUnitsHolder));
        u->r0=malloc((u->npairs)*sizeof(*u->r0));
        assert(u->pair!=NULL && pairHolder!=NULL && u->k!=NULL && u->energyIOUnits!=NULL
               && energyIOUnitsHolder!=NULL && u->r0!=NULL /*malloc*/);
        memset(energyIOUnitsHolder,'\0',u->npairs * maxUnitsLett);
        //line-by-line load
        char** words;
        char* ptr=energyIOUnitsHolder;
        int nwords[1]={0};
        for(int i=0;i<(u->npairs);++i) {
                u->pair[i]=(pairHolder+2*i);
                words=split(lines[i+1]," ",nwords);
                u->pair[i][0]=(s+atoi(words[1]));
                u->pair[i][1]=(s+atoi(words[2]));
                u->k[i]=atof(words[3]);
                strcpy(ptr,words[4]);
                u->energyIOUnits[i]=ptr;
                ptr+=strlen(ptr)+1;
                if(strcmp(u->energyIOUnits[i],"kT")==0) u->k[i]*=p->kT; //pump up by factor of kT internally
                else if(strcmp(u->energyIOUnits[i],"raw")==0);
                else {
                        fprintf(myerr,"loadUmbrellas: received energyIOUnits energy of '%s'; acceptable inputs are 'kT' or 'raw'\n",u->energyIOUnits[i]);
                        exit(1);
                }
                u->r0[i]=atof(words[5]);
                
                free(*words);   free(words);    words=NULL;
        }
        //setup header for myumbout
        if(myumbout!=NULL) {
                fprintf(myumbout,"#umbN\tsuA<->suB\tk\tr0");
                for(int i=0;i<(u->npairs);++i) {
                        fprintf(myumbout,"#%d\t%d<->%d\t%f\t%f\n",i,(int)(u->pair[i][0]-s),(int)(u->pair[i][1]-s),u->k[i],u->r0[i]);
                }
                fprintf(myumbout,"#ts\t");
                for(int i=0;i<(u->npairs);++i) fprintf(myumbout,"umb%d r\t",i);
                fprintf(myumbout,"\n");
        }

        fprintf(myerr,"loadUmbrellas: load success. dumping for confirmation:\n");
        dumpUmbrellas(u,s,p,myerr,myerr);
}

//dump umbrellas to .para file
void dumpUmbrellas(umb u, su s, para p, FILE* f, FILE* myerr) {
        extern int d;
        fprintf(f,"\nUmbrellas\n");
        for(int i=0;i<u->npairs;++i) {
                const unsigned int s_i0=(unsigned int)(u->pair[i][0]-s);
                const unsigned int s_i1=(unsigned int)(u->pair[i][1]-s);
                double k=u->k[i];
                const char* energyIOUnits=u->energyIOUnits[i];
                if(strcmp(energyIOUnits,"kT")==0) k/=p->kT; //pumped up by factor of kT internally
                else if(strcmp(energyIOUnits,"raw")==0);
                else {
                        fprintf(myerr,"dumpUmbrellas: received energyIOUnits energy of '%s'; acceptable inputs are 'kT' or 'raw'. This should have been caught at loadUmbrellas\n",energyIOUnits);
                }
                const double r0=u->r0[i];
                fprintf(f,"%d %u %u %f %s %f\n",i,s_i0,s_i1,k,energyIOUnits,r0);
        }
        fprintf(f,"End Umbrellas\n");
        fflush(f);
}

//read lines in from file f to lines
int readLines(char* header, char*** lines, char** inLines, int ninLines, FILE* f, FILE* myerr) {
        const int lh=strlen(header);
        char* headern=NULL; //ensure header ends with '\n'
        char headernHolder[lh+2];
        headern=headernHolder;
        memset(headern,'\0',lh+2);
        strncpy(headern,header,lh);
        if(headern[lh-1]!='\n') {
                headern[lh]='\n';
                headern[lh+1]='\0';
        }

        char endern[strlen("End ")+strlen(headern)+4];
        memset(endern,'\0',strlen("End ")+strlen(headern)+4);
        strncpy(endern,"End ",strlen("End "));
        strncat(endern,headern,strlen(headern));
        if(inLines!=NULL) { //read from inLines
                int i=0;
                while(strcmp(inLines[i],headern)!=0) {
                        i+=1;
                        if(i==ninLines) {
                                fprintf(myerr,"readLines: didn't find header (%s) in str mode\n",header);
                                return 0;
                        }
                }
                int j=i;
                while(strcmp(inLines[j],endern)!=0) {
                        j+=1;
                        if(j==ninLines) {
                                fprintf(myerr,"readLines: didn't find ender (End %s) in str mode\n",header);
                                break;
                        }
                }
                (*lines)=malloc((j-i)*sizeof(**lines));
                assert(*lines!=NULL /*malloc*/);
                for(int k=0;k<(j-i);++k) (*lines)[k]=inLines[i+k];

                return (j-i);
        }
        else if(f!=NULL) { //read from f
                const int lbuf=32768;
                char buf[lbuf];
                int fseekret;
                memset(buf,'\0',lbuf);
                long int ph=0l;
                while(fgets(buf,lbuf,f) != NULL) { //find beginning of header sect.
                        if(strcmp(buf,headern)!=0) ph=ftell(f);
                        else                       break;
                }
                assert(ph>0l /*didn't find header*/);
                fseekret=fseek(f,ph,SEEK_SET);
                assert(fseekret==0);
                int nlines=0, totchars=0;
                while(fgets(buf,lbuf,f) != NULL) { //find length of header sect.
                        if(strcmp(buf,endern) != 0) {
                                const int l=strlen(buf);
                                totchars+=l;
                                if(buf[l-1]=='\n') ++nlines;
                        }
                        else    break;
                }

                fseekret=fseek(f,ph,SEEK_SET);
                assert(fseekret==0);
                (*lines)=malloc(nlines*sizeof(**lines));
                char* linesStore=malloc(totchars+nlines); //extra nlines for str sep.s '\0'
                assert((*lines)!=NULL && linesStore!=NULL /*malloc*/);
                linesStore[0]='\0'; //so first call of strncat() <-> strncpy()
                int i=0,l=0,prevstrlen=-1;
                while(fgets(buf,lbuf,f) != NULL) { //read in header sect. including header
                        if(i>=nlines)    break;
                        (*lines)[i]=linesStore+l;
                        int tmp=strlen(buf);
                        strncat((*lines)[i],buf,tmp);
                        if(i!=0) (*lines)[i-1][prevstrlen-1]='\0';
                        prevstrlen=tmp;
                        l+=tmp;
                        ++i;
                }
                if((*lines)[i-1][prevstrlen-1]=='\n') (*lines)[i-1][prevstrlen-1]='\0';
                fprintf(myerr,"readLines(): returning with nlines=%d, from f mode lines:\n",nlines);
                #if TEST_BUILD
                for(int i=0;i<nlines;++i) fprintf(myerr,"%s\n",(*lines)[i]);
                #endif
                return nlines;
        }
        else {
                fprintf(myerr,"readLines: input method (char** or file) not given\n");
                exit(1);
        }
}

//Output in LAMMPS format for visualization in VMD;
//currently allows visualization for 2D & 3D lattices.
//Uses ctplus & ctminus because LAMMPS wants particles
//which do not change type, ie. if particle 'm' is
//declared as type 't', it must remain so throughout
//the simulation. But the lattice definition for this
//code has cells which remain at constant positions
//with spins that change orientation. ctplus&ctminus
//are a hack to make the output type compatible with
//this code's interal representation of the lattice.
void dumpLammps(lattice c, su s, int ts, int coreNum, FILE* f, FILE* myerr) {
        extern int N,d,Nsu;
        extern int *nPlusSvSites,*nMinusSvSites,*nZeroSvSites;
        int ctminus=0,ctplus=nMinusSvSites[coreNum],ctzero=nMinusSvSites[coreNum]+nPlusSvSites[coreNum];
        int ctsu[Nsu];
        int tmp=0;
        for(int i=0;i<Nsu;++i) {
                ctsu[i]=nPlusSvSites[coreNum]+nMinusSvSites[coreNum]+nZeroSvSites[coreNum]+tmp;
                tmp+=s[i].nsites;
        }

        extern int* L;
        extern int* bc;
        //LAMMPS style header
        fprintf(f,"ITEM: TIMESTEP\n%d\n",ts);
        fprintf(f,"ITEM: NUMBER OF ATOMS\n%d\n",N);
        fprintf(f,"ITEM: BOX BOUNDS ");
        for(int e=0;e<d;++e) {
                if(bc[2*e]==0 || bc[2*e+1]==0) fprintf(f,"ff ");
                else fprintf(f,"pp ");
        }
        if(d==2) fprintf(f,"pp ");
        fprintf(f,"\n");
        for(int e=0;e<d;++e) fprintf(f,"0 %d\n",L[e]);
        if(d==2) fprintf(f,"0 1\n");
        fprintf(f,"ITEM: ATOMS id type ");
        if(d==3 || d==2) fprintf(f,"x y z\n");
        else {
                for(int e=0;e<d;++e) fprintf(f,"e%d ",e);
                fprintf(f,"\n");
        }
        //cell/config
        const int nTypesBeforeSu = (nZeroSvSites==0) ? 2 : 3;
        for(int i=0;i<N;++i) {
                if(c[i].su!=-1)          fprintf(f,"%d %d ",ctsu[c[i].su]++,nTypesBeforeSu+c[i].su);
                else if(*c[i].val<0.0)   fprintf(f,"%d 0 ",ctminus++);
                else if(*c[i].val>0.0)   fprintf(f,"%d 1 ",ctplus++);
                else if(*c[i].val==0.0)  fprintf(f,"%d 2 ",ctzero++);
                else                     fprintf(myerr,"dumpLammps: got unexpected val: %f\n",*c[i].val);
                int pos[d];
                getPos(i,pos);
                for(int e=0;e<d;++e) fprintf(f,"%d ",pos[e]);
                if(d==2) fprintf(f,"0 "); //make 2D plane in 3D so VMD can visualize in LAMMPS format
                fprintf(f,"\n");
        }
        if(ctminus!=nMinusSvSites[coreNum] || ctplus!=nMinusSvSites[coreNum]+nPlusSvSites[coreNum] || ctzero!=nMinusSvSites[coreNum]+nPlusSvSites[coreNum]+nZeroSvSites[coreNum]) {
                fprintf(myerr,"dumpLammps: ts %d final ct.s\n\texpected\tactual\n"
                              "zero\t%d\t\t%d\nplus\t%d\t\t%d\nminus\t%d\t\t%d\n",
                              ts,nMinusSvSites[coreNum]+nPlusSvSites[coreNum]+nZeroSvSites[coreNum],
                              ctzero,nMinusSvSites[coreNum]+nPlusSvSites[coreNum],ctplus,nMinusSvSites[coreNum],ctminus);
                tmp=0;
                for(int i=0;i<Nsu;++i) {
                        tmp+=s[i].nsites;
                        fprintf(myerr,"su[%d]\t%d\t\t%d\t\t(type:%d)\n",i,N-nPlusSvSites[coreNum]-nMinusSvSites[coreNum]-nZeroSvSites[coreNum]+tmp,ctsu[i],(int)lround(s[i].surfType));
                }
        }
}

//dumps in format similar to LAMMPS,
//but with more information, namely
//the type as well as which su for
//each site. for a given site, format
//is:
//site# type    su#     x       y       z
//
//currently supported types are:
//internal      output
//-             0
//+             1
//0             2
//hydrophobic   3
//
//su# is -1 if site is not part of
//a su, and in [0,Nsu) otherwise
//
//file extensions to be canonically
//called .fitrj
void dumpFI(lattice c, su s, int ts, FILE* f, FILE* myerr) {
        extern int N,d,Nsu;
        extern int* L;
        extern int* bc;
        //LAMMPS style header
        fprintf(f,"ITEM: TIMESTEP\n%d\n",ts);
        fprintf(f,"ITEM: NUMBER OF ATOMS\n%d\n",N);
        fprintf(f,"ITEM: BOX BOUNDS ");
        for(int e=0;e<d;++e) {
                if(bc[2*e]==0 || bc[2*e+1]==0) fprintf(f,"ff ");
                else fprintf(f,"pp ");
        }
        if(d==2) fprintf(f,"pp ");
        fprintf(f,"\n");
        for(int e=0;e<d;++e) fprintf(f,"0 %d\n",L[e]);
        if(d==2) fprintf(f,"0 1\n");
        fprintf(f,"ITEM: # type su ");
        if(d==3 || d==2) fprintf(f,"x y z\n");
        else {
                for(int e=0;e<d;++e) fprintf(f,"e%d ",e);
                fprintf(f,"\n");
        }
        //cell/config
        for(int i=0;i<N;++i) {
                fprintf(f,"%d ",i);
                if(c[i].hydrophobic==true) fprintf(f,"3 ");
                else if(*c[i].val<0.0)     fprintf(f,"0 ");
                else if(*c[i].val>0.0)     fprintf(f,"1 ");
                else if(*c[i].val==0.0)    fprintf(f,"2 ");
                else                       fprintf(myerr,"dumpLammps: got unexpected val: %f\n",*c[i].val);
                fprintf(f,"%d ",c[i].su);
                int pos[d];
                getPos(i,pos);
                for(int e=0;e<d;++e) fprintf(f,"%d ",pos[e]);
                if(d==2) fprintf(f,"0 "); //make 2D plane in 3D so VMD can visualize in LAMMPS format
                fprintf(f,"\n");
        }
}
