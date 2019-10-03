//******************
//sections (191001)
//INPUTS
//SOLUTE SETUP
//LATTICE SETUP
//SPOOF?
//SWAP MOVES
//CLUSTER MOVES
//******************

#include "headerBundle.h"

//DEFINING EXTERNAL VAR.S
double pi=0;
double* totSvPlus=NULL;
double* totSvMinus=NULL;
double* totSuPlus=NULL;
double* totSuMinus=NULL;
double plusCharge=1.0;
double minusCharge=-1.0;
int N=0;
int Nsu=0;
int Ntotsusites=0;
int Nc=0;
int d=0;
int shlL=0;
int* nPlusSvSites=NULL;
int* nMinusSvSites=NULL;
int* nZeroSvSites=NULL;
int* nPlusSuSites=NULL;
int* nMinusSuSites=NULL;
int* nHydrSuSites=NULL;
int types=2; //default: -,+; read in from .para
int inshlopt=0;
bool hydrophobicExist=false;
int* L=NULL; 
int* bc=NULL;
int** posArr=NULL;

//DONE DEFINING EXT VAR.S

int main(int argc, char* argv[]) {
        pi=acos(-1);
        //**************************************************************
        //********************** INPUTS ********************************
        //**************************************************************
        //get inputs & assign parameters
        fprintf(stderr,"args:\n");
        for(int i=0;i<argc;++i) fprintf(stderr,"%s ",*(argv+i));
        fprintf(stderr,"\n");
        const int nargs=21;
        if(argc != nargs) {
                fprintf(stderr,"incorrect num. args (%d instead of %d); usage:\n",argc,nargs);
                fprintf(stderr,"./MCIsing verbose?(0/1) spoof?(0/in.lammpstrj) neutralOverallOnSoluteInsertion?(0/1) rotation?(0/1) N L0,L1,L2 sigma parametersIn0,1,.. dataIn0,1,..(or'none') sweepsMult0,1,..(or0.0) +kTeff_su0,1,..(or0.0) enFreq dumpFreq dataFreq suComFreq umbrFreq equiSteps prodSteps rngSeed,rngSeq confDump0,1,..\n");
                return 1;
        }

        fprintf(stderr,"args:\n");
        for(int i=0;i<argc;++i) fprintf(stderr,"%s ",*(argv+i));
        fprintf(stderr,"\n");

        //********************************************************************
        //pre-startup:
        //-start timer
        //-setup pointer tracker for easy freeing at end of run
        //-setup myout & myerr to sidestep need to write to stdout & stderr;
        //     myout & myerr names based on user input dumpout name
        //********************************************************************
        clock_t progStart=time(NULL);

        const int totExtPtrs=12;
        void* extPtrs[totExtPtrs];
        int nExtPtrs=0;

        int nCoresO;
        char* fnames=argv[argc-1];
        char** fname=split(fnames,",",&nCoresO);

        //done w/ pre-startup

        int verbose=atoi(*(++argv));
        const int lSpoofIn=strlen(*(++argv));
        char* spoofIn=malloc(lSpoofIn+1);
        assert(spoofIn!=NULL /*malloc*/);
        extPtrs[nExtPtrs++]=(void*)spoofIn;
        memset(spoofIn,'\0',lSpoofIn+1);
        strncpy(spoofIn,*argv,lSpoofIn);
        const int neutralOverall=atoi(*(++argv));
        bool rotationOpt=(atoi(*(++argv))==0) ? false : true;
        N=atoi(*(++argv)); //extern

        d=3; //extern hardcode; todo: generalize this hardcode
        assert(d==3 /*hardcode*/);

        L=malloc(d*sizeof(*L)); //extern
        bc=malloc(2*d*sizeof(*bc)); //extern
        assert(L!=NULL && bc!=NULL /*malloc*/);
        char* lInputStr=*(++argv);
        char** lIn=NULL;
        int lDim=0;
        lIn=split(lInputStr,",",&lDim);
        for(int e=0;e<d;++e) L[e]=atoi(lIn[e]);
        posArr=setupPos();

        extPtrs[nExtPtrs++]=L;
        extPtrs[nExtPtrs++]=bc;
        extPtrs[nExtPtrs++]=lIn;
        extPtrs[nExtPtrs++]=*posArr;
        extPtrs[nExtPtrs++]=posArr;

        const double sigma=atof(*(++argv));

        //load para&data filenames, su move modifiers: private to each thread
        int nCoresP,nCoresD,nCoresS,nCoresK;
        char* paraInputStr=*(++argv);
        char** parasIn=split(paraInputStr,",",&nCoresP);

        char* dataInputStr=*(++argv);
        char** datasIn=NULL;
        if(strcmp(dataInputStr,"none")!=0) {
                datasIn=split(dataInputStr,",",&nCoresD);
                if(nCoresD!=nCoresP) {
                        fprintf(stderr,"nCoresD!=nCoresP: %d!=%d. exiting\n",nCoresD,nCoresP);
                }
                else {
                        fprintf(stderr,"data inputs:\n");
                        for(int i=0;i<nCoresD;++i) {
                                fprintf(stderr,"%s;",datasIn[i]);
                        }
                        fprintf(stderr,"\n");
                }
        }
        else nCoresD=nCoresP;

        allocGlobalNSiteArrays(&nPlusSvSites,&nMinusSvSites,&nZeroSvSites,
                               &nPlusSuSites,&nMinusSuSites,&nHydrSuSites,
                               &totSvPlus,&totSvMinus,&totSuPlus,&totSuMinus,nCoresD);

        //compute Ewald sum (omp)
        double* ewald=NULL;
        #if FFT_ON
        if(verbose==1) {
                fprintf(stderr,"before ewald\n");
                fflush(stderr);
        }
        ewald=setupEwald(L,sigma,nCoresP);
        assert(ewald!=NULL);
        extPtrs[nExtPtrs++]=ewald;
        if(verbose==1) fprintf(stderr,"ewald sucessfully loaded\n");
        #endif

        char* sweepsMultInputStr=*(++argv);
        char** sweepsMultIn=NULL;
        if(strcmp(sweepsMultInputStr,"0.0")!=0) {
                sweepsMultIn=split(sweepsMultInputStr,",",&nCoresS);
                if(nCoresS!=nCoresP) {
                        fprintf(stderr,"nCoresS!=nCoresP: %d!=%d. exiting\n",nCoresS,nCoresP);
                }
        }
        else nCoresS=nCoresP;

        char* kTeff_suInputStr=*(++argv);
        char** kTeff_suIn=NULL;
        if(strcmp(kTeff_suInputStr,"0.0")!=0) {
                kTeff_suIn=split(kTeff_suInputStr,",",&nCoresK);
                if(nCoresK!=nCoresP) {
                        fprintf(stderr,"nCoresK!=nCoresP: %d!=%d. exiting\n",nCoresK,nCoresP);
                }
        }
        else nCoresK=nCoresP;

        assert(nCoresD==nCoresP && nCoresS==nCoresP && nCoresK==nCoresP);
        fprintf(stderr,"got nCores %d\n",nCoresP);
        extPtrs[nExtPtrs++]=parasIn;
        extPtrs[nExtPtrs++]=datasIn;
        extPtrs[nExtPtrs++]=sweepsMultIn;
        extPtrs[nExtPtrs++]=kTeff_suIn;

        int enFreq=atoi(*(++argv));
        int dumpFreq=atoi(*(++argv));
        int dataFreq=atoi(*(++argv));
        int suComFreq=atoi(*(++argv));
        int umbrFreq=atoi(*(++argv));
        int equiSteps=atoi(*(++argv));
        int prodSteps=atoi(*(++argv));
        //start RNG: get seeds; init
        int extltmp;
        char** rngSeeds=split(*(++argv), ",", &extltmp);
        uint64_t rngseed=strtoul(*rngSeeds,NULL,10);
        uint64_t rngseq=strtoul(*(rngSeeds+1),NULL,10);
        free(*rngSeeds);   free(rngSeeds);    rngSeeds=NULL;
        pcg32_random_t extRandomNumberGenerator;
        pcg32_random_t* extRng=&extRandomNumberGenerator;
        pcg32_srandom_r(extRng,rngseed,rngseq);
        uint32_t threadRngSeeds[2*nCoresP];
        for(int ii=0;ii<2*nCoresP;++ii) threadRngSeeds[ii]=pcg32_random_r(extRng);

	//split into threads; continue loading before launching sim.s on each thread
        #pragma omp parallel for num_threads(nCoresP)
        for(int ii=0;ii<nCoresP;++ii) {
		//thread pre-startup
                const int totIntPtrs=39;
                void* intPtrs[totIntPtrs];
                int nIntPtrs=0;

                const int l=(int)strlen(fname[ii]);
                char* prefix=malloc(l+1);
                assert(prefix!=NULL /*malloc*/);
                intPtrs[nIntPtrs++]=(void*)prefix;
                memset(prefix,'\0',l+1);
                int ltmp;
                char** words=split(fname[ii], ".", &ltmp);
                assert(words!=NULL /*split*/);
                for(int i=0;i<ltmp-1;++i) {
                        strncat(prefix,words[i],strlen(words[i]));
                        if(i!=ltmp-2) strncat(prefix,".",strlen("."));
                }
                bool outputStyleFI;
                if(strcmp(words[ltmp-1],"fitrj")==0)          outputStyleFI=true;
                else if(strcmp(words[ltmp-1],"lammpstrj")==0) outputStyleFI=false;
                else {
                        fprintf(stderr,"unexpected dump file suffix:\n.%s\n"
                                       "please use either .lammpstrj for lammps-style"
                                       "output or .fitrj for FI-style output\n",
                                       words[ltmp-1]);
                        exit(1);
                }
                free(*words);   free(words);    words=NULL;
                char* myoutname=(char*)malloc(strlen(prefix)+4+1);
                char* myerrname=(char*)malloc(strlen(prefix)+4+1);
                char* mvstatsname=(char*)malloc(strlen(prefix)+8+4+1);
                assert(myoutname!=NULL && myerrname!=NULL && mvstatsname!=NULL /*malloc*/);
                intPtrs[nIntPtrs++]=(void*)myoutname;
                intPtrs[nIntPtrs++]=(void*)myerrname;
                intPtrs[nIntPtrs++]=(void*)mvstatsname;
                memset(myoutname,'\0',strlen(prefix)+4+1);
                memset(myerrname,'\0',strlen(prefix)+4+1);
                memset(mvstatsname,'\0',strlen(prefix)+8+4+1);
                strncat(myoutname,prefix,strlen(prefix));
                strncat(myoutname,".out",strlen(".out"));
                strncat(myerrname,prefix,strlen(prefix));
                strncat(myerrname,".err",strlen(".err"));
                strncat(mvstatsname,prefix,strlen(prefix));
                strncat(mvstatsname,"_mvstats.out",strlen("_mvstats.out"));
                fprintf(stderr,"thread %d reporting. fname: %s\nprefix: %s\nmyoutname: %s\nmyerrname: %s\nmvstatsname: %s\n",
                        ii,fname[ii],prefix,myoutname,myerrname,mvstatsname);
                FILE* myout=fopen(myoutname,"w");
                FILE* myerr=fopen(myerrname,"w");
                FILE* mvstats=fopen(mvstatsname,"w");
                assert(myout!=NULL /*fopen*/);
                assert(myerr!=NULL /*fopen*/);
                assert(mvstats!=NULL /*fopen*/);

                fprintf(myerr,"args:\n");
                for(int i=0;i<argc;++i) fprintf(myerr,"%s ",*(argv-(argc-2)+i));
                fprintf(myerr,"\n");

                //done w/ thread pre-startup

                FILE* readIn=fopen(parasIn[ii],"r");        //readIn parameters

                assert(readIn!=NULL /*fopen*/);
                if(verbose==1) {
                        fprintf(myerr,"main: opened parameter in file at %s\n",parasIn[ii]);
                        if(ii==0) {
                                fprintf(stderr,"thread %d: opened parameter in file at %s\n",ii,parasIn[ii]);
                        }
                }
                //setup holders for parameters
                struct parameters parameterHolder;
                para p=&parameterHolder;
                initPara(p); //sets all parameters to 0 by default

                char*** linesByType=NULL;
                int* nLinesByType=readSetParameters(p,&linesByType,ii,readIn,myerr);
                if(p->sigma!=sigma) {
                        fprintf(stderr,"thread %d sigma differs from cmd line arg sigma: %f != %f\n",ii,p->sigma,sigma);
                        exit(1);
                }
                intPtrs[nIntPtrs++]=(void*)linesByType[1]; //hardcode: 1-> solutes
                intPtrs[nIntPtrs++]=(void*)linesByType[2]; //hardcode: 2-> umbrellas
                int isFluidNeutral=willParaLatBeChargeNeut(ii);
                if(isFluidNeutral==0) {
                        fprintf(myerr,"WARNING: pure fluid cannot be neutral with given inputs\n");
                }
                else if(isFluidNeutral==1) {
                        fprintf(myerr,"system is neutral with given inputs\n");
                }
                else exit(1);
                #if FFT_ON      //check runtime and compiletime options agree
                if(p->fftOn==0) {
                        fprintf(myerr,"code built for FFT, but passed p->fftOn=0. Use with FFT or rebuild w FFT_ON switch set to 0\n");
                        exit(1);
                }
                #endif
                #if !FFT_ON
                if(p->fftOn==1) {
                        fprintf(myerr,"code built for non-FFT, but passed p->fftOn=1. Use without FFT or rebuild w FFT_ON switch set to 1\n");
                        exit(1);
                }
                #endif
                #if ISI_ON      //check runtime and compiletime options agree
                if(p->J==0.0) {
                        fprintf(myerr,"code built w Ising on, but passed J=0. can get more efficiency w recompile\n");
                }
                #endif
                #if !ISI_ON
                if(p->J!=0.0) {
                        fprintf(myerr,"code built w Ising off, but passed J!=0. recompile with ISI_ON 1. exiting\n");
                        exit(1);
                }
                #endif
                #if !ISI_ON
                fprintf(myerr,"ISI_ON 0 not yet implemented!\n");
                exit(1);
                #endif
                if(nLinesByType[0]>=1) assert(nLinesByType[1]!=-1 /*readSolutes fail code*/);
                assert(fclose(readIn)==0);
                readIn=NULL;

                int datain,ndatalines;
                char** datalines=NULL;   //used if reading in .data
                if(datasIn!=NULL) {
                        datain=1;
                        readIn=fopen(datasIn[ii],"r");    //readIn .data
                        assert(readIn!=NULL /*fopen*/);
                        if(verbose==1) {
                                fprintf(myerr,"main: opened data in file at %s\n",datasIn[ii]);
                                if(ii==0) {
                                        fprintf(stderr,"thread %d: opened data in file at %s\n",ii,datasIn[ii]);
                                }
                        }
                        ndatalines=readData(&datalines, readIn);
                        if(nLinesByType[0]>=1 && nLinesByType[1]==0) { //Solutes info in data lines
                                nLinesByType[1]=readLines("Solutes\n",&(linesByType[1]),datalines,ndatalines,(FILE*)NULL,myerr);
                        }
                        assert(fclose(readIn)==0);
                }
                else    datain=0;

                double sweepsMult=0.0,kTeff_su=0.0;
                if(sweepsMultIn!=NULL) sweepsMult=atof(sweepsMultIn[ii]);
                if(kTeff_suIn!=NULL) kTeff_su=atof(kTeff_suIn[ii]);

                //**************************************************************
                //********************** SOLUTE SETUP **************************
                //**************************************************************

                //functions expect inputs for solutes, etc.; if no
                //solutes, still pass them a NULL pointer
                su solutes=NULL;
                su trialsolutes=NULL;
                char* shapeNameHolder=NULL;
                double* suComHolder=NULL;
                double* unwrappedSuComHolder=NULL;
                double* orientationHolderH=NULL;
                double** orientationHolder=NULL;
                double* suRelPosBlock=NULL;
                double** suRelPosPtrBlock=NULL;
                double* inShlRelPosBlock=NULL;
                double** inShlRelPosPtrBlock=NULL;
                lattice* suCurrPosBlock=NULL;
                lattice* inShlCurrPosBlock=NULL;
                bool* hydrH=NULL;
                if(Nsu!=0) {
                        loadSu(linesByType[1],nLinesByType[1],&solutes,&trialsolutes,&shapeNameHolder,        \
                                    &suComHolder,&unwrappedSuComHolder,&orientationHolderH,&orientationHolder, \
                                    &suRelPosBlock,&suRelPosPtrBlock,&inShlRelPosBlock,&inShlRelPosPtrBlock,    \
                                    &suCurrPosBlock,&inShlCurrPosBlock,&hydrH,myerr                              );
                        intPtrs[nIntPtrs++]=(void*)solutes;
                        intPtrs[nIntPtrs++]=(void*)trialsolutes;
                        intPtrs[nIntPtrs++]=(void*)shapeNameHolder;
                        intPtrs[nIntPtrs++]=(void*)suComHolder;
                        intPtrs[nIntPtrs++]=(void*)unwrappedSuComHolder;
                        intPtrs[nIntPtrs++]=(void*)orientationHolderH;
                        intPtrs[nIntPtrs++]=(void*)orientationHolder;
                        intPtrs[nIntPtrs++]=(void*)suRelPosBlock;
                        intPtrs[nIntPtrs++]=(void*)suRelPosPtrBlock;
                        intPtrs[nIntPtrs++]=(void*)inShlRelPosBlock;
                        intPtrs[nIntPtrs++]=(void*)inShlRelPosPtrBlock;
                        intPtrs[nIntPtrs++]=(void*)suCurrPosBlock;
                        intPtrs[nIntPtrs++]=(void*)inShlCurrPosBlock;
                        intPtrs[nIntPtrs++]=(void*)hydrH;
                }
                if(ii==0) {
                        fprintf(stderr,"thread %d: loaded su\n",ii);
                }

                //umbrella sampling?
                //setup umbout if umbrFreq != 0
                char* myumboutname=NULL;
                FILE* myumbout=NULL;
                if(umbrFreq!=0) {
                        myumboutname=(char*)malloc(strlen(prefix)+4+2+1);
                        assert(myumboutname!=NULL /*malloc*/);
                        intPtrs[nIntPtrs++]=(void*)myumboutname;
                        memset(myumboutname,'\0',strlen(prefix)+4+2+1);
                        strncat(myumboutname,prefix,strlen(prefix));
                        strncat(myumboutname,"_u.out",strlen("_u.out"));
                        myumbout=fopen(myumboutname,"w");
                        assert(myumbout!=NULL /*fopen*/);
                }
                umb umbr=NULL;
                umb tumbr=NULL;
                struct umbrellas umbrHolder;
                struct umbrellas tumbrHolder;
                if(nLinesByType[0]>=2 && nLinesByType[2]>0) {
                        umbr=&umbrHolder;
                        tumbr=&tumbrHolder;
                        loadUmbrellas(linesByType[2],nLinesByType[2],umbr,solutes,p,myumbout,myerr);
                        loadUmbrellas(linesByType[2],nLinesByType[2],tumbr,trialsolutes,p,myumbout,myerr);
                        intPtrs[nIntPtrs++]=(void*)*umbr->pair;
                        intPtrs[nIntPtrs++]=(void*)umbr->pair;
                        intPtrs[nIntPtrs++]=(void*)umbr->k;
                        intPtrs[nIntPtrs++]=(void*)umbr->r0;
                        intPtrs[nIntPtrs++]=(void*)*tumbr->pair;
                        intPtrs[nIntPtrs++]=(void*)tumbr->pair;
                        intPtrs[nIntPtrs++]=(void*)tumbr->k;
                        intPtrs[nIntPtrs++]=(void*)tumbr->r0;
                }
                free(linesByType);
                free(nLinesByType);
                linesByType=NULL;
                nLinesByType=NULL;
                if(ii==0) {
                        fprintf(stderr,"thread %d: loaded umb\n",ii);
                }

                //setup suComOut if suComFreq != 0
                char* mySuComOutName=NULL;
                FILE* mySuComOut=NULL;
                if(suComFreq!=0) {
                        mySuComOutName=(char*)malloc(strlen(prefix)+4+6+1);
                        assert(mySuComOutName!=NULL /*malloc*/);
                        intPtrs[nIntPtrs++]=(void*)mySuComOutName;
                        memset(mySuComOutName,'\0',strlen(prefix)+4+6+1);
                        strncat(mySuComOutName,prefix,strlen(prefix));
                        strncat(mySuComOutName,"_suCom.out",strlen("_suCom.out"));
                        mySuComOut=fopen(mySuComOutName,"w");
                        assert(mySuComOut!=NULL /*fopen*/);
                }

                if(verbose==1) {
                        fprintf(myerr,"local var.s:\nJ:%f\tQ:%f\tbeta:%f\tsigma:%f\t"        \
                                      "dumpFreq:%d\tequiSteps:%d\tprodSteps:%d\nother param:" \
                                      "N:%d d:%d L,bc: ",p->J,p->Q,p->beta,p->sigma,           \
                                                                dumpFreq,equiSteps,prodSteps,N,d);
                        for(int e=0;e<d;++e) fprintf(myerr,"L[%d]:%d bc[%d]:%d,%d ",e,L[e],e,bc[2*e],bc[2*e+1]);
                        fprintf(myerr,"\n");
                }

                //**************************************************************
                //********************** LATTICE SETUP *************************
                //**************************************************************
                fprintf(myerr,"\nrngseed,seq: %lu,%lu\n",threadRngSeeds[ii],threadRngSeeds[ii+nCoresP]);
                pcg32_random_t randomNumberGenerator;
                pcg32_random_t* rng=&randomNumberGenerator;
                pcg32_srandom_r(rng,threadRngSeeds[ii],threadRngSeeds[ii+nCoresP]);
                fprintf(myerr,"before valBlock\n");
                fflush(myerr);
                //initialize lattice
                lattice lat=malloc(N*sizeof(*lat));
                double* valBlock=fftw_malloc(2*N*sizeof(*valBlock));
                assert(lat!=NULL && valBlock!=NULL /*malloc*/);
                intPtrs[nIntPtrs++]=(void*)lat;
                for(int i=0;i<N;++i) lat[i].val=valBlock+i;
                setLatInitSu(lat);
                if(verbose==1) fprintf(myerr,"got past lat mallocs\n");
                fprintf(myerr,"before setLatNeigh\n");
                fflush(myerr);
                if(ii==0) {
                        fprintf(stderr,"thread %d: before setting neighbors\n",ii);
                }
                void** tmpptrs=setLatNeigh_coul_twoway(lat,p,myerr);
                assert(tmpptrs!=NULL /*setLatNeigh_coul_twoway fail*/);
                for(int i=0;i<2;++i) intPtrs[nIntPtrs++]=tmpptrs[i];
                free(tmpptrs);  tmpptrs=NULL;
                if(ii==0) {
                        fprintf(stderr,"thread %d: finished setting neighbors\n",ii);
                }

                double* chg=NULL;
                if(datain==1) {
                        chg=setLatVal_data(lat,datalines,ndatalines,ii,myerr);
                        free(*datalines);   free(datalines);    datalines=NULL;
                        fprintf(myerr,"\n\nloaded config:\n");
                        if(ii==0) {
                                fprintf(stderr,"\n\nthread %d: loaded config\n",ii);
                        }
                }
                else {
                        chg=setLatVal_rand(lat,solutes,ii,rng,myerr);
                        fprintf(myerr,"\n\nbuilt config:\n");
                        if(ii==0) {
                                fprintf(stderr,"\n\nthread %d: built config\n",ii);
                        }
                }
                free(chg);
                //if(outputStyleFI==true) dumpFI(lat,solutes,0,myerr,myerr);
                //else                    dumpLammps(lat,solutes,0,myerr,myerr);
                //fprintf(myerr,"\n\n");

                //if(strcmp(spoofIn,"0")==0) {
                //        if(datain==1) {
                //                double* chg=setLatVal_data(lat, datalines, ndatalines,myerr);
                //                free(*datalines);   free(datalines);    datalines=NULL;
                //                free(chg);
                //                fprintf(myerr,"\n\nloaded config:\n");
                //                if(outputStyleFI==true) dumpFI(lat,solutes,0,myerr,myerr);
                //                else                    dumpLammps(lat,solutes,0,myerr,myerr);
                //                fprintf(myerr,"\n\n");
                //        }
                //        else {
                //                double* chg=setLatVal_rand(lat,rng,myerr);
                //                free(chg);
                //        }
                //}

                lattice triallat=malloc(N*sizeof(*triallat));
                assert(triallat!=NULL /*malloc*/);
                intPtrs[nIntPtrs++]=(void*)triallat;
                for(int i=0;i<N;++i) {
                        triallat[i].val=valBlock+i+N;
                        *triallat[i].val=*lat[i].val;
                }
                setLatInitSu(triallat); //set lat su values to -1: no su
                tmpptrs=setLatNeigh_coul_twoway(triallat,p,myerr);
                assert(tmpptrs != NULL /*setLatNeigh_coul_twoway fail*/);
                for(int i=0;i<2;++i) intPtrs[nIntPtrs++]=tmpptrs[i];
                if(Nsu!=0) {
                        const uint64_t rngStateAtSuRelPosStart=rng->state;
                        const uint64_t rngIncAtSuRelPosStart=rng->inc; //ensure c,tc same
                        int setSuRelPosRet;
                        setSuRelPosRet=setSuRelPos(lat,solutes,inshlopt,neutralOverall,ii,rng,myerr);
                        assert(setSuRelPosRet==0);
                        if(verbose==1) fprintf(myerr,"setSuRelPosRet for c: %d\n",setSuRelPosRet);
                        rng->state=rngStateAtSuRelPosStart;
                        rng->inc=rngIncAtSuRelPosStart; //ensure c,tc same
                        setSuRelPosRet=setSuRelPos(triallat,trialsolutes,inshlopt,neutralOverall,ii,rng,myerr);
                        assert(setSuRelPosRet==0);
                        if(verbose==1) {
                                fprintf(myerr,"setSuRelPosRet for tc: %d\n",setSuRelPosRet);
                                fflush(myerr);
                        }
                }
                else {
                        for(int i=0;i<N;++i) {
                                lat[i].hydrophobic=false;
                                triallat[i].hydrophobic=false;
                        }
                }
                hydrophobicExist=setHydrophobicExist(solutes);
                if(verbose==1) {
                        fprintf(myerr,"hydrophobic sites exist on lattice: ");
                        if(hydrophobicExist==true) fprintf(myerr,"true\n");
                        else                       fprintf(myerr,"false\n");
                }
                if(strcmp(spoofIn,"0")==0) {
                        if(verbose==1) {
                                const int sumSv=nPlusSvSites[ii]+nMinusSvSites[ii]+nZeroSvSites[ii];
                                const double totSvChg=totSvPlus[ii]+totSvMinus[ii];
                                const double totSuChg=totSuPlus[ii]+totSuMinus[ii];
                                fprintf(myerr,"main: init lat val.s set on lattice.\n"
                                              "nPlusSites,Minus,Zero: %d,%d,%d of %d solvent, %d sites.\n"
                                              "\tor %f,%f,%f per solvent\n"
                                              "\tor %f,%f,%f per SITE\n"
                                              "totSvChg,totSuChg,sum: %.1f,%.1f,%.1f\n",
                                        nPlusSvSites[ii],nMinusSvSites[ii],nZeroSvSites[ii],sumSv,N,(double)(nPlusSvSites[ii])/sumSv,
                                        (double)(nMinusSvSites[ii])/sumSv,(double)(nZeroSvSites[ii])/sumSv,(double)(nPlusSvSites[ii])/N,
                                        (double)(nMinusSvSites[ii])/N,(double)(nZeroSvSites[ii])/N,totSvChg,totSuChg,totSvChg+totSuChg);
                                fprintf(myout,"####################################\n");
                                fflush(myerr);
                        }
                        setCharge(lat,ii);
                        if(verbose==1) {
                                const int sumSv=nPlusSvSites[ii]+nMinusSvSites[ii]+nZeroSvSites[ii];
                                const double totSvChg=totSvPlus[ii]+totSvMinus[ii];
                                const double totSuChg=totSuPlus[ii]+totSuMinus[ii];
                                fprintf(myerr,"main: init lat val.s set on lattice.\n"
                                              "nPlusSites,Minus,Zero: %d,%d,%d of %d solvent, %d sites.\n"
                                              "\tor %f,%f,%f per solvent\n"
                                              "\tor %f,%f,%f per SITE\n"
                                              "totSvChg,totSuChg,sum: %.1f,%.1f,%.1f\n",
                                        nPlusSvSites[ii],nMinusSvSites[ii],nZeroSvSites[ii],sumSv,N,(double)(nPlusSvSites[ii])/sumSv,
                                        (double)(nMinusSvSites[ii])/sumSv,(double)(nZeroSvSites[ii])/sumSv,(double)(nPlusSvSites[ii])/N,
                                        (double)(nMinusSvSites[ii])/N,(double)(nZeroSvSites[ii])/N,totSvChg,totSuChg,totSvChg+totSuChg);

                                fprintf(myout,"####################################\n");
                                fflush(myerr);
                        }
                }
                fprintf(myerr,"post-su placement\n");
                if(outputStyleFI==true) dumpFI(lat,solutes,0,myerr,myerr);
                else                    dumpLammps(lat,solutes,0,ii,myerr,myerr);
                fprintf(myerr,"\n\n");

                ////set rng to same seed for testing
                ////COMMENT THIS OUT FOR NORMAL RUNS
                //const int testSeed=1;
                //const int testSeq=2;
                //pcg32_srandom_r(rng,testSeed,testSeq);

                //**************************************************************
                //*********************** SPOOF? *******************************
                //**************************************************************
                //"spoof" is an option to, rather than running a fresh simulation,
                //read in the trajectory info. from a previously run simulation
                //and output the energies of those config.s according to the
                //parameters input rather than the parameters according to which
                //they were generated. This is an entirely deterministic process
                //and should give the same energies for the same parameter sets.
                if(strcmp(spoofIn,"0")!=0) { //compute E.s from config.s @ spoofIn
                        fprintf(myerr,"entering spoof mode using trajectory at %s\n",spoofIn);
                        FILE* trjIn=fopen(spoofIn,"r");
                        assert(trjIn!=NULL /*fopen*/);
                        runSpoof(lat,solutes,ewald,p,umbr,trjIn,myout);
                        assert(fclose(trjIn)==0);
                        freeAll(intPtrs,nIntPtrs);
                        fprintf(myerr,"spoof complete. exiting\n");
                        continue;
                }

                //**************************************************************
                //********************** SWAP MOVES ****************************
                //**************************************************************
                
                double energyHolder[ENERGY_LENGTH]; //hardcode
                double* en=energyHolder;
                energy_full(lat,solutes,en,ewald,p);

                fprintf(myerr,"pre-run energy\n");
                const double eUmb=energy_umbr(umbr);
                const double eSelf=-p->e_self;
                fprintf(myerr,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",        \
                        0,en[0]+en[1]+en[2]+en[3]+en[4]+en[5]+en[6]+en[7]+en[8]+en[9]+eUmb+eSelf, \
                        en[0],en[1],en[2],en[3],en[4],en[5],en[6],en[7],en[8],en[9],eUmb,eSelf     );

                const int mvtypes=2*3;  //hardcode
                int nmoves[mvtypes]; //hardcode
                for(int i=0;i<mvtypes;++i) nmoves[i]=0;
                int* mvtmp=NULL;
                time_t equiStart = time(NULL);
                FILE* dumpout=NULL;

                //skip unnecessary code if equisteps == 0
                if(equiSteps > 0) {
                        fprintf(myerr,"fname: %s\tlen: %d\n",fname[ii],l);
                        char* feq=(char*)malloc(l+6);
                        assert(feq!=NULL /*malloc*/);
                        memset(feq,'\0',l+6);

                        strncat(feq,prefix,strlen(prefix)); //found @ start of main()
                        if(verbose==1) fprintf(myerr,"feq after cating prefix: %s\n",feq);
                        if(outputStyleFI==false) {
                                strncat(feq,"_swap.lammpstrj",strlen("_swap.lammpstrj"));
                                if(verbose==1) fprintf(myerr,"feq after cating '_swap.lammpstrj': %s\n",feq);
                        }
                        else {
                                strncat(feq,"_swap.fitrj",strlen("_swap.fitrj"));
                                if(verbose==1) fprintf(myerr,"feq after cating '_swap.fitrj': %s\n",feq);
                        }



                        if(dumpout==NULL) dumpout=fopen(feq,"w");
                        assert(dumpout!=NULL /*fopen*/);
                        if(verbose==1)     fprintf(myerr,"main: opened dump file at %s\n",feq);
                        free(feq);      feq=NULL;
                        fprintf(myout,"##################equi########################\n" \
                                      "#step\tE_tot\tE_isvsv\tE_isvsu\tE_isusu\tE_csrsvsv\tE_csrsvsu\tE_csrsusu\tE_clr\tE_ihh\tE_isvh\tE_isuh\tE_umb\tE_self\n");
                        for(int i=0;i<equiSteps;++i) {
                                #if EQUI_INNER_DUMP
                                mvtmp=MCstep_flipswap(lat,triallat,en,solutes,trialsolutes,umbr,tumbr,ewald,inshlopt,rotationOpt,sweepsMult,\
                                                      kTeff_su,p,ii,rng,enFreq,dumpFreq,suComFreq,umbrFreq,dataFreq,dumpout,prefix,i,\
                                                      outputStyleFI,myout,myumbout,mySuComOut,myerr);
                                #endif
                                #if !EQUI_INNER_DUMP
                                mvtmp=MCstep_flipswap(lat,triallat,en,solutes,trialsolutes,umbr,tumbr,ewald,inshlopt,rotationOpt,sweepsMult,\
                                                      kTeff_su,p,ii,rng,myout,myumbout,mySuComOut,myerr);
                                #endif
                                for(int j=0;j<mvtypes;++j) nmoves[j]+=mvtmp[j];
                                free(mvtmp);    mvtmp=NULL;
                                #if !EQUI_INNER_DUMP
                                if(i%enFreq==0) {
                                        const double eUmb=energy_umbr(umbr);
                                        const double eSelf=-p->e_self;
                                        //energy_full(lat,solutes,en,ewald,p);
                                        fprintf(myout,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",        \
                                                i,en[0]+en[1]+en[2]+en[3]+en[4]+en[5]+en[6]+en[7]+en[8]+en[9]+eUmb+eSelf, \
                                                en[0],en[1],en[2],en[3],en[4],en[5],en[6],en[7],en[8],en[9],eUmb,eSelf     );
                                }
                                if(dumpFreq!=0 && i%dumpFreq==0) {
                                        if(outputStyleFI==true) dumpFI(lat,solutes,i,dumpout,myerr);
                                        else                    dumpLammps(lat,solutes,i,ii,dumpout,myerr);
                                }
                                if(solutes!=NULL && suComFreq!=0 && i%suComFreq==0) dumpAllSuComs(solutes,mySuComOut);
                                if(umbr!=NULL && umbrFreq!=0 && i%umbrFreq==0) {
                                        fprintf(myumbout,"%d\t",i);
                                        for(int j=0;j<(umbr->npairs);++j) {
                                                const double r_scal=dist_vect((umbr->pair[j][0])->com,(umbr->pair[j][1])->com);
                                                fprintf(myumbout,"%f\t",r_scal);
                                        }
                                        fprintf(myumbout,"\n");
                                }
                                if(dataFreq!=0 && i%dataFreq==0) {
                                        const uint16_t lprefix=(uint16_t)strlen(prefix); //strlen ret pos of '\0'
                                        const uint16_t ldotdata=(uint16_t)strlen(".data");
                                        const uint16_t ndig=ndigits((long int)i);
                                        char* datafstr=(char*)malloc(lprefix+1u+ndig+ldotdata);
                                        strncpy(datafstr,prefix,lprefix);
                                        const int snpfret=snprintf(datafstr+lprefix,ndig+2u,"_%u",i);
                                        assert(snpfret==ndig+1 /*snprintf*/);
                                        strncat(datafstr,".data",ldotdata);
                                        FILE* dataout=fopen(datafstr,"w");
                                        assert(dataout!=NULL /*fopen*/);
                                        dumpData(lat,solutes,dataout);
                                        assert(fclose(dataout)==0);
                                        free(datafstr); datafstr=NULL;
                                }
                                #endif
                        }
                        fprintf(myout,"Move type\t\tAccepted\t\tAttempted\t\tFraction accepted\n");
                        fprintf(myout,"Solvent flip/swap\t%d\t\t%d\t\t%f\n",nmoves[0],nmoves[1],(double)nmoves[0]/nmoves[1]);
                        fprintf(myout,"Solvent cluster  \t%d\t\t%d\t\t%f\n",nmoves[4],nmoves[5],(double)nmoves[4]/nmoves[5]);
                        if(Nsu!=0) fprintf(myout,"Solute  \t\t%d\t\t%d\t\t%f\n",nmoves[2],nmoves[3],(double)nmoves[2]/nmoves[3]);

                        assert(fclose(dumpout)==0);
                        dumpout=NULL;
                }
                //**************************************************************
                //*********************** CLUSTER MOVES ************************
                //**************************************************************

                time_t prodStart = time(NULL);
                if(prodSteps>0) {
                        if(nZeroSvSites[ii]>0) {
                                fprintf(myerr,"cluster moves not currently implemented for lattice with defects.\n"
                                              "please either:\n"
                                              "-set nZeroSvSites=0 & run with cluster moves\n"
                                              "-set nZeroSvSites>0 & run with ONLY swap moves\n"
                                              "-implement (or convince Nick to implement) cluster moves with defects\n"
                                              "ending run.\n");
                                exit(1);
                        }
                        fprintf(myout,"##################prod########################\n" \
                                      "#step\tE_tot\tE_isvsv\tE_isvsu\tE_isusu\tE_csrsvsv\tE_csrsvsu\tE_csrsusu\tE_clr\tE_ihh\tE_isvh\tE_isuh\tE_umb\tE_self\n");
                        if(dumpout==NULL) dumpout=fopen(fname[ii],"w");
                        assert(dumpout!=NULL /*fopen*/);
                        if(verbose==1) {
                                fprintf(myerr,"main: opened dump file at %s\n"\
                                              "starting prod. run\n",fname[ii]);
                        }
                        for(int j=0;j<mvtypes;++j) nmoves[j]=0;
                        for(int i=0;i<prodSteps;++i) {
                                mvtmp=MCstep_cluster(lat,triallat,en,solutes,trialsolutes,umbr,tumbr,ewald,\
                                                     inshlopt,rotationOpt,sweepsMult,kTeff_su,p,ii,rng,mvstats,myerr);
                                for(int j=0;j<mvtypes;++j) nmoves[j]+=mvtmp[j]; // hardcode
                                free(mvtmp);    mvtmp=NULL;
                                if(i%enFreq==0) {
                                        const double eUmb=energy_umbr(umbr);
                                        const double eSelf=-p->e_self;
                                        //energy_full(lat,solutes,en,ewald,p);
                                        fprintf(myout,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",        \
                                                i,en[0]+en[1]+en[2]+en[3]+en[4]+en[5]+en[6]+en[7]+en[8]+en[9]+eUmb+eSelf, \
                                                en[0],en[1],en[2],en[3],en[4],en[5],en[6],en[7],en[8],en[9],eUmb,eSelf     );
                                }
                                if(dumpFreq!=0 && i%dumpFreq==0) {
                                        if(outputStyleFI==true) dumpFI(lat,solutes,i,dumpout,myerr);
                                        else                    dumpLammps(lat,solutes,i,ii,dumpout,myerr);
                                }
                                if(solutes!=NULL && suComFreq!=0 && i%suComFreq==0) dumpAllSuComs(solutes,mySuComOut);
                                if(umbr!=NULL && umbrFreq!=0 && i%umbrFreq==0) {
                                        fprintf(myumbout,"%d\t",i);
                                        for(int j=0;j<(umbr->npairs);++j) {
                                                const double r_scal=dist_vect((umbr->pair[j][0])->com,(umbr->pair[j][1])->com);
                                                fprintf(myumbout,"%f\t",r_scal);
                                        }
                                        fprintf(myumbout,"\n");
                                }
                                if(dataFreq!=0 && i%dataFreq==0) {
                                        const uint16_t lprefix=(uint16_t)strlen(prefix); //strlen ret pos of '\0'
                                        const uint16_t ldotdata=(uint16_t)strlen(".data");
                                        const uint16_t ndig=ndigits((long int)(i+equiSteps));
                                        char* datafstr=malloc((lprefix+ndig+ldotdata+2u)*sizeof(*datafstr));
                                        assert(datafstr!=NULL /*malloc*/);
                                        const int snpfret=snprintf(datafstr,lprefix+ndig+ldotdata+2u,"%s_%u.data",prefix,(i+equiSteps));
                                        FILE* dataout=fopen(datafstr,"w");
                                        assert(dataout!=NULL /*fopen*/);
                                        dumpData(lat,solutes,dataout);
                                        assert(fclose(dataout)==0);
                                        fprintf(stderr,"datafstr: %s\n",datafstr);
                                        fflush(stderr);
                                        free(datafstr);
                                }
                        }
                }
                fprintf(myout,"Move type\t\tAccepted\t\tAttempted\t\tFraction accepted\n");
                fprintf(myout,"Solvent flip/swap\t%d\t\t%d\t\t%f\n",nmoves[0],nmoves[1],(double)nmoves[0]/nmoves[1]);
                fprintf(myout,"Solvent cluster\t%d\t\t%d\t\t%f\n",nmoves[4],nmoves[5],(double)nmoves[4]/nmoves[5]);
                if(Nsu!=0) fprintf(myout,"Solute\t\t%d\t\t%d\t\t%f\n",nmoves[2],nmoves[3],(double)nmoves[2]/(double)nmoves[3]);

                if(verbose==1) {
                        time_t end = time(NULL);
                        double totTime = difftime(end,progStart);
                        double setupTime = difftime(equiStart,progStart);
                        double equiTime = difftime(prodStart,equiStart);
                        double prodTime = difftime(end,prodStart);
                        fprintf(myout,"total CPU time: %f s = %f mins = %f hrs = %f days\n",totTime,totTime/60.0,totTime/3600.0,totTime/3600.0/24.0);
                        fprintf(myout,"setup CPU time: %f s = %f mins = %f hrs = %f days\n",setupTime,setupTime/60.0,setupTime/3600.0,setupTime/3600.0/24.0);
                        fprintf(myout,"equilibration run CPU time: %f s = %f mins = %f hrs = %f days; that is ~ %f ms per MC cycle/step\n",equiTime,equiTime/60.0,equiTime/3600.0,equiTime/3600.0/24.0,equiTime/equiSteps*1000.0);
                        fprintf(myout,"production run CPU time: %f s = %f mins = %f hrs = %f days; that is ~ %f ms per MC cycle/step\n",prodTime,prodTime/60.0,prodTime/3600.0,prodTime/3600.0/24.0,prodTime/prodSteps*1000.0);
                }

                //cleanup & final output
                if(myout!=NULL) {
                        const int filecloseout=fclose(myout);
                        assert(filecloseout==0 /*fclose*/);
                        myout=NULL;
                }
                if(myerr!=NULL) {
                        const int filecloseout=fclose(myerr);
                        assert(filecloseout==0 /*fclose*/);
                        myerr=NULL;
                }
                if(dumpout!=NULL) {
                        const int filecloseout=fclose(dumpout);
                        assert(filecloseout==0 /*fclose*/);
                        dumpout=NULL;
                }
                if(mySuComOut!=NULL) {
                        const int filecloseout=fclose(mySuComOut);
                        assert(filecloseout==0 /*fclose*/);
                        mySuComOut=NULL;
                }
                if(myumbout!=NULL) {
                        const int filecloseout=fclose(myumbout);
                        assert(filecloseout==0 /*fclose*/);
                        myumbout=NULL;
                }
                if(mvstats!=NULL) {
                        const int filecloseout=fclose(mvstats);
                        assert(filecloseout==0 /*fclose*/);
                        mvstats=NULL;
                }

                char* paraname=(char*)malloc(strlen(prefix)+5+1);
                char* dataname=(char*)malloc(strlen(prefix)+5+1);
                assert(paraname!=NULL && dataname!=NULL /*malloc*/);
                intPtrs[nIntPtrs++]=(void*)paraname;
                intPtrs[nIntPtrs++]=(void*)dataname;
                memset(paraname,'\0',strlen(prefix)+5+1);
                memset(dataname,'\0',strlen(prefix)+5+1);
                strncat(paraname,prefix,strlen(prefix));
                strncat(paraname,".para",strlen(".para"));
                strncat(dataname,prefix,strlen(prefix));
                strncat(dataname,".data",strlen(".data"));

                dumpout=fopen(paraname,"w");
                assert(dumpout!=NULL /*fopen*/);
                dumpParameters(p,solutes,umbr,ii,dumpout,myerr);
                assert(fclose(dumpout)==0);
                dumpout=fopen(dataname,"w");
                assert(dumpout!=NULL /*fopen*/);
                dumpData(lat,solutes,dumpout);
                assert(fclose(dumpout)==0);

                freeAll(intPtrs,nIntPtrs);
                fftw_free(valBlock);
        }
        freeAll(extPtrs,nExtPtrs);

        return 0;
}
