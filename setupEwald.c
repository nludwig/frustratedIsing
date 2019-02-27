#include "setupEwald.h"
//compile w/ icc for ~50x speedup

void getPos(int n, int* pos, int* L) {
        const int d=3;
        int prod=1;
        n+=1; //req.d to work
        for(int e=0;e<d;++e) prod*=L[e];
        for(int e=0;e<d;++e) {
                prod/=L[e];
                int i=0;
                while(i*prod < n) ++i;
                --i;    //go one too far
                pos[e]=i;
                n-=i*prod;
        }
}

double pbc(double dst, int e, int* L) {
        while(dst<-L[e]/2) dst+=L[e];
        while(dst>=L[e]/2) dst-=L[e];
        return dst;
}

//dr must have d entries
void diff_lat(double* dr, unsigned long i, unsigned long j, int* L) {
        const int d=3;
        int u[d];
        int v[d];
        getPos(i,u,L);
        getPos(j,v,L);
        for(int e=0;e<d;++e) {
                dr[e]=(double)(u[e]-v[e]);
                dr[e]=pbc(dr[e],e,L);
        }
}

//use data N(N-1)/2 + 1 , cos's, not e^ikx's so only sum over +k's
void setupEwald(double*** ewald, int* L, para p) {
        const int d=3;
        assert(d==3 /*setupEwald assumes...*/);
        const int N=L[0]*L[1]*L[2];
        const double pi=acos(-1);
        const double s2=p->sigma * p->sigma;
        double* ewaldHolder=malloc((N*(N-1)/2+1)*sizeof(*ewaldHolder));
        *ewald=malloc(((N-1)+1)*sizeof(**ewald));
        double complex* cexpValsHolder=malloc(L[0]*L[0]*sizeof(*cexpValsHolder)); //add precompute of vectors
        double complex** cexpVals=malloc(L[0]*sizeof(*cexpVals));
        assert(ewaldHolder!=NULL && *ewald!=NULL && cexpValsHolder!=NULL && cexpVals!=NULL /*malloc*/);
        int i=0;
        for(int k=0;k<L[0];++k) {
                cexpVals[k]=cexpValsHolder+i;
                for(int x=-L[0]/2;x<L[0]/2;++x) {
                        cexpVals[k][x+L[0]/2]=cexp(2.0*pi*I/L[0]*k*x);
                        ++i;
                }
        }
        double cx,cy,cz;
        i=0;
        for(int j=0;j<N-1;++j) {
                (*ewald)[j]=ewaldHolder+i;
                for(int l=j+1;l<N;++l) {
                        double rjl[d];
                        diff_lat(rjl,j,l,L);

                        double sum=0.0;
                        for(int kx=0;kx<L[0];++kx) {
                                const double kx2=(double)(kx*kx)/(L[0]*L[0]);
                                const double xkx=rjl[0]*kx/L[0];
                                if(kx!=0) cx=1.0;
                                else      cx=0.5;
                                for(int ky=0;ky<L[1];++ky) {
                                        const double ky2=(double)(ky*ky)/(L[1]*L[1]);
                                        const double yky=rjl[1]*ky/L[1];
                                        if(ky!=0) cy=1.0;
                                        else      cy=0.5;
                                        for(int kz=0;kz<L[2];++kz) {
                                                if(kx!=0 || ky!=0 || kz!=0) {
                                                        const double kz2=(double)(kz*kz)/(L[2]*L[2]);
                                                        const double zkz=rjl[2]*kz/L[2];
                                                        if(kz!=0) cz=1.0;
                                                        else      cz=0.5;

                                                        const double k2=4.0*pi*pi*(kx2+ky2+kz2);
                                                        //use binary numbers to repr permutations
                                                        //0-> +1
                                                        //1-> -1
                                                        const int x=(int)rjl[0];
                                                        const int y=(int)rjl[1];
                                                        const int z=(int)rjl[2];
                                                        const double complex ekr000=cexpVals[kx][x+L[0]/2]*
                                                                cexpVals[ky][y+L[0]/2]*cexpVals[kz][z+L[0]/2];
                                                        const double complex ekr001=cexpVals[kx][x+L[0]/2]*
                                                                cexpVals[ky][y+L[0]/2]*conj(cexpVals[kz][z+L[0]/2]);
                                                        const double complex ekr010=cexpVals[kx][x+L[0]/2]*
                                                                conj(cexpVals[ky][y+L[0]/2])*cexpVals[kz][z+L[0]/2];
                                                        const double complex ekr011=cexpVals[kx][x+L[0]/2]*
                                                                conj(cexpVals[ky][y+L[0]/2])*conj(cexpVals[kz][z+L[0]/2]);

                                                        //const double kr000=2.0*pi*(xkx+yky+zkz);
                                                        //const double kr001=2.0*pi*(xkx+yky-zkz);
                                                        //const double kr010=2.0*pi*(xkx-yky+zkz);
                                                        //const double kr011=2.0*pi*(xkx-yky-zkz);
                                                        //fprintf(myerr,"kx,ky,kz %d,%d,%d\tx,y,z %d,%d,%d\tkr000 %f\tekr000:%f+I%f\tcexp(Ikr):%f+I%f\n",kx,ky,kz,x,y,z,kr000,creal(ekr000),cimag(ekr000),creal(cexp(I*kr000)),cimag(cexp(I*kr000)));
                                                        //fflush(myerr);
                                                        const double Ek=exp(-k2*s2/4.0)/k2*creal(ekr000+ekr001+ekr010+ekr011);
                                                        ////use binary numbers to repr permutations
                                                        ////0-> +1
                                                        ////1-> -1
                                                        //const double kr000=2.0*pi*(xkx+yky+zkz);
                                                        //const double kr001=2.0*pi*(xkx+yky-zkz);
                                                        //const double kr010=2.0*pi*(xkx-yky+zkz);
                                                        //const double kr011=2.0*pi*(xkx-yky-zkz);
                                                        //const double Ek=exp(-k2*s2/4.0)/k2*
                                                        //        creal(cexp(I*kr000)+cexp(I*kr001)+cexp(I*kr010)+cexp(I*kr011));
                                                        sum+=2.0*cx*cy*cz*Ek;
                                                }
                                        }
                                }
                        }
                        const int lmj=l-j-1;
                        (*ewald)[j][lmj]=sum;
                        ++i;
                }
        }

        //diag element: 'self E'
        double rjl[d];
        rjl[0]=0.0;
        rjl[1]=0.0;
        rjl[2]=0.0;
        double sum=0.0;
        for(int kx=0;kx<L[0];++kx) {
                const double kx2=(double)(kx*kx)/(L[0]*L[0]);
                const double xkx=rjl[0]*kx/L[0];
                if(kx!=0) cx=1.0;
                else      cx=0.5;
                for(int ky=0;ky<L[1];++ky) {
                        const double ky2=(double)(ky*ky)/(L[1]*L[1]);
                        const double yky=rjl[1]*ky/L[1];
                        if(ky!=0) cy=1.0;
                        else      cy=0.5;
                        for(int kz=0;kz<L[2];++kz) {
                                if(kx!=0 || ky!=0 || kz!=0) {
                                        const double kz2=(double)(kz*kz)/(L[2]*L[2]);
                                        const double zkz=rjl[2]*kz/L[2];
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
                                        sum+=2.0*cx*cy*cz*Ek;
                                }
                        }
                }
        }

        (*ewald)[N-1]=ewaldHolder+i;
        (*ewald)[N-1][0]=sum;
        free(cexpValsHolder); free(cexpVals);
}

////use data N(N-1)/2 + 1 , cos's, not e^ikx's so only sum over +k's
//void setupEwald_omp(double*** ewald, int* L, para p) {
//        fprintf(stderr,"setupEwald_omp: fn called\n");
//        fflush(stderr);
//        const int d=3;
//        assert(d==3 /*setupEwald assumes...*/);
//        const int N=L[0]*L[1]*L[2];
//        //const double err=10e-10; //hardcode
//        const double pi=acos(-1);
//        const double s2=p->sigma * p->sigma;
//        double* ewaldHolder=malloc(((size_t)N*((size_t)N-(size_t)1)/(size_t)2+(size_t)1)*sizeof(*ewaldHolder)); //ensure no wrapping
//        assert(ewaldHolder!=NULL);
//        *ewald=malloc(((N-1)+1)*sizeof(**ewald));
//        assert(*ewald!=NULL);
//        double complex* cexpValsHolderX=malloc(L[0]*L[0]*sizeof(*cexpValsHolderX));
//        assert(cexpValsHolderX!=NULL);
//        double complex* cexpValsHolderY=malloc(L[1]*L[1]*sizeof(*cexpValsHolderY));
//        assert(cexpValsHolderY!=NULL);
//        double complex* cexpValsHolderZ=malloc(L[2]*L[2]*sizeof(*cexpValsHolderZ));
//        assert(cexpValsHolderZ!=NULL);
//        double complex** cexpValsX=malloc(L[0]*sizeof(*cexpValsX)); //add precompute of vectors
//        assert(cexpValsX!=NULL);
//        double complex** cexpValsY=malloc(L[1]*sizeof(*cexpValsY));
//        assert(cexpValsY!=NULL);
//        double complex** cexpValsZ=malloc(L[2]*sizeof(*cexpValsZ));
//        assert(cexpValsZ!=NULL);
//        //assert(ewaldHolder!=NULL && *ewald!=NULL && 
//        //       cexpValsHolderX!=NULL && cexpValsX!=NULL &&
//        //       cexpValsHolderY!=NULL && cexpValsY!=NULL &&
//        //       cexpValsHolderZ!=NULL && cexpValsZ!=NULL /*malloc*/);
//        int i=0;
//        for(int k=0;k<L[0];++k) {
//                cexpValsX[k]=cexpValsHolderX+i;
//                for(int x=-L[0]/2;x<L[0]/2;++x) {
//                        cexpValsX[k][x+L[0]/2]=cexp(2.0*pi*I/L[0]*k*x);
//                        ++i;
//                }
//        }
//        i=0;
//        for(int k=0;k<L[1];++k) {
//                cexpValsY[k]=cexpValsHolderY+i;
//                for(int y=-L[1]/2;y<L[1]/2;++y) {
//                        cexpValsY[k][y+L[1]/2]=cexp(2.0*pi*I/L[1]*k*y);
//                        ++i;
//                }
//        }
//        i=0;
//        for(int k=0;k<L[2];++k) {
//                cexpValsZ[k]=cexpValsHolderZ+i;
//                for(int z=-L[2]/2;z<L[2]/2;++z) {
//                        cexpValsZ[k][z+L[2]/2]=cexp(2.0*pi*I/L[2]*k*z);
//                        ++i;
//                }
//        }
//        i=0;
//        for(int j=0;j<N-1;++j) {
//                (*ewald)[j]=ewaldHolder+i;
//                i+=N-j-1;
//        }
//        (*ewald)[N-1]=ewaldHolder+i;
//        fprintf(stderr,"setupEwald_omp: pre-omp loop\n");
//        fflush(stderr);
//        #pragma omp parallel for schedule(dynamic)
//        for(int j=0;j<N-1;++j) {
//                for(int l=j+1;l<N;++l) {
//                        double rjl[d];
//                        diff_lat(rjl,j,l,L);
//                        double sum=0.0;
//                        double cx,cy,cz;
//                        for(int kx=0;kx<L[0];++kx) {
//                                const double kx2=(double)(kx*kx)/(L[0]*L[0]);
//                                if(kx!=0) cx=1.0;
//                                else      cx=0.5;
//                                for(int ky=0;ky<L[1];++ky) {
//                                        const double ky2=(double)(ky*ky)/(L[1]*L[1]);
//                                        if(ky!=0) cy=1.0;
//                                        else      cy=0.5;
//                                        for(int kz=0;kz<L[2];++kz) {
//                                                if(kx!=0 || ky!=0 || kz!=0) {
//                                                        const double kz2=(double)(kz*kz)/(L[2]*L[2]);
//                                                        if(kz!=0) cz=1.0;
//                                                        else      cz=0.5;
//                                                        const double k2=4.0*pi*pi*(kx2+ky2+kz2);
//                                                        //use binary numbers to repr permutations
//                                                        //0-> +1
//                                                        //1-> -1
//                                                        const int x=(int)rjl[0]+L[0]/2;
//                                                        const int y=(int)rjl[1]+L[1]/2;
//                                                        const int z=(int)rjl[2]+L[2]/2;
//                                                        const double complex ekr000=cexpValsX[kx][x]*
//                                                                                    cexpValsY[ky][y]*
//                                                                                    cexpValsZ[kz][z];
//                                                        const double complex ekr001=cexpValsX[kx][x]*
//                                                                                    cexpValsY[ky][y]*
//                                                                                    conj(cexpValsZ[kz][z]);
//                                                        const double complex ekr010=cexpValsX[kx][x]*
//                                                                                    conj(cexpValsY[ky][y])*
//                                                                                    cexpValsZ[kz][z];
//                                                        const double complex ekr011=cexpValsX[kx][x]*
//                                                                                    conj(cexpValsY[ky][y])*
//                                                                                    conj(cexpValsZ[kz][z]);
//                                                        const double Ek=exp(-k2*s2/4.0)/k2*creal(ekr000+ekr001+ekr010+ekr011);
//                                                        sum+=2.0*cx*cy*cz*Ek;
//                                                }
//                                        }
//                                }
//                        }
//                        const int lmj=l-j-1;
//                        (*ewald)[j][lmj]=sum;
//                }
//        }
//        fprintf(stderr,"setupEwald_omp: post-omp loop\n");
//        fflush(stderr);
//        //diag element: 'self E'
//        double rjl[d];
//        rjl[0]=0.0;
//        rjl[1]=0.0;
//        rjl[2]=0.0;
//        double sum=0.0;
//        double cx,cy,cz;
//        for(int kx=0;kx<L[0];++kx) {
//                const double kx2=(double)(kx*kx)/(L[0]*L[0]);
//                const double xkx=rjl[0]*kx/L[0];
//                if(kx!=0) cx=1.0;
//                else      cx=0.5;
//                for(int ky=0;ky<L[1];++ky) {
//                        const double ky2=(double)(ky*ky)/(L[1]*L[1]);
//                        const double yky=rjl[1]*ky/L[1];
//                        if(ky!=0) cy=1.0;
//                        else      cy=0.5;
//                        for(int kz=0;kz<L[2];++kz) {
//                                if(kx!=0 || ky!=0 || kz!=0) {
//                                        const double kz2=(double)(kz*kz)/(L[2]*L[2]);
//                                        const double zkz=rjl[2]*kz/L[2];
//                                        if(kz!=0) cz=1.0;
//                                        else      cz=0.5;
//                                        const double k2=4.0*pi*pi*(kx2+ky2+kz2);
//                                        //use binary numbers to repr permutations
//                                        //0-> +1
//                                        //1-> -1
//                                        const double kr000=2.0*pi*(xkx+yky+zkz);
//                                        const double kr001=2.0*pi*(xkx+yky-zkz);
//                                        const double kr010=2.0*pi*(xkx-yky+zkz);
//                                        const double kr011=2.0*pi*(xkx-yky-zkz);
//                                        const double Ek=exp(-k2*s2/4.0)/k2*
//                                                creal(cexp(I*kr000)+cexp(I*kr001)+cexp(I*kr010)+cexp(I*kr011));
//                                        sum+=2.0*cx*cy*cz*Ek;
//                                }
//                        }
//                }
//        }
//        (*ewald)[N-1][0]=sum;
//        fprintf(stderr,"setupEwald_omp: work done, pre-cleanup\n");
//        fflush(stderr);
//        free(cexpValsHolderX); free(cexpValsX);
//        free(cexpValsHolderY); free(cexpValsY);
//        free(cexpValsHolderZ); free(cexpValsZ);
//        fprintf(stderr,"setupEwald_omp: post-cleanup. exiting fn\n");
//        fflush(stderr);
//}

//use data N(N-1)/2 + 1 , cos's, not e^ikx's so only sum over +k's
void setupEwald_omp(double*** ewald, int* L, para p) {
        fprintf(stderr,"setupEwald_omp: fn called\n");
        fflush(stderr);
        const int d=3;
        assert(d==3 /*setupEwald assumes...*/);
        const unsigned long N=L[0]*L[1]*L[2];
        //const double err=10e-10; //hardcode
        const double pi=acos(-1);
        const double s2=p->sigma * p->sigma;
        double* ewaldHolder=malloc(((size_t)N*(size_t)(N-1)/2+(size_t)1)*sizeof(*ewaldHolder)); //ensure no wrapping
        assert(ewaldHolder!=NULL);
        *ewald=malloc(((N-1)+1)*sizeof(**ewald));
        assert(*ewald!=NULL);
        double complex* cexpValsHolderX=malloc(L[0]*L[0]*sizeof(*cexpValsHolderX));
        assert(cexpValsHolderX!=NULL);
        double complex* cexpValsHolderY=malloc(L[1]*L[1]*sizeof(*cexpValsHolderY));
        assert(cexpValsHolderY!=NULL);
        double complex* cexpValsHolderZ=malloc(L[2]*L[2]*sizeof(*cexpValsHolderZ));
        assert(cexpValsHolderZ!=NULL);
        double complex** cexpValsX=malloc(L[0]*sizeof(*cexpValsX)); //add precompute of vectors
        assert(cexpValsX!=NULL);
        double complex** cexpValsY=malloc(L[1]*sizeof(*cexpValsY));
        assert(cexpValsY!=NULL);
        double complex** cexpValsZ=malloc(L[2]*sizeof(*cexpValsZ));
        assert(cexpValsZ!=NULL);
        //assert(ewaldHolder!=NULL && *ewald!=NULL && 
        //       cexpValsHolderX!=NULL && cexpValsX!=NULL &&
        //       cexpValsHolderY!=NULL && cexpValsY!=NULL &&
        //       cexpValsHolderZ!=NULL && cexpValsZ!=NULL /*malloc*/);
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
        i=0;
        for(unsigned long j=0;j<N-1;++j) {
                (*ewald)[j]=ewaldHolder+i;
                i+=N-j-1;
        }
        (*ewald)[N-1]=ewaldHolder+i;
        const unsigned long iter=(unsigned long)N*(unsigned long)(N-1)/2;
        fprintf(stderr,"setupEwald_omp: pre-omp loop. L0,L1,L2,N,iter,N*(N-1)/2: %d,%d,%d,%d,%lu,%lu\n",
                L[0],L[1],L[2],N,iter,(unsigned long)N*(unsigned long)(N-1)/2);
        fflush(stderr);
        #pragma omp parallel for schedule(static)
        for(unsigned long m=0;m<iter;++m) {
                unsigned long j=m/N;
                unsigned long l=m%N;
                if(l>j) {
                        ;
                }
                else {
                        j=N-2-j;
                        l=j+l+1;
                }
                const unsigned long lmj=l-j-1;
                (*ewald)[j][lmj]=0.0;
                double rjl[d];
                diff_lat(rjl,j,l,L);
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
                                                const int x=(int)rjl[0]+L[0]/2;
                                                const int y=(int)rjl[1]+L[1]/2;
                                                const int z=(int)rjl[2]+L[2]/2;
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
                                                (*ewald)[j][lmj]+=2.0*cx*cy*cz*Ek;
                                        }
                                }
                        }
                }
        }
        fprintf(stderr,"setupEwald_omp: post-omp loop\n");
        fflush(stderr);
        //diag element: 'self E'
        double rjl[d];
        rjl[0]=0.0;
        rjl[1]=0.0;
        rjl[2]=0.0;
        double sum=0.0;
        double cx,cy,cz;
        for(int kx=0;kx<L[0];++kx) {
                const double kx2=(double)(kx*kx)/(L[0]*L[0]);
                const double xkx=rjl[0]*kx/L[0];
                if(kx!=0) cx=1.0;
                else      cx=0.5;
                for(int ky=0;ky<L[1];++ky) {
                        const double ky2=(double)(ky*ky)/(L[1]*L[1]);
                        const double yky=rjl[1]*ky/L[1];
                        if(ky!=0) cy=1.0;
                        else      cy=0.5;
                        for(int kz=0;kz<L[2];++kz) {
                                if(kx!=0 || ky!=0 || kz!=0) {
                                        const double kz2=(double)(kz*kz)/(L[2]*L[2]);
                                        const double zkz=rjl[2]*kz/L[2];
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
                                        sum+=2.0*cx*cy*cz*Ek;
                                }
                        }
                }
        }
        (*ewald)[N-1][0]=sum;
        fprintf(stderr,"setupEwald_omp: work done, pre-cleanup\n");
        fflush(stderr);
        free(cexpValsHolderX); free(cexpValsX);
        free(cexpValsHolderY); free(cexpValsY);
        free(cexpValsHolderZ); free(cexpValsZ);
        fprintf(stderr,"setupEwald_omp: post-cleanup. exiting fn\n");
        fflush(stderr);
}

void dumpEwald(double** ewald, int* L, FILE* f) {
        fprintf(stderr,"dumpEwald: fn called\n");
        fflush(stderr);
        const int d=3;
        assert(d==3);
        const int N=L[0]*L[1]*L[2];
        const int nDig=DBL_DECIMAL_DIG; //from float.h
        for(int i=0;i<N-1;++i) {
                for(int j=i+1;j<N;++j) {
                        const int jmi=j-i-1;
                        fprintf(f,"%.*e\n",nDig,ewald[i][jmi]);
                }
        }
        fprintf(stderr,"dumpEwald: post-loop\n");
        fflush(stderr);
        fprintf(f,"%.*e\n",nDig,ewald[N-1][0]);
        fprintf(stderr,"dumpEwald: work done. exiting fn\n");
        fflush(stderr);
}

void loadEwald(double*** ewald, int N, char** lines, int nlines)
{
        double* ewaldHolder=malloc((N*(N-1)/2+1)*sizeof(*ewaldHolder));
        *ewald=malloc(((N-1)+1)*sizeof(**ewald));
        assert(ewaldHolder!=NULL && *ewald!=NULL /*malloc*/);
        int i=0;
        for(int j=0;j<N-1;++j) {
                (*ewald)[j]=ewaldHolder+i;
                for(int k=j+1;k<N;++k) {
                        const int kmj=k-j-1;
                        (*ewald)[j][kmj]=atof(lines[i]);
                        ++i;
                }
        }
        (*ewald)[N-1]=ewaldHolder+i;
        (*ewald)[N-1][0]=atof(lines[i]);
}

int main(int argc, char* argv[]) {
        fprintf(stderr,"args:\n");
        for(int i=0;i<argc;++i) fprintf(stderr,"%s ",argv[i]);
        fprintf(stderr,"\n");
        const int nargs=6;
        if(argc!=nargs) {
                fprintf(stderr,"usage: ./setupEwald L0 L1 L2 sigma f\n");
                fprintf(stderr,"nargs!=number args given. (nargs %d\tn given %d) exiting\n",nargs,argc);
                exit(1);
        }
        const int L0=atoi(*(++argv));
        const int L1=atoi(*(++argv));
        const int L2=atoi(*(++argv));
        const int d=3; //assumes d==3
        int L[3];
        L[0]=L0;
        L[1]=L1;
        L[2]=L2;
        fprintf(stderr,"L: %d,%d,%d\n",L[0],L[1],L[2]);
        struct parameters pHolder;
        para p=&pHolder;
        p->sigma=atof(*(++argv));

        fprintf(stderr,"main: calling setupEwald_omp\n");
        fflush(stderr);
        double** ewald=NULL;
        //setupEwald(&ewald,L,p);
        setupEwald_omp(&ewald,L,p);
        fprintf(stderr,"main: setupEwald_omp returned\n");
        fflush(stderr);

        fprintf(stderr,"main: calling dumpEwald\n");
        fflush(stderr);
        FILE* f=fopen(*(++argv),"w");
        assert(f!=NULL /*fopen*/);
        dumpEwald(ewald,L,f);
        fprintf(stderr,"main: dumpEwald returned\n");
        fflush(stderr);

        fprintf(stderr,"main: exiting successfully\n");
        fflush(stderr);
        return 0;
}
