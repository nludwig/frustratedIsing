#include "prebuildEwald.h"

int N=0;
int d=0;
int* L=NULL;
int** posArr=NULL;

double* setupEwald(int* L, double sigma) {
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

void dumpEwald(double* ewald, int* L, FILE* f) {
        fprintf(stderr,"dumpEwald: fn called\n");
        fflush(stderr);
        const int d=3;
        assert(d==3);
        const int N=L[0]*L[1]*L[2];
        const int nDig=DBL_DECIMAL_DIG; //from float.h
        for(int i=0;i<N;++i) {
                fprintf(f,"%.*e\n",nDig,ewald[i]);
        }
        fprintf(stderr,"dumpEwald: post-loop\n");
        fprintf(stderr,"dumpEwald: work done. exiting fn\n");
        fflush(stderr);
}

int** setupPos() {
        extern int N,d;
        assert(d==3);
        extern int* L;
        int* posArrHolder=malloc(N*d*sizeof(*posArrHolder));
        int** posArray=malloc(N*sizeof(*posArray));
        assert(posArrHolder!=NULL && posArray!=NULL);
        for(int i=0;i<N;++i) {
                posArray[i]=posArrHolder+i*d;
        }
        int i=0;
        for(int x=0;x<L[0];++x) {
                for(int y=0;y<L[1];++y) {
                        for(int z=0;z<L[2];++z) {
                                posArray[i][0]=x;
                                posArray[i][1]=y;
                                posArray[i][2]=z;
                                ++i;
                        }
                }
        }
        return posArray;
}

void getPos(int n, int* pos) {
        extern int d;
        extern int** posArr;
        for(int e=0;e<d;++e) {
                pos[e]=posArr[n][e];
        }
}

void minImageVect(int* pos, int* L_in) {
        if(L_in==NULL) {
                extern int* L;
                L_in=L;
        }
        extern int d;
        extern int* bc;
        for(int e=0;e<d;++e) {
                //assert(bc[2*e]==1 || bc[2*e+1]==1 /*wrapIntoL: bc.s must be set*/);
                while(pos[e]<-L_in[e]/2) pos[e]+=L_in[e];
                while(pos[e]>=L_in[e]/2) pos[e]-=L_in[e];
        }
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
        d=3; //assumes d==3
        int Lh[3];
        Lh[0]=L0;
        Lh[1]=L1;
        Lh[2]=L2;
        L=Lh;
        fprintf(stderr,"L: %d,%d,%d\n",L[0],L[1],L[2]);
        N=L[0]*L[1]*L[2];
        const double sigma=atof(*(++argv));
        posArr=setupPos();

        fprintf(stderr,"main: calling setupEwald\n");
        fflush(stderr);
        double* ewald=setupEwald(L,sigma);
        fprintf(stderr,"main: setupEwald returned\n");
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
