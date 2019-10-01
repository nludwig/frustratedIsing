#include "shapes.h"

int** buildCube(int l) {
        extern int d;
        int nsites=1;
        for(int e=0;e<d;++e) nsites*=l;
        int** sites=malloc(nsites*sizeof(*sites));
        int* points=malloc(nsites*d*sizeof(*points));
        assert(sites!=NULL && points!=NULL /*malloc*/);
        for(int i=0;i<nsites;++i) sites[i]=points+d*i;
        int i=0;
        for(int x=0;x<l;++x) { //row major
                for(int y=0;y<l;++y) {
                        for(int z=0;z<l;++z) {
                                sites[i][0]=x;
                                sites[i][1]=y;
                                sites[i][2]=z;
                                ++i;
                        }
                }
        }
        return sites;
}

int** buildRectangle(int* l) {
        extern int d;
        int nsites=1;
        for(int e=0;e<d;++e) nsites*=l[e];
        int** sites=malloc(nsites*sizeof(*sites));
        int* points=malloc(nsites*d*sizeof(*points));
        assert(sites!=NULL && points!=NULL /*malloc*/);
        for(int i=0;i<nsites;++i) sites[i]=points+d*i;
        int i=0;
        for(int x=0;x<l[0];++x) { //row major
                for(int y=0;y<l[1];++y) {
                        for(int z=0;z<l[2];++z) {
                                sites[i][0]=x;
                                sites[i][1]=y;
                                sites[i][2]=z;
                                ++i;
                        }
                }
        }
        return sites;
}

int* getShapeInterface(int** shape, int nsites, int* ninterface, int** Lcube, int* border) {
        extern int d;
        *Lcube=embedShapeInCube(shape,nsites,border);
        const int nNN=nsites*d*2;
        int* nnList=malloc(nNN*sizeof(*nnList));
        int* cmpList=malloc(nsites*sizeof(*cmpList));
        assert(nnList!=NULL && cmpList!=NULL /*malloc*/);
        int j=0;
        for(int i=0;i<nsites;++i) {
                cmpList[i]=getWrappedSite(shape[i],*Lcube);
                for(int e=0;e<d;++e) {
                        shape[i][e]-=1; //'-1' neigh.s
                        nnList[j++]=getWrappedSite(shape[i],*Lcube);
                        shape[i][e]+=1;

                        shape[i][e]+=1; //'+1' neigh.s
                        nnList[j++]=getWrappedSite(shape[i],*Lcube);
                        shape[i][e]-=1;
                }
        }
        qsort(nnList,nNN,sizeof(int),(int (*)(const void*,const void*))iCmp);
        const int nCubeShapePlusInterface=uniqueint(&nnList,nNN);
        *ninterface=diff(&nnList,nCubeShapePlusInterface,&cmpList,nsites);
        return nnList;
}

int getSA(int** shape, int nsites, int** Lcube, int* border) {
        int SA;
        int* interface=getShapeInterface(shape,nsites,&SA,Lcube,border);
        free(interface);
        return SA;
}

int** getShapeInterfacePos(int** shape, int nsites, int* ninterface, int** Lcube, int* border) {
        int* interface=getShapeInterface(shape,nsites,ninterface,Lcube,border);
        extern int d;
        int** interShape=malloc(*ninterface*sizeof(*interShape));
        int* interShapePos=malloc(*ninterface*d*sizeof(*interShapePos));
        assert(interShape!=NULL && interShapePos!=NULL /*malloc*/);
        for(int i=0;i<*ninterface;++i) {
                interShape[i]=interShapePos+i*d;
                getPos(interface[i],interShape[i]);
        }
        if(*ninterface>0) {
                free(interface);
        }
        return interShape;
}


//assumes *arr is sorted
int uniqueint(int** arr, int larr) {
        int i=0,j=1;
        while(i<larr) {
                while(j<larr && iCmp(*arr+i,*arr+j)==0) ++j;
                if(j>=larr) break;
                (*arr)[++i]=(*arr)[j++];
        }
        return i+1;
}

int uniqueuint32_t(uint32_t** arr, int larr) {
        int i=0,j=1;
        while(i<larr) {
                while(j<larr && uint32_tCmp(*arr+i,*arr+j)==0) ++j;
                if(j>=larr) break;
                (*arr)[++i]=(*arr)[j++];
        }
        return i+1;
}

int uniquelattice(lattice** arr, int larr) {
        int i=0,j=1;
        while(i<larr) {
                while(j<larr && lattCmp(*arr+i,*arr+j)==0) ++j;
                if(j>=larr) break;
                (*arr)[++i]=(*arr)[j++];
        }
        if((*arr)[0]==NULL) {
                for(int k=1;k<i+1;++k) (*arr)[k-1]=(*arr)[k];
                i-=1;
        }
        return i+1;
}

int* enclosingCubeSize(int** shape, int nsites)
        //find maxes
        extern int d;
        int min[d];
        int max[d];
        for(int e=0;e<d;++e) {
                min[e]=INT_MAX;
                max[e]=0;
        }
        for(int i=0;i<nsites;++i) {
                for(int e=0;e<d;++e) {
                        if(shape[i][e]<min[e]) min[e]=shape[i][e];
                        if(shape[i][e]>max[e]) max[e]=shape[i][e];
                }
        }
        int* shapeMinMax=malloc(d*2*sizeof(*shapeMinMax));
        assert(shapeMinMax!=NULL /*malloc*/);
        for(int e=0;e<d;++e) {
                shapeMinMax[e]=min[e];
                shapeMinMax[e+d]=max[e];
        }
        return shapeMinMax;
}

int* embedShapeInCube(int** shape, int nsites, int* border)
        extern int d;
        int* minMax=enclosingCubeSize(shape,nsites);
        int* LwBorder=malloc(d*sizeof(*LwBorder));
        assert(LwBorder!=NULL /*malloc*/);
        for(int e=0;e<d;++e) LwBorder[e]=minMax[e+d]-minMax[e]+2*border[e]+1;
        for(int i=0;i<nsites;++i) {
                for(int e=0;e<d;++e) shape[i][e]+=border[e]-minMax[e];
        }
        free(minMax);
        return LwBorder;
}

int** buildNeighCube(double Rc) {
        const int iRc=(int)ceil(Rc);
        const int l=2*iRc+1;
        int** cube=buildCube(l);
        extern int d;
        int Ncube=1;
        for(int e=0;e<d;++e) Ncube*=l;
        for(int i=0;i<Ncube;++i) for(int e=0;e<d;++e) cube[i][e]-=iRc;
        return cube;
}

//assumptions:
//1) a,b sorted
//2) a,b unique (no repeats w/in a or w/in b)
//3) la>=lb; more specfically, la contains lb
int diff(int** a, int la, int** b, int lb) {
        assert(la>=lb /*diff assump violated*/);
        int i=0,j=0;
        while(i<la && lb>0) {
                while(j<lb) {
                        if((*a)[i]==(*b)[j]) {
                                --la;
                                --lb;
                                for(int k=i;k<la;++k) (*a)[k]=(*a)[k+1];
                                for(int k=j;k<lb;++k) (*b)[k]=(*b)[k+1];
                        }
                        else ++j;
                }
                ++i;
                j=0;
        }
        free(*b); *b=NULL;
        if(la==0) {
                free(*a); *a=NULL;
                return 0;
        }
        else {
                *a=realloc(*a,la*sizeof(*a));
                assert(*a!=NULL /*realloc*/);
                return la;
        }
}

double* matrixOnVector(double** M, double* v) {
        extern int d;
        double* u=malloc(d*sizeof(*u));
        assert(u!=NULL /*malloc*/);
        for(int i=0;i<d;++i) u[i]=0.0;
        for(int i=0;i<d;++i) {
                for(int j=0;j<d;++j) {
                       u[i]+=M[i][j]*v[j];
                }
        }
        return u;
}

void matrixOnVectorInPlace(double** M, double* v) {
        double* u=matrixOnVector(M,v);
        extern int d;
        for(int i=0;i<d;++i) v[i]=u[i];
        free(u);
}

double** matrixOnMatrix(double** M1, double** M2) {
        extern int d;
        double* Ahold=malloc(d*d*sizeof(*Ahold));
        double** A=malloc(d*sizeof(*A));
        assert(Ahold!=NULL && A!=NULL /*malloc*/);
        for(int i=0;i<d;++i) {
                A[i]=Ahold+i*d;
                for(int j=0;j<d;++j) {
                        A[i][j]=0.0;
                }
        }
        for(int i=0;i<d;++i) {
                for(int j=0;j<d;++j) {
                        for(int k=0;k<d;++k) {
                                A[i][k]+=M1[i][j]*M2[j][k];
                        }
                }
        }
        return A;
}

void matrixOnMatrixInPlace(double** M, double** Mout) {
        double** A=matrixOnMatrix(M,Mout);
        extern int d;
        for(int i=0;i<d;++i) {
                for(int j=0;j<d;++j) {
                        Mout[i][j]=A[i][j];
                }
        }
        free(*A); free(A);
}

void copyMatrixInPlace(double** M, double** Mout, int d) {
        for(int i=0;i<d;++i) {
                for(int j=0;j<d;++j) {
                        Mout[i][j]=M[i][j];
                }
        }
}

double** copyMatrix(double** M, int d) {
        double* Ah=malloc(d*d*sizeof(*Ah));
        double** A=malloc(d*sizeof(*A));
        assert(Ah!=NULL && A!=NULL);
        for(int i=0;i<d;++i) A[i]=Ah+i*d;
        copyMatrixInPlace(M,A,d);
        return A;
}

double*** buildPiO2RotationMatricies(int d) {
        double* Rhh=malloc(d*d*d*sizeof(*Rhh));
        double** Rh=malloc(d*d*sizeof(*Rh));
        double*** R=malloc(d*sizeof(*R));
        assert(Rhh!=NULL && Rh!=NULL && R!=NULL);
        for(int i=0;i<d;++i) {
                R[i]=Rh+i*d;
                for(int j=0;j<d;++j) {
                        Rh[j+i*d]=Rhh+(j+i*d)*d;
                        for(int k=0;k<d;++k) {
                                if(j==i || k==i || j==k) R[i][j][k]=0.0;
                                else {
                                        int modj=j-i, modk=k-i;
                                        if(modj<0) modj+=d;
                                        if(modk<0) modk+=d;
                                        modj%=d;
                                        modk%=d;
                                        R[i][j][k]=(double)(modk-modj);
                                }
                        }
                }
                R[i][i][i]=1.0;
        }
        return R;
}
