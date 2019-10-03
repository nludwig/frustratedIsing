#include "forRandOrder.h"

//generate uniformly psuedo-random ordered "for" loop
//(ie. non-repeating sequence in [0,l))

int* forRandOrder(int l, pcg32_random_t* rng) {
        int* order=malloc(l*sizeof(*order));
        bool* taken=malloc(l*sizeof(*taken));
        assert(order!=NULL && taken!=NULL /*malloc*/);
        for(int i=0;i<l;++i) taken=false;
        int i=0;
        while(!allTrue(taken,l)) {
                const int rn=pcg32_boundedrand_r(rng,l);
                if(taken[rn]==false) {
                        taken[rn]=true;
                        order[i]=rn;
                        ++i;
                }
        }
        free(taken);
        return order;
}

bool allTrue(bool* a, int la) {
        for(int i=0;i<la;++i) {
                if(a[i]==false) {
                        return false;
                }
        }
        return true;
}
