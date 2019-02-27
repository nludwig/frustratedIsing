#pragma once

#define ISI_ON 1        //0 -> Ising NN off;    1 -> Ising NN on
#define COUL_ON 1       //0 -> pure Ising;      1 -> Frustrated Ising
#define FFT_ON 1        //0 -> only real space interactions;    1 -> LR part of Coulomb via FFT (Ewald sum)
#define SHIFT_COUL 1    //0 -> E_c truncated @ R_c=coulCutoff*sigma;   1 -> E_c truncated & shifted to 0 @..
#define EXPAND_E 1      //1 -> spoof outputs separated E_i, E_c (and E_c.fft) for each step
#define EQUI_INNER_DUMP 1
#define ISI_USE_SIGN 0  //0 -> use q+, q- for Ising energies. 1-> use sign val (q+ -> +1, q- -> -1) for Ising energies
#define ISI_SPIN_UP 1.0  //to generalize, use these macros to choose the ising spin val that is used
#define ISI_SPIN_DN -1.0 //for spin up/down instead of just hardcoding in +1,-1; no effect w/ ISI_USE_SIGN=0
#define TEST_BUILD 0
#define ENERGY_LENGTH 10
