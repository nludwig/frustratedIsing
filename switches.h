#pragma once

#define ISI_ON 1                //0 -> Ising NN off; 1 -> Ising NN on
#define COUL_ON 1               //0 -> pure Ising; 1 -> Frustrated Ising
#define FFT_ON 1                //0 -> only real space interactions; 1 -> LR part of Coulomb via FFT (Ewald sum)
#define SHIFT_COUL 1            //0 -> E_c truncated @ R_c=coulCutoff*sigma; 1 -> E_c truncated & shifted to 0 @..
#define EQUI_INNER_DUMP 1       //0 -> output test every N attempted moves; 1-> output test every move
#define ISI_USE_SIGN 0          //0 -> use q+, q- for Ising energies; 1-> use sign val (q+ -> +1, q- -> -1) for Ising energies
#define ISI_SPIN_UP 1.0         //to generalize, use these macros to choose the ising spin val that is used
#define ISI_SPIN_DN -1.0        //for spin up/down instead of just hardcoding in +1,-1; no effect w/ ISI_USE_SIGN=0
#define TEST_BUILD 0            //enables extra checks & outputs, such as full vs local E comparisons
#define ENERGY_LENGTH 10        //isvsv, isvsu, isusu, cSRsvsv, cSRsvsu, cSRsusu, cLR, ihh, isvh, isuh
