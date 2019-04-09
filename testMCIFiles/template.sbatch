#!/bin/bash
#SBATCH --job-name=NAME
#SBATCH --output=./NAME.out
#SBATCH --error=./NAME.err
#SBATCH --time=36:00:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END
#SBATCH --mail-user=nickludwig@uchicago.edu

#./MCIsing verbose?(0/1) parametersIn dataIn(or'none') enFreq dumpFreq equiSteps prodSteps rngSeed,rngSeq confDump
