#!/bin/bash

# Script to run package mcdiff
#
# how to run:
#    chmod u+x example.sh        # to make the file executable
#    ./example.sh

export PYTHONPATH=$HOME/lib/python/:$PYTHONPATH

# manage dirs/files
#------------------
# working with replicate A
outdir=output/HexdWat/A/out.ncosf10d6.run1/
if [ ! -d "output" ];  then mkdir output ;  fi
if [ ! -d "output/HexdWat" ];  then mkdir output/HexdWat ;  fi
if [ ! -d "output/HexdWat/A" ];  then mkdir output/HexdWat/A ;  fi
if [ ! -d "$outdir" ]; then mkdir $outdir ; fi 

echo outdir  $outdir

# define the location of the transition matrices
#------------------------------------------------
line="data/HexdWat/transitions.nbins100.20.pbc.A.dat"
echo LINE $line

# run the Monte Carlo (Bayesian Analysis)
#----------------------------------------
nmc=5000    # test run

run-mcdiff --nmc $nmc -T 1 --Tend 1 \
    $line \
   --model CosinusModel --ncosF 10 --ncosD 6 \
   --nmc_update 1000 \
   --dv=0.05 --dw=0.05 \
   -o $outdir/out.nbins100.20.pbc.dat \
    > $outdir/out.nbins100.20.pbc.log

# in this example:
# using 10 for the cosine basis set of F
# using 6 for the cosine basis set of D (name 'ncosf10d6')
# "Bayesian analysis" when parameter temperature is set to T=1, Tend=1
# nmc is number of Monte Carlo steps
# both F profile and D profile are updated (dv>0, dw>0)
# steps in F and D parameters space are updated every 1000 Monte Carlo steps

# plot the results
#-------------------
plotresults $outdir/out.nbins100.20.pbc.dat -o $outdir/finalprofile
#plotresults $outdir/out.nbins100.20.pbc.dat -o $outdir/averageprofile
