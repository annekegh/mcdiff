#!/bin/bash

# Script to run package mcdiff
#
# how to run:
#    chmod u+x example.sh        # to make the file executable
#    ./example.sh

export PYTHONPATH=$HOME/lib/python/:$PYTHONPATH

# manage dirs/files
#------------------
initdir=output/HexdWat/A/out.ncosf10d6.run1/
if [ ! -d "output" ];  then
  echo initdir is $initdir
  echo initdir does not exist
  exit
fi
outdir=output/HexdWat/A/out.ncosf10d6.run2/
if [ ! -d "$outdir" ]; then mkdir $outdir ; fi 

echo outdir  $outdir

# define the location of the transition matrices
#------------------------------------------------
line="data/HexdWat/transitions.nbins100.20.pbc.A.dat"
echo LINE $line

# run the Monte Carlo (Bayesian Analysis)
#----------------------------------------
nmc=100000    # test run

run-mcdiff --nmc $nmc -T 1 --Tend 1 \
    $line \
   --model CosinusModel --ncosF 10 --ncosD 6 \
   --nmc_update 1000 \
   --dv=0.05 --dw=0.05 \
   --initf $initdir/out.nbins100.20.pbc.dat \
   -o $outdir/out.nbins100.20.pbc.dat \
    > $outdir/out.nbins100.20.pbc.log

# in this example:
# using 10 for the cosine basis set of F
# using 6 for the cosine basis set of D (name 'ncosf10d6')
# "Bayesian analysis" when parameter temperature is set to T=1, Tend=1
# nmc is number of Monte Carlo steps
# both F profile and D profile are updated (dv>0, dw>0)
# steps in F and D parameters space are updated every 1000 Monte Carlo steps

#   --initf $initdir/out.nbins100.20.pbc.dat \
# note about this line:
# this is used as the initialization profile from the previous run
# to be used when you restart from the previous profile

# plot the results
#-------------------
plotresults $outdir/out.nbins100.20.pbc.dat -o $outdir/finalprofile
plotresults $outdir/out.nbins100.20.pbc.dat.pic -o $outdir/averageprofile-errorbars --pic
plotresults $outdir/out.nbins100.20.pbc.dat.pic.ave.dat -o $outdir/averageprofile 
