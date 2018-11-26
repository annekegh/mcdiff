#!/bin/bash

# Script to run package mcdiff
# $1 is the name of the system (to find directories)
# $2 is the run
# $3 is the lag time
#
# how to run:
#    chmod u+x example.sh        # to make the file executable
#    ./example.sh HexdWat 1 20   # system=HexdWat, run=1, lagtime=20
#    ./example.sh HexdWat 2 20   # continue the previous Monte Carlo run

export PYTHONPATH=$HOME/lib/python/:$PYTHONPATH

#nmc=2000    # test run
#nmc=50000   # equilibration of parameters, more could be needed
nmc=100000  # production run

system=$1       # the system, here HexdWat
run=$2          # run  (if you do multiple after each other)
lt=$3           # lag time

nbins=100       # number of bins
fk=10    # basis set for F
dk=6     # basis set for D

# manage dirs/files
#------------------
outdir=output/$system/out.ncosf${fk}d${dk}.run$run/
if [ ! -d "output" ];  then mkdir output ;  fi
if [ ! -d "output/$system" ];  then mkdir output/$system ;  fi
if [ ! -d "$outdir" ]; then mkdir $outdir ; fi 

if [ "$run" -gt 1 ]; then
  # when continuing previous run, point to the previous directory
  initrun=`echo "$run-1" | bc`   # compute the previous run
  initdir=output/$system/out.ncosf${fk}d${dk}.run$initrun/
fi

# print some info
#----------------
echo run     $run
echo outdir  $outdir
echo initrun $initrun
echo initdir $initdir


# define the location of the transition matrices
#------------------------------------------------
## use all replicas
#line=""
#for letter in A B C D
#do
#  #extend the line
#  dir="data/$system/"
#  line="$line $dir/transitions.nbins$nbins.$lt.pbc.$letter.dat"
#done

# ... or, just use one file with letter  # this is PREFERRED!
# after having executed merge.py in directory manipulate-transition-matrices/
line="manipulate-transition-matrices/transitions.merge.dat"
echo LINE $line
# this line will be fed to the Monte Carlo

# advice: Merge the transition matrices of the different replicates first,
# to a new merged transition matrix, and then run the Monte Carlo.
# This is faster + more accurate than running Monte Carlo with
# a set of transition matrices.

# run the Monte Carlo (Bayesian Analysis)
#----------------------------------------
if [ "$run" -gt 1 ]; then
  run-mcdiff --nmc $nmc -T 1 --Tend 1 \
    $line \
   --model CosinusModel --ncosD $dk --ncosF $fk \
   --nmc_update 1000 \
   --dv=0.05 --dw=0.05 \
   --initf $initdir/out.nbins$nbins.$lt.pbc.dat \
   -o $outdir/out.nbins$nbins.$lt.pbc.dat \
    > $outdir/out.nbins$nbins.$lt.pbc.log

else
  run-mcdiff --nmc $nmc -T 1 --Tend 1 \
    $line \
   --model CosinusModel --ncosD $dk --ncosF $fk \
   --nmc_update 1000 \
   --dv=0.05 --dw=0.05 \
   -o $outdir/out.nbins$nbins.$lt.pbc.dat \
    > $outdir/out.nbins$nbins.$lt.pbc.log
fi

# in this example:
# using 10 for the cosine basis set of F
# using 6 for the cosine basis set of D (name 'ncosf10d6')
# "Bayesian analysis" when parameter temperature is set to T=1, Tend=1
# nmc is number of Monte Carlo steps
# both F profile and D profile are updated (dv>0, dw>0)
# steps in F and D parameters space are updated every 1000 Monte Carlo steps

#   --initf $initdir/out.nbins$nbins.$lt.pbc.dat \
# note about this line:
# to be used when you start from the previous profile, so in run2, run3, ...

