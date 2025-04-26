#!/bin/bash

# Run the small world pc calculations for alpha = 0..6

gtype=7
#nsites=1048576
nsites=16384
mcull=3
#nsamples=100
nsamples=1
coord=4

#./HIABP $gtype $nsites $mcull $nsamples 0 $coord
#./ABPCDist_te $gtype $nsites $mcull $nsamples 0 $coord
./MaxClusterProbability_p $gtype $nsites $mcull $nsamples 0 $coord
