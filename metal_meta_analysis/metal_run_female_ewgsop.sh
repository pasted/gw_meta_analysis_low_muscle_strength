#!/bin/bash
#PBS -V   # export all environment variables to the batch job
#PBS -d . # set working directory
#PBS -A Research_Project-MRC158833
#PBS -q sq # change to mrchq for high-memory queue / mrcq for mrc queue
#PBS -l walltime=90:00:00,ncpus=1,mem=10G


/gpfs/mrc0/projects/Research_Project-MRC158833/programs/metal/metal < _syntax.2020-03-13.gc.female.ewgsop.METAL.se.txt