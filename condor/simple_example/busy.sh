#!/bin/bash
# file name: busy.sh

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/setup.sh` 
/data/condor_builds/users/blaufuss/umd_computing_examples/condor/simple_example/busy.py $1
