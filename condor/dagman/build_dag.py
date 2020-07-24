#!/usr/bin/env python

#Let's make a dag.
# ./build_dag.py > newdag.dag
# then hand it to condor
import numpy as np

jobbase = 'stacked_sens_test'
script = '/data/condor_builds/users/blaufuss/umd_computing_examples/jupyter/stacked_sensitivity_refactor.py'

counter = 0

time_windows = np.logspace(1,6,6)
spectra = np.linspace(-3.5,-1,11)

for twin in time_windows:
    for spec_ind in spectra:
        #twin = 100
        #spec_ind = -2.0
        command = f'python {script} {twin} {spec_ind} {counter}'
        job_id = f'{jobbase}_{counter}'
        print(f'JOB {job_id} /data/condor_builds/users/blaufuss/umd_computing_examples/condor/dagman/submit.sub')
        print(f'VARS {job_id} JOBNAME="{job_id}" command="{command}"')
        counter += 1
