import sources, time
import os, sys, glob, numpy as np, matplotlib, scipy, healpy as hp
from matplotlib import pyplot as plt, colors
from scipy import stats, interpolate, optimize
import numpy.lib.recfunctions as rf

# #### Plotting - Let's read the output from disk and make pretty pictures.

OutDir = '/data/condor_builds/users/blaufuss/stacked_outputs/'
import glob
filenames = glob.glob(OutDir + '*.npy')

#time_windows = np.logspace(1,6,6)
#spectra = np.linspace(-3.5,-1,6)

loaded = []
for file in sorted(filenames):
    loaded.append(np.load(file, allow_pickle = True))
s = []
for entry in loaded:
    # each entry has [index, timewindow, flux norm]
    print(entry[0], entry[1], entry[2]*1e10*entry[1])
    s.append([entry[0], entry[1], entry[2]*1e10*entry[1]])
sarray = np.array(s)
print(sarray)
# Make a list of unique spectral indices
idxs = np.unique(sarray[:,0])
# split out the TW and time int. flux limits for each
foo = np.array( [list(sarray[sarray[:,0]==i,1:]) for i in idxs] )
# zip with spectral indices fo easy plotting
IndexSplitTWLimits = np.array(list(zip(idxs, foo)))
print(IndexSplitTWLimits)
fig, ax = plt.subplots(figsize = (12,12))
import itertools
marker = itertools.cycle(('x', 's', 'd', 'o', '*', 'p', '3'))
# plot a line per spectral index
for entry in IndexSplitTWLimits:
    idx = entry[0]
    tw_limits = np.array(entry[1])
    np.ndarray.sort(tw_limits,axis=0)
    print(idx, tw_limits[:,0],tw_limits[:,1])
    labeled = 'Index: ' + str(idx)
    plt.plot(np.log10(tw_limits[:,0]),tw_limits[:,1],
	label=labeled, marker=next(marker), markersize = 9)
ax.set_yscale('log')
ax.set_title('Sensitivity for a stacked analysis of unrelated sources', fontsize = 18)
ax.set_xlabel('Log10(time window [s])', fontsize = 16)
ax.set_ylabel(f'fluence ($E_0^2F$)', fontsize = 16)
leg = ax.legend();
plt.savefig(OutDir + 'results.pdf')
