#!/usr/bin/env python

# A do-nothing example
import numpy as np
import sys

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} length')
    exit(-1)

length = int(sys.argv[1])

sizer = 300000


for idx in range(0,length):
    tester1 = []
    tester2 = []
    for idx2 in range(0,sizer):
        tester1.append(float(idx2))
        tester2.append(tester1[idx2]/(tester1[idx2]+idx+0.01))
print('All done')
