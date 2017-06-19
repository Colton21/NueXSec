# Author: Colton Hill
# This file should briefly examine the optical information

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import root_numpy as rnp
from tqdm import tqdm


def Pandafy(fileName, tree):
    df = pd.DataFrame(rnp.root2array(fileName, tree))
    return df

infile = 'nue_xsec_extraction.root'
infile_tree_name = 'NueXsec/optical_tree'

df = Pandafy(infile, infile_tree_name)

nEvents = len(df.index)
# check if data frame is empty
if(nEvents == 0):
    print >> sys.stderr, 'Data Frame is Null!'
    exit(1)

print 'Number of Events: %s' % (nEvents)

# Filter for events with fewer than 50 PE and out-of-time
pe_threshold = 150
start_time_window = 5.0
end_time_window = 16.0
df = df.drop(df[(df.OpFlashPE <= pe_threshold)].index)
nPassed_PE = len(df.index)
print 'Number of Events Post PE Cut of %s : %s' % (pe_threshold, nPassed_PE)
df = df.drop(df[(df.OpFlashTime >= end_time_window) |
                (df.OpFlashTime <= start_time_window)].index)
nPassed_Time = len(df.index)
print 'Number of Events Post Time Cut from %s to %s : %s' % (start_time_window, end_time_window, nPassed_Time)
passing_fraction = (float(nPassed_Time) / float(nEvents)) * 100.
print 'Total Passing Fraction: ', passing_fraction

photoelectrons = df.OpFlashPE
timing = df.OpFlashTime
plt.figure()
photoelectrons.plot.hist(alpha=0.5)
plt.show()
timing.plot.hist(alpha=0.5)
plt.show()
