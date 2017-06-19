# Author: Colton Hill
# This file should briefly examine the reco-true information

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
infile_tree_name = 'TrueRecoMon/pandora'

df = Pandafy(infile, infile_tree_name)

nEvents = len(df.index)
# check if data frame is empty
if(nEvents == 0):
    print >> sys.stderr, 'Data Frame is Null!'
    exit(1)
# for guidance of data product names
print df.head(1)

print 'Number of Events: %s' % (nEvents)

# let's trim the sample down to only events in the TPC
tpc_x1 = 0
tpc_x2 = 256.35
tpc_y1 = -116.5
tpc_y2 = 116.5
tpc_z1 = 0
tpc_z2 = 1036.8
df = df.drop(df[(df.mcVtxX < tpc_x1) | (df.mcVtxX > tpc_x2)].index)
df = df.drop(df[(df.mcVtxY < tpc_y1) | (df.mcVtxY > tpc_y2)].index)
df = df.drop(df[(df.mcVtxZ < tpc_z1) | (df.mcVtxZ > tpc_z2)].index)
nEvents_inTPC = len(df.index)
print 'Events in TPC: ', nEvents_inTPC

# let's try and remove all of the "bad" events first, before looking at reco
df = df.drop(df[df.nMCParticles == 0].index)
nEvents_MCParticle = len(df.index)
print 'Number of Events w/ MCParticles: ', nEvents_MCParticle

dummy_df = df
dummy_df = dummy_df.drop(dummy_df[dummy_df.mcEnergy != 0].index)
energy = dummy_df.mcEnergy
print 'Particles with 0.0 MCEnergy: ', len(energy.index)
# plt.figure()
# energy.plot.hist(alpha=0.5)
# plt.show()

print energy
