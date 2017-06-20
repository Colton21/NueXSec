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
import timeit

# Start function definitions

# define a function to get a dataframe using pandas


def Pandafy(fileName, tree):
    df = pd.DataFrame(rnp.root2array(fileName, tree))
    return df

# let's trim the sample down to only events in the TPC


def inTPC(dataframe):
    tpc_x1 = 0
    tpc_x2 = 256.35
    tpc_y1 = -116.5
    tpc_y2 = 116.5
    tpc_z1 = 0
    tpc_z2 = 1036.8
    dataframe = dataframe.drop(
        dataframe[(dataframe.mcVtxX < tpc_x1) | (dataframe.mcVtxX > tpc_x2)].index)
    dataframe = dataframe.drop(
        dataframe[(dataframe.mcVtxY < tpc_y1) | (dataframe.mcVtxY > tpc_y2)].index)
    dataframe = dataframe.drop(
        dataframe[(dataframe.mcVtxZ < tpc_z1) | (dataframe.mcVtxZ > tpc_z2)].index)
    return dataframe

# some events have no MC particles - these are uninteresting!


def removeZeroMCP(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.nMCParticles == 0].index)
    return dataframe

# some mc particles have zero MC energy - we should remove these


def removeZeroMCEng(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.mcEnergy == 0.0].index)
    return dataframe


# Register `pandas.progress_apply` and `pandas.Series.map_apply` with `tqdm`
# (can use `tqdm_gui`, `tqdm_notebook`, optional kwargs, etc.)
tqdm.pandas(desc="Drop Time")

# Begin timer to measure execution duration
start_time = timeit.default_timer()
print 'Begin Reco True Analysis: ', start_time

# open file and tree root table
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

# let's try and remove all of the "bad" events first, before looking at reco
df = removeZeroMCP(df)
nEvents_MCParticle = len(df.index)
print 'Number of Events w/ MCParticles: ', nEvents_MCParticle

df = removeZeroMCEng(df)
nEvents_MCEnergy = len(df.index)
print 'Number of Events w/ Non-Zero MCEng: ', nEvents_MCEnergy

good_fraction = (float(nEvents_MCEnergy) / float(nEvents)) * 100.
print 'Good Fraction: ', good_fraction


# filter events to be inTPC
df = inTPC(df)
nEvents_inTPC = len(df.index)
print 'Events in TPC: ', nEvents_inTPC


elapsed = timeit.default_timer() - start_time
print 'Time Elapsed: ', elapsed
