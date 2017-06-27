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
from reco_true_functions import *

# Start function definitions

# define a function to get a dataframe using pandas


def Pandafy(fileName, tree):
    df = pd.DataFrame(rnp.root2array(fileName, tree))
    return df


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
print 'Number of MC Particles: %s' % (nEvents)

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
print 'MC Particles in TPC: ', nEvents_inTPC


# let's look at some characteristics of the dataset
df_zero = zeroIndex(df)
df_nues = df_zero.drop(df_zero[df_zero.mcPdg != 12].index)
nMCNues = len(df_nues.index)
print 'Num MC Nues: ', nMCNues

# number of MCParticles
n_mcpart = df_zero.nMCParticles
fig_mcpart = plt.figure()
ax = fig_mcpart.add_subplot(111)
n_mcpart.plot.hist(bins=18, alpha=0.5)
ax.set_xlabel('MC Particle per event')

# number of PFParticles
df_valid_pfpart = validPFParticle(df)

# let's get just the PFP Neutrino candidates
df_pfp_neutrinos = getPfpNeutrino(df_valid_pfpart)
# what is the number of pfp neutrinos per event?
df_pfp_neutrino_pdgs = df_zero.nNeutrinoPfos
fig_pfpnu = plt.figure()
ax = fig_pfpnu.add_subplot(111)
df_pfp_neutrino_pdgs.plot.hist(alpha=0.5)
ax.set_xlabel('PFP Neutrino per event')

# how many 'primary' pfps
# does this include the neutrino in the count for the pfo?
# it should not include the neutrino, however this counter takes into account
# pfps which may have been reconstructed without a valid pdg code?
# fig_pfps = plt.figure()
# ax = fig_pfps.add_subplot(111)
# df_pfp_neutrinos.nPrimaryPfos.plot.hist(alpha=0.5)
# ax.set_xlabel('Primary PFPs per event')
# print df_pfp_neutrinos.nPrimaryPfos
df_mcprimary = getMCPrimary(df)
mc_shower_counter = 0
mc_track_counter = 0
mc_shower_array = []
mc_track_array = []
mc_array = []
last_num = 0
for mcpart in df_mcprimary.index:
    # reset counters
    this_num = int(df_mcprimary.loc[mcpart, 'index'])
    if(this_num < last_num):
        mc_shower_array.append(mc_shower_counter)
        mc_track_array.append(mc_track_counter)
        mc_array.append(mc_shower_counter + mc_track_counter)
        mc_track_counter = 0
        mc_shower_counter = 0
    # shower events
    if(df_mcprimary.loc[mcpart, 'mcPdg'] == 11 or
       df_mcprimary.loc[mcpart, 'mcPdg'] == 22 or
       df_mcprimary.loc[mcpart, 'mcPdg'] == -11):
        mc_shower_counter = mc_shower_counter + 1
    # track events
    if(df_mcprimary.loc[mcpart, 'mcPdg'] == 2212 or
       df_mcprimary.loc[mcpart, 'mcPdg'] == 211 or
       df_mcprimary.loc[mcpart, 'mcPdg'] == -211 or
       df_mcprimary.loc[mcpart, 'mcPdg'] == 13 or
       df_mcprimary.loc[mcpart, 'mcPdg'] == -13):
        mc_track_counter = mc_track_counter + 1
    last_num = this_num

# histogram of mc tracks and showers
fig_mc_trackshower = plt.figure()
ax = fig_mc_trackshower.add_subplot(111)
mult = [mc_array, mc_shower_array, mc_track_array]
_ = plt.hist(mult, 10, (0, 10), histtype='bar', color=[
             'wheat', 'goldenrod', 'darkmagenta'], label=['Track+Shower', 'Showers', 'Tracks'])
ax.set_xlabel('Primary MC per event')
plt.legend()
plt.show()

df_pfps = getAllPfps(df_valid_pfpart)
pfp_shower_counter = 0
pfp_track_counter = 0
first_event = True
pfp_shower_array = []
pfp_track_array = []
pfp_array = []
for pfp in df_pfps.index:
    isPrimary = df_pfps.loc[pfp, 'pfoIsPrimary']
    # if(isPrimary == True):
    #     print df_pfps.loc[pfp, 'pfoPdg']
    if(df_pfps.loc[pfp, 'index'] == 0 and first_event == False):  # new event, restart counter
        pfp_shower_array.append(pfp_shower_counter)
        pfp_track_array.append(pfp_track_counter)
        pfp_array.append(pfp_shower_counter + pfp_track_counter)
        pfp_shower_counter = 0
        pfp_track_counter = 0

    first_event = False
    # shower-like events
    if(df_pfps.loc[pfp, 'pfoPdg'] == 11):
        pfp_shower_counter = pfp_shower_counter + 1
    # track-like events
    if(df_pfps.loc[pfp, 'pfoPdg'] == 13):
        pfp_track_counter = pfp_track_counter + 1

# histogram of pfp tracks and showers
fig_pfp_trackshower = plt.figure()
ax = fig_pfp_trackshower.add_subplot(111)
mult = [pfp_array, pfp_shower_array, pfp_track_array]
_ = plt.hist(mult, 8, (0, 8), histtype='bar', color=[
             'wheat', 'goldenrod', 'darkmagenta'], label=['Track+Shower', 'Showers', 'Tracks'])
ax.set_xlabel('PFP per event')
plt.legend()
plt.show()

# let's also look at the pfp accuracy

# how often are showers actually showers?
df_pfp_showers = getPfpShower(df_valid_pfpart)
mcPdg_shower = df_pfp_showers.mcPdg
pfpPdg_shower = df_pfp_showers.pfoPdg
for showers in mcPdg_shower.index:
    if(mcPdg_shower[showers] != pfpPdg_shower[showers]):
        print 'MC PDG and PFP PDG (showers) do not match!'

# how often are nues actually nues?
df_pfp_nues = getPfpNue(df_valid_pfpart)
nPfpNues = len(df_pfp_nues.index)
print 'Number of PFP Nues: ', nPfpNues
mcPdg_nue = df_pfp_nues.mcPdg
pfpPdg_nue = df_pfp_nues.pfoPdg
n_available_hits_nue = []
available_hits_nue = df_pfp_nues.nAvailableHits
for nues in mcPdg_nue.index:
    if(mcPdg_nue[nues] != pfpPdg_nue[nues]):
        print 'MC PDG and PFP PDG do not match!'
        continue
    n_available_hits_nue.append(available_hits_nue[nues])

print n_available_hits_nue
nue_efficiency = float(nPfpNues) / float(nMCNues) * 100.
print 'Nue Reco Efficiency: ', nue_efficiency

elapsed = timeit.default_timer() - start_time
print 'Time Elapsed: ', elapsed
