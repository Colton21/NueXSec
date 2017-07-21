import os
import sys
import json
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import root_numpy as rnp
from tqdm import tqdm
import timeit
from reco_true_functions import *
from cross_section_function import *


def Pandafy(fileName, tree):
    df = pd.DataFrame(rnp.root2array(fileName, tree))
    return df

# How are my neutrinos doing in these cuts? - look at index 0 events to
# get a count - make sure these are true reco and not others, how many others?

# Begin timer to measure execution duration
start_time = timeit.default_timer()
print 'Begin Selection - Time: ', start_time

# open all relevant dataframes
#infile = 'nue_matching.root'
infile = sys.argv[1]  # command line input file
infile_tree_name = 'TrueRecoMon/pandora'
infile_tree_name2 = 'NueXsec/optical_tree'
df = Pandafy(infile, infile_tree_name)
df_opt = Pandafy(infile, infile_tree_name2)

# Number of MC Particles
nEvents = len(df.index)
nEvents_opt = len(df_opt.index)
print 'Number of Entries: ', nEvents

# check if data frame is empty
if(nEvents == 0):
    print >> sys.stderr, 'Data Frame is Empty!'
    exit(1)
if(nEvents == 0):
    print >> sys.stderr, 'Data Frame Opt is Empty!'
    exit(1)


############################################################################
# let's try and remove all of the "bad" events first, before looking at reco
############################################################################

df = removeZeroMCP(df)
nEvents_MCParticle = len(df.index)
print 'Number of MC Particles: ', nEvents_MCParticle

df = removeZeroMCEng(df)
nEvents_MCEnergy = len(df.index)
print 'Number of MCParticles w/ Non-Zero MCEng: ', nEvents_MCEnergy

good_fraction = (float(nEvents_MCEnergy) / float(nEvents)) * 100.
print 'Good Fraction: ', good_fraction

##########################################################
# How many true neutrino events do we have in the sample?
#########################################################
df_mc_nu = getMCNeutrino(df)
df_mc_cc_nue = getMCCCNue(df_mc_nu)
num_mc_cc_nue = len(df_mc_cc_nue.index)

#############################
# filter events to be in time
# ############################
start_time_window = 5.0
end_time_window = 16.0
flashThreshold = 50
failEvents = inTimeList(df_opt, start_time_window,
                        end_time_window, flashThreshold)
df_opt = inTimeDataFrame(df_opt, failEvents)
df = inTimeDataFrame(df, failEvents)
print 'MC Particles within time: ', len(df.index)


###########################
# filter events to be inTPC
###########################
# df = inTPC(df)
# nEvents_inTPC = len(df.index)
# print 'MC Particles in TPC: ', nEvents_inTPC

######################################
# filter events where no reco particle
######################################
df = validPFParticle(df)
nEvents_validPFP = len(df.index)
print 'MC Particles reco as PFPs: ', nEvents_validPFP

x_min = 10
x_max = 10
y_min = 10
y_max = 50
z_min = 10
z_max = 10
df = pfpInFV(df, x_min, x_max, y_min, y_max, z_min, z_max)
nEvents_pfpInTPC = len(df.index)
print 'PFP in FV: ', nEvents_pfpInTPC

###################################
# vertex to flash information
##################################
distance = 100
df = flashRecoVtxDist(df, df_opt, distance)
print 'Passed Flash-Vtx Cut: ', len(df.index)


# some other selection cut ideas - reco energy, reco hits, shower to
# neutrino vertex, at least one reco shower, angles?

######################################
# final number of neutrinos selected
#####################################
df_nu = getPfpNeutrino(df)
# print df_nu
print 'Number of True CC Nue Neutrinos: ', num_mc_cc_nue
print 'Number of Reco Nue Neutrinos   : ', len(df_nu.index)
info_list = mcPartBreakdown(df_nu)
num_nue_cosmic = cosmicBreakdown(df)
num_pfp_cc_nue = info_list[1]
print 'Efficiency: ', float(num_pfp_cc_nue) / float(num_mc_cc_nue) * 100.
purity = float(num_pfp_cc_nue) / \
    (float(info_list[0]) + float(num_nue_cosmic)) * 100.
print 'Purity: ', purity

elapsed = timeit.default_timer() - start_time
print 'Time Elapsed: ', elapsed
