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

# Begin timer to measure execution duration
start_time = timeit.default_timer()
print ' ============================ '
print ' ==== Starting Selection ===='
print ' ============================ '
###############################
# open all relevant dataframes
###############################
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
if(df.empty):
    print >> sys.stderr, 'Data Frame DF is Empty!'
    exit(1)
if(df_opt.empty):
    print >> sys.stderr, 'Data Frame DF_OPT is Empty!'
    exit(1)


############################################################################
# let's try and remove all of the "bad" events first, before looking at reco
############################################################################
df_mc_nu = getMCNeutrino(df)
df_mc_cc_nue = getMCCCNue(df_mc_nu)
num_mc_cc_nue = len(df_mc_cc_nue.index)
printInfo(df, num_mc_cc_nue)
print 'Removing Bad Events First...'

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
df_mc_nc_nue = getMCNCNue(df_mc_nu)
num_mc_cc_nue = len(df_mc_cc_nue.index)
num_mc_nc_nue = len(df_mc_nc_nue.index)
print 'True Nue CC: ', num_mc_cc_nue
print 'True Nue NC: ', num_mc_nc_nue
##print statements##
printInfo(df, num_mc_cc_nue)

##############################
# filter events to be in time
##############################
print ' ============================ '
print ' ======= In-Time Cut: ======='
print ' ============================ '
start_time_window = 5.0
end_time_window = 16.0
flashThreshold = 50
failEvents = inTimeList(df_opt, start_time_window,
                        end_time_window, flashThreshold)
df_opt = inTimeDataFrame(df_opt, failEvents)
df = inTimeDataFrame(df, failEvents)
print 'MC Particles within time: ', len(df.index)
##print statements##
printInfo(df, num_mc_cc_nue)

######################################
# filter events where no reco particle
######################################
print ' ============================ '
print ' ====== Reco-Particle ======='
print ' ============================ '
df = validPFParticle(df)
nEvents_validPFP = len(df.index)
print 'MC Particles reco as PFPs: ', nEvents_validPFP
##print statements##
printInfo(df, num_mc_cc_nue)

###################################
# Minimum hits threshold
##################################
print ' ============================= '
print ' ====== Hits Threshold ======='
print ' ============================= '
hit_threshold = 90
df = pfpHitThreshold(df, hit_threshold)
print 'Hit Threshold, ', hit_threshold
printInfo(df, num_mc_cc_nue)

######################
# Fiducial Volume Cut
######################
print ' ============================ '
print ' ==== Fiducial Volume Cut: =='
print ' ============================ '
x_min = 0
x_max = 0
y_min = 0
y_max = 50  # this one is the top
z_min = 0
z_max = 0
df = pfpInFV(df, x_min, x_max, y_min, y_max, z_min, z_max)
nEvents_pfpInTPC = len(df.index)
print 'PFP in FV: ', nEvents_pfpInTPC
##print statements##
printInfo(df, num_mc_cc_nue)

################################
# Distance between vertex to nue
################################
print ' ============================ '
print ' ==== Shower - Nue Cut: ====='
print ' ============================ '
vtx_to_nu_distance = 50  # cm
print 'Shower to Nue Distance: ', vtx_to_nu_distance, ' cm '
df = vtxNuDistance(df, vtx_to_nu_distance)
printInfo(df, num_mc_cc_nue)


###################################
# vertex to flash information
##################################
print ' ============================ '
print ' ===== Vertex-To-Flash: ====='
print ' ============================ '
distance = 100
df = flashRecoVtxDist(df, df_opt, distance)
print 'Vertex to Flash: ', distance, ' cm '
print 'Passed Flash-Vtx Cut: ', len(df.index)


# some other selection cut ideas - reco energy, reco hits, shower to
# neutrino vertex, at least one reco shower, angles?

######################################
# final number of neutrinos selected
#####################################
print ' ============================ '
print ' ======= Final Values ======='
print ' ============================ '
printInfo(df, num_mc_cc_nue)


###############################
# cross section calculation
###############################
# flux = 3.6 * pow(10, 20) # from summing input file
# flux = 4.68 * pow(10, 21) #from 9500 * 2.4e17 nue cc / POT
# flux = 9.4 * pow(10, 20) #from 3941 nue cc
flux = 1.24 * pow(10, 21)  # from counting num nue cc and scaling
printXsec(df, num_mc_cc_nue, flux)


################
#### timer #####
################
elapsed = timeit.default_timer() - start_time
print 'Time Elapsed: ', elapsed
