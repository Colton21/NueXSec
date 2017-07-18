# Author: Colton Hill
# This file should briefly examine the reco-true information

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import pandas as pd
from collections import Counter
import root_numpy as rnp
from tqdm import tqdm
import timeit
from reco_true_functions import *

####################################################################
# Start adding the muon stuff! Also add a variable to differentiate
# between cosmic and beam file configs
# finish adding savefig's
####################################################################

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
# infile = 'nue_xsec_extraction_2.root'
# infile = 'nue_matching.root'
# infile = 'nue_cosmic_matching.root'
infile = sys.argv[1]  # give path to file
cosmic_file = sys.argv[2]  # give true or false if cosmics in datafile
infile_tree_name = 'TrueRecoMon/pandora'
df = Pandafy(infile, infile_tree_name)

# Number of MC Particles
nEvents = len(df.index)
print 'Number of MC Particles: %s' % (nEvents)

# check if data frame is empty
if(nEvents == 0):
    print >> sys.stderr, 'Data Frame is Null!'
    exit(1)
# for guidance of data product names
print df.head(1)

############################################################################
# let's try and remove all of the "bad" events first, before looking at reco
############################################################################
df = removeZeroMCP(df)
nEvents_MCParticle = len(df.index)
print 'Number of MCParticles: ', nEvents_MCParticle

df = removeZeroMCEng(df)
nEvents_MCEnergy = len(df.index)
print 'Number of MCParticles w/ Non-Zero MCEng: ', nEvents_MCEnergy

good_fraction = (float(nEvents_MCEnergy) / float(nEvents)) * 100.
print 'Good Fraction: ', good_fraction

###########################
# filter events to be inTPC
###########################
# df = inTPC(df)
nEvents_inTPC = len(df.index)
print 'MC Particles in TPC: ', nEvents_inTPC
###################################################
# let's look at some characteristics of the dataset
###################################################
# zeroIndex should return true MC neutrino objects
# as they are all stored at index 0 per event
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
print 'Dropping all non-reconstructed all objects: '
print 'Remaining MC Particles: ', len(df_valid_pfpart)


# let's get just the PFP Neutrino candidates
df_pfp_neutrinos = getPfpNeutrino(df_valid_pfpart)
# what is the number of pfp neutrinos per event?
df_pfp_neutrino_pdgs = df_zero.nNeutrinoPfos
# fig_pfpnu = plt.figure()
# ax = fig_pfpnu.add_subplot(111)
df_pfp_neutrino_pdgs.plot.hist(alpha=0.5)
# ax.set_xlabel('PFP Neutrino per event')

####################################
# compare MC showers/tracks per event
###################################
df_mcprimary = getMCPrimary(df)
mc_shower_counter = 0
mc_track_counter = 0
mc_shower_array = []
mc_track_array = []
mc_array = []
last_num = 0
mc_electron_counter = 0
for mcpart in tqdm(df_mcprimary.index):
    # reset counters
    this_num = int(df_mcprimary.loc[mcpart, 'index'])
    if(this_num < last_num):
        mc_shower_array.append(mc_shower_counter)
        mc_track_array.append(mc_track_counter)
        mc_array.append(mc_shower_counter + mc_track_counter)
        mc_track_counter = 0
        mc_shower_counter = 0
    # shower events
    if(df_mcprimary.loc[mcpart, 'mcPdg'] == 11):
        mc_electron_counter = mc_electron_counter + 1
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
fig_mc_trackshower.savefig('mc_per_event_track_shower.pdf')
plt.close()

########################################
# compare pfp showers / tracks per event
########################################
df_pfps = getAllPfps(df_valid_pfpart)
pfp_shower_counter = 0
pfp_track_counter = 0
first_event = True
pfp_shower_array = []
pfp_track_array = []
pfp_array = []
for pfp in tqdm(df_pfps.index):
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
fig_pfp_trackshower.savefig('pfp_per_event_track_shower.pdf')
plt.close()

#####################################
# let's also look at the pfp accuracy
####################################
# how often are showers actually showers?
nPfpShowers = 0
df_pfp_showers = getPfpShower(df_valid_pfpart)
mcPdg_showers = df_pfp_showers.mcPdg
pfpPdg_showers = df_pfp_showers.pfoPdg
available_hits_shower = df_pfp_showers.nAvailableHits
pfo_hits_shower = df_pfp_showers.nPfoHits
mc_hits_shower = df_pfp_showers.nMCHits
false_pfp_counter = 0
true_pfp_counter = 0
photon_pfp_counter = 0
cosmic_pfp_counter = 0


# energy
energy_mc_shwr = df_pfp_showers.mcMomentum
energy_pfp_shwr = df_pfp_showers.pfoEnergy
mcEnergy_shwr = []
mcEnergy_elec = []
mcEnergy_notElec = []
mcEnergy_prot = []
mcEnergy_gamma = []
mcEnergy_pion = []
mcEnergy_neut = []
mcEnergy_cosmic = []
pfpEnergy_shwr = []
pfpEnergy_elec = []
pfpEnergy_notElec = []
pfpEnergy_prot = []
pfpEnergy_gamma = []
pfpEnergy_pion = []
pfpEnergy_neut = []
pfpEnergy_cosmic = []
diff_energy_shwr = []
diff_energy_elec = []
diff_energy_notElec = []
diff_energy_prot = []
diff_energy_pion = []
diff_energy_neut = []
diff_energy_gamma = []
diff_energy_cosmic = []
# opening angle
openangle_pfp_shwr = df_pfp_showers.pfoOpenAngle
pfpOpenAngle_shwr = []
pfpOpenAngle_elec = []
pfpOpenAngle_notElec = []
# shower Length
pfp_shower_length = df_pfp_showers.pfoLength
mc_shower_length = df_pfp_showers.mcLength
mcLength_shwr = []
mcLength_elec = []
mcLength_prot = []
mcLength_gamma = []
mcLength_neut = []
mcLength_pion = []
mcLength_cosmic = []
pfpLength_shwr = []
pfpLength_elec = []
pfpLength_prot = []
pfpLength_gamma = []
pfpLength_neut = []
pfpLength_pion = []
pfpLength_cosmic = []
# vertex
mcVtxX_shwr = df_pfp_showers.mcVtxX
mcVtxY_shwr = df_pfp_showers.mcVtxY
mcVtxZ_shwr = df_pfp_showers.mcVtxZ
pfpVtxX_shwr = df_pfp_showers.pfoVtxX
pfpVtxY_shwr = df_pfp_showers.pfoVtxY
pfpVtxZ_shwr = df_pfp_showers.pfoVtxZ
vtx_diff_X_shwr = []
vtx_diff_Y_shwr = []
vtx_diff_Z_shwr = []
vtx_diff_shwr = []
vtx_diff_X_elec = []
vtx_diff_Y_elec = []
vtx_diff_Z_elec = []
vtx_diff_elec = []
vtx_Z_pfp_shwr = []
vtx_X_pfp_elec = []
vtx_X_pfp_prot = []
vtx_X_pfp_gamma = []
vtx_X_pfp_cosmic = []
vtx_Y_pfp_elec = []
vtx_Y_pfp_prot = []
vtx_Y_pfp_gamma = []
vtx_Y_pfp_cosmic = []
vtx_Z_pfp_elec = []
vtx_Z_pfp_prot = []
vtx_Z_pfp_gamma = []
vtx_Z_pfp_cosmic = []
vtx_Z_mc_shwr = []
vtx_Z_mc_elec = []
vtx_Z_mc_prot = []
vtx_Z_mc_gamma = []
# direction
mcDirX_shwr = df_pfp_showers.mcDirX
mcDirY_shwr = df_pfp_showers.mcDirY
mcDirZ_shwr = df_pfp_showers.mcDirZ
pfpDirX_shwr = df_pfp_showers.pfoDirX
pfpDirY_shwr = df_pfp_showers.pfoDirY
pfpDirZ_shwr = df_pfp_showers.pfoDirZ
dir_diff_X_shwr = []
dir_diff_Y_shwr = []
dir_diff_Z_shwr = []
dir_diff_shwr = []
dir_diff_X_elec = []
dir_diff_Y_elec = []
dir_diff_Z_elec = []
dir_diff_elec = []
# interaciton modes
mcEnergy_ccqe = []
mcEnergy_res = []
mcEnergy_dis = []
mcEnergy_coh = []
pfpEnergy_ccqe = []
pfpEnergy_res = []
pfpEnergy_dis = []
pfpEnergy_coh = []
interaction_list = []
mcEnergy_ccqe_elec = []
mcEnergy_res_elec = []
mcEnergy_dis_elec = []
mcEnergy_coh_elec = []
pfpEnergy_ccqe_elec = []
pfpEnergy_res_elec = []
pfpEnergy_dis_elec = []
pfpEnergy_coh_elec = []
interaction_list_elec = []
# hits
particle_list = []
pion_pfp_hits = []
pion_mc_hits = []
pion_ratio_hits = []
neutron_pfp_hits = []
neutron_mc_hits = []
neutron_ratio_hits = []
proton_pfp_hits = []
proton_mc_hits = []
proton_ratio_hits = []
photon_pfp_hits = []
photon_mc_hits = []
photon_ratio_hits = []
electron_pfp_hits = []
electron_mc_hits = []
electron_ratio_hits = []
cosmic_pfp_hits = []
cosmic_mc_hits = []
cosmic_ratio_hits = []
# completeness
# nMatchedHits / nMCHits
pfp_completeness = df_pfp_showers.completeness
pion_completeness = []
neutron_completeness = []
proton_completeness = []
photon_completeness = []
electron_completeness = []
cosmic_completeness = []
ccqe_completeness = []
res_completeness = []
dis_completeness = []
coh_completeness = []
# Begin looping showers
for showers in tqdm(mcPdg_showers.index):
    nPfpShowers = nPfpShowers + 1
    mcPdg_shower = mcPdg_showers[showers]
    pfpPdg_shower = pfpPdg_showers[showers]
    pfpHits = pfo_hits_shower[showers]
    mcHits = mc_hits_shower[showers]
    mcMode = df_pfp_showers.mcMode[showers]
    this_completeness = pfp_completeness[showers]
    # let's do some true mc comparisons
    # vertex
    x_diff = pfpVtxX_shwr[showers] - mcVtxX_shwr[showers]
    y_diff = pfpVtxY_shwr[showers] - mcVtxY_shwr[showers]
    z_diff = pfpVtxZ_shwr[showers] - mcVtxZ_shwr[showers]
    total_diff = np.sqrt((x_diff * x_diff) +
                         (y_diff * y_diff) + (z_diff * z_diff))
    vtx_diff_X_shwr.append(x_diff)
    vtx_diff_Y_shwr.append(y_diff)
    vtx_diff_Z_shwr.append(z_diff)
    vtx_diff_shwr.append(total_diff)
    vtx_Z_pfp_shwr.append(pfpVtxZ_shwr[showers])
    vtx_Z_mc_shwr.append(mcVtxZ_shwr[showers])
    # direction
    dir_x_diff_shwr = pfpDirX_shwr[showers] - mcDirX_shwr[showers]
    dir_y_diff_shwr = pfpDirY_shwr[showers] - mcDirY_shwr[showers]
    dir_z_diff_shwr = pfpDirZ_shwr[showers] - mcDirZ_shwr[showers]
    dir_dot_shwr = np.arccos(np.dot([pfpDirX_shwr[showers], pfpDirY_shwr[showers], pfpDirZ_shwr[showers]], [
        mcDirX_shwr[showers], mcDirY_shwr[showers], mcDirZ_shwr[showers]])) * (180 / 3.1415)
    # print 'PFP', pfpDirX_shwr[showers], pfpDirY_shwr[showers], pfpDirZ_shwr[showers]
    # print 'MC ', mcDirX_shwr[showers], mcDirY_shwr[showers],
    dir_diff_X_shwr.append(dir_x_diff_shwr)
    dir_diff_Y_shwr.append(dir_y_diff_shwr)
    dir_diff_Z_shwr.append(dir_z_diff_shwr)
    dir_diff_shwr.append(dir_dot_shwr)
    # Length
    mc_length = mc_shower_length[showers]
    pfp_length = pfp_shower_length[showers]
    mcLength_shwr.append(mc_length)
    pfpLength_shwr.append(pfp_length)

    # energy
    mc_energy_shower = energy_mc_shwr[showers]
    pfp_energy_shower = energy_pfp_shwr[showers]
    if(pfp_energy_shower < 0):
        pfp_energy_shower = 0
    mcEnergy_shwr.append(mc_energy_shower)
    pfpEnergy_shwr.append(pfp_energy_shower)
    diff_energy_shwr.append(pfp_energy_shower - mc_energy_shower)
    if(mcPdg_shower != pfpPdg_shower and mcPdg_shower != 22):
        false_pfp_counter = false_pfp_counter + 1
        mcEnergy_notElec.append(mc_energy_shower)
        pfpEnergy_notElec.append(pfp_energy_shower)
        diff_energy_notElec.append(pfp_energy_shower - mc_energy_shower)
        if(mcPdg_shower == 211 or mcPdg_shower == -211):
            particle_list.append('Pion')
            pion_pfp_hits.append(pfpHits)
            pion_mc_hits.append(mcHits)
            pion_ratio_hits.append(float(pfpHits) / float(mcHits))
            mcEnergy_pion.append(mc_energy_shower)
            pfpEnergy_pion.append(pfp_energy_shower)
            diff_energy_pion.append(pfp_energy_shower - mc_energy_shower)
            pion_completeness.append(this_completeness)
            pfpLength_pion.append(pfp_length)
            mcLength_pion.append(mc_length)
        if(mcPdg_shower == 2112):
            particle_list.append('Neut')
            neutron_pfp_hits.append(pfpHits)
            neutron_mc_hits.append(mcHits)
            neutron_ratio_hits.append(float(pfpHits) / float(mcHits))
            mcEnergy_neut.append(mc_energy_shower)
            pfpEnergy_neut.append(pfp_energy_shower)
            diff_energy_neut.append(pfp_energy_shower - mc_energy_shower)
            neutron_completeness.append(this_completeness)
            pfpLength_neut.append(pfp_length)
            mcLength_neut.append(mc_length)
        if(mcPdg_shower == 2212):
            particle_list.append('Proton')
            proton_pfp_hits.append(pfpHits)
            proton_mc_hits.append(mcHits)
            proton_ratio_hits.append(float(pfpHits) / float(mcHits))
            mcEnergy_prot.append(mc_energy_shower)
            pfpEnergy_prot.append(pfp_energy_shower)
            diff_energy_prot.append(pfp_energy_shower - mc_energy_shower)
            vtx_X_pfp_prot.append(pfpVtxX_shwr[showers])
            vtx_Y_pfp_prot.append(pfpVtxY_shwr[showers])
            vtx_Z_pfp_prot.append(pfpVtxZ_shwr[showers])
            proton_completeness.append(this_completeness)
            pfpLength_prot.append(pfp_length)
            mcLength_prot.append(mc_length)

    # if mcNuPdg == 0 then it's not from a neutrino!
    if(cosmic_file == 'True'):
        if(df_pfp_showers.mcNuPdg[showers] == 0):
            particle_list.append('Cosmic')
            cosmic_pfp_counter = cosmic_pfp_counter + 1
            mcEnergy_cosmic.append(mc_energy_shower)
            pfpEnergy_cosmic.append(pfp_energy_shower)
            diff_energy_cosmic.append(pfp_energy_shower - mc_energy_shower)
            cosmic_pfp_hits.append(pfpHits)
            cosmic_mc_hits.append(mcHits)
            cosmic_ratio_hits.append(float(pfpHits) / float(mcHits))
            vtx_X_pfp_cosmic.append(pfpVtxX_shwr[showers])
            vtx_Y_pfp_cosmic.append(pfpVtxY_shwr[showers])
            vtx_Z_pfp_cosmic.append(pfpVtxZ_shwr[showers])
            cosmic_completeness.append(this_completeness)
            pfpLength_cosmic.append(pfp_length)
            mcLength_cosmic.append(mc_length)
    if(mcPdg_shower == 22):
        photon_pfp_counter = photon_pfp_counter + 1
        mcEnergy_notElec.append(mc_energy_shower)
        pfpEnergy_notElec.append(pfp_energy_shower)
        diff_energy_notElec.append(pfp_energy_shower - mc_energy_shower)
        particle_list.append('Photon')
        photon_pfp_hits.append(pfpHits)
        photon_mc_hits.append(mcHits)
        photon_ratio_hits.append(float(pfpHits) / float(mcHits))
        mcEnergy_gamma.append(mc_energy_shower)
        pfpEnergy_gamma.append(pfp_energy_shower)
        diff_energy_gamma.append(pfp_energy_shower - mc_energy_shower)
        vtx_X_pfp_gamma.append(pfpVtxX_shwr[showers])
        vtx_Y_pfp_gamma.append(pfpVtxY_shwr[showers])
        vtx_Z_pfp_gamma.append(pfpVtxZ_shwr[showers])
        photon_completeness.append(this_completeness)
        pfpLength_gamma.append(pfp_length)
        mcLength_gamma.append(mc_length)

    if(mcPdg_shower == pfpPdg_shower):
        # print 'Electron Num Unmatched Hits: ', available_hits_shower[showers]
        # print 'Electron: ', pfo_hits_shower[showers]
        true_pfp_counter = true_pfp_counter + 1
        particle_list.append('Elec')
        electron_pfp_hits.append(pfpHits)
        electron_mc_hits.append(mcHits)
        electron_ratio_hits.append(float(pfpHits) / float(mcHits))
        elec_x_diff = pfpVtxX_shwr[showers] - mcVtxX_shwr[showers]
        elec_y_diff = pfpVtxY_shwr[showers] - mcVtxY_shwr[showers]
        elec_z_diff = pfpVtxZ_shwr[showers] - mcVtxZ_shwr[showers]
        elec_total_diff = np.sqrt((elec_x_diff * elec_x_diff) +
                                  (elec_y_diff * elec_y_diff) + (elec_z_diff * elec_z_diff))
        vtx_diff_X_elec.append(elec_x_diff)
        vtx_diff_Y_elec.append(elec_y_diff)
        vtx_diff_Z_elec.append(elec_z_diff)
        vtx_diff_elec.append(elec_total_diff)
        vtx_X_pfp_elec.append(pfpVtxX_shwr[showers])
        vtx_Y_pfp_elec.append(pfpVtxY_shwr[showers])
        vtx_Z_pfp_elec.append(pfpVtxZ_shwr[showers])
        vtx_Z_mc_elec.append(mcVtxZ_shwr[showers])
        mcEnergy_elec.append(mc_energy_shower)
        pfpEnergy_elec.append(pfp_energy_shower)
        diff_energy_elec.append(pfp_energy_shower - mc_energy_shower)
        dir_diff_X_elec.append(dir_x_diff_shwr)
        dir_diff_Y_elec.append(dir_y_diff_shwr)
        dir_diff_Z_elec.append(dir_z_diff_shwr)
        dir_diff_elec.append(dir_dot_shwr)
        electron_completeness.append(this_completeness)
        # Length
        mcLength_elec.append(mc_length)
        pfpLength_elec.append(pfp_length)

        # what do the interaction modes look like for the true electrons?
        # ccqe /cc0pi - is it actually 0 pi?
        if(mcMode == 0):
            interaction_list_elec.append('CCQE')
            mcEnergy_ccqe_elec.append(mc_energy_shower)
            pfpEnergy_ccqe_elec.append(pfp_energy_shower)
        # res
        if(mcMode == 1):
            interaction_list_elec.append('Res')
            mcEnergy_res_elec.append(mc_energy_shower)
            pfpEnergy_res_elec.append(pfp_energy_shower)
        # dis
        if(mcMode == 2):
            interaction_list_elec.append('DIS')
            mcEnergy_dis_elec.append(mc_energy_shower)
            pfpEnergy_dis_elec.append(pfp_energy_shower)
        # coh
        if(mcMode == 3):
            interaction_list_elec.append('Coh')
            mcEnergy_coh_elec.append(mc_energy_shower)
            pfpEnergy_coh_elec.append(pfp_energy_shower)

    # I'm also going to compare some of the interaction modes in the loop
    # ccqe /cc0pi - is it actually 0 pi?
    if(mcMode == 0):
        interaction_list.append('CCQE')
        mcEnergy_ccqe.append(mc_energy_shower)
        pfpEnergy_ccqe.append(pfp_energy_shower)
        ccqe_completeness.append(this_completeness)
    # res
    if(mcMode == 1):
        interaction_list.append('Res')
        mcEnergy_res.append(mc_energy_shower)
        pfpEnergy_res.append(pfp_energy_shower)
        res_completeness.append(this_completeness)
    # dis
    if(mcMode == 2):
        interaction_list.append('DIS')
        mcEnergy_dis.append(mc_energy_shower)
        pfpEnergy_dis.append(pfp_energy_shower)
        dis_completeness.append(this_completeness)
    # coh
    if(mcMode == 3):
        interaction_list.append('Coh')
        mcEnergy_coh.append(mc_energy_shower)
        pfpEnergy_coh.append(pfp_energy_shower)
        coh_completeness.append(this_completeness)


print 'Showers Reconstructed: ', len(mcPdg_showers.index)
print 'True Electron Showers: ', true_pfp_counter
print 'Not Showers          : ', false_pfp_counter
print 'True Photon Showers  : ', photon_pfp_counter
print 'True Cosmic Showers  : ', cosmic_pfp_counter
# what are the particles most reconstructed as a shower?
particle_counts = Counter(particle_list)
temp_plot_df = pd.DataFrame.from_dict(particle_counts, orient='index')
# fig_pfp_showers = plt.figure()
# ax = fig_pfp_showers.add_subplot(111)
temp_plot_df.plot(kind='bar', alpha=0.5, legend=False)
# plt.show()
# ax.set_xlabel('True Particle Reconstructed as PFP Shower')

interaction_counts = Counter(interaction_list)
temp_plot_df = pd.DataFrame.from_dict(interaction_counts, orient='index')
# fig_pfp_showers = plt.figure()
# ax = fig_pfp_showers.add_subplot(111)
temp_plot_df.plot(kind='bar', alpha=0.5, legend=False)
# plt.show()
# ax.set_xlabel('Interaction Mode - Reco PFP Shower')

interaction_counts_elec = Counter(interaction_list_elec)
temp_plot_df = pd.DataFrame.from_dict(interaction_counts_elec, orient='index')
# fig_pfp_showers = plt.figure()
# ax = fig_pfp_showers.add_subplot(111)
temp_plot_df.plot(kind='bar', alpha=0.5, legend=False)
# plt.show()
# ax.set_xlabel('Interaction Mode - Reco PFP Electrons')

# histogram of mc shower hits
fig_mc_shower_hits = plt.figure()
ax = fig_mc_shower_hits.add_subplot(111)
if(cosmic_file == 'False'):
    mult_mc = [pion_mc_hits, neutron_mc_hits, proton_mc_hits,
               photon_mc_hits, electron_mc_hits]
    _ = plt.hist(mult_mc, 30, (0, 3500), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_mc = [pion_mc_hits, neutron_mc_hits, proton_mc_hits,
               photon_mc_hits, electron_mc_hits, cosmic_mc_hits]
    _ = plt.hist(mult_mc, 30, (0, 3500), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('True MC Hits')
plt.legend()
fig_mc_shower_hits.savefig('true_hits_type.pdf')
plt.close()

# histogram of pfp shower hits
fig_pfp_shower_hits = plt.figure()
ax = fig_pfp_shower_hits.add_subplot(111)
if(cosmic_file == 'False'):
    mult_pfp = [pion_pfp_hits, neutron_pfp_hits, proton_pfp_hits,
                photon_pfp_hits, electron_pfp_hits]
    _ = plt.hist(mult_pfp, 30, (0, 3500), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_pfp = [pion_pfp_hits, neutron_pfp_hits, proton_pfp_hits,
                photon_pfp_hits, electron_pfp_hits, cosmic_pfp_hits]
    _ = plt.hist(mult_pfp, 30, (0, 3500), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('Reconstructed Hits')
plt.legend()
fig_pfp_shower_hits.savefig('reco_hits_type.pdf')
plt.close()

# histogram of ratio pfp/MC shower hits
fig_pfp_shower_hits = plt.figure()
ax = fig_pfp_shower_hits.add_subplot(111)
if(cosmic_file == 'False'):
    mult_ratio = [pion_ratio_hits, neutron_ratio_hits,
                  proton_ratio_hits, photon_ratio_hits, electron_ratio_hits]
    _ = plt.hist(mult_ratio, 30, (0, 4), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_ratio = [pion_ratio_hits, neutron_ratio_hits,
                  proton_ratio_hits, photon_ratio_hits, electron_ratio_hits, cosmic_ratio_hits]
    _ = plt.hist(mult_ratio, 30, (0, 4), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('Ratio Reco / MC Hits')
plt.legend()
fig_pfp_shower_hits.savefig('ratio_hits_type.pdf')
plt.close()

# energy per particle type vs num showers per event

####################################
# 2d histogram of ratio vs true hits
####################################
fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(electron_mc_hits, electron_ratio_hits,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Electron True MC Hits')
ax.set_ylabel('Electon Ratio Reco / MC Hits')
plt.colorbar()
fig_shower_hits_2d.savefig('elec_ratio_true_hits.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(pion_mc_hits, pion_ratio_hits,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Pion True MC Hits')
ax.set_ylabel('Pion Ratio Reco / MC Hits')
plt.colorbar()
fig_shower_hits_2d.savefig('pion_ratio_true_hits.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(neutron_mc_hits, neutron_ratio_hits,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Neutron True MC Hits')
ax.set_ylabel('Neutron Ratio Reco / MC Hits')
plt.colorbar()
fig_shower_hits_2d.savefig('neutron_ratio_true_hits.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(proton_mc_hits, proton_ratio_hits,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Proton True MC Hits')
ax.set_ylabel('Proton Ratio Reco / MC Hits')
plt.colorbar()
fig_shower_hits_2d.savefig('proton_ratio_hits_true_hits.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(photon_mc_hits, photon_ratio_hits,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Photon True MC Hits')
ax.set_ylabel('Photon Ratio Reco / MC Hits')
plt.colorbar()
fig_shower_hits_2d.savefig('photon_ratio_hits_true_hits.pdf')
plt.close()

if(cosmic_file == 'True'):
    fig_shower_hits_2d = plt.figure()
    ax = fig_shower_hits_2d.add_subplot(111)
    _ = plt.hist2d(cosmic_mc_hits, cosmic_ratio_hits,
                   20, cmap=cm.summer, norm=LogNorm())
    ax.set_xlabel('Cosmic True MC Hits')
    ax.set_ylabel('Cosmic Ratio Reco / MC Hits')
    plt.colorbar()
    fig_shower_hits_2d.savefig('cosmic_ratio_hits_true_hits.pdf')
    plt.close()

# 2d histogram of ratio vs true hits
fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(electron_pfp_hits, electron_ratio_hits,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Electron Reco Hits')
ax.set_ylabel('Electron Ratio Reco / MC Hits')
plt.colorbar()
fig_shower_hits_2d.savefig('electron_ratio_hits_reco_hits.pdf')
plt.close()

# histograms for energy
fig_energy_mc = plt.figure()
ax = fig_energy_mc.add_subplot(111)
if(cosmic_file == 'False'):
    mult_eng_mc = [mcEnergy_pion, mcEnergy_neut,
                   mcEnergy_prot, mcEnergy_gamma, mcEnergy_elec]
    _ = plt.hist(mult_eng_mc, 80, (0, 2), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_eng_mc = [mcEnergy_pion, mcEnergy_neut,
                   mcEnergy_prot, mcEnergy_gamma, mcEnergy_elec, mcEnergy_cosmic]
    _ = plt.hist(mult_eng_mc, 80, (0, 2), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('True MC Shower Momentum [GeV]')
plt.legend()
fig_energy_mc.savefig('true_shower_momentum_type.pdf')
plt.close()

fig_energy_pfp = plt.figure()
ax = fig_energy_pfp.add_subplot(111)
if(cosmic_file == 'False'):
    mult_eng_pfp = [pfpEnergy_pion, pfpEnergy_neut,
                    pfpEnergy_prot, pfpEnergy_gamma, pfpEnergy_elec]
    _ = plt.hist(mult_eng_pfp, 40, (0, 2), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_eng_pfp = [pfpEnergy_pion, pfpEnergy_neut,
                    pfpEnergy_prot, pfpEnergy_gamma, pfpEnergy_elec, pfpEnergy_cosmic]
    _ = plt.hist(mult_eng_pfp, 40, (0, 2), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('Reco Shower Momentum [GeV]')
plt.legend()
fig_energy_pfp.savefig('reco_shower_momentum_type.pdf')
plt.close()

fig_energy_diff = plt.figure()
ax = fig_energy_diff.add_subplot(111)
if(cosmic_file == 'False'):
    mult_eng_diff = [diff_energy_pion, diff_energy_neut,
                     diff_energy_prot, diff_energy_gamma, diff_energy_elec]
    _ = plt.hist(mult_eng_diff, 40, (-2, 2), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_eng_diff = [diff_energy_pion, diff_energy_neut,
                     diff_energy_prot, diff_energy_gamma, diff_energy_elec, diff_energy_cosmic]
    _ = plt.hist(mult_eng_diff, 40, (-2, 2), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('Reco - True Shower Momentum [GeV]')
plt.legend()
fig_energy_diff.savefig('reco-true_shower_momentum_type.pdf')
plt.close()

fig_energy_diff = plt.figure()
ax = fig_energy_diff.add_subplot(111)
if(cosmic_file == 'False'):
    mult_eng_diff = [pfpLength_pion, pfpLength_neut,
                     pfpLength_prot, pfpLength_gamma, pfpLength_elec]
    _ = plt.hist(mult_eng_diff, 40, (0, 100), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_eng_diff = [pfpLength_pion, pfpLength_neut,
                     pfpLength_prot, pfpLength_gamma, pfpLength_elec, pfpLength_cosmic]
    _ = plt.hist(mult_eng_diff, 40, (0, 100), histtype='bar', fill=True, stacked=True,
                 color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('Reco Shower Length [cm]')
plt.legend()
fig_energy_diff.savefig('reco_shower_length_type.pdf')
plt.close()

fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(mcEnergy_elec, diff_energy_elec,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('True Electron Momentum [GeV]')
ax.set_ylabel('Reco - True Electron Momentum [Gev]')
plt.colorbar()
fig_shower_eng_2d.savefig('reco-true_electron_energy_true_energy.pdf')
plt.close()

fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(mcEnergy_shwr, diff_energy_shwr,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('True Shower Momentum [GeV]')
ax.set_ylabel('Reco - True Shower Momentum [Gev]')
plt.colorbar()
fig_shower_eng_2d.savefig('reco-true_shower_energy_true_energy.pdf')
plt.close()
#####
fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(pfpLength_elec, mcEnergy_elec,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Electron Reco Shower Length [cm]')
ax.set_ylabel('True Electron Momentum [Gev]')
plt.colorbar()
fig_shower_eng_2d.savefig('true_electron_momentum_reco_shower_length.pdf')
plt.close()

fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(pfpLength_shwr, mcEnergy_shwr,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('All Reco Shower Length [cm]')
ax.set_ylabel('True Shower Momentum [GeV]')
plt.colorbar()
fig_shower_eng_2d.savefig('true_shower_momentum_reco_shower_length.pdf')
plt.close()

fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(pfpLength_elec, pfpEnergy_elec,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Electron Reco Shower Length [cm]')
ax.set_ylabel('Reco Electron Momentum [Gev]')
plt.colorbar()
fig_shower_eng_2d.savefig('reco_electron_momentum_reco_shower_length.pdf')
plt.close()

fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(pfpLength_shwr, pfpEnergy_shwr,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('All Reco Shower Length [cm]')
ax.set_ylabel('Reco Shower Momentum [GeV]')
plt.colorbar()
fig_shower_eng_2d.savefig('reco_shower_momentum_reco_shower_length.pdf')
plt.close()

fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(pfpLength_elec, diff_energy_elec,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Electron Reco-True Shower Length [cm]')
ax.set_ylabel('Reco Electron Momentum [Gev]')
plt.colorbar()
fig_shower_eng_2d.savefig('reco-true_electron_momentum_reco_shower_length.pdf')
plt.close()

fig_shower_eng_2d = plt.figure()
ax = fig_shower_eng_2d.add_subplot(111)
_ = plt.hist2d(pfpLength_shwr, diff_energy_shwr,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('All Reco-True Shower Length [cm]')
ax.set_ylabel('Reco Shower Momentum [GeV]')
plt.colorbar()
fig_shower_eng_2d.savefig('reco-true_shower_momentum_reco_shower_length.pdf')
plt.close()

# histograms for energy in terms of interaction modes
fig_energy_mode = plt.figure()
ax = fig_energy_mode.add_subplot(111)
mult_eng_mc_mode = [mcEnergy_ccqe, mcEnergy_res,
                    mcEnergy_dis, mcEnergy_coh]
_ = plt.hist(mult_eng_mc_mode, 40, (0, 2), histtype='bar', fill=True, stacked=True,
             color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue'], label=['CCQE', 'Resonant', 'DIS', 'Coherent'])
ax.set_xlabel('True Shower Momentum [GeV]')
plt.legend()
fig_energy_mode.savefig('true_shower_momentum_mode.pdf')
plt.close()

fig_energy_mode = plt.figure()
ax = fig_energy_mode.add_subplot(111)
mult_eng_pfp_mode = [pfpEnergy_ccqe, pfpEnergy_res,
                     pfpEnergy_dis, pfpEnergy_coh]
_ = plt.hist(mult_eng_pfp_mode, 40, (0, 2), histtype='bar', fill=True, stacked=True,
             color=['wheat', 'goldenrod', 'darkmagenta', 'skyblue'], label=['CCQE', 'Resonant', 'DIS', 'Coherent'])
ax.set_xlabel('Reco Shower Momentum [GeV]')
plt.legend()
fig_energy_mode.savefig('reco_shower_momentum_mode.pdf')
plt.close()

# vtx and energy histograms
# fig_shower_hits_2d = plt.figure()
# ax = fig_shower_hits_2d.add_subplot(111)
# _ = plt.hist2d(pfpEnergy_shwr, vtx_diff_shwr,
#                20, cmap=cm.summer, norm=LogNorm())
# ax.set_xlabel('Reco Shower Energy')
# ax.set_ylabel('Reco - True Shower Vertex')
# plt.colorbar()
# fig_shower_hits_2d.savefig('reco-true_shower_vtx_reco_energy.pdf')
# plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(mcEnergy_shwr, vtx_diff_shwr,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('True Shower Energy')
ax.set_ylabel('Reco - True Shower Vertex')
plt.colorbar()
fig_shower_hits_2d.savefig('reco-true_shower_vtx_true_energy.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(vtx_Z_pfp_shwr, vtx_diff_shwr,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Reco Shower Vtx Z')
ax.set_ylabel('Reco - True Shower Vertex')
plt.colorbar()
fig_shower_hits_2d.savefig('reco-true_shower_vtx_reco_vtxZ.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(vtx_Z_mc_shwr, vtx_diff_shwr,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('True Shower Vtx Z')
ax.set_ylabel('Reco - True Shower Vertex')
plt.colorbar()
fig_shower_hits_2d.savefig('reco-true_shower_vtx_true_vtxZ.pdf')
plt.close()

# fig_shower_hits_2d = plt.figure()
# ax = fig_shower_hits_2d.add_subplot(111)
# _ = plt.hist2d(pfpEnergy_elec, vtx_diff_elec,
#                20, cmap=cm.summer, norm=LogNorm())
# ax.set_xlabel('Reco Electron Energy')
# ax.set_ylabel('Reco - True Electron Vertex')
# plt.colorbar()
# fig_shower_hits_2d.savefig('reco-true_electron_vtx_reco_energy.pdf')
# plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(mcEnergy_elec, vtx_diff_elec,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('True Electron Energy')
ax.set_ylabel('Reco - True Electron Vertex')
plt.colorbar()
fig_shower_hits_2d.savefig('reco-true_electron_vtx_true_energy.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(vtx_Z_pfp_elec, vtx_diff_elec,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('Reco Electron Vtx Z')
ax.set_ylabel('Reco - True Electron Vertex')
plt.colorbar()
fig_shower_hits_2d.savefig('reco-true_electron_vtx_reco_vtxZ.pdf')
plt.close()

fig_shower_hits_2d = plt.figure()
ax = fig_shower_hits_2d.add_subplot(111)
_ = plt.hist2d(vtx_Z_mc_elec, vtx_diff_elec,
               20, cmap=cm.summer, norm=LogNorm())
ax.set_xlabel('True Electron Vtx Z')
ax.set_ylabel('Reco - True Electron Vertex')
plt.colorbar()
fig_shower_hits_2d.savefig('reco-true_electron_vtx_true_vtxZ.pdf')
plt.close()

# histograms for vtx
fig_vtx_diff_shwr = plt.figure()
ax = fig_vtx_diff_shwr.add_subplot(111)
mult_vtx = [vtx_diff_X_shwr, vtx_diff_Y_shwr,
            vtx_diff_Z_shwr]
_ = plt.hist(mult_vtx, 80, (-80, 80), histtype='bar', fill=True, stacked=True, color=[
             'tomato', 'goldenrod', 'darkmagenta'], label=['X', 'Y', 'Z'])
ax.set_xlabel('Reco - True Shower Vertex [cm]')
plt.legend()
fig_vtx_diff_shwr.savefig('reco-true_shower_vtx.pdf')
plt.close()

# histograms for vtx, true electron only
fig_vtx_diff_elec = plt.figure()
ax = fig_vtx_diff_elec.add_subplot(111)
mult_vtx = [vtx_diff_X_elec, vtx_diff_Y_elec,
            vtx_diff_Z_elec]
_ = plt.hist(mult_vtx, 80, (-80, 80), histtype='bar', fill=True, stacked=True, color=[
             'tomato', 'goldenrod', 'darkmagenta'], label=['X', 'Y', 'Z'])
ax.set_xlabel('Reco - True Shower (Electron) Vertex [cm]')
plt.legend()
fig_vtx_diff_elec.savefig('reco-true_electron_vtx.pdf')
plt.close()

# histograms for vtx - total distances
fig_vtx_diff_shwr_total = plt.figure()
ax = fig_vtx_diff_shwr_total.add_subplot(111)
mult_vtx = [vtx_diff_shwr, vtx_diff_elec]
_ = plt.hist(mult_vtx, 80, (0, 80), histtype='bar', fill=True, stacked=False, color=[
             'goldenrod', 'darkmagenta'], label=['All Showers', 'True Electron'])
ax.set_xlabel('Reco - True Shower Vertex Total Distance[cm]')
plt.legend()
fig_vtx_diff_shwr_total.savefig('reco-true_shower_vtx_total.pdf')
plt.close()

# histograms for direction
fig_dir_shwr_diff = plt.figure()
ax = fig_dir_shwr_diff.add_subplot(111)
mult_dir_shwr = [dir_diff_X_shwr, dir_diff_Y_shwr,
                 dir_diff_Z_shwr]
_ = plt.hist(mult_dir_shwr, 20, (-2, 2), histtype='bar', fill=True, stacked=True, color=[
             'tomato', 'goldenrod', 'darkmagenta'], label=['X', 'Y', 'Z'])
ax.set_xlabel('Reco - True Shower Direction')
plt.legend()
fig_dir_shwr_diff.savefig('shower_reco-true_direction.pdf')
plt.close()

fig_dir_shwr_dot = plt.figure()
ax = fig_dir_shwr_dot.add_subplot(111)
_ = plt.hist(dir_diff_shwr, 40, (0, 180), histtype='bar',
             fill=True, stacked=True, color='tomato', label='Dot Product')
ax.set_xlabel('Shower Reco:True Theta')
plt.legend()
fig_dir_shwr_dot.savefig('shower_reco_true_dot_dir.pdf')
plt.close()

fig_dir_elec_diff = plt.figure()
ax = fig_dir_elec_diff.add_subplot(111)
mult_dir_elec = [dir_diff_X_elec, dir_diff_Y_elec,
                 dir_diff_Z_elec]
_ = plt.hist(mult_dir_elec, 20, (-2, 2), histtype='bar', fill=True, stacked=True, color=[
             'tomato', 'goldenrod', 'darkmagenta'], label=['X', 'Y', 'Z'])
ax.set_xlabel('Reco - True Electron Direction')
plt.legend()
fig_dir_elec_diff.savefig('reco-true_electron_direction.pdf')
plt.close()

fig_dir_elec_dot = plt.figure()
ax = fig_dir_elec_dot.add_subplot(111)
_ = plt.hist(dir_diff_elec, 40, (0, 180), histtype='bar',
             fill=True, stacked=True, color='tomato', label='Dot Product')
ax.set_xlabel('Electron Reco:True Theta')
plt.legend()
fig_dir_elec_dot.savefig('electron_reco_true_dot_dir.pdf')
plt.close()

# histograms completeness
fig_energy_pfp = plt.figure()
ax = fig_energy_pfp.add_subplot(111)
if(cosmic_file == 'False'):
    mult_eng_pfp = [pion_completeness, neutron_completeness,
                    proton_completeness, photon_completeness, electron_completeness]
    _ = plt.hist(mult_eng_pfp, 40, (0, 1), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron'])
if(cosmic_file == 'True'):
    mult_eng_pfp = [pion_completeness, neutron_completeness,
                    proton_completeness, photon_completeness, electron_completeness, cosmic_completeness]
    _ = plt.hist(mult_eng_pfp, 40, (0, 1), histtype='bar', fill=True, stacked=True, color=[
        'wheat', 'goldenrod', 'darkmagenta', 'skyblue', 'tomato', 'darkslategray'], label=['Pion', 'Neutron', 'Proton', 'Photon', 'Electron', 'Cosmics'])
ax.set_xlabel('Completeness - Matched / MC Hits')
plt.legend(loc='upper left')
fig_energy_pfp.savefig('shower_completeness_type.pdf')
plt.close()

fig_completeness_mode = plt.figure()
ax = fig_completeness_mode.add_subplot(111)
mult_dir = [ccqe_completeness, res_completeness,
            dis_completeness, coh_completeness]
_ = plt.hist(mult_dir, 40, (0, 1), histtype='bar', fill=True, stacked=True, color=[
             'wheat', 'goldenrod', 'darkmagenta', 'skyblue'], label=['CCQE', 'Resonant', 'DIS', 'Coherent'])
ax.set_xlabel('Completeness - Matched / MC Hits')
plt.legend(loc='upper left')
fig_completeness_mode.savefig('shower_completeness_mode.pdf')
plt.close()

####################################
# how often are nues actually nues?
####################################
df_pfp_nues = getPfpNue(df_valid_pfpart)
nPfpNues = len(df_pfp_nues.index)
print 'Number of PFP Nues: ', nPfpNues
mcPdg_nue = df_pfp_nues.mcPdg
pfpPdg_nue = df_pfp_nues.pfoPdg
n_available_hits_nue = []
available_hits_nue = df_pfp_nues.nAvailableHits
mcVtxX = df_pfp_nues.mcVtxX
mcVtxY = df_pfp_nues.mcVtxY
mcVtxZ = df_pfp_nues.mcVtxZ
pfpVtxX = df_pfp_nues.pfoVtxX
pfpVtxY = df_pfp_nues.pfoVtxY
pfpVtxZ = df_pfp_nues.pfoVtxZ
vtx_diff_X = []
vtx_diff_Y = []
vtx_diff_Z = []
vtx_diff = []
nue_pfp_vtxX = []
nue_pfp_vtxY = []
nue_pfp_vtxZ = []
mcDirX = df_pfp_nues.mcDirX
mcDirY = df_pfp_nues.mcDirY
mcDirZ = df_pfp_nues.mcDirZ
pfpDirX = df_pfp_nues.pfoDirX
pfpDirY = df_pfp_nues.pfoDirY
pfpDirZ = df_pfp_nues.pfoDirZ
dir_diff_X = []
dir_diff_Y = []
dir_diff_Z = []
dir_diff = []
pfp_nues_counter = 0
all_pfp_nues_counter = 0
for nues in tqdm(mcPdg_nue.index):
    if(mcPdg_nue[nues] != pfpPdg_nue[nues]):
        tqdm.write('Nue: MC PDG and PFP PDG do not match!')
        continue
    n_available_hits_nue.append(available_hits_nue[nues])
    if(mcPdg_nue[nues] == 12):
        pfp_nues_counter = pfp_nues_counter + 1

    # let's do some true mc comparisons
    # vertex
    x_diff = pfpVtxX[nues] - mcVtxX[nues]
    y_diff = pfpVtxY[nues] - mcVtxY[nues]
    z_diff = pfpVtxZ[nues] - mcVtxZ[nues]
    total_diff = np.sqrt((x_diff * x_diff) +
                         (y_diff * y_diff) + (z_diff * z_diff))
    vtx_diff_X.append(x_diff)
    vtx_diff_Y.append(y_diff)
    vtx_diff_Z.append(z_diff)
    vtx_diff.append(total_diff)
    nue_pfp_vtxX.append(pfpVtxX[nues])
    nue_pfp_vtxY.append(pfpVtxY[nues])
    nue_pfp_vtxZ.append(pfpVtxZ[nues])
    # direction
    dir_x_diff = pfpDirX[nues] - mcDirX[nues]
    dir_y_diff = pfpDirY[nues] - mcDirY[nues]
    dir_z_diff = pfpDirZ[nues] - mcDirZ[nues]
    dir_dot = np.arccos(np.dot([pfpDirX[nues], pfpDirY[nues], pfpDirZ[nues]], [
        mcDirX[nues], mcDirY[nues], mcDirZ[nues]])) * (180 / 3.1415)
    # print 'PFP', pfpDirX[nues], pfpDirY[nues], pfpDirZ[nues]
    # print 'MC ', mcDirX[nues], mcDirY[nues], mcDirZ[nues]
    #! PFP Direction for neutrinos doesn't exist - it's zero! -- this needs fixing!?! !#
    dir_diff_X.append(dir_x_diff)
    dir_diff_Y.append(dir_y_diff)
    dir_diff_Z.append(dir_z_diff)
    dir_diff.append(dir_dot)

# histograms for vtx
fig_vtx_diff = plt.figure()
ax = fig_vtx_diff.add_subplot(111)
mult_vtx = [vtx_diff_X, vtx_diff_Y,
            vtx_diff_Z]
_ = plt.hist(mult_vtx, 80, (-80, 80), histtype='bar', fill=True, stacked=True, color=[
             'tomato', 'goldenrod', 'darkmagenta'], label=['X', 'Y', 'Z'])
ax.set_xlabel('Reco - True Nue Vertex [cm]')
plt.legend()
fig_vtx_diff.savefig('reco-true_nue_vtx.pdf')
plt.close()

fig_vtx_diff = plt.figure()
ax = fig_vtx_diff.add_subplot(111)
mult_vtx = [vtx_diff]
_ = plt.hist(mult_vtx, 80, (0, 160), histtype='bar', fill=True, stacked=True, color=[
             'skyblue'], label=['Total'])
ax.set_xlabel('Reco - True Nue Vertex [cm]')
plt.legend()
fig_vtx_diff.savefig('reco-true_nue_vtx_total.pdf')
plt.close()

################################################
# direction isn't being filled for pfp properly...
# histograms for direction
# fig_dir_diff = plt.figure()
# ax = fig_dir_diff.add_subplot(111)
# mult_dir = [dir_diff_X, dir_diff_Y,
#             dir_diff_Z]
# _ = plt.hist(mult_dir, 20, (-1, 1), histtype='bar', fill=True, stacked=True, color=[
#              'tomato', 'goldenrod', 'darkmagenta'], label=['X', 'Y', 'Z'])
# ax.set_xlabel('Reco - True Direction')
# plt.legend()
# plt.show()
#
# fig_dir_dot = plt.figure()
# ax = fig_dir_dot.add_subplot(111)
# _ = plt.hist(dir_diff, 40, (0, 180), histtype='bar',
#              fill=True, stacked=True, color='tomato', label='Dot Product')
# ax.set_xlabel('Cosine(Theta)')
# plt.legend()
# plt.show()

##########################################
# some comparisons between nue and showers
##########################################
diff_nue_shower_vtx_elec = []
diff_nue_shower_vtx_gamma = []
diff_nue_shower_vtx_prot = []
for nues in tqdm(df_pfp_nues.index):
    if(mcPdg_nue[nues] != pfpPdg_nue[nues]):
        tqdm.write('Nue: MC PDG and PFP PDG do not match!')
        continue

    temp_df1 = df_pfp_showers.drop(
        df_pfp_showers[df_pfp_showers.event != df_pfp_nues.event[nues]].index)
    for showers in (temp_df1.index):
        # print 'Nue Vtx X: ', pfpVtxX[nues]
        # print 'Particle: ', showers, ' has a Vtx X: ', pfpVtxX_shwr[showers]
        if(temp_df1.pfoPdg[showers] == 11):
            vtxX_shwr = pfpVtxX_shwr[showers]
            vtxY_shwr = pfpVtxY_shwr[showers]
            vtxZ_shwr = pfpVtxZ_shwr[showers]
            _diff_nue_shower_vtx_elec = np.sqrt(
                ((vtxX_shwr - pfpVtxX[nues]) * (vtxX_shwr - pfpVtxX[nues])) +
                ((vtxY_shwr - pfpVtxY[nues]) * (vtxY_shwr - pfpVtxY[nues])) +
                ((vtxZ_shwr - pfpVtxZ[nues]) * (vtxZ_shwr - pfpVtxZ[nues])))
            diff_nue_shower_vtx_elec.append(_diff_nue_shower_vtx_elec)

        if(temp_df1.pfoPdg[showers] == 22):
            vtxX_shwr = pfpVtxX_shwr[showers]
            vtxY_shwr = pfpVtxY_shwr[showers]
            vtxZ_shwr = pfpVtxZ_shwr[showers]
            _diff_nue_shower_vtx_gamma = np.sqrt(
                ((vtxX_shwr - pfpVtxX[nues]) * (vtxX_shwr - pfpVtxX[nues]))
                + ((vtxY_shwr - pfpVtxY[nues]) * (vtxY_shwr - pfpVtxY[nues]))
                + ((vtxZ_shwr - pfpVtxZ[nues]) * (vtxZ_shwr - pfpVtxZ[nues])))
            diff_nue_shower_vtx_gamma.append(_diff_nue_shower_vtx_gamma)

        if(temp_df1.pfoPdg[showers] == 2212):
            vtxX_shwr = pfpVtxX_shwr[showers]
            vtxY_shwr = pfpVtxY_shwr[showers]
            vtxZ_shwr = pfpVtxZ_shwr[showers]
            _diff_nue_shower_vtx_prot = np.sqrt(
                ((vtxX_shwr - pfpVtxX[nues]) * (vtxX_shwr - pfpVtxX[nues]))
                + ((vtxY_shwr - pfpVtxY[nues]) * (vtxY_shwr - pfpVtxY[nues]))
                + ((vtxZ_shwr - pfpVtxZ[nues]) * (vtxZ_shwr - pfpVtxZ[nues])))
            diff_nue_shower_vtx_prot.append(_diff_nue_shower_vtx_prot)

# histogram of nue vtx to shower vtx
fig_vtx_diff = plt.figure()
ax = fig_vtx_diff.add_subplot(111)
mult_vtx = [diff_nue_shower_vtx_elec, diff_nue_shower_vtx_gamma,
            diff_nue_shower_vtx_prot]
_ = plt.hist(mult_vtx, 80, (0, 80), histtype='bar', fill=True, stacked=False, color=[
    'tomato', 'skyblue', 'darkmagenta'], label=['Electron', 'Photon', 'Proton'])
ax.set_xlabel('Reco: Shower - Nue Vertex [cm]')
plt.legend()
fig_vtx_diff.savefig('distance_nue_shower_vtx.pdf')
plt.close()

####################
# Print statements
####################
print 'Nues: '
nue_efficiency = float(nPfpNues) / float(nMCNues) * 100.
print 'All  Nue Reco Efficiency: ', nue_efficiency, '(', float(nPfpNues), ' / ', float(nMCNues), ')'

true_nue_efficiency = float(pfp_nues_counter) / float(nMCNues) * 100.
print 'True Nue Reco Efficiency: ', true_nue_efficiency, '(', float(pfp_nues_counter), ' / ', float(nMCNues), ')'

nue_purity = float(pfp_nues_counter) / float(nPfpNues) * 100.
print 'True Nue Reco Purity', nue_purity, '(', float(pfp_nues_counter), ' / ', float(nPfpNues), ')'

# for showers too
print 'Showers: '
shower_efficiency = float(true_pfp_counter) / float(mc_electron_counter) * 100.
print 'True Electron Reco Efficiency: ', shower_efficiency, '(', float(true_pfp_counter), ' / ', float(mc_electron_counter), ')'

shower_purity = float(true_pfp_counter) / float(nPfpShowers) * 100.
print 'Electron Reco Purity', shower_purity, '(', float(true_pfp_counter), ' / ', float(nPfpShowers), ')'

elapsed = timeit.default_timer() - start_time
print 'Time Elapsed: ', elapsed
