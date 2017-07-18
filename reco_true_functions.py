# Author: Colton Hill
# This file should create functions used to compare the reco-true information

# let's only have events in our time window
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import sys


def inTimeList(dataframe, startTime, endTime, flashThreshold):
    failEvents = []
    # more than one entry in some events -> need to drop duplicates
    temp_df_evt = dataframe.drop_duplicates('event')
    for i in tqdm(temp_df_evt.event):
        eventPass = False
        temp_df1 = dataframe.drop(
            dataframe[dataframe.event != i].index)
        for flash in temp_df1.index:
            if(temp_df1.OpFlashPE[flash] >= flashThreshold):
                if(temp_df1.OpFlashTime[flash] <= endTime and temp_df1.OpFlashTime[flash] >= startTime):
                    # then this whole event passes!
                    eventPass = True
        if(eventPass == False):
            failEvents.append(i)

    return failEvents


def inTimeDataFrame(dataframe, failEvents):
    for events in failEvents:
        dataframe = dataframe.drop(
            dataframe[dataframe.event == events].index)
    return dataframe


# let's work on the flash to reco vertex distance
def flashRecoVtxDist(dataframe_reco, dataframe_opt, distance):

    failEvents = []
    yz_dist_list = []
    yz_dist_list_cosmic = []
    shwr_to_vtx_dist = distance
    # make sure we're only looking at pfp showers
    dataframe_showers = dataframe_reco.drop(
        dataframe_reco[(dataframe_reco.pfoPdg != 11)].index)

    # I want to compare the vertex of each shower to that of the largest flash
    # vertex in YZ PER EVENT
    # need to first make a unique list of event - some duplicates, as more
    # than one event in an event!
    temp_df_evt = dataframe_showers.drop_duplicates('event')
    for i in tqdm(temp_df_evt.event):
        # construct dataframes which only have information for these events
        temp_df1 = dataframe_showers.drop(
            dataframe_showers[dataframe_showers.event != i].index)
        temp_df2 = dataframe_opt.drop(
            dataframe_opt[dataframe_opt.event != i].index)
        if(temp_df2.empty == True):
            #tqdm.write('Data Frame Empty - maybe something is wrong!')
            failEvents.append(i)
            continue
        # find flash with largest PE
        lrgFlash_index = temp_df2.OpFlashPE.argmax()
        event_IsBad = True
        opt_vtxY = temp_df2.get_value(lrgFlash_index, 'OpFlashCenterY')
        opt_vtxZ = temp_df2.get_value(lrgFlash_index, 'OpFlashCenterZ')
        for shower in temp_df1.index:
            truePdg = dataframe_showers.mcNuPdg[shower]
            shower_vtxY = dataframe_showers.pfoVtxY[shower]
            shower_vtxZ = dataframe_showers.pfoVtxZ[shower]
            yz_dist = np.sqrt(((opt_vtxY - shower_vtxY) * (opt_vtxY - shower_vtxY)) +
                              ((opt_vtxZ - shower_vtxZ) * (opt_vtxZ - shower_vtxZ)))
            if(truePdg == 0):
                yz_dist_list_cosmic.append(yz_dist)
            if(truePdg != 0):
                yz_dist_list.append(yz_dist)
            if(yz_dist <= shwr_to_vtx_dist):
                event_IsBad = False
            if(event_IsBad == False):
                break

        if(event_IsBad == True):
            failEvents.append(i)

    for events in failEvents:
        dataframe_reco = dataframe_reco.drop(
            dataframe_reco[dataframe_reco.event == events].index)

    # let's just do the plotting here for now...
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mult = [yz_dist_list, yz_dist_list_cosmic]
    _ = plt.hist(mult, 50, (0, 300), histtype='bar',
                 fill=True, stacked=True, color=['tomato', 'darkslategray'], label=['Nue', 'Cosmics'])
    ax.set_xlabel('Dist. Reco Shower Vertex to Largest Flash YZ [cm]')
    plt.legend()
    plt.show()

    return dataframe_reco

# let's also make a list of dropped events


def outTimeWindow(dataframe, startTime, endTime):
    dataframe = dataframe.drop(dataframe[(dataframe.OpFlashTime <= endTime) & (
        dataframe.OpFlashTime >= startTime)].index)
    return dataframe

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

# let's trim the sample down to only events in the TPC


def pfpInFV(dataframe, x_min, x_max, y_min, y_max, z_min, z_max):
    tpc_x1 = 0
    tpc_x2 = 256.35
    tpc_y1 = -116.5
    tpc_y2 = 116.5
    tpc_z1 = 0
    tpc_z2 = 1036.8
    dataframe = dataframe.drop(
        dataframe[(dataframe.pfoVtxX < tpc_x1 + x_min) | (dataframe.pfoVtxX > tpc_x2 - x_max)].index)
    dataframe = dataframe.drop(
        dataframe[(dataframe.pfoVtxY < tpc_y1 + y_min) | (dataframe.pfoVtxY > tpc_y2 - y_max)].index)
    dataframe = dataframe.drop(
        dataframe[(dataframe.pfoVtxZ < tpc_z1 + z_min) | (dataframe.pfoVtxZ > tpc_z2 - z_max)].index)
    return dataframe

# some events have no MC particles - these are uninteresting!


def removeZeroMCP(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.nMCParticles == 0].index)
    return dataframe

# some mc particles have zero MC energy - we should remove these


def removeZeroMCEng(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.mcEnergy == 0.0].index)
    return dataframe

# grab the 0 index particles - index is overloaded so have to use iloc to get
# the fourth column


def zeroIndex(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.iloc[:, 3] != 0].index)
    return dataframe

# construct dataframe of only valid PFParticles]


def validPFParticle(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.pfoPdg == 0].index)
    dataframe = dataframe.drop(dataframe[
        (dataframe.pfoVtxX == 0.0) &
        (dataframe.pfoVtxY == 0.0) &
        (dataframe.pfoVtxZ == 0.0)].index)
    #dataframe = dataframe.drop(dataframe[dataframe.pfoEnergy == 0.0].index)
    return dataframe


# get the mc neutrino objects
def getMCNeutrino(dataframe):
    dataframe = dataframe.drop(
        dataframe[(dataframe.mcIsNeutrino == False)].index)
    return dataframe

# get the primary mc particles


def getMCPrimary(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.mcIsPrimary == False].index)
    return dataframe

# get the pfp neutrino objects


def getPfpNeutrino(dataframe):
    dataframe = dataframe.drop(
        dataframe[(dataframe.pfoPdg != 12) & (dataframe.pfoPdg != 14)].index)
    return dataframe

# get the pfp nue objects


def getPfpNue(dataframe):
    dataframe = dataframe.drop(
        dataframe[(dataframe.pfoPdg != 12)].index)
    return dataframe

# get all pfp objects


def getAllPfps(dataframe):
    dataframe = dataframe.drop(dataframe[(dataframe.pfoPdg == 0)].index)
    return dataframe

# get all pfp showers


def getPfpShower(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.pfoPdg != 11].index)
    return dataframe

# get only CC MC events


def getCConly(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.mcIsCC == False].index)
    return dataframe
