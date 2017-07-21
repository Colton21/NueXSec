# Author: Colton Hill
# This file should create functions used to compare the reco-true information

# let's only have events in our time window
from tqdm import tqdm
from shapely import geometry
from shapely.geometry import Polygon
from shapely.geometry import Point
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
            # tqdm.write('Data Frame Empty - maybe something is wrong!')
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

# construct a plane


def constructPlane(point1, point2, point3):
    # determine the vectors
    v_ab_i = (point2[0] - point1[0])
    v_ab_j = (point2[1] - point1[1])
    v_ab_k = (point2[2] - point1[2])
    v_ab = [v_ab_i, v_ab_j, v_ab_k]

    v_ac_i = (point3[0] - point1[0])
    v_ac_j = (point3[1] - point1[1])
    v_ac_k = (point3[2] - point1[2])
    v_ac = [v_ac_i, v_ac_j, v_ac_k]

    x_ab_ac = np.cross(v_ab, v_ac)
    return x_ab_ac

# distance to nearest wall calculation


def distToWall(start_x, start_y, start_z):
    tpc_x1 = 0
    tpc_x2 = 256.35
    tpc_y1 = -116.5
    tpc_y2 = 116.5
    tpc_z1 = 0
    tpc_z2 = 1036.8

    tpc_x_mid = (float(tpc_x2) - float(tpc_x1)) / 2.
    tpc_y_mid = (float(tpc_y2) - float(abs(tpc_y1))) / 2.
    tpc_z_mid = (float(tpc_z2) - float(tpc_z1)) / 2.

    nearest_x = 0
    nearest_y = 0
    nearest_z = 0

    point1_yz = []
    point2_yz = []
    point3_yz = []
    point1_xz = []
    point2_xz = []
    point3_xz = []
    point1_xy = []
    point2_xy = []
    point3_xy = []

    using_x = []
    using_y = []
    using_z = []

    if(start_x <= tpc_x_mid):
        # yz plane at x1
        point1_yz = (tpc_x1, tpc_y1, tpc_z1)
        point2_yz = (tpc_x1, tpc_y2, tpc_z1)
        point3_yz = (tpc_x1, tpc_y1, tpc_z2)
        using_x = tpc_x1
    if(start_x > tpc_x_mid):
        # yz plane at x2
        point1_yz = (tpc_x2, tpc_y1, tpc_z1)
        point2_yz = (tpc_x2, tpc_y2, tpc_z1)
        point3_yz = (tpc_x2, tpc_y1, tpc_z2)
        using_x = tpc_x2

    if(start_y <= tpc_y_mid):
        # xz plane at y1
        point1_xz = (tpc_x1, tpc_y1, tpc_z1)
        point2_xz = (tpc_x2, tpc_y1, tpc_z1)
        point3_xz = (tpc_x1, tpc_y1, tpc_z2)
        using_y = tpc_y1
    if(start_y > tpc_y_mid):
        # xz plane at y2
        point1_xz = (tpc_x1, tpc_y2, tpc_z1)
        point2_xz = (tpc_x2, tpc_y2, tpc_z1)
        point3_xz = (tpc_x1, tpc_y2, tpc_z2)
        using_y = tpc_y2

    if(start_z <= tpc_z_mid):
        # xy plane at z1
        point1_xy = (tpc_x1, tpc_y1, tpc_z1)
        point2_xy = (tpc_x2, tpc_y1, tpc_z1)
        point3_xy = (tpc_x1, tpc_y2, tpc_z1)
        using_z = tpc_z1
    if(start_z > tpc_z_mid):
        # xy plane at z2
        point1_xy = (tpc_x1, tpc_y1, tpc_z2)
        point2_xy = (tpc_x2, tpc_y1, tpc_z2)
        point3_xy = (tpc_x1, tpc_y2, tpc_z2)
        using_z = tpc_z2

    # now I need to construct 3 planes and check the minimum distance
    plane_xy = constructPlane(point1_xy, point2_xy, point3_xy)
    # ax+by+cz+d=0
    d_xy = -1 * ((plane_xy[0] * using_x) + (plane_xy[1]
                                            * using_y) + (plane_xy[2] * using_z))
    mod_xy = np.sqrt((plane_xy[0] * plane_xy[0]) + (plane_xy[1] *
                                                    plane_xy[1]) + (plane_xy[2] * plane_xy[2]))
    plane_xz = constructPlane(point1_xz, point2_xz, point3_xz)
    d_xz = -1 * ((plane_xz[0] * using_x) + (plane_xz[1]
                                            * using_y) + (plane_xz[2] * using_z))
    mod_xz = np.sqrt((plane_xz[0] * plane_xz[0]) + (plane_xz[1] *
                                                    plane_xz[1]) + (plane_xz[2] * plane_xz[2]))
    plane_yz = constructPlane(point1_yz, point2_yz, point3_yz)
    d_yz = -1 * ((plane_yz[0] * using_x) + (plane_yz[1]
                                            * using_y) + (plane_yz[2] * using_z))
    mod_yz = np.sqrt((plane_yz[0] * plane_yz[0]) + (plane_yz[1] *
                                                    plane_yz[1]) + (plane_yz[2] * plane_yz[2]))

    # calculate the 3D distance to the nearest wall
    start_point = (start_x, start_y, start_z)
    dist_xy = abs(((np.dot(start_point, plane_xy) + d_xy) / mod_xy))
    dist_xz = abs(((np.dot(start_point, plane_xz) + d_xz) / mod_xz))
    dist_yz = abs(((np.dot(start_point, plane_yz) + d_yz) / mod_yz))

    # find the minimum
    # print 'XY: ', dist_xy
    # print 'XZ: ', dist_xz
    # print 'YZ: ', dist_yz
    if(dist_xy <= dist_xz):
        if(dist_xy <= dist_yz):
            return dist_xy
    if(dist_xz <= dist_xy):
        if(dist_xz <= dist_yz):
            return dist_xz
    if(dist_yz <= dist_xy):
        if(dist_yz <= dist_xz):
            return dist_yz

    print 'Ooops you are an idiot!'
    return 0


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
    # dataframe = dataframe.drop(dataframe[dataframe.pfoEnergy == 0.0].index)
    return dataframe


# get the mc neutrino objects
def getMCNeutrino(dataframe):
    dataframe = dataframe.drop(
        dataframe[(dataframe.mcIsNeutrino == False)].index)
    return dataframe

# get the nue cc objects


def getMCCCNue(dataframe):
    dataframe = dataframe.drop(
        dataframe[(dataframe.mcIsCC != True) & (dataframe.mcPdg != 12)].index)
    return dataframe

# use this to breakdown mc Pdg of remaining reco nus


def mcPartBreakdown(dataframe):
    if(dataframe.empty):
        print 'Dataframe is Empty!'
        exit(1)
    num_all = len(dataframe.index)
    info_list = []

    df_nue_cc = dataframe.drop(
        dataframe[(dataframe.mcIsCC == False) | (dataframe.mcPdg != 12)].index)
    num_nue_cc = len(df_nue_cc.index)
    df_nue_nc = dataframe.drop(
        dataframe[(dataframe.mcIsCC == True) | (dataframe.mcPdg != 12)].index)
    num_nue_nc = len(df_nue_nc.index)
    df_numu = dataframe.drop(
        dataframe[(dataframe.mcPdg != 14)].index)
    num_numu = len(df_numu.index)

    # calculate the purity of the sample
    #purity = float(num_nue_cc) / float(num_all) * 100.

    info_list.append(num_all)
    info_list.append(num_nue_cc)
    info_list.append(num_nue_nc)
    info_list.append(num_numu)
    # info_list.append(purity)

    print 'Number of Nue CC: ', num_nue_cc
    print 'Number of Nue NC: ', num_nue_nc
    print 'Number of Numu  : ', num_numu
    # print 'Purity          : ', purity

    return info_list

# cosmic breakdown - most conservative


def cosmicBreakdown(dataframe):
    df_cosmic = dataframe.drop(
        dataframe[(dataframe.mcNuPdg != 0)].index)
    df_cosmic = df_cosmic.drop(
        df_cosmic[(df_cosmic.pfoNuPdg != 12)].index)
    num_cosmic = len(df_cosmic.index)

    print 'Number of Cosmic: ', num_cosmic

    return num_cosmic


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
        dataframe[dataframe.pfoNuPdg != 12].index)
    return dataframe

# get the pfp cosmic nue objects


def getCosmicPfpNue(dataframe):
    dataframe = dataframe.drop(
        dataframe[(dataframe.mcNuPdg != 0)].index)
    dataframe = dataframe.drop(
        dataframe[(dataframe.pfoNuPdg != 12)].index)
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
