# Author: Colton Hill
# This file should create functions used to compare the reco-true information

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

# grab the 0 index particles - index is overloaded so have to use iloc to get
# the third column


def zeroIndex(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.iloc[:, 2] != 0].index)
    return dataframe

# construct dataframe of only valid PFParticles]


def validPFParticle(dataframe):
    dataframe = dataframe.drop(dataframe[dataframe.pfoPdg == 0].index)
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
