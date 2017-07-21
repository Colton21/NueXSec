import math
from reco_true_functions import *

# Calculating the cross section requires a few things
# xsec = N_total - N_bkg / (flux * ar_nucleons * efficiency)

# 5822 nu events for 3.8322x10^19 POT w/ 225 nu_e, 167 nu_e CC
# NuMIFlux.root results for 6.6x10^20 POT/cm^2
# ar_nucleons = 4.76*pow(10,31); //This is #nucleons - /40. to get # Ar
# molecules
# total_measured = 167./(3.86375*pow(10,19)); /// total # nu_e CC
# scale factor: 2409 nue cc per 6.0e20 POT (1 year of running)


def calcXSec(n_total, n_bkg, flux, efficiency):
    n_events = n_total - n_bkg
    # scale flux and events
    scale_factor = 2.4 * math.pow(10, 17)
    # number of nucleons
    num_nucleons = 4.76 * math.pow(10, 31)

    xsec_cc = (n_events * scale_factor) / (flux * num_nucleons * efficiency)

    return xsec_cc


def printInfo(df, num_mc_cc_nue):
    ##print statements##
    df_nu = getPfpNeutrino(df)
    info_list = mcPartBreakdown(df_nu)
    num_nue_cosmic = cosmicBreakdown(df)
    num_pfp_cc_nue = info_list[1]
    print 'Efficiency: ', float(num_pfp_cc_nue) / float(num_mc_cc_nue) * 100.
    purity = float(num_pfp_cc_nue) / \
        (float(info_list[0]) + float(num_nue_cosmic)) * 100.
    print 'Purity: ', purity
