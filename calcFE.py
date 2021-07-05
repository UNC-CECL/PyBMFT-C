"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _5 July 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import math
import bisect


def calcFE(bfoc, bfop, elevation, yr, organic_dep_autoch, organic_dep_alloch, mineral_dep, rhos):
    """Function to calculate the flux of organic matter (FE_org) and the flux of mineral sediment (FE_min) from the marsh to the bay,using the fetch for the current year (bfoc)
    the fetch for the previous year (bfop) and the stragtigraphy of organic and mineral deposition """

    # Calculate OM eroded from the marsh platform
    organic_dep = organic_dep_autoch + organic_dep_alloch
    pyr = yr - 1  # Previous year
    E = bfoc - bfop  # Amount of erosion b/t the previous yr and the current yr

    x_m1 = math.ceil(bfop) - 1  # First marsh cell to erode

    x_m2 = math.ceil(bfoc) - 1  # Last  marsh cell to erode

    if E <= 0:
        FE_org = 0  # [g] If there is no erosion, no OM is eroded
        FE_min = 0  # [g] If there is no erosion, no MM is eroded
    elif x_m1 == x_m2:  # If actively eroding marsh edge has not changed since previous year
        if elevation[pyr, 0] < elevation[0, x_m1]:  # If depth of erosion is below the lowest marsh deposit
            us = elevation[0, x_m1] - elevation[yr-1, 0]  # [m] Depth of underlying stratigraphy
            usmass = us * rhos * 1000  # [g] Mass of pure mineral sediment underlying marsh at marsh edge
            FE_org = np.sum(organic_dep[0: pyr + 1, x_m1]) * E  # [g] OM eroded is equal to the total amount of OM in the eroding marsh edge (including both initial deposit and OM
            # deposited since the model run began) times the fraction of the marsh edge cell that is eroded
            FE_min = np.sum(mineral_dep[0: pyr + 1, x_m1]) * E + usmass  # [g] MIN eroded is equal to the total amount of OM in the eroding marsh edge (including both initial
            # deposit and OM deposited since the model run began) times the fraction of the marsh edge cell that is eroded
        else:  # If depth of erosion is less than marsh deposits
            boundyr = bisect.bisect_right(elevation[:, x_m1], elevation[0, pyr])  # Year at which deposit is above the boundary for the depth of erosion
            FE_org = np.sum(organic_dep[boundyr: pyr + 1, x_m1]) * E  # [g] Total mass of OM deposited in the marsh edge cell
            FE_min = np.sum(mineral_dep[boundyr: pyr + 1, x_m1]) * E  # [g] Total mass of MIN deposited in the marsh edge cell
    else:
        ecells = x_m2 - x_m1  # Number of cells eroded
        Hfrac_ero = np.ones(ecells)
        Hfrac_ero[0] = x_m1 - (bfop - 1)  # Horizontal fraction of previous marsh edge that is eroded
        Hfrac_ero[-1] = bfoc - math.floor(bfoc)  # Horizontal fraction of current marsh edge that is eroded

        FE_org = 0
        FE_min = 0
        i = 0
        for x_m in range(x_m1, x_m2 + 1):
            if elevation[pyr, 0] < elevation[0, x_m]:  # If depth of erosion is below the lowest marsh deposit
                us = elevation[0, x_m] - elevation[yr - 1, 0]  # [m] Depth of underlying stratigraphy
                usmass = us * rhos * 1000  # [g] Mass of pure mineral sediment underlying marsh at marsh edge
                FE_org = FE_org + np.sum(organic_dep[0: pyr + 1, x_m]) * Hfrac_ero[i]  # [g] OM eroded from previous marsh edge cell
                FE_min = FE_min + np.sum(mineral_dep[0: pyr + 1, x_m]) * Hfrac_ero[i] + usmass  # [g] MIN eroded from previous marsh edge cell
            else:  # If depth of erosion is less than marsh deposit
                boundyr = bisect.bisect_right(elevation[:, x_m], elevation[0, pyr])  # Year at which deposit is above the boundary for the depth of erosion
                FE_org = np.sum(organic_dep[boundyr: pyr + 1, x_m]) * Hfrac_ero[i]  # [g] Total mass of OM deposited in the marsh edge cell
                FE_min = np.sum(mineral_dep[boundyr: pyr + 1, x_m]) * Hfrac_ero[i]  # [g] Total mass of MIN deposited in the marsh edge cell

    return FE_org, FE_min





