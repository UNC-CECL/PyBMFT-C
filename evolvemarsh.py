"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _10 June 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import math


def evolvemarsh(
        marshelevation,
        msl,
        C_e,
        OCb,
        tr,
        numiterations,  # Check naming
        P,
        dt,  # Check naming
        ws,
        timestep,  # Check naming
        BMax,
        Dmin,  # Add this back to transect.py
        Dmax,  # Add this back to transect.py
        rhoo,
        rhos,
):
    """Calculates biomass and mineral and organic deposition for each cell in the marsh as a function of flooding
    frequency; calculates the total flux of sediment onto the marsh from the bay."""

    L = len(marshelevation)
    time_submerged = np.zeros([numiterations, L])
    sedimentcycle = np.zeros([L])
    vegtype = np.zeros([L])
    bgb = np.zeros([L])
    depth = np.zeros([numiterations, L])
    C = np.zeros([L])
    d = 0


    # Loop through a tidal cycle to determine flooding duration for each point in the marsh
    for i in range(1, numiterations):
        depth_iteration = 0.5 * tr * math.sin(2 * math.pi * (i * dt / P)) + (
                    msl - marshelevation)  # [m] Depth at each position in the marsh
        depth[i, :] = depth_iteration  # Store depths in array
        time_submerged[i, :][depth_iteration > 0] = dt  # Inundation for a single flood cycle recorded for each cell

    # -------------------------
    # Belowground Productivity

    # Creates a biomass curve (Mariotti & Carr, 2014) where peak biomass occurs at a depth halfway between the maximum
    # depth for vegetation and the minimum (here, mean high water level)

    dm = msl + tr / 2 - marshelevation  # [m] Depth of the marsh surface below HWL at any given point

    bgb = np.zeros([L])  # Belowground biomass
    agb = np.zeros([L])  # Aboveground biomass
    organic_autoch = np.zeros([L])  # Autochthonous organic material

    for ii in range(L):
        if dm[ii] > Dmax:  # If depth is below vegetation maximum, there is no production
            agb[ii] = 0  # [g] Mudflat
            bgb[ii] = 0  # [g] Mudflat
            organic_autoch[ii] = 0  # [g] No autochthonous material stored in the soil; we do not multiply by lingin content (as in earlier
            # version of the model) because we subtract mass due to decomposition in the 'decompose' function
        elif dm[ii] <= Dmin:  # If depth is above vegetation minimum, there is very little belowground productivity
            if ii > 6000:
                agb[ii] = 100  # [g] Forest aboveground - constant for now, should depend on elevation
                bgb[ii] = 0.00001  # [g] Forest
                organic_autoch[ii] = bgb[ii]  # [g] Forest organic matter stored in soil
            else:
                agb[ii] = 1 * (BMax * (dm[ii] - Dmax) * (dm[ii] - Dmin) / (0.25 * (-Dmin - Dmax) * (Dmax - 3 * Dmin)))
                bgb[ii] = 0.1  # [g] "Forest" - really high marsh
                organic_autoch[ii] = bgb[ii]  # [g] Forest organic matter stored in soil
        else:
            agb[ii] = 1 * (BMax * (dm[ii] - Dmax) * (dm[ii] - Dmin) / (0.25 * (-Dmin - Dmax) * (Dmax - 3 * Dmin)))  # [g] Marsh
            bgb[ii] = BMax * (dm[ii] - Dmax) * (dm[ii] - Dmin) / (0.25 * (-Dmin - Dmax) * (Dmax - 3 * Dmin))  # [g] Marsh
            organic_autoch[ii] = bgb[ii]  # [g] As mentioned above, no longer multiply by lingin content

    # Marsh width only includes cells with belowground biomass
    #marshwidth = ... doesnt appear this needs to be calculated

    # -------------------------
    # Mineral Deposition

    distance = 0  # [m] Initialize, distance from marsh edge?

    for xx in range(L):
        if bgb[xx] > 0:
            distance += 1  # [m]
            C[xx] = C_e * math.exp(-0.0031 * distance)  # [kg/m3] coefficient of -0.0031 is a fitted parameter for realistic marsh topography
        else:
            distance = 1  # [m]
            C_e = C_e * 0.9  # [kg/m3] Decrease concentration at the new "marsh edge" by 10% with each subsequent pond formation
            C[xx] = C_e * math.exp(-0.0031 * distance)  # [kg/m3] coefficient of -0.0031 is a fitted parameter for realistic marsh topography

    floodfraction = np.sum(time_submerged, axis=1) / P  # Portion of the tidal cycle that each point is submerged

    print(floodfraction)  # IR 10Jun21: Paused here; running but not working properly

    organic_alloch = 0
    mineral = 0
    Fm_min = 0
    Fm_org = 0
    accretion = 0


    return marshelevation, organic_autoch, organic_alloch, mineral, Fm_min, Fm_org, bgb, accretion, agb
