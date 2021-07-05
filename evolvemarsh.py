"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _5 July 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import math
import matplotlib.pyplot as plt


def evolvemarsh(
        marshelevation,
        msl,
        C_e,
        OCb,
        tr,
        numiterations,
        P,
        dt,
        ws,
        timestep,
        BMax,
        Dmin,
        Dmax,
        rhoo,
        rhos,
        plot
):
    """Calculates biomass and mineral and organic deposition for each cell in the marsh as a function of flooding
    frequency; calculates the total flux of sediment onto the marsh from the bay."""

    L = len(marshelevation)
    time_submerged = np.zeros([numiterations, L])
    sedimentcycle = np.zeros([numiterations, L])
    depth = np.zeros([numiterations, L])
    C = np.zeros([L])

    # Loop through a tidal cycle to determine flooding duration for each point in the marsh
    for i in range(1, numiterations):
        depth_iteration = 0.5 * tr * math.sin(2 * math.pi * ((i + 1) * dt / P)) + (
                msl - marshelevation[0: L + 1])  # [m] Depth at each position in the marsh
        depth[i, :] = depth_iteration  # Store depths in array
        temp = np.zeros([L])
        temp[depth_iteration > 0] = dt
        time_submerged[i, :] = temp  # Inundation for a single flood cycle recorded for each cell

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

    # -------------------------
    # Mineral Deposition

    distance = 0  # [m] Initialize, distance from marsh edge
    for xx in range(L):
        if bgb[xx] > 0:
            distance += 1  # [m]
            C[xx] = C_e * math.exp(-0.0031 * distance)  # [kg/m3] coefficient of -0.0031 is a fitted parameter for realistic marsh
            # topography
        else:
            distance = 1  # [m]
            C_e = C_e * 0.9  # [kg/m3] Decrease concentration at the new "marsh edge" by 10% with each subsequent pond formation
            C[xx] = C_e * math.exp(-0.0031 * distance)  # [kg/m3] coefficient of -0.0031 is a fitted parameter for realistic marsh
            # topography

    floodfraction = np.sum(time_submerged, axis=0) / P  # Portion of the tidal cycle that each point is submerged

    for i in range(1, numiterations):
        tempdepth = depth[i, :]
        tempy = np.zeros([L])
        tempy[tempdepth > 0] = C[tempdepth > 0] * ws * dt  # [kg] mass of mineral sediment deposited, where depth > 0
        sedimentcycle[i, :] = tempy

    susp_dep = np.sum(sedimentcycle, axis=0) * timestep * 1000  # [g/yr] suspended sediment deposition in an entire year, in each cell

    if plot:
        plt.figure()
        plt.subplot(3, 1, 1)
        plt.plot(np.transpose(sedimentcycle))
        plt.xlabel("Distance Across Marsh [m]")
        plt.ylabel("Mineral Sediment Deposited Per Cycle [g/cycle]")
        plt.subplot(3, 1, 2)
        plt.plot(C)
        plt.xlabel("Distance Across Marsh [m]")
        plt.ylabel("Suspended Sediment Concentration [kg/m3]")
        plt.subplot(3, 1, 3)
        plt.plot(susp_dep)
        plt.xlabel("Distance Across Marsh [m]")
        plt.ylabel("Suspended Sediment Deposition [g/yr]")
        plt.show()

    mineral = susp_dep * (1 - OCb)  # [g] Mineral deposition of suspended sediment in a given year is equal determined by the organic content of the bay sediment
    organic_alloch = susp_dep * OCb  # [g] Organic deposition of suspended sediment is equal determined by the organic content of bay sed

    # -------------------------
    # Calculations & Conversions

    Fm_min = np.sum(mineral) / 1000  # [kg/yr] Flux of mineral sediment from the bay
    Fm_org = np.sum(organic_alloch) / 1000  # [kg/yr] Flux of organic sediment from the bay

    # Calculate thickness of new sediment (mineral+organic) based off LOI and its effect on density
    loi = (organic_autoch + organic_alloch) / (mineral + organic_autoch + organic_alloch)
    density = 1 / ((loi / rhoo) + ((1 - loi) / rhos)) * 1000  # [g/m3] Bulk density is calculated according to Morris et al. (2016)
    density[np.isnan(density)] = 1  # If there is zero deposition, loi calculation will divide by zero and make density nan. Set density
    # equal to one in this case, so that accretion is zero, instead of nan.

    accretion = (mineral + organic_autoch + organic_alloch) / density  # [m] accretion in a given year

    # Update elevation
    marshelevation += accretion

    return marshelevation, organic_autoch, organic_alloch, mineral, Fm_min, Fm_org, bgb, accretion, agb
