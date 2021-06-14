"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _9 June 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import os
import scipy.io
import math
import matplotlib.pyplot as plt


def buildtransect(R, C, slope, mwo, elev_25, amp, wind, bfo, endyear, plot):
    """Creates model domain and initial morphology of the bay, marsh, and upland slope. Marsh and bay depth are set
    to values close to equilibrium for the given sea level rise rate and suspended sediment concentration."""

    directory = "Calibrate/Fetch" + str(bfo) + "_Wind" + str(wind)  # Bay depth lookup table as function of fetch & wind

    spindur = np.size(elev_25, axis=0)

    # Determine initial bay depth, such that change in depth will be small
    if os.path.isdir(directory) is False:
        dfo = 2
        print("Warning: Initial conditions have not been calibrated for these parameters")
        # raise ValueError("Initial conditions have not been calibrated for these parameters")
    else:
        db_eq_dict = scipy.io.loadmat(directory + "/Equilibrium Bay Depth.mat")
        db_eq = db_eq_dict["db_eq"]
        if C / 10 > np.size(db_eq, axis=0) or R > np.size(db_eq, axis=1):
            dfo = 2
            print("Warning: Initial conditions have not been calibrated for these parameters")
            # raise ValueError("Initial conditions have not been calibrated for these parameters")
        elif C / 10 < 0.5 or R < 0.5:
            dfo = 2
            print("Warning: Initial conditions have not been calibrated for these parameters")
            # raise ValueError("Initial conditions have not been calibrated for these parameters")
        else:
            dfo = db_eq[(round(C / 10)) - 1, round(R) - 1]

    x_m = bfo + 1  # First marsh cell?

    maxY = R / 1000 * (endyear + 1) + amp + 0.5  # Maximum sea-level excursion
    max_potentialwidth = math.ceil(maxY / slope)
    upland_width = math.ceil(maxY / slope)

    if max_potentialwidth > upland_width:
        raise ValueError("Slope/sea-level rise conditions are such that the model domain is too small for the "
                         "scenario. Adjust the upland width accordingly.")

    B = bfo + mwo + upland_width  # [m] Total domain width, and also number of cells in domain each with 1 m width
    x = np.linspace(0, B - 1, num=B)  # x-position of each cell in model domain
    elevation = np.zeros([endyear + 1, B])
    elevation[:spindur, :x_m - 1] = amp - dfo  # Bay depth for first 25 (?) years is at equilibrium
    elevation[:spindur, x_m: x_m + mwo] = elev_25 - (spindur * (1 / 1000))  # Marsh elevation comes from model spinup, adjusted to modern sea level

    modernslope = slope * (np.linspace(1, upland_width, num=upland_width)) + elevation[spindur - 1, x_m + mwo - 1]

    # Plot initial stratigraphy
    if plot:
        plt.figure(figsize=(12, 8))

    for i in range(spindur):
        elevation[i, x_m + mwo - 1: B] = modernslope - (0.025 * (spindur - i - 1))
        if plot:  # Plot surfaces
            plt.plot(x, elevation[i, :], c='tan')

    if plot:
        plt.plot(x, elevation[spindur - 1, :], c='black')
        plt.scatter(x[x_m], elevation[spindur - 1, x_m], c='green')
        plt.scatter(x[x_m + mwo], elevation[spindur - 1, x_m + mwo - 1], c='green')
        plt.xlabel("Distance [m]")
        plt.ylabel("Elevation Relative to Initial Sea Level [m]")
        plt.show()

        print("Initial Conditions:")
        print("RSLR = ", R, " mm/yr")
        print("Co = ", C, " kg/m")
        print("Wind = ", wind, " m/s")
        print("Bay Depth = ", round(100 * dfo) / 100, " m")
        print("db = ", dfo)
        print("Fetch = ", bfo, " m")
        print("Platform = ", mwo, " m")
        print("Upland Slope = ", slope)
        print("B = ", B)
        print("Upland Width = ", upland_width)
        print("mwo = ", mwo)
        print("slope = ", slope)
        print("maxY = ", maxY)


    return B, dfo, elevation
