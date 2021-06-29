"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _23 June 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import math


def decompose(
        x_m,
        x_f,
        yr,
        organic_dep_autoch,
        elevation,
        B,
        mui,
        mki,
        rhoo,
):
    """Decomposes all of the organic sediment within the marsh soil profile at a rate determined by depth."""

    compaction = np.zeros([B])
    Fd = 0

    # IR 26Jun21 00:20 Disconnects with Matlab version here

    # Decompose the marsh sediment
    for x in range(x_m, x_f):  # Loop through each marsh and upland cell in the domain
        decomp = np.zeros([yr + 1])
        for tempyr in range(yr, 0, -1):  # Loop through each pocket of sediment in each cell, starting at the most recently deposited
            # packet of sediment at the surface
            depth = elevation[yr, x] - elevation[tempyr, x]  # Depth of sediment pocket below the surface
            print()
            if depth > mui:  # Maximum depth at which decomposition occurs
                print()
                decomp[tempyr] = 0
                break
            else:
                decomp[tempyr] = organic_dep_autoch[tempyr, x] * (mki * math.exp(-depth / mui))  # [g] Mass of organic material
                # decomposed from a given "pocket" of sediment
                organic_dep_autoch[tempyr, x] -= decomp[tempyr]  # [g] Autochthanous organic material in a
                # given "pocket" of sediment updated for deomposition
                print()
        compaction[x] = np.sum(decomp) / 1000 / rhoo  # [m] Total compaction in a given cell is a result of the sum of all decomposition
        # in that cell
        Fd += np.sum(decomp)  # [kg] Flux of organic matter out of the marsh due to decomposition

    return compaction, Fd
