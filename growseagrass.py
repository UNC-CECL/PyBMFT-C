"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _20 August 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""


def growseagrass(
        seagrass,
        seagrass_density_table,
        elevation,
        x_m,
        msl,
        amp,
        yr,
):
    """Grows seagrass in specified bay cells. Loops through each cell in bay, finds depth, grows seagrass using lookup table, and returns seagrass array with each bay cell
    containing a shoot density. Based off of Reeves et al. (2020)."""

    bay = elevation[yr - 1, :x_m]

    # Define spatial limit of bay (i.e., buffer between seagrass meadow and shorelines) - Percent Bay Cover (PBC) in Reeves et al. (2020)
    barrier_limit = 0
    mainland_limit = len(bay)

    if len(bay) >= 10:  # Do not grow seagrass if bay is less than 10 m wide
        for i in range(x_m):  # Loop through each cell in the bay

            # Find depth of cell
            depth = msl[yr] + amp - bay[i]  # [m]
            depth = round(depth / 0.05) * 0.05  # [m] Round to nearest 0.05 m

            # Grow seagrass only if bay cell is within depth and spatial range for seagrass growth
            if barrier_limit <= i <= mainland_limit and 0 < depth < 4:
                if seagrass[yr - 1, i] > 0:  # If seagrass was present in cell in previous year
                    row = int(depth / 0.05)
                    shootdensity = seagrass_density_table[row, 1]  # Use lookup table to find density at specific depth for prior seagrass
                else:  # If cell did not contain seagrass in prior timestep
                    row = int(depth / 0.05)
                    shootdensity = seagrass_density_table[row, 2]  # Use lookup table to find density at specific depth for prior bare bed

                seagrass[yr, i] = shootdensity

    return seagrass
