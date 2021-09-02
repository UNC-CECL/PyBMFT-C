"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _1 September 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import time
import numpy as np
import matplotlib.pyplot as plt
import warnings

from bmftc import Bmftc

warnings.filterwarnings("ignore")

# ==================================================================================================================================================================================
# Create an instance of the BMI class

# Mainland shoreline
bmftc_ML = Bmftc(
            name="default",
            time_step=1,
            time_step_count=50,
            relative_sea_level_rise=3,
            reference_concentration=50,
            slope_upland=0.005,
            bay_fetch_initial=3000,
            wind_speed=5,
            seagrass_on=False,
            critical_shear_mudflat=0.2,
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
)

# Back-barrier shoreline
bmftc_BB = Bmftc(
            name="default",
            time_step=1,
            time_step_count=50,
            relative_sea_level_rise=3,
            reference_concentration=50,
            slope_upland=0.001,
            bay_fetch_initial=3000,
            marsh_width_initial=500,
            wind_speed=5,
            seagrass_on=False,
            critical_shear_mudflat=0.2,
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_500_all_RSLR1_CO50.mat",
)


# ==================================================================================================================================================================================
# Run the PyBMFT-C model

# Record start time
Time = time.time()

start_yr_ML = bmftc_ML.startyear
start_yr_BB = bmftc_BB.startyear

# Loop through time
for time_step in range(int(bmftc_ML.dur)):

    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

    # Run time step
    bmftc_ML.update()
    bmftc_BB.update()

    # Set fetch, depth
    delta_fetch_ML = bmftc_ML.fetch[start_yr_ML + time_step] - bmftc_ML.fetch[start_yr_ML + time_step - 1]  # [m] Calculate change in fetch from erosion of mainland marsh
    delta_fetch_BB = bmftc_BB.fetch[start_yr_BB + time_step] - bmftc_BB.fetch[start_yr_BB + time_step - 1]  # [m] Calculate change in fetch from erosion of back-barrier marsh

    bmftc_ML._bfo = bmftc_ML.bfo + delta_fetch_BB
    bmftc_ML.fetch[start_yr_ML + time_step] = bmftc_ML.bfo
    bmftc_BB._bfo = bmftc_BB.bfo + delta_fetch_ML
    bmftc_BB.fetch[start_yr_BB + time_step] = bmftc_BB.bfo

    print("Fetch diff:", bmftc_ML.bfo - bmftc_BB.bfo)

# Print elapsed time of simulation
print()
SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")


# ==================================================================================================================================================================================
# Plot

# ===========================
# Elevation Profile
plt.figure()
fig = plt.gcf()
fig.set_size_inches(12, 4)
plt.rcParams.update({"font.size": 12})
ML_profile = bmftc_ML.elevation[bmftc_ML.endyear - 1, :]
BB_profile = np.flip(bmftc_BB.elevation[bmftc_BB.endyear - 1, :])
whole_profile = np.append(BB_profile, ML_profile)
plt.plot(whole_profile)
plt.xlabel("Distance")
plt.ylabel("Elevation [m MSL]")


# ===========================
# Seagrass: bay depth and shoot density
plt.figure()
fig = plt.gcf()
fig.set_size_inches(12, 14)
plt.rcParams.update({"font.size": 12})

plt.subplot(3, 2, 1)
plt.plot(bmftc_ML.Bay_depth[bmftc_ML.startyear:] * -1)
plt.xlabel("Year")
plt.ylabel("Bay Depth [m MSL]")
plt.title("Mainland")

plt.subplot(3, 2, 3)
plt.plot(bmftc_ML.seagrass[bmftc_ML.startyear: , 1])
plt.xlabel("Year")
plt.ylabel("Shoot Density [shoots/m^2]")

plt.subplot(3, 2, 5)
plt.plot(bmftc_ML.Marsh_edge[bmftc_ML.startyear:])
plt.xlabel("Year")
plt.ylabel("Marsh Edge Location [m cross-shore]")

plt.subplot(3, 2, 2)
plt.plot(bmftc_BB.Bay_depth[bmftc_BB.startyear:] * -1)
plt.xlabel("Year")
plt.ylabel("Bay Depth [m MSL]")
plt.title("Back-Barrier")

plt.subplot(3, 2, 4)
plt.plot(bmftc_BB.seagrass[bmftc_BB.startyear: , 1])
plt.xlabel("Year")
plt.ylabel("Shoot Density [shoots/m^2]")

plt.subplot(3, 2, 6)
plt.plot(bmftc_BB.Marsh_edge[bmftc_BB.startyear:])
plt.xlabel("Year")
plt.ylabel("Marsh Edge Location [m cross-shore]")

# ===========================
plt.show()
