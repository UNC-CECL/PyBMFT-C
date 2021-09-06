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
            time_step_count=30,
            relative_sea_level_rise=3,
            reference_concentration=50,
            slope_upland=0.005,
            bay_fetch_initial=3000,
            wind_speed=5,
            seagrass_on=False,
            forest_on=True,
            critical_shear_mudflat=0.2,
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
)

# Back-barrier shoreline
bmftc_BB = Bmftc(
            # name="default",
            # time_step=1,
            # time_step_count=30,
            # relative_sea_level_rise=3,
            # reference_concentration=50,
            # slope_upland=0.0005,
            # bay_fetch_initial=3000,
            # marsh_width_initial=500,  # Different initial widths will make the depths slightly (i.e. 1 * 10**-3 to 10**-4 m) off in spinup
            # wind_speed=5,
            # seagrass_on=False,
            # critical_shear_mudflat=0.2,
            # filename_marshspinup="Input/PyBMFT-C/MarshStrat_500_all_RSLR1_CO50.mat",
            name="default",
            time_step=1,
            time_step_count=30,
            relative_sea_level_rise=3,
            reference_concentration=50,
            slope_upland=0.005,
            bay_fetch_initial=3000,
            wind_speed=5,
            seagrass_on=False,
            forest_on=False,
            critical_shear_mudflat=0.2,
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
)


# ==================================================================================================================================================================================
# Run the PyBMFT-C model

# Record start time
Time = time.time()

x_b_TS_ML = np.zeros([bmftc_ML.dur])
x_b_TS_BB = np.zeros([bmftc_BB.dur])

sum_delta_fetch_ML = 0
sum_delta_fetch_BB = 0

# Loop through time
for time_step in range(int(bmftc_ML.dur)):

    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

    # ========================================
    # Run time step
    bmftc_ML.update()
    bmftc_BB.update()

    # ========================================
    # Set fetch, depth
    delta_fetch_ML = bmftc_ML.bfo - bmftc_ML.fetch[bmftc_ML.startyear + time_step - 1]  # [m] Calculate change in fetch from erosion of mainland marsh
    delta_fetch_BB = bmftc_BB.bfo - bmftc_BB.fetch[bmftc_BB.startyear + time_step - 1]  # [m] Calculate change in fetch from erosion of back-barrier marsh

    # Temp sum
    sum_delta_fetch_ML += delta_fetch_ML
    sum_delta_fetch_BB += delta_fetch_BB
    sumchange = delta_fetch_ML + delta_fetch_BB

    # Determine change in x_b location
    bmftc_ML._x_b = bmftc_ML.x_b - delta_fetch_BB
    bmftc_BB._x_b = bmftc_BB.x_b - delta_fetch_ML
    x_b_TS_ML[time_step] = bmftc_ML.x_b  # Save to array
    x_b_TS_BB[time_step] = bmftc_BB.x_b  # Save to array

    # Determine new fetch based on change in opposite marsh - both fetches should be exactly the same!
    bmftc_ML._bfo = bmftc_ML.bfo + delta_fetch_BB
    bmftc_BB._bfo = bmftc_BB.bfo + delta_fetch_ML
    bmftc_ML.fetch[bmftc_ML.startyear + time_step] = bmftc_ML.bfo  # Save to array
    bmftc_BB.fetch[bmftc_BB.startyear + time_step] = bmftc_BB.bfo  # Save to array
    # print("bfo2:", bmftc_ML.bfo)
    # print("bfo2-delta:", bmftc_ML.bfo - delta_fetch_BB)
    # print("-- sumchange:", sumchange)
    # print("-- NEW_delta_fetch_ML:", bmftc_ML.bfo - bmftc_ML.fetch[bmftc_ML.startyear + time_step - 1])

    # ========================================
    # Check consistency
    if bmftc_ML.db != bmftc_BB.db:
        print("Depths unequal!", bmftc_ML.db - bmftc_BB.db)
    if bmftc_ML.bfo != bmftc_BB.bfo:
        print("Depths unequal!", bmftc_ML.bfo - bmftc_BB.bfo)


# Print elapsed time of simulation
print()
SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")

print("sum_ML:", sum_delta_fetch_ML)
print("sum_BB:", sum_delta_fetch_BB)

# ==================================================================================================================================================================================
# Plot

# ===========================
# Elevation Profile
plt.figure()
fig = plt.gcf()
fig.set_size_inches(12, 4)
plt.rcParams.update({"font.size": 9})
ML_profile = bmftc_ML.elevation[bmftc_ML.endyear - 1, bmftc_ML.x_m:]
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
plt.rcParams.update({"font.size": 9})

plt.subplot(5, 2, 1)
plt.plot(bmftc_ML.Bay_depth[bmftc_ML.startyear:] * -1)
plt.ylabel("Bay Depth [m MSL]")
plt.title("Mainland")

plt.subplot(5, 2, 3)
plt.plot(bmftc_ML.seagrass[bmftc_ML.startyear: , 1])
plt.ylabel("Shoot Density [shoots/m^2]")

plt.subplot(5, 2, 5)
plt.plot(bmftc_ML.Marsh_edge[bmftc_ML.startyear:])
plt.ylabel("Marsh Edge Location [m cross-shore]")

plt.subplot(5, 2, 7)
plt.plot(bmftc_ML.fetch[bmftc_ML.startyear:])
plt.ylabel("Fetch")

plt.subplot(5, 2, 9)
plt.plot(x_b_TS_ML)
plt.xlabel("Year")
plt.ylabel("x_b Location")

plt.subplot(5, 2, 2)
plt.plot(bmftc_BB.Bay_depth[bmftc_BB.startyear:] * -1)
plt.ylabel("Bay Depth [m MSL]")
plt.title("Back-Barrier")

plt.subplot(5, 2, 4)
plt.plot(bmftc_BB.seagrass[bmftc_BB.startyear: , 1])
plt.ylabel("Shoot Density [shoots/m^2]")

plt.subplot(5, 2, 6)
plt.plot(bmftc_BB.Marsh_edge[bmftc_BB.startyear:])
plt.ylabel("Marsh Edge Location [m cross-shore]")

plt.subplot(5, 2, 8)
plt.plot(bmftc_BB.fetch[bmftc_BB.startyear:])
plt.ylabel("Fetch")

plt.subplot(5, 2, 10)
plt.plot(x_b_TS_BB)
plt.xlabel("Year")
plt.ylabel("x_b Location")

# ===========================
plt.show()
