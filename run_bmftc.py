"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _20 August 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import time
import numpy as np
import matplotlib.pyplot as plt

from bmftc import Bmftc

# ==================================================================================================================================================================================
# Create an instance of the BMI class
bmftc = Bmftc(
            name="default",
            time_step=1,
            time_step_count=50,
            relative_sea_level_rise=1,
            reference_concentration=50,
            slope_upland=0.005,
            bay_fetch_initial=3000,
            wind_speed=4,
            seagrass_on=True,
)

# ==================================================================================================================================================================================
# Run the PyBMFT-C model

# Initialize in-run figure plotting
fig1 = plt.figure()
ax1 = fig1.add_subplot(2, 2, 1)
plt.xlabel("Distance [m]")
plt.ylabel("Organic Deposition Autoch [g]")

ax2 = fig1.add_subplot(2, 2, 3)
plt.xlabel("Distance [m]")
plt.ylabel("Mineral Deposition [g]")

ax3 = fig1.add_subplot(2, 2, 4)
plt.xlabel("Distance [m]")
plt.ylabel("Elevation [m]")

ax4 = fig1.add_subplot(2, 2, 2)
plt.xlabel("Distance [m]")
plt.ylabel("Organic Deposition Alloch [g]")

# Record start time
Time = time.time()

# Loop through time
for time_step in range(int(bmftc.dur)):

    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

    # Run time step
    bmftc.update()

    # Plot each time step
    ax1.plot(bmftc.organic_dep_autoch[550 + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))
    ax2.plot(bmftc.mineral_dep[550 + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))
    ax3.plot(bmftc.elevation[550 + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))
    ax4.plot(bmftc.organic_dep_alloch[550 + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))

# Print elapsed time of simulation
print()
SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")


# ==================================================================================================================================================================================
# Sum Major Fluxes and Output Variables for Analysis

# Organic matter deposited in the marsh over the past 30 years [g]
organic_dep_last30yrs = bmftc.organic_dep_autoch[bmftc.endyear - 31: bmftc.endyear, bmftc.x_m: bmftc.x_f + 1] + bmftc.organic_dep_alloch[bmftc.endyear - 31: bmftc.endyear, bmftc.x_m: bmftc.x_f + 1]
# Mineral matter deposited in the marsh over the past 30 years [g]
mineral_dep_last30yrs = bmftc.mineral_dep[bmftc.endyear - 31: bmftc.endyear, bmftc.x_m: bmftc.x_f + 1]
# # Percent organic matter [%]
loi_last30yrs = organic_dep_last30yrs / (organic_dep_last30yrs + mineral_dep_last30yrs) * 100
# Organic carbon content (%) from Craft et al (1991)
OCP_last30yrs = 0.4 * loi_last30yrs + 0.0025 * loi_last30yrs ** 2
# Organic Carbon deposited in the marsh over the past 30 years [g]
OC_last30yrs = OCP_last30yrs / 100 * (organic_dep_last30yrs + mineral_dep_last30yrs)
# Average Organic Carbon accumulation rate over the last 30 years [g C / m2 / yr]
OC_avg_last30yrs = np.mean(OC_last30yrs, axis=1)

# Total mass of organic matter in the marsh at the end of the simulation [kg]
marshOM_final = np.sum(np.sum(bmftc.organic_dep_autoch[:, bmftc.x_m: bmftc.x_f + 1]) + np.sum(np.sum(bmftc.organic_dep_alloch[:, bmftc.x_m: bmftc.x_f + 1]))) / 1000
# Total mass of mineral matter in the marsh at the end of the simulation [kg]
marshMM_final = np.sum(np.sum(bmftc.mineral_dep[:, bmftc.x_m: bmftc.x_f + 1])) / 1000
# Average loi of the marsh across the marsh platform at the end of the simulation [%]
marshLOI_final = marshOM_final / (marshOM_final + marshMM_final) * 100
# Organic carbon content (%) from Craft et al (1991)
marshOCP_final = 0.4 * marshLOI_final + 0.0025 * marshLOI_final ** 2
# Mass of organic carbon stored in the marsh at the end of the simulation [kg]
marshOC_final = marshOCP_final / 100 * (marshOM_final + marshMM_final)


# ==================================================================================================================================================================================
# Plot

plt.figure()
plt.plot(organic_dep_last30yrs)
plt.xlabel("Year (previous 30)")
plt.ylabel("Organic Deposition [g]")

plt.figure()
plt.plot(mineral_dep_last30yrs)
plt.xlabel("Year (previous 30)")
plt.ylabel("Mineral Deposition [g]")

# plt.figure()
# plt.plot(bmftc.organic_dep_autoch[bmftc.endyear - 30: bmftc.endyear + 1, bmftc.x_m: bmftc.x_f + 1])
# plt.xlabel("Distance")
# plt.ylabel("Organic Deposition Autoch [g]")

# plt.figure()
# plt.plot(bmftc.organic_dep_alloch[bmftc.endyear - 30: bmftc.endyear + 1, bmftc.x_m: bmftc.x_f + 1])
# plt.xlabel("Distance")
# plt.ylabel("Organic Deposition Alloch [g]")

plt.figure()
plt.plot(bmftc.elevation[bmftc.endyear - 1, :])
plt.xlabel("Distance")
plt.ylabel("Elevation [m MSL]")

# Seagrass: bay depth and shoot density
plt.figure()
fig = plt.gcf()
fig.set_size_inches(12, 14)
plt.rcParams.update({"font.size": 12})

plt.subplot(2,1,1)
plt.plot(bmftc.Bay_depth[bmftc.startyear:] * -1)
plt.xlabel("Year")
plt.ylabel("Bay Depth [m MSL]")

plt.subplot(2,1,2)
plt.plot(bmftc.seagrass[bmftc.startyear: , 1])
plt.xlabel("Year")
plt.ylabel("Shoot Density [shoots/m^2]")

plt.show()
