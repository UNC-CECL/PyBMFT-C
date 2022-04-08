"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _20 January 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import time
import numpy as np
import matplotlib.pyplot as plt
import warnings

from bmftc import Bmftc

warnings.filterwarnings("ignore")

# ==================================================================================================================================================================================
# Create an instance of the BMI class
bmftc = Bmftc(
            name="default",
            time_step_count=125,
            relative_sea_level_rise=8,
            reference_concentration=50,
            slope_upland=0.005,
            bay_fetch_initial=5000,
            wind_speed=6,
            seagrass_on=False,
            forest_on=True,
            # filename_equilbaydepth="Input/PyBMFT-C/EquilibriumBayDepth_f3000_w5.mat",
            filename_equilbaydepth="Input/PyBMFT-C/Equilibrium Bay Depth.mat",
            # filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50_width500.mat",
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
            marsh_width_initial=1000
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
    ax1.plot(bmftc.organic_dep_autoch[bmftc.startyear + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))
    ax2.plot(bmftc.mineral_dep[bmftc.startyear + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))
    ax3.plot(bmftc.elevation[bmftc.startyear + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))
    ax4.plot(bmftc.organic_dep_alloch[bmftc.startyear + time_step, bmftc.x_m: bmftc.x_f + 1], label=str(time_step))

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
for t in range(bmftc.startyear, bmftc.endyear, 10):
    plt.plot(bmftc.elevation[t, :])
plt.xlabel("Distance")
plt.ylabel("Elevation [m MSL]")

# Seagrass: bay depth and shoot density
plt.figure()
fig = plt.gcf()
fig.set_size_inches(12, 14)
plt.rcParams.update({"font.size": 12})

plt.subplot(3,1,1)
plt.plot(bmftc.Bay_depth[bmftc.startyear:] * -1)
plt.xlabel("Year")
plt.ylabel("Bay Depth [m MSL]")

plt.subplot(3,1,2)
plt.plot(bmftc.seagrass[bmftc.startyear: , 1])
plt.xlabel("Year")
plt.ylabel("Shoot Density [shoots/m^2]")

plt.subplot(3,1,3)
plt.plot(bmftc.Marsh_edge[bmftc.startyear:])
plt.xlabel("Year")
plt.ylabel("Marsh Edge Location [m cross-shore]")

# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(7, 15)

plt.subplot(4, 1, 1)
plt.plot(bmftc.OCb[bmftc.startyear: bmftc.endyear + 1])
plt.xlabel("Time [yr]")
plt.ylabel("Bay Org Cont")

plt.subplot(4, 1, 2)
plt.plot(bmftc.BaySedDensity[:bmftc.endyear + 1])
plt.xlabel("Time [yr]")
plt.ylabel("Bay Sed Dens (rhob)")

plt.subplot(4, 1, 3)
plt.plot(bmftc.rhomt[:bmftc.endyear + 1])
plt.xlabel("Time [yr]")
plt.ylabel("Mar Edg Dens (rhom)")

plt.subplot(4, 1, 4)
plt.plot(bmftc.Bay_depth[bmftc.startyear: bmftc.endyear + 1])
plt.xlabel("Time [yr]")
plt.ylabel("Bay Depth")

# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.rcParams.update({"font.size": 12})

# Fe_min
plt.subplot(2, 4, 1)
plt.plot(bmftc.fluxes[0, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fe_min: marsh to bay")

# Fe_org
plt.subplot(2, 4, 2)
plt.plot(bmftc.fluxes[1, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fe_org: marsh to bay")

# Fm_min
plt.subplot(2, 4, 3)
plt.plot(bmftc.fluxes[2, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fm_min: bay to marsh")

# Fm_org
plt.subplot(2, 4, 4)
plt.plot(bmftc.fluxes[3, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fm_org: bay to marsh")

# Fc_min
plt.subplot(2, 4, 5)
plt.plot(bmftc.fluxes[4, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fc_min: external to bay")

# Fc_org
plt.subplot(2, 4, 6)
plt.plot(bmftc.fluxes[5, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fc_org: external to bay")

# Fb_min
plt.subplot(2, 4, 7)
plt.plot(bmftc.fluxes[6, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fb_min: net flux into bay")

# Fb_org
plt.subplot(2, 4, 8)
plt.plot(bmftc.fluxes[7, bmftc.startyear: bmftc.endyear + 1])
plt.ylabel("Fb_org: net flux into bay")
plt.tight_layout()

plt.show()
