"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _15 June 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""


import time
import numpy as np
import matplotlib.pyplot as plt

from bmftc import Bmftc

# Create an instance of the BMI class
bmftc = Bmftc(
            name="default",
            time_step=1,
            time_step_count=10,
            relative_sea_level_rise=1,
            reference_concentration=10,
            slope_upland=0.005,
)

# ==================================================================================================================================================================================
# Run the PyBMFT-C model
Time = time.time()  # Record start time
for time_step in range(int(bmftc.dur)):
    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

    # Run time step
    bmftc.update()

# Print elapsed time of simulation
print()
SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")

# ==================================================================================================================================================================================
# Sum major fluxes and output variables for analysis

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

print()
print("SUM ORG", np.sum(organic_dep_last30yrs))
print("SUM MIN", np.sum(mineral_dep_last30yrs))

plt.figure()
plt.plot(organic_dep_last30yrs)
plt.xlabel("Year (previous 30)")
plt.ylabel("Organic Deposition [g]")
plt.show()

plt.figure()
plt.plot(mineral_dep_last30yrs)
plt.xlabel("Year (previous 30)")
plt.ylabel("Mineral Deposition [g]")
plt.show()

# plt.figure()
# plt.plot(bmftc.organic_dep_autoch[bmftc.endyear - 30: bmftc.endyear + 1, bmftc.x_m: bmftc.x_f + 1])
# plt.xlabel("Distance")
# plt.ylabel("Mineral Deposition Autoch [g]")
# plt.show()

# print()
# print("SUM ORG", np.sum(bmftc.organic_dep_autoch[bmftc.endyear - 30: bmftc.endyear + 1, bmftc.x_m: bmftc.x_f + 1]))
# print("SUM MIN", np.sum(bmftc.organic_dep_alloch[bmftc.endyear - 30: bmftc.endyear + 1, bmftc.x_m: bmftc.x_f + 1]))

# plt.figure()
# plt.plot(bmftc.organic_dep_alloch[bmftc.endyear - 30: bmftc.endyear + 1, bmftc.x_m: bmftc.x_f + 1])
# plt.xlabel("Distance")
# plt.ylabel("Mineral Deposition Alloch [g]")
# plt.show()

plt.figure()
plt.plot(bmftc.elevation[bmftc.endyear - 1, :])
plt.xlabel("Distance")
plt.ylabel("Elevation [m MSL]")
plt.show()
