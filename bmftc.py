"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _14 June 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import scipy.io
import math
import bisect

from buildtransect import buildtransect
from evolvemarsh import evolvemarsh


# from .functions import (
#     decompose,
# )

class Bmftc:
    def __init__(
            self,
            name="default",
            time_step=1,
            time_step_count=100,
            relative_sea_level_rise=4,
            reference_concentration=10,
            slope_upland=0.005,

            marsh_width_initial=1000,
            bay_fetch_initial=5000,
            forest_age_initial=60,

            # KV Organic
            bulk_density_mineral=2000,
            bulk_density_organic=2000,
            tidal_period=12.5 * 3600 * 1,
            settling_velocity_effective=0.05 * 10 ** (-3),
            settling_velocity_mudflat=0.5 * 10 ** (-3),
            critical_shear_mudflat=0.1,
            wind_speed=6,
            tidal_amplitude=1.4 / 2,
            marsh_progradation_coeff=2,
            marsh_erosion_coeff=0.16 / (365 * 24 * 3600),
            mudflat_erodibility_coeff=0.0001,
            dist_marsh_bank=10,
            tide_cycles_yearly=365 * (24 / 12.5),

            # Vegetation
            maximum_biomass_marsh=2500,
            veg_minimum_depth=0,
            maximum_biomass_forest=5000,
            forest_aboveground_unknown_a=4,
            forest_aboveground_unknown_b=2,
            forest_aboveground_unknown_f0=0.0001,
            forest_aboveground_unknown_fwet=5,
            forest_aboveground_unknown_fgrow=2,
            zero_decomposition_depth_marsh=0.4,
            decomposition_coefficient_marsh=0.1,

            # Bay/marsh
            tidal_iterations=500,
            unknown_fm_min=0,
            unknown_fm_org=0,
            unknown_fm_flood=0,
            sed_flux_pond=0,

    ):
        """Bay-Marsh-Forest Transect Carbon Model (Python version)

        Parameters
        ----------
        name: string, optional
            Name of simulation
        time_step: float, optional
            Time step of the numerical model [yr]

        Examples
        --------
        >>> from bmftc import Bmftc
        >>> model = Bmftc()
        """

        self._name = name
        self._RSLRi = relative_sea_level_rise  # mm/yr
        self._RSLR = relative_sea_level_rise * 10 ** (-3) / (3600 * 24 * 365)  # Convert from mm/yr to m/s
        self._time_index = 0
        self._dt = time_step
        self._dur = time_step_count
        self._Coi = reference_concentration  # mg/L
        self._Co = reference_concentration / 1000  # Convert to kg/m3
        self._slope = slope_upland

        self._mwo = marsh_width_initial
        self._bfo = bay_fetch_initial
        self._startforestage = forest_age_initial

        self._rhos = bulk_density_mineral
        self._rhoo = bulk_density_organic
        self._P = tidal_period
        self._ws = settling_velocity_effective
        self._wsf = settling_velocity_mudflat
        self._tcr = critical_shear_mudflat
        self._wind = wind_speed
        self._amp = tidal_amplitude
        self._Ba = marsh_progradation_coeff
        self._Be = marsh_erosion_coeff
        self._lamda = mudflat_erodibility_coeff
        self._dist = dist_marsh_bank
        self._cyclestep = tide_cycles_yearly

        self._BMax = maximum_biomass_marsh
        self._Dmin = veg_minimum_depth
        self._Bmax_forest = maximum_biomass_forest
        self._a = forest_aboveground_unknown_a
        self._b = forest_aboveground_unknown_b
        self._f0 = forest_aboveground_unknown_f0
        self._fwet = forest_aboveground_unknown_fwet
        self._fgrow = forest_aboveground_unknown_fgrow
        self._mui = zero_decomposition_depth_marsh  # [m] Depth below which decomposition goes to zero in the marsh
        self._mki = decomposition_coefficient_marsh  # Coefficient of decomposition in the marsh
        self._numiterations = tidal_iterations

        self._Fm_min = unknown_fm_min  # Unknown variable
        self._Fm_org = unknown_fm_org  # Unknown variable
        self._Fm_flood = unknown_fm_flood  # Unknown variable
        self._Fp_sum = sed_flux_pond  # Amount of sediment taken from ponds to recharge sedimentation to drowning interior marsh

        # Calculate additional variables
        self._SLR = self._RSLR * (3600 * 24 * 365)  # Convert to m/yr
        self._rhob = self._rhos  # [kg/m3] Bulk density of bay, which is initially all mineral
        self._tr = self._amp * 2  # [m] Tidal range
        self._Dmax = 0.7167 * 2 * self._amp - 0.483  # [m] Maximum depth below high water that marsh veg can grow

        # Load MarshStrat spin up file
        filename_spinup = "Input/MarshStrat_all_RSLR1_CO50.mat"  # IR: Need to make easily changeable
        marsh_spinup = scipy.io.loadmat(filename_spinup)
        self._elev25 = marsh_spinup["elev_25"]
        self._min_25 = marsh_spinup["min_25"]
        self._orgAL_25 = marsh_spinup["orgAL_25"]
        self._orgAT_25 = marsh_spinup["orgAT_25"]

        # Load Forest Organic Profile files: Look-up table with soil organic matter for forest based on age and depth
        directory_fop = "Input/Forest_Organic_Profile"
        self._forestOM = scipy.io.loadmat(directory_fop + "/forestOM.mat")  # [g] Table with forest organic matter profile stored in 25 depth increments of 2.5cm (rows) for
        # forests of different ages (columns) from 1 to 80 years
        self._forestMIN = scipy.io.loadmat(directory_fop + "/forestMIN.mat")  # [g] Table with forest mineral matter profile stored in 25 depth increments of 2.5cm (rows) for
        # forests of different ages (columns) from 1 to 80 years
        self._B_rts = scipy.io.loadmat(directory_fop + "/B_rts.mat")

        # Continue variable initializations
        self._startyear = np.size(self._elev25, axis=0)
        self._endyear = self._dur + self._startyear

        self._msl = np.zeros([self._endyear + 1])
        self._msl[self._startyear:self._endyear] = np.linspace(1, self._dur, num=self._dur) * self._SLR  # [m] Mean sea level over time relative to start

        # Time
        self._to = np.linspace(1, 3600 * 24 * 365 * 1, 2)

        # Initialize marsh and forest edge variables
        self._x_m = math.ceil(self._bfo) + 1  # First marsh cell
        self._Marsh_edge = np.zeros([self._endyear + 1])
        self._Marsh_edge[:self._startyear] = self._x_m
        self._Forest_edge = np.zeros(self._endyear + 1)
        self._fetch = np.zeros([self._endyear + 1])
        self._fetch[:self._startyear] = self._bfo

        self._tidal_dt = self._P / self._numiterations  # Inundation time?
        self._OCb = np.zeros(self._endyear + 1)  # Organic content of uppermost layer of bay sediment, which determines the organic content of suspended material deposited onto the
        # marsh. Initially set to zero.
        self._OCb[:551] = 0.15  # IR 6/8: Appears hardwired; need to fix
        self._edge_flood = np.zeros(self._endyear + 1)  # IR 6/8: Undefined variable

        self._marshOM_initial = (np.sum(np.sum(self._orgAL_25)) + np.sum(np.sum(self._orgAT_25))) / 1000  # [kg] Total mass of organic matter in the marsh at the beginning of
        # the simulation (both alloch and autoch)
        self._marshMM_initial = np.sum(np.sum(self._min_25)) / 1000  # [kg] Total mass of mineral matter in the marsh at the beginning of the simulation
        self._marshLOI_initial = self._marshOM_initial / (self._marshOM_initial + self._marshMM_initial) * 100  # [%] LOI of the initial marsh deposit
        self._marshOCP_initial = 0.4 * self._marshLOI_initial + 0.0025 * self._marshLOI_initial ** 2  # [%] Organic carbon content from Craft et al. (1991)
        self._marshOC_initial = self._marshOCP_initial / 100 * (self._marshOM_initial + self._marshMM_initial)  # [kg] Organic carbon deposited in the marsh over the past 25 (?
        # 30?) years

        # Build starting transect
        self._B, self._db, self._elevation = buildtransect(self._RSLRi, self._Coi, self._slope, self._mwo, self._elev25, self._amp, self._wind, self._bfo, self._endyear,
                                                           plot=False)

        # Set up vectors for deposition
        self._organic_dep_alloch = np.zeros([self._endyear + 1, self._B])
        self._organic_dep_autoch = np.zeros([self._endyear + 1, self._B])
        self._mineral_dep = np.zeros([self._endyear + 1, self._B])
        self._organic_dep_alloch[:self._startyear, self._x_m: self._x_m + self._mwo] = self._orgAL_25  # Set the first 25[?] years to be the spin up values for deposition
        self._organic_dep_autoch[:self._startyear, self._x_m: self._x_m + self._mwo] = self._orgAT_25
        self._mineral_dep[:self._startyear, self._x_m: self._x_m + self._mwo] = self._min_25

        # Run ODE
        # OPT = odeset...

        # Calculate where elevation is right for the forest to start
        self._Forest_edge[self._startyear - 1] = bisect.bisect_left(self._elevation[self._startyear - 1, :], self._msl[self._startyear - 1] + self._amp + self._Dmin)
        self._forestage = self._startforestage

        self._Bay_depth = np.zeros([self._endyear + 1])
        self._Bay_depth[:self._startyear] = self._db
        self._dmo = self._elevation[self._startyear - 1, self._x_m]  # Set marsh edge depth to the elevation of the marsh edge at year 25[?]

        # Initialize additional data storage arrays
        self._mortality = np.zeros([self._endyear + 1, self._B])
        self._BayExport = np.zeros([self._endyear + 1, 2])
        self._BayOM = np.zeros([self._endyear + 1])
        self._BayMM = np.zeros([self._endyear + 1])
        self._fluxes = np.zeros([8, self._endyear + 1])
        self._bgb_sum = np.zeros([self._endyear + 1])  # [g] Sum of organic matter deposited across the marsh platform in a given year
        self._Fd = np.zeros([self._endyear + 1])  # [kg] Flux of organic matter out of the marsh due to decomposition
        self._avg_accretion = np.zeros([self._endyear + 1])  # [m/yr] Annual accretion rate averaged across the marsh platform

    def update(self):
        """Update Bmftc by a single time step"""

        # Increase time
        self._time_index += 1

        # Calculate the density of the marsh edge cell
        boundyr = bisect.bisect_left(self._elevation[:, self._x_m], self._elevation[self._time_index - 1, 0])
        if boundyr == 0:
            us = self._elevation[0, self._x_m] - self._elevation[self._time_index - 1, 0]  # [m] Depth of underlying stratigraphy
            usmass = us * self._rhos  # [kg] Mass of pure mineral sediment underlying marsh at marsh edge
        else:
            usmass = 0  # [kg] Mass of pure mineral sediment underlying marsh at marsh edge

        # Mass of sediment to be eroded at the current marsh edge above the depth of erosion [kg]
        massm = np.sum(self._organic_dep_autoch[:, self._x_m]) / 1000 + np.sum(self._organic_dep_alloch[:, self._x_m]) / 1000 + np.sum(self._mineral_dep[:, self._x_m]) / 1000 + usmass
        # Volume of sediment to be eroded at the current marsh edge above the depth of erosion [m3]
        volm = self._elevation[self._time_index - 1, self._x_m] - self._elevation[self._time_index - 1, 0]

        # IR 14Jun21: Returning wrong values for volm, has to do with indexing of year
        if self._time_index == 1:
            print("boundyr = ", boundyr)
            print("us = ", us)
            print("usmass = ", usmass)
            print("massm = ", massm)
            print("volm = ", volm)
            print()
            print("1: ", self._elevation[self._time_index - 1, self._x_m])
            print("2: ", self._elevation[self._time_index - 1, 0])
            print()

    @property
    def time_index(self):
        return self._time_index

    @property
    def dur(self):
        return self._dur
