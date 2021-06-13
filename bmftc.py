import numpy as np
import scipy.io
import math


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
        self._RSLR = relative_sea_level_rise * 10 ** (-3) / (3600 * 24 * 365)  # Convert from mm/yr to m/s
        self._time_index = 0
        self._dt = time_step
        self._dur = time_step_count
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
        self._mui = zero_decomposition_depth_marsh
        self._mki = decomposition_coefficient_marsh
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
        filename_spinup = "/SpinUp/MarshStrat_all_RSLR1_CO50.mat"  # IR: Need to make easily changeable
        marsh_spinup = scipy.io.loadmat(filename_spinup)
        self._elev25 = marsh_spinup["elev_25"]
        self._min_25 = marsh_spinup["min_25"]
        self._orgAL_25 = marsh_spinup["orgAL_25"]
        self._orgAT_25 = marsh_spinup["orgAT_25"]

        # Continue variable initializations
        self._startyear = np.size(self._elev25, axis=0) + 1
        self._endyear = self._dur + self._startyear

        self._msl = np.zeros([self._endyear])
        self._msl[self._startyear:self._endyear] = np.linspace(1, self._dur, num=self._dur) * self._SLR  # [m] Mean sea level over time relative to start

        # Time
        self._to = np.linspace(1, 3600*24*365*1, 2)

        # Initialize marsh and forest edge variables
        self._x_m = math.ceil(self._bfo) + 1  # First marsh cell
        self._Marsh_edge = np.zeros([self._endyear])
        self._Marsh_edge[:self._startyear] = self._x_m
        self._Forest_edge = np.zeros(self._endyear)
        self._fetch = np.zeros([self._endyear])
        self._fetch[:self._startyear] = self._bfo

        self._tidal_dt = self._P / self._numiterations  # Inundation time?
        self._OCb = np.zeros(self._endyear)  # Organic content of uppermost layer of bay sediment, which determines the OC of suspended material deposited onto the marsh. Initially set to zero.
        self._OCb[:551] = 0.15  # IR 6/8: Appears hardwired; need to fix
        self._edge_flood = np.zeros(self._endyear)  # IR 6/8: Undefined variable

        self._marshOM_initial = (np.sum(np.sum(self._orgAL_25)) + np.sum(np.sum(self._orgAT_25))) / 1000  # [kg] Total mass of organic matter in the marsh at the beginning of the simulation (both alloch and autoch)
        self._marshMM_initial = np.sum(np.sum(self._min_25)) / 1000  # [kg] Total mass of mineral matter in the marsh at the beginning of the simulation
        self._marshLOI_initial = self._marshOM_initial / (self._marshOM_initial + self._marshMM_initial) * 100  # [%] LOI of the initial marsh deposit
        self._marshOCP_initial = 0.4 * self._marshLOI_initial + 0.0025 * self._marshLOI_initial ** 2  # [%] Organic carbon content from Craft et al. (1991)
        self._marshOC_initial = self._marshOCP_initial / 100 * (self._marshOM_initial + self._marshMM_initial)  # [kg] Organic carbon deposited in the marsh over the past 25 (? 30?) years


    def update(self):
        """Update Bmftc by a single time step"""

        self._time_index += 1

    @property
    def time_index(self):
        return self._time_index
