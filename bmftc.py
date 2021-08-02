"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _5 July 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import scipy.io
from scipy.integrate import solve_ivp
import math
import bisect

from buildtransect import buildtransect
from funBAY import funBAY
from funBAY import POOLstopp5
from calcFE import calcFE
from evolvemarsh import evolvemarsh
from decompose import decompose


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
            bulk_density_organic=85,
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
        self._dur = time_step_count + 1
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
        self._elev25 = marsh_spinup["elev_25"]  # IR 29Jun21: very slightly off from Matlab version when importing (rounding error)
        self._min_25 = marsh_spinup["min_25"]
        self._orgAL_25 = marsh_spinup["orgAL_25"]
        self._orgAT_25 = marsh_spinup["orgAT_25"]

        # Load Forest Organic Profile files: Look-up table with soil organic matter for forest based on age and depth
        directory_fop = "Input/Forest_Organic_Profile"
        file_forestOM = scipy.io.loadmat(directory_fop + "/forestOM.mat")  # [g] Table with forest organic matter profile stored in 25 depth increments of 2.5cm (rows) for
        # forests of different ages (columns) from 1 to 80 years
        self._forestOM = file_forestOM["forestOM"]
        file_forestMIN = scipy.io.loadmat(directory_fop + "/forestMIN.mat")  # [g] Table with forest mineral matter profile stored in 25 depth increments of 2.5cm (rows) for
        self._forestMIN = file_forestMIN["forestMIN"]
        # forests of different ages (columns) from 1 to 80 years
        file_B_rts = scipy.io.loadmat(directory_fop + "/B_rts.mat")
        self._B_rts = file_B_rts["B_rts"]

        # Continue variable initializations
        self._startyear = np.size(self._elev25, axis=0)
        self._endyear = self._dur + self._startyear

        self._msl = np.zeros([self._endyear])
        self._msl[self._startyear:self._endyear] = np.linspace(1, self._dur, num=self._dur) * self._SLR  # [m] Mean sea level over time relative to start

        # Time
        self._to = np.linspace(0, 3600 * 24 * 365 * 1, 2)
        self._timestep = 365 * (24 / 12.5)  # [tidal cycles per year] number to multiply accretion simulated over a tidal cycle by

        # Initialize marsh and forest edge variables
        self._x_m = math.ceil(self._bfo)  # First marsh cell
        self._x_f = None  # First forest cell
        self._Marsh_edge = np.zeros([self._endyear])
        self._Marsh_edge[:self._startyear] = self._x_m
        self._Forest_edge = np.zeros(self._endyear)
        self._fetch = np.zeros([self._endyear])
        self._fetch[:self._startyear] = self._bfo

        self._tidal_dt = self._P / self._numiterations  # Inundation time?
        self._OCb = np.zeros(self._endyear)  # Organic content of uppermost layer of bay sediment, which determines the organic content of suspended material deposited onto the
        # marsh. Initially set to zero.
        self._OCb[:552] = 0.15  # IR 6/8: Appears hardwired; need to fix
        self._edge_flood = np.zeros(self._endyear)  # IR 6/8: Undefined variable
        self._Edge_ht = np.zeros(self._endyear)  # IR 6/24: Undefined variable

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
        self._organic_dep_alloch = np.zeros([self._endyear, self._B])
        self._organic_dep_autoch = np.zeros([self._endyear, self._B])
        self._mineral_dep = np.zeros([self._endyear, self._B])
        self._organic_dep_alloch[:self._startyear, self._x_m: self._x_m + self._mwo] = self._orgAL_25  # Set the first 25[?] years to be the spin up values for deposition
        self._organic_dep_autoch[:self._startyear, self._x_m: self._x_m + self._mwo] = self._orgAT_25
        self._mineral_dep[:self._startyear, self._x_m: self._x_m + self._mwo] = self._min_25

        # Set options for ODE solver
        POOLstopp5.terminal = True
        # POOLstopp5.direction = 1  # IR 15Jun21: Is this correct? Trouble translating from Matlab

        # Calculate where elevation is right for the forest to start
        self._Forest_edge[self._startyear - 1] = bisect.bisect_left(self._elevation[self._startyear - 1, :], self._msl[self._startyear - 1] + self._amp + self._Dmin)
        self._forestage = self._startforestage

        self._Bay_depth = np.zeros([self._endyear])
        self._Bay_depth[:self._startyear] = self._db
        self._dmo = self._elevation[self._startyear - 1, self._x_m]  # Set marsh edge depth to the elevation of the marsh edge at year 25[?]

        # Initialize
        self._C_e_ODE = []
        self._Fc_ODE = []

        # Initialize additional data storage arrays
        self._mortality = np.zeros([self._endyear, self._B])
        self._BayExport = np.zeros([self._endyear, 2])
        self._BayOM = np.zeros([self._endyear])
        self._BayMM = np.zeros([self._endyear])
        self._fluxes = np.zeros([8, self._endyear])
        self._bgb_sum = np.zeros([self._endyear])  # [g] Sum of organic matter deposited across the marsh platform in a given year
        self._Fd = np.zeros([self._endyear])  # [kg] Flux of organic matter out of the marsh due to decomposition
        self._avg_accretion = np.zeros([self._endyear])  # [m/yr] Annual accretion rate averaged across the marsh platform
        self._rhomt = np.zeros([self._dur])
        self._C_e = np.zeros([self._endyear])
        self._aboveground_forest = np.zeros([self._endyear, self._B])  # Forest aboveground biomass
        self._OM_sum_au = np.zeros([self._endyear, self._B])
        self._OM_sum_al = np.zeros([self._endyear, self._B])

    def update(self):
        """Update Bmftc by a single time step"""

        # Year including spinup
        yr = self._time_index + self._startyear

        # Calculate the density of the marsh edge cell
        boundyr = bisect.bisect_left(self._elevation[:, self._x_m], self._elevation[yr - 1, 0])
        if boundyr == 0:
            us = self._elevation[0, self._x_m] - self._elevation[yr - 1, 0]  # [m] Depth of underlying stratigraphy
            usmass = us * self._rhos  # [kg] Mass of pure mineral sediment underlying marsh at marsh edge
        else:
            usmass = 0  # [kg] Mass of pure mineral sediment underlying marsh at marsh edge

        # Mass of sediment to be eroded at the current marsh edge above the depth of erosion [kg]
        massm = np.sum(self._organic_dep_autoch[:, self._x_m]) / 1000 + np.sum(self._organic_dep_alloch[:, self._x_m]) / 1000 + np.sum(
            self._mineral_dep[:, self._x_m]) / 1000 + usmass
        # Volume of sediment to be eroded at the current marsh edge above the depth of erosion [m3]
        volm = self._elevation[yr - 1, self._x_m] - self._elevation[yr - 1, 0]

        rhom = massm / volm  # [kg/m3] Bulk density of marsh edge
        self._rhomt[self._time_index] = rhom

        Fm = (self._Fm_min + self._Fm_org) / (3600 * 24 * 365)  # [kg/s] Mass flux of both mineral and organic sediment from the bay to the marsh

        # Parameters to feed into ODE
        PAR = [
            self._rhos,
            self._P,
            self._B,
            self._wsf,
            self._tcr,
            self._Co,
            self._wind,
            self._Ba,
            self._Be,
            self._amp,
            self._RSLR,
            Fm,
            self._lamda,
            self._dist,
            self._dmo,
            self._rhob,
            rhom,
            self,
        ]

        # ODE solves for change in bay depth and width
        # IR 5July21: Small deviations in the solved values from the Matlab version (on the order of ~ 10^-4 to 10^-5)
        ode = solve_ivp(funBAY,
                        t_span=self._to,
                        y0=[self._bfo, self._db],
                        atol=10 ** (-6),
                        rtol=10 ** (-6),
                        method='BDF',
                        args=PAR,
                        )

        fetch_ODE = ode.y[0, :]
        db_ODE = ode.y[1, :]

        self._db = db_ODE[-1]  # Set initial bay depth of the bay to final depth from funBAY
        self._fetch[yr] = fetch_ODE[-1]  # Set initial bay width of the bay to final width from funBAY
        self._bfo = self._fetch[yr]  # Set initial bay width of the bay to final width from funBAY
        self._C_e[yr] = self._C_e_ODE[-1]

        Fc = self._Fc_ODE[-1] * 3600 * 24 * 365  # [kg/yr] Annual net flux of sediment out of/into the bay from outside the system
        Fc_org = Fc * self._OCb[yr - 1]  # [kg/yr] Annual net flux of organic sediment out of/into the bay from outside the system
        Fc_min = Fc * (1 - self._OCb[yr - 1])  # [kg/yr] Annual net flux of mineral sediment out of/into the bay from outside the system

        # Calculate the flux of organic and mineral sediment to the bay from erosion of the marsh
        Fe_org, Fe_min = calcFE(self._bfo, self._fetch[yr - 1], self._elevation, yr, self._organic_dep_autoch, self._organic_dep_alloch, self._mineral_dep, self._rhos)
        Fe_org /= 1000  # [kg/yr] Annual net flux of organic sediment to the bay due to erosion
        Fe_min /= 1000  # [kg/yr] Annual net flux of mineral sediment to the bay due to erosion

        Fb_org = Fe_org - self._Fm_org - Fc_org  # [kg/yr] Net flux of organic sediment into (or out of, if negative) the bay
        Fb_min = Fe_min - self._Fm_min - Fc_min  # [kg/yr] Net flux of mineral sediment into (or out of, if negative) the bay

        self._BayExport[yr, :] = [Fc_org, Fc_min]  # [kg/yr] Mass of organic and mineral sediment exported from the bay each year
        self._BayOM[yr] = Fb_org  # [kg/yr] Mass of organic sediment stored in the bay in each year
        self._BayMM[yr] = Fb_min  # [kg/yr] Mass of mineral sediment stored in the bay in each year

        if Fb_org > 0 and Fb_min > 0:
            self._OCb[yr] = Fb_org / (Fb_org + Fb_min) + 0.05  # BIG CHANGE HERE
        elif Fb_org > 0:
            self._OCb[yr] = 1
        elif Fb_min > 0:
            self._OCb[yr] = 0
        else:
            self._OCb[yr] = self._OCb[yr - 1]

        # If bay has eroded down to depth below initial bay bottom, there is only mineral sediment remaining
        if self._db < self._Bay_depth[0]:
            self._OCb[yr] = 0.05

        self._rhob = 1 / ((1 - self._OCb[yr]) / self._rhos + self._OCb[yr] / self._rhoo)  # [kg/m3] Density of bay sediment

        if int(self._bfo) <= 0:
            print("Marsh has eroded completely away")
            self._endyear = yr
            return  # Exit program

        self._x_m = math.ceil(self._bfo)  # New first marsh cell
        self._x_f = bisect.bisect_left(self._elevation[yr - 1, :], self._msl[yr] + self._amp - self._Dmin)  # New first forest cell

        tempelevation = self._elevation[yr - 1, self._x_m: self._x_f + 1]
        Dcells = self._Marsh_edge[yr - 1] - self._x_m  # Gives the change in the number of marsh cells

        if Dcells > 0:  # Prograde the marsh, with new marsh cells having the same elevation as the previous marsh edge
            tempelevation[0: Dcells] = self._elevation[yr - 1, self._Marsh_edge[yr - 1]]
            # Account for mineral and organic material deposited in new marsh cells  # IR 21Jun21: IS THIS TO-DO??

        self._msl[yr] = self._msl[yr - 1] + self._SLR
        self._elevation[yr, 0: self._x_m] = self._msl[yr] + self._amp - self._db  # All bay cells have the same depth

        # Mineral and organic marsh deposition
        (
            tempelevation,
            temporg_autoch,
            temporg_alloch,
            tempmin,
            self._Fm_min,
            self._Fm_org,
            tempbgb,
            accretion,
            tempagb,
        ) = evolvemarsh(
            tempelevation,
            self._msl[yr],
            self._C_e[yr],
            self._OCb[yr],
            self._tr,
            self._numiterations,
            self._P,
            self._tidal_dt,
            self._ws,
            self._timestep,
            self._BMax,
            self._Dmin,
            self._Dmax,
            self._rhoo,
            self._rhos,
            plot=False
        )

        self._elevation[yr, self._x_m: self._x_f + 1] = tempelevation  # [m] Set new elevation to current year
        self._elevation[yr, self._x_f + 1: self._B] = self._elevation[yr - 1, self._x_f + 1: self._B]  # Forest elevation remains unchanged
        self._mineral_dep[yr, self._x_m: self._x_f + 1] = tempmin  # [g] Mineral sediment deposited in a given year
        self._organic_dep_autoch[yr, self._x_m: self._x_f + 1] = temporg_autoch  # [g] Belowground plant material deposited in a given year
        self._mortality[yr, self._x_m: self._x_f + 1] = temporg_autoch  # [g] Belowground plant material deposited in a given year, for keeping track of without decomposition
        self._organic_dep_alloch[yr, self._x_m: self._x_f + 1] = temporg_alloch  # [g] Allochthonous organic material deposited in a given year
        self._bgb_sum[yr] = np.sum(tempbgb)  # [g] Belowground biomass deposition summed across the marsh platform. Saved through time without decomposition for analysis

        avg_accretion = np.mean(accretion)  # [m/yr] Accretion rate for a given year averaged across the marsh platform
        self._x_f = bisect.bisect_left(self._elevation[yr - 1, :], self._msl[yr] + self._amp - self._Dmin)  # New first forest cell

        # Update forest soil organic matter
        self._forestage += 1  # Age the forest
        for x in range(int(self._Forest_edge[yr - 1]), self._x_f + 1):
            if self._forestage < 80:
                self._organic_dep_autoch[self._startyear - 25: self._startyear, x] = self._forestOM[:, yr - 525] + self._B_rts[:, yr - 525]
            else:
                self._organic_dep_autoch[self._startyear - 25: self._startyear, x] = self._forestOM[:, 79] + self._B_rts[:, 79]
        for x in range(self._x_f, self._B):
            if self._forestage < 80:
                self._organic_dep_autoch[self._startyear - 25: self._startyear, x] = self._forestOM[:, yr - 525]
                self._mineral_dep[self._startyear - 25: self._startyear, x] = self._forestMIN[:, yr - 525]
            else:
                self._organic_dep_autoch[self._startyear - 25: self._startyear, x] = self._forestOM[:, 79]
                self._mineral_dep[self._startyear - 25: self._startyear, x] = self._forestMIN[:, 79]

        df = -self._msl[yr] + self._elevation[yr, self._x_f: self._B]

        self._organic_dep_autoch[yr, self._x_f: self._B] = self._f0 + self._fwet * np.exp(-self._fgrow * df)
        self._mineral_dep[yr, self._x_f: self._B] = self._forestMIN[0, 79]

        # Update forest aboveground biomass
        self._aboveground_forest[yr, self._x_f: self._B] = self._Bmax_forest / (1 + self._a * np.exp(-self._b * df))

        (
            compaction,
            tempFd,
        ) = decompose(
            self._x_m,
            self._x_f,
            yr,
            self._organic_dep_autoch,
            self._elevation,
            self._B,
            self._mui,
            self._mki,
            self._rhoo,
        )

        self._Fd[yr] = tempFd  # [kg] Flux of organic matter out of the marsh due to decomposition

        # Adjust marsh and forest elevation due to compaction from decomposition
        self._elevation[yr, self._x_m: self._B] -= compaction[self._x_m: self._B]
        self._OM_sum_au[yr, :len(self._elevation) + 1] = np.sum(self._organic_dep_autoch[:yr + 1, :])
        self._OM_sum_al[yr, :len(self._elevation) + 1] = np.sum(self._organic_dep_alloch[:yr + 1, :])

        F = 0
        while self._x_m < self._B:
            if self._organic_dep_autoch[yr, self._x_m] > 0:  # If organic deposition is greater than zero, marsh is no longer growing
                break
            else:  # Otherwise, the marsh has drowned, and will be eroded to form new bay
                F = 1
                self._edge_flood[yr] += 1  # Count that cell as a flooded cell
                self._bfo += 1  # Increase the bay fetch by one cell
                self._x_m += 1  # Update the new location of the marsh edge

        if F == 1:  # If flooding occurred, adjust marsh flux
            # Calculate the amount of organic and mineral sediment liberated from the flooded cells
            FF_org, FF_min = calcFE(self._bfo, self._fetch[yr - 1], self._elevation, yr, self._organic_dep_autoch, self._organic_dep_alloch, self._mineral_dep, self._rhos)
            # Adjust flux of mineral sediment to the marsh
            self._Fm_min -= FF_min
            # Adjust flux of organic sediment to the marsh
            self._Fm_org -= FF_org
            # Change the drowned marsh cell to z bay cell
            self._elevation[yr, :self._x_m + 1] = self._elevation[yr, 0]

        self._fluxes[:, yr] = [
            Fe_min,
            Fe_org,
            self._Fm_min,
            self._Fm_org,
            Fc_min,
            Fc_org,
            Fb_min,
            Fb_org,
        ]

        # Update inputs for marsh edge
        self._Marsh_edge[yr] = self._x_m
        self._Forest_edge[yr] = self._x_f
        self._Bay_depth[yr] = self._db

        if 0 < self._x_m < self._B:
            self._dmo = self._msl[yr] + self._amp - self._elevation[yr, self._x_m]
            self._Edge_ht[yr] = self._dmo
        elif self._x_m <= 0:  # Condition for if the marsh has expanded to fill the basin
            print("Marsh has expanded to fill the basin.")
            self._endyear = yr
            self._edge_flood[self._endyear + 1:] = []  # ?
            return  # Exit program
        elif self._x_m >= self._B:  # Condition for if the marsh has eroded completely away
            print("Marsh has retreated. Basin is completely flooded.")
            self._endyear = yr
            self._edge_flood[self._endyear + 1:] = []  # ?
            return  # Exit program

        if self._db < 0.3:  # Condition for if the bay gets very shallow. Should this number be calculated within the code?
            print("Bay has filled in to form marsh.")
            self._endyear = yr
            self._edge_flood[self._endyear + 1:] = []  # ?
            return  # Exit program

        self._fetch[yr] = self._bfo  # Save change in bay fetch through time

        self._Fc_ODE = []
        self._C_e_ODE = []

        # Increase time
        self._time_index += 1

        # TIME STEP COMPLETE
        # ==========================================================================================================================================================================

    @property
    def time_index(self):
        return self._time_index

    @property
    def dur(self):
        return self._dur

    @property
    def organic_dep_autoch(self):
        return self._organic_dep_autoch

    @property
    def x_m(self):
        return self._x_m

    @property
    def x_f(self):
        return self._x_f

    @property
    def organic_dep_alloch(self):
        return self._organic_dep_alloch

    @property
    def endyear(self):
        return self._endyear

    @property
    def mineral_dep(self):
        return self._mineral_dep

    @property
    def elevation(self):
        return self._elevation

    @property
    def B(self):
        return self._B
