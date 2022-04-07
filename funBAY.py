"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _20 August 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import numpy as np
import math


def funBAY(t,
           X,
           rhos,
           P,
           B,
           wsf,
           tcr,
           Co,
           wind,
           Ba,
           Be,
           amp,
           RSLR,
           Fm2,
           lamda,
           dist,
           dmo,
           rhob,
           rhom,
           seagrass_on,
           seagrass,
           self
           ):
    """Determines change in bay depth and width by solving mass balance between fluxes of sediment into and out of the bay from marsh edge erosion, tidal exchange with the
    outside sediment source, and sediment deposited onto the marsh surface """

    # TEMP SEAGRASS; to do: move variables to main file
    seagrass_max_density = 667  # [shoots/m] Maximum shoot density a cell can achieve
    seagrass_meadow_width = np.count_nonzero(seagrass)  # [m] Width of seagrass meadow
    if seagrass_meadow_width == 0:
        max_density_pct = 0
    else:
        seagrass_meadow_density = np.mean(seagrass[seagrass > 0])  # Average shoot density of seagrass meadow
        max_density_pct = seagrass_meadow_density / seagrass_max_density  # Percent of max shoot density
    max_decay_coeff = 0.01  # Maximum wave decay coefficient
    effective_decay_coeff = max_decay_coeff * max_density_pct  # Effective way decay coefficient (after adjusting for seagrass meadows with shoot density below maximum); decay coeff varies roughly +1:1 with shoot density (e.g. Manca et al., 2012)
    max_attenuation_pct = 0.45  # [%]  Maximum percent of wave height that can be attenuated - 45-70% from Hansen & Reidenbach (2012)

    # Dynamic Variable X
    fetch = X[0]  # Mudflat width
    df = X[1]  # Mudflat depth

    fac = min(1, df / (2 * amp))  # Proportion of tide that the bay is flooded, dimensionless
    Df = (df + (df - fac * 2 * amp)) / 2  # [m] Average bay depth over tidal cycle
    dm = dmo  # [m] Marsh edge depth

    tw = wavetau(fetch, wind, Df, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct)  # Calculates wave bed shear stress [Pa]

    tau = max((tw - tcr) / tcr, 0) * lamda  # Excess shear stress, dimensionless
    Cr = rhos * tau / (1 + tau)  # Reference suspended sediment concentration in the basin [kg/m3]

    try:
        hb = dm + (df - dm) * (1 - math.exp(-dist * 0.1 / df))  # [m] scarp height at a fixed distance from the marsh according to semi-empirical shoaling profile
    except OverflowError:
        print("  <-- hb error, dm:", dm, ", df:", df)
        df = self.db
        hb = dm + (df - dm) * (1 - math.exp(-dist * 0.1 / df))

    try:
        W = waveTRNS(amp, wind, fetch, hb, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct)  # [W] Wave power density at the marsh boundary
    except OverflowError:
        print("  <-- W overflow error, dm:", dm, ", df:", df)
        df = self.db
        hb = dm + (df - dm) * (1 - math.exp(-dist * 0.1 / df))
        W = waveTRNS(amp, wind, fetch, hb, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct)

    E = (Be * W / (hb - dm) - Ba * Cr * wsf / rhom)  # (m2/s) Net flux of sediment eroded from/deposited to the marsh edge

    Fc = (Cr - Co) * (fac * 2 * amp) / P / rhob  # (m2/s) Net flux of sediment lost or gained through tidal exchange with external sediment supply/sink

    # IR 5July21: Global variables! Could be improved by passing the information in a non-global way
    self._Fc_ODE.append(Fc * rhob * fetch)  # Save Fc as a mass flux (kg/s) for each iteration of the ODE
    self._C_e_ODE.append(Cr)  # Save C_e (SSC at marsh edge, kg/m3) for each iteration of the ODE to use in marsh model

    dX = np.zeros([2])
    dX[0] = E  # [m2/s, or m/s if integrated over 1m transect width] Change in bay width due to erosion
    dX[1] = -E * (df - dm) / fetch + Fm2 / fetch / rhos + Fc + RSLR  # [m/s] Change in bay depth due to mass balance between fluxes into and out of bay

    return dX


def wavetau(fetch, wind, Df, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct):
    """Calculates wave bed shear stress, as a function of fetch, wind speed, and bay depth. From Mariotti and Fagherazzi (2013)."""

    (Hs, Tp) = YeV(fetch, wind, Df, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct)  # Significant wave height (Hs, [m]) and wave period (Tp, [m])
    kk = wavek(1 / Tp, Df)  # Calculates wave number [m^-1]
    Um = math.pi * Hs / Tp / math.sinh(kk * Df)  # Term in equation for shear stress [m/s]
    aw = Tp * Um / (2 * math.pi)  # Term in equation for shear stress[m]
    ko = 0.001  # Roughness [m] (Mariotti and Fagherazzi, 2013)
    fw = 0.4 * (aw / ko) ** (-0.75)  # Friction factor, dimensionlesss
    tw = 1 / 2 * 1020 * fw * Um ** 2  # Wave bed shear stress [Pa]

    return tw


def YeV(fetch, wind, h, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct):
    """Calculates wave height (Hs) and period (Tp) for a given fetch, wind speed, and depth; based on a set of semi-empirical equations from Young and Verhagen (1996). Modified for eefect of seagrass. """

    g = 9.8  # [m/s2] Acceleration due to gravity
    delta = h * g / wind ** 2  # Dimensionless coefficient
    chi = fetch * g / wind ** 2  # Dimensionless coefficient
    epsilon = 3.64 * 10 ** (-3) * (
            math.tanh(0.493 * delta ** 0.75) * math.tanh(3.13 * 10 ** (-3) * chi ** 0.57 / math.tanh(0.493 * delta ** 0.75))) ** 1.74  # Dimensionless coefficient
    ni = 0.133 * (math.tanh(0.331 * delta ** 1.01) * math.tanh(5.215 * 10 ** (-4) * chi ** 0.73 / math.tanh(0.331 * delta ** 1.01))) ** (-0.37)  # Dimensionless coefficient
    Tp = wind / ni / g  # [s] Wave period
    Hs = 4 * math.sqrt(wind ** 4 * epsilon / g ** 2)  # [m] Wave height (without seagrass attenuation)
    if seagrass_on and seagrass_meadow_width > 0:
        min_Hs = Hs * (1 - max_attenuation_pct)  # Minimum wave weight allowed, prevents wave height from approaching zero
        effective_Hs = Hs * math.exp(-effective_decay_coeff * seagrass_meadow_width)  # Seagrass wave height attenuation
        if effective_Hs > min_Hs:
            Hs = effective_Hs
        else:  # Enforce minimum wave height
            Hs = min_Hs

    return Hs, Tp


def wavek(F, H):
    """Computes wave number via dispersion relationship. Copyright (C) 2001, Lee Gordon, NortekUSA LLC

    This routine use an approximate equation, then sharpens the result with one interpolation. The result is good to around 1 part in 10^-6.
    The equation came out of a textbook, but I have long since forgotten which one. If you know, please tell me! lgordon@nortekusa.com"""

    g = 9.80171

    e1 = 4 * math.pi ** 2 * F ** 2 * H / g
    e2 = 1 + 0.6666666 * e1 + 0.355555555 * e1 ** 2 + 0.1608465608 * e1 ** 3 + 0.0632098765 * e1 ** 4 + 0.0217540484 * e1 ** 5 + 0.0065407983 * e1 ** 6
    e3 = e1 ** 2 + e1 / e2
    K1 = math.sqrt(e3) / H

    # Compute error as basis for interpolation
    o1 = math.sqrt(g * K1 * math.tanh(K1 * H))
    e1 = o1 ** 2 * H / g
    e2 = 1 + 0.6666666 * e1 + 0.355555555 * e1 ** 2 + 0.1608465608 * e1 ** 3 + 0.0632098765 * e1 ** 4 + 0.0217540484 * e1 ** 5 + 0.0065407983 * e1 ** 6
    e3 = e1 ** 2 + e1 / e2
    K2 = math.sqrt(e3) / H

    # Interpolate
    K = 2 * K1 - K2

    return K


def waveTRNS(amp, wind, fetch, hb, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct):
    """Calculates wave power density at marsh boundary. From Mariotti and Carr (2014)."""

    depth = hb  # Scarp height
    fac = min(1, depth / (2 * amp))  # Proportion of tide that the bay is flooded, dimensionless
    D = (depth + (depth - fac * 2 * amp)) / 2  # [m] average bay depth over tidal cycle
    (Hs, Tp) = YeV(fetch, wind, D, seagrass_on, effective_decay_coeff, seagrass_meadow_width, max_attenuation_pct)  # Solves for wave height (Hs, [m]) and wave period (Tp, [s])
    kk = wavek(1 / Tp, D)  # Solves for wave number [m^-1]
    cg = 2 * math.pi / kk / Tp * 0.5 * (1 + 2 * kk * D / (math.sinh(2 * kk * D)))  # Wave group celerity at marsh edge [m/s]
    W = cg * 9800 / 16 * abs(Hs) ** 2  # [W] Wave power density at the marsh boundary

    return W


def POOLstopp5(t, X):
    B = 10000 * 0.999
    isterminal = [1, 1, 1, 1]  # Stop the integration when value = 0  # IR 15Jun21: commented out because input into Python ODE solver not presently working
    direction = [0, 1, 1, 0]  # 0 = negative direction  # IR 15Jun21: commented out because input into Python ODE solver not presently working
    value = [X[0] - 1, X[0] - B, X[1] - 0.00, X[1] - 0.00]

    return value  # , isterminal, direction
