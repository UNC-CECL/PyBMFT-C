"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _9 June 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import scipy.io
import numpy as np
import math
import matplotlib.pyplot as plt
import bisect

from buildtransect import buildtransect
from evolvemarsh import evolvemarsh

# TEMP VARIABLE DEF buildtransect
C = 10
R = 1
amp = 0.7
bfo = 5000
endyear = 561
mwo = 1000
slope = 0.005
wind = 6

# TEMP VARIABLE DEF evolvemarsh
x_m = 5005
yr = 550
msl = 1 * 10 ** (-3)
C_e = 0.0114
OCb = 0.15
tr = 1.4
numiterations = 500
P = 45000
dt = 90
ws = 5 * 10 ** (-5)
timestep = 700.8
BMax = 2500
Dmin = 0
Dmax = 0.5204
rhoo = 85
rhos = 2000

# Test buildtransect
stratfile = "Input/MarshStrat_all_RSLR1_CO50.mat"
mat = scipy.io.loadmat(stratfile)
elev_25 = mat["elev_25"]

B, dfo, elevation = buildtransect(R, C, slope, mwo, elev_25, amp, wind, bfo, endyear, False)

# Test evolvemarsh
x_f = bisect.bisect_left(elevation[yr - 1], amp)  # First forest cell. Technically amp + dmin + msl[yr]
temp_elevation_marsh = elevation[yr - 1, x_m: x_f + 1]  # values are slightly different here

(marshelevation,
 organic_autoch,
 organic_alloch,
 mineral,
 Fm_min,
 Fm_org,
 bgb,
 accretion,
 agb
 ) = evolvemarsh(temp_elevation_marsh,
                 msl,
                 C_e,
                 OCb,
                 tr,
                 numiterations,
                 P,
                 dt,
                 ws,
                 timestep,
                 BMax,
                 Dmin,
                 Dmax,
                 rhoo,
                 rhos)
