"""----------------------------------------------------------------------
PyBMFT-C: Bay-Marsh-Forest Transect Carbon Model (Python version)

Last updated _9 June 2021_ by _IRB Reeves_
----------------------------------------------------------------------"""

import scipy.io
import numpy

from buildtransect import buildtransect

# TEMP VARIABLE DEF
C = 10
R = 1
amp = 0.7
bfo = 5000
endyear = 561
mwo = 1000
slope = 0.005
wind = 6

stratfile = "C:/Users/ceclmoore2/PycharmProjects/PyBMFT-C/SpinUp/MarshStrat_all_RSLR1_CO50.mat"
mat = scipy.io.loadmat(stratfile)
elev_25 = mat["elev_25"]

B, dfo, elevation = buildtransect(R, C, slope, mwo, elev_25, amp, wind, bfo, endyear, True)
