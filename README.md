# PyBMFT-C
Bay-Marsh-Forest Transect with Carbon model

## About
**This is a Python version of the [Coastal Landscape Transect model (*CoLT*) from Valentine et al. (2023)](https://github.com/csdms-contrib/colt).** *CoLT* is originally 
written in Matlab, and expands upon the Bay Marsh Forest Transect model (*BMFT*) from Kirwan et al. (2016) by integrating 
carbon cycling dynamics. See Valentine et al. (2023) for a full description of the model. 

This Python version of *CoLT* - *PyBMFT-C* - is used in the [*BarrierBMFT* coupled model framework](https://github.com/UNC-CECL/BarrierBMFT) (Reeves et al., in review) and 
differs from the original Matlab version of *CoLT* with the following ALTERATIONS/NEW ADDITIONS:

1) *PyBMFT-C* does not contain the option to disperse autochthonous carbon with depth, and therefore only deposits 
autochthonous carbon on the surface of the marsh
2) Sediment supply to the marsh platform is allowed past the first drowned/ponded cell (`evolvemarsh.py`)
3) The organic content of bay sediment is held constant at 5% (`bmftc.py`)
4) Marsh edge progradation and bay accretion are recorded to depositional record, and marsh progradation is included in flux of 
sediment from bay to marsh; eroded sediment is removed from the depositional record (`bmftc.py`)
5) The density of the marsh edge cell is calculated using the **most recent** year where the elevation of the marsh is above the depth
of erosion, rather than the first year (`bmftc.py`)

## Installation

To download the source code for *PyBMFT-C*, either clone the repository with git:

    git clone git@github.com/UNC-CECL/pybmft-c

or download the zip file:

    https://github.com/UNC-CECL/PyBMFT-C/archive/refs/heads/main.zip

Then, run the following from the top-level folder (the one that contains `setup.py`) to install *PyBMFT-C* into the current environment:

    pip install -e .

## Input Files & Parameters

Input files for running a **standalone** *PyBMFT-C* simulation (i.e., not within *BarrierBMFT*) are located in the `PyBMFT-C/Input` directory.
See documentation in the [CoLT repository](https://github.com/csdms-contrib/colt) for additional information on input files.
*PyBMFT-C* parameter values can be set in the initialization of the *PyBMFT-C* classes within a run script (see example below).

## Example Simulation

The following describes the approach for running a **standalone** *PyBMFT-C* simulation. For a more comprehensive example run 
script, see `run_bmftc.py`. 

To run *PyBMFT-C*, first set the working direcory to the main `PyBMFT-C` directory.

Then, import *PyBMFT-C* and dependencies:

    from bmftc import Bmftc  
    import numpy as np
    import matplotlib.pyplot as plt

Next, initialize an instance of the *PyBMFT-C* class. Here, you can set certain parameter values such as the simulation duration (time step count; 
yrs), RSLR rate (mm/yr), or the external suspended sediment supply (reference concentration; mg/L). Parameters values not set here will revert
to their default values:

    bmftc = Bmftc(
                name="PyBMFT-C Example Simulation",
                time_step_count=50,
                relative_sea_level_rise=4,
                reference_concentration=50,
    )
    print(bmftc.name)

Next, run the simulation by progressing through time:

    for time_step in range(int(bmftc.dur)):
    
        # Print time step to screen
        print("\r", "Time Step: ", time_step, end="")
    
        # Run time step
        bmftc.update()

Once the simulation finishes, we can plot some results. For example, the following plots the elevation profile for every 10<sup>th</sup> year:
    
    plt.figure()
    for t in range(bmftc.startyear, bmftc.endyear, 10):
        plt.plot(bmftc.elevation[t, :])
    plt.xlabel("Distance")
    plt.ylabel("Elevation [m MSL]")
    plt.title(bmftc.name)


## References

#### CoLT (Coastal Landscape Transect model)
    Code: https://github.com/csdms-contrib/colt
    Wiki: https://csdms.colorado.edu/wiki/Model:Coastal_Landscape_Transect_Model_(CoLT)
    DOI: 10.5281/zenodo.7625873

    Valentine, K., Herbert, E. R., Walters, D. C., Chen, Y., Smith, A. J., & Kirwan, M. L. (2023). Climate-driven tradeoffs 
    between landscape connectivity and the maintenance of the coastal carbon sink. Nature Communications, 14, 1137. 
    https://doi.org/10.1038/s41467-023-36803-7.
#### BMFT (Bay-Marsh-Forest Transect model)
    Kirwan, M. L., Walters, D. C., Reay, W. G., & Carr, J. A. (2016). Sea level driven marsh expansion in a coupled model 
    of marsh erosion and migration. Geophysical Research Letters, 43, 4366–4373. https://doi.org/10.1002/2016GL068507
#### BarrierBMFT (Barrier Bay-Marsh-Forest Transect model)
    Code: https://github.com/UNC-CECL/BarrierBMFT
    Wiki: TBD
    DOI: TBD

    Reeves, I.R.B., Moore, L.J., Valentine, K., Fagherazzi, S., & Kirwan, M.L. (in review). Sediment exchange across coastal 
    barrier landscapes alters ecosystem extents.



