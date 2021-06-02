import numpy as np


# from .functions import (
#     decompose,
# )

class Bmftc:
    def __init__(
            self,
            name="default",
            relative_sea_level_rise=0.004,
    ):
        """Bay Marsh Forest Transects with Carbon

        Parameters
        ----------
        name: string, optional
            Name of simulation
        total_number_of_agents: int, optional
            Total number of agents in simulation

        Examples
        --------
        >>> from bmftc import Bmftc
        >>> model = Bmftc()
        """

        self._name = name
        self._rslr = relative_sea_level_rise
        self._time_index = 0

    def update(self):
        """Update Bmftc by a single time step"""

        self._time_index += 1

    @property
    def time_index(self):
        return self._time_index
