import logging
import numpy as np

log = logging.getLogger(__name__) #Nothing will prob log here?




class EquilibriumCondition:

    @staticmethod
    def checkInternalPressureStable( list_of_internal_pressure, tol= 1e-3,):
        """
        Used to check if the difference in  internal pressure is below threshold
        return:
            true if the difference is below threshold
            false if the difference is above threshold
        """
        mean_pressure = np.mean(list_of_internal_pressure)
        delta = np.abs(np.max(list_of_internal_pressure) - np.min(list_of_internal_pressure))
        return (delta <= tol) and (np.isclose(mean_pressure, 0, atol=tol) )


    @staticmethod
    def checkStable(list_of_values, tol = 1e-5):
        """
        Used to check if energy oscillates around certain value, then we probably found the equilibrium,
        meant to input either a list of MSD values or energy, which unlike internal pressure don't have to oscillate
        around 0.

        return:
            true if the difference is below threshold
            false if the difference is above threshold
        """
        print(1)
        delta = np.abs(np.max(list_of_values) - np.min(list_of_values))
        print(delta)
        return delta <= tol


