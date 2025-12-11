import logging
import copy
from MDClasses import EquilibriumRun, SampleRun, StretchRun

logger = logging.getLogger(__name__)


class MDManager:
    def __init__(self, sim_list, atomic_structure):
        """
        Manages the order of MD simulations to be run, with refernce to the settings file.
        The logic behind why the quantities are distributed as they are, among these 4 big dict class variables:
                    If a quantity, to be accurately computed, NEEDS combination(s) of specific ensemble(s) with specific simulation type(s)
        ----------->The quantity then appears in the class variable corresponding to the ensemble, and as an item to the key that
                    corresponds to the simulation type.

                    If, for example, a quantity needs an acceptable conclusion to an equilibration simulation to grab reference values from,
                    but the actual quantity itself is computed in the sample simulation, the quantity will appear as an item in both the key
                    corresponding to the equilibrating simulation and the sample simulation. (In the class variable corrresponding to the correct ensemble ofc)
        """
        self.sim_list = sim_list
        self.atomic_structure = atomic_structure

    def run(self, store_traj=False):
        for simulation, settings in self.sim_list:
            match simulation:
                case "Equilibrium":
                    equilibrium_run = EquilibriumRun(settings=settings)
                    logger.info("Running equilibrium simulation...")
                    self.atomic_structure = equilibrium_run.run(atomic_structure=self.atomic_structure,
                                                                num_steps=settings.num_steps, store_traj=store_traj)
                case "Stretch":
                    stretch_run = StretchRun(settings=settings)
                    logger.info("Running stretching simulation...")
                    stretch_run.run(atomic_structure=self.atomic_structure)
                case "Sample":
                    sample_run = SampleRun(settings=settings)
                    logger.info("Running sample simulation...")
                    sample_run.run(atomic_structure=self.atomic_structure, num_steps=settings.num_steps,
                                   store_traj=store_traj)

        logger.info("MD done")
        logger.info(f"Stored results in {self.atomic_structure.label}/Outputfiles")
