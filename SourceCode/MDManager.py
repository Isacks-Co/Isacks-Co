import logging
import copy
from MDClasses import EquilibriumRun, SampleRun, StretchRun

logger = logging.getLogger(__name__)


class MDManager:
    """
    Run a configured sequence of MD simulation stages.

    The manager iterates through a simulation list produced by preprocessing
    and dispatches to the corresponding run class:
    - ``"Equilibrium"`` -> :class:`MDClasses.EquilibriumRun`
    - ``"Stretch"`` -> :class:`MDClasses.StretchRun`
    - ``"Sample"`` -> :class:`MDClasses.SampleRun`

    Parameters
    ----------
    sim_list : iterable
        Iterable of ``(simulation_name, settings)`` pairs describing what to run.
    atomic_structure
        Initial structure object passed into the first stage and updated by
        stages that return a modified structure.
    """
    def __init__(self, sim_list, atomic_structure):
        """
        Initialize the MD simulation manager.

        Parameters
        ----------
        sim_list : iterable
            Iterable of ``(simulation_name, settings)`` pairs. The simulation name
            controls which run stage is executed (e.g. ``"Equilibrium"``).
        atomic_structure
            Structure object to be used as the starting point for the workflow.

        Notes
        -----
        The simulation list is expected to already reflect any dependencies between
        computed quantities and simulation stages (e.g., quantities that require
        equilibration before sampling should be scheduled appropriately by the
        preprocessing step).
        """
        self.sim_list = sim_list
        self.atomic_structure = atomic_structure

    def run(self, store_traj=False):
        """
        Execute the configured MD workflow.

        Iterates through `sim_list` and runs each stage in order. The internal
        `atomic_structure` is updated when stages return an updated structure
        (e.g., equilibration).

        Parameters
        ----------
        store_traj : bool, optional
            If True, enable writing ASE trajectory files for stages that support it.

        Returns
        -------
        None
        """
        for simulation, settings in self.sim_list:
            # Dispatch based on stage name produced by preprocessing
            match simulation:
                case "Equilibrium":
                    # Equilibration updates the structure reference for subsequent stages
                    equilibrium_run = EquilibriumRun(settings=settings)
                    logger.info("Running equilibrium simulation...")
                    self.atomic_structure = equilibrium_run.run(
                        atomic_structure=self.atomic_structure,
                        num_steps=settings.num_steps,
                        store_traj=store_traj,
                        check_conv=True,
                    )

                case "Stretch":
                    # Stretching computes response properties but does not replace the structure
                    stretch_run = StretchRun(settings=settings)
                    logger.info("Running stretching simulation...")
                    stretch_run.run(atomic_structure=self.atomic_structure)

                case "Sample":
                    # Sampling stores per-step quantities; may also write trajectory output
                    sample_run = SampleRun(settings=settings)
                    logger.info("Running sample simulation...")
                    sample_run.run(
                        atomic_structure=self.atomic_structure,
                        num_steps=settings.num_steps,
                        store_traj=store_traj,
                    )

        logger.info("MD done")
        logger.info(f"Stored results in {self.atomic_structure.label}/Outputfiles")
