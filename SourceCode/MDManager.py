import logging
import copy
from MDClasses import EquilibriumRun, SampleRun, StrecthRun

logger = logging.getLogger(__name__)

class MDManager:
    def __init__(self, sim_list, quants):
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
        self.quants = quants
        q_sims_nvt_deps = {"equil": ["p_i", "b", "g", "cvt", "debye"], "sample": ["cvt"], "stretch": ["p_i", "b"], "any": []}
        q_sims_npt_deps ={"equil": [], "sample": ["lat_const"], "stretch": [], "any": []}
        q_sims_indep_of_ensemble = {"equil": ["d", "e_coh", "l_crit", "lat_const"], "sample": ["e_coh"], "stretch": [] , "any": ["msd", "l_crit"]}
        self.catergorized_compatibility = [q_sims_nvt_deps, q_sims_npt_deps, q_sims_indep_of_ensemble]
        self.order_of_operations = _simulations_to_run(sim_list, self.catergorized_compatibility, quants)


    def run(self, atomic_structure):
        for sim in self.order_of_operations:
            logger.debug(f"Running {sim}")
            for ensemble, variations in sim.items():
                logger.debug(f"Variations: {variations}")
                for var in variations:
                    self.simulate(ensemble=ensemble, variation=var, atomic_structure=atomic_structure)

    def simulate(self, ensemble, variation, atomic_structure):
        if variation == "equil":
            equil_settings = self.sim_list[0]
            equil_MD = EquilibriumRun(settings=equil_settings)
            logger.info("Relaxing structure")
            equil_struct = equil_MD.run(atomic_structure, equil_settings.num_steps, init_vel=True)

        if variation == "sample" or variation == "any":
            if ensemble == "nvt":
                sample_settings = self.sim_list[1]
            else:
                sample_settings = self.sim_list[2]
            sample_MD = SampleRun(settings=self.sim_list[1])
            logger.info("Sampling structure")
            sample_data = sample_MD.run(equil_struct, sample_settings.num_steps)

        if variation == "stretch":
            stretch_MD = StrecthRun(settings=self.sim_list[1])
            logger.info("Running stretch sequence")
            sample_data.storeTxtFile()
            C_matrix = stretch_MD.run(equil_struct)

        logger.info("MD done")
        logger.info(f"Stored results in {equil_struct.label}/Outputfiles")


def _simulations_to_run(sim_list, cat_resp, quants = None):
    order_of_operations = []
    candidate_run_variations = []
    categorized_compatibility = cat_resp

    if sim_list[0]:
        order_of_operations.append({"npt": ["equil"]})

    if quants:
        for deps in categorized_compatibility:
            for q in quants:
                q = q.lower()
                logger.debug(f"q: {q}")
                """
                if q in deps["equil"]:
                    if "equil" not in candidate_run_variations:
                        candidate_run_variations.append("equil")
                """
                if q in deps["sample"]:
                    if "sample" not in candidate_run_variations:
                        candidate_run_variations.append("sample")
                if q in deps["stretch"] :
                    if "stretch" not in candidate_run_variations:
                        candidate_run_variations.append("stretch")
                """
                if q in deps["any"] :
                    if "any" not in candidate_run_variations:
                        candidate_run_variations.append("any")
                """

            if candidate_run_variations:
                if categorized_compatibility.index(deps) == 0:
                    order_of_operations.append({"nvt": candidate_run_variations})
                elif categorized_compatibility.index(deps) == 1:
                    order_of_operations.append({"npt": candidate_run_variations})
                elif categorized_compatibility.index(deps) == 2:
                    order_of_operations.append({"indep": candidate_run_variations})
            candidate_run_variations = []

    logger.debug(f"Order of operations: {order_of_operations}")

    return order_of_operations