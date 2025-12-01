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
        self.quants_by_variation = {"equil": [], "sample": [], "stretch": [], "any":[]}
        q_sims_nvt_deps = {"equil": ["p_i", "b", "g", "cvt", "debye"], "sample": ["cvt", "t", "e_tot", "e_kin", "e_pot", "v"], "stretch": ["p_i", "b"], "any": []}
        q_sims_npt_deps ={"equil": [], "sample": ["lat_const", "t", "e_tot", "e_kin", "e_pot", "v"], "stretch": [], "any": []}
        q_sims_indep_of_ensemble = {"equil": ["d", "e_coh", "l_crit", "lat_const"], "sample": ["e_coh", "t", "e_tot", "e_kin", "e_pot", "v"], "stretch": [] , "any": ["msd", "l_crit"]}
        self.catergorized_compatibility = [q_sims_nvt_deps, q_sims_npt_deps, q_sims_indep_of_ensemble]
        self.order_of_operations = self._simulations_to_run(sim_list, self.catergorized_compatibility, quants)
        self.equil_struct = None


    def run(self, atomic_structure, init_vel=False, store_traj=False):
        for sim in self.order_of_operations:
            logger.debug(f"Running {sim}")
            for ensemble, variations in sim.items():
                logger.debug(f"Variations: {variations}")
                for var in variations:
                    if var != "any":
                        self.simulate(ensemble=ensemble, variation=var, atomic_structure=atomic_structure, init_vel=init_vel, store_traj=store_traj)

    def simulate(self, ensemble, variation, atomic_structure, init_vel, store_traj):
        if variation == "equil":
            equil_settings = self.sim_list[0]
            equil_MD = EquilibriumRun(settings=equil_settings)
            logger.info("Relaxing structure")
            self.equil_struct = equil_MD.run(atomic_structure=atomic_structure, num_steps=equil_settings.num_steps, init_vel=init_vel, store_traj=store_traj)

        elif variation == "sample":
            if self.quants_by_variation["any"]:
                sample_quants = self.quants_by_variation["sample"] + self.quants_by_variation["any"]
            else:
                sample_quants = self.quants_by_variation["sample"]
            if ensemble == "nvt" or ensemble == "indep":
                sample_settings = self.sim_list[1]
            else:
                sample_settings = self.sim_list[2]
            sample_MD = SampleRun(settings=sample_settings, sample_data=sample_quants)
            logger.info("Sampling structure")
            if self.sim_list[0]:
                sample_data = sample_MD.run(atomic_structure=self.equil_struct, num_steps=sample_settings.num_steps, store_traj=store_traj)
            else:
                sample_data = sample_MD.run(atomic_structure=atomic_structure, num_steps=sample_settings.num_steps, store_traj=store_traj)
            sample_data.storeTxtFile()

        elif variation == "stretch":
            stretch_settings = self.sim_list[1]
            stretch_MD = StrecthRun(settings=stretch_settings)
            logger.info("Running stretch sequence")
            if self.sim_list[0]:
                C_matrix = stretch_MD.run(atomic_structure=self.equil_struct)
            else:
                C_matrix = stretch_MD.run(atomic_structure=atomic_structure)

        logger.info("MD done")
        logger.info(f"Stored results in {equil_struct.label}/Outputfiles")


    def _simulations_to_run(self, sim_list, cat_resp, quants = None):
        order_of_operations = []
        candidate_run_variations = []
        categorized_compatibility = cat_resp

        if sim_list[0]:
            order_of_operations.append({"npt": ["equil"]})

        if quants:
            for deps in categorized_compatibility:
                quants_by_deps = []
                for q in quants:
                    q = q.lower()
                    if q in deps["any"] :
                        if q not in quants_by_deps:
                            quants_by_deps.append(q)
                            self.quants_by_variation["any"].append(q)
                        if "any" not in candidate_run_variations:
                            candidate_run_variations.append("any")

                    if q in deps["sample"]:
                        if q not in quants_by_deps:
                            self.quants_by_variation["sample"].append(q)
                            quants_by_deps.append(q)
                        if "sample" not in candidate_run_variations:
                            candidate_run_variations.append("sample")

                    if q in deps["stretch"]:
                        if q not in quants_by_deps:
                            self.quants_by_variation["stretch"].append(q)
                            quants_by_deps.append(q)
                        if "stretch" not in candidate_run_variations:
                            candidate_run_variations.append("stretch")

                if candidate_run_variations:
                    if categorized_compatibility.index(deps) == 0:
                        order_of_operations.append({"nvt": candidate_run_variations})
                    elif categorized_compatibility.index(deps) == 1:
                        order_of_operations.append({"npt": candidate_run_variations})
                    elif categorized_compatibility.index(deps) == 2:
                        order_of_operations.append({"indep": candidate_run_variations})
                candidate_run_variations = []

        return order_of_operations