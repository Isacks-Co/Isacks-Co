import logging
import copy
from MDBase import MDBase

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
        q_sims_nve_deps = {"equil": ["cve"], "sample": ["cve"], "stretch1": [], "stretch2": [], "stretch3": [], "any": [] }
        q_sims_nvt_deps = {"equil": ["p_i", "b1", "b2", "b3", "g2", "g3", "cvt", "debye"], "sample": ["cvt"], "stretch1": ["p_i", "b1"], "stretch2": ["b2", "g2", "debye"], "stretch3": ["b3", "g3", "debye"], "any": []}
        q_sims_npt_deps ={"equil": [], "sample": ["lat_const"], "stretch1": [], "stretch2": [], "stretch3": [], "any": []}
        q_sims_indep_of_ensemble = {"equil": ["d", "e_coh", "l_crit", "lat_const"], "sample": ["e_coh"], "stretch1": [], "stretch2": [], "stretch3": [] , "any": ["msd", "l_crit"]}
        self.catergorized_compatibility = [q_sims_nve_deps, q_sims_nvt_deps, q_sims_npt_deps, q_sims_indep_of_ensemble]
        self.order_of_operations = _simulations_to_run(sim_list, self.catergorized_compatibility, quants)


    def run(self, atoms):
        for sim in self.order_of_operations:
            logger.info(f"Running {sim}")
            for ensemble, variations in sim.items():
                match ensemble:
                    case "npt":
                        for var in variations:
                            MD = MDBase(self.sim_list[0], ensemble, var)
                            MD.simulate(atoms)
                    case "nve":
                        for var in variations:
                            MD = MDBase(self.sim_list[1], ensemble, var)
                            MD.simulate(atoms)
                    case "nvt":
                        for var in variations:
                            MD = MDBase(self.sim_list[2], ensemble, var)
                            MD.simulate(atoms)
                    case "indep":
                        for var in variations:
                            MD = MDBase(self.sim_list[2], ensemble, var)
                            MD.simulate(atoms)


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
                if q in deps["stretch1"] :
                    if "stretch1" not in candidate_run_variations:
                        candidate_run_variations.append("stretch1")
                if q in deps["stretch2"] :
                    if "stretch2" not in candidate_run_variations:
                        candidate_run_variations.append("stretch2")
                if q in deps["stretch3"] :
                    if "stretch3" not in candidate_run_variations:
                        candidate_run_variations.append("stretch3")
                """
                if q in deps["any"] :
                    if "any" not in candidate_run_variations:
                        candidate_run_variations.append("any")
                """

            if candidate_run_variations:
                if categorized_compatibility.index(deps) == 0:
                    order_of_operations.append({"nve": candidate_run_variations})
                elif categorized_compatibility.index(deps) == 1:
                    order_of_operations.append({"nvt": candidate_run_variations})
                elif categorized_compatibility.index(deps) == 2:
                    order_of_operations.append({"npt": candidate_run_variations})
                elif categorized_compatibility.index(deps) == 3:
                    order_of_operations.append({"indep": candidate_run_variations})
            candidate_run_variations = []

    logger.debug(f"Order of operations: {order_of_operations}")

    return order_of_operations