import numpy as np
import logging
from ase.io.trajectory import Trajectory
from collections import defaultdict

logger = logging.getLogger(__name__)
from Utils.unitConversions import auToGPascal, evToJ
from Utils.plotting import secondOrderNumericalDerivative

def calculateCMatrix(strech_sequence): #TODO Move so it is computed on the fly.
    
    betas = [[], [], [], [], [], []]
    for frame in strech_sequence:
        # Add the corresponding values to the right beta (strain direction)
        betas[frame.info["beta"]].append([frame.info["strain"], frame.info["stress"]])
    beta_arrays = [np.array(beta, dtype=object) for beta in betas]

    # Create dictionaries, one for each beta, and store the arrays of matrices with epsilon as key
    beta_dicts = [defaultdict(list) for i in range(6)]
    for i in range(6):
        for epsilon, matrix in beta_arrays[i]:
            beta_dicts[i][epsilon].append(matrix)

    averages = []
    for beta in beta_dicts:
        avg_data = []

        for epsilon, matrices in beta.items():
            # Take elementwise average over all the matrices for each epsilon
            stacked = np.stack(matrices)
            avg_matrix = stacked.mean(axis=0)
            avg_data.append(np.array((epsilon, avg_matrix), dtype=object))

        avg_data = sorted(avg_data, key=lambda x: x[0])
        averages.append(avg_data)
    C = np.zeros((6, 6))
    for i in range(6):
        # Line fit epsilon vs sigma to find each c_ij
        epsilons = np.array([x[0] for x in averages[i]], dtype=float)
        for j in range(6):
            sigmas = np.array([x[1][j] for x in averages[i]], dtype=float)
            C[j, i] = np.polyfit(epsilons, sigmas, 1)[0]
    C *= 160.21766208  # Convert to GPa
    return C




def _numericalC(stretch_sequence): # TODO Move to do on the fly
    """
    Calculates the elastic constants C11, C22, C33, C12, C44

    Input
        ---

    Output
        C_from_U: matrix, where C11, C12, C44 = C[0,0], C[0,1], C[3,3]
    """
    betas = [[], [], [], [], [], []]
    # Prefer reconstructing 1D slices from 2D trajectory to avoid stale 1D data
    
    used_2d = False
    try:
        
        tol = 1e-12
        for frame in stretch_sequence:
            info = getattr(frame, 'info', {})
            b1 = info.get('beta1')
            b2 = info.get('beta2')
            if b1 is None or b2 is None:
                continue
            if int(b1) == int(b2) and np.isclose(float(info.get('strain1', 0.0)), 0.0, atol=tol):                   #Looks at difference between strains when creating 2D data
                i = int(b1)
                eps = float(info.get('strain2', 0.0))
                energy = float(info.get('total_energy', np.nan))
                stress = np.array(info.get('stress'), dtype=float)
                betas[i].append([eps, stress, energy])
        used_2d = any(len(b) > 0 for b in betas)
    except Exception as e2:
        logger.info(f"Failed")
    if not used_2d:
        # Fallback to legacy 1D if available
        try:
            stretch_trajectory = Trajectory(self.settings.output_file + "_stretch_data.traj")
            for frame in stretch_trajectory:
                energy = frame.info.get('total_energy')
                betas[int(frame.info['beta'])].append([float(frame.info['strain']), np.array(frame.info['stress'],
                                                                                    dtype=float), float(energy)])
            logger.info("Used legacy 1D stretch trajectory to compute elastic constants.")
        except Exception as e:
            logger.info(f"No usable stretch data found (2D or 1D): {e}")
    beta_arrays = [np.array(beta, dtype=object) for beta in betas]
    beta_dicts = [defaultdict(list) for i in range(6)]
    averages = []

    for i in range(6):
        for eps, matrix, energy in beta_arrays[i]:
            beta_dicts[i][eps].append((matrix, energy))

    for beta in beta_dicts:
        avg_data = []

        for eps, matrices_and_energies in beta.items():
            # Take elementwise average for each epsilon
            matrices = np.stack([me[0] for me in matrices_and_energies])
            energies = np.array([me[1] for me in matrices_and_energies], dtype=float)
            avg_matrix = matrices.mean(axis=0)
            avg_energy = np.nanmean(energies)  # ignore possible NaNs
            avg_data.append(np.array((eps, avg_matrix, avg_energy), dtype=object))

        avg_data = sorted(avg_data, key=lambda x: x[0])
        averages.append(avg_data)

    second_deriv = np.zeros((6,6))
    

    twoD_energies = stretch_sequence[-1].info["2D Energies"]
    strains_axis = stretch_sequence[-1].info["Strains axis"]
    number_of_pairs = stretch_sequence[-1].info["Number of pairs"]         # Should usually be 6

    for i in range(int(np.sqrt(number_of_pairs))):
        energy_1 = np.array([x[2] for x in averages[i]], dtype=float)
        for j in range(int(np.sqrt(number_of_pairs))):
            energy_2 = np.array([x[2] for x in averages[j]], dtype=float)
            if i == j:
                second_deriv[i, j] = secondOrderNumericalDerivative(strains_axis, [energy_1, energy_2])
                continue
            else:
                try:
                    second_deriv[i, j] = secondOrderNumericalDerivative(strains_axis, twoD_energies[i][j])
                except Exception as e:
                    logger.info(f"2D stretch calc failed for ({i},{j}), WHY THE FRICK???!: {e}")
                    second_deriv[i, j] = 0.0

    C_from_U = second_deriv / stretch_sequence[0].get_volume()
    logger.debug(f"C_from_U = \n {C_from_U * auToGPascal(1)} \n")
    logger.info(f" \n C_11 = {auToGPascal(C_from_U[0,0])} \n C_12 = {auToGPascal(C_from_U[0,1])} \n C_44 = {auToGPascal(C_from_U[3,3])}")

    return C_from_U









