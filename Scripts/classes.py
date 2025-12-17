#!/usr/bin/env python
#
# Description:
#    This file was created as part of the Automatic Defect Analysis and Qualification (ADAQ) software packages.
#    For additional information about ADAQ, visit: https://defects.anyterial.se/
#
# License:
#    This file is part of the project distributed under the MIT License.
#    Copyright (c) 2025 Joel Davidsson
#
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

from __future__ import print_function, division

try:
    import httk
except Exception:
    import sys, os.path, inspect
    _realpath = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
    sys.path.insert(1, os.path.join(_realpath, '../..'))
    import httk

from httk.atomistic import Structure, StructureTag
from httk.core.vectors import FracVector

# This class decide if you want a cubic structure or copies of your structure
class SupercellShape(object):
    # cubic, tolerance
    # copy, [x,y,z]

    def __init__(self, desc, data, dim):
        self.desc = desc
        self.data = data
        self.dim = dim

    # describes the shape copy or cubic
    def desc(self):
        return self.desc

    # data can be how many copies in x,y,z or tolerance for cubic
    def data(self):
        return self.data

    # dimension of the material 2D or 3D
    def dim(self):
        return self.dim

# This class stores the version of adaq
class ADAQ(httk.HttkObject):

    @httk.httk_typed_init({'git_commit': str})
    def __init__(self, git_commit):
        self.git_commit = git_commit

# This class for the host material properties
class HostUnitCellResult(httk.Result):

    @httk.httk_typed_init({'host_name': str, 'computation': httk.Computation,  'structure': Structure, 'formation_entalpy': float, 'mp_id': str}, index=['host_name'])
    def __init__(self, host_name, computation, structure, mp_id, formation_entalpy):
        self.host_name = host_name
        self.computation = computation
        self.structure = structure
        self.mp_id = mp_id
        self.formation_entalpy = formation_entalpy # eV/atom

# Class to store hand hanle the dielectric tensor
class DielectricTensor(httk.Result):
    @httk.httk_typed_init({'host_name': str, 'effective_epsilon': float, 'madelung_potential': float, 'dielectric_tensor_electric': [float], 'dielectric_tensor_ionic': [float], 'dielectric_tensor_total_supercell': [float], 'refractive_index': float}, index=['host_name'])
    def __init__(self, host_name, effective_epsilon, madelung_potential, dielectric_tensor_electric, dielectric_tensor_ionic, dielectric_tensor_total_supercell, refractive_index):
        self.host_name = host_name
        self.effective_epsilon = effective_epsilon # dimensionless
        self.madelung_potential = madelung_potential
        self.dielectric_tensor_electric = dielectric_tensor_electric #dimensionless, unitcell cartesian coordinates
        self.dielectric_tensor_ionic = dielectric_tensor_ionic #dimensionless, unitcell cartesian coordinates
        self.dielectric_tensor_total_supercell = dielectric_tensor_total_supercell #dimensionless, supercell cartesian coordinates
        self.refractive_index = refractive_index #dimensionless


# This class for the Host material
class HostSuperCell(httk.HttkObject):

    @httk.httk_typed_init({'unit_cell': Structure, 'host_supercell': Structure, 'description': str, 'data': [int], 'dimension': str, 'material': str, 'functional': str})
    def __init__(self, unit_cell, host_supercell, description, data, dimension, material, functional):
        self.unit_cell = unit_cell
        self.host_supercell = host_supercell
        self.description = description
        self.data = data
        self.dimension = dimension
        # which material is it
        self.material = material
        # which functional was used to determine the lattice parameters
        self.functional = functional

    def name(self):
        return self.material+'_'+self.functional

class HostSuperCellResult(httk.Result):

    @httk.httk_typed_init({'host_name': str, 'computation': httk.Computation, 'screen_energy': float, 'screen_valance_band_max': float, 'screen_conduction_band_min': float, 'total_energy_coarse': float, 'coarse_valance_band_max': float, 'coarse_conduction_band_min': float, 'total_energy': float, 'full_valance_band_max': float, 'full_conduction_band_min': float, 'madelung_potential': float}, index=['host_name'])
    def __init__(self, host_name, computation, screen_energy, screen_valance_band_max, screen_conduction_band_min, total_energy_coarse, coarse_valance_band_max, coarse_conduction_band_min, total_energy, full_valance_band_max, full_conduction_band_min, madelung_potential):
        self.host_name = host_name
        self.computation = computation
        # screen
        self.screen_energy = screen_energy
        self.screen_valance_band_max = screen_valance_band_max
        self.screen_conduction_band_min = screen_conduction_band_min
        # full
        self.total_energy_coarse = total_energy_coarse
        self.coarse_valance_band_max = coarse_valance_band_max
        self.coarse_conduction_band_min = coarse_conduction_band_min
        self.total_energy = total_energy
        self.full_valance_band_max = full_valance_band_max
        self.full_conduction_band_min = full_conduction_band_min

        # charge correction
        self.madelung_potential = madelung_potential

# These two classes contain the info about the defect:
# defectinfo -  contain overview info, fast class to search
# defectcell - contain the structure
class DefectInfo(httk.HttkObject):

    @httk.httk_typed_init({'key': int, 'host_name': str, 'defect_name': str, 'defect_stoichiometry': str, 'defect_type': str, 'configuration': str, 'defect_size': int, 'vacancy': bool, 'substitutional': bool, 'interstitial': bool, 'S': bool, 'P': bool, 'D': bool, 'F': bool}, index=['key', 'defect_name', 'defect_stoichiometry', 'defect_type', 'defect_size', 'vacancy', 'substitutional', 'interstitial', 'S', 'P', 'D', 'F'])
    def __init__(self, key, host_name, defect_name, defect_stoichiometry, defect_type, configuration, defect_size, vacancy, substitutional, interstitial, S, P, D, F):
        self.key = key
        self.host_name = host_name
        # unique name for the defect
        self.defect_name = defect_name
        # what stoichiometry the defect has
        self.defect_stoichiometry = defect_stoichiometry
        # and where in the supercell it is localated
        self.defect_type = defect_type
        # unique hash for the defect
        self.configuration = configuration
        self.defect_size = defect_size
        # defect info
        self.vacancy = vacancy
        self.substitutional = substitutional
        self.interstitial = interstitial
        self.S = S
        self.P = P
        self.D = D
        self.F = F

class DefectCell(httk.HttkObject):

    @httk.httk_typed_init({'host_name': str, 'defect_structure': Structure, 'defect_name': str, 'defect_types': [str], 'defect_positions': (FracVector, 0, 3), 'key': int, 'priority': int, 'description': str}, index=['defect_name', 'key', 'priority'])
    def __init__(self, host_name, defect_structure, defect_name, defect_types, defect_positions, key, priority, description):
        self.host_name = host_name
        self.defect_structure = defect_structure
        # unique name for the defect
        self.defect_name = defect_name
        # what defect is it
        self.defect_types = defect_types
        # and where in the supercell it is localated
        self.defect_positions = defect_positions
        # unique hash for the defect
        self.key = key
        # priority for the defect
        self.priority = priority
        # general description of the defect
        self.description = description


# This class is for the band data
# Only ipr data for ground state, the excited ipr are the same as the ground state
class BandData(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'spin': float, 'defect_folder_name':str, 'workflow': str, 'state': str, 'name': str, 'bands': [int], 'eigenvalue_spin1': [float], 'occupation_spin1': [float], 'ipr_spin1': [float], 'eigenvalue_spin2': [float], 'occupation_spin2': [float], 'ipr_spin2': [float]}, index=['defect_key', 'workflow', 'state'])
    def __init__(self, defect_key, computation, charge, spin, defect_folder_name, workflow, state, name, bands, eigenvalue_spin1, occupation_spin1, ipr_spin1, eigenvalue_spin2, occupation_spin2, ipr_spin2):
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        # data
        self.bands = bands
        # spin 1
        self.eigenvalue_spin1 = eigenvalue_spin1
        self.occupation_spin1 = occupation_spin1
        self.ipr_spin1 = ipr_spin1
        # spin 2
        self.eigenvalue_spin2 = eigenvalue_spin2
        self.occupation_spin2 = occupation_spin2
        self.ipr_spin2 = ipr_spin2


# This class is for the adaq band data
# it stores the states used by ADAQ and any lone states.
class ADAQBandData(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'spin': float, 'defect_folder_name':str, 'workflow': str, 'state': str, 'name': str, 'ADAQ_state_spin1': [int], 'ADAQ_state_spin2': [int], 'lone_states_spin1': [int], 'lone_states_spin2': [int]}, index=['defect_key', 'workflow', 'state'])
    def __init__(self, defect_key, computation, charge, spin, defect_folder_name, workflow, state, name, ADAQ_state_spin1, ADAQ_state_spin2, lone_states_spin1, lone_states_spin2):
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        # spin 1
        self.ADAQ_state_spin1 = ADAQ_state_spin1
        self.lone_states_spin1 = lone_states_spin1
        # spin 2
        self.ADAQ_state_spin2 = ADAQ_state_spin2
        self.lone_states_spin2 = lone_states_spin2

# This class has all the charge transition levels
class ChargeTransitionLevels(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'description': [str], 'from_charge': [int], 'to_charge': [int], 'from_spin': [int], 'to_spin': [int], 'level': [float]}, index=['defect_key'])
    def __init__(self, defect_key, description, from_charge, to_charge, from_spin, to_spin, level):
        self.defect_key = defect_key
        # screen or full
        self.description = description
        # transition levels
        self.from_charge = from_charge
        self.to_charge = to_charge
        # spins
        self.from_spin = from_spin
        self.to_spin = to_spin
        # relative to host VBM
        self.level = level

# These two classes contain the screening results for a defect:
# ScreenResult - just the results from the screening workflow, fast class to search
# ScreenCell - structures and relaxation from starting geometry
class ScreenResult(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'spin': float, 'defect_folder_name': str, 'total_energy_coarse': float, 'smallest_abs_estimate': float, 'smallest_abs_transition_intensity': float, 'smallest_abs_partial_density_difference': float, 'smallest_abs_polarization_x': float, 'smallest_abs_polarization_y': float, 'smallest_abs_polarization_z': float, 'excitation': str, 'abs_estimate': float, 'abs_transition_intensity': float, 'abs_partial_density_difference': float, 'abs_polarization_x': float, 'abs_polarization_y': float, 'abs_polarization_z': float, 'abs_radiative_lifetime': float, 'ZPL_converged': bool, 'ZPL_estimate': float, 'ZPL_transition_intensity': float, 'ZPL_partial_density_difference': float, 'ZPL_polarization_x': float, 'ZPL_polarization_y': float, 'ZPL_polarization_z': float, 'ZPL_radiative_lifetime': float, 'max_relaxation': float, 'average_relaxation': float, 'delta_R': float, 'delta_Q': float, 'huang_rhys': float, 'debye_waller': float, 'omega': float, 'time': float}, index=['defect_key'])
    def __init__(self, defect_key, computation, charge, spin, defect_folder_name, total_energy_coarse, smallest_abs_estimate, smallest_abs_transition_intensity, smallest_abs_partial_density_difference, smallest_abs_polarization_x, smallest_abs_polarization_y, smallest_abs_polarization_z, excitation, abs_estimate, abs_transition_intensity, abs_partial_density_difference, abs_polarization_x, abs_polarization_y, abs_polarization_z, abs_radiative_lifetime, ZPL_converged, ZPL_estimate, ZPL_transition_intensity, ZPL_partial_density_difference, ZPL_polarization_x, ZPL_polarization_y, ZPL_polarization_z, ZPL_radiative_lifetime, max_relaxation, average_relaxation, delta_R, delta_Q, huang_rhys, debye_waller, omega, time):
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.total_energy_coarse = total_energy_coarse
        # smallest optical estimates, from ground state only
        self.smallest_abs_estimate = smallest_abs_estimate
        self.smallest_abs_transition_intensity = smallest_abs_transition_intensity
        self.smallest_abs_partial_density_difference = smallest_abs_partial_density_difference
        self.smallest_abs_polarization_x = smallest_abs_polarization_x
        self.smallest_abs_polarization_y = smallest_abs_polarization_y
        self.smallest_abs_polarization_z = smallest_abs_polarization_z
        # optical estimates from ground and excited state
        self.excitation = excitation
        # abs
        self.abs_estimate = abs_estimate
        self.abs_transition_intensity = abs_transition_intensity
        self.abs_partial_density_difference = abs_partial_density_difference
        self.abs_polarization_x = abs_polarization_x
        self.abs_polarization_y = abs_polarization_y
        self.abs_polarization_z = abs_polarization_z
        self.abs_radiative_lifetime = abs_radiative_lifetime
        # zpl
        self.ZPL_converged = ZPL_converged
        self.ZPL_estimate = ZPL_estimate
        self.ZPL_transition_intensity = ZPL_transition_intensity
        self.ZPL_partial_density_difference = ZPL_partial_density_difference
        self.ZPL_polarization_x = ZPL_polarization_x
        self.ZPL_polarization_y = ZPL_polarization_y
        self.ZPL_polarization_z = ZPL_polarization_z
        self.ZPL_radiative_lifetime = ZPL_radiative_lifetime
        # relaxation from ground to excited geometry
        self.max_relaxation = max_relaxation
        self.average_relaxation = average_relaxation
        self.delta_R = delta_R
        self.delta_Q = delta_Q
        # one-phonon approximation
        self.huang_rhys = huang_rhys
        self.debye_waller = debye_waller
        self.omega = omega # E(Q) = 1/2 Omega^2 Q^2, saved in eV
        # calculation time
        self.time = time

class ScreenCell(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'spin': float, 'defect_folder_name': str, 'workflow': str, 'state': str, 'name': str , 'structure': Structure}, index=['defect_key', 'workflow', 'state'])
    def __init__(self, defect_key, computation, charge, spin, defect_folder_name, workflow, state, name, structure):
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        self.structure = structure

# this class contain info from the full characterization workflow
# similiar to ScreenResult
# contain all info for webpage visualaizion, check! all converged ZPL, like transition script
# add lifetime here and to screen.
class FullResult(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'spin': float, 'total_energy_coarse': float, 'smallest_ZPL_estimate': float, 'smallest_transition_intensity': float, 'smallest_partial_density_difference': float, 'excitation': str, 'converged': bool, 'ZPL_estimate': float, 'transition_intensity': float, 'partial_density_difference': float, 'max_relaxation': float, 'average_relaxation': float, 'time': float}, index=['defect_key'])
    def __init__(self, defect_key, computation, charge, spin, total_energy_coarse, smallest_ZPL_estimate, smallest_transition_intensity, smallest_partial_density_difference, excitation, converged, ZPL_estimate, transition_intensity, partial_density_difference, max_relaxation, average_relaxation, time):
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.spin = spin
        self.total_energy_coarse = total_energy_coarse
        # smallest optical estimates, from ground state only
        self.smallest_ZPL_estimate = smallest_ZPL_estimate
        self.smallest_transition_intensity = smallest_transition_intensity
        self.smallest_partial_density_difference = smallest_partial_density_difference
        # optical estimates from ground and excited state
        self.excitation = excitation
        self.converged = converged
        self.ZPL_estimate = ZPL_estimate
        self.transition_intensity = transition_intensity
        self.partial_density_difference = partial_density_difference
        # relaxation from start to present geometry
        self.max_relaxation = max_relaxation
        self.average_relaxation = average_relaxation
        # calculation time
        self.time = time

# This class is for the spin density
class SpinDensity(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'charge': int, 'spin': float, 'defect_folder_name': str, 'workflow': str, 'state': str, 'name': str , 'center': [float], 'ipr': float}, index=['defect_key'])
    def __init__(self, defect_key, charge, spin, defect_folder_name, workflow, state, name, center, ipr):
        # info
        self.defect_key = defect_key
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        # center of spin density, stored in fractional coordinates
        self.center = center
        # The localization of spin density
        self.ipr = ipr

# This class is for the D-tensor
class ZeroFieldSplitting(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'charge': int, 'spin': float, 'defect_folder_name': str, 'workflow': str, 'state': str, 'name': str , 'tensor': [float], 'eigenvalues': [float], 'eigenvectors': [float], 'D': float, 'E': float}, index=['defect_key'])
    def __init__(self, defect_key, charge, spin, defect_folder_name, workflow, state, name, tensor, eigenvalues, eigenvectors, D , E):
        # info
        self.defect_key = defect_key
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        # tensor in Mhz order D_xx, D_yy, D_zz, D_xy, D_xz, D_yz, as calculated
        self.tensor = tensor
        # eigenvalues in Mhz (e1, e2, e3), diagonalized values
        self.eigenvalues = eigenvalues
        # eigenvectors (e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z), corresponding eigenvectors
        self.eigenvectors = eigenvectors
        # D=D_zz-(D_yy+D_xx)/2.0) [MHz]
        self.D = D
        # E=(D_yy-D_xx)/2.0 [MHz]
        self.E = E

# This class is for the band data for ZFS
# it stores the states used by ADAQ to calculate the D tensor.
class DBandData(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'spin': float, 'defect_folder_name':str, 'workflow': str, 'state': str, 'name': str, 'ADAQ_D_state_spin1': [int], 'ADAQ_D_state_spin2': [int]}, index=['defect_key', 'workflow', 'state'])
    def __init__(self, defect_key, computation, charge, spin, defect_folder_name, workflow, state, name, ADAQ_D_state_spin1, ADAQ_D_state_spin2):
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        # spin 1
        self.ADAQ_D_state_spin1 = ADAQ_D_state_spin1
        # spin 2
        self.ADAQ_D_state_spin2 = ADAQ_D_state_spin2

# This class is for hyperfine tensor
# update properties in this class
class HyperFineTensor(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'charge': int, 'spin': float, 'defect_folder_name': str, 'workflow': str, 'state': str, 'name': str , 'fermi_contact_point': [float], 'A_xx': [float], 'A_yy': [float], 'A_zz': [float], 'A_xy': [float], 'A_xz': [float], 'A_yz': [float]}, index=['defect_key'])
    def __init__(self, defect_key, charge, spin, defect_folder_name, workflow, state, name, fermi_contact_point, A_xx, A_yy, A_zz, A_xy, A_xz, A_yz):
        # info
        self.defect_key = defect_key
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        # Fermi contact point in MHz, needs to be multiplied with the gyromag in Mhz T^-1 for selected isotope
        self.fermi_contact_point = fermi_contact_point
        # A tensor in Mhz, needs to be multiplied with the gyromag in Mhz T^-1 for selected isotope
        self.A_xx = A_xx
        self.A_yy = A_yy
        self.A_zz = A_zz
        self.A_xy = A_xy
        self.A_xz = A_xz
        self.A_yz = A_yz

# This class is for Quadruple tensor
# update properties in this class
class QuadrupleTensor(httk.HttkObject):

    @httk.httk_typed_init({'defect_key': int, 'charge': int, 'spin': float, 'defect_folder_name': str, 'workflow': str, 'state': str, 'name': str , 'Q_xx': [float], 'Q_yy': [float], 'Q_zz': [float], 'Q_xy': [float], 'Q_xz': [float], 'Q_yz': [float]}, index=['defect_key'])
    def __init__(self, defect_key, charge, spin, defect_folder_name, workflow, state, name, Q_xx, Q_yy, Q_zz, Q_xy, Q_xz, Q_yz):
        # info
        self.defect_key = defect_key
        self.charge = charge
        self.spin = spin
        self.defect_folder_name = defect_folder_name
        self.workflow = workflow
        self.state = state
        self.name = name
        # Q tensor in Mhz, needs to be multiplied with the quadruple moment for selected isotope
        self.Q_xx = Q_xx
        self.Q_yy = Q_yy
        self.Q_zz = Q_zz
        self.Q_xy = Q_xy
        self.Q_xz = Q_xz
        self.Q_yz = Q_yz

# this class is for a ground state of a defect run
class GroundStateResult(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'structure': Structure, 'charge': int, 'spin': float, 'total_energy_coarse': float, 'total_energy': float}, index=['defect_key'])
    def __init__(self, defect_key, computation, structure, charge, spin, total_energy_coarse, total_energy):
        self.defect_key = defect_key
        self.computation = computation
        self.structure = structure
        self.charge = charge
        self.spin = spin
        self.total_energy_coarse = total_energy_coarse
        self.total_energy = total_energy


# this class is for a excited state of a defect run
# update
class ExcitedStateResult(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'name': str, 'charge': int, 'spin': float, 'excitation': str, 'absorption': bool, 'total_energy_absorption': float, 'excited_coarse': bool, 'total_energy_coarse': float, 'excited': bool, 'total_energy': float, 'emission': bool, 'total_energy_emission': float}, index=['defect_key', 'excited'])
    def __init__(self, defect_key, computation, name, charge, spin, absorption, total_energy_absorption, excited_coarse, total_energy_coarse, excited, total_energy, emission, total_energy_emission):
        self.defect_key = defect_key
        self.computation = computation
        self.name = name
        self.charge = charge
        self.spin = spin
        self.absorption = absorption
        self.total_energy_absorption = total_energy_absorption
        self.excited_coarse = excited_coarse
        self.total_energy_coarse = total_energy_coarse
        self.excited = excited
        self.total_energy = total_energy
        self.emission = emission
        self.total_energy_emission = total_energy_emission


# this class is for a ZPL of a defect run
# update for all transitions, clearer refernce to other data
# look over indexe
class TransitionResult(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'defect_spin': float, 'transition': str, 'step': str, 'name': str, 'spin_channel': int, 'from_band': int, 'to_band': int, 'transition_energy': float, 'polarization_x': float, 'polarization_y': float, 'polarization_z': float, 'transition_intensity': float, 'delta_Q': float}, index=['defect_key'])
    def __init__(self, defect_key, computation, charge, defect_spin, transition, step, name, spin_channel, from_band, to_band, transition_energy, polarization_x, polarization_y, polarization_z, transition_intensity, delta_Q):
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.defect_spin = defect_spin
        self.transition = transition
        self.step = step
        self.name = name
        self.spin_channel = spin_channel
        self.from_band = from_band
        self.to_band = to_band
        self.transition_energy = transition_energy
        self.polarization_x = polarization_x
        self.polarization_y = polarization_y
        self.polarization_z = polarization_z
        self.transition_intensity = transition_intensity
        self.delta_Q = delta_Q

# this class is for the chemical potential for each material
class ChemicalPotential(httk.Result):

    @httk.httk_typed_init({'material': str, 'computation': httk.Computation, 'chemical_potential_coarse': float, 'chemical_potential': float, 'functional': str}, index=['material'])
    def __init__(self, material, computation, chemical_potential_coarse, chemical_potential, functional):
        self.material = material
        self.computation = computation
        self.chemical_potential_coarse = chemical_potential_coarse
        self.chemical_potential = chemical_potential
        self.functional = functional

# this class stores the only defects on the defect hull
class DefectHull(httk.HttkObject):

    @httk.httk_typed_init({'host_name': str, 'defect_stoichiometry': str, 'defect_key': int, 'defect_charge': int, 'defect_spin': float, 'fermi_levels': [float], 'formation_energies': [float]}, index=['host_name', 'defect_stoichiometry', 'defect_key'])
    def __init__(self, host_name, defect_stoichiometry, defect_key, defect_charge, defect_spin, fermi_levels, formation_energies):
        # host material and defect stoichiometry
        self.host_name = host_name
        self.defect_stoichiometry = defect_stoichiometry
        # defect data
        self.defect_key = defect_key
        self.defect_charge = defect_charge
        self.defect_spin = defect_spin
        # formation energy data
        self.fermi_levels = fermi_levels
        self.formation_energies = formation_energies

# this class stores the the distance to the defect hull for all defects
class HullDistance(httk.HttkObject):

    @httk.httk_typed_init({'host_name': str, 'defect_stoichiometry': str, 'defect_key': int, 'defect_charge': int, 'defect_spin': float, 'formation_energy_const': float, 'min_distance': float, 'interval': [float], 'span': float}, index=['host_name', 'defect_stoichiometry', 'defect_key'])
    def __init__(self, host_name, defect_stoichiometry, defect_key, defect_charge, defect_spin, formation_energy_const, min_distance, interval, span):
        # host material and defect stoichiometry
        self.host_name = host_name
        self.defect_stoichiometry = defect_stoichiometry
        # defect data
        self.defect_key = defect_key
        self.defect_charge = defect_charge
        self.defect_spin = defect_spin
        self.formation_energy_const = formation_energy_const
        # distance data
        self.min_distance = min_distance
        self.interval = interval
        self.span = span

# class for binding energy?

# this class stores the HSE screen results
# 'max_relaxation': float, 'average_relaxation': float, 'delta_Q': float,vv
class ScreenHSEResult(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation, 'charge': int, 'spin': float, 'ground_bool': bool, 'total_energy_pbe': float,  'total_energy_hse': float, 'excitation': str, 'excited_bool': bool, 'excited_spin': str, 'excited_converged': bool, 'ZPL_estimate_pbe': float, 'transition_intensity_pbe': float, 'partial_density_difference_pbe': float, 'polarization_pbe_abs_x': float, 'polarization_pbe_abs_y': float, 'polarization_pbe_abs_z': float, 'polarization_pbe_zpl_x': float, 'polarization_pbe_zpl_y': float, 'polarization_pbe_zpl_z': float, 'delta_R': float, 'delta_Q': float, 'ex_bool_hse': float, 'ZPL_estimate_hse': float, 'transition_intensity_hse': float, 'partial_density_difference_hse': float, 'polarization_hse_abs_x': float, 'polarization_hse_abs_y': float, 'polarization_hse_abs_z': float, 'polarization_hse_zpl_x': float, 'polarization_hse_zpl_y': float, 'polarization_hse_zpl_z': float, 'time': float}, index=['defect_key'])
    def __init__(self, defect_key, computation, charge, spin, ground_bool, total_energy_pbe, total_energy_hse, excitation, excited_bool, excited_spin, excited_converged, ZPL_estimate_pbe, transition_intensity_pbe, partial_density_difference_pbe, polarization_pbe_abs_x, polarization_pbe_abs_y, polarization_pbe_abs_z, polarization_pbe_zpl_x, polarization_pbe_zpl_y, polarization_pbe_zpl_z, delta_R, delta_Q, ex_bool_hse, ZPL_estimate_hse, transition_intensity_hse, partial_density_difference_hse, polarization_hse_abs_x, polarization_hse_abs_y, polarization_hse_abs_z, polarization_hse_zpl_x, polarization_hse_zpl_y, polarization_hse_zpl_z, time):
# max_relaxation, average_relaxation, delta_Q,
        self.defect_key = defect_key
        self.computation = computation
        self.charge = charge
        self.spin = spin
        # ground
        self.ground_bool = ground_bool
        self.total_energy_pbe = total_energy_pbe        
        self.total_energy_hse = total_energy_hse
        # excited
        self.excited_bool = excited_bool
        self.excited_spin = excited_spin
        self.excitation = excitation
        self.excited_converged = excited_converged
        # pbe
        self.ZPL_estimate_pbe = ZPL_estimate_pbe
        self.transition_intensity_pbe = transition_intensity_pbe
        self.partial_density_difference_pbe = partial_density_difference_pbe
        self.polarization_pbe_abs_x = polarization_pbe_abs_x
        self.polarization_pbe_abs_y = polarization_pbe_abs_y
        self.polarization_pbe_abs_z = polarization_pbe_abs_z
        self.polarization_pbe_zpl_x = polarization_pbe_zpl_x
        self.polarization_pbe_zpl_y = polarization_pbe_zpl_y
        self.polarization_pbe_zpl_z = polarization_pbe_zpl_z
        #self.max_relaxation = max_relaxation
        #self.average_relaxation = average_relaxation
        self.delta_R = delta_R # only for PBE
        self.delta_Q = delta_Q # only for PBE
        # hse
        self.ex_bool_hse = ex_bool_hse
        self.ZPL_estimate_hse = ZPL_estimate_hse
        self.transition_intensity_hse = transition_intensity_hse
        self.partial_density_difference_hse = partial_density_difference_hse
        self.polarization_hse_abs_x = polarization_hse_abs_x
        self.polarization_hse_abs_y = polarization_hse_abs_y
        self.polarization_hse_abs_z = polarization_hse_abs_z
        self.polarization_hse_zpl_x = polarization_hse_zpl_x
        self.polarization_hse_zpl_y = polarization_hse_zpl_y
        self.polarization_hse_zpl_z = polarization_hse_zpl_z
        # calculation time
        self.time = time


# This class is for meta data about the projects
class Project(httk.HttkObject):

    @httk.httk_typed_init({'title': str, 'date': str, 'producer': str, 'contributors': str, 'citations': str, 'counter': int}, index=['counter'])
    def __init__(self, title, date, producer, contributors, citations, counter):
        self.title = title
        self.date = date
        self.producer = producer
        self.contributors = contributors
        self.citations = citations
        self.counter = counter
        
        
