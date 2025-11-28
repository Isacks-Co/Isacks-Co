import sys

sys.path.append("../SourceCode")
import logging
import pytest
import unittest
import numpy as np
from TestBase import TestBase
from QuantityConstants import E_TOT_SEQ_NVT, E_TOT_SEQ_NVE, MSD_SEQ, C_MATRIX
from Utils.unitConversions import GPascalToAu, auToGPascal
from quantityCalculator import QuantityCalculator
from PreProcessing import PreProcessing
from Utils.unitConversions import auToGPascal
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.cubic import SimpleCubic
from ase.visualize import view
from ase.io.trajectory import Trajectory

logger = logging.getLogger(__name__)


class TestQuantityCalculator(TestBase):
    """Tests that we get close the real quantities when calculating with mostly real data.
    All tests currently use the copper as the reference. Some data can be found in QuantityConstants.py
    that was too big to implement in this file. The sequence data comes from simulation done by the program
    No tests implemented for nearestNeighbourMean and computeBulkModulus as they are expected to change."""

    def setUp(self):
        """Sets up the physical values that will be used. Simulates fcc copper, 5x5x5 supercell"""
        super().setUp()
        self.T = 300
        self.B = GPascalToAu(138)
        self.G = GPascalToAu(48)
        self.E = GPascalToAu(130)
        self.T_D = 343
        self.E_tot_seq_NVT = E_TOT_SEQ_NVT
        self.E_tot_seq_NVE = E_TOT_SEQ_NVE
        self.msd_seq = MSD_SEQ
        self.sample_spacing = 1
        self.total_atoms = 216
        self.total_mass_amu = self.total_atoms * 63.546
        self.total_volume = self.total_atoms * 11.777
        self.msd = 8.016646992112521e-33
        self.nearest_neighbour = 2.556
        self.C_matrix = np.asarray(C_MATRIX)

    def testSpecificHeatNVT(self):
        """Test to check that we get reasonable numbers from computeSpecificHeatNVT
        Expected 0.385 J/g*K or 3.9902E-6 eV/au*K"""
        heat_capacity = QuantityCalculator.computeSpecificHeatNVT(self.E_tot_seq_NVT, self.total_mass_amu, self.T)
        logger.debug(f"This is specific heat capacity from NVT data : {heat_capacity}")
        self.assertTrue(abs(heat_capacity - 3.9902E-6) < 4E-6)
        self.assertGreater(heat_capacity, 1E-7)

    @pytest.mark.skip("Not working as intended, on hold")
    def testSpecificHeatNVE(self):
        """Test to check that we get reasonable numbers from computeSpecificHeatNVE
        Expected 0.385 J/g*K or 3.9902E-6 eV/au*K"""
        heat_capacity = QuantityCalculator.computeSpecificHeatNVE(self.E_tot_seq_NVE, self.total_mass_amu, self.T)
        logger.debug(f"This is specific heat capacity from NVE data : {heat_capacity}")
        self.assertTrue(abs(heat_capacity - 3.9902E-6) < 4E-6)
        self.assertGreater(heat_capacity, 1E-7)

    def testSelfDiffusionCoefficient(self):
        """Compute self diffusion coefficient. Can be expanded by data from a melted example as well"""
        diffusion_coeff = QuantityCalculator.computeSelfDiffusionCoefficient(self.msd_seq, self.sample_spacing)
        logger.debug(f"This is the self diffusion coefficient : {diffusion_coeff}")
        self.assertAlmostEqual(diffusion_coeff, 8E-8, delta=1E6,
                               msg="Self diffusion differed to much from expected 8E-8 Å^2/fs")

    def testDebyeTemperature(self):
        """Checks that we get within a 5 degree margin of the true debye temperature for copper"""
        debye_temp = QuantityCalculator.computeDebyeTemperature(
            self.total_volume, self.total_mass_amu, self.total_atoms, self.G, self.B)
        logger.debug(f"This is the debye temperature : {debye_temp}")
        self.assertAlmostEqual(debye_temp, 343, delta=5, msg="Calculation of debye temperature differed to much")

    def testLindemannIndex(self):
        """Test for melted phase missing, should be implemented"""
        lindemann_index = QuantityCalculator.computeLindemannIndex(self.msd, self.nearest_neighbour)
        logger.debug(f"This is the lindemann_index : {lindemann_index}")
        self.assertAlmostEqual(lindemann_index, 0, 4, msg="Lindemann index to high for solid state data")

    def testCalculateModuli(self):
        """Test that we get the right moduli from a C matrix"""
        B, G, E = QuantityCalculator.calculateModuli(self.C_matrix)
        B, G, E = auToGPascal(B), auToGPascal(G), auToGPascal(E)
        logger.debug(f"This is bulk modulus (B) : {B}")
        logger.debug(f"This is the shear modulus (G) : {G}")
        logger.debug(f"This is the young modulus (E) : {E}")
        print(f"B, G, E = {B}, {G}, {E}")
        self.assertAlmostEqual(B, auToGPascal(self.B), delta=15,
                               msg="Bulk modulus from C matrix differed to much from real value")
        self.assertAlmostEqual(G, auToGPascal(self.G), delta=15,
                               msg="Shear modulus from C matrix differed to much from real value")
        self.assertAlmostEqual(E, auToGPascal(self.E), delta=15,
                               msg="Young modulus from C matrix differed to much from real value")
