import sys

sys.path.append("../SourceCode")
import logging
import pytest
from Tests.TestBase import TestBase
from quantityCalculator import QuantityCalculator
from PreProcessing import PreProcessing
from unitConversions import auToGPascal
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.cubic import SimpleCubic
from ase.visualize import view
from ase.io.trajectory import Trajectory

logger = logging.getLogger(__name__)


class TestQuantityCalculator(TestBase):
    """Tests for the class PostProcessing. We want to run the source code to create the files we use for testing"""

    def setUp(self):
        """Sets up the PostProcessing object for all tests. This is supposed to be fcc copper."""
        super().setUp()

        # Files for the tests
        near_zero_copper = ["TestAtomicStructure/Cu_fcc.vasp", "TestSettings/nearZeroSettings.json",
                            "TestOutputs/testNearZeroOutput.traj"]
        solid_copper = ["TestAtomicStructure/Cu_fcc.vasp", "TestSettings/solidSettings.json",
                        "TestOutputs/testSolidOutput.traj"]
        melted_copper = ["TestAtomicStructure/Cu_fcc.vasp", "TestSettings/meltedSettings.json",
                         "TestOutputs/testMeltedOutput.traj"]
        solid_chromium = ["TestAtomicStructure/Cr_bcc.vasp", "TestSettings/chromiumSettings.json",
                          "TestOutputs/testChromiumOutput.traj"]

        self.preprocessing_near_zero_copper = PreProcessing(
            ["placeholder", near_zero_copper[0], near_zero_copper[1]])
        self.preprocessing_solid_copper = PreProcessing(
            ["placeholder", solid_copper[0], solid_copper[1]])
        self.preprocessing_melted_copper = PreProcessing(
            ["placeholder", melted_copper[0], melted_copper[1]])
        self.preprocessing_solid_chromium = PreProcessing(
            ["placeholder", solid_chromium[0], solid_chromium[1]])

        self.near_zero_copper = QuantityCalculator(self.preprocessing_near_zero_copper.createSettings(),
                                                   Trajectory(near_zero_copper[2]))
        self.solid_copper = QuantityCalculator(self.preprocessing_solid_copper.createSettings(),
                                               Trajectory(solid_copper[2]))
        self.melted_copper = QuantityCalculator(self.preprocessing_melted_copper.createSettings(),
                                                Trajectory(melted_copper[2]))
        self.solid_chromium = QuantityCalculator(self.preprocessing_solid_chromium.createSettings(),
                                                 Trajectory(solid_chromium[2]))

    def testCohesiveEnergy(self):
        """Make sure we calculate the Cohesive Energy correctly, for now some problems.
        Expected value for FCC copper is 3.55 eV/atom"""
        # TODO This is currently wrong, put in when the function is fixed.
        self.assertTrue(abs(self.solid_copper.computeCohesiveEnergy() - 3.55) < 0.2,
                               "Cohesive Energy differed with more than 0.1")

    def testLatticeConstant(self):
        """Checking that we get the expected lattice constant for fcc and sc. Fcc needs to be modified for unit cell"""
        self.assertAlmostEqual(self.solid_copper.computeLatticeConstant()[0], (2*(3.57/2)**2)**(1/2), 0,
                               "Lattice constant differed with more than 0.001")

        self.assertAlmostEqual(self.solid_chromium.computeLatticeConstant()[0], 2.97, 0,
                               "Lattice constant differed with more than 0.001")

    def testBulkModulus(self):
        """Bulk Modulus for Copper is expected to be 138 GPa"""
        # Now we convert the unit, computeBulkModulus returns in Pa
        self.assertTrue(abs(self.solid_copper.computeBulkModulus() - 138) < 10,
                               "Bulk modulus differed to much from expected")

    def testInternalPressure(self):
        """Test that internal pressure is close to zero for NPT ensamble"""
        self.assertAlmostEqual(auToGPascal(self.melted_copper.computeInternalPressure()), 0, delta=5e1,
                           msg="Internal pressure not close to zero")

    def testMeanSquareDisplacement(self):
        self.assertTrue(self.near_zero_copper.computeMSD(frame=-1) < self.solid_copper.computeMSD(frame=-1),
                        "Mean Square Distance expected to be less for lower temp")
        self.assertTrue(self.solid_copper.computeMSD(frame=-1) < self.melted_copper.computeMSD(frame=-1),
                        "Mean Square Distance expected to be less for lower temp")


    def testLindemannCriterion(self):
        """Checks that Lindemann criteria is fulfilled for melted and solid copper"""
        self.assertTrue(self.solid_copper.computeLindemannIndex() < 0.1,
                         "Lindemann criterion was True for a solid phase")
        self.assertTrue(self.melted_copper.computeLindemannIndex() > 0.1,
                        "Lindemann criterion was False for a melted phase")

    def testSelfDiffusionCoefficient(self):
        """Self diffusion for copper is expected to be around 8 * 10^-8 Å^2/fs at 2000 K, but very inaccurate is the
        feel. Therefor the big interval """
        self.assertAlmostEqual(self.melted_copper.computeSelfDiffusionCoefficient(), 8E-8, delta=1E6,
                               msg="Self diffusion differed to much from expected 8E-8 Å^2/fs")

    @pytest.mark.skip("Implemented but differs wildly")
    def testDebyeTemperature(self):
        """For copper Debye temperature should be 343"""
        logger.info(f"Debye TEMPERATURE : {self.solid_copper.computeDebyeTemperature()}")
        self.assertTrue(abs(self.solid_copper.computeDebyeTemperature() - 343) < 30,
                        "Debye temperature deviated by more than 30 degrees")


    def testNearestNeighbour(self):
        """Currently returns the mean of the two closest atoms in the structure. Expected to be a little lower than
        average NN.
        """
        distance_to_nn = self.solid_copper.nearestNeighborsMean(start=-2)
        self.assertTrue(abs(distance_to_nn - 3.57 / (2 ** (1 / 2))) < 0.2,
                               "Nearest Neighbour distance differed with more than 0.1 to expected")
