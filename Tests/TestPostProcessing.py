import sys

sys.path.append("..")
import logging
import pytest
from Tests.TestBase import TestBase
from SourceCode.PostProcessing import PostProcessing
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.cubic import SimpleCubic
from ase.visualize import view

logger = logging.getLogger(__name__)


class TestPostProcessing(TestBase):
    """Tests for the class PostProcessing. We want to run the source code to create the files we use for testing"""

    def setUp(self):
        """Sets up the PostProcessing object for all tests. This is supposed to be fcc copper."""
        super().setUp()
        self.five_K_copper = PostProcessing("TestSettings/nearZeroSettings.json", "TestOutputs/testNearZeroOutput.traj")
        self.solid_copper = PostProcessing("TestSettings/solidSettings.json", "TestOutputs/testSolidOutput.traj")
        self.melted_copper = PostProcessing("TestSettings/meltedSettings.json", "TestOutputs/testMeltedOutput.traj")
        self.solid_chromium = PostProcessing("TestSettings/chromiumSettings.json",
                                             "TestOutputs/testChromiumOutput.traj")
        self.solid_polonium = PostProcessing("TestSettings/poloniumSettings.json",
                                             "TestOutputs/testPoloniumOutput.traj")

    def testCohesiveEnergy(self):
        """Make sure we calculate the Cohesive Energy correctly, for now some problems.
        Expected value for FCC copper is 3.55 eV/atom"""
        # TODO This is currently wrong, put in when the function is fixed.
        self.assertAlmostEqual(self.solid_copper.computeCohesiveEnergy(return_unit="eV_per_atom"), 3.55, 1,
                               "Cohesive Energy differed with more than 0.1")

    def testLatticeConstant(self):
        """Checking that we get the expected lattice constant for fcc, bcc and sc"""
        self.assertAlmostEqual(self.solid_copper.computeLatticeConstant(), 3.57, 0,
                               "Lattice constant differed with more than 0.001")

        self.assertAlmostEqual(self.solid_chromium.computeLatticeConstant(), 2.97, 0,
                               "Lattice constant differed with more than 0.001")

        # This is assumed to be a BBC, but actually SC
        self.assertAlmostEqual(self.solid_polonium.computeLatticeConstant(), 3.35, 0,
                               "Lattice constant differed with more than 0.001")

    def testBulkModulus(self):
        """Bulk Modulus for Copper is expected to be 138 GPa, currently we get 134. Is allowed to diff with +/- 0.5"""
        # Now we convert the unit, computeBulkModulus returns in Pa
        self.assertAlmostEqual(self.solid_copper.computeBulkModulus() * 1e-9, 134, 0,
                               "Bulk modulus differed to much from expected")

    def testInternalPressure(self):
        """Internal Pressure for Copper is expected to be somewhere around 20-50 GPa (I think)
        TODO Expand this test, to week"""
        self.assertGreater(self.solid_copper.computeInternalPressure() * 1e-9, 20, "Internal pressure lower than 20 GPa")
        self.assertLess(self.solid_copper.computeInternalPressure() * 1e-9, 50, "Internal pressure higher than 50 GPa")

    def testMeanSquareDistance(self):
        logger.info(f"MEAN SQUARE DISTANCE {self.solid_copper.computeMSD()}")

        self.assertTrue(self.five_K_copper.computeMSD() < self.solid_copper.computeMSD(),
                        "Mean Square Distance expected to be less for lower temp")
        self.assertTrue(self.solid_copper.computeMSD() < self.melted_copper.computeMSD(),
                        "Mean Square Distance expected to be less for lower temp")

        # TODO Want to implement one that is stable and one that explodes. Waiting for test library.

    def testLindemannCriterion(self):
        # TODO Fix test sample directory when merge with test library is finished. Remember to fix all inits.
        # Maybe initialize many in the setUp function could be an idea for cleaner code.

        self.assertFalse(self.solid_copper.computeLindemannCriterion(),
                         "Lindemann criterion was True for a solid phase")
        self.assertTrue(self.melted_copper.computeLindemannCriterion(),
                        "Lindemann criterion was False for a melted phase")

    def testSelfDiffusionCoefficient(self):
        # Seems like we multiply by 1e-8 for some reason I don't understand. Needs to be discussed.
        # ON HOLD
        logger.info(f"Diffusion Coeff --- {self.solid_copper.computeSelfDiffusionCoefficient()}")

    def testDebyeTemperature(self):
        """Debye model is good for low temperatures. For copper, Debye temperature should be 343
        Currently we observe 386 which is pretty good. The test is made for our current result, with deviation of
        50 degrees K  allowed"""
        self.assertTrue(abs(self.five_K_copper.computeDebyeTemperature() - 343) < 50,
                        "Debye temperature deviated by more than 10 degrees")

    #def testShearModulus(self):
        # TODO Not in postprocessing anymore, needs to be disused

    def testNearestNeighbour(self):
        """Want to make sure we calculate NN correctly. Currently, the function is a bit messy, I believe it can be made
        alot better.
        TODO When implemented, change this test to test for multiple structures
        """
        distance_to_nn = self.solid_copper.nearestNeightborsMean()
        self.assertAlmostEqual(distance_to_nn, 3.57 / (2 ** (1 / 2)), 1,
                               "Nearest Neighbour distance differed with more than 0.1 to expected")
