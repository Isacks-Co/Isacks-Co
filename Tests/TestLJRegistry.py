import sys

sys.path.append("..")
from Tests.TestBase import TestBase
from SourceCode.LJRegistry import  LJParams, LJ_DB
import logging
logger = logging.getLogger(__name__)


class TestLJParams(TestBase):

    def testNonEmptyArgs(self):
        logger.info("Testing setting dictionary with empty sigma and epsilon parameters")
        known_material = LJParams(material="Ar")
        self.assertIsInstance(known_material, dict)
        self.assertGreater(known_material["epsilon_eV"], 0.0)
        self.assertGreater(known_material["sigma_A"], 0.0)


    def testRadiusArgs(self):
        known_material = LJParams(material="Ar")
        self.assertLess(known_material["ro_A"], known_material["rc_A"])


    def testCaseInsensitive(self):
        logger.info("Testing case insesitivity")
        self.assertEqual(LJParams("ar"), LJParams("Ar"))
        self.assertEqual(LJParams("AR"), LJParams("aR"))
        self.assertEqual(LJParams("ar"), LJParams("AR"))


    def testOverride(self):
        logger.info("Checking that the override values are used instead of the ones from database")
        known = LJParams(material="Ar")
        override = LJParams(material="Ar", epsilon_eV= 0.02, sigma_A= 3.5)
        self.assertNotEqual(known["epsilon_eV"], override["epsilon_eV"])
        self.assertNotEqual(known["sigma_A"], override["sigma_A"])




if __name__ == "__main__":
    import unittest
    unittest.main(verbosity=2)


