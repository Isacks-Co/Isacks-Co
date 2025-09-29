import unittest
import logging

logger = logging.getLogger(__name__)


class TestBase(unittest.TestCase):
    """This is a testclass that is inherited by all other tests. For now limited functionality but can be used setup
    tests
    """

    def setUp(self):
        self.method_name = unittest.TestCase.id(self)
        logger.info(f"Running test: {self.method_name}")

    def tearDown(self):
        logger.info(f"Finished test: {self.method_name}")
