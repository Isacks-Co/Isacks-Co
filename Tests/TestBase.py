import unittest
import logging

logger = logging.getLogger(__name__)

class TestBase(unittest.TestCase):
    def setUp(self):
        self.method_name = unittest.TestCase.id(self)
        logger.info(f"Running test: {self.method_name}")

    def tearDown(self):
        logger.info(f"Finished test: {self.method_name}")