import sys

sys.path.append("../SourceCode")
from Tests.TestBase import TestBase
from SourceCode.PostProcessing import PostProcessing

class TestPostProcessing(TestBase):
    """Tests for the class PostProcessing"""
    def testInitilizer(self):
        """Only test for now, but more can be implemented"""
        PostViz = PostProcessing("testsettings.json", "testoutput.traj")