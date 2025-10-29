import sys

sys.path.append("..")
from Tests.TestBase import TestBase
from SourceCode.newPostProcessing import PostProcessing

class TestPostProcessing(TestBase):
    """Tests for the class PostProcessing"""
    def testInitilizer(self):
        """Only test for now, but more can be implemented"""
        PostViz = PostProcessing("testsettings.json", "testoutput.traj")