import unittest
import sys

sys.path.insert(0, "..")  # noqa
import pyEPR as epr


class TestProjectInfo(unittest.TestCase):
    """Test pyEPR.project_info.py"""

    def setUp(self):
        path_to_project = r"..\_example_files"
        try:
            self.pinfo = epr.ProjectInfo(
                project_path=path_to_project,
                project_name="pyEPR_tutorial1",
                design_name="1. single_transmon",
            )
        except:
            assert ConnectionError("Failed to connect to HFSS. Opening it manually")

    def test_dissipative(self):
        """Test change of _Dissipative from a class to a dict with deprecation warnings"""
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__getattr__,
            "mot_exist",
            msg="Failed calling non-existing attr",
        )
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__getitem__,
            "not_exist",
            msg="Failed calling non-existing item",
        )
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__setattr__,
            "seams",
            1,
            msg="Failed setting invalid attr",
        )
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__setitem__,
            "seams",
            1,
            msg="Failed setting invalid item",
        )
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__setitem__,
            "seams",
            ["a"],
            msg="Failed setting item to non-existing HFSS obj",
        )
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__setattr__,
            "seams",
            ["a"],
            msg="Failed setting attr to non-existing HFSS obj",
        )
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__setattr__,
            "not_exist",
            1,
            msg="Failed setting invalid value and attr",
        )
        self.assertRaises(
            Exception,
            self.pinfo.dissipative.__setitem__,
            "not_exist",
            1,
            msg="Failed setting invalid value and key",
        )
        self.pinfo.dissipative["seams"]
        self.pinfo.dissipative["seams"] = []
        self.pinfo.dissipative["seams"] = ["substrate"]
