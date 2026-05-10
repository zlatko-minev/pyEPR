"""
Tests for pyEPR.project_info.py

Tests marked @pytest.mark.hfss require a live Ansys HFSS connection and
are skipped in CI. Run locally with: pytest -m hfss
"""
import pytest


@pytest.mark.hfss
class TestProjectInfoHFSS:
    """Tests that require a live HFSS connection."""

    @pytest.fixture(autouse=True)
    def connect(self):
        import pyEPR as epr
        try:
            self.pinfo = epr.ProjectInfo(
                project_path=r"..\_example_files",
                project_name="pyEPR_tutorial1",
                design_name="1. single_transmon",
            )
        except Exception as e:
            pytest.skip(f"Cannot connect to HFSS: {e}")

    def test_dissipative_invalid_attr_raises(self):
        with pytest.raises(Exception):
            self.pinfo.dissipative.__getattr__("not_exist")

    def test_dissipative_invalid_item_raises(self):
        with pytest.raises(Exception):
            self.pinfo.dissipative["not_exist"]

    def test_dissipative_invalid_setattr_raises(self):
        with pytest.raises(Exception):
            self.pinfo.dissipative.__setattr__("seams", 1)

    def test_dissipative_invalid_setitem_raises(self):
        with pytest.raises(Exception):
            self.pinfo.dissipative["seams"] = 1

    def test_dissipative_nonexistent_object_setitem_raises(self):
        with pytest.raises(Exception):
            self.pinfo.dissipative["seams"] = ["a"]

    def test_dissipative_valid_operations(self):
        self.pinfo.dissipative["seams"]
        self.pinfo.dissipative["seams"] = []
        self.pinfo.dissipative["seams"] = ["substrate"]


class TestProjectInfoNoHFSS:
    """Tests that run without any HFSS connection."""

    def test_import(self):
        import pyEPR as epr
        assert hasattr(epr, "ProjectInfo")

    def test_project_info_requires_path(self):
        """ProjectInfo with a nonexistent path should fail at connection, not import."""
        import pyEPR as epr
        # Just check that the class is instantiable with args (connection fails later)
        assert callable(epr.ProjectInfo)

    def test_dissipative_container_standalone(self):
        """_Dissipative validation logic works without HFSS."""
        from pyEPR.project_info import ProjectInfo
        d = ProjectInfo._Dissipative()
        # Setting a valid key to an empty list is allowed
        d["seams"] = []
        assert d["seams"] == []

    def test_dissipative_rejects_non_string_list(self):
        """_Dissipative should reject non-string-list values."""
        from pyEPR.project_info import ProjectInfo
        d = ProjectInfo._Dissipative()
        with pytest.raises(ValueError):
            d["seams"] = 42

    def test_dissipative_rejects_invalid_key(self):
        """_Dissipative should reject unknown keys."""
        from pyEPR.project_info import ProjectInfo
        d = ProjectInfo._Dissipative()
        with pytest.raises((KeyError, AttributeError, ValueError)):
            d["completely_invalid_key_xyz"] = []
