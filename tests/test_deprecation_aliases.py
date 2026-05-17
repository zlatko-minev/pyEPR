"""
Tests for DeprecationWarning on the old pyEPR class name aliases:
  - Project_Info       → ProjectInfo
  - pyEPR_HFSSAnalysis → DistributedAnalysis
  - pyEPR_Analysis     → QuantumAnalysis
"""
import warnings
import pytest


class TestTopLevelAliases:
    """pyEPR.Project_Info etc. should warn and return the right class."""

    def test_project_info_alias_warns(self):
        import pyEPR
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cls = pyEPR.Project_Info
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "Project_Info" in str(w[0].message)
        assert "ProjectInfo" in str(w[0].message)

    def test_project_info_alias_returns_correct_class(self):
        import pyEPR
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            cls = pyEPR.Project_Info
        assert cls is pyEPR.ProjectInfo

    def test_hfss_analysis_alias_warns(self):
        import pyEPR
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cls = pyEPR.pyEPR_HFSSAnalysis
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "pyEPR_HFSSAnalysis" in str(w[0].message)
        assert "DistributedAnalysis" in str(w[0].message)

    def test_hfss_analysis_alias_returns_correct_class(self):
        import pyEPR
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            cls = pyEPR.pyEPR_HFSSAnalysis
        assert cls is pyEPR.DistributedAnalysis

    def test_quantum_analysis_alias_warns(self):
        import pyEPR
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cls = pyEPR.pyEPR_Analysis
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "pyEPR_Analysis" in str(w[0].message)
        assert "QuantumAnalysis" in str(w[0].message)

    def test_quantum_analysis_alias_returns_correct_class(self):
        import pyEPR
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            cls = pyEPR.pyEPR_Analysis
        assert cls is pyEPR.QuantumAnalysis

    def test_unknown_attribute_raises(self):
        import pyEPR
        with pytest.raises(AttributeError):
            _ = pyEPR.nonexistent_name_xyz


class TestCoreModuleAliases:
    """pyEPR.core.Project_Info etc. should warn and return the right class."""

    def test_project_info_alias_warns(self):
        from pyEPR import core
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cls = core.Project_Info
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "Project_Info" in str(w[0].message)

    def test_project_info_alias_returns_correct_class(self):
        from pyEPR import core
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            cls = core.Project_Info
        assert cls is core.ProjectInfo

    def test_hfss_analysis_alias_returns_correct_class(self):
        from pyEPR import core
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            cls = core.pyEPR_HFSSAnalysis
        assert cls is core.DistributedAnalysis

    def test_quantum_analysis_alias_returns_correct_class(self):
        from pyEPR import core
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            cls = core.pyEPR_Analysis
        assert cls is core.QuantumAnalysis

    def test_unknown_attribute_raises(self):
        from pyEPR import core
        with pytest.raises(AttributeError):
            _ = core.nonexistent_name_xyz


class TestCanonicalNamesNoWarning:
    """Canonical names must not emit any warnings."""

    def test_project_info_no_warning(self):
        import pyEPR
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _ = pyEPR.ProjectInfo
        dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
        assert dep_warnings == []

    def test_distributed_analysis_no_warning(self):
        import pyEPR
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _ = pyEPR.DistributedAnalysis
        dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
        assert dep_warnings == []

    def test_quantum_analysis_no_warning(self):
        import pyEPR
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _ = pyEPR.QuantumAnalysis
        dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
        assert dep_warnings == []
