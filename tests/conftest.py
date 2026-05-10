import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "hfss: marks tests that require a live Ansys HFSS connection (skipped in CI)"
    )
