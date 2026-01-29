# Unit Tests Module (work-in-progress)
This module is used for development unit testing, to be expanded in the future.

Tests are written using the Python [built-in unittest module](https://docs.python.org/3/library/unittest.html), but are executed using [pytest](https://pytest.org/) for better test discovery and reporting.

## Running Tests

**Using uv (recommended):**
```sh
uv sync --group test
uv run pytest
```

**Using pytest directly:**
```sh
pytest
```

**Running specific tests:**
```sh
pytest tests/test_project_info.py           # Run all tests in a file
pytest tests/test_project_info.py -k test_dissipative  # Run a specific test
```

**Using unittest (legacy):**
```sh
python -m unittest                           # Run all tests
python -m unittest tests.test_project_info   # Run tests from a specific module
```

Note: Tests require Python >=3.10,<3.13 as specified in the project requirements.
