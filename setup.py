"""
Python (py) Energy-Participation-Ratio (EPR) package
pyEPR is an open source, BSD-licensed library providing high-efficiency,
easy-to-use analysis functions and automation for the design of quantum
chips based on superconducting quantum circuits, both distributed and lumped.
pyEPR interfaces the classical distributed microwave analysis with that of
quantum structures and Hamiltonians. It is chiefly based on the energy participation
ratio approach; however, it has since v0.4 extended to cover a broad range of
design approaches. pyEPR straddles the analysis from Maxwell’s to Schrodinger’s
equations, and converts the solutions of distributed microwave (typically eigenmode
simulations) to a fully diagonalized spectrum of the energy levels, couplings,
and key parameters of a many-body quantum Hamiltonian.

Read the docs: https://pyepr-docs.readthedocs.io/en/latest/
Github page: https://github.com/zlatko-minev/pyEPR
"""

from pathlib import Path
from setuptools import setup, find_packages

here = Path(__file__).parent.absolute()

# Get the long description from the README file
with open(here / "README.md", encoding="utf-8") as f:
    long_description = f.read()

with open(here / "requirements.txt", encoding="utf-8") as f:
    requirements = f.read().splitlines()

doclines = __doc__.split('\n')

setup(
    name='pyEPR-quantum',
    version='0.8.5.7',
    description=doclines[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Zlatko K. Minev',
    packages=find_packages(),
    author_email='zlatko.minev@aya.yale.edu',
    maintainer='Zlatko Minev, pyEPR team',
    license='BSD-3-Clause',
    url=r'https://github.com/zlatko-minev/pyEPR',
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS", "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering", "Environment :: Console",
        "License :: OSI Approved :: Apache Software License"
    ],
    python_requires=">=3.5, <4",
    # install_requires=['numpy','pandas','pint','matplotlib','addict','sympy','IPython'],
    install_requires=requirements)
