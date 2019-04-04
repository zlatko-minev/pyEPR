"""
      Python (py) Energy-Participation-Ratio (EPR) package
"""

from setuptools import setup, find_packages

doclines = __doc__.split('\n')

setup(name='pyEPR',
      version='0.6',
      description = doclines[0],
      long_description = '\n'.join(doclines[2:]),
      author='Zlatko Minev',
      packages=['pyEPR'],
      author_email='zminev@gmail.com',
      license=' BSD-3-Clause',
      install_requires=['numpy','pandas','pint','matplotlib']
      )
