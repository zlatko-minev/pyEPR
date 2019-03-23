"""Experiment control and data acquisition software for the Parity measurment.
"""

from setuptools import setup, find_packages

doclines = __doc__.split('\n')

setup(name='pyEPR',
      version='0.0',
      description = doclines[0],
      long_description = '\n'.join(doclines[2:]),
      author='Zlatko Minev',
      packages=['pyEPR'],
      author_email='???',
      license='???',
      install_requires=['numpy','pandas','pint','matplotlib']
      )