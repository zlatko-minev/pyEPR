"""
pyEPR constants and convenience definitions.

@author: Zlatko Minev
"""
# pylint: disable=invalid-name

from scipy.constants import Planck, elementary_charge, epsilon_0, pi  # pylint: disable=unused-import

# Pi
π = pi

# Reduced Planks constant
ħ = hbar = Planck/(2*pi)

# Reduced Flux Quantum  (3.29105976 × 10-16 Webers)
ϕ0 = fluxQ = ħ / (2*elementary_charge)

# Magnitude of the electric charge carried by a single electron
e_el = elementary_charge
