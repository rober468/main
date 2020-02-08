# k_TST.py
# Calculates the transition state theory rate constant analytically.

# Load modules
import math as ma
import numpy as np

# Load constants
a = 3.				# Quartic potential parameter
b = 2.				# Quartic potential parameter
xi = 4.				# Kondo coupling constant
omega_max = 3.			# Cutoff frequency
beta = 6.			# Inverse temperature parameter

# Calculate derived quantities
Ra = -ma.sqrt(b/a+xi/a*(1.-np.exp(-omega_max)))
Udp = 3.*a*Ra**2 - ( b+xi*(1.-np.exp(-omega_max)) )
W_Ra = a/4*Ra**4-b/2*Ra**2-(1./2.)*xi*(1.-np.exp(-omega_max))*Ra**2
n_A = ma.sqrt(2*ma.pi/(beta*Udp))*np.exp(-beta*W_Ra)

# Calculate k_TST
k_TST = ( ma.sqrt(Udp)/ (2.*ma.pi) ) * np.exp(beta*W_Ra)
print k_TST
