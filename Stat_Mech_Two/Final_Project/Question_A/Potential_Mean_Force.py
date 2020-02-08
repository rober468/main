# Potential_Mean_Force.py
# Plots the potential of mean force for the quartic, bistable
# system coupled to a harmonic bath.
# Potential of mean force given by:
# W(R0) = 0.25*a*R0**4 - 0.5*b*R0**2 - 0.5*xi*(1-exp[-omega_max])*R0**2

# Load modules
import numpy as np
import matplotlib.pyplot as plt
import math as ma

# Load constants
a = 3.				# Quartic potential parameter
b = np.array([-2,2])		# Quartic potential parameter
omega_max = 3.			# Cutoff frequency
xi = 1.				# Kondo coupling constant

# Plot parameters
xmin = -2.0			# Minimum value for R0
xmax = 2.0			# Maximum value for R0
ymin = -2.			# Below minimum for W(R_0)
ymax = 16.			# Above maximum for W(R_0)
xfrac = 100.			# Number of intervals
xdiff = (xmax-xmin)/xfrac	# Interval length
figx = 8.0			# Figure size, x
figy = 5.0			# Figure size, y

# Potential calculation
R0 = np.arange(xmin,xmax+xdiff,xdiff)
W_Osc = np.zeros([len(R0),len(b)])
W_No_Osc = np.zeros([len(R0),len(b)])
for i in np.arange(0,2,1):
    W_Osc[:,i] = 0.25*a*R0[:]**4 - 0.5*b[i]*R0[:]**2 - (1./2.)*xi*(1.-np.exp(-omega_max))*R0[:]**2
    W_No_Osc[:,i] = 0.25*a*R0[:]**4 - 0.5*b[i]*R0[:]**2

# Plot potentials
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.plot(R0,W_Osc[:,0],'c-',label=r"Osc.: $\tilde{a}=3 , \tilde{b}=-2$")
ax.plot(R0,W_Osc[:,1],'b-',label=r"Osc.: $\tilde{a}=3 , \tilde{b}=2$")
ax.plot(R0,W_No_Osc[:,0],'c--',label=r"No Osc.: $\tilde{a}=3 , \tilde{b}=-2$")
ax.plot(R0,W_No_Osc[:,1],'b--',label=r"No Osc.: $\tilde{a}=3 , \tilde{b}=2$")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r"$\tilde{R}_{0}$")
ax.set_ylabel(r"$W(\tilde{R}_{0})$")
ax.set_title("Figure 1. Potential of mean force as a function of reaction coordinate.",fontsize=14)
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
plt.legend(loc="upper center",fontsize=12)
plt.tight_layout(pad=0.1)
plt.savefig("Potential_Mean_Force.pdf",format="pdf")
