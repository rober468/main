# Poisson_Dichotomous.py
# Simulates a particles undergoing dynamics with Poisson
# dichotomous noise h(t)=+/-delta with equation of motion:
# dx/dt=-kx+h(t)
# where jumps between states are exponentially distributed
# with characteristic time tau. Calculation of steady state
# of position is also calculated.

# Load modules
import numpy as np
import matplotlib.pyplot as plt

# Constants
traj = 20000 	                   	# Number of trajectories
time = 2500.				# Total time
timestep = 0.1				# Timestep
tau = 2.0				# Characteristic time
delta = 0.5				# Dichotomous noise value
k = 1.0					# Potential parameter
fpos = np.zeros(traj)			# Final position
figx = 4.0                      	# x dim. for fig. size
figy = 3.25                     	# y dim. for fig. size

# Seed MT PRNG and make array of random numbers
num = 2798048
np.random.seed(num)

# Time evolution
for i in np.arange(traj):
    print i
    ctime=0.
    pos = (2.*np.random.random_sample() -1.)*10.
    if (np.random.random_sample() < 0.5):
        ostate = 0
    else:
        ostate = 1
    while ( ctime < time ):
        transt = -tau*np.log(np.random.random_sample())
        steps = int(transt/timestep)
        for j in np.arange(steps):
            ctime = ctime + timestep
            if (ostate == 0):
                pos = pos + (-k*pos+delta)*timestep
            if (ostate == 1):
                pos = pos + (-k*pos-delta)*timestep
        if (ostate == 0):
            nstate = 1
        if (ostate == 1):
            nstate = 0
        ostate = nstate
    fpos[i] = pos

# Plot histogram for steady state
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
#ax.set_xlim(xmin,xmax)
#ax.set_ylim(0,1)
ax.set_xlabel(r"$x$",fontsize=10)
ax.set_ylabel(r"$p_{s}(x)$",fontsize=10)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
weights = np.ones_like(fpos)/len(fpos)
ax.hist(fpos[:],bins=250,weights=weights)
ax.set_title('Problem 2c - Histogram of $p_{s}(x)$',fontsize=10)
plt.tight_layout(pad=0.1)
plt.savefig('SS_Histogram_Poisson_Dichotomous.pdf',dpi=600)
