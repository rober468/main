# Rate_Constant.py
# Calculates the trajectory over time for the first ten runs for R_0(t).

# Load modules
import numpy as np
import matplotlib.pyplot as plt

# Load constants
timestep = 0.005
figx = 8.0
figy = 5.0

# Load data
a = np.loadtxt('trajectory_R_0.dat')
b = np.arange(timestep,len(a[:,0])*timestep+timestep,timestep)


# Plot data
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.plot(b,a[:,0],'k-')
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$\tilde{R_0}(t)$")
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
ax.set_title(r"Figure 6. Example of trajectory with $\tilde{\beta}=$6 and 1000 oscillators.",fontsize=12)
plt.tight_layout(pad=0.1)
plt.savefig("trajectory_R_0.pdf",format="pdf",dpi=600)
