# Rate_Constant_Temp_Dependence.py
# Calculates time-dependent rate-constant from simulation data.

# Load modules
import numpy as np
import matplotlib.pyplot as plt

# Load constants
timestep = 0.005		# MD timestep
figx = 8.0			# Figure x dim.
figy = 5.0			# Figure y dim.	

# Load data
a = np.loadtxt('Time_Dep_Rate_3_100.dat')
b = np.loadtxt('Time_Dep_Rate_6_100.dat')
c = np.loadtxt('Time_Dep_Rate_9_100.dat')
d = np.loadtxt('Time_Dep_Rate_12_100.dat')
e = np.loadtxt('Time_Dep_Rate_15_100.dat')
f = np.arange(timestep,len(a[:,0])*timestep+timestep,timestep)

# Plot data
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.plot(f,a[:,0],'k-',label=r"$\tilde{\beta}$=3")
ax.plot(f,b[:,0],'b-',label=r"$\tilde{\beta}$=6")
ax.plot(f,c[:,0],'r-',label=r"$\tilde{\beta}$=9")
ax.plot(f,a[:,0]-1.96*a[:,1],'k--')
ax.plot(f,a[:,0]+1.96*a[:,1],'k--')
ax.plot(f,b[:,0]-1.96*b[:,1],'b--')
ax.plot(f,b[:,0]+1.96*b[:,1],'b--')
ax.plot(f,c[:,0]-1.96*c[:,1],'r--')
ax.plot(f,c[:,0]+1.96*c[:,1],'r--')
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$k_{f}(t)$")
ax.set_title(r"Figure 2. Temperature dependence of rate constant over time for low values of $\tilde{\beta}$, $N{osc}=$100.",fontsize=10)
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
plt.tight_layout(pad=0.1)
plt.legend(loc='upper right',fontsize=12)
plt.savefig("Plots/Rate_Constant_Temp_Dep_Low.pdf",format="pdf",dpi=600)
plt.clf()
ax = fig.add_subplot(111)
ax.plot(f,c[:,0],'k-',label=r"$\tilde{\beta}$=9")
ax.plot(f,d[:,0],'b-',label=r"$\tilde{\beta}$=12")
ax.plot(f,e[:,0],'r-',label=r"$\tilde{\beta}$=15")
ax.plot(f,c[:,0]-1.96*c[:,1],'k--')
ax.plot(f,c[:,0]+1.96*c[:,1],'k--')
ax.plot(f,d[:,0]-1.96*d[:,1],'b--')
ax.plot(f,d[:,0]+1.96*d[:,1],'b--')
ax.plot(f,e[:,0]-1.96*e[:,1],'r--')
ax.plot(f,e[:,0]+1.96*e[:,1],'r--')
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$k_{f}(t)$")
ax.set_title(r"Figure 3. Temperature dependence of rate constant over time for high values of $\tilde{\beta}$, $N{osc}=$100.",fontsize=10)
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
plt.tight_layout(pad=0.1)
plt.legend(loc='upper right',fontsize=12)
plt.savefig("Plots/Rate_Constant_Temp_Dep_High.pdf",format="pdf",dpi=600)
