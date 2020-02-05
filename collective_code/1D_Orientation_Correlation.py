# 1D_Orientation_Correlation.py
#
# Calculates the average of the orientation for atom i with other atom at distance j,
# over all particles.

# Load libraries
import h5py
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np

# Load data
data = h5py.File("data.h5","r")
dimer = data["Dimer"]
nuc = dimer["nucvec"]

# Load constants
freq = dimer.attrs["io_frequency"]
MDt = data.attrs["MD_time"]
figx = 8.0
figy = 4.0
xmin = 0.
yminp1 = -1.
ymaxp1 = 1.
yminp2 = -0.5
ymaxp2 = 1.

# Derived quantities
time = np.arange(0,freq*int(len(nuc[:,0,0]))*MDt,freq*MDt)
u = np.ones((len(nuc[:,0,0]),len(nuc[0,0,:]),len(nuc[0,0,:])))
if ( len(nuc[0,0,:])%2 == 0 ):
    P1 = np.zeros((len(nuc[:,0,0]),int(0.5*len(nuc[0,0,:])+1)))
    P2 = np.zeros((len(nuc[:,0,0]),int(0.5*len(nuc[0,0,:])+1)))
if ( len(nuc[0,0,:])%2 == 1 ):
    P1 = np.zeros((len(nuc[:,0,0]),int(0.5*(len(nuc[0,0,:]+1)))))
    P2 = np.zeros((len(nuc[:,0,0]),int(0.5*(len(nuc[0,0,:]+1)))))
xmax = freq*MDt*(len(nuc[:,0,0])-1)

# Calculation
## Individual correlations (no averaging)
for j in np.arange(len(nuc[0,0,:])):
    for i in np.arange(len(nuc[0,0,:])):
        u[:,j,i] = nuc[:,0,i]*nuc[:,0,(i+j)%len(nuc[0,0,:])]+\
                   nuc[:,1,i]*nuc[:,1,(i+j)%len(nuc[0,0,:])]+\
                   nuc[:,2,i]*nuc[:,2,(i+j)%len(nuc[0,0,:])]
## Averaging using first and second Legendre polynomials , <P1(u)>(j,t) and <P2(u)>(j,t)
for j in np.arange(len(P1[0,:])):
    for i in np.arange(len(nuc[0,0,:])):
        P1[:,j] = P1[:,j] + u[:,j,i] + u[:,j,(i-j)%len(nuc[0,0,:])]
        P2[:,j] = P2[:,j] + 0.5*(3.*u[:,j,i]**2-1) + 0.5*(3.*u[:,j,(i-j)%len(nuc[0,0,:])]**2-1)
    P1[:,j] = P1[:,j] / (2.*len(nuc[0,0,:]))
    P2[:,j] = P2[:,j] / (2.*len(nuc[0,0,:]))

# Plot
colors = np.linspace(0,1,len(P1[0,:]))
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax.set_xlim(xmin,xmax)
ax.set_ylim(yminp1,ymaxp1)
ax.set_xlabel(r"$t$",fontsize=10)
ax.set_ylabel(r"$C_{1}(t,r)$",fontsize=10)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
for j in np.arange(len(P1[0,:])):
    c = [float(j)/len(P1[0,:]),0.,float(j)/len(P1[0,:])]
    plt.plot(time,P1[:,j],color=c,linewidth=1.,label=r"$C_{1}(t,$"+str(j+1)+"$)$")
#ax.legend(loc='best',fontsize=8)
plt.tight_layout(pad=0.1)
plt.savefig("Orientation_Correlation_P1.eps",format='eps',dpi=300)
plt.clf()

fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax.set_xlim(xmin,xmax)
ax.set_ylim(yminp2,ymaxp2)
ax.set_xlabel(r"$t$",fontsize=10)
ax.set_ylabel(r"$C_{2}(t,r)$",fontsize=10)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
for j in np.arange(len(P1[0,:])):
    c = [float(j)/len(P2[0,:]),0.,float(j)/len(P2[0,:])]
    plt.plot(time,P2[:,j],color=c,linewidth=1.,label=r"$C_{2}(t,$"+str(j+1)+"$)$")
#ax.legend(loc='best',fontsize=8)
plt.tight_layout(pad=0.1)
plt.savefig("Orientation_Correlation_P2.eps",format='eps',dpi=300)
plt.close()
