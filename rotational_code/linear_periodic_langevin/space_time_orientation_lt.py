# space_time_orientation_lt.py
#
# Calculates the 1D concentration field as a function of x; used when
# initial conditions are not homogeneous on x axis.

# Import module
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
import matplotlib.animation as animation

# Open HDF5 file
datafile = h5py.File("data.h5","r")
# Open groups
Angle = datafile["angle"]

# Read attributes
Rc = datafile.attrs["C_sphere_radius"][(0)]
Rn = datafile.attrs["N_sphere_radius"][(0)]
d_CN = datafile.attrs["d_CN"][(0)]
sep = datafile.attrs["separation"][(0)]
numrotors = datafile.attrs["number_rotors"][(0)]
timestep = datafile.attrs["timestep"][(0)]
freq = Angle.attrs["io_frequency"][(0)]

# Read dataset
angle = Angle["angle"]

# Derived quantities
len_rotor = 2.**(1./6.)*(Rc+Rn) + d_CN
sim_box_len = numrotors*(sep+len_rotor)
time = np.arange(0,freq*(len(angle[:,0])-1)*timestep+freq*timestep,freq*timestep)
motor_ang = np.arange(0,sim_box_len+sim_box_len/numrotors,sim_box_len/numrotors)

# Space-time plot
figx = 4.0
figy = 2.0
xmin = 0.
xmax = freq*timestep*(len(angle[:,0])-1)
ymin = 0.
ymax = sim_box_len
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r"$t$",fontsize=10)
ax.set_ylabel(r"$x$",fontsize=10)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
cax = ax.pcolormesh(time,motor_ang,np.transpose(angle[:,:]),cmap=mpl.cm.coolwarm)
cbar = fig.colorbar(cax,orientation='vertical')
#ax.set_xticks(np.arange(0,len(conc2d[0,0,:,i])+1,len(conc2d[0,0,:,i])/int(box_x)*xtix))
#ax.set_yticks(np.arange(0,len(conc2d[0,:,0,i])+1,len(conc2d[0,:,0,i])/int(box_y)*ytix))
#ax.set_xticklabels(np.arange(0,box_x+1,xtix))
#ax.set_yticklabels(np.arange(0,box_y+1,ytix))
plt.savefig("space_time_orientation_lt.png",format="png",dpi=300)
plt.close()
