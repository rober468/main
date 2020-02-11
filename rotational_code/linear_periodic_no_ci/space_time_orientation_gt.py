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
dimer = datafile["Dimer"]
concentration = datafile["Concentration"]
orient = dimer["nucvec"][()]

# Read attributes
box_x = datafile.attrs["simulation_box_length"][0][0]
numdims = datafile.attrs["number_of_dimers"]
MDtime = datafile.attrs["MD_time"][0]
freq = concentration.attrs["io_frequency"][0]

# Read dataset
time = np.arange(0,freq*(len(orient[:,0,0])-1)*MDtime+freq*MDtime,freq*MDtime)
motor_pos = np.arange(0,box_x+box_x/numdims,box_x/numdims)

# Space-time plot
figx = 8.0
figy = 4.0
xmin = 0.
xmax = freq*MDtime*(len(orient[:,0,0])-1)
ymin = 0.
ymax = box_x
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
cax = ax.pcolormesh(time,motor_pos,np.transpose(np.arccos(orient[:,1,:])),cmap=mpl.cm.coolwarm)
cbar = fig.colorbar(cax,orientation='vertical')
#ax.set_xticks(np.arange(0,len(conc2d[0,0,:,i])+1,len(conc2d[0,0,:,i])/int(box_x)*xtix))
#ax.set_yticks(np.arange(0,len(conc2d[0,:,0,i])+1,len(conc2d[0,:,0,i])/int(box_y)*ytix))
#ax.set_xticklabels(np.arange(0,box_x+1,xtix))
#ax.set_yticklabels(np.arange(0,box_y+1,ytix))
plt.savefig("space_time_orientation_gt.pdf",format="pdf",dpi=300)
plt.close()
