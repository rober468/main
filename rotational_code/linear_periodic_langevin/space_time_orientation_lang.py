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
import math as ma

# begin overall plot
figx = 5.
figy = 6.
fig = plt.figure(figsize=(figx,figy))

# Open HDF5 file
datafile_0_lt = h5py.File("k3/sep_3563/eps_B_lt_A/Run_1/data.h5","r")
datafile_0_gt = h5py.File("k3/sep_3563/eps_B_gt_A/Run_3/data.h5","r")
datafile_3_lt = h5py.File("k3/sep_33102/eps_B_lt_A/Run_2/data.h5","r")
datafile_3_gt = h5py.File("k3/sep_33102/eps_B_gt_A/Run_1/data.h5","r")
# Open groups
# 0,lt
Angle_0_lt = datafile_0_lt["angle"]
# 0 ,gt
Angle_0_gt = datafile_0_gt["angle"]
# 3,lt
Angle_3_lt = datafile_3_lt["angle"]
# 3 ,gt
Angle_3_gt = datafile_3_gt["angle"]

# Read attributes
Rc = datafile_0_lt.attrs["C_sphere_radius"][(0)]
Rn = datafile_0_lt.attrs["N_sphere_radius"][(0)]
d_CN = datafile_0_lt.attrs["d_CN"][(0)]
sep_0 = datafile_0_lt.attrs["separation"][(0)]
sep_3 = datafile_3_lt.attrs["separation"][(0)]
numrotors_0 = datafile_0_lt.attrs["number_rotors"][(0)]
numrotors_3 = datafile_3_lt.attrs["number_rotors"][(0)]
timestep = datafile_0_lt.attrs["timestep"][(0)]
freq = Angle_0_lt.attrs["io_frequency"][(0)]

# Read dataset
angle_0_lt = Angle_0_lt["angle"]
angle_0_gt = Angle_0_gt["angle"]
angle_3_lt = Angle_3_lt["angle"]
angle_3_gt = Angle_3_gt["angle"]

# Switch angles
angle_0_lt_mod = np.zeros((len(angle_0_lt[:,0]),len(angle_0_lt[0,:])))
angle_3_lt_mod = np.zeros((len(angle_3_lt[:,0]),len(angle_3_lt[0,:])))
angle_0_gt_mod = np.zeros((len(angle_0_gt[:,0]),len(angle_0_gt[0,:])))
angle_3_gt_mod = np.zeros((len(angle_3_gt[:,0]),len(angle_3_gt[0,:])))
for i in np.arange(len(angle_0_lt[:,0])):
    for j in np.arange(len(angle_0_lt[0,:])):
        if ( angle_0_lt[i,j] > ma.pi ):
            angle_0_lt_mod[i,j] = 2.*ma.pi - angle_0_lt[i,j]
        else:
            angle_0_lt_mod[i,j] = angle_0_lt[i,j]
for i in np.arange(len(angle_3_lt[:,0])):
    for j in np.arange(len(angle_3_lt[0,:])):
        if ( angle_3_lt[i,j] > ma.pi ):
            angle_3_lt_mod[i,j] = 2.*ma.pi - angle_3_lt[i,j]
        else:
            angle_3_lt_mod[i,j] = angle_3_lt[i,j]

for i in np.arange(len(angle_0_gt[:,0])):
    for j in np.arange(len(angle_0_gt[0,:])):
        angle_0_gt_mod[i,j] = angle_0_gt[i,j] - 0.5*ma.pi
        if (angle_0_gt_mod[i,j] < 0. ):
            angle_0_gt_mod[i,j] = angle_0_gt_mod[i,j] + 2.*ma.pi
        if ( angle_0_gt_mod[i,j] > ma.pi ):
            angle_0_gt_mod[i,j] = 2.*ma.pi - angle_0_gt_mod[i,j]
#            angle_0_gt_mod[i,j] = (2.*ma.pi - angle_0_gt[i,j]) + ma.pi
#        else:
#            angle_0_gt_mod[i,j] = angle_0_gt[i,j]
#        angle_0_gt_mod[i,j] = angle_0_gt_mod[i,j] - ma.pi/2.
for i in np.arange(len(angle_3_gt[:,0])):
    for j in np.arange(len(angle_3_gt[0,:])):
        angle_3_gt_mod[i,j] = angle_3_gt[i,j] - 0.5*ma.pi
        if (angle_3_gt_mod[i,j] < 0. ):
            angle_3_gt_mod[i,j] = angle_3_gt_mod[i,j] + 2.*ma.pi
        if ( angle_3_gt_mod[i,j] > ma.pi ):
            angle_3_gt_mod[i,j] = 2.*ma.pi - angle_3_gt_mod[i,j]

# Derived quantities
skip=20
len_rotor = 2.**(1./6.)*(Rc+Rn) + d_CN
sim_box_len_0 = numrotors_0*(sep_0+len_rotor)
sim_box_len_3 = numrotors_3*(sep_3+len_rotor)
print sim_box_len_0,numrotors_0
print sim_box_len_3,numrotors_3
time = np.arange(0,freq*(len(angle_0_lt[:,0])-1)*timestep+freq*timestep,skip*freq*timestep)
motor_ang_0 = np.arange(0,sim_box_len_0+sim_box_len_0/numrotors_0,sim_box_len_0/numrotors_0)
motor_ang_3 = np.arange(0,sim_box_len_3+sim_box_len_3/numrotors_3,sim_box_len_3/numrotors_3)

# Space-time plot
figx = 4.0
figy = 5.0
xmin = 0.
xmax = freq*timestep*(len(angle_0_lt[:,0])-1)
ymin = 0.
ymax = 256.
xtix = xmax/2.
ytix = 256.

fig = plt.figure(figsize=(figx,figy))
gs = mpl.gridspec.GridSpec(4,1,height_ratios=[1,1,1,1])

# 0,lt
ax_0_lt = fig.add_subplot(gs[0])
plt.setp(ax_0_lt.get_xticklabels(),visible=False)
ax_0_lt.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.0f'))
ax_0_lt.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.0f'))
ax_0_lt.set_xlim(xmin,xmax)
ax_0_lt.set_ylim(ymin,ymax)
#ax_0_lt.set_xlabel(r"$t$",fontsize=6)
#ax_0_lt.set_ylabel(r"$x$",fontsize=10)
ax_0_lt.xaxis.set_tick_params(labelsize=8)
ax_0_lt.yaxis.set_tick_params(labelsize=8)
ax_0_lt.set_xticks(np.arange(0,xmax+1,xtix))
ax_0_lt.set_yticks(np.arange(0,ymax+1,ytix))
ax_0_lt.set_xticklabels(np.arange(0,int(xmax+1),int(xtix)))
#ax_0_lt.set_yticklabels(np.arange(0,int(ymax+1),int(ytix)))
ax_0_lt.set_yticklabels(())
cax_0_lt = ax_0_lt.pcolormesh(time,motor_ang_0,np.transpose(angle_0_lt_mod[::skip,:]),cmap=mpl.cm.coolwarm)
ax_0_lt.text(500,200,"a)",fontsize=10,bbox={'facecolor':'white','edgecolor':'none'})
# 3,lt
ax_3_lt = fig.add_subplot(gs[1],sharex=ax_0_lt)
plt.setp(ax_3_lt.get_xticklabels(),visible=False)
#ax_3_lt.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
#ax_3_lt.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax_3_lt.set_xlim(xmin,xmax)
ax_3_lt.set_ylim(ymin,ymax)
#ax_3_lt.set_xlabel(r"$t$",fontsize=6)
#ax_3_lt.set_ylabel(r"$x$",fontsize=10)
ax_3_lt.xaxis.set_tick_params(labelsize=8)
ax_3_lt.yaxis.set_tick_params(labelsize=8)
ax_3_lt.set_yticks(np.arange(0,ymax+1,ytix))
#ax_3_lt.set_yticklabels(np.arange(0,int(ymax+1),int(ytix)))
ax_3_lt.set_yticklabels(())
cax_3_lt = ax_3_lt.pcolormesh(time,motor_ang_3,np.transpose(angle_3_lt_mod[::skip,:]),cmap=mpl.cm.coolwarm)
ax_3_lt.text(500,200,"b)",fontsize=10,bbox={'facecolor':'white','edgecolor':'none'})
# 0,gt
ax_0_gt = fig.add_subplot(gs[2],sharex=ax_0_lt)
plt.setp(ax_0_gt.get_xticklabels(),visible=False)
#ax_0_gt.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
#ax_0_gt.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax_0_gt.set_xlim(xmin,xmax)
ax_0_gt.set_ylim(ymin,ymax)
#ax_0_gt.set_xlabel(r"$t$",fontsize=6)
#ax_0_gt.set_ylabel(r"$x$",fontsize=10)
ax_0_gt.xaxis.set_tick_params(labelsize=8)
ax_0_gt.yaxis.set_tick_params(labelsize=8)
ax_0_gt.set_yticks(np.arange(0,ymax+1,ytix))
#ax_0_gt.set_yticklabels(np.arange(0,int(ymax+1),int(ytix)))
ax_0_gt.set_yticklabels(())
cax_0_gt = ax_0_gt.pcolormesh(time,motor_ang_0,np.transpose(angle_0_gt_mod[::skip,:]),cmap=mpl.cm.coolwarm)
ax_0_gt.text(500,200,"c)",fontsize=10,bbox={'facecolor':'white','edgecolor':'none'})
# 3,gt
ax_3_gt = fig.add_subplot(gs[3],sharex=ax_0_lt)
#ax_3_gt.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
#ax_3_gt.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax_3_gt.set_xlim(xmin,xmax)
ax_3_gt.set_ylim(ymin,ymax)
ax_3_gt.set_xlabel(r"$t$",fontsize=10)
#ax_3_gt.set_ylabel(r"$x$",fontsize=10)
ax_3_gt.xaxis.set_tick_params(labelsize=8)
ax_3_gt.yaxis.set_tick_params(labelsize=8)
ax_3_gt.set_yticks(np.arange(0,ymax+1,ytix))
#ax_3_gt.set_yticklabels(np.arange(0,int(ymax+1),int(ytix)))
ax_3_gt.set_yticklabels(())
cax_3_gt = ax_3_gt.pcolormesh(time,motor_ang_3,np.transpose(angle_3_gt_mod[::skip,:]),cmap=mpl.cm.coolwarm)
ax_3_gt.text(500,200,"d)",fontsize=10,bbox={'facecolor':'white','edgecolor':'none'})

plt.tight_layout(pad=0.1)
plt.subplots_adjust(hspace=0.15)
fig.subplots_adjust(top=0.9)
cbar_ax = fig.add_axes([0.01,0.97,0.945,0.02])
cbar = fig.colorbar(cax_0_lt,cax=cbar_ax,orientation='horizontal')
cbar.ax.tick_params(labelsize=8)
plt.setp(cbar.ax.xaxis.get_ticklabels()[0],visible=False)
plt.savefig("space_time_orientation_lang.eps",format="eps",dpi=300)
plt.close()
