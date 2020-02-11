# 2d_concentration_field.py
#
# Calculates the 1D concentration field as a function of x; used when
# initial conditions are not homogeneous on x axis.

# Import module
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import h5py
import matplotlib.animation as animation
import matplotlib.pyplot as plt

# Constants
Rc = 1
xmin = 0.
xmax = 256.
ymin = 0.
ymax = 64.
xtix = 32.
ytix = 8.
datatime = -1
space = 1
figx = 2.0*xmax/xtix
figy = 2.0*ymax/xtix

# Open HDF5 file
datafile = h5py.File("data.h5","r")
# Open groups
dimer = datafile["Dimer"]
concentration = datafile["Concentration"]
# Read dataset
dimer_x = dimer["position_x"]
dimer_y = dimer["position_y"]
# Read attributes
box_x = datafile.attrs["simulation_box_length"][0][0]
box_y = datafile.attrs["simulation_box_length"][0][1]
MDtime = datafile.attrs["MD_time"][0]
freq = concentration.attrs["io_frequency"][0]

# Plot at specific time
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r"$x$",fontsize=10)
ax.set_ylabel(r"$y$",fontsize=10)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.set_xticks(np.arange(0,xmax+1,xtix))
ax.set_yticks(np.arange(0,ymax+1,ytix))
ax.set_xticklabels(np.arange(0,xmax+1,xtix))
ax.set_yticklabels(np.arange(0,ymax+1,ytix))
circle1 = []
circle2 = []
nuc = []
for i in np.arange(len(dimer_x[0,0,:])):
    for j in np.arange(len(dimer_x[0,:,0])):
        if ( j == 0 ):
            circle1.append(plt.Circle((dimer_x[-1,j,i],dimer_y[-1,j,i]),Rc,color='r'))
            ax.add_patch(circle1[i])
        if ( j == 1 ):
            circle2.append(plt.Circle((dimer_x[-1,j,i],dimer_y[-1,j,i]),Rc,color='b'))
            ax.add_patch(circle2[i])
    nuc.append(mpl.lines.Line2D([dimer_x[-1,0,i],dimer_x[-1,1,i]],[dimer_y[-1,0,i],dimer_y[-1,1,i]],color='k',linewidth=2.0))
    ax.add_line(nuc[i])
plt.tight_layout(pad=0.1)
plt.savefig("motors_instantaneous.pdf",format="pdf",dpi=300)
plt.close(fig)

# Animation of data over time
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r"$x$",fontsize=10)
ax.set_ylabel(r"$y$",fontsize=10)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
time_text = ax.text(xmax-20.,ymin+5. , '')
dim_x = dimer_x[::space,:,:]
dim_y = dimer_y[::space,:,:]
circle1 = []
circle2 = []
nuc = []
for i in np.arange(len(dim_x[0,0,:])):
    for j in np.arange(len(dim_x[0,:,0])):
        if ( j == 0 ):
            circle1.append(plt.Circle((dim_x[0,j,i],dim_y[0,j,i]),Rc,color='r'))
        if ( j == 1 ):
            circle2.append(plt.Circle((dim_x[0,j,i],dim_y[0,j,i]),Rc,color='b'))
    nuc.append(mpl.lines.Line2D([dim_x[0,0,i],dim_x[0,1,i]],[dim_y[0,0,i],dim_y[0,1,i]],color='k',linewidth=2.0))
plt.tight_layout(pad=0.1)

def init():
    for i in np.arange(len(dim_x[0,0,:])):
        for j in np.arange(len(dim_x[0,:,0])):
            if ( j == 0 ):
                ax.add_patch(circle1[i])
            if ( j == 1 ):
                ax.add_patch(circle2[i])
        ax.add_line(nuc[i])
    time_text.set_text('')
    return ax,time_text

def animate(k):
    for i in np.arange(len(dim_x[0,0,:])):
        for j in np.arange(len(dim_x[0,:,0])):
            if ( j == 0 ):
                circle1[i].center = (dim_x[k,j,i],dim_y[k,j,i])
            if ( j == 1 ):
                circle2[i].center = (dim_x[k,j,i],dim_y[k,j,i])
        nuc[i].set_data([dim_x[k,0,i],dim_x[k,1,i]],[dim_y[k,0,i],dim_y[k,1,i]])
    time_text.set_text(k*freq*MDtime*space)
    return time_text,ax
    
#anim = animation.FuncAnimation(fig,animate,init_func=init,frames=len(dim_x[:,0,0]),interval=10,blit=True)
#anim.save("motors_dynamics.mp4", fps=int(len(dim_x[:,0,0])*0.5), codec='libx264',extra_args=['-vcodec', 'libx264'],dpi=300)
plt.close(fig)

# Close file
datafile.close()
