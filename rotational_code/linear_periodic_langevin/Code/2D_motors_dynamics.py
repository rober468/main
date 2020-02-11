# 2d_concentration_field.py
#
# Calculates the 1D concentration field as a function of x; used when
# initial conditions are not homogeneous on x axis.

# Import module
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import h5py
import matplotlib.animation as anim
import matplotlib.pyplot as plt

# Constants
Rc = 2.
xmin = 0.
ymin = 0.
ymax = 32.
xtix = 16.
ytix = 16.
datatime = -1
space = 1
figy = ymax/ytix
box_y = ymax
skip = 10

# Open HDF5 file
datafile = h5py.File("data.h5","r")
# Open groups
Angle = datafile["angle"]
Restart = datafile["restart"]
# Read dataset
angle = Angle["angle"]
# Read attributes
Rc = datafile.attrs["C_sphere_radius"][(0)]
Rn = datafile.attrs["N_sphere_radius"][(0)]
d_CN = datafile.attrs["d_CN"][(0)]
sep = datafile.attrs["separation"][(0)]
numrotors = datafile.attrs["number_rotors"][(0)]
timestep = datafile.attrs["timestep"][0]
rotor_CM = datafile.attrs["rotor_CM"][:,:]
freq = Angle.attrs["io_frequency"][0]
cstep = Restart["cstep"][(0)]

# Derived quantities
len_rotor = 2.**(1./6.)*(Rc+Rn) + d_CN
box_x = numrotors*(sep+len_rotor)
xmax = box_x
figx = xmax/xtix
r_CN = 0.5*d_CN
dimer_x = np.zeros((len(angle[:,0]),len(angle[0,:]),2))
dimer_y = np.zeros((len(angle[:,0]),len(angle[0,:]),2))
for i in np.arange(len(dimer_x[:,0,0])):
    for j in np.arange(len(dimer_x[0,:,0])):
        dimer_x[i,j,0] = rotor_CM[0,j] + r_CN * np.cos(angle[i,j])
        dimer_x[i,j,1] = rotor_CM[0,j] - r_CN * np.cos(angle[i,j])
        dimer_y[i,j,0] = rotor_CM[1,j] + r_CN * np.sin(angle[i,j]) + 16.
        dimer_y[i,j,1] = rotor_CM[1,j] - r_CN * np.sin(angle[i,j]) + 16.

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
for i in np.arange(len(dimer_x[0,:,0])):
    for j in np.arange(len(dimer_x[0,0,:])):
        if ( j == 0 ):
            circle1.append(plt.Circle((dimer_x[-1,i,j],dimer_y[-1,i,j]),Rc,color='r'))
            ax.add_patch(circle1[i])
        if ( j == 1 ):
            circle2.append(plt.Circle((dimer_x[-1,i,j],dimer_y[-1,i,j]),Rc,color='b'))
            ax.add_patch(circle2[i])
    nuc.append(mpl.lines.Line2D([dimer_x[-1,i,0],dimer_x[-1,i,1]],[dimer_y[-1,i,0],dimer_y[-1,i,1]],color='k',linewidth=2.0))
    ax.add_line(nuc[i])
plt.tight_layout(pad=0.1)
plt.savefig("motors_instantaneous.pdf",format="pdf",dpi=300)
plt.close(fig)

# Animation of data over time
metadata = dict(title="Langevin,sep=0.5,B<A")
ffwriter = anim.FFMpegWriter(metadata=metadata)
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
time_text = ax.text(xmax-20.,ymin+5. , '')
circle1 = []
circle2 = []
nuc = []
for i in np.arange(len(dimer_x[0,:,0])):
    for j in np.arange(len(dimer_x[0,0,:])):
        if ( j == 0 ):
            circle1.append(plt.Circle((dimer_x[0,i,j],dimer_y[0,i,j]),Rc,color='b'))
        if ( j == 1 ):
            circle2.append(plt.Circle((dimer_x[0,i,j],dimer_y[0,i,j]),Rc,color='r'))
    nuc.append(mpl.lines.Line2D([dimer_x[0,i,0],dimer_x[0,i,1]],[dimer_y[0,i,0],dimer_y[0,i,1]],color='k',linewidth=2.0))


def init():
    for i in np.arange(len(dimer_x[0,:,0])):
        for j in np.arange(len(dimer_x[0,0,:])):
            if ( j == 0 ):
                ax.add_patch(circle1[i])
            if ( j == 1 ):
                ax.add_patch(circle2[i])
        ax.add_line(nuc[i])
    time_text.set_text('')
    return ax,time_text

def animate(k):
    for i in np.arange(0,len(dimer_x[0,:,0]),):
        for j in np.arange(0,len(dimer_x[0,0,:])):
            if ( j == 0 ):
                circle1[i].center = (dimer_x[k*skip,i,j],dimer_y[k*skip,i,j])
            if ( j == 1 ):
                circle2[i].center = (dimer_x[k*skip,i,j],dimer_y[k*skip,i,j])
        nuc[i].set_data([dimer_x[k*skip,i,0],dimer_x[k*skip,i,1]],[dimer_y[k*skip,i,0],dimer_y[k*skip,i,1]])
    time_text.set_text(k*freq*timestep*skip)
    print k*skip
    return time_text,ax
    
animation = anim.FuncAnimation(fig,animate,init_func=init,frames=int((len(dimer_x[:,0,0])-1)/skip+1),interval=5,blit=True)
#anim.save("motors_dynamics.mp4", fps=100*cstep/(freq*int(len(dimer_x[:,0,0])-1)), codec='libx264',extra_args=['-vcodec', 'libx264'],dpi=300)
animation.save("motors_dynamics.mp4", writer=ffwriter,fps=10*cstep/(freq*skip*int(len(dimer_x[:,0,0])-1)), dpi=300)#codec='libx264',extra_args=['-vcodec', 'libx264'],dpi=300)
plt.close(fig)

# Close file
datafile.close()
