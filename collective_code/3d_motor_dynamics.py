from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np
import h5py as h5py
import matplotlib.image as mgimg

# Constants
figx = 4.0
figy = 4.0
figz = 2.0
xmin = 0.
xmax = 64.
ymin = 0.
ymax = 64.
zmin = 0.
zmax = 32.
xtix = 8.
ytix = 8.
ztix = 4.
datatime = -1
int = 1./45
numdims = 1
Rc = 2.
Rn = 2.
space = 100

# Open HDF5 file
datafile = h5py.File("data.h5","r")
# Open groups
dimer = datafile["Dimer"]
# Read dataset
pos_x = dimer["position_x"].value
pos_y = dimer["position_y"].value
pos_z = dimer["position_z"].value
# Read attributes
box_x = datafile.attrs["simulation_box_length"][0][0]
box_y = datafile.attrs["simulation_box_length"][0][1]
box_z = datafile.attrs["simulation_box_length"][0][2]
MDtime = datafile.attrs["MD_time"][0]
freq = dimer.attrs["io_frequency"][0]

# Make spheres
u = np.arange(0,2*np.pi+int,int)
v = np.arange(0,np.pi+int,int)
x_c = np.zeros((len(pos_x[0,0,:]),len(u),len(v)))
y_c = np.zeros((len(pos_y[0,0,:]),len(u),len(v)))
z_c = np.zeros((len(pos_z[0,0,:]),len(u),len(v)))
x_n = np.zeros((len(pos_x[0,1,:]),len(u),len(v)))
y_n = np.zeros((len(pos_y[0,1,:]),len(u),len(v)))
z_n = np.zeros((len(pos_z[0,1,:]),len(u),len(v)))

for i in np.arange(0,len(pos_x[:,0,0]),space):
    for j in np.arange(0,numdims):
        x_c[j,:,:] = pos_x[i,0,j] + Rc * np.outer(np.cos(u),np.sin(v))
        y_c[j,:,:] = pos_y[i,0,j] + Rc * np.outer(np.sin(u),np.sin(v))
        z_c[j,:,:] = pos_z[i,0,j] + Rc * np.outer(np.ones(np.size(u)),np.cos(v))
        x_n[j,:,:] = pos_x[i,1,j] + Rn * np.outer(np.cos(u),np.sin(v))
        y_n[j,:,:] = pos_y[i,1,j] + Rn * np.outer(np.sin(u),np.sin(v))
        z_n[j,:,:] = pos_z[i,1,j] + Rn * np.outer(np.ones(np.size(u)),np.cos(v))

    # Plot surfaces
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111,projection="3d")
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_zlim(zmin,zmax)
    for j in np.arange(numdims):
        ax.plot_surface(x_c[j,:,:],y_c[j,:,:],z_c[j,:,:],color="r",linewidth=0,edgecolors="none")
        ax.plot_surface(x_n[j,:,:],y_n[j,:,:],z_n[j,:,:],color="b",linewidth=0,edgecolors="none")
        ax.plot_surface(0,y_c[j,:,:],z_c[j,:,:],color="grey",linewidth=0,edgecolors="none")
        ax.plot_surface(0,y_n[j,:,:],z_n[j,:,:],color="grey",linewidth=0,edgecolors="none")
        ax.plot_surface(x_c[j,:,:],ymax,z_c[j,:,:],color="grey",linewidth=0,edgecolors="none")
        ax.plot_surface(x_n[j,:,:],ymax,z_n[j,:,:],color="grey",linewidth=0,edgecolors="none")
        ax.plot_surface(x_c[j,:,:],y_c[j,:,:],0,color="grey",linewidth=0,edgecolors="none")
        ax.plot_surface(x_n[j,:,:],y_n[j,:,:],0,color="grey",linewidth=0,edgecolors="none")
    plt.savefig("3d_motor_dynamics_"+str(i/space)+".png",format="png",dpi=300)
    plt.close()

images = []

fig = plt.figure()

for frame in np.arange(0,len(pos_x[:,0,0]),space):
    fname = "3d_motor_dynamics_"+str(frame/space)+".png" %frame
    img = mgimg.imread(fname)
    imgplot = plt.imshow(img)
    images.append([imgplot])

metadata = dict(title="chemical_patterns",artist="Robertson_Kapral",genre="video")
ffwriter = anim.FFMpegWriter(metadata=metadata)


animation = anim.ArtistAnimation(fig,images,interval=100,blit=True)
animation.save("3d_motor_dynamics.mp4", writer=ffwriter,dpi=300)#fps=int(len(dim_x_lt[:,0,0]))/5,dpi=300)
