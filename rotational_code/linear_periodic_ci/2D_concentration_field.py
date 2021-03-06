# 2d_concentration_field.py
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

# Constants
figx = 8.0
figy = 4.0
xmin = 0.
xmax = 64.
ymin = 0.
ymax = 64.
xtix = 8.
ytix = 8.
datatime = -1

# Open HDF5 file
datafile = h5py.File("data.h5","r")
# Open groups
concentration = datafile["Concentration"]
# Read dataset
conc2d = np.zeros((len(concentration["2D_concentration_A"].value[:,0,0]),len(concentration["2D_concentration_A"].value[0,:,0]),len(concentration["2D_concentration_A"].value[0,0,:]),2))
conc2dname = ["A","B"]
conc2d[:,:,:,0] = concentration["2D_concentration_A"].value
conc2d[:,:,:,1] = concentration["2D_concentration_B"].value
# Read attributes
box_x = datafile.attrs["simulation_box_length"][0][0]
box_y = datafile.attrs["simulation_box_length"][0][1]
MDtime = datafile.attrs["MD_time"][0]
freq = concentration.attrs["io_frequency"][0]

for i in np.arange(2):
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
    cax = ax.imshow(conc2d[datatime,:,:,i],interpolation='none',cmap=mpl.cm.coolwarm)
    cbar = fig.colorbar(cax,orientation='horizontal')
    ax.set_xticks(np.arange(0,len(conc2d[0,0,:,i])+1,len(conc2d[0,0,:,i])/int(box_x)*xtix))
    ax.set_yticks(np.arange(0,len(conc2d[0,:,0,i])+1,len(conc2d[0,:,0,i])/int(box_y)*ytix))
    ax.set_xticklabels(np.arange(0,box_x+1,xtix))
    ax.set_yticklabels(np.arange(0,box_y+1,ytix))
    plt.savefig("2D_concentration_"+conc2dname[i]+"_field.pdf",format="pdf",dpi=300)
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
    cax = plt.imshow(conc2d[-1,:,:,i],interpolation='none',cmap=mpl.cm.jet)
    cbar = fig.colorbar(cax,orientation='horizontal')
    ax.set_xticks(np.arange(0,len(conc2d[0,0,:,i])+1,len(conc2d[0,0,:,i])/int(box_x)*xtix))
    ax.set_yticks(np.arange(0,len(conc2d[0,:,0,i])+1,len(conc2d[0,:,0,i])/int(box_y)*ytix))
    ax.set_xticklabels(np.arange(0,box_x+1,xtix))
    ax.set_yticklabels(np.arange(0,box_y+1,ytix))
    time_text = ax.text(xmax-20.,ymin+5. , '')

    def init():
        cax.set_data(np.zeros((len(conc2d[0,0,:,i]),len(conc2d[0,:,0,i]))))
        time_text.set_text('')
        return cax,time_text

    def animate(j):
        cax.set_array(conc2d[j,:,:,i])
        time_text.set_text(j*freq*MDtime)
        return cax,time_text
    
    anim = animation.FuncAnimation(fig,animate,init_func=init,frames=len(conc2d[:,0,0,i]),interval=10,blit=True)
    anim.save("2d_concentration_"+conc2dname[i]+"_field.mp4", fps=int(len(conc2d[:,0,0,i])*0.5), codec='libx264',extra_args=['-vcodec', 'libx264'],dpi=300)
    plt.close(fig)

# Close file
datafile.close()
