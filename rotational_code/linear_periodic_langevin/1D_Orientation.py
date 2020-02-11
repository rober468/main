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
import math as ma

figx = 4.
figy = 4.

# k2
in_file = h5py.File("k2/averaged_data.h5","r")
mean_1_sep_0_k2 = in_file["mean_1_sep_0"][:,:]
mean_2_sep_0_k2 = in_file["mean_2_sep_0"][:,:]
mean_1_sep_1_k2 = in_file["mean_1_sep_1"][:,:]
mean_2_sep_1_k2 = in_file["mean_2_sep_1"][:,:]
mean_1_sep_3_k2 = in_file["mean_1_sep_3"][:,:]
mean_2_sep_3_k2 = in_file["mean_2_sep_3"][:,:]
mean_1_sep_6_k2 = in_file["mean_1_sep_6"][:,:]
mean_2_sep_6_k2 = in_file["mean_2_sep_6"][:,:]
# k3
in_file = h5py.File("k3/averaged_data.h5","r")
mean_1_sep_0_k3 = in_file["mean_1_sep_0"][:,:]
mean_2_sep_0_k3 = in_file["mean_2_sep_0"][:,:]
mean_1_sep_1_k3 = in_file["mean_1_sep_1"][:,:]
mean_2_sep_1_k3 = in_file["mean_2_sep_1"][:,:]
mean_1_sep_3_k3 = in_file["mean_1_sep_3"][:,:]
mean_2_sep_3_k3 = in_file["mean_2_sep_3"][:,:]
mean_1_sep_6_k3 = in_file["mean_1_sep_6"][:,:]
mean_2_sep_6_k3 = in_file["mean_2_sep_6"][:,:]
# k4
in_file = h5py.File("k4/averaged_data.h5","r")
in_file_1 = h5py.File("k4/averaged_data_1.h5","r")
mean_1_sep_0_k4_lt = in_file["mean_1_sep_0"][0,:]
mean_2_sep_0_k4_lt = in_file["mean_2_sep_0"][0,:]
mean_1_sep_0_k4_gt = in_file_1["mean_1_sep_0"][1,:]
mean_2_sep_0_k4_gt = in_file_1["mean_2_sep_0"][1,:]
mean_1_sep_1_k4 = in_file["mean_1_sep_1"][:,:]
mean_2_sep_1_k4 = in_file["mean_2_sep_1"][:,:]
mean_1_sep_3_k4 = in_file["mean_1_sep_3"][:,:]
mean_2_sep_3_k4 = in_file["mean_2_sep_3"][:,:]
mean_1_sep_6_k4 = in_file["mean_1_sep_6"][:,:]
mean_2_sep_6_k4 = in_file["mean_2_sep_6"][:,:]

mean_1_sep_0_k4 = np.zeros((2,len(mean_1_sep_0_k4_gt)))
mean_2_sep_0_k4 = np.zeros((2,len(mean_2_sep_0_k4_gt)))
mean_1_sep_0_k4[0,:] = mean_1_sep_0_k4_lt[:]
mean_1_sep_0_k4[1,:] = mean_1_sep_0_k4_gt[:]
mean_2_sep_0_k4[0,:] = mean_2_sep_0_k4_lt[:]
mean_2_sep_0_k4[1,:] = mean_2_sep_0_k4_gt[:]

# Plot P1
for k in np.arange(2):
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
    ax.set_xlabel(r"$\sigma_{sep}$",fontsize=10)
    ax.set_ylabel(r"$<\hat{e}_{i}\cdot{}\hat{e}_{i+j}>$",fontsize=10)
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    sep = [0.3563,1.1768,3.3102,5.569]
    line = np.zeros((4,3))
    uncorr_sep = np.linspace(0,6,10)
    uncorr = np.zeros(len(uncorr_sep))
    line[0,0] = mean_1_sep_0_k2[k,1]
    line[1,0] = mean_1_sep_1_k2[k,1]
    line[2,0] = mean_1_sep_3_k2[k,1]
    line[3,0] = mean_1_sep_6_k2[k,1]
    line[0,1] = mean_1_sep_0_k3[k,1]
    line[1,1] = mean_1_sep_1_k3[k,1]
    line[2,1] = mean_1_sep_3_k3[k,1]
    line[3,1] = mean_1_sep_6_k3[k,1]
    line[0,2] = mean_1_sep_0_k4[k,1]
    line[1,2] = mean_1_sep_1_k4[k,1]
    line[2,2] = mean_1_sep_3_k4[k,1]
    line[3,2] = mean_1_sep_6_k4[k,1]
    plt.plot(sep,line[:,0],marker="^",color="k",linestyle="--",label=r"$10^{-2}$")
    plt.plot(sep,line[:,1],marker="D",color="k",linestyle="--",label=r"$10^{-3}$")
    plt.plot(sep,line[:,2],marker="o",color="k",linestyle="--",label=r"$10^{-4}$")
    plt.plot(uncorr_sep,uncorr,color="k",linestyle=":")
    ax.legend(loc='best',fontsize=10)
    plt.tight_layout(pad=0.1)
    if ( k == 0 ):
        ax.set_xlim(0,6)
        ax.set_ylim(-0.1,0.5)
        plt.savefig("orientation_1_lt.png",format='pdf',dpi=300)
    if ( k == 1 ):
        ax.set_xlim(0,6)
        ax.set_ylim(-0.1,0.5)
        plt.savefig("orientation_1_gt.png",format='pdf',dpi=300)
    plt.clf()

# Plot P1
fig = plt.figure(figsize=(figx,2.*figy))
gs = mpl.gridspec.GridSpec(2,1,height_ratios=[1,1])
for k in np.arange(2):
    if ( k == 0 ):
        ax0 = fig.add_subplot(gs[k])
        ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
        ax0.set_ylabel(r"$<\hat{e}_{i}\cdot{}\hat{e}_{i+j}>$",fontsize=10)
        ax0.yaxis.set_tick_params(labelsize=10)
    if ( k == 0 ):
        ax1 = fig.add_subplot(gs[1],sharex=ax0)
        plt.setp(ax0.get_xticklabels(),visible=False)
        ax1.set_xlabel(r"$\sigma_{sep}$",fontsize=10)
        ax0.set_ylabel(r"$<\hat{e}_{i}\cdot{}\hat{e}_{i+j}>$",fontsize=10)
        ax1.xaxis.set_tick_params(labelsize=10)
        ax1.yaxis.set_tick_params(labelsize=10)
    sep = [0.3563,1.1768,3.3102,5.569]
    line = np.zeros((4,3))
    uncorr_sep = np.linspace(0,6,10)
    uncorr = np.zeros(len(uncorr_sep))
    line[0,0] = mean_1_sep_0_k2[k,1]
    line[1,0] = mean_1_sep_1_k2[k,1]
    line[2,0] = mean_1_sep_3_k2[k,1]
    line[3,0] = mean_1_sep_6_k2[k,1]
    line[0,1] = mean_1_sep_0_k3[k,1]
    line[1,1] = mean_1_sep_1_k3[k,1]
    line[2,1] = mean_1_sep_3_k3[k,1]
    line[3,1] = mean_1_sep_6_k3[k,1]
    line[0,2] = mean_1_sep_0_k4[k,1]
    line[1,2] = mean_1_sep_1_k4[k,1]
    line[2,2] = mean_1_sep_3_k4[k,1]
    line[3,2] = mean_1_sep_6_k4[k,1]
    if ( k == 0 ):
        ax0.plot(sep,line[:,0],marker="^",color="k",linestyle="--",label=r"$10^{-2}$")
        ax0.plot(sep,line[:,1],marker="D",color="k",linestyle="--",label=r"$10^{-3}$")
        ax0.plot(sep,line[:,2],marker="o",color="k",linestyle="--",label=r"$10^{-4}$")
        ax0.plot(uncorr_sep,uncorr,color="k",linestyle=":")
        ax0.set_ylim(-0.1,0.5)
    if ( k == 1 ):
        ax1.plot(sep,line[:,0],marker="^",color="k",linestyle="--",label=r"$10^{-2}$")
        ax1.plot(sep,line[:,1],marker="D",color="k",linestyle="--",label=r"$10^{-3}$")
        ax1.plot(sep,line[:,2],marker="o",color="k",linestyle="--",label=r"$10^{-4}$")
        ax1.plot(uncorr_sep,uncorr,color="k",linestyle=":")
        ax1.set_xlim(0,6)
        ax1.set_ylim(-0.1,0.2)
        plt.setp(ax1.get_yticklabels()[-1],visible=False)
        ax0.legend(loc='upper right',fontsize=10)
plt.tight_layout(pad=0.1)
plt.subplots_adjust(hspace=0.001)
plt.savefig("orientation_1.png",format='png',dpi=300)
plt.clf()

# Plot P2
for k in np.arange(2):
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
    #ax.set_xlim(xmin,xmax)
    #ax.set_ylim(yminp1,ymaxp1)
    ax.set_xlabel(r"$\sigma_{sep}$",fontsize=10)
    ax.set_ylabel(r"$\frac{1}{2}(3<(\hat{e}_{i}\cdot{}\hat{e}_{i+j})-1)>$",fontsize=10)
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    uncorr_sep = np.linspace(0,6,10)
    uncorr = 0.25*np.ones(len(uncorr_sep))
    sep = [0.3563,1.1768,3.3102,5.569]
    line = np.zeros((4,3))
    line[0,0] = mean_2_sep_0_k2[k,1]
    line[1,0] = mean_2_sep_1_k2[k,1]
    line[2,0] = mean_2_sep_3_k2[k,1]
    line[3,0] = mean_2_sep_6_k2[k,1]
    line[0,1] = mean_2_sep_0_k3[k,1]
    line[1,1] = mean_2_sep_1_k3[k,1]
    line[2,1] = mean_2_sep_3_k3[k,1]
    line[3,1] = mean_2_sep_6_k3[k,1]
    line[0,2] = mean_2_sep_0_k4[k,1]
    line[1,2] = mean_2_sep_1_k4[k,1]
    line[2,2] = mean_2_sep_3_k4[k,1]
    line[3,2] = mean_2_sep_6_k4[k,1]
    plt.plot(sep,line[:,0],marker="^",color="k",linestyle="--",label=r"$10^{-2}$")
    plt.plot(sep,line[:,1],marker="D",color="k",linestyle="--",label=r"$10^{-3}$")
    plt.plot(sep,line[:,2],marker="o",color="k",linestyle="--",label=r"$10^{-4}$")
    plt.plot(uncorr_sep,uncorr,color="k",linestyle=":")
    ax.legend(loc='best',fontsize=10)
    plt.tight_layout(pad=0.1)
    if ( k == 0 ):
        ax.set_xlim(0,6)
        ax.set_ylim(0.2,0.7)
        plt.savefig("orientation_2_lt.png",format='png',dpi=300)
    if ( k == 1 ):
        ax.set_xlim(0,6)
        ax.set_ylim(0.2,0.7)
        plt.savefig("orientation_2_gt.png",format='png',dpi=300)
    plt.clf()

# Plot P2
fig = plt.figure(figsize=(figx,1.5*figy))
gs = mpl.gridspec.GridSpec(2,1,height_ratios=[1,1])
for k in np.arange(2):
    if ( k == 0 ):
        ax0 = fig.add_subplot(gs[k])
        ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
        ax0.set_ylabel(r"$S_{2}$",fontsize=10)
        ax0.yaxis.set_tick_params(labelsize=10)
    if ( k == 0 ):
        ax1 = fig.add_subplot(gs[1],sharex=ax0)
        plt.setp(ax0.get_xticklabels(),visible=False)
        ax1.set_xlabel(r"$\sigma_{sep}$",fontsize=10)
        ax1.set_ylabel(r"$S_{2}$",fontsize=10)
        ax1.xaxis.set_tick_params(labelsize=10)
        ax1.yaxis.set_tick_params(labelsize=10)
    sep = [0.3563,1.1768,3.3102,5.569]
    line = np.zeros((4,3))
    uncorr_sep = np.linspace(0,6,10)
    uncorr = 0.25*np.ones(len(uncorr_sep))
    line[0,0] = mean_2_sep_0_k2[k,1]
    line[1,0] = mean_2_sep_1_k2[k,1]
    line[2,0] = mean_2_sep_3_k2[k,1]
    line[3,0] = mean_2_sep_6_k2[k,1]
    line[0,1] = mean_2_sep_0_k3[k,1]
    line[1,1] = mean_2_sep_1_k3[k,1]
    line[2,1] = mean_2_sep_3_k3[k,1]
    line[3,1] = mean_2_sep_6_k3[k,1]
    line[0,2] = mean_2_sep_0_k4[k,1]
    line[1,2] = mean_2_sep_1_k4[k,1]
    line[2,2] = mean_2_sep_3_k4[k,1]
    line[3,2] = mean_2_sep_6_k4[k,1]
    if ( k == 0 ):
        ax0.plot(sep,line[:,0],marker="^",color="k",linestyle="--",label=r"$10^{-2}$")
        ax0.plot(sep,line[:,1],marker="D",color="k",linestyle="--",label=r"$10^{-3}$")
        ax0.plot(sep,line[:,2],marker="o",color="k",linestyle="--",label=r"$10^{-4}$")
        ax0.plot(uncorr_sep,uncorr,color="k",linestyle=":")
        ax0.set_ylim(0.2,0.7)
    if ( k == 1 ):
        ax1.plot(sep,line[:,0],marker="^",color="k",linestyle="--",label=r"$10^{-2}$")
        ax1.plot(sep,line[:,1],marker="D",color="k",linestyle="--",label=r"$10^{-3}$")
        ax1.plot(sep,line[:,2],marker="o",color="k",linestyle="--",label=r"$10^{-4}$")
        ax1.plot(uncorr_sep,uncorr,color="k",linestyle=":")
        ax1.set_xlim(0,6)
        ax1.set_ylim(0.2,0.5)
        plt.setp(ax1.get_yticklabels()[-1],visible=False)
        ax0.legend(loc='upper right',fontsize=10)
plt.tight_layout(pad=0.1)
plt.subplots_adjust(hspace=0.001)
plt.savefig("orientation_2.png",format='png',dpi=300)
plt.clf()

# Plot for P1 and P2
figx = 5.
figy = 5.
fig = plt.figure(figsize=(figx,figy))
gs = mpl.gridspec.GridSpec(2,2)#height_ratios=[1,1,1,1])
# left side
ax00 = fig.add_subplot(gs[0,0])
ax10 = fig.add_subplot(gs[1,0],sharex=ax00)
ax00.set_ylabel(r"$S_{2}(1)$",fontsize=10)
ax10.set_ylabel(r"$S_{1}(1)$",fontsize=10)
ax10.set_xlabel(r"$\ell_{s}$",fontsize=10)
ax00.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax10.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax00.yaxis.set_tick_params(labelsize=10)
ax10.yaxis.set_tick_params(labelsize=10)
ax10.xaxis.set_tick_params(labelsize=10)
# right side
ax01 = fig.add_subplot(gs[0,1],sharey=ax00)
ax11 = fig.add_subplot(gs[1,1],sharex=ax01,sharey=ax10)
plt.setp(ax11.get_yticklabels(),visible=False)
ax11.set_xlabel(r"$\ell_{s}$",fontsize=10)
ax11.xaxis.set_tick_params(labelsize=10)

#sep = [0.3563,1.1768,3.3102,5.569]
sep = [9.846,10.667,12.8,15.059]
line = np.zeros((2,2,4,3))
uncorr_sep = np.linspace(9,16,10)
uncorr_1 = np.zeros(len(uncorr_sep))
uncorr_2 = 0.25*np.ones(len(uncorr_sep))
for k in np.arange(2):
    line[0,k,0,0] = mean_1_sep_0_k2[k,1]
    line[0,k,1,0] = mean_1_sep_1_k2[k,1]
    line[0,k,2,0] = mean_1_sep_3_k2[k,1]
    line[0,k,3,0] = mean_1_sep_6_k2[k,1]
    line[0,k,0,1] = mean_1_sep_0_k3[k,1]
    line[0,k,1,1] = mean_1_sep_1_k3[k,1]
    line[0,k,2,1] = mean_1_sep_3_k3[k,1]
    line[0,k,3,1] = mean_1_sep_6_k3[k,1]
    line[0,k,0,2] = mean_1_sep_0_k4[k,1]
    line[0,k,1,2] = mean_1_sep_1_k4[k,1]
    line[0,k,2,2] = mean_1_sep_3_k4[k,1]
    line[0,k,3,2] = mean_1_sep_6_k4[k,1]
    line[1,k,0,0] = mean_2_sep_0_k2[k,1]
    line[1,k,1,0] = mean_2_sep_1_k2[k,1]
    line[1,k,2,0] = mean_2_sep_3_k2[k,1]
    line[1,k,3,0] = mean_2_sep_6_k2[k,1]
    line[1,k,0,1] = mean_2_sep_0_k3[k,1]
    line[1,k,1,1] = mean_2_sep_1_k3[k,1]
    line[1,k,2,1] = mean_2_sep_3_k3[k,1]
    line[1,k,3,1] = mean_2_sep_6_k3[k,1]
    line[1,k,0,2] = mean_2_sep_0_k4[k,1]
    line[1,k,1,2] = mean_2_sep_1_k4[k,1]
    line[1,k,2,2] = mean_2_sep_3_k4[k,1]
    line[1,k,3,2] = mean_2_sep_6_k4[k,1]
# left side
ax00.plot(sep,line[1,0,:,0],marker="^",color="k",linestyle="--")
ax00.plot(sep,line[1,0,:,1],marker="D",color="k",linestyle="--")
ax00.plot(sep,line[1,0,:,2],marker="o",color="k",linestyle="--")
ax00.plot(uncorr_sep,uncorr_2,color="k",linestyle=":")
ax10.plot(sep,line[0,0,:,0],marker="^",color="k",linestyle="--")
ax10.plot(sep,line[0,0,:,1],marker="D",color="k",linestyle="--")
ax10.plot(sep,line[0,0,:,2],marker="o",color="k",linestyle="--")
ax10.plot(uncorr_sep,uncorr_1,color="k",linestyle=":")
ax00.set_ylim(0.2,0.7)
ax10.set_ylim(-0.1,0.5)
ax10.set_xlim(9.5,15.5)
ax00.set_yticks(np.arange(0.25,0.71,0.15))
ax00.set_yticklabels(np.arange(0.25,0.71,0.15))
ax10.set_yticks(np.arange(-0.05,0.5,0.15))
ax10.set_yticklabels(np.arange(-0.05,0.5,0.15))
ax10.set_xticks(np.arange(10,15,2))
ax10.set_xticklabels(np.arange(10,15,2))
plt.setp(ax00.get_xticklabels(),visible=False)
#plt.setp(ax10.get_xticklabels()[-1],visible=False)
#plt.setp(ax10.get_yticklabels()[-1],visible=False)
# right side
ax01.plot(sep,line[1,1,:,0],marker="^",color="k",linestyle="--",label=r"$\ell_{c}=1.53$")
ax01.plot(sep,line[1,1,:,1],marker="D",color="k",linestyle="--",label=r"$\ell_{c}=4.83$")
ax01.plot(sep,line[1,1,:,2],marker="o",color="k",linestyle="--",label=r"$\ell_{c}=15.27$")
ax01.plot(uncorr_sep,uncorr_2,color="k",linestyle=":")
ax11.plot(sep,line[0,1,:,0],marker="^",color="k",linestyle="--")
ax11.plot(sep,line[0,1,:,1],marker="D",color="k",linestyle="--")
ax11.plot(sep,line[0,1,:,2],marker="o",color="k",linestyle="--")
ax11.plot(uncorr_sep,uncorr_1,color="k",linestyle=":")
ax11.set_xlim(9.5,15.5)
ax11.set_xticks(np.arange(10,15,2))
ax11.set_xticklabels(np.arange(10,15,2))
plt.setp(ax01.get_xticklabels(),visible=False)
plt.setp(ax01.get_yticklabels(),visible=False)
plt.setp(ax11.get_yticklabels(),visible=False)
ax01.legend(loc='upper right',fontsize=10)
plt.tight_layout(pad=0.1)
plt.subplots_adjust(wspace=0.001,hspace=0.001)
ax00.text(10,0.66,"a)",fontsize=10)
ax01.text(10,0.66,"b)",fontsize=10)
ax10.text(10,0.46,"c)",fontsize=10)
ax11.text(10,0.46,"d)",fontsize=10)

plt.savefig("orientation_1_2.png",format='png',dpi=300)
plt.clf()

# j-th motor away
figx = 4.
figy = 2.5
fig = plt.figure(figsize=(figx,figy))
gs = mpl.gridspec.GridSpec(1,2)#height_ratios=[1,1,1,1])
# left side
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1],sharey=ax0)
ax0.set_ylabel(r"$S_{2}(j)$",fontsize=10)
ax0.set_xlabel(r"$j$",fontsize=10)
ax1.set_xlabel(r"$j$",fontsize=10)
ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax0.yaxis.set_tick_params(labelsize=10)
ax0.xaxis.set_tick_params(labelsize=10)
ax1.xaxis.set_tick_params(labelsize=10)

jth = np.arange(1,8)
line = np.zeros((2,3,7))
for k in np.arange(2):
    line[k,0,0] = mean_2_sep_0_k2[k,1]
    line[k,0,1] = mean_2_sep_0_k2[k,2]
    line[k,0,2] = mean_2_sep_0_k2[k,3]
    line[k,0,3] = mean_2_sep_0_k2[k,4]
    line[k,0,4] = mean_2_sep_0_k2[k,5]
    line[k,0,5] = mean_2_sep_0_k2[k,6]
    line[k,0,6] = mean_2_sep_0_k2[k,7]
    line[k,1,0] = mean_2_sep_0_k3[k,1]
    line[k,1,1] = mean_2_sep_0_k3[k,2]
    line[k,1,2] = mean_2_sep_0_k3[k,3]
    line[k,1,3] = mean_2_sep_0_k3[k,4]
    line[k,1,4] = mean_2_sep_0_k3[k,5]
    line[k,1,5] = mean_2_sep_0_k3[k,6]
    line[k,1,6] = mean_2_sep_0_k3[k,7]
    line[k,2,0] = mean_2_sep_0_k4[k,1]
    line[k,2,1] = mean_2_sep_0_k4[k,2]
    line[k,2,2] = mean_2_sep_0_k4[k,3]
    line[k,2,3] = mean_2_sep_0_k4[k,4]
    line[k,2,4] = mean_2_sep_0_k4[k,5]
    line[k,2,5] = mean_2_sep_0_k4[k,6]
    line[k,2,6] = mean_2_sep_0_k4[k,7]
uncorr_j = np.linspace(0,8,10)
uncorr_2 = 0.25*np.ones(len(uncorr_j))
ax0.plot(uncorr_j,uncorr_2,color="k",linestyle=":")
ax1.plot(uncorr_j,uncorr_2,color="k",linestyle=":")
# left side
ax0.plot(jth,line[0,0,:],marker="^",color="k",linestyle="--")
ax0.plot(jth,line[0,1,:],marker="D",color="k",linestyle="--")
ax0.plot(jth,line[0,2,:],marker="o",color="k",linestyle="--")
ax0.set_ylim(0.2,0.75)
ax0.set_xlim(0,8)
ax0.set_xticks(np.arange(1,8,2))
ax0.set_xticklabels(np.arange(1,8,2))
#plt.setp(ax0.get_xticklabels()[-1],visible=False)
ax0.set_yticks(np.arange(0.2,0.76,0.15))
ax0.set_yticklabels(np.arange(0.2,0.76,0.15))
#ax0.set_xticks(np.arange(1,6,2))
#ax0.set_xticklabels(np.arange(1,6,2))
# right side
ax1.plot(jth,line[1,0,:],marker="^",color="k",linestyle="--",label=r"$\ell_{c}=1.53$")
ax1.plot(jth,line[1,1,:],marker="D",color="k",linestyle="--",label=r"$\ell_{c}=4.83$")
ax1.plot(jth,line[1,2,:],marker="o",color="k",linestyle="--",label=r"$\ell_{c}=15.27$")
ax1.set_ylim(0.2,0.75)
ax1.set_xlim(0,8)
ax1.set_xticks(np.arange(1,8,2))
ax1.set_xticklabels(np.arange(1,8,2))
#plt.setp(ax1.get_xticklabels()[-1],visible=False)
ax1.set_yticks(np.arange(0.2,0.76,0.15))
ax1.set_yticklabels(np.arange(0.2,0.76,0.15))
plt.setp(ax1.get_yticklabels(),visible=False)

ax1.legend(loc='upper right',fontsize=10)
plt.tight_layout(pad=0.1)
plt.subplots_adjust(wspace=0.001)
ax0.text(0.5,0.71,"a)",fontsize=10)
ax1.text(0.5,0.71,"b)",fontsize=10)

plt.savefig("s_2_vs_j_sep_0.eps",format='eps',dpi=600)
plt.clf()

# j-th motor away
figx = 4.
figy = 2.5
fig = plt.figure(figsize=(figx,figy))
gs = mpl.gridspec.GridSpec(1,2)#height_ratios=[1,1,1,1])
# left side
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1],sharey=ax0)
ax0.set_ylabel(r"$S_{1}(j)$",fontsize=10)
ax0.set_xlabel(r"$j$",fontsize=10)
ax1.set_xlabel(r"$j$",fontsize=10)
ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
ax0.yaxis.set_tick_params(labelsize=10)
ax0.xaxis.set_tick_params(labelsize=10)
ax1.xaxis.set_tick_params(labelsize=10)

jth = np.arange(1,8,3)
line = np.zeros((2,3,3))
for k in np.arange(2):
    line[k,0,0] = mean_1_sep_0_k2[k,1]
    line[k,0,1] = mean_1_sep_0_k2[k,4]
    line[k,0,2] = mean_1_sep_0_k2[k,7]
    line[k,1,0] = mean_1_sep_0_k3[k,1]
    line[k,1,1] = mean_1_sep_0_k3[k,4]
    line[k,1,2] = mean_1_sep_0_k3[k,7]
    line[k,2,0] = mean_1_sep_0_k4[k,1]
    line[k,2,1] = mean_1_sep_0_k4[k,4]
    line[k,2,2] = mean_1_sep_0_k4[k,7]
uncorr_j = np.linspace(0,8,10)
uncorr_1 = np.ones(len(uncorr_j))
ax0.plot(uncorr_j,uncorr_1,color="k",linestyle=":")
ax1.plot(uncorr_j,uncorr_1,color="k",linestyle=":")
# left side
ax0.plot(jth,line[0,0,:],marker="^",color="k",linestyle="--")
ax0.plot(jth,line[0,1,:],marker="D",color="k",linestyle="--")
ax0.plot(jth,line[0,2,:],marker="o",color="k",linestyle="--")
#ax0.set_ylim(0.2,0.7)
#ax0.set_xlim(0,8)
#ax0.set_xticks(np.arange(1,8,3))
#ax0.set_xticklabels(np.arange(1,8,3))
#plt.setp(ax0.get_xticklabels()[-1],visible=False)
#ax0.set_yticks(np.arange(0.2,0.71,0.15))
#ax0.set_yticklabels(np.arange(0.2,0.71,0.15))
#ax0.set_xticks(np.arange(1,6,2))
#ax0.set_xticklabels(np.arange(1,6,2))
# right side
ax1.plot(jth,line[1,0,:],marker="^",color="k",linestyle="--",label=r"$\ell_{c}=1.53$")
ax1.plot(jth,line[1,1,:],marker="D",color="k",linestyle="--",label=r"$\ell_{c}=4.83$")
ax1.plot(jth,line[1,2,:],marker="o",color="k",linestyle="--",label=r"$\ell_{c}=15.27$")
#ax1.set_ylim(0.2,0.7)
#ax1.set_xlim(0,8)
#ax1.set_xticks(np.arange(1,8,3))
#ax1.set_xticklabels(np.arange(1,8,3))
#plt.setp(ax1.get_xticklabels()[-1],visible=False)
#ax1.set_yticks(np.arange(0.2,0.71,0.15))
#ax1.set_yticklabels(np.arange(0.2,0.71,0.15))
#plt.setp(ax1.get_yticklabels(),visible=False)

ax1.legend(loc='upper right',fontsize=10)
plt.tight_layout(pad=0.1)
plt.subplots_adjust(wspace=0.001)
plt.savefig("s_1_vs_j_sep_0.png",format='png',dpi=600)
plt.clf()
