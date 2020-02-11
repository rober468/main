# subjobs_list.py
#
# Make a list of subjobs for gnu-parallel
#
#
#
# inport libraries
import numpy as np

numrun = 50

# input parameters
sep = ["3563","11768","33102","55690"]
eps = ["eps_B_lt_A","eps_B_gt_A"]
runs = np.arange(1,numrun+1,1)

# open file
file = open("subjobs.txt","w")

# write subjobs to file
for j in np.arange(len(sep)):
    for k in np.arange(len(eps)):
        for l in np.arange(len(runs)):
            file.write("cd sep_"+sep[j]+"/"+eps[k]+"/Run_"+str(runs[l])+"/ ; ./main_langevin_rotors\n")

# close file
file.close()
