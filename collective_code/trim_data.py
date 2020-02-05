# trim_data.py

import h5py as h5py


# Open original file
file = h5py.File("data.h5","r")

# Open new file
newfile = h5py.File("newdata.h5","w")

# Open groups
conc = file["Concentration"]
dimer = file["Dimer"]
solv = file["Solvent"]
sys = file["System"]

# Create groups in new file
concnew = newfile.create_group("Concentration")
dimernew = newfile.create_group("Dimer")
solvnew = newfile.create_group("Solvent")
sysnew = newfile.create_group("System")

# Read attributes
MD_time = file.attrs["MD_time"]
MPC_time = file.attrs["MPC_time"]
number_of_dimers = file.attrs["number_of_dimers"]
simulation_box_length = file.attrs["simulation_box_length"]
solvent_temperature = file.attrs["solvent_temperature"]
total_particles = file.attrs["total_#_particles"]

# Write attributes to new file
newfile.attrs["MD_time"] = file.attrs["MD_time"]
newfile.attrs["MPC_time"] = file.attrs["MPC_time"]
newfile.attrs["number_of_dimers"] = file.attrs["number_of_dimers"]
newfile.attrs["simulation_box_length"] = file.attrs["simulation_box_length"]
newfile.attrs["solvent_temperature"] = file.attrs["solvent_temperature"]
newfile.attrs["total_#_particles"] = file.attrs["total_#_particles"]

# Read attributes and relevant data for concentration
conc_io_frequency = conc.attrs["io_frequency"]
conc_A = conc["concentration_A"]
conc_B = conc["concentration_B"]
total_A = conc["total_A"]
total_B = conc["total_B"]

# Write attributes and relevant data to new file for concentration
concnew.attrs["io_frequency"] = conc.attrs["io_frequency"]
conc_A_new = concnew.create_dataset("concentration_A",data=conc_A)
conc_B_new = concnew.create_dataset("concentration_B",data=conc_B)
total_A_new = concnew.create_dataset("total_A",data=total_A)
total_B_new = concnew.create_dataset("total_B",data=total_B)

# Read attributes and data for dimer
dimer_io_frequency = dimer.attrs["io_frequency"]
nucvec = dimer["nucvec"]
position_x = dimer["position_x"]
position_y = dimer["position_y"]
position_z = dimer["position_z"]
velocity_x = dimer["velocity_x"]
velocity_y = dimer["velocity_y"]
velocity_z = dimer["velocity_z"]

# Write attributes and data to new file for dimer
dimernew.attrs["io_frequency"] = dimer.attrs["io_frequency"]
nucvec_new = dimernew.create_dataset("nucvec",data=nucvec)
position_x_new = dimernew.create_dataset("position_x",data=position_x)
position_y_new = dimernew.create_dataset("position_y",data=position_y)
position_z_new = dimernew.create_dataset("position_z",data=position_z)
velocity_x_new = dimernew.create_dataset("velocity_x",data=velocity_x)
velocity_y_new = dimernew.create_dataset("velocity_y",data=velocity_y)
velocity_z_new = dimernew.create_dataset("velocity_z",data=velocity_z)

# Read attributes and data for solvent
solv_io_frequency = solv.attrs["io_frequency"]
flag = solv["flag"]
position_x_sol = solv["position_x"]
position_y_sol = solv["position_y"]
position_z_sol = solv["position_z"]
velocity_x_sol = solv["velocity_x"]
velocity_y_sol = solv["velocity_y"]
velocity_z_sol = solv["velocity_z"]

# Write attributes and data to new file for solvent
solvnew.attrs["io_frequency"] = solv.attrs["io_frequency"]
flag_new = solvnew.create_dataset("flag",data=flag)
position_x_new = solvnew.create_dataset("position_x",data=position_x_sol)
position_y_new = solvnew.create_dataset("position_y",data=position_y_sol)
position_z_new = solvnew.create_dataset("position_z",data=position_z_sol)
velocity_x_new = solvnew.create_dataset("velocity_x",data=velocity_x_sol)
velocity_y_new = solvnew.create_dataset("velocity_y",data=velocity_y_sol)
velocity_z_new = solvnew.create_dataset("velocity_z",data=velocity_z_sol)

# Read attributes and data for system
sys_io_frequency = sys.attrs["io_frequency"]
dimer_temp = sys["dimer_temperature"]
pxsol = sys["pxsol"]
pysol = sys["pysol"]
pzsol = sys["pzsol"]
solv_temp = sys["solvent_temperature"]
total_energy = sys["total_energy"]

# Write attributes and data to new file for solvent
sysnew.attrs["io_frequency"] = sys.attrs["io_frequency"]
dimer_temp_new = sysnew.create_dataset("dimer_temperature",data=dimer_temp)
pxsol_new = sysnew.create_dataset("pxsol",data=pxsol)
pysol_new = sysnew.create_dataset("pysol",data=pysol)
pzsol_new = sysnew.create_dataset("pzsol",data=pzsol)
solv_temp_new = sysnew.create_dataset("solvent_temperature",data=solv_temp)
total_energy_new = sysnew.create_dataset("total_energy",data=total_energy)

file.close()
newfile.close()
