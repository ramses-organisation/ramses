import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import csv
import shutil
import visu_ramses

## Default params

fig = plt.figure(figsize=(8, 12))
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)
ax5 = ax1.twinx()

# Load RAMSES output
data = visu_ramses.load_snapshot(2)

order = data["data"]["x"].argsort()
x_sim = data["data"]["x"][order]
amrlev = data["data"]["level"][order]
rho_sim = data["data"]["density"][order]
u_sim = data["data"]["velocity_x"][order]
p_sim = data["data"]["pressure"][order]

# Density
ax1.plot(x_sim, rho_sim, "o", color="black", markerfacecolor="none")
ax1.set_xlabel("Distance (cm)")
ax1.set_ylabel("Density (g/cm3)")
ax5.plot(x_sim, amrlev, color="black", ls="dotted")
ax5.set_ylabel("AMR Level")

# Velocity
ax2.plot(x_sim, u_sim, "o", color="black", markerfacecolor="none")
ax2.set_xlabel("Distance (cm)")
ax2.set_ylabel("Velocity (cm/s)")

# Pressure
ax3.plot(x_sim, p_sim, "o", color="black", markerfacecolor="none")
ax3.set_xlabel("Distance (cm)")
ax3.set_ylabel("Pressure (g/cm/s2)")

# Read analytical solution
data_ana = np.loadtxt("sod-tube-ana.dat")
x_ana = data_ana[:, 1]
rho_ana = data_ana[:, 3]
u_ana = data_ana[:, 2]
p_ana = data_ana[:, 4]
ax1.plot(x_ana, rho_ana, color="red")
ax2.plot(x_ana, u_ana, color="red")
ax3.plot(x_ana, p_ana, color="red")

fig.subplots_adjust(wspace=0.3)
fig.savefig("sod-tube.pdf", bbox_inches="tight")


# Interpolate the analytical solution to the simulation grid
rho_ana_interp = np.interp(x_sim, x_ana, rho_ana)
u_ana_interp = np.interp(x_sim, x_ana, u_ana)
p_ana_interp = np.interp(x_sim, x_ana, p_ana)

tolerance={}
for var, sim, ana in zip(
    ["density", "velocity", "pressure"], [rho_sim, u_sim, p_sim], [rho_ana_interp, u_ana_interp, p_ana_interp]
):
    abs_error = np.abs(sim - ana)
    rel_error = np.zeros(len(abs_error))
    mask = ana != 0
    rel_error[mask] = abs_error[mask] / np.abs(ana[mask])
    rel_error[np.logical_not(mask)] = abs_error[np.logical_not(mask)]

    data["data"][f"{var}_med_error"] = np.median(rel_error)
    data["data"][f"{var}_avg_error"] = np.mean(rel_error)
    tolerance[f"{var}_med_error"] = 1e-8


## Parameter study

# 1. Read error computed with other parameters
parameter_study_errors = {}
with open("sod-tube-parameter-study.csv", "r", encoding="utf-8") as csvfile:
    reader = csv.DictReader(csvfile)          # reads the header row automatically
    for row in reader:
        combo = row["test_name_combi"]
        err   = float(row["error"])           # cast if numeric
        parameter_study_errors[combo] = err

# Overwrite reference solution
overwrite=False
if overwrite:
    shutil.copyfile("sod-tube-parameter-study.csv", "sod-tube-parameter-study-ref.csv")

# 2. Read ref
parameter_study_ref = {}
with open("sod-tube-parameter-study-ref.csv", "r", encoding="utf-8") as csvfile:
    reader = csv.DictReader(csvfile)          # reads the header row automatically
    for row in reader:
        combo = row["test_name_combi"]
        err   = float(row["error"])           # cast if numeric
        parameter_study_ref[combo] = err


# 3. Compare
tolerance_param_study = 1e-8
num_wrong = 0

for key in parameter_study_errors:

    if key not in parameter_study_ref:
        print(f"Error, {key} not in parameter study reference.\n")
        num_wrong += 1
    else:
        rel_diff = (parameter_study_errors[key] - parameter_study_ref[key]) / parameter_study_ref[key]
        ok_key = np.all(rel_diff < tolerance_param_study)
        if not ok_key:
            print(f"Error in parameter combinaison {key}.\n")
            print(f"Now {parameter_study_errors[key]} | Ref {parameter_study_ref[key]}.\n")
            num_wrong += 1

data["data"]["nb_failed_param_combinaison"] = num_wrong

## Check results against reference solution

visu_ramses.check_solution(data["data"], "sod-tube", tolerance=tolerance)
