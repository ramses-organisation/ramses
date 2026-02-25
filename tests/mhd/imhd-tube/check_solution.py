import numpy as np
import scipy.interpolate as interp
import visu_ramses

def check_solution():
    # Load RAMSES output
    data = visu_ramses.load_snapshot(2)

    order = data["data"]["x"].argsort()
    x_sim = data["data"]["x"][order] - 1.5
    amrlev = data["data"]["level"][order]
    rho_sim = data["data"]["density"][order]
    u_sim = data["data"]["velocity_x"][order]
    p_sim = data["data"]["pressure"][order]
    By_sim     = 0.5*(data["data"]["B_y_left"][order] + data["data"]["B_y_right"][order])

    # Read analytical solution
    data_ana = np.loadtxt("imhd-tube-ana.dat")
    x_ana    = data_ana[:,0]
    rho_ana  = data_ana[:,1]
    u_ana    = data_ana[:,2]
    p_ana    = data_ana[:,8]
    By_ana   = data_ana[:,6]

    errors = {}

    # Interpolate the analytical solution to the simulation grid
    rho_ana_interp = np.interp(x_sim, x_ana, rho_ana)
    u_ana_interp = np.interp(x_sim, x_ana, u_ana)
    p_ana_interp = np.interp(x_sim, x_ana, p_ana)
    By_ana_interp = np.interp(x_sim, x_ana, By_ana)

    for var, sim, ana in zip(
        ["density", "velocity", "pressure", "mag_field"],
        [rho_sim, u_sim, p_sim, By_sim],
        [rho_ana_interp, u_ana_interp, p_ana_interp, By_ana_interp]
   ):
        abs_error = np.abs(sim - ana)
        rel_error = np.zeros(len(abs_error))
        mask = ana != 0
        rel_error[mask] = abs_error[mask] / np.abs(ana[mask])
        rel_error[np.logical_not(mask)] = abs_error[np.logical_not(mask)]

        #errors[f"{var}_med_error"] = np.median(rel_error)
        errors[f"{var}_avg_error"] = np.mean(rel_error)

    return errors["velocity_avg_error"]
