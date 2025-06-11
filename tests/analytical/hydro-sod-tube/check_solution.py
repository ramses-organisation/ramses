import numpy as np
import scipy.interpolate as interp
import visu_ramses

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
order    = data["data"]["x"].argsort()
x_sim    = data["data"]["x"][order]
rho_sim  = data["data"]["density"][order]
u_sim    = data["data"]["velocity_x"][order]
p_sim    = data["data"]["pressure"][order]

# Read analytical solution, which is given on 1024 points (lvl 10)
data_ana = np.loadtxt('hydro-sod-tube-ana.dat')
x_ana    = data_ana[:,1]
rho_ana  = data_ana[:,3]
u_ana    = data_ana[:,2]
p_ana    = data_ana[:,4]

# Match to closest analytical point
closest_indices = np.array([np.abs(x_ana - x).argmin() for x in x_sim])
x_ana    = x_ana[closest_indices]
rho_ana  = rho_ana[closest_indices]
u_ana    = u_ana[closest_indices]
p_ana    = p_ana[closest_indices]

# calculate the error for density, velocity and pressure
errors = {}
for var, sim, ana in zip(['density', 'velocity', 'pressure'],
                         [rho_sim, u_sim, p_sim],
                         [rho_ana, u_ana, p_ana]):
    abs_error = np.abs(sim - ana)
    rel_error = np.zeros(len(abs_error))
    myfilter = (ana>0)
    rel_error[myfilter] = abs_error[myfilter] / ana[myfilter]
    myfilter = (ana<=0)&(sim>0)
    rel_error[myfilter] = abs_error[myfilter] / sim[myfilter]
    # we take the max error I guess
    errors[var] = np.max(rel_error)


    error = np.abs(sim_y_interp - analytical_y)  # Absolute error
    rel_error = np.abs(error / (analytical_y + 1e-8))  # Relative error (avoid divide by zero)
    
    return np.max(error), np.mean(error), np.max(rel_error), np.mean(rel_error)

# Check results against reference solution
visu_ramses.check_solution(data["data"],'sod-tube')
