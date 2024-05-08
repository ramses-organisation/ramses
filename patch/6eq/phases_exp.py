# String searching library
import re
# Standard libraries
import os
import numpy as np
import pandas as pd
from functools import partial

from scipy import constants
from scipy.io import FortranFile
from scipy.optimize import newton

# Declaring functions
# Newton iteration functions to find the saturation pressure/corresponding temperature
def psat_eq2(P_sat,T,C_vl,C_vg,C_pl,C_pg,pinf_l,pinf_g,ql,qg,qql,qqg):
    A = (C_pl - C_pg + qqg - qql)/(C_pg - C_vg)
    B = (ql - qg)/(C_pg - C_vg)
    C = (C_pg - C_pl)/(C_pg - C_vg)
    D = (C_pl - C_vl)/(C_pg - C_vg)
    return np.log(P_sat + pinf_g) - A - B/T - C*np.log(T) - D*np.log(P_sat+pinf_l)
def psat_eq3(T,P_sat,C_vl,C_vg,C_pl,C_pg,pinf_l,pinf_g,ql,qg,qql,qqg):
    A = (C_pl - C_pg + qqg - qql)/(C_pg - C_vg)
    B = (ql - qg)/(C_pg - C_vg)
    C = (C_pg - C_pl)/(C_pg - C_vg)
    D = (C_pl - C_vl)/(C_pg - C_vg)
    return np.log(P_sat + pinf_g) - A - B/T - C*np.log(T) - D*np.log(P_sat+pinf_l)

# Functions for specific volume
def v_l(p_sat,T,C_pl,C_vl,pinf_l):
    return (C_pl - C_vl)*T/ (p_sat + pinf_l)
def v_g(p_sat,T,C_pg,C_vg,pinf_g):
    return (C_pg - C_vg)*T/ (p_sat + pinf_g)

# Reverse eos calculation: Given P --> Find e
def e_sg(P, T, pinf, q, rho, smallgamma):
    return ((P + smallgamma*pinf)/(smallgamma-1) + rho*q)

# Setting code units by choosing G=1
t0=1e-3
l0=1
r0=constants.G/t0**2
p0=r0*l0**2/t0**2
v0= l0/t0


# Experimental reference values
T0, T1, T2, T3 = 298, 473, 439, 588
psat0, psat1, psat2, psat3 = 3166, 15.551e5, 7.152e5, 105.3e5
hl0, hl1 = 104.74e3, 851.6e3
hg0, hg1 = 2473.42e3, 2733.669e3
vg0, vg1 = 42.41, 0.124
vl0, vl1 = 1.0756e-3, 1.4267e-3

# Trying to reproduce the constants given in the paper
Cpexp_l = (hl1 - hl0)/(T1-T0)
Cpexp_g = (hg1 - hg0)/(T1-T0)
pinf_expl = (vl0*T3*psat2 - vl1*T2*psat3)/(vl1*T2 - vl0*T3)
pinf_expg = (vg0*T1*psat0 - vg1*T0*psat1)/(vg1*T0 - vg0*T1)
Cvexp_l = Cpexp_l - vl0/T2 * (psat2 + pinf_expl)
Cvexp_g = Cpexp_g - vg0/T0 * (psat0 + pinf_expg)
gamma_expl = Cpexp_l/Cvexp_l
gamma_expg = Cpexp_g/Cvexp_g
ql_exp = hl0 - Cpexp_l*T0
qg_exp = hg0 - Cpexp_g*T0

# Print cavitation case parameters in code units
print(" t0=%.3E \n l0=%.3E \n r0=%.3E \n p0=%.3E \n v0=%.3E \n"%(t0,l0,r0,p0,v0))
print(" 1 atm in code units: %.3E"%(1e5/p0))
print("\n SG EOS parameters in code units: \n pinf_l=%E \n pinf_g=%E \n q_l=%E \n qg=%E \n"%(pinf_expl/p0,
                                                                                               pinf_expg/p0,
                                                                                               ql_exp*r0/p0,
                                                                                               qg_exp*r0/p0))
# Create empty dataframe
cols = ['P_sat','d_l','d_v','e_l','e_v']
results = pd.DataFrame(columns=cols)

# Calculate saturation pressures and corresponding vapor/liquid densities for linearly spaced temperatures
T_int = np.linspace(298,647,10000)
temp_guess = 1000
qql = 0.
qqg_exp = -23.95e3 # Free parameter adjusted for the curve fit
data = []
for T_temp in T_int:
    psat_temp = partial(psat_eq2,
                        T=T_temp,
                        C_vl=Cvexp_l,
                        C_vg=Cvexp_g,
                        C_pl=Cpexp_l,
                        C_pg=Cpexp_g,
                        pinf_l=pinf_expl,
                        pinf_g=pinf_expg,
                        ql=ql_exp,
                        qg=qg_exp,
                        qql=qql,
                        qqg=qqg_exp)
    # print("Temperature:%.1f ; Guess: %.2f"%(T_temp, temp_guess))
    temp_guess = newton(psat_temp, temp_guess,disp=False)
    # print("Results: %.2f"%(temp_guess))
    # Append (P_sat, rho_l, rho_g) to the array
    data.append([temp_guess,
                 1/v_l(temp_guess,T_temp,Cpexp_l,Cvexp_l,pinf_expl),
                 1/v_g(temp_guess,T_temp,Cpexp_g,Cvexp_g,pinf_expg),
                 e_sg(temp_guess,T_temp, pinf_expl, ql_exp, 1/v_l(temp_guess,T_temp,Cpexp_l,Cvexp_l,pinf_expl), gamma_expl),
                 e_sg(temp_guess,T_temp, pinf_expg, qg_exp, 1/v_g(temp_guess,T_temp,Cpexp_g,Cvexp_g,pinf_expg), gamma_expg)])

# Add (P_{sat,0}, d_{l,0}, d_{v,0}, e_{l,0}, e_{v,0}) = (0,d_{l,1},0, e_{l,1},0)
init = [[0.0, data[1][1], 0.0, data[3][1], 0.0]]
init.extend(data)
data = init
# Add the critical point to the end
T_crit,p_crit, d_crit = 647.14, 220.64e5, 1/3.01e-3
data.append([p_crit,
             d_crit,
             d_crit,
             e_sg(p_crit,T_crit, pinf_expl, ql_exp, d_crit, gamma_expl),
             e_sg(p_crit,T_crit, pinf_expg, qg_exp, d_crit, gamma_expg)])
T_int = np.append(T_int,T_crit)
T_int = np.insert(T_int, 0, 0)

# Fill dataframe with calculated values
for i in range(len(cols)):
    results[cols[i]] = np.asarray(data)[:,i]


# Calculate the bins from previously calculated pressures
ser, bins = pd.cut(results['P_sat'], bins=9000, retbins=True, labels=False)

# New dataframe
results2 = pd.DataFrame(columns=cols)

# Calculate saturation pressures and corresponding vapor/liquid densities with pressures binned
temp_guess = 300
T_int2 = []
data2 = []
for p in bins[1:]:
    psat_temp = partial(psat_eq3,
                        P_sat=p,
                        C_vl=Cvexp_l,
                        C_vg=Cvexp_g,
                        C_pl=Cpexp_l,
                        C_pg=Cpexp_g,
                        pinf_l=pinf_expl,
                        pinf_g=pinf_expg,
                        ql=ql_exp,
                        qg=qg_exp,
                        qql=qql,
                        qqg=qqg_exp)
    # print("Pressure:%.1f ; Guess: %.2f"%(p, temp_guess))
    temp_guess = newton(psat_temp, temp_guess,disp=False)
    T_int2.append(temp_guess)
    # print("Results: %.2f"%(temp_guess))
    # Append (P_sat, rho_l, rho_g) to the array
    data2.append([p/p0,
                 (1/v_l(p,temp_guess,Cpexp_l,Cvexp_l,pinf_expl))/r0,
                 (1/v_g(p,temp_guess,Cpexp_g,Cvexp_g,pinf_expg))/r0,
                 e_sg(p,temp_guess, pinf_expl, ql_exp, 1/v_l(p,temp_guess,Cpexp_l,Cvexp_l,pinf_expl), gamma_expl)/p0,
                 e_sg(p,temp_guess, pinf_expg, qg_exp, 1/v_g(p,temp_guess,Cpexp_g,Cvexp_g,pinf_expg), gamma_expg)/p0])

# Remove the last element because it has the wrong values for P_crit, d_crit
data2.pop(-1)
# Add the (P_{sat,0}, d_{l,0}, d_{v,0}, e_{l,0}, e_{v,0}) = (0,d_{l,1},0, e_{l,1},0)
init = [[0.0, data2[0][1], 0.0, data2[0][3], 0.0]]
init.extend(data2)
data2 = init
T_int2.insert(0, 0)
# Add the correct critical point entry
T_crit,p_crit, d_crit = 647.14, 220.64e5, 1/3.01e-3
data2.append([p_crit/p0,
              d_crit/r0,
              d_crit/r0,
              e_sg(p_crit,T_crit, pinf_expl, ql_exp, d_crit, gamma_expl)/p0,
              e_sg(p_crit,T_crit, pinf_expg, qg_exp, d_crit, gamma_expg)/p0])
T_int2.append(T_crit)
for i in range(len(cols)):
    results2[cols[i]] = np.asarray(data2)[:,i]

results2.to_csv('psat.csv', index=False, encoding='ascii',  header=False)
