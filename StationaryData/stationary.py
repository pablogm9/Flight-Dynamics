#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:58:44 2020

@author: pablo
"""


import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import os
import subprocess


# ------------- A/C Parameters -------------
BEM = 13600*0.453592  #Basic Empty Mass [kg]
S = 30.00 #[m^2]
b = 15.911 #[m]
A = (b**2)/S #[-]

# ------------- Constants -------------
p_0 = 101325. #[Pa]
rho_0 = 1.225 #[kg/m^3]
gamma = 1.4 #[-]
R = 287.058 #[J/(kg K)]
T_a = -6.5/1000 #[deg C / m]
T_0 = 273.15+15 #[K]
g_0 = 9.80665 #[m/s^2]



# ------------- Read EXCEL file -------------

datasheet = pd.read_excel('../Data/REFERENCE_Post_Flight_Datasheet.xlsx')
datasheet.columns = ['A','B','C','D','E','F','G','H','I','J','K','L','M']
datasheet.index = range(2,85)


# ------------- Extract parameters -------------

# --- Names, masses, blockfuel ---
passenger_names = datasheet.loc[[8,9,10,11,12,13,14,15,16],'D']
passenger_masses = datasheet.loc[[8,9,10,11,12,13,14,15,16],'H'] #[kg]
block_fuel = datasheet.loc[18,'D']*0.453592 #[kg]


# --- CL_CD SERIES ---
series = datasheet.loc[[28,29,30,31,32,33],['B','C','D','E','F','G','H','I','J']]
time = series['B'] #[hh:mm:ss]
et = series['C'] #[sec]
h_p = series['D'] #[ft]
IAS = series['E'] #[kts]
a = series['F'] #[deg]
FFL = series['G'] #[lbs/hr]
FFR = series['H'] #[lbs/hr]
F_used = np.array(series['I'])*0.453592 #[kg]
TAT = series['J'] #[deg C]



# --- ELEVATOR TRIM CURVE ---
elevator_series = datasheet.loc[[59,60,61,62,63,64,65],['B','C','D','E','F','G','H','I','J','K','L','M']]
time_elevator = elevator_series['B'] #[hh:mm:ss]
et_elevator = elevator_series['C'] #[sec]
h_p_elevator = elevator_series['D'] #[ft]
IAS_elevator  = elevator_series['E'] #[kts]
a_elevator = elevator_series['F'] #[deg]
delta_e = elevator_series['G'] #[deg]
delta_e_t = elevator_series['H'] #[deg]
F_e = elevator_series['I'] #[N]
FFL_elevator = elevator_series['J'] #[lbs/hr]
FFR_elevator = elevator_series['K'] #[lbs/hr]
F_used_elevator = elevator_series['L'] #[lbs]
TAT_elevator = elevator_series['M'] #[deg C]


# --- CG SHIFT ---
cg_series = datasheet.loc[[75,76],['B','C','D','E','F','G','H','I','J','K','L','M']]
time_cg = cg_series['B'] #[hh:mm:ss]
et_cg = cg_series['C'] #[sec]
h_p_cg = cg_series['D'] #[ft]
IAS_cg  = cg_series['E'] #[kts]
a_cg = cg_series['F'] #[deg]
delta_e_cg = cg_series['G'] #[deg]
delta_e_t_cg = cg_series['H'] #[deg]
F_e = cg_series['I'] #[N]
FFL_cg = cg_series['J'] #[lbs/hr]
FFR_cg = cg_series['K'] #[lbs/hr]
F_used_cg = cg_series['L'] #[lbs]
TAT_cg = cg_series['M'] #[deg C]


# --- TIME OF START OF EIGENMOTIONS ---
phugoid_time = datasheet.loc[[83],['D']] #[hh:mm:ss]
short_period_time = datasheet.loc[[84],['D']] #[hh:mm:ss]
dutch_roll_time = datasheet.loc[[83],['G']] #[hh:mm:ss]
dutch_roll_YD_time = datasheet.loc[[84],['G']] #[hh:mm:ss]
aper_roll_time = datasheet.loc[[83],['J']] #[hh:mm:ss]
spiral_time = datasheet.loc[[84],['J']] #[hh:mm:ss]



# ------------- Compute Weight -------------

#BlockFuel = Block_fuel #Block Fuel mass in kg [kg]
BlockFuel = block_fuel

#Total Payload mass in kg [kg]
Payload = np.sum(passenger_masses)

# Ramp weight [kg]
rampmass=(BEM+Payload+BlockFuel)

# Weights [N]
W = np.ones(6)

for i in range(len(W)):
    W[i]=(rampmass - F_used[i])*g_0


# ------------- Compute Thrust -------------

# Pressure altitudes in [m]
hp = np.array(h_p * 0.3048)

# Calibrated airspeed in [m/s]
# V_C = V_IAS (neglecting instrument and position errors)
V_C = IAS * 0.514444


# Static Pressures [Pa]
p = p_0*((1+T_a*hp/T_0)**(-g_0/(T_a*R)))

# Mach numbers [-]
M = np.sqrt(np.array(((2/(gamma-1))*(((1+(p_0/p)*(((1+((gamma-1)/(2*gamma))*(rho_0/p_0)*(V_C**2))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma))-1)),dtype= np.float))       

# T_ISA in [K], for temperature difference
T_ISA = (hp*T_a + T_0-273.15)+273.15

# Total temperature in [K]
TAT = np.array(TAT,dtype=np.float)+273.15

# Static temperature in [K]
# Obtained by correcting the TAT for ram rise
T = TAT/(1+((gamma-1)/2)*(M**2))

# Temperature differences
delta_T = np.array(T - T_ISA)

# Left engine fuel flows [kg/s]
FFL = np.array(FFL*0.000125998)

# Right engine fuel flows [kg/s]
FFR  = np.array(FFR*0.000125998)

'''
# Write .dat file
lines =[]

for i in range(len(W)):
    line = str(hp[i])+' '+str(M[i])+' '+str(delta_T[i])+' '+str(FFL[i])+' '+str(FFR[i])+'\n'
    lines.append(line)

input_file = open('matlab.dat', 'w')
input_file.writelines(lines)
input_file.close()

#Call thrust.exe
subprocess.call(['thrust.exe'])

# Extract computed thrust values from new .dat file
output_file = open('thrust.dat','r')
thrust_lines = output_file.readlines()
output_file.close()

# Delete created files, both input and output
os.remove('matlab.dat')
os.remove('thrust.dat')

# Create thrust arrays
thrust_L=[]
thrust_R=[]

for line in thrust_lines:
    values = line.replace('\n','').split('\t')
    thrust_L.append(float(values[0]))
    thrust_R.append(float(values[1]))

thrust_L = np.array(thrust_L)
thrust_R = np.array(thrust_R)
'''


thrust_L=np.array([3665.28, 2995.58, 2399.84, 1863.53, 1892.36, 2208.98])
thrust_R=np.array([3771.2,3057.48,2526.29,2016.03,2070.92,2405.45])

thrust_total=thrust_L+thrust_R


# ------------- Compute C_L -------------

# True airspeed [m/s]
V_TAS = M*np.sqrt(np.array((gamma*R*T),dtype=np.float))

# Density [kg/m^3]
rho = p/(R*T)

C_L = np.array(W/(0.5*rho*(V_TAS**2)*S),dtype=np.float)


# ------------- Compute C_D -------------
C_D = np.array(thrust_total/(0.5*rho*(V_TAS**2)*S),dtype=np.float)


# ------------- Calculation of C_D_0 and e -------------

# Linear regression of C_D vs C_L^2
C_L_squared = C_L**2
z = np.polyfit(C_L_squared,C_D,1) #Coefficients, highest power first

# C_D_0 and e
e = 1/(np.pi*A*z[0])
C_D_0 = z[1]


# Plot regression as sanity check
z_p = np.poly1d(z)

C_L_range = np.linspace(C_L_squared[0],C_L_squared[-1],1000)
z_line = z_p(C_L_range)

plt.plot(C_L_squared,C_D,'rD')
plt.plot(C_L_range,z_line,'b-')
plt.show()



# ------------- Linear interpolation of C_L vs alpha -------------




# ------------- Quadratic interpolation of C_D vs alpha -------------














