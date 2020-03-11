#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 00:39:25 2020

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


# -----------------------------------------------
# ---------- FIRST STATIC MEASUREMENTS ----------
# -----------------------------------------------



# ------------- Compute Weight -------------

#BlockFuel = Block_fuel #Block Fuel mass in [kg]
BlockFuel = block_fuel

#Total Payload mass in [kg]
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

'''
# Plot regression as sanity check
z_p = np.poly1d(z)

C_L_range = np.linspace(C_L_squared[0],C_L_squared[-1],1000)
z_line = z_p(C_L_range)

plt.plot(C_L_squared,C_D,'rD')
plt.plot(C_L_range,z_line,'b-')
plt.show()
'''


# ------------- Plot C_D vs C_L using C_D_0 (and e?) -------------












# -----------------------------------------------
# ---------- SECOND STATIC MEASUREMENTS ---------
# -----------------------------------------------


# ------------- Compute Weight -------------

# Weights [N]
W_elevator = np.ones(7)

for i in range(len(W_elevator)):
    W_elevator[i]=(rampmass - F_used_elevator[i])*g_0


# ------------- Compute Thrust -------------

# Pressure altitudes in [m]
hp_elevator = np.array(h_p_elevator * 0.3048)

# Calibrated airspeed in [m/s]
# V_C = V_IAS (neglecting instrument and position errors)
V_C_elevator = IAS_elevator * 0.514444

# Static Pressures [Pa]
p_elevator = p_0*((1+T_a*hp_elevator/T_0)**(-g_0/(T_a*R)))

# Mach numbers [-]
M_elevator = np.sqrt(np.array(((2/(gamma-1))*(((1+(p_0/p_elevator)*(((1+((gamma-1)/(2*gamma))*(rho_0/p_0)*(V_C_elevator**2))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma))-1)),dtype= np.float))       

# T_ISA in [K], for temperature difference
T_ISA_elevator = (hp_elevator*T_a + T_0-273.15)+273.15

# Total temperature in [K]
TAT_elevayor = np.array(TAT_elevator,dtype=np.float)+273.15

# Static temperature in [K]
# Obtained by correcting the TAT for ram rise
T_elevator = TAT_elevator/(1+((gamma-1)/2)*(M_elevator**2))

# Temperature differences
delta_T_elevator = np.array(T_elevator - T_ISA_elevator)

# Left engine fuel flows [kg/s]
FFL_elevator = np.array(FFL_elevator*0.000125998)

# Right engine fuel flows [kg/s]
FFR_elevayor  = np.array(FFR_elevator*0.000125998)


# Write .dat file
lines_elevator =[]

for i in range(len(W)):
    line = str(hp_elevator[i])+' '+str(M_elevator[i])+' '+str(delta_T_elevator[i])+' '+str(FFL_elevator[i])+' '+str(FFR_elevator[i])+'\n'
    lines_elevator.append(line)

input_file = open('matlab.dat', 'w')
input_file.writelines(lines_elevator)
input_file.close()

#Call thrust.exe
subprocess.call(['thrust.exe'])

# Extract computed thrust values from new .dat file
output_file = open('thrust.dat','r')
thrust_lines_elevator = output_file.readlines()
output_file.close()

# Delete created files, both input and output
os.remove('matlab.dat')
os.remove('thrust.dat')

# Create thrust arrays 
thrust_L_elevator=[]
thrust_R_elevator=[]

for line in thrust_lines_elevator:
    values = line.replace('\n','').split('\t')
    thrust_L_elevator.append(float(values[0]))
    thrust_R_elevator.append(float(values[1]))

thrust_L_elevator = np.array(thrust_L_elevator)
thrust_R_elevator = np.array(thrust_R_elevator)


thrust_total_elevator = thrust_L_elevator+thrust_R_elevator


# ------------- True airspeed [m/s] -------------
V_TAS_elevator = M_elevator*np.sqrt(np.array((gamma*R*T_elevator),dtype=np.float))

# ------------- Density [kg/m^3] -------------
rho_elevator = p_elevator/(R*T_elevator)



# ------------- Shift in c.g. -------------








# ------------- Elevator trim curve -------------

# --- Elevator effectiveness ---
# This value should be changed to the one found from the c.g. shift 
C_m_delta = -1.1642 #[-]

# --- Dimensionless thrust moment arm (from Table C.2) ---
C_m_T_c = -0.0064  #[-]


# --- Engine inlet diameter ---
# (arbitrary value since apparently it doesn't matter, see page 8)
D = 1.

# --- Thrust coefficient --- 
T_c = np.array([])

for i in range(len(thrust_total)):
    T_c_i = thrust_total_elevator[i]/(rho_elevator[i]*(V_TAS_elevator[i]**2)*(D**2))
    T_c = np.append(T_c,T_c_i)


# --- Standard thrust coefficients --- 
# To calculate it, first calculate the standard thrust T_s 
# by rerunning thrust.exe using the value of m_f_s instead of FFL and FFR

# Standard engine fuel flow (from Table B.1)
m_f_s = 0.048 #[kg/s]


# Rerunning thrust.exe

standard_lines =[]

for i in range(len(W)):
    line = str(hp_elevator[i])+' '+str(M_elevator[i])+' '+str(delta_T_elevator[i])+' '+str(m_f_s)+' '+str(m_f_s)+'\n'
    standard_lines.append(line)

input_file = open('matlab.dat', 'w')
input_file.writelines(standard_lines)
input_file.close()

subprocess.call(['thrust.exe'])

output_file = open('thrust.dat','r')
standard_thrust_lines = output_file.readlines()
output_file.close()

os.remove('matlab.dat')
os.remove('thrust.dat')

standard_thrust_L = []
standard_thrust_R = []

for line in standard_thrust_lines:
    values = line.replace('\n','').split('\t')
    standard_thrust_L.append(float(values[0]))
    standard_thrust_R.append(float(values[1]))

standard_thrust_L = np.array(standard_thrust_L)
standard_thrust_R = np.array(standard_thrust_R)

standard_thrust = standard_thrust_L + standard_thrust_R


# Calculating standard thrust coefficients
T_c_s = np.array([])

for i in range(len(standard_thrust)):
    T_c_s_i = standard_thrust[i]/(rho[i]*(V_TAS[i]**2)*(D**2))
    T_c_s = np.append(T_c_s,T_c_s_i)


# --- Reduced elevator deflection --- 
    
d_e_reduced = np.array([])
    
for i in range(len(standard_thrust)): 
    a = delta_e[i]-((1/C_m_delta)*C_m_T_c *(T_c_s[i]-T_c[i]))
    d_e_reduced = np.append(d_e_reduced,a)


# --- Reduced EAS --- 
    
# Standard A/C weight (from Table B.1)
W_s = 60500 #[N]

# EAS
V_EAS_elevator = V_TAS_elevator*np.sqrt(rho_elevator/rho_0)

# Reduced EAS
V_EAS_reduced = V_EAS_elevator*np.sqrt(W_s/W_elevator)


'''
# --- Plot elevator trim curve ---
plt.plot(V_EAS_reduced,delta_e_reduced)
plt.xlabel('Reduced equivalent airspeed [m/s]')
plt.ylabel('Reduced elevator deflection for equilibrium [deg]')
plt.title('Elevator trim curve')
plt.grid()
plt.show()
'''


# ------------- Elevator control force curve -------------


# Reduced elevator force
F_e_reduced = np.array([])
    
for i in range(len(standard_thrust)): 
    a = F_e*(W_s/W_elevator[i])
    F_e_reduced = np.append(F_e_reduced,a)


'''
# --- Plot elevator control force curve ---
plt.plot(V_EAS_reduced,F_e_reduced)
plt.xlabel('Reduced equivalent airspeed [m/s]')
plt.ylabel('Reduced elevator control force [N]')
plt.title('Elevator control force curve')
plt.grid()
plt.show()
'''





