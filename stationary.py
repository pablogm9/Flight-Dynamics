#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:58:44 2020

@author: pablo
"""


import numpy as np
import pandas as pd
import os
import subprocess


# ------------- A/C Parameters -------------
BEM = 13600*0.453592  #Basic Empty Mass [kg]
S = 30.00 #[m^2]

# ------------- Constants -------------
p_0 = 101325. #[Pa]
rho_0 = 1.225 #[kg/m^3]
gamma = 1.4 #[-]
R = 287.058 #[J/(kg K)]
a = -6.5/1000 #[deg C / m]
T_0 = 273.15+15 #[K]
g_0 = 9.80665 #[m/s^2]



# ------------- Read EXCEL file -------------

datasheet = pd.read_excel('Data/TESTFLIGHT2_Post_Flight_Datasheet.xlsx')
datasheet.columns = ['A','B','C','D','E','F','G','H','I','J','K','L','M']
datasheet.index = range(2,85)


# ------------- Extract parameters -------------


# --- Names, masses, blockfuel ---
passenger_names = datasheet.loc[[8,9,10,11,12,13,14,15,16],'D']
passenger_masses = datasheet.loc[[8,9,10,11,12,13,14,15,16],'H'] #[kg]
block_fuel = datasheet.loc[18,'D']*0.453592 #[kg]


# --- SERIES 1 ---
series1 = datasheet.loc[[28,29,30,31,32,33,34],['B','C','D','E','F','G','H','I','J']]
time_1 = series1['B'] #[hh:mm:ss]
et_1 = series1['C'] #[sec]
h_p_1 = series1['D'] #[ft]
IAS_1 = series1['E'] #[kts]
a_1 = series1['F'] #[deg]
FFL_1 = series1['G'] #[lbs/hr]
FFR_1 = series1['H'] #[lbs/hr]
F_used_1 = series1['I']*0.453592 #[kg]
TAT_1 = series1['J'] #[deg C]


# --- SERIES 2 ---
series2 = datasheet.loc[[44,45,46,47,48,49,50],['B','C','D','E','F','G','H','I','J']]
time_2 = series2['B'] #[hh:mm:ss]
et_2 = series2['C'] #[sec]
h_p_2 = series2['D'] #[ft]
IAS_2 = series2['E'] #[kts]
a_2 = series2['F'] #[deg]
FFL_2 = series2['G'] #[lbs/hr]
FFR_2 = series2['H'] #[lbs/hr]
F_used_2 = series2['I']*0.453592 #[kg]
TAT_2 = series2['J'] #[deg C]


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

# Used fuel [kg]
F_used_1 = np.array(F_used_1)
F_used_2 = np.array(F_used_2)

# Weights [N]
W_1 = np.ones(7)
W_2 = np.ones(7)


for i in range(len(W_1)):
    W_1[i]=(rampmass - F_used_1[i])*g_0
    W_2[i]=(rampmass - F_used_2[i])*g_0


# ------------- Compute Thrust -------------

# Pressure altitudes in [m]
hp_1 = np.array(h_p_1 * 0.3048)
hp_2 = np.array(h_p_2 * 0.3048)


# Calibrated airspeed in [m/s]
# V_C = V_IAS (neglecting instrument and position errors)
V_C_1 = IAS_1 * 0.514444
V_C_2 = IAS_2 * 0.514444


# Static Pressures [Pa]
p_1 = p_0*((1+a*hp_1/T_0)**(-g_0/(a*R)))
p_2 = p_0*((1+a*hp_2/T_0)**(-g_0/(a*R)))

# Mach numbers [-]
M_1 = np.sqrt(np.array(((2/(gamma-1))*(((1+(p_0/p_1)*(((1+((gamma-1)/(2*gamma))*(rho_0/p_0)*(V_C_1**2))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma))-1)),dtype= np.float))       
M_2 = np.sqrt(np.array(((2/(gamma-1))*(((1+(p_0/p_2)*(((1+((gamma-1)/(2*gamma))*(rho_0/p_0)*(V_C_2**2))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma))-1)),dtype= np.float))         

# T_ISA in [K], for temperature difference
T_ISA_1 = hp_1*a + T_0
T_ISA_2 = hp_2*a + T_0

# Total temperature in [K]
TAT_1 = TAT_1+273.15
TAT_2 = TAT_2+273.15

# Static temperature in [K]
# Obtained by correcting the TAT for ram rise
T_1 = TAT_1/(1+((gamma-1)/2)*(M_1**2))
T_2 = TAT_2/(1+((gamma-1)/2)*(M_2**2))

# Temperature differences
delta_T_1 = np.array(T_1 - T_ISA_1)
delta_T_2 = np.array(T_2 - T_ISA_2)

# Left engine fuel flows [kg/s]
FFL_1 = np.array(FFL_1*0.000125998)
FFL_2 = np.array(FFL_2*0.000125998)

# Right engine fuel flows [kg/s]
FFR_1 = np.array(FFR_1*0.000125998)
FFR_2 = np.array(FFR_2*0.000125998)


# Write .dat file
lines_1 =[]
lines_2 = []

for i in range(len(W_1)):
    line_1 = str(hp_1[i])+' '+str(M_1[i])+' '+str(delta_T_1[i])+' '+str(FFL_1[i])+' '+str(FFR_1[i])+'\n'
    line_2 = str(hp_2[i])+' '+str(M_2[i])+' '+str(delta_T_2[i])+' '+str(FFL_2[i])+' '+str(FFR_2[i])+'\n'
    lines_1.append(line_1)
    lines_2.append(line_2)

all_lines =  lines_1 + lines_2

input_file = open('StationaryData/matlab.dat','w')
input_file.writelines(all_lines)
input_file.close()


# Call thrust.exe



# Extract computed thrust values from new .dat file
'''
output_file = open('StationaryData/thrust.dat','r')
thrust = output_file.readlines()
output_file.close()
'''

# Delete created files, both input and output
os.remove('StationaryData/matlab.dat')
#os.remove('StationaryData/thrust.dat')




# ------------- Compute C_L -------------

# True airspeed [m/s]
V_TAS_1 = M_1*np.sqrt(np.array((gamma*R*T_1),dtype=np.float))
V_TAS_2 = M_2*np.sqrt(np.array((gamma*R*T_2),dtype=np.float))

# Density [kg/m^3]
rho_1 = p_1/(R*T_1)
rho_2 = p_2/(R*T_2)

C_L_1 = W_1/(0.5*rho_1*(V_TAS_1**2)*S)
C_L_2 = W_2/(0.5*rho_2*(V_TAS_2**2)*S)


# ------------- Compute C_D -------------




# ------------- Return C_L, C_D, Weight -------------


