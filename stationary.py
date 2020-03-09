#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:58:44 2020

@author: pablo
"""


import numpy as np
import pandas as pd



# ------------- Read EXCEL file -------------

datasheet = pd.read_excel('Data/TESTFLIGHT2_Post_Flight_Datasheet.xlsx')
datasheet.columns = ['A','B','C','D','E','F','G','H','I','J','K','L','M']
datasheet.index = range(2,85)


# ------------- Extract parameters -------------

# --- Names, masses, blockfuel ---
passenger_names = datasheet.loc[[8,9,10,11,12,13,14,15,16],'D']
passenger_masses = datasheet.loc[[8,9,10,11,12,13,14,15,16],'H'] #[kg]
block_fuel = datasheet.loc[18,'D']*0.453592


# --- SERIES 1 ---
series1 = datasheet.loc[[28,29,39,31,32,33,34],['B','C','D','E','F','G','H','I','J']]
time_1 = series1['B'] #[hh:mm:ss]
et_1 = series1['C'] #[sec]
h_p_1 = series1['D'] #[ft]
IAS_1 = series1['E'] #[kts]
a_1 = series1['F'] #[deg]
FFL_1 = series1['G'] #[lbs/hr]
FFR_1 = series1['H'] #[lbs/hr]
F_used_1 = series1['I'] #[lbs]
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
F_used_2 = series2['I'] #[lbs]
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





# ------------- Compute Thrust -------------

# Pressure altitudes in [m]
hp_1 = h_p_1 * 0.3048
hp_2 = h_p_2 * 0.3048


# Calibrated airspeed in [m/s]
# V_C = V_IAS (neglecting instrument and position errors)
V_C_1 = IAS_1 * 0.514444
V_C_2 = IAS_2 * 0.514444


# Constants
p_0 = 101325. #[Pa]
rho_0 = 1.225 #[kg/m^3]
gamma = 1.4 #[-]
g_0 = 9.80665  #[m/s^2]
R = 287.058 #[J/(kg K)]
a = -6.5/1000 #[deg C / m]
T_0 = 273.15+15 #[K]

# Pressures
p_1 = p_0*((1+a*hp_1/T_0)**(-g_0/(a*R)))
p_2 = p_0*((1+a*hp_2/T_0)**(-g_0/(a*R)))

# Mach numbers [-]
M_1 = np.sqrt((2/(gamma-1))*(((1+(p_0/p_1)*(((1+((gamma-1)/(2*gamma))*(rho_0/p_0)*(V_C_1**2))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma))-1))         
M_2 = np.sqrt((2/(gamma-1))*(((1+(p_0/p_2)*(((1+((gamma-1)/(2*gamma))*(rho_0/p_0)*(V_C_2**2))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma))-1))         

# T_ISA in [K], for temperature difference
T_ISA_1 = hp_1*a + T_0
T_ISA_2 = hp_2*a + T_0

# Total temperature in [K], corrected for ram rise
TAT_1 = TAT_1+273.15
TAT_2 = TAT_2+273.15

T_1 = TAT_1/(1+((gamma-1)/2)*(M_1**2))
T_2 = TAT_2/(1+((gamma-1)/2)*(M_2**2))

# Temperature differences
delta_T_1 = T_1 - T_ISA_1
delta_T_2 = T_2 - T_ISA_2

# Left engine fuel flows [kg/s]
FFL_1 = FFL_1*0.000125998
FFL_2 = FFL_2*0.000125998

# Right engine fuel flows [kg/s]
FFR_1 = FFR_1*0.000125998
FFR_2 = FFR_2*0.000125998

# Write .mat file

# Call thrust.exe


# Extract computed thrust values



'''
# Pressure altitude [m]
h_p_index = np.where(reference_headers=='Dadc1_alt [ft]')[0].flat[0]
h_p = np.array(reference_data[reference_headers[h_p_index]])*0.3048

# Mach Number [-]
M_index = np.where(reference_headers=='Dadc1_mach [mach]')[0].flat[0]
M = np.array(reference_data[reference_headers[M_index]])

# Left engine fuel flow [kg/s]
left_fuelflow_index = np.where(reference_headers=='lh_engine_FMF [lbs/hr]')[0].flat[0]
left_fuelflow = np.array(reference_data[reference_headers[left_fuelflow_index]])

# Right engine fuel flow [kg/s]
right_fuelflow_index = np.where(reference_headers=='rh_engine_FMF [lbs/hr]')[0].flat[0]
right_fuelflow = np.array(reference_data[reference_headers[right_fuelflow_index]])

# Static temperature converted to [K]
T_static_index = np.where(reference_headers=='Dadc1_sat [deg C]')[0].flat[0]
T_static = np.array(reference_data[reference_headers[T_static_index]])+273.15

# ISA temperature [K]
T0 = 273.15+15 #[K]
a = -6.5/1000 #[deg C / m]
T_ISA = h_p*a + T0

# Temperature difference
delta_T = T_static - T_ISA

# Write parameters to .mat file


# Call thrust.exe 
'''


# ------------- Compute C_L -------------




# ------------- Compute C_D -------------




# ------------- Return C_L, C_D, Weight -------------


