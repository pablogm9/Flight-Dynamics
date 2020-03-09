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
hp_1 = series1['D'] #[ft]
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
hp_2 = series2['D'] #[ft]
IAS_2 = series2['E'] #[kts]
a_2 = series2['F'] #[deg]
FFL_2 = series2['G'] #[lbs/hr]
FFR_2 = series2['H'] #[lbs/hr]
F_used_2 = series2['I'] #[lbs]
TAT_2 = series2['J'] #[deg C]



# --- ELEVATOR TRIM CURVE ---
elevator_series = datasheet.loc[[44,45,46,47,48,49,50],['B','C','D','E','F','G','H','I','J']]
time_2 = series2['B'] #[hh:mm:ss]
et_2 = series2['C'] #[sec]
hp_2 = series2['D'] #[ft]
IAS_2 = series2['E'] #[kts]
a_2 = series2['F'] #[deg]
FFL_2 = series2['G'] #[lbs/hr]
FFR_2 = series2['H'] #[lbs/hr]
F_used_2 = series2['I'] #[lbs]
TAT_2 = series2['J'] #[deg C]



# --- CG SHIFT ---



# --- EIGENMOTIONS ---



# --- Compute weight ---




# --- Compute thrust ---
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


# --- Compute C_L ---




# --- Compute C_D ---





# --- Return C_L, C_D, Weight ---


