#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:58:44 2020

@author: pablo
"""


import numpy as np
import pandas as pd





# --- Read EXCEL file ---




# --- Compute weight ---



a = 'sku'




# --- Compute thrust ---

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




# --- Compute C_L ---




# --- Compute C_D ---





# --- Return C_L, C_D, Weight ---


