#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 21:40:44 2020

@author: pablo
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Read





# ---------------- GETTING REFERENCE DATA ----------------

# Parameters --- Included for selecting data
parameters =  np.array(['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt', 'rh_engine_itt', 'lh_engine_OP', 'rh_engine_OP', 'lh_engine_fan_N1', 'lh_engine_turbine_N2', 'rh_engine_fan_N1', 'rh_engine_turbine_N2',
                'lh_engine_FU', 'rh_engine_FU', 'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll', 'Ahrs1_Pitch', 'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate', 'Ahrs1_bYawRate', 
                'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc', 'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas', 'Dadc1_tas', 
                'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state', 'display_active_screen', 'time'])

reference_data,reference_headers,reference_descriptions = Read.get_data('ref_data')
flight_data,flight_headers,flight_descriptions = Read.get_data('testflight2')



# ---------------- SAMPLE PLOTS ----------------

# Sample plot, AOA vs t
alpha_index = np.where(parameters=='vane_AOA')[0].flat[0]
time_index = np.where(parameters=='time')[0].flat[0]


time = np.array(reference_data[[reference_headers[time_index]]])
reference_alpha = np.array(reference_data[[reference_headers[alpha_index]]])
flight_alpha = np.array(flight_data[[flight_headers[alpha_index]]])


f1 = plt.figure(1)
plt.plot(time,reference_alpha,'r',label='reference data')
plt.plot(time,flight_alpha,'b',label='flight data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Angle of Attack [deg]')
plt.title('Angle of Attack over time',pad=10)


# Sample plot, delta_e vs t
delta_e_index = np.where(parameters=='delta_e')[0].flat[0]

reference_delta_e = np.array(reference_data[flight_headers[delta_e_index]])
flight_delta_e = np.array(flight_data[flight_headers[delta_e_index]])

f2 = plt.figure(2)
plt.plot(time,reference_delta_e,'r',label='reference data')
plt.plot(time,flight_delta_e,'b',label='flight data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Delta_e [deg]')
plt.title('Elevator deflection over time',pad=10)


plt.show()






