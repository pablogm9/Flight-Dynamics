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
from math import pi

c      = 2.0569
b      = 15.911
#-------------------Definitions-------------------------
def statvalues(FlightType, Motion):
    parameters = np.array(
        ['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt', 'rh_engine_itt',
         'lh_engine_OP', 'rh_engine_OP', 'lh_engine_fan_N1', 'lh_engine_turbine_N2', 'rh_engine_fan_N1',
         'rh_engine_turbine_N2',
         'lh_engine_FU', 'rh_engine_FU', 'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll',
         'Ahrs1_Pitch', 'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate',
         'Ahrs1_bYawRate',
         'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc',
         'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas', 'Dadc1_tas',
         'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state', 'display_active_screen',
         'time'])

    reference_data, reference_headers, reference_descriptions = Read.get_data('ref_data')
    flight_data, flight_headers, flight_descriptions = Read.get_data('testflight')

    alpha_a_index = np.where(parameters == 'vane_AOA')[0].flat[0]
    velocity_index = np.where(parameters == 'Dadc1_tas')[0].flat[0]
    hp0_index = np.where(parameters == 'Dadc1_bcAlt')[0].flat[0]
    th0_index = np.where(parameters == 'Ahrs1_Pitch')[0].flat[0]
    time_index = np.where(parameters == 'time')[0].flat[0]

    time = np.array(reference_data[[reference_headers[time_index]]])
    reference_V0 = np.array(reference_data[reference_headers[velocity_index]])
    reference_hp0 = np.array(reference_data[reference_headers[hp0_index]])
    reference_th0 = np.array(reference_data[reference_headers[th0_index]])

    flight_V0 = np.array(flight_data[flight_headers[velocity_index]])
    flight_hp0 = np.array(flight_data[flight_headers[hp0_index]])
    flight_th0 = np.array(flight_data[flight_headers[th0_index]])
    flight_a = np.array(flight_data[flight_headers[alpha_a_index]])

    if Motion == 1:
        if FlightType == 1:
            time_ini = 3634
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0] * 0.51444
            hp0 = reference_hp0[cell_hp0] * 0.3048
            th0 = reference_th0[cell_th0] * pi / 180
            a01 = flight_a[cell_th0] * pi / 180

        elif FlightType == 2:
            time_ini = 3157
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            cell_a0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0] * 0.51444
            hp0 = flight_hp0[cell_hp0] * 0.3048
            th0 = flight_th0[cell_th0] * pi / 180
            a01 = flight_a[cell_a0] * pi / 180

    elif Motion == 2:
        if FlightType == 1:
            time_ini = 3236
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0] * 0.51444
            hp0 = reference_hp0[cell_hp0] * 0.3048
            th0 = reference_th0[cell_th0] * pi / 180
        elif FlightType == 2:
            time_ini = 3228
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            cell_a0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0] * 0.51444
            hp0 = flight_hp0[cell_hp0] * 0.3048
            th0 = flight_th0[cell_th0] * pi / 180
            a01 = flight_a[cell_a0] * pi / 180

    elif Motion == 3:
        if FlightType == 1:
            time_ini = 3717
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            cell_a0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0] * 0.51444
            hp0 = reference_hp0[cell_hp0] * 0.3048
            th0 = reference_th0[cell_th0] * pi / 180
        elif FlightType == 2:
            time_ini = 3479
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            cell_a0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0] * 0.51444
            hp0 = flight_hp0[cell_hp0] * 0.3048
            th0 = flight_th0[cell_th0] * pi / 180
            a01 = flight_a[cell_a0] * pi / 180
    elif Motion == 4:
        if FlightType == 1:
            time_ini = 3550
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0] * 0.51444
            hp0 = reference_hp0[cell_hp0] * 0.3048
            th0 = reference_th0[cell_th0] * pi / 180
        elif FlightType == 2:
            time_ini = 3607
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            cell_a0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0] * 0.51444
            hp0 = flight_hp0[cell_hp0] * 0.3048
            th0 = (flight_th0[cell_th0]) * pi / 180
            a01 = flight_a[cell_a0] * pi / 180
    elif Motion == 5:
        if FlightType == 1:
            time_ini = 3920
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0] * 0.51444
            hp0 = reference_hp0[cell_hp0] * 0.3048
            th0 = reference_th0[cell_th0] * pi / 180
        elif FlightType == 2:
            time_ini = 3675
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            cell_a0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0] * 0.51444
            hp0 = flight_hp0[cell_hp0] * 0.3048
            th0 = flight_th0[cell_th0] * pi / 180
            a01 = flight_a[cell_a0] * pi / 180

    else:
        print("Not a valid input")
        sys.exit()

    return a01, V0, th0


a0_1, V0_1, th0_1 = statvalues(2,1)
a0_2, V0_2, th0_2 = statvalues(2,2)
a0_3, V0_3, th0_3 = statvalues(2,3)
a0_4, V0_4, th0_4 = statvalues(2,4)
a0_5, V0_5, th0_5 = statvalues(2,5)



# ---------------- GETTING REFERENCE DATA ----------------

# Parameters --- Included for selecting data
parameters =  np.array(['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt', 'rh_engine_itt', 'lh_engine_OP', 'rh_engine_OP', 'lh_engine_fan_N1', 'lh_engine_turbine_N2', 'rh_engine_fan_N1', 'rh_engine_turbine_N2',
                'lh_engine_FU', 'rh_engine_FU', 'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll', 'Ahrs1_Pitch', 'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate', 'Ahrs1_bYawRate', 
                'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc', 'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas', 'Dadc1_tas', 
                'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state', 'display_active_screen', 'time'])

reference_data,reference_headers,reference_descriptions = Read.get_data('ref_data')
flight_data,flight_headers,flight_descriptions = Read.get_data('testflight')



#-----------------------Extracting Validation arrays----------------------------

time_index = np.where(parameters=='time')[0].flat[0]
alpha_a_index = np.where(parameters == 'vane_AOA')[0].flat[0]
theta_index =np.where(parameters == 'Ahrs1_Pitch')[0].flat[0]
q_index = np.where(parameters == 'Ahrs1_bPitchRate')[0].flat[0]
V_index = np.where(parameters == 'Dadc1_tas')[0].flat[0]
r_index = np.where(parameters == 'Dadc1_tas')[0].flat[0]



time = np.array(flight_data[[flight_headers[time_index]]])

alpha = np.array(flight_data[[flight_headers[alpha_a_index]]])*pi/180
theta =  np.array(flight_data[[flight_headers[theta_index]]])*pi/180
Pitch_rate_q =  np.array(flight_data[[flight_headers[q_index]]])*pi/180
V_tas = np.array(flight_data[[flight_headers[V_index]]])*0.51444

roll_rate_r = np.array(flight_data[[flight_headers[r_index]]])*pi/180

alpha_1 = np.array([i- a0_1 for i in alpha])
theta_1 = np.array([i-th0_1 for i in theta])
u_1 = np.array([(i-V0_1)/V0_1 for i in V_tas])
q_1 = np.array([i*c/V0_1 for i in Pitch_rate_q])

alpha_2 = np.array([i- a0_2 for i in alpha])
theta_2 = np.array([i-th0_2 for i in theta])
u_2 = np.array([(i-V0_2)/V0_2 for i in V_tas])
q_2 = np.array([i*c/V0_2 for i in Pitch_rate_q])

r_3 = np.array([i*b/(2*V0_3) for i in roll_rate_r])


# ---------------- SAMPLE PLOTS ----------------
'''

# Sample plot, AOA vs t
delta_r_index = np.where(parameters=='delta_r')[0].flat[0]
time_index = np.where(parameters=='time')[0].flat[0]
delta_a_index = np.where(parameters=='delta_a')[0].flat[0]


time = np.array(reference_data[[reference_headers[time_index]]])
reference_delta_r = np.array(reference_data[[reference_headers[delta_r_index]]])
reference_delta_a = np.array(reference_data[reference_headers[delta_a_index]])


f1 = plt.figure(1)
plt.plot(time,reference_delta_r,'r',label='delta R')
plt.plot(time,reference_delta_a,'b',label='delta A')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Delta_r [deg]')
plt.title('Delta R over time',pad=10)


# Sample plot, delta_e vs t
delta_e_index = np.where(parameters=='delta_e')[0].flat[0]
alpha_index = np.where(parameters=='vane_AOA')[0].flat[0]

reference_delta_e = np.array(reference_data[reference_headers[delta_e_index]])
reference_alpha_a = np.array(reference_data[reference_headers[alpha_index]])

f2 = plt.figure(2)
plt.plot(time,reference_delta_e,'r',label='delta e')
plt.plot(time,reference_alpha_a,'b',label='AOA')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Delta_e [deg]')
plt.title('Elevator deflection over time',pad=10)
plt.show()


# Sample plot, elevator_dte vs t
time_index = np.where(parameters=='time')[0].flat[0]
time = np.array(reference_data[[reference_headers[time_index]]])

trim_index = np.where(parameters=='elevator_dte')[0].flat[0]

reference_trim = np.array(reference_data[flight_headers[trim_index]])
flight_trim = np.array(flight_data[flight_headers[trim_index]])

f2 = plt.figure(2)
plt.plot(time,reference_trim,'r',label='Arbitrary data (sine function)')
plt.plot(time,flight_trim,'b',label='Flight data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Elevator trim tab deflection [deg]')
plt.title('Trim tab deflection over time',pad=10)

plt.show()
'''





