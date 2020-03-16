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
BEM = 9165*0.453592  #Basic Empty Mass [kg]
S = 30.00 #[m^2]
b = 15.911 #[m]
A = (b**2)/S #[-]
c = 2.0569 #[m]


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
delta_e = np.array(elevator_series['G']) #[deg]
delta_e_t = elevator_series['H'] #[deg]
F_e_elevator = np.array(elevator_series['I'],dtype=np.float) #[N]
FFL_elevator = elevator_series['J'] #[lbs/hr]
FFR_elevator = elevator_series['K'] #[lbs/hr]
F_used_elevator = np.array(elevator_series['L'])*0.453592 #[kg]
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
F_used_cg = np.array(cg_series['L']) #[lbs]
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

#thrust_total = thrust_L + thrust_R

thrust_total = np.array([7436.48, 6053.06, 4926.13, 3879.56, 3963.28, 4614.43])


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

# ------------- Calculating C_L_alpha, alpha_O -------------


C_L_alpha_coefficients = np.polyfit(np.array(a,dtype=np.float),C_L,1)
C_L_alpha_polynomial = np.poly1d(C_L_alpha_coefficients)

alpha_range = np.linspace(-3,12,10000)
C_L_range = C_L_alpha_polynomial(alpha_range)

# Slope of C_L vs alpha curve
C_L_alpha = C_L_alpha_coefficients[0]

# Alpha_0 
idx = np.where(C_L_range==np.abs(C_L_range).min())
alpha_0 = alpha_range[idx][0]


# ------------- Calculating C_D_alpha -------------

C_D_range = C_D_0 + (1/(np.pi*A*e))*(C_L_range**2)

# Slope of C_D vs alpha curve







# -----------------------------------------------
# ---------- SECOND STATIC MEASUREMENTS ---------
# -----------------------------------------------


# ------------- Compute Weight -------------

# Weights [N]
W_elevator = np.ones(7)

for i in range(len(W_elevator)):
    W_elevator[i] = (rampmass - F_used_elevator[i]) * g_0

# ------------- Compute Thrust -------------

# Pressure altitudes in [m]
hp_elevator = np.array(h_p_elevator * 0.3048)

# Calibrated airspeed in [m/s]
# V_C = V_IAS (neglecting instrument and position errors)
V_C_elevator = IAS_elevator * 0.514444

# Static Pressures [Pa]
p_elevator = p_0 * ((1 + T_a * hp_elevator / T_0) ** (-g_0 / (T_a * R)))

# Mach numbers [-]
M_elevator = np.sqrt(np.array(((2 / (gamma - 1)) * (((1 + (p_0 / p_elevator) * (((1 + ((gamma - 1) / (2 * gamma)) * (
            rho_0 / p_0) * (V_C_elevator ** 2)) ** (gamma / (gamma - 1))) - 1)) ** ((gamma - 1) / gamma)) - 1)),
                              dtype=np.float))

# T_ISA in [K], for temperature difference
T_ISA_elevator = (hp_elevator * T_a + T_0 - 273.15) + 273.15

# Total temperature in [K]
TAT_elevator = np.array(TAT_elevator, dtype=np.float) + 273.15

# Static temperature in [K]
# Obtained by correcting the TAT for ram rise
T_elevator = TAT_elevator / (1 + ((gamma - 1) / 2) * (M_elevator ** 2))

# Temperature differences
delta_T_elevator = np.array(T_elevator - T_ISA_elevator)

# Left engine fuel flows [kg/s]
FFL_elevator = np.array(FFL_elevator * 0.000125998)

# Right engine fuel flows [kg/s]
FFR_elevator = np.array(FFR_elevator * 0.000125998)

'''
# Write .dat file
lines_elevator = []

for i in range(len(W_elevator)):
    line = str(hp_elevator[i]) + ' ' + str(M_elevator[i]) + ' ' + str(delta_T_elevator[i]) + ' ' + str(
        FFL_elevator[i]) + ' ' + str(FFR_elevator[i]) + '\n'
    lines_elevator.append(line)

input_file = open('matlab.dat', 'w')
input_file.writelines(lines_elevator)
input_file.close()

# Call thrust.exe
subprocess.call(['thrust.exe'])

# Extract computed thrust values from new .dat file
output_file = open('thrust.dat', 'r')
thrust_lines_elevator = output_file.readlines()
output_file.close()

# Delete created files, both input and output
os.remove('matlab.dat')
os.remove('thrust.dat')

# Create thrust arrays
thrust_L_elevator = []
thrust_R_elevator = []

for line in thrust_lines_elevator:
    values = line.replace('\n', '').split('\t')
    thrust_L_elevator.append(float(values[0]))
    thrust_R_elevator.append(float(values[1]))

thrust_L_elevator = np.array(thrust_L_elevator)
thrust_R_elevator = np.array(thrust_R_elevator)
'''

#thrust_total_elevator = thrust_L_elevator + thrust_R_elevator
thrust_total_elevator = np.array([4009.72, 4103.43, 4173.47, 4256.74, 3943.49, 3941.25, 3882.84])

# ------------- True airspeed [m/s] -------------
V_TAS_elevator = M_elevator * np.sqrt(np.array((gamma * R * T_elevator), dtype=np.float))

# ------------- Density [kg/m^3] -------------
rho_elevator = p_elevator / (R * T_elevator)

# ------------- Shift in c.g. -------------

#--- Calculating C N ----

#Total Mass in [kg]
W_cg = rampmass-F_used_cg*0.453592

# Pressure altitudes in [m]
hp_cg = np.array(h_p_cg * 0.3048)

# Calibrated airspeed in [m/s]
V_C_cg = IAS_cg* 0.514444

# Static Pressures [Pa]
p_cg = p_0*((1+T_a*hp_cg/T_0)**(-g_0/(T_a*R)))

# Mach numbers [-]
M_cg = np.sqrt(np.array(((2/(gamma-1))*(((1+(p_0/p_cg)*(((1+((gamma-1)/(2*gamma))*(rho_0/p_0)*(V_C_cg**2))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma))-1)),dtype= np.float))

# Total temperature in [K]
TAT_cg= np.array(TAT_cg,dtype=np.float)+273.15

# Static temperature in [K] - Obtained by correcting the TAT for ram rise
T_cg = TAT_cg/(1+((gamma-1)/2)*(M_cg**2))

# True airspeed [m/s]
V_TAS_cg = M_cg*np.sqrt(np.array((gamma*R*T_cg),dtype=np.float))

# Density [kg/m^3]
rho_cg = p_cg/(R*T_cg)

# Coefficient of normal Force:C N
C_N_cg = np.array((W_cg*g_0)/(0.5*rho_cg*(V_TAS_cg**2)*S),dtype=np.float)

# Averaging out the C_N at both points:
C_N_cg_av=np.mean(C_N_cg)


# --- Shift in CG Calculations ----
#Value of BEM Moment in [lbs*inch]
moment_BEM = 2672972.25

#Array of x_cg of each passenger in [inches]
x_cg_payload = np.array([131,131,170,214,214,251,251,288,288])

#Array of moment of each passenger in [lbs*inch]
moment_payload = x_cg_payload*np.array(passenger_masses)*2.20462

#Value of total moment of all passengers in [lbs*inch]
moment_total_payload = np.sum(moment_payload)

#Array of 2 fuel loads in [lbs]
fuel_load_cg = np.array([block_fuel,block_fuel])*2.20462-F_used_cg

#Array of 2 moment of fuel loads [lbs*inch]
moment_fuel_load = np.array([9036.2144,8953.344]) # Choose this Value according to the value of fuel_load_cg from E.2 [Inch*pounds]

#Total moment of aircraft [lbs*inch]
moment_total_aircraft = np.array([moment_BEM,moment_BEM])+moment_fuel_load+np.array([moment_total_payload,moment_total_payload])

#Array of 2 cg locations [inches] from the datum line
x_cg_max = (moment_BEM+moment_total_payload+11451.85)/(rampmass*2.20462) #Once again choose the last value according to your Fuel Load
x_cg = moment_total_aircraft/(W_cg*2.20462)


# --- Calculating C_m_delta ---

# Change in Elevator deflection
Delta_e_cg = np.array(delta_e_cg)[1]-np.array(delta_e_cg)[0]

# C_m_delta
C_m_delta = (1/Delta_e_cg)*C_N_cg_av*((x_cg[1]-x_cg[0])/c)


# ------------- Elevator trim curve -------------

# --- Elevator effectiveness ---
# This value should be changed to the one found from the c.g. shift
#C_m_delta = -1.1642  # [-]

# --- Dimensionless thrust moment arm (from Table C.2) ---
C_m_T_c = -0.0064  # [-]

# --- Engine inlet diameter ---
# Engine model: Pratt & Whitney Canada JT15D-4
# Value taken from all-aero.com
D = 0.686  # [m]

# --- Thrust coefficient ---
T_c = np.array([])

for i in range(len(thrust_total_elevator)):
    T_c_i = thrust_total_elevator[i] / (rho_elevator[i] * (V_TAS_elevator[i] ** 2) * (D ** 2))
    T_c = np.append(T_c, T_c_i)

# --- Standard thrust coefficients ---
# To calculate it, first calculate the standard thrust T_s
# by rerunning thrust.exe using the value of m_f_s instead of FFL and FFR

# Standard engine fuel flow (from Table B.1)
m_f_s = 0.048  # [kg/s]

'''
# Rerunning thrust.exe
standard_lines = []

for i in range(len(W_elevator)):
    line = str(hp_elevator[i]) + ' ' + str(M_elevator[i]) + ' ' + str(delta_T_elevator[i]) + ' ' + str(m_f_s) + ' ' + str(m_f_s) + '\n'
    standard_lines.append(line)

input_file = open('matlab.dat', 'w')
input_file.writelines(standard_lines)
input_file.close()

subprocess.call(['thrust.exe'])

output_file = open('thrust.dat', 'r')
standard_thrust_lines = output_file.readlines()
output_file.close()

os.remove('matlab.dat')
os.remove('thrust.dat')

standard_thrust_L = []
standard_thrust_R = []

for line in standard_thrust_lines:
    values = line.replace('\n', '').split('\t')
    standard_thrust_L.append(float(values[0]))
    standard_thrust_R.append(float(values[1]))

standard_thrust_L = np.array(standard_thrust_L)
standard_thrust_R = np.array(standard_thrust_R)
'''

#standard_thrust = standard_thrust_L + standard_thrust_R
standard_thrust = np.array([2678.12, 2801.7 , 2914.08, 3040.54, 2593., 2508.48, 2355.64])


# Calculating standard thrust coefficients
T_c_s = np.array([])

for i in range(len(standard_thrust)):
    T_c_s_i = standard_thrust[i] / (rho_elevator[i] * (V_TAS_elevator[i] ** 2) * (D ** 2))
    T_c_s = np.append(T_c_s, T_c_s_i)

# --- Reduced elevator deflection ---

delta_e_reduced = np.array([])

for i in range(len(standard_thrust)):
    b = delta_e[i] - ((1 / C_m_delta) * C_m_T_c * (T_c_s[i] - T_c[i]))
    delta_e_reduced = np.append(delta_e_reduced, b)

# --- Reduced EAS ---

# Standard A/C weight (from Table B.1)
W_s = 60500  # [N]

# EAS
V_EAS_elevator = V_TAS_elevator * np.sqrt(np.array(rho_elevator,dtype=np.float) / rho_0)

# Reduced EAS
V_EAS_reduced = V_EAS_elevator * np.sqrt(W_s / W_elevator)

'''
# --- Plot elevator trim curve (VELOCITY) ---

# C_m_0 (from Table C.2)
C_m_0 = 0.0297


plt.scatter(V_EAS_reduced,delta_e_reduced)
plt.xlabel('Reduced equivalent airspeed [m/s]')
plt.ylabel('Reduced elevator deflection for equilibrium [deg]')
plt.ylim((np.max(delta_e_reduced)*1.5,np.min(delta_e_reduced)*1.1))
plt.hlines(-C_m_0/C_m_delta,np.min(V_EAS_reduced),np.max(V_EAS_reduced),colors='r',linestyles='dashed',label='delta_e_eq for neutral stability')
plt.title('Elevator trim curve (vs V_EAS')
plt.legend()
plt.grid()
plt.show()
'''


# ------------- Compute C_m_alpha -------------
# (from trim curve slope and C_m_delta)

# Linear regression of delta_e_reduced vs a_elevator
v = np.polyfit(np.array(a_elevator,dtype=np.float),np.array(delta_e_reduced,dtype=np.float),1) #Coefficients, highest power first

elevator_alpha_slope = v[0]
elevator_alpha_intercept = v[1]
C_m_alpha = elevator_alpha_slope*-C_m_delta


# ------------- Elevator control force curve -------------


# Reduced elevator force
F_e_reduced = np.array([])

for i in range(len(standard_thrust)):
    b = F_e_elevator[i] * (W_s / W_elevator[i])
    F_e_reduced = np.append(F_e_reduced, b)




# -----------------------------------------------
# ------------------- PLOTS ---------------------
# -----------------------------------------------

# ------------- C_L vs Alpha -------------
'''
angle_of_attack=np.array(a,dtype = np.float)

plt.plot(alpha_range,C_L_range,'r--',label=r'Constructed C$_L$ vs. $\alpha$ curve',zorder=1)
plt.scatter(angle_of_attack,C_L,label='Measurement points',marker='D',color='g',zorder=2,s=75,facecolors='None')
plt.scatter(alpha_0,C_L_range[idx],marker='o',label=r'Zero-lift angle of attack, $\alpha_0$',zorder=3,s=75,color='b')
plt.axvline(0,color='k')
plt.axhline(0,color='k')
plt.xlabel(r'Angle of Attack, $\alpha$ [deg]',fontsize=13)
plt.ylabel(r'Coefficient of Lift, C$_L$ [-]',fontsize=13)
plt.xlim((-3,12))
plt.title(r'C$_L$ vs. $\alpha$ (linear part)',fontsize=17,pad=25)
plt.legend()
plt.grid()
plt.show()
'''


# ------------- C_D vs Alpha -------------
'''
plt.plot(alpha_range,C_D_range,'r--',label=r'Constructed C$_D$ vs. $\alpha$ curve',zorder=1)
plt.scatter(angle_of_attack,C_D,label='Measurement points',marker='D',color='g',zorder=2,s=75,facecolors='None')
plt.scatter(alpha_0,C_D_range[idx],marker='o',label=r'Zero-lift drag coefficient, $C_{D_0}$',zorder=3,s=75,color='b')
plt.axvline(0,color='k')
plt.axhline(0,color='k')
plt.xlabel(r'Angle of Attack, $\alpha$ [deg]',fontsize=13)
plt.ylabel(r'Coefficient of Drag, C$_D$ [-]',fontsize=13)
plt.xlim((-3,12))
plt.title(r'C$_D$ vs. $\alpha$',fontsize=17,pad=25)
plt.legend()
plt.grid()
plt.show()
'''


# ------------- Drag polar, C_L vs C_D -------------
'''
plt.plot(C_D_range,C_L_range,'r--',label=r'Constructed C$_L$ vs. C$_D$ curve',zorder=1)
plt.scatter(C_D,C_L,label='Measurement points',marker='D',color='g',zorder=2,s=75,facecolors='None')
plt.scatter(C_D_range[idx],C_L_range[idx],marker='o',label=r'Zero-lift drag coefficient, $C_{D_0}$',zorder=3,s=75,color='b')
plt.axvline(0,color='k')
plt.axhline(0,color='k')
plt.xlabel(r'Coefficient of Drag, C$_D$ [-]',fontsize=13)
plt.ylabel(r'Coefficient of Lift, C$_L$ [-]',fontsize=13)
#plt.xlim((-3,12))
plt.title(r'Drag polar (C$_L$ vs. C$_D$)',fontsize=17,pad=25)
plt.legend(loc='upper left')
plt.grid()
plt.show()
'''

# --- Plot elevator trim curve (AOA) ---

# C_m_0 (from Table C.2)
C_m_0 = 0.0297

'''
delta_range = elevator_alpha_intercept + elevator_alpha_slope*alpha_range

plt.plot(alpha_range,delta_range,'r--',label=r'Constructed trim curve',zorder=1)
plt.scatter(np.array(a_elevator,dtype=np.float),delta_e_reduced,label='Measurement points',marker='D',color='g',zorder=2,s=75,facecolors='None')

plt.ylim((-5,5))
plt.gca().invert_yaxis()

plt.axvline(0,color='k')
plt.axhline(0,color='k')
plt.axhline(-C_m_0/C_m_delta,color='k',linestyle='--')

plt.xlim((0,10))

plt.xlabel(r'Angle of Attack, $\alpha$ [deg]',fontsize=13)
plt.ylabel(r'Reduced elevator deflection, $\delta^*_{e_{eq}}$ [deg]',fontsize=13)
plt.title(r'Elevator trim curve (vs. $\alpha$)',fontsize=17,pad=25)

plt.grid()
plt.legend()
plt.show()
'''


# --- Plot elevator trim curve (V_EAS_reduced) ---


# Stall speed in clean cruise configuration
# From all-aero.com, 95 kts. Could be ISA but it does NOT matter, 
# since it is only used as an indication in the graph

V_stall = 48.87 #[m/s] 


'''
w = np.polyfit((1/(V_EAS_reduced**2)),delta_e_reduced,1)

V_range = np.linspace(V_stall,100,10000)
delta_V_range = (1/(V_range**2))*w[0] + w[1]

plt.plot(V_range,delta_V_range,'r--',label=r'Constructed trim curve',zorder=1)
plt.scatter(V_EAS_reduced,delta_e_reduced,label='Measurement points',marker='D',color='g',zorder=2,s=75,facecolors='None')
plt.ylim((-3,3))

plt.gca().invert_yaxis()

plt.axhline(0,color='k')
plt.axvline(V_stall,color='b',linestyle='--',label='Stall speed (clean configuration)')

plt.xlabel(r'Reduced equivalent airspeed, $\tilde{V}_{EAS}$ [m/s]',fontsize=13)
plt.ylabel(r'Reduced elevator deflection, $\delta^*_{e_{eq}}$ [deg]',fontsize=13)
plt.title(r'Elevator trim curve (vs. $\tilde{V}_{EAS}$)',fontsize=17,pad=25)

plt.grid()
plt.legend(loc='lower center')
plt.show()

'''


# --- Plot elevator force curve (V_EAS_reduced) ---

'''
t = np.polyfit(V_EAS_reduced**2,F_e_reduced,1)
V_range = np.linspace(V_stall,100,10000)
F_range = (V_range**2)*t[0] + t[1]


plt.plot(V_range,F_range,'r--',label=r'Constructed control force curve',zorder=1)
plt.scatter(V_EAS_reduced,np.array(F_e_reduced,dtype=np.float),label='Measurement points',marker='D',color='g',zorder=2,s=75,facecolors='None')


plt.axhline(0,color='k')
plt.axvline(V_stall,color='b',linestyle='--',label='Stall speed (clean configuration)')

plt.ylim((-60,100))
plt.gca().invert_yaxis()


plt.xlabel(r'Reduced equivalent airspeed, $\tilde{V}_{EAS}$ [m/s]',fontsize=13)
plt.ylabel(r'Reduced stick force, $F^*_e$ [N]',fontsize=13)
plt.title(r'Elevator control force curve (vs. $\tilde{V}_{EAS}$)',fontsize=17,pad=25)

plt.grid()
plt.legend(loc='lower center')
plt.show()

'''







