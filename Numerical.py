from Resources.Cit_par import *
import scipy as np
import control.matlab as ml
import matplotlib.pyplot as plt
import control as ctr
import Read

#####DEFINITIONS######

def inputcr(delta_values, time_value, t_array, time_ini, time_fin): #(deflections values, time of flight, time of simulation, start time of input, end time of input)
     cell_ini = np.where(time_value == time_ini)[0][0]
     cell_fin = np.where(time_value == time_fin)[0][0]
     values = delta_values[cell_ini:cell_fin]
     delta_values_a = np.array([values])
     delta_values_aa = np.array([i * 0.017455 for i in delta_values_a])
     missing = len(t_array) - len(values)
     zeros = np.zeros((1, missing))
     delta_values_aaa = np.hstack((delta_values_aa, zeros))
     delta_values_aaaa = np.transpose(delta_values_aaa)

     return delta_values_aaaa


####GETTING DATA FOR INPUTS###########

parameters =  np.array(['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt', 'rh_engine_itt', 'lh_engine_OP', 'rh_engine_OP', 'lh_engine_fan_N1', 'lh_engine_turbine_N2', 'rh_engine_fan_N1', 'rh_engine_turbine_N2',
                'lh_engine_FU', 'rh_engine_FU', 'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll', 'Ahrs1_Pitch', 'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate', 'Ahrs1_bYawRate',
                'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc', 'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas', 'Dadc1_tas',
                'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state', 'display_active_screen', 'time'])

reference_data,reference_headers,reference_descriptions = Read.get_data('ref_data')
flight_data,flight_headers,flight_descriptions = Read.get_data('testflight2')

delta_e_index = np.where(parameters=='delta_e')[0].flat[0]
delta_r_index = np.where(parameters=='delta_r')[0].flat[0]
delta_a_index = np.where(parameters=='delta_a')[0].flat[0]
time_index = np.where(parameters=='time')[0].flat[0]

reference_delta_e = np.array(reference_data[reference_headers[delta_e_index]])
reference_delta_r = np.array(reference_data[reference_headers[delta_r_index]])
reference_delta_a = np.array(reference_data[reference_headers[delta_a_index]])
time = np.array(reference_data[[reference_headers[time_index]]])


########## ASYMMETRIC EQUATIONS OF MOTION IN STATE-SPACE FORM ##########

P1 = np.array([[(CYb - 2 * mub) * (c/V0), 0, 0, 0],
               [0, -(1/2), 0, 0],
               [0, 0, - 4 * mub * KX2 * (c/V0), 4 * mub * KXZ * (c/V0)],
               [Cnbdot * (c/V0), 0, 4 * mub * KXZ * (c/V0), - 4 * mub * KZ2 * (c/V0)]])

Q1 = np.array([[-CYb, -CL, -CYp, 4 * mub - CYr],
               [0, 0, -1, 0],
               [-Clb, 0, -Clp, -Clr],
               [-Cnb, 0, -Cnp, -Cnr]])

P11 = np.linalg.inv(P1)

R1 = np.array([[-CYda, -CYdr],
               [0, 0],
               [-Clda, -Cldr],
               [-Cnda, -Cndr]])

A1 = np.matmul(P11, Q1)

B1 = np.matmul(P11, R1)

C1 = np.identity(4)

D1 = np.array([[0, 0], [0, 0], [0, 0], [0, 0]])


########## SYMMETRIC EQUATIONS OF MOTION IN STATE-SPACE FORM ##########

P2 = np.array([[- 2 * muc * (c/V0), 0, 0, 0],
              [0, (CZadot - 2 * muc) * (c/V0), 0, 0],
              [0, 0, - (c/V0), 0],
              [0, Cmadot * (c/V0), 0, - 2 * muc * KY2 * (c/V0)]])

Q2 = np.array([[-CXu, -CXa, -CZ0, -CXq],
               [-CZu, -CZa, CX0, -(CZadot + 2 * muc)],
               [0, 0, 0, -1],
               [-Cmu, -Cma, 0, -Cmq]])

P21 = np.linalg.inv(P2)

R2 = np.array([[-CXde],
               [-CZde],
               [0],
               [-Cmde]])


#print(np.shape(P21))
#print(np.shape(Q2))
#print(np.shape(R2))


A2 = np.matmul(P21, Q2)

B2 = np.matmul(P21, R2)

C2 = np.identity(4)

D2 = np.array([[0], [0], [0], [0]])

#print(A2)
#print(B2)
#print(C2)
#print(D2)

########## Eigenvalues and Eigenvectors ##########

[e1, v1] = np.linalg.eig(A1)
[e2, v2] = np.linalg.eig(A2)

#print(e1, v1)
#print(e2, v2)

########## Continuous time-state-space Model ##########

#t = np.arange(0 , 100, 1)

sys1 = ml.ss(A1, B1, C1, D1)

sys2 = ml.ss(A2, B2, C2, D2)


ml.pzmap(sys1) #eigenvalues plot
ml.pzmap(sys2)


########### Simulation of Eigen Motions########
dt= 0.1
t = np.arange(0 , 100, dt)

#Short Period

t_delta = int(10/dt) #time the elevator is deflected
s_def = -0.031 #angle of deflection[rad]
s_e = np.ones((t_delta,1))
s_en = np.array([i*s_def for i in s_e])
zeros = np.zeros((900,1))
s_e = np.vstack((s_en,zeros))

y , T, x = ml.lsim(sys2, s_e, t)


#plt.plot(T, y[:,1])
#plt.show()

# Phugoid
t2 = np.arange(0 , 1000, dt)

u2 = inputcr(reference_delta_e, time, t2, 3237, 3247)


y , T, x = ml.lsim(sys2, u2, t2)
#plt.plot(T, y[:,0])
#plt.show()

#Dutch Roll
t3 = np.arange(0 , 100, dt)

u3_t = inputcr(reference_delta_r, time, t3, 3717, 3718.8)

u3 = np.hstack((u3_t, np.zeros((len(t3), 1))))

y , T, x = ml.lsim(sys1, u3, t3)
#plt.plot(T, y[:,0])
#plt.show()

#Aperiodic roll
t4 = np.arange(0 , 12, dt)

u4_t = inputcr(reference_delta_a, time, t4, 3550, 3551)

u4 = np.hstack((np.zeros((len(t4), 1)), u4_t))

y , T, x = ml.lsim(sys1, u4, t4)
#plt.plot(T, y[:,0])
#plt.show()

#Spiral
t5 = np.arange(0, 100,dt)

u5_t = inputcr(reference_delta_a, time, t5, 3912, 3920)


u5 = np.hstack((np.zeros((len(t5), 1)), u5_t))

y , T, x = ml.lsim(sys1, u5, t5)
#plt.plot(T, y[:,1])
#plt.show()





########## Print Commands ##########


#print(e1, v2)

#print(e2, v2)