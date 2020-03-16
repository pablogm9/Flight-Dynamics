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

########## EIGENVALUES ##########

[e1, v1] = np.linalg.eig(A1)
[e2, v2] = np.linalg.eig(A2)

#print(e1, v1)
#print(e2, v2)

########## STATE SPACE MODEL ##########

sys1 = ml.ss(A1, B1, C1, D1)

sys2 = ml.ss(A2, B2, C2, D2)


#ml.pzmap(sys1) #eigenvalues plot
#ml.pzmap(sys2)


########### SIMULATION OF EIGENMOTIONS########
dt= 0.1


#Short Period
t1 = np.arange(0 , 10, dt)

u1 = inputcr(reference_delta_e, time, t1, 3634, 3636)

y1_r , T1_r, x1_r = ml.lsim(sys2, u1, t1)


# Phugoid
t2 = np.arange(0 , 1000, dt)

u2 = inputcr(reference_delta_e, time, t2, 3237, 3247)

y2_r , T2_r, x2_r = ml.lsim(sys2, u2, t2)


#Dutch Roll
t3 = np.arange(0 , 100, dt)

u3_t = inputcr(reference_delta_r, time, t3, 3717, 3718.8)

u3 = np.hstack((u3_t, np.zeros((len(t3), 1))))

y3_r , T3_r, x3_r = ml.lsim(sys1, u3, t3)


#Aperiodic roll
t4 = np.arange(0 , 12, dt)

u4_t = inputcr(reference_delta_a, time, t4, 3550, 3551)

u4 = np.hstack((np.zeros((len(t4), 1)), u4_t))

y4_r , T4_r, x4_r = ml.lsim(sys1, u4, t4)

#Spiral
t5 = np.arange(0, 100,dt)

u5_t = inputcr(reference_delta_a, time, t5, 3912, 3920)

u5 = np.hstack((np.zeros((len(t5), 1)), u5_t))

y5_r , T5_r, x5_r = ml.lsim(sys1, u5, t5)



##########PLOTS OF RESPONSES###########
#Short Period
f1 = plt.figure(1)
plt.plot(T1_r,y1_r[:,0],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Absoulte Velocity')
plt.title('Short Period Response - velocity',pad=10)

f2 = plt.figure(2)
plt.plot(T1_r,y1_r[:,1],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Angle of Attack[-]')
plt.title('Short Period Response - Angle of Attack',pad=10)

#Phugoid
f3 = plt.figure(3)
plt.plot(T2_r,y2_r[:,0],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Absoulte Velocity')
plt.title('Phugoid Response - velocity',pad=10)

f4 = plt.figure(4)
plt.plot(T2_r,y2_r[:,1],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Angle of Attack[-]')
plt.title('Phugoid Response - Angle of Attack',pad=10)

#Dutch Roll
f5 = plt.figure(5)
plt.plot(T3_r,y3_r[:,0],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Side Slip [-]')
plt.title('Dutch Roll Repsonse - Sideslip',pad=10)

f6 = plt.figure(6)
plt.plot(T3_r,y3_r[:,1],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Roll Angle [-]')
plt.title('Dutch Roll Repsonse - Roll Angle',pad=10)

#Aperiodic Roll
f7 = plt.figure(7)
plt.plot(T4_r,y4_r[:,0],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Side Slip [-]')
plt.title('Aperiodic Roll Repsonse - Sideslip',pad=10)

f8 = plt.figure(8)
plt.plot(T4_r,y4_r[:,1],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Roll Angle [-]')
plt.title('Aperiodic Roll Repsonse - Roll Angle',pad=10)

#Spiral

f9 = plt.figure(8)
plt.plot(T5_r,y5_r[:,0],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Side Slip [-]')
plt.title('Spiral Repsonse - Sideslip',pad=10)

f10 = plt.figure(10)
plt.plot(T5_r,y5_r[:,1],'r',label='Reference Data')
#plt.plot(T1_f,y1_r,'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Roll Angle [-]')
plt.title('Spiral Repsonse - Roll Angle',pad=10)
########## Print Commands ##########



#print(e1, v2)

#print(e2, v2)