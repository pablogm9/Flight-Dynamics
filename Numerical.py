from Resources.Cit_par import *
import scipy as np
import control.matlab as ml
import matplotlib.pyplot as plt
import control as ctr
import Read
####GETTING DATA FOR INPUTS###########

parameters =  np.array(['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt', 'rh_engine_itt', 'lh_engine_OP', 'rh_engine_OP', 'lh_engine_fan_N1', 'lh_engine_turbine_N2', 'rh_engine_fan_N1', 'rh_engine_turbine_N2',
                'lh_engine_FU', 'rh_engine_FU', 'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll', 'Ahrs1_Pitch', 'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate', 'Ahrs1_bYawRate',
                'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc', 'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas', 'Dadc1_tas',
                'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state', 'display_active_screen', 'time'])

reference_data,reference_headers,reference_descriptions = Read.get_data('ref_data')
flight_data,flight_headers,flight_descriptions = Read.get_data('testflight2')


########## ASYMMETRIC EQUATIONS OF MOTION IN STATE-SPACE FORM ##########

y = [(V0/b) * (CYb/(2 * mub)),
     (V0/b) * (CL/(2 * mub)),
     (V0/b) * (CYp/(2 * mub)),
     (V0/b) * ((CYr - 4 * mub)/(2 * mub)),
     (V0/b) * (CYda/(2 * mub)),
     (V0/b) * (CYdr/(2 * mub))]

l = [(V0/b) * ((Clb * KZ2 + Cnb * KXZ)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     0,
     (V0/b) * ((Clp * KZ2 + Cnp * KXZ)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     (V0/b) * ((Clr * KZ2 + Cnr * KXZ)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     (V0/b) * ((Clda * KZ2 + Cnda * KXZ)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     (V0/b) * ((Cldr * KZ2 + Cndr * KXZ)/(4 * mub * (KX2 * KZ2 - KXZ**2)))]

n = [(V0/b) * ((Clb * KXZ + Cnb * KX2)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     0,
     (V0/b) * ((Clp * KXZ + Cnp * KX2)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     (V0/b) * ((Clr * KXZ + Cnr * KX2)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     (V0/b) * ((Clda * KXZ + Cnda * KX2)/(4 * mub * (KX2 * KZ2 - KXZ**2))),
     (V0/b) * ((Cldr * KXZ + Cndr * KX2)/(4 * mub * (KX2 * KZ2 - KXZ**2)))]

A1 = np.array([[y[0], y[1], y[2], y[3]], [0, 0, 2 * (V0/b), 0], [l[0], l[1], l[2], l[3]], [n[0], n[1], n[2], n[3]]])

B1 = np.array([[0, y[5]], [0, 0], [l[4], l[5]], [n[4], n[5]]])

C1 = np.identity(4)

D1 = np.array([[0,0], [0,0], [0,0], [0,0]])




########## SYMMETRIC EQUATIONS OF MOTION IN STATE-SPACE FORM ##########

x = [(V0/c) * (CXu/(2 * muc)),
     (V0/c) * (CXa/(2 * muc)),
     (V0/c) * (CZ0/(2 * muc)),
     (V0/c) * (CXq/(2 * muc)),
     (V0/c) * (CXde/(2 * muc))]

z = [(V0/c) * (CZu/(2 * muc - CZadot)),
     (V0/c) * (CZa/(2 * muc - CZadot)),
     (V0/c) * (CX0/(2 * muc - CZadot)),
     (V0/c) * ((2 * muc + CZq)/(2 * muc - CZadot)),
     (V0/c) * (CZde/(2 * muc - CZadot))]

m = [(V0/c) * ((Cmu + CZu * ((Cmadot)/(2 * muc - CZadot)))/(2 * muc * KY2)),
     (V0/c) * ((Cma + CZa * ((Cmadot)/(2 * muc - CZadot)))/(2 * muc * KY2)),
     - (V0/c) * ((CX0 * ((Cmadot)/(2 * muc - CZadot)))/(2 * muc * KY2)),
     (V0/c) * ((Cmq + Cmadot * ((2 * muc + CZq)/(2 * muc - CZadot)))/(2 * muc * KY2)),
     (V0/c) * ((Cmde + CZde * ((Cmadot)/(2 * muc - CZadot)))/(2 * muc * KY2))]

A2 = np.array([[x[0], x[1], x[2], 0], [z[0], z[1], z[2], z[3]], [0, 0, 0, V0/c], [m[0], m[1], m[2], m[3]]])

B2 = np.array([[x[4]], [z[4]], [0], [m[4]]])

C2 = np.identity(4)

D2 = np.array([[0], [0], [0], [0]])

########## Eigenvalues and Eigenvectors ##########

[e2, v2] = np.linalg.eig(A2)
[e1, v1] = np.linalg.eig(A1)

########## Continuous time-state-space Model ##########



sys1 = ml.ss(A1, B1, C1, D1)

sys2 = ml.ss(A2, B2, C2, D2)


#ml.pzmap(sys1) #eigenvalues plot
#ml.pzmap(sys2)


########### Simulation of Eigen Motions########
dt= 0.1
t = np.arange(0 , 10, dt)

#Short Period

t_delta = int(0.1/dt) #time the elevator is deflected
s_def = -0.031 #angle of deflection[rad]
s_e = np.ones((t_delta,1))
s_en = np.array([i*s_def for i in s_e])
zeros = np.zeros((99,1))
s_e = np.vstack((s_en,zeros))

y , T, x = ml.lsim(sys2, s_e, t)


plt.plot(T, y[:,3])
plt.show()

# Phugoid
t = np.arange(0 , 7000, dt)

#getting data
delta_e_index = np.where(parameters=='delta_e')[0].flat[0]

flight_delta_e = np.array(flight_data[flight_headers[delta_e_index]])
value = flight_delta_e[32280:33400]
s_e_phu = np.array([flight_delta_e[32280:33400]])

s_e_phun = np.array([i*0.017455 for i in s_e_phu])
missing = len(t) - len(value)
zeros_phu = np.zeros((1,missing))
s_e_phu_n = np.hstack((s_e_phun,zeros_phu))
s_e_phu_n = np.transpose(s_e_phu_n)

#y , T, x = ml.lsim(sys2, s_e_phu_n, t)
#plt.plot(T[1:300], y[:,0][1:300])
#plt.show()

####Dutch Roll###########

u = [-0.130899,0] #deflection

T, y = ctr.impulse_response(sys2,t)
#plt.plot(T[1:300], y[1][1:300])

########## Print Commands ##########


#print(e1, v2)

#print(e2, v2)