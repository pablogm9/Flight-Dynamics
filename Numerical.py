from Resources.Cit_par import Cit_par_Values
import scipy as np
import control.matlab as ml
import matplotlib.pyplot as plt
import control as ctr
import Read
Cl = 1.3

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
flight_data,flight_headers,flight_descriptions = Read.get_data('testflight')

delta_e_index = np.where(parameters=='delta_e')[0].flat[0]
delta_r_index = np.where(parameters=='delta_r')[0].flat[0]
delta_a_index = np.where(parameters=='delta_a')[0].flat[0]
time_index = np.where(parameters=='time')[0].flat[0]

reference_delta_e = np.array(reference_data[reference_headers[delta_e_index]])
reference_delta_r = np.array(reference_data[reference_headers[delta_r_index]])
reference_delta_a = np.array(reference_data[reference_headers[delta_a_index]])
time = np.array(reference_data[[reference_headers[time_index]]])

flight_delta_e = np.array(flight_data[reference_headers[delta_e_index]])
flight_delta_r = np.array(reference_data[reference_headers[delta_r_index]])
flight_delta_a = np.array(reference_data[reference_headers[delta_a_index]])
time = np.array(reference_data[[reference_headers[time_index]]])

#class ABCDmat():
    #from Resources.Cit_par import Cit_par_class
    #def __init__(self):
        #Cit_par_Values(self)

def ABCD(flight,motion):
    
    rho, m, Cma, CZ0, Cl, hp0, V0, th0, Cmde, S, Sh, Sh_S, lh, c, lh_c, b, bh, A, Ah, Vh_V, ih, rho0, lamda, Temp0, R, g, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNha, depsda, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb, CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr = Cit_par_Values(flight,motion)

    ########## ASYMMETRIC EQUATIONS OF MOTION IN STATE-SPACE FORM ##########
    A1 = np.array([[(V0/b)*(CYb/2/mub), (V0/b)*(Cl/2/mub), (V0/b)*(CYp/2/mub),(V0/b)*(CYr - 4*mub)/2/mub],
                   [   0                     ,     0           ,        2*V0/b      ,           0   ],
                   [(V0/b)*((Clb*KZ2+Cnb*KXZ)/(4*mub*(KX2*KZ2-KXZ**2))), 0 , (V0/b)*((Clp*KZ2+Cnp*KXZ)/(4*mub*(KX2*KZ2-KXZ**2))) , (V0/b)*((Clr*KZ2+Cnr*KXZ)/(4*mub*(KX2*KZ2-KXZ**2)))],
                   [(V0/b)*((Clb*KXZ+Cnb*KX2)/(4*mub*(KX2*KZ2-KXZ**2))), 0 , (V0/b)*((Clp*KXZ+Cnp*KX2)/(4*mub*(KX2*KZ2-KXZ**2))) , (V0/b)*((Clr*KXZ+Cnr*KX2)/(4*mub*(KX2*KZ2-KXZ**2)))]])

    B1 = np.array([[0          ,    (V0/b) * (CYdr/(2*mub))],
                   [0          ,             0],
                   [(V0/b)*((Clda*KZ2+Cnda*KXZ)/(4*mub*(KX2*KZ2-KXZ**2)))  ,  (V0/b)*((Cldr*KZ2+Cndr*KXZ)/(4*mub*(KX2*KZ2-KXZ**2)))],
                   [(V0/b)*((Clda*KXZ+Cnda*KX2)/(4*mub*(KX2*KZ2-KXZ**2)))  ,  (V0/b)*((Cldr*KXZ+Cndr*KX2)/(4*mub*(KX2*KZ2-KXZ**2)))]])


    C1 = np.identity(4)

    D1 = np.array([[0, 0], [0, 0], [0, 0], [0, 0]])


    ########## SYMMETRIC EQUATIONS OF MOTION IN STATE-SPACE FORM ##########
    A2 = np.array([[(V0/c)*CXu/(2*muc)   , (V0/c)*CXa/(2*muc) , (V0/c)*CZ0/(2*muc) , (V0/c)*CXq/(2*muc)],
                    [(V0/c)*CZu/(2*muc-CZadot)   , (V0/c)*CZa/(2*muc-CZadot)  , -(V0/c)*CX0/(2*muc-CZadot) , (V0/c)*(2*muc + CZq)/(2*muc-CZadot)],
                   [        0                  ,                0           ,           0               ,           V0/c                     ],
                   [(V0/c)*(Cmu+CZu*(Cmadot/(2*muc-CZadot)))/(2*muc*KY2), (V0/c)*(Cma + CZa*(Cmadot/(2*muc-CZadot)))/(2*muc*KY2), -(V0/c)*CX0*(Cmadot/(2*muc-CZadot))/(2*muc*KY2) , (V0/c)* (Cmq + Cmadot*((2*muc + CZq)/(2*muc - CZadot)))/(2*muc*KY2)]])

    B2 = np.array([[(V0/c)*CXde/(2*muc)],
                 [(V0/c)*CZde/(2*muc -CZadot)],
                  [0],
                  [(V0/c)*(Cmde + CZde*(Cmadot/(2*muc-CZadot)))/(2*muc*KY2)]])

    C2 = np.identity(4)

    D2 = np.array([[0], [0], [0], [0]])

    #print(A2)
    #print(B2)
    #print(C2)
    #print(D2)


    return A1, B1, C1, D1, A2, B2, C2, D2



A1, B1, C1, D1, A2, B2, C2, D2 = ABCD()


########## EIGENVALUES ##########

[e1, v1] = np.linalg.eig(A1)
[e2, v2] = np.linalg.eig(A2)

print(e1)
print(e2)

########## STATE SPACE MODEL ##########

sys1 = ml.ss(A1, B1, C1, D1)

sys2 = ml.ss(A2, B2, C2, D2)

a = [-507.87304253328404, -1.1211704696917422, -5468.616677123462, -468.2605031612816, 0.17942559940696895]
#ml.pzmap(sys1) #eigenvalues plot
#ml.pzmap(sys2)


########### SIMULATION OF EIGENMOTIONS########
dt= 0.1


#Short Period
t1 = np.arange(0 , 10, dt)
t_ini = 3634
t_fin = 3636

u1 = inputcr(reference_delta_e, time, t1, t_ini , t_fin)

y1_r , T1_r, x1_r = ml.lsim(sys2, u1, t1)


u1_f = inputcr(flight_delta_e, time, t1, 3158 , 3160)

y1_f , T1_f, x1_f = ml.lsim(sys2, u1_f, t1)

# Phugoid
t2 = np.arange(0 , 1000, dt)

u2 = inputcr(reference_delta_e, time, t2, 3237, 3247)

y2_r , T2_r, x2_r = ml.lsim(sys2, u2, t2)


u2_f = inputcr(flight_delta_e, time, t2, 3230 , 3240)

y2_f , T2_f, x2_f = ml.lsim(sys2, u2_f, t2)

#Dutch Roll
t3 = np.arange(0 , 100, dt)

u3_t = inputcr(reference_delta_r, time, t3, 3717, 3718.8)

u3 = np.hstack((np.zeros((len(t3), 1)), u3_t))

y3_r , T3_r, x3_r = ml.lsim(sys1, u3, t3)


u3_t_f = inputcr(flight_delta_r, time, t3, 3479, 3480.1)

u3_f = np.hstack((np.zeros((len(t3), 1)), u3_t_f))

y3_f , T3_f, x3_f = ml.lsim(sys1, u3_f, t3)

#Aperiodic roll
t4 = np.arange(0 , 12, dt)

u4_t = inputcr(reference_delta_a, time, t4, 3550, 3551)

u4 = np.hstack((u4_t, np.zeros((len(t4), 1))))

y4_r , T4_r, x4_r = ml.lsim(sys1, u4, t4)


u4_t_f = inputcr(flight_delta_a, time, t4, 3607, 3608)

u4_f = np.hstack((u4_t, np.zeros((len(t4), 1))))

y4_f , T4_f, x4_f = ml.lsim(sys1, u4_f, t4)


#Spiral
t5 = np.arange(0, 100,dt)

u5_t = inputcr(reference_delta_a, time, t5, 3912, 3920)

u5 = np.hstack((u5_t, np.zeros((len(t5), 1))))

y5_r , T5_r, x5_r = ml.lsim(sys1, u5, t5)

u5_t_f = inputcr(flight_delta_a, time, t5, 3675, 3682)

u5_f = np.hstack((u5_t, np.zeros((len(t5), 1))))

y5_f , T5_f, x5_f = ml.lsim(sys1, u5_f, t5)



##########PLOTS OF RESPONSES###########
#Short Period
f1 = plt.figure(1)
plt.plot(T1_r,y1_r[:,0],'r',label='Reference Data')
plt.plot(T1_f,y1_f[:,0],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Absoulte Velocity')
plt.title('Short Period Response - velocity',pad=10)

f2 = plt.figure(2)
plt.plot(T1_r,y1_r[:,1],'r',label='Reference Data')
plt.plot(T1_f,y1_f[:,1],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Angle of Attack[-]')
plt.title('Short Period Response - Angle of Attack',pad=10)

f3 = plt.figure(3)
plt.plot(T1_r,y1_r[:,2],'r',label='Reference Data')
plt.plot(T1_f,y1_f[:,2],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Theta[-]')
plt.title('Short Period Response - Theta',pad=10)

f4 = plt.figure(4)
plt.plot(T1_r,y1_r[:,3],'r',label='Reference Data')
plt.plot(T1_f,y1_f[:,3],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('qc/V[-]')
plt.title('Short Period Response - qc/V',pad=10)


#Phugoid
f5 = plt.figure(5)
plt.plot(T2_r,y2_r[:,0],'r',label='Reference Data')
plt.plot(T2_f,y2_f[:,0],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Absoulte Velocity')
plt.title('Phugoid Response - velocity',pad=10)

f6 = plt.figure(6)
plt.plot(T2_r,y2_r[:,1],'r',label='Reference Data')
plt.plot(T2_f,y2_f[:,1],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Angle of Attack[-]')
plt.title('Phugoid Response - Angle of Attack',pad=10)

f7 = plt.figure(7)
plt.plot(T2_r,y2_r[:,2],'r',label='Reference Data')
plt.plot(T2_f,y2_f[:,2],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Theta[-]')
plt.title('Phugoid Response- Theta',pad=10)

f8 = plt.figure(8)
plt.plot(T2_r,y2_r[:,3],'r',label='Reference Data')
plt.plot(T2_f,y2_f[:,3],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('qc/V[-]')
plt.title('Phugoid Response - qc/V',pad=10)

#Dutch Roll
f9 = plt.figure(9)
plt.plot(T3_r,y3_r[:,0],'r',label='Reference Data')
plt.plot(T3_f,y3_f[:,0],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Side Slip [-]')
plt.title('Dutch Roll Repsonse - Sideslip',pad=10)

f10 = plt.figure(10)
plt.plot(T3_r,y3_r[:,1],'r',label='Reference Data')
plt.plot(T3_f,y3_f[:,1],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Roll Angle [-]')
plt.title('Dutch Roll Repsonse - Roll Angle',pad=10)

f11 = plt.figure(11)
plt.plot(T3_r,y3_r[:,2],'r',label='Reference Data')
plt.plot(T3_f,y3_f[:,2],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('pb/V [-]')
plt.title('Dutch Roll Repsonse - pb/V',pad=10)

f12 = plt.figure(12)
plt.plot(T3_r,y3_r[:,3],'r',label='Reference Data')
plt.plot(T3_f,y3_f[:,3],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('rb/V [-]')
plt.title('Dutch Roll Repsonse - rb/V',pad=10)

#Aperiodic Roll
f13 = plt.figure(13)
plt.plot(T4_r,y4_r[:,0],'r',label='Reference Data')
plt.plot(T4_f,y4_f[:,0],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Side Slip [-]')
plt.title('Aperiodic Roll Repsonse - Sideslip',pad=10)

f14 = plt.figure(14)
plt.plot(T4_r,y4_r[:,1],'r',label='Reference Data')
plt.plot(T4_f,y4_f[:,1],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Roll Angle [-]')
plt.title('Aperiodic Roll Repsonse - Roll Angle',pad=10)

f15 = plt.figure(15)
plt.plot(T4_r,y4_r[:,2],'r',label='Reference Data')
plt.plot(T4_f,y4_f[:,2],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('pb/V [-]')
plt.title('Aperiodic Roll Repsonse - pb/V',pad=10)

f16 = plt.figure(16)
plt.plot(T4_r,y4_r[:,3],'r',label='Reference Data')
plt.plot(T4_f,y4_f[:,3],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('rb/V [-]')
plt.title('Aperiodic Roll Repsonse - rb/V',pad=10)

#Spiral

f17 = plt.figure(17)
plt.plot(T5_r,y5_r[:,0],'r',label='Reference Data')
plt.plot(T5_f,y5_f[:,0],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Side Slip [-]')
plt.title('Spiral Repsonse - Sideslip',pad=10)

f18 = plt.figure(18)
plt.plot(T5_r,y5_r[:,1],'r',label='Reference Data')
plt.plot(T5_f,y5_f[:,1],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('Roll Angle [-]')
plt.title('Spiral Repsonse - Roll Angle',pad=10)

f19 = plt.figure(19)
plt.plot(T5_r,y5_r[:,2],'r',label='Reference Data')
plt.plot(T5_f,y5_f[:,2],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('pb/V [-]')
plt.title('Spiral Repsonse - pb/V',pad=10)

f20 = plt.figure(20)
plt.plot(T5_r,y5_r[:,3],'r',label='Reference Data')
plt.plot(T5_f,y5_f[:,3],'b',label='Flight Data')

plt.legend()
plt.grid()
plt.xlabel('Time [sec]')
plt.ylabel('rb/V [-]')
plt.title('Spiral Repsonse - rb/V',pad=10)



########## Print Commands ##########

plt.show()

#print(e1, v2)

#print(e2, v2)