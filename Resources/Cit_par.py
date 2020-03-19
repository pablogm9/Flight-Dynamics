# Citation 550 - Linear simulation
from math import *
import sys
import Read
import numpy as np

#class Cit_par_class:

def Cit_par_Values(FlightType,Motion):              #FlightType: Reference Data = 1, Flight Test Data = 2
                                                    #Motion: Short Period = 1, Phugoid = 2, Dutch Roll = 3, Aperiodic Roll = 4, Spiral = 5
    # self.xcg = 0.25 * c
    rho = 1.225
    m = 12000
    Cma = -0.8
    CZ0 = 0.002
    Cl = 1.2

    # Stationary flight condition
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

    velocity_index = np.where(parameters == 'Dadc1_tas')[0].flat[0]
    hp0_index = np.where(parameters == 'Dadc1_bcAlt')[0].flat[0]
    th0_index = np.where(parameters == 'Ahrs1_Pitch')[0].flat[0]
    time_index = np.where(parameters=='time')[0].flat[0]

    time = np.array(reference_data[[reference_headers[time_index]]])

    if Motion == 1:
        if FlightType == 1:
            time_ini = 3635
            reference_V0 = np.array(reference_data[reference_headers[velocity_index]])
            reference_hp0 = np.array(reference_data[reference_headers[hp0_index]])
            reference_th0 = np.array(reference_data[reference_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0]
            hp0 = reference_hp0[cell_hp0]
            th0 = reference_th0[cell_th0]
        elif FlightType == 2:
            time_ini = 3158
            flight_V0 = np.array(flight_data[flight_headers[velocity_index]])
            flight_hp0 = np.array(flight_data[flight_headers[hp0_index]])
            flight_th0 = np.array(flight_data[flight_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0]
            hp0 = flight_hp0[cell_hp0]
            th0 = flight_th0[cell_th0]

    elif Motion == 2:
        if FlightType == 1:
            time_ini = 3237
            reference_V0 = np.array(reference_data[reference_headers[velocity_index]])
            reference_hp0 = np.array(reference_data[reference_headers[hp0_index]])
            reference_th0 = np.array(reference_data[reference_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0]
            hp0 = reference_hp0[cell_hp0]
            th0 = reference_th0[cell_th0]
        elif FlightType == 2:
            time_ini = 3230
            flight_V0 = np.array(flight_data[flight_headers[velocity_index]])
            flight_hp0 = np.array(flight_data[flight_headers[hp0_index]])
            flight_th0 = np.array(flight_data[flight_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0]
            hp0 = flight_hp0[cell_hp0]
            th0 = flight_th0[cell_th0]
    elif Motion == 3:
        if FlightType == 1:
            time_ini = 3717
            reference_V0 = np.array(reference_data[reference_headers[velocity_index]])
            reference_hp0 = np.array(reference_data[reference_headers[hp0_index]])
            reference_th0 = np.array(reference_data[reference_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0]
            hp0 = reference_hp0[cell_hp0]
            th0 = reference_th0[cell_th0]
        elif FlightType == 2:
            time_ini = 3479
            flight_V0 = np.array(flight_data[flight_headers[velocity_index]])
            flight_hp0 = np.array(flight_data[flight_headers[hp0_index]])
            flight_th0 = np.array(flight_data[flight_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0]
            hp0 = flight_hp0[cell_hp0]
            th0 = flight_th0[cell_th0]
    elif Motion == 4:
        if FlightType == 1:
            time_ini = 3550
            reference_V0 = np.array(reference_data[reference_headers[velocity_index]])
            reference_hp0 = np.array(reference_data[reference_headers[hp0_index]])
            reference_th0 = np.array(reference_data[reference_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0]
            hp0 = reference_hp0[cell_hp0]
            th0 = reference_th0[cell_th0]
        elif FlightType == 2:
            time_ini = 3607
            flight_V0 = np.array(flight_data[flight_headers[velocity_index]])
            flight_hp0 = np.array(flight_data[flight_headers[hp0_index]])
            flight_th0 = np.array(flight_data[flight_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0]
            hp0 = flight_hp0[cell_hp0]
            th0 = flight_th0[cell_th0]
    elif Motion == 5:
        if FlightType == 1:
            time_ini = 3920
            reference_V0 = np.array(reference_data[reference_headers[velocity_index]])
            reference_hp0 = np.array(reference_data[reference_headers[hp0_index]])
            reference_th0 = np.array(reference_data[reference_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = reference_V0[cell_V0]
            hp0 = reference_hp0[cell_hp0]
            th0 = reference_th0[cell_th0]
        elif FlightType == 2:
            time_ini = 3675
            flight_V0 = np.array(flight_data[flight_headers[velocity_index]])
            flight_hp0 = np.array(flight_data[flight_headers[hp0_index]])
            flight_th0 = np.array(flight_data[flight_headers[th0_index]])
            cell_V0 = np.where(time == time_ini)[0][0]
            cell_hp0 = np.where(time == time_ini)[0][0]
            cell_th0 = np.where(time == time_ini)[0][0]
            V0 = flight_V0[cell_V0]
            hp0 = flight_hp0[cell_hp0]
            th0 = flight_th0[cell_th0]
    else:
        print("Not a valid input")
        sys.exit()

    #hp0 = pressure altitude in the stationary flight condition [m]
    #V0 = true airspeed in the stationary flight condition [m/sec]
    #alpha0 = angle of attack in the stationary flight condition [rad]
    #th0    =  pitch angle in the stationary flight condition [rad]

    # Aircraft mass
    #m      =             # mass [kg]

    # aerodynamic properties
    #e      =             # Oswald factor [ ]
    #CD0    =             # Zero lift drag coefficient [ ]
    #CLa    =             # Slope of CL-alpha curve [ ]

    # Longitudinal stability
    if FlightType == 1:
        from StationaryData.stationary_reference_call import C_m_alpha, C_m_delta
        Cma    = C_m_alpha           # longitudinal stabilty [ ]
        Cmde   = C_m_delta      # elevator effectiveness [ ]
    elif FlightType == 2:
        from StationaryData.stationary_testflight_call import C_m_alpha, C_m_delta
        Cma = C_m_alpha
        Cmde = C_m_delta
    else:
        print("Not a valid input")
        sys.exit()

    # Aircraft geometry

    S      = 30.00	          # wing area [m^2]
    Sh     = 0.2 * S         # stabiliser area [m^2]
    Sh_S   = Sh / S	          # [ ]
    lh     = 0.71 * 5.968    # tail length [m]
    c      = 2.0569	          # mean aerodynamic cord [m]
    lh_c   = lh / c	          # [ ]
    b      = 15.911	          # wing span [m]
    bh     = 5.791	          # stabilser span [m]
    A      = b ** 2 / S      # wing aspect ratio [ ]
    Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
    Vh_V   = 1	          # [ ]
    ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

    # Constant values concerning atmosphere and gravity

    rho0   = 1.2250          # air density at sea level [kg/m^3]
    lamda = -0.0065         # temperature gradient in ISA [K/m]
    Temp0  = 288.15          # temperature at sea level in ISA [K]
    R      = 287.05          # specific gas constant [m^2/sec^2K]
    g      = 9.81            # [m/sec^2] (gravity constant)

    # air density [kg/m^3]
    #rho    = rho0 **( ((1+(lamda * hp0 / Temp0))), (-((g / (lamda*R)) + 1)))
    W      = m * g            # [N]       (aircraft weight)

    # Constant values concerning aircraft inertia

    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
    KX2    = 0.019
    KZ2    = 0.042
    KXZ    = 0.002
    KY2    = 1.25 * 1.114

    # Aerodynamic constants
    Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
    #CNwa   = CLa                    # Wing normal force slope [ ]
    CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
    depsda = 4 / (A + 2)            # Downwash gradient [ ]

    # Lift and drag coefficient

    #CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
    #CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

    # Stabiblity derivatives

    CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CXu    = -0.02792
    CXa    = +0.47966		# Positive! (has been erroneously negative since 1993)
    CXadot = +0.08330
    CXq    = -0.28170
    CXde   = -0.03728

    CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu    = -0.37616
    CZa    = -5.74340
    CZadot = -0.00350
    CZq    = -5.66290
    CZde   = -0.69612

    Cmu    = +0.06990
    Cmadot = +0.17800
    Cmq    = -8.79415

    CYb    = -0.7500
    CYbdot =  0
    CYp    = -0.0304
    CYr    = +0.8495
    CYda   = -0.0400
    CYdr   = +0.2300

    Clb    = -0.10260
    Clp    = -0.71085
    Clr    = +0.23760
    Clda   = -0.23088
    Cldr   = +0.03440

    Cnb    =  +0.1348
    Cnbdot =   0
    Cnp    =  -0.0602
    Cnr    =  -0.2061
    Cnda   =  -0.0120
    Cndr   =  -0.0939

    return rho, m, Cma, CZ0, Cl, hp0, V0, th0, Cmde, S, Sh, Sh_S, lh, c, lh_c, b, bh, A, Ah, Vh_V, ih, rho0, lamda, Temp0, R, g, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNha, depsda, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb, CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr