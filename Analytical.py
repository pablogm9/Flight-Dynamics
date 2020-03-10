from Resources.Cit_par import *
from math import *

Cl=1
import numpy as np



#empy lists to fill with values
lambdas     = []
lambdas_real = []
lambdas_imag = []
t_half      = []
taus        = []
periods     = []
no_per_half = []
log_decr    = []
damp_ratio  = []


#========================SYMMETRIC-MOTIONS=============================

#--------------------------SHORT-PERIOD----------------------------
A_sp = 2 * muc * KY2**2 * (2 * muc - CZadot)
B_sp = - 2 * muc * KY2**2 * CZa - (2 * muc + CZq) * Cmadot - (2 * muc - CZadot) * Cmq
C_sp = CZa * Cmq - (2 * muc + CZq) * Cma

Eig_val_short_period_real = (-B_sp) / (2 * A_sp)
Eig_val_short_period_imag = sqrt(4 * A_sp * C_sp - B_sp**2) / (2 * A_sp)
Eig_val_short_period_1 = complex(Eig_val_short_period_real, Eig_val_short_period_imag)
Eig_val_short_period_2 = complex(Eig_val_short_period_real, -1*Eig_val_short_period_imag)
lambdas.append(Eig_val_short_period_1)
lambdas.append(Eig_val_short_period_2)
lambdas_real.append(Eig_val_short_period_real)
lambdas_imag.append(Eig_val_short_period_imag)



#--------------------------PHUGOID---------------------------------------------------------
A_ph = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
B_ph = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq * (CZu * CXa - CXu * CZa)
C_ph = CZ0 * (Cmu * CZa - CZu * Cma)

Eig_val_phugoid_real = (-B_ph) / (2 * A_ph)
Eig_val_phugoid_img = cmath.sqrt(4 * A_ph * C_ph - B_ph**2) / (2 * A_ph)
Eig_val_phugoid_1 = complex(Eig_val_phugoid_real, Eig_val_phugoid_img)
Eig_val_phugoid_2 = complex(Eig_val_phugoid_real, -1*Eig_val_phugoid_img)
lambdas.append(Eig_val_phugoid_1)
lambdas.append(Eig_val_phugoid_2)
lambdas_real.append(Eig_val_phugoid_real)
lambdas_imag.append(Eig_val_phugoid_img)




#============================ASYMMETRIC- MOTIONS========================================



#--------------------Aperiodic roll----------------------
lambda_aroll = Clp/(4*mub*KX2)
print('Eigenvalue corresponding to aperiodic roll=', lambda_aroll)
lambdas.append(lambda_aroll)
lambdas_real.append(lambda_aroll)



#---------------------Dutch roll--------------------------
#lambda_droll_1=(2*(Cnr+2*KZ2*CYb)+ sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)**2))/(16*mub*KZ2)
#lambda_droll_2= (2*(Cnr+2*KZ2*CYb)- sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)**2))/(16*mub*KZ2)

#lambdas.append(lambda_droll_1)
#lambdas.append(lambda_droll_2)
A_droll= 8*mub**2*KZ2
B_droll = -2*mub*(Cnr+ 2*KZ2*CYb)
C_droll = 4*mub*Cnb+CYb*Cnr
disc_droll = B_droll**2-4*A_droll*C_droll
lambda_droll_real = -B_droll/(2*A_droll)
lambda_droll_imag = cmath.sqrt(disc_droll)/(2*A_droll)
print('Eigenvalues corresponding to dutch roll=', [lambda_droll_real+lambda_droll_imag])
lambdas.append(lambda_droll_real+lambda_droll_imag)
lambdas_real.append(lambda_droll_real)
#----------------------Spiral motion----------------------------------
lambda_spiral= (2*Cl*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))
print('Eigenvalues corresponding to spiral motion=', lambda_spiral)
lambdas.append(lambda_spiral)
lambdas_real.append(lambda_spiral)
#--------------Time to half the amplitude---------------


def tau(y):
    return -1*c/(Vh_V*y)
#def per(z):
 #   return 2*pi*c/()

for i in lambdas:
    taus.append(tau(i))
print('Times to reach half amplitude=',t_half)
print('The time constants=', taus)

# #----------------------Characteristics------------------------------------------------------------


#---------For short period------------

#------------Half time-----------------
def thalf(x):
    return log(0.5, e) * c / (x * Vh_V) ### The velocity needs to be changed for every eigenmotion, the speed is taken at the beginning of the motion
for i in lambdas_real:
    t_half.append(thalf(i))

#----------Period---------------------








