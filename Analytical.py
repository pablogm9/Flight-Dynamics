from Resources.Cit_par import *
from math import *

Cl=1
import numpy as np



#Calculation of eigenvectors and values for symmetric and asymmetric
#eigenmotions

#List of available parameters for aymmetric flight
lambdas=[]
t_half=[]
taus= []
periods = []

#--------------------Aperiodic roll----------------------
lambda_aroll = Clp/(4*mub*KX2)
print('Eigenvalue corresponding to aperiodic roll=', lambda_aroll)
lambdas.append(lambda_aroll)



#---------------------Dutch roll--------------------------
lambda_droll_1=(2*(Cnr+2*KZ2*CYb)+ sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)**2))/(16*mub*KZ2)
lambda_droll_2= (2*(Cnr+2*KZ2*CYb)- sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)**2))/(16*mub*KZ2)
print('Eigenvalues corresponding to dutch roll=', [lambda_droll_1,lambda_droll_2])
lambdas.append(lambda_droll_1)
lambdas.append(lambda_droll_2)
#----------------------Spiral motion----------------------------------
lambda_spiral= (2*Cl*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))
print('Eigenvalues corresponding to spiral motion=', lambda_spiral)
lambdas.append(lambda_spiral)
#--------------Time to half the amplitude---------------
def thalf(x):
    return log(0.5, e) * c / (x * Vh_V)
print(lambdas)
def tau(y):
    return -1*c/(Vh_V*y)
def per(z):
    return 2*pi*c/()
for i in lambdas:
    t_half.append(thalf(i))
for i in lambdas:
    taus.append(tau(i))
print('Times to reach half amplitude=',t_half)
print('The time constants=', taus)






