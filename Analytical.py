from Resources.Cit_par import *
import math as m
import numpy as np


#Calculation of eigenvectors and values for symmetric and asymmetric
#eigenmotions

#List of available parameters for aymmetric flight
# to use a variable from Cit_par > print(Cit_par.CXu)

#Aperiodic roll
lambda_aroll = Clp/(4*mub*KX2)
#Dutch roll
lambda_droll_1=(2*(Cnr+2*KZ2*CYb)+ m.sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)**2))/(16*mub*KZ2)
lambda_droll_2= (2*(Cnr+2*KZ2*CYb)- m.sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)**2))/(16*mub*KZ2)
#Spiral motion
lambda_spiral= (2*Cl*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))





