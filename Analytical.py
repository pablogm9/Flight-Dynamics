from Resources.Cit_par import *
import numpy as np

#--------------------SYMMETRIC-EoM---------------------------------

#------------SHORT-PERIOD----------------------------
A_sp = 2 * muc * KY2**2 * (2 * muc - CZadot)
B_sp = - 2 * muc * KY2**2 * CZa - (2 * muc + CZq) * Cmadot - (2 * muc - CZadot) * Cmq
C_sp = CZa * Cmq - (2 * muc + CZq) * Cma

Eig_val_short_period_1 = (-B_sp + sqrt(4 * A_sp * C_sp - B_sp**2)) / (2 * A_sp)
Eig_val_short_period_2 = (-B_sp - sqrt(4 * A_sp * C_sp - B_sp**2)) / (2 * A_sp)

#--------------------PHUGOID-------------------------
A_ph = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
B_ph = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq * (CZu * CXa - CXu * CZa)
C_ph = CZ0 * (Cmu * Cza - CZu * Cma)

Eig_val_phugoid_1 = (-B_ph + sqrt(4 * A_ph * C_ph - B_ph**2)) / (2 * A_ph)
Eig_val_phugoid_2 = (-B_ph - sqrt(4 * A_ph * C_ph - B_ph**2)) / (2 * A_ph)


#Calculation of eigenvectors and values for symmetric and asymmetric
#eigenmotions


#List of available parameters for aymmetric flight


