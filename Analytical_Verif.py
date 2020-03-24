from Resources.Cit_par import *
from math import *
import cmath
#import matplotlib.pyplot as plt

rho, m, Cma, CZ0, Cl, hp0, V0, th0, Cmde, S, Sh, Sh_S, lh, c, lh_c, b, bh, A, Ah, Vh_V, ih, rho0, lamda, Temp0, R, g, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNha, depsda, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb, CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr= Cit_par_Values(1,1)



#empty lists to fill with values
lambdas     = []
lambdas_real = []
lambdas_imag = []
t_half_sym      = []
t_half_asym      = []
taus        = []
periods     = []
no_per_half = []
log_decr    = []
damp_ratio  = []


#========================SYMMETRIC-MOTIONS=============================

#--------------------------SHORT-PERIOD----------------------------
A_sp = 2 * muc * KY2 * (2 * muc - CZadot)
B_sp = - 2 * muc * KY2 * CZa - (2 * muc + CZq) * Cmadot - (2 * muc - CZadot) * Cmq
C_sp = CZa * Cmq - (2 * muc + CZq) * Cma

Eig_val_short_period_real = (-B_sp) / (2 * A_sp)
Eig_val_short_period_imag = cmath.sqrt(- 4 * A_sp * C_sp + B_sp**2) / (2 * A_sp)
Eig_val_short_period_1 = (Eig_val_short_period_real + Eig_val_short_period_imag) * V0/c
Eig_val_short_period_2 = (Eig_val_short_period_real - Eig_val_short_period_imag) * V0/c

lambdas.append(Eig_val_short_period_1)
lambdas.append(Eig_val_short_period_2)
lambdas_real.append(Eig_val_short_period_real)
lambdas_imag.append(abs(Eig_val_short_period_imag))
print('Eigenvalues corresponding to short period=', Eig_val_short_period_1, Eig_val_short_period_2)


#--------------------------PHUGOID---------------------------------------------------------
A_ph = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
B_ph = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq * (CZu * CXa - CXu * CZa)
C_ph = CZ0 * (Cmu * CZa - CZu * Cma)
Eig_val_phugoid_real = (-B_ph) / (2 * A_ph)
Eig_val_phugoid_img = cmath.sqrt(- 4 * A_ph * C_ph + B_ph**2) / (2 * A_ph)
Eig_val_phugoid_1 = (Eig_val_phugoid_real + Eig_val_phugoid_img) * V0/c
Eig_val_phugoid_2 = (Eig_val_phugoid_real - Eig_val_phugoid_img) * V0/c

lambdas.append(Eig_val_phugoid_1)
lambdas.append(Eig_val_phugoid_2)
lambdas_real.append(Eig_val_phugoid_real)
lambdas_imag.append(abs(Eig_val_phugoid_img))
print('Eigenvalues corresponding to phugoid=', Eig_val_phugoid_1, Eig_val_phugoid_2)




#============================ASYMMETRIC- MOTIONS========================================



#--------------------Aperiodic roll----------------------
lambda_aroll = Clp/(4*mub*KX2)* V0/b
print('Eigenvalue corresponding to aperiodic roll=', lambda_aroll)
lambdas.append(lambda_aroll)
lambdas_real.append(lambda_aroll)



#---------------------Dutch roll--------------------------
A_droll=  8 * mub**2 * KZ2
B_droll = -2 * mub * (Cnr + 2 * KZ2 * CYb)
C_droll = 4 * mub * Cnb + CYb * Cnr
disc_droll = B_droll**2 - 4 * A_droll * C_droll
lambda_droll_real = - B_droll / (2 * A_droll)
lambda_droll_imag = cmath.sqrt(disc_droll)/(2*A_droll)
Eig_val_droll_1 = (lambda_droll_real + lambda_droll_imag) * V0/b
Eig_val_droll_2 = (lambda_droll_real - lambda_droll_imag) * V0/b
print('Eigenvalues corresponding to dutch roll=', Eig_val_droll_1, Eig_val_droll_2)
lambdas.append(Eig_val_droll_1)
lambdas.append(Eig_val_droll_2)
lambdas_real.append(lambda_droll_real)
lambdas_imag.append(abs(lambda_droll_imag))




#----------------------Spiral motion----------------------------------
lambda_spiral= ((2 * Cl * (Clb * Cnr - Cnb * Clr))/(Clp * (CYb * Cnr + 4*mub * Cnb) - Cnp * (CYb * Clr + 4*mub * Clb))) * V0/b
disc_spiral = Cl * (Clb * Cnr - Cnb * Clr)
print('Eigenvalues corresponding to spiral motion=', lambda_spiral)
lambdas.append(lambda_spiral)
lambdas_real.append(lambda_spiral)


x = [i.real for i in lambdas]
y = [i.imag for i in lambdas]
# plt.scatter(x, y)
# plt.show()



#--------------Time to half the amplitude---------------

#for symmetric motions
def tau(y):
    return -1*c/(Vh_V*y)
#def per(z):
 #   return 2*pi*c/()

for i in lambdas:
    taus.append(tau(i))

# print('The time constants=', taus)

#----------------------Characteristics-------------------------

#------------Half time-----------------
#applies to all eigenvalues
def thalf_sym(x):
    return log(0.5, e) * c / (x * Vh_V) ### The velocity needs to be changed for every eigenmotion, the speed is taken at the beginning of the motion
def thalf_asym(x):
    return log(0.5, e) * b / (x * Vh_V) ### The velocity needs to be changed for every eigenmotion, the speed is taken at the beginning of the motion



for k in range(2):
        t_half_sym.append(thalf_sym(lambdas_real[k]))
for k1 in range(2,5):
        t_half_asym.append(thalf_asym(lambdas_real[k1]))
# print('Times to reach half amplitude for asymmetric motions=',t_half_asym)
#----------Period---------------------
#applies to complex eigvalues only

def per(z):
    return 2*pi*c/(z*Vh_V)#again the velocity needs to be appended
for j in lambdas_imag:
    periods.append(per(j))

#-----------Time constant------------------
#only applies to real eigenvalues
tau_aperiod = -1*c/lambda_aroll
tau_spiral = -1*c/lambda_spiral
taus.append(tau_aperiod)
taus.append(tau_spiral)

#----------C half values--------------
#only applies to motions that have an imaginary part
c_half_shawty = t_half_sym[0]/periods[0]
c_half_phugoid= t_half_sym[1]/periods[1]
c_half_droll = t_half_asym[1]/periods[2]
no_per_half.append(c_half_shawty)
no_per_half.append(c_half_phugoid)
no_per_half.append(c_half_droll)

#-----------Logarithmic decrement---------
#only applies to motions that have an imaginary part
dec_shawty = lambdas_real[0]*Vh_V*periods[0]/c
dec_phugoid = lambdas_real[1]*Vh_V*periods[1]/c
dec_droll = lambdas_real[3]*Vh_V*periods[2]/c
log_decr.append(dec_shawty)
log_decr.append(dec_phugoid)
log_decr.append(dec_droll)

#-----------Damping ratio------------------
#only applies to motions that have an imaginary part
damp_shawty = -1*lambdas[0].real/abs(lambdas[0])
damp_phugoid = -1*lambdas[2].real/abs(lambdas[2])
damp_droll = -1*lambda_droll_real/abs(lambda_droll_real+lambda_droll_imag)
damp_ratio.append(damp_shawty)
damp_ratio.append(damp_phugoid)
damp_ratio.append(damp_droll)



#============================COMPLETE SOLUTION SYMMETRICAL======================
A_s = 4 * muc**2 * KY2 * (CZadot - 2*muc)
B_s = Cmadot * 2 * muc * (CZq + 2 * muc) - Cmq * 2 * muc * (CZadot - 2 * muc) - 2*muc * KY2 * (CXu * (CZadot - 2*muc) - 2*muc*CZa)
C_s = Cma * 2*muc * (CZq + 2*muc) - Cmadot * (2*muc * CX0 + CXu * (CZq + 2*muc)) + Cmq * (CXu * (CZadot - 2*muc) - 2*muc * CZa) + 2*muc * KY2 * (CXa * CZu - CZa * CXu)
D_s = Cmu * ( CXa * ( CZq + 2*muc) - CZ0 * (CZadot - 2*muc)) - Cma * (2*muc * CX0 + CXu * (CZq + 2*muc)) + Cmadot * (CX0 * CXu - CZ0 * CZu) + Cmq *  (CXu * CZa - CZu * CXa)
E_s = -Cmu * (CX0 * CXa + CZ0 * CZa) + Cma * (CX0 * CXu + CZ0 * CZu)
p_s = np.array([A_s, B_s, C_s, D_s, E_s])

roots_s = np.roots(p_s)
root1_s = roots_s[0]*V0/c
root2_s = roots_s[1]*V0/c
root3_s = roots_s[2]*V0/c
root4_s = roots_s[3]*V0/c


#===========================COMPLETE SOLUTION ASYMMETRICAL==========================
A_a = 16 * mub**3 * (KX2*KZ2 - KXZ**2)
B_a = -4*mub**2 * ( 2*CYb * (KX2*KZ2 - KXZ**2) + Cnr * KX2 + Clp * KZ2 + (Clr + Cnp) * KXZ)
C_a = 2*mub * ((CYb * Cnr - CYr * Cnb) * KX2 + (CYb * Clp - Clb * CYp) * KZ2 + ((CYb * Cnp - Cnb * CYp) + (CYb * Clr - Clb * CYr)) * KXZ + 4*mub * Cnb * KX2 + 4*mub * Clb * KXZ + 0.5 * (Clp * Cnr - Cnp * Clr))
D_a = -4*mub * Cl * (Clb * KZ2 + Cnb * KXZ) + 2*mub * (Clb * Cnp - Cnb * Clp) + 0.5 * CYb * (Clr * Cnp - Cnr * Clp) + 0.5 * CYp * (Clb * Cnr - Cnb * Clr) + 0.5 * CYr * (Clp * Cnb - Cnp * Clb)
E_a = Cl * (Clb * Cnr - Cnb * Clr)
p_a = np.array([A_a, B_a, C_a, D_a, E_a])

roots_a = np.roots(p_a)
root1_a = roots_a[0]*V0/b
root2_a = roots_a[1]*V0/b
root3_a = roots_a[2]*V0/b
root4_a = roots_a[3]*V0/b
print(root1_a, root2_a, root3_a, root4_a)