from Resources.Cit_par import *
import numpy as np

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

D1 = np.array([[0], [0], [0], [0]])


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

#print(y)
#print(l)
#print(n)
#print(A)
#print(B)