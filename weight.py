import numpy as np


BEM = 13600*0.453592 #Basic Empty Mass  [kg]
#Payload = np.sum(payload_array) #Total Payload mass in kg [kg]
Payload = np.sum([102,80,79,79,69,69,78,80,96])
#BlockFuel = Block_fuel #Block Fuel mass in kg [kg]
BlockFuel= 2438*0.454592
g=9.81 #gravitaional acceleration

#F_used= f_used_data
F_used=np.array([261,308,336,352,377,403,420])
rampmass=(BEM+Payload+BlockFuel)*g


W = np.ones(7)

for i in range(len(W)):
    W[i]=rampmass - F_used[i]
