import numpy as np
import matplotlib.pyplot as plt

a = np.linspace(0, 10, 1000)
x = np.cos(a)

plt.plot(a, x)
plt.grid()
plt.show()
