#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 21:40:44 2020

@author: pablo
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Read



# Getting reference data
reference_data,reference_headers,reference_descriptions = Read.get_data('ref_data')


# Sample plot, AOA vs t
time = np.array(reference_data[[reference_headers[-1]]])
alpha = np.array(reference_data[[reference_headers[0]]])
plt.plot(time,alpha)
plt.xlabel('Time [sec]')
plt.ylabel('Angle of Attack [deg]')
plt.title('Angle of Attack over time',pad=10)
plt.show()
