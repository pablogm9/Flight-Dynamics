#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:58:44 2020

@author: pablo
"""

import numpy as np
import pandas as pd
from scipy.io import loadmat


def get_data(filename):
    
    '''
    Reads data from the provided MATLAB files. The files must be 
    located in a folder called Data inside the current working directory.
    
    INPUT:
        - file_name: Name of the file WITHOUT the file extension.
        
    OUTPUTS:
        - flightdata_df: Pandas dataframe containing the data. The header 
                         of each column is the name of the corresponding parameter
                         and the corresponding unit for that parameter, as
                         provided in the MATLAB file.
                         
        - headers: Numpy array containing the name of each of the columns in the 
                   dataframe. Can be used to select specific columns from the
                   dataframe.
                         
        - descriptions: Numpy array containing the description of each of the 
                        parameters. Included just in case you need to check
                        what a specific parameter actually is.
    
    '''
    
    # Extract data from .mat file
    mat = loadmat('Data/'+str(filename)+'.mat')
    flight_data = mat['flightdata']
    
    # List of parameters to loop through
    parameters =  ['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt', 'rh_engine_itt', 'lh_engine_OP', 'rh_engine_OP', 'lh_engine_fan_N1', 'lh_engine_turbine_N2', 'rh_engine_fan_N1', 'rh_engine_turbine_N2',
                'lh_engine_FU', 'rh_engine_FU', 'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll', 'Ahrs1_Pitch', 'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate', 'Ahrs1_bYawRate', 
                'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc', 'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas', 'Dadc1_tas', 
                'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state', 'display_active_screen', 'time']
    
    # Generate arrays where data will be stored
    flightdata = np.empty((48,48321,1))
    units = np.array([])
    descriptions = np.array([]) # Description of each parameter
    headers = np.array([])
    
    
    # Loop through parameters
    for i in range(len(parameters)):
    
        # Extract unit
        # Exception to deal with case where unit does not exist
        try:
            unit = str(flight_data[0][0][parameters[i]]['units'].flat[0].flat[0].flat[0])
        except IndexError:
            unit = 'N/A'
        
        # Extract and reshape data
        data = flight_data[0][0][parameters[i]]['data'].flat[0]
        data = np.resize(data, (48321,1))
        
        # Extract description
        description = str(flight_data[0][0][parameters[i]]['description'].flat[0].flat[0].flat[0])
        
        # Generate header for pandas DataFrame
        header = parameters[i]+' ['+unit+']'
        
        # Append
        flightdata[i] = data
        units = np.append(units,unit)
        descriptions = np.append(descriptions,description)
        headers = np.append(headers,header)
        
    
    # Convert data to pandas dataframe
    flightdata_df = pd.DataFrame(columns=headers)
    
    for i in range(len(parameters)): 
        column_df = pd.DataFrame(data=flightdata[i])
        flightdata_df[headers[i]] = np.resize(flightdata[i],(48321,))
    

    return flightdata_df,headers,descriptions








    