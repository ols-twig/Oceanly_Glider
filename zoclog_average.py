# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 13:03:10 2025

@author: OllieTwigge
"""
import re
import pandas as pd

def zoclog_average(input_file, phrase):
    input_file = r"C:\temp\GliderLog-250305_01.log"

    clog = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if phrase in line:
                clog.append(line[-8:-1])
                
                

    clog = [float(''.join(re.findall(r'-?\d*\.?\d+', item))) for item in clog]
    clog = pd.Series(clog)
    
    print(phrase)
    print(clog.describe())
#%%
input_file = r"C:\temp\GliderLog-250305_01.log"
zoclog_average(input_file, 'science wrote:sci_rbrctd_temperature_00(degC)')
zoclog_average(input_file, 'science wrote:sci_rbrctd_salinity_00(psu)')

