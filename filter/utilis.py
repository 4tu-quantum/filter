#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 08:54:50 2025

@author: matteocamagna
"""

import pandas as pd
from math import pi

# Physical constants
R = 8.31446261815324 # ideal gas constant [J / (K mol)]
P = 101325 # atmospheric pressure [Pa]
sigma = 5.670*10**-8 # Stefan-Boltzmann constant for radiative heat transfer [W / (m^2 K^4)]
g  = 9.81 # gravitational accelleration [m^2 / s]
T0 = 273.15 # 0 Celsius degrees in Kelvin
Tg_default = 550 # Default value for the temperature of the hot gases inside the stove [Celsius degress]
Tair_default = 25 # Default value for the temperature of the ambient air [Celsius degrees]

# average 2025 exchange rate Euros to United States Dollars
eur_to_usd = 1.08

def read_excel_to_dict(file_path, sheet_name=0):
    """
    Reads an Excel file and returns a dictionary with header cells as keys 
    and the corresponding column values as lists.
    
    Parameters:
        file_path (str): Path to the Excel file.
        sheet_name (str or int, optional): Sheet name or index to read. 
                                           Defaults to 0 (the first sheet).
    
    Returns:
        dict: A dictionary with keys from the header row and corresponding 
              lists of column values.
    """
    # Read the Excel file into a DataFrame
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    
    # Convert the DataFrame to a dictionary with lists as values
    data_dict = df.to_dict(orient='list')
    return data_dict

def calculate_mass(z, D1, t1, D2, t2, density1, density2, t_top, t_bottom, density_top, density_bottom):
    """
    Calculates a rough estimate of the total mass of the stove
    """
    Vc1 = volume_hollow_cylinder(D1, t1, z)
    Vc2 = volume_hollow_cylinder(D2, t2, z)
    V_top = pi * (D2 + 2*t2)**2 / 4 * t_top
    V_bottom = pi * (D2 + 2*t2)**2 / 4 * t_bottom
    mass = Vc1 * density1 + Vc2 * density2 + V_top * density_top + V_bottom * density_bottom
    return mass

def volume_hollow_cylinder(diameter_in, thickness, height):
    """
    Calculates volume of a hollow cylinder given inner diameter, thickness and height
    """
    return pi * ((diameter_in + 2*thickness)**2 - diameter_in**2) / 4 * height
    
def calculate_price(z, D1, t1, D2, t2, density1, density2, t_top, t_bottom, density_top, density_bottom, price1, price2, price_top, price_bottom):
    """
    Calculates total price of the stove
    """
    Vc1 = volume_hollow_cylinder(D1, t1, z)
    Vc2 = volume_hollow_cylinder(D2, t2, z)
    V_top = pi * (D2 + 2*t2)**2 / 4 * t_top
    V_bottom = pi * (D2 + 2*t2)**2 / 4 * t_bottom
    return Vc1 * density1 * price1 + Vc2 * density2 * price2 + V_top * density_top * price_top + V_bottom * density_bottom * price_bottom
    
    