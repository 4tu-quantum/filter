#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 29 21:37:04 2025

@author: matteocamagna
"""
import heat_transfer_model as htm
import math
from utilis import T0
# This is a test
Tw1i = 800              # Temperature of the inner wall of the inner cylinder [C].
Tw1o = 790              # Temperature of the outer wall of the inner cylinder [C].
Tw2i = 100              # Temperature of the inner wall of the outer cylinder [C].
Tw2o = 90               # Temperature of the outer wall of the outer cylinder [C].
Tg = 1000               # Temperature of combustion gases [C].
Tair = 25               # Temperature of the outer air [C].
D1 = 0.3                # Inner diameter of the inner cylinder [m].
D2 = 0.4                # Inner diameter of the outer cylinder [m].
t1 = 0.01               # thickness of the inner cylinder [m].
t2 = 0.01               # thickness of the outer cylinder [m].
z = 0.4                 # height of the cylinder [m].
k1 = 16                 # conductivity of the inner cylinder [W / (m K)].
k2 = k1                 # conductivity of the outer cylinder [W / (m K)].
epsi1 = 0.7             # emissivity of the inner cylinder [-].
epsi2 = epsi1           # emissivity of the outer cylinder [-].

residuals = htm.heat_balance(Tw1i + T0, Tw1o + T0, Tw2i + T0, Tw2o + T0, Tg + T0, Tair + T0, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2)

print(residuals)

from scipy.optimize import fsolve

hb = lambda Temps : htm.heat_balance(Temps[0], Temps[1], Temps[2], Temps[3], Tg + T0, Tair + T0, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2)

sol, info, ier, mesg = fsolve(hb, (Tw1i+T0, Tw1o+T0, Tw2i+T0, Tw2o+T0), full_output=True)

res = htm.heat_balance(sol[0], sol[1], sol[2], sol[3], Tg, Tair, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2)

print(sol-273.15)

Tw1i = sol[0]
Tw1o = sol[1]
Tw2i = sol[2]
Tw2o = sol[3]

A1i = math.pi * D1 * z
print((htm.q_rad(Tg, Tw1i, D1/2, z) + htm.q_conv_in(Tg, Tair, Tw1i, D1, z)) * A1i / (Tg - Tw1i))

#%% Really basic filter

def technical_filter(material, tresholds, reflect = True):
    
    suitable = False
    
    if reflect == True:
        melt, therm_diff, strength, reflect = material
        tr_melt, tr_therm_diff, tr_strength, tr_reflect = tresholds
        if all([melt >= tr_melt, therm_diff <= tr_therm_diff, strength >= tr_strength, reflect >= tr_reflect]):
            suitable = True
    else:
        melt, therm_diff, strength = material
        tr_melt, tr_therm_diff, tr_strength = tresholds
        if all([melt >= tr_melt, therm_diff <= tr_therm_diff, strength >= tr_strength]):
            suitable = True
    
    return suitable


material_1 = [300, 2, 60, 0.9]

tresholds_1 = [200, 5, 40, 0.85]

suitable_1 = technical_filter(material_1, tresholds_1)

print('material 1', suitable_1)

material_2 = [300, 7, 60]

tresholds_2 = [200, 5, 40]

suitable_2 = technical_filter(material_2, tresholds_2, reflect=False)

print('material 2', suitable_2)

from technical import calculate_temperatures
print(calculate_temperatures(Tg, Tair, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2) - 273.15)

from utilis import calculate_mass

print(calculate_mass(z, D1, t1, D2, t2, 1000, 1000, 0.005, 0.005, 1000, 1000))
