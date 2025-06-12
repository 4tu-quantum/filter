#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 21:09:31 2025

@author: matteocamagna
"""

import pandas as pd
import filtering_framework as ff

material_data = pd.read_excel('Data/Case study materials.xlsx', sheet_name=0, keep_default_na=False)

stove_geometry = {
    'Inner diameter inner cylinder' : 0.30,
    'Thickness inner cylinder' : 0.002,
    'Inner diameter outer cylinder' : 0.35,
    'Thickness outer cylinder' : 0.002,
    'Height' : 0.5,
    'Thickness top disk' : 0.002,
    'Thickness bottom disk' : 0.004
    }

thresholds = {
    'Temperature Safety Margin' : 10,
    'Cost' : 120,
    'Acceptable outer wall temperature' : 250,
    'Environmental score' : 1.5,
    'Social score' : 1.5,
    }

suitable_materials, unsuitable_materials = ff.filtering_framework(stove_geometry, material_data, thresholds)

print('normal')
for m in suitable_materials:
    print(m['Material'])
if len(suitable_materials) == 0:
    print('No suitable material')
    
for k, i in unsuitable_materials.items():
    print(k)
    print([i[j]['Material'] for j in range(len(i))])
    
suitable_materials, unsuitable_materials = ff.filtering_framework(stove_geometry, material_data, thresholds, change_wall_thickness=True)

print('changing wall thickness')
for m in suitable_materials:
    print(m['Material'])
if len(suitable_materials) == 0:
    print('No suitable material')
    
for k, i in unsuitable_materials.items():
    print(k)
    print([i[j]['Material'] for j in range(len(i))])
    
