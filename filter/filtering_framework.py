#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 16:27:40 2025

@author: matteocamagna
"""
import technical
import economic
import environmental
import social
import pandas as pd
from utilis import Tg_default, Tair_default

def filtering_framework(stove_geometry, materials_parameters, thresholds, change_wall_thickness=False):
    
    not_technical = []
    not_economic = []
    not_environmental= []
    not_social = []
    suitable_materials = []
    
    thicknesses = (stove_geometry['Thickness inner cylinder'], 
                   stove_geometry['Thickness outer cylinder'])
    
    materials_parameters = pd.DataFrame(materials_parameters)
        
    for _, material in materials_parameters.iterrows():
        
        if change_wall_thickness:
            stove_geometry['Thickness inner cylinder'], stove_geometry['Thickness outer cylinder'] = thicknesses
            stove_geometry['Thickness inner cylinder'], stove_geometry['Thickness outer cylinder'] =\
                wall_thick_to_threshold(stove_geometry, material, thresholds)
            material['Thickness inner cylinder'] = stove_geometry['Thickness inner cylinder']
            material['Thickness outer cylinder'] = stove_geometry['Thickness outer cylinder']
        
        if not technical_filter(stove_geometry, material, thresholds):
            not_technical.append(material)
            continue
        
        if not economic_filter(stove_geometry, material, thresholds):
            not_economic.append(material)
            continue
        
        if not environmental_filter(material, thresholds):
            not_environmental.append(material)
            continue
        
        if not social_filter(material, thresholds):
            not_social.append(material)
            continue
        
        suitable_materials.append(material)
    
    unsuitable_materials = {
        "Didn't satisfy technical requirements" : not_technical,
        "Didn't satisfy economic requirements" : not_economic,
        "Didn't satisfy environmental requirements" : not_environmental,
        "Didn't satisfy social requirements" : not_social
        } 
    
    return suitable_materials, unsuitable_materials

def technical_filter(stove_geometry, material_params, thresholds, Tg=Tg_default, Tair=Tair_default):
    
    temperatures = technical.calculate_temperatures(
        Tg, 
        Tair, 
        stove_geometry['Inner diameter inner cylinder'],
        stove_geometry['Thickness inner cylinder'],
        stove_geometry['Inner diameter outer cylinder'],
        stove_geometry['Thickness outer cylinder'],
        stove_geometry['Height'],
        material_params['Thermal Conductivity (W/m*C)'],
        material_params['Thermal Conductivity (W/m*C)'],
        material_params['Emissivity'],
        material_params['Emissivity'])
    
    temperatures = {
        'Inner wall inner cylinder' : temperatures[0],
        'Outer wall inner cylinder' : temperatures[1],
        'Inner wall outer cylinder' : temperatures[2],
        'Outer wall outer cylinder' : temperatures[3]
        }
    
    degradation_temp = technical.check_degradation_temperature(
        temperatures['Inner wall inner cylinder'],
        material_params['Melting Point (C)'] - thresholds['Temperature Safety Margin']
        )
    
    stress_failure = technical.check_thermal_stress_failure(
        material_params['Yield Strength (Pa) (Elastic Limit)'], 
        material_params['Youngs Modulus (Pa)'], 
        material_params['Thermal Expansion Coefficient (strain/C)'], 
        temperatures['Inner wall inner cylinder'] - temperatures['Outer wall inner cylinder']
        )
    
    outer_wall_temp = technical.check_outer_wall_temperature(
        thresholds['Acceptable outer wall temperature'],
        temperatures['Outer wall outer cylinder']
        )
    
    return all([degradation_temp, stress_failure, outer_wall_temp])

def economic_filter(stove_geometry, material_params, thresholds):
    
    cost = economic.estimate_cost(
        stove_geometry['Height'],
        stove_geometry['Inner diameter inner cylinder'],
        stove_geometry['Thickness inner cylinder'],
        stove_geometry['Inner diameter outer cylinder'],
        stove_geometry['Thickness outer cylinder'], 
        material_params['Density kg/m3'],
        material_params['Density kg/m3'],
        stove_geometry['Thickness top disk'], 
        stove_geometry['Thickness bottom disk'],
        material_params['Density kg/m3'],
        material_params['Density kg/m3'],
        material_params['Price ($/kg)'],
        material_params['Price ($/kg)'],
        material_params['Price ($/kg)'],
        material_params['Price ($/kg)'],
        material_params['Change manufacturing process needed?']
        )
    
    return economic.check_price(cost, thresholds['Cost'])
    
def environmental_filter(material_params, thresholds):
    
    score = environmental.calculate_score(material_params[[
        'Recyclability',
        'Degradability',
        'Toxicity in Use/Disposal',
        'Hazard Classification',
        'Use-Phase Emissions',
        'Repairability'
        ]])
    
    return environmental.check_environmental_sustainability(score[0], thresholds['Environmental score'])

def social_filter(material_params, thresholds):
    
    score = social.calculate_score(material_params[[
        'Labor Ethics',
        'Local Jobs',
        'Cooking Fit',
        'Supply Transparency',
        'Local Availability',
        ]])
    
    return social.check_social_sustainability(score[0], thresholds['Social score'])

def wall_thick_to_threshold(stove_geometry, material_params, thresholds, Tg=Tg_default, Tair=Tair_default):
    
    temperatures = technical.calculate_temperatures(
        Tg, 
        Tair, 
        stove_geometry['Inner diameter inner cylinder'],
        stove_geometry['Thickness inner cylinder'],
        stove_geometry['Inner diameter outer cylinder'],
        stove_geometry['Thickness outer cylinder'],
        stove_geometry['Height'],
        material_params['Thermal Conductivity (W/m*C)'],
        material_params['Thermal Conductivity (W/m*C)'],
        material_params['Emissivity'],
        material_params['Emissivity'])
    
    temperatures = {
        'Inner wall inner cylinder' : temperatures[0],
        'Outer wall inner cylinder' : temperatures[1],
        'Inner wall outer cylinder' : temperatures[2],
        'Outer wall outer cylinder' : temperatures[3]
        }
    
    while not technical.check_outer_wall_temperature(
        thresholds['Acceptable outer wall temperature'],
        temperatures['Outer wall outer cylinder']
        ):
        
        stove_geometry['Thickness inner cylinder'] += 0.002
        stove_geometry['Thickness outer cylinder'] += 0.002
           
        temperatures = technical.calculate_temperatures(
            Tg, 
            Tair, 
            stove_geometry['Inner diameter inner cylinder'],
            stove_geometry['Thickness inner cylinder'],
            stove_geometry['Inner diameter outer cylinder'],
            stove_geometry['Thickness outer cylinder'],
            stove_geometry['Height'],
            material_params['Thermal Conductivity (W/m*C)'],
            material_params['Thermal Conductivity (W/m*C)'],
            material_params['Emissivity'],
            material_params['Emissivity'])
        
        temperatures = {
            'Inner wall inner cylinder' : temperatures[0],
            'Outer wall inner cylinder' : temperatures[1],
            'Inner wall outer cylinder' : temperatures[2],
            'Outer wall outer cylinder' : temperatures[3]
            }
    
    return stove_geometry['Thickness inner cylinder'], stove_geometry['Thickness outer cylinder'] 
  