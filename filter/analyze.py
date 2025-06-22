#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 14:34:06 2025

@author: matteocamagna
"""
from utilis import Tg_default, Tair_default
import technical
import economic
import social
import environmental

def evaluate_criteria(material_params, stove_geometry, thresholds):
    """
    Evaluates the material performance for each criteria,

    Parameters
    ----------
    material_params : dict
        contains material paramters.
    stove_geometry : dict
        contains geometrical dimensions of the stove.
    thresholds : dict
        contains the user-defined thresholds.

    Returns
    -------
    criteria : dict
        Dictionary containing the maximum/minimum value of a property
        and the effective value of that property for each criteria.

    """
    criteria = {}
    
    temperatures = technical.calculate_temperatures(
        Tg_default, 
        Tair_default, 
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

    criteria['maximum temperature'] = {
        'maximum service temperature' : material_params['Melting Point (C)'],
        'highest material temperature' : temperatures['Inner wall inner cylinder']
        }
    
    thermal_stress = (material_params['Youngs Modulus (Pa)'] * material_params['Thermal Expansion Coefficient (strain/C)'] *\
         (temperatures['Inner wall inner cylinder'] - temperatures['Outer wall inner cylinder']))          
            
    criteria['stress failure'] = {
        'yield strength' : material_params['Yield Strength (Pa) (Elastic Limit)'],
        'thermal stress' : thermal_stress
        }
    
    criteria['outer temperature'] = {
        'maximum outer wall temperature' : thresholds['Acceptable outer wall temperature'],
        'outer wall temperature' : temperatures['Outer wall outer cylinder']
        }

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
    
    criteria['economic'] = {
        'maximum acceptable cost' : thresholds['Maximum Cost'],
        'cost' : cost}

    environmental_score = environmental.calculate_score(material_params[[
        'Recyclability',
        'Degradability',
        'Toxicity in Use/Disposal',
        'Hazard Classification',
        'Use-Phase Emissions',
        'Repairability'
        ]])[0]
    
    criteria['environmental'] = {
        'minimum accceptable score' : thresholds['Acceptable Environmental score'],
        'score' : environmental_score
        }

    social_score = social.calculate_score(material_params[[
        'Labor Ethics',
        'Local Jobs',
        'Cooking Fit',
        'Supply Transparency',
        'Local Availability',
        ]])[0]
    
    criteria['social'] = {
        'minimum accceptable score' : thresholds['Acceptable Social score'],
        'score' : social_score
        }
    
    return criteria

def normalize_criteria(criteria):
    """
    Normalizes the scores in such a way that the material is suitable if the score
    is larger than 1.

    Parameters
    ----------
    criteria : dictionary of dictionaries
        Dictionary containing the maximum/minimum value of a property
        and the effective value of that property for each criteria.
        The dictionary has the form:
            {
                'maximum temperature' : {
                    'maximum service temperature' : , 'highest material temperature' : ,
                        },
                'stress failure' : {
                    'yield strength' : , 'thermal stress' : ,
                    },
                'outer temperature' : {
                    'maximum outer wall temperature' : , 'outer wall temperature' : ,
                    },
                'economic' : {
                    'maximum acceptable cost' : , 'cost' : ,
                    },
                'environmental' : {
                    'minimum accceptable score' : , 'score' : ,
                    },
                'social' : {
                    'minimum accceptable score' : , 'score' : ,
                    },
                    }

    Returns
    -------
    norm : list
        list containing the normalized scores.

    """    
    norm = []
    
    norm[0] = criteria['maximum temperature']['maximum service temperature'] /\
        criteria['maximum temperature']['highest material temperature']
    
    norm[1] = criteria['stress failure']['yield strength'] /\
        criteria['stress failure']['thermal stress']

    norm[2] = criteria['outer temperature']['maximum outer wall temperature'] /\
        criteria['outer temperature']['outer wall temperature']
        
    norm[3] = criteria['economic']['maximum acceptable cost'] / criteria['economic']['cost']
    
    norm[4] = criteria['environmental']['score'] /\
        criteria['environmental']['minimum accceptable score']
        
    norm[5] = criteria['social']['score'] /\
        criteria['social']['minimum accceptable score']
    
    return norm
