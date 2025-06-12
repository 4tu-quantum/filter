#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 12:54:12 2025

@author: matteocamagna
"""

import pandas as pd

def check_environmental_sustainability(score, threshold_score):
    """
    Checks whether the material is enviornmentally sustainable

    Parameters
    ----------
    score : ndarray
        Average environmental score attributed to the stove.
    threshold_score : float
        Threshold enviornmental score.

    Returns
    -------
    bool
        Whether the material respects the threshold.

    """
    return score > threshold_score

def calculate_score(parameters):
    """
    

    Parameters
    ----------
    parameters : dict or pandas DataFrame
        Qualitative evaluation of 'Recyclability', 'Degradability', 'Toxicity in Use/Disposal'
        'Hazard Classification', 'Use-Phase Emissions', and 'Repairability' of the material.

    Returns
    -------
    score : list
        scores of the sustainability of each material.

    """
    parameters = pd.DataFrame(parameters)
    _, c = parameters.shape
    score = []
    
    for i in range(c):
        score.append(analyze_sustainability(parameters.iloc[:, i]))
    
    return score

def analyze_sustainability(parameters):
    """
    Converts the qualitative evaluations of each sustainability parameter
    into a score

    Parameters
    ----------
    parameters : dict or pandas DataFrame
        Qualitative evaluation of 'Recyclability', 'Degradability', 'Toxicity in Use/Disposal'
        'Hazard Classification', 'Use-Phase Emissions', and 'Repairability' of the material.

    Returns
    -------
    float
        score assigned to the sustainability of the material.

    """
    score = 0
    par_count = 6
    
    rec = parameters['Recyclability'].lower()
    deg = parameters['Degradability'].lower()
    tox = parameters['Toxicity in Use/Disposal'].lower()
    haz = parameters['Hazard Classification'].lower()
    emis = parameters['Use-Phase Emissions'].lower()
    rep = parameters['Repairability'].lower()
    
    if rec == 'low':
        score += 1
    elif rec == 'moderate':
        score += 2
    elif rec == 'high':
        score += 3
    elif rec == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Recyclability' must be either 'low', 'moderate', 'high', or 'unknown'")
    
    if deg == 'low':
        score += 1
    elif deg == 'moderate':
        score += 2
    elif deg == 'high':
        score += 3
    elif deg == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Degradability' must be either 'low', 'moderate', 'high', or 'unknown'")
    
    if tox == 'high':
        score += 1
    elif tox == 'moderate':
        score += 2
    elif tox == 'low':
        score += 3
    elif tox == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Toxicity in Use/Disposal' must be either 'low', 'moderate', 'high', or 'unknown'")
    
    if haz == 'monitored':
        score += 3
    elif haz == 'restricted':
        score += 1
    elif haz == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Hazard Classification' must be either 'monitored', 'restricted', or 'unknown'")
    
    if emis == 'none':
        score += 3
    elif emis == 'low':
        score += 2
    elif emis == 'moderate':
        score += 1
    elif emis == 'high':
        score += 0
    elif emis == 'unknown':
        par_count -= 1   
    else:
        raise Exception("'Use-Phase Emissions' must be either 'none', 'low', 'moderate', 'high' or 'unknown'")
     
    if rep == 'easy':
        score += 3
    elif rep == 'moderate':
        score += 2
    elif rep == 'difficult':
        score += 1
    elif rep == 'unknown':
        par_count -= 1   
    else:
        raise Exception("'Repairability' must be either 'easy, 'moderate', 'difficult', or 'unknown'")
    
    if par_count == 0:
        raise Exception("At least one parameter is required to assess the sustainability of the material")
    else:
        return score / par_count

