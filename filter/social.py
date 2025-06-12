#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 15:23:32 2025

@author: matteocamagna
"""
import pandas as pd

def check_social_sustainability(score, threshold_score):
    """
    Checks whether the material is enviornmentally sustainable

    Parameters
    ----------
    score : ndarray
        Average social score attributed to the stove.
    threshold_score : float
        Threshold social score.

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
        Qualitative evaluation of 'Labor Ethics', 'Local Jobs', 'Cooking Fit'
        'Supply Transparency', and 'Local Availability' of the material.

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
        Qualitative evaluation of 'Labor Ethics', 'Local Jobs', 'Cooking Fit'
        'Supply Transparency', and 'Local Availability'.

    Returns
    -------
    float
        score assigned to the sustainability of the material.

    """
    score = 0
    par_count = 5
    
    ethics = parameters['Labor Ethics'].lower()
    jobs = parameters['Local Jobs'].lower()
    fit = parameters['Cooking Fit'].lower()
    transp = parameters['Supply Transparency'].lower()
    avail = parameters['Local Availability'].lower()
    
    if ethics == 'yes':
        score += 3
    elif ethics == 'partial':
        score += 2
    elif ethics == 'no':
        score += 1
    elif ethics == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Labor Ethics' must be either 'yes', 'partial', 'no', or 'unknown'")
        
    if jobs == 'high':
        score += 3
    elif jobs == 'partial':
        score += 2
    elif jobs == 'low':
        score += 1
    elif jobs == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Local Jobs' must be either 'high', 'partial', 'low', or 'unknown'")
    
    if fit == 'yes':
        score += 3
    elif fit == 'partial':
        score += 2
    elif fit == 'no':
        score += 1
    elif fit == 'unknown':
         par_count -= 1
    else:
        raise Exception("'Cooking Fit' must be either 'yes', 'partial', 'no', or 'unknown'")
    
    if transp == 'high':
        score += 3
    elif transp == 'partial':
        score += 2
    elif transp == 'low':
        score += 1
    elif transp == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Supply Transparency' must be either 'high', 'partial', 'low', or 'unknown'")
    
    if avail == 'yes':
        score += 3
    elif avail == 'partial':
        score += 2
    elif avail == 'no':
        score += 1
    elif avail == 'unknown':
        par_count -= 1
    else:
        raise Exception("'Local Availability' must be either 'yes', 'partial', 'no', or 'unknown'")
    
    if par_count == 0:
        raise Exception("At least one parameter is required to assess the sustainability of the material")
    else:
        return score / par_count
    
    