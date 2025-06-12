#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 12:03:04 2025

@author: matteocamagna
"""
from utilis import calculate_price, eur_to_usd

def check_price(estimated_cost, threshold_cost):
    """
    Checks whether the estimated price is higher than the 
    threshold price

    Parameters
    ----------
    estimated_price : ndarray
        Estimated price of the stove [2025 USD].
    threshold_price : float
        Threshold price of the stove [2025 USD].

    Returns
    -------
    bool
        Whether the material respects the threshold.

    """
    return estimated_cost < threshold_cost

def estimate_cost(z, D1, t1, D2, t2, density1, density2, t_top, t_bottom, density_top,
                  density_bottom, price1, price2, price_top, price_bottom, change_manufacturing_process):
    """
    Estimates the cost of production of the stove
    
    Parameters
    ----------
    z : float
        Height of the stove [m].
    D1 : float
        inner diamter of the inner cylinder [m].
    t1 : float
        thickness of the inner cylinder [m].
    D2 : float
        inner diamter of the outer cylinder [m].
    t2 : float
        thickness of the outer cylinder [m].
    density1 : float
        density of the material of the inner cylinder [kg/m^3].
    density2 : float
        density of the material of the outer cylinder [kg/m^3].
    t_top : float
        thickness of the top of the stove [m].
    t_bottom : float
        thickness of the bottom of the stove [m].
    density_top : float
        density of the top of the stove [kg/m^3].
    density_bottom : float
        density of the top of the stove [kg/m^3].
    price1 : float
        price of the material of the inner cylinder [USD/kg].
    price2 : float
        price of the material of the outer cylinder [USD/kg].
    price_top : float
        price of the material of the top [USD/kg].
    price_bottom : float
        price of the material of the bottom [USD/kg].
    change_manufacturing_process : string
        whether the material considered requires a different manufacturing process
        can be 'yes' or 'no'.

    Raises
    ------
    Exception
        Error if 'change_manufacturing_process' is not 'yes' or 'no'.

    Returns
    -------
    float
        estimated cost of the stove [2025 USD].

    """
    # Cost of the material in Euros
    material_cost = calculate_price(z, D1, t1, D2, t2, density1, density2, t_top,
                                    t_bottom, density_top, density_bottom, price1,
                                    price2, price_top, price_bottom)
    
    ## convert the price to USD
    #material_cost *= eur_to_usd
    
    # assume 5% of the material is lost during the manufacturing process
    material_cost =  material_cost * 1.05
    
    if change_manufacturing_process.lower() == 'yes':
        return material_cost + 40
    elif change_manufacturing_process.lower() == 'no':
        return material_cost + 100
    else:
        raise Exception("'Change manufacturing process' can be either 'yes' or 'no'")
        
        