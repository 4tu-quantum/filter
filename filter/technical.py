#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 15:02:27 2025

@author: matteocamagna
"""

import heat_transfer_model as htm
from scipy.optimize import fsolve, root
from utilis import T0

def check_degradation_temperature(T, T_deg):
    """
    Checks whether the melting temperature of the material is higher than the 
    temperature inside the stove

    Parameters
    ----------
    T : float
        Temperature inside the stove.
    T_deg : ndarray
        Degradation temperature of the material.

    Returns
    -------
    bool
        Whether the condition is satisfied.

    """
    return T_deg > T

def check_thermal_stress_failure(sigma_f, E, alpha, DeltaT):
    """
    Checks whether the material is able to withstand thermal stresses

    Parameters
    ----------
    sigma_f : ndarray
        Failure stress of the material (GPa).
    E : ndarray
        Yung's modulus of the material (GPa).
    alpha : ndarray
        Thermal expansion coefficient (1/C).
    DeltaT : 
        Temperature difference across the wall of the stove.

    Returns
    -------
    bool
        Whether the material satisfies the condition.

    """
    return sigma_f > E * alpha * DeltaT

def check_outer_wall_temperature(T_treshold, T_outer_wall):
    """
    Checks whether the calculated heat loss is lower than the threshold value

    Parameters
    ----------
    T_treshold : float
        Treshold temperature of outer wall of the stove (C).
    q : ndarray
        Temperature of the outer wall of the stove (C).

    Returns
    -------
    bool
        Whether the temperature satisfies the condition.

    """
    
    return T_outer_wall < T_treshold

def calculate_temperatures(Tg, Tair, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2):
    """
    Calculates the wall temperatures of the stove
    Parameters
    ----------
    Tg : float
        Temperature of the gas [C].
    Tair : float
        Temperature of the air [C].
    D1 : float
        Inner diameter of the inner cylinder [m].
    t1 : float
        Thickness of the inner cylinder [m].
    D2 : float
        Inner diameter of the outer cylinder [m].
    t2 : float
        Thickness of the outer cylinder [m].
    z : float
        Height of the stove [m].
    k1 : float
        Thermal conductivity of the inner cylinder [W / (m K)].
    k2 : float
        Thermal conductivity of the outer cylinder [W / (m K)].
    epsi1 : float
        emissivity of the inner cylinder [-].
    epsi2 : float
        emissivity of the outer cylinder [-].

    Raises
    ------
    Exception
        Error if the solver doesn't converge.

    Returns
    -------
    numpy.ndarray
        [Inner wall temperature of inner cylinder, Outer wall temperature of inner cylinder,
         Inner wall temperature of outer cylinder, Outer wall temperature of outer cylinder] [K].

    """
    # First guesses
    Tw1i = Tg - 100     # Temperature of the inner wall of the inner cylinder [C].
    Tw1o = Tw1i - 1     # Temperature of the outer wall of the inner cylinder [C].
    Tw2o = Tair + 50    # Temperature of the outer wall of the outer cylinder [C].
    Tw2i = Tw2o + 1     # Temperature of the inner wall of the outer cylinder [C].

    def heat_balance(Temps):
        return htm.heat_balance(Temps[0], Temps[1], Temps[2], Temps[3], Tg + T0, Tair + T0, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2)
    
    #sol, info, ier, mesg = fsolve(heat_balance, (Tw1i + T0, Tw1o + T0, Tw2i + T0, Tw2o + T0), full_output=True, xtol=1e-9)
    sol = root(heat_balance, (Tw1i + T0, Tw1o + T0, Tw2i + T0, Tw2o + T0), method='lm', tol=1e-9)
    
    #if ier == 1:
    if sol.success:
        return sol.x - T0

    else:
       # raise Exception(mesg)
       raise Exception(sol.message)

    