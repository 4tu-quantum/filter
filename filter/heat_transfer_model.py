#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 29 10:22:49 2025

@author: matteocamagna
"""
import math
from utilis import sigma, g

########## Model developed by DA Quan ###############

#air properties
def mu_air(T):
    """
    Returns dynamic viscosity of air at temperature T [K]
    """
    return(0.0447*10**-5*T**0.7775) # kg / (m s)

def k_air(T):
    """
    Returns air conductivity at temperature T [K]
    """
    return(0.00031847*T**0.7775) # W / (m K)

def Cp_air(T):
    """
    Returns specific heat capacity at constant pressure of air at temperature T [K]
    """
    return(0.9362 +0.0002*T) #J / (g K)

def nu_air(T):
    """
    Returns kinematic viscosity of air at temperature T [K]
    """
    return((0.0000644*T**2+0.0631*T-9.54)*10**-6) # 

def heat_balance(Tw1i, Tw1o, Tw2i, Tw2o, Tg, Tair, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2):
    """
    When the heat balance is closed the four values returned by this code should be
    equal to zero
    
    Tw1i : float
        Temperature of the inner wall of the inner cylinder [K].
    Tw1o : float
        Temperature of the outer wall of the inner cylinder [K].
    Tw2i : float
        Temperature of the inner wall of the outer cylinder [K].
    Tw2o : float
        temperature of the outer wall of the outer cylinder [K].
    Tg : float
        Temperature of combustion gases [K].
    Tair : float
        Temperature of the outer air [K].
    D1 : float
        Inner diameter of the inner cylinder [m].
    D2 : float
        Inner diameter of the outer cylinder [m].
    t1 : float
        thickness of the inner cylinder [m].
    t2 : float
        thickness of the outer cylinder [m].
    z : float
        height of the cylinder [m].
    k1 : float
        conductivity of the inner cylinder [W / (m K)].
    k2 : float
        conductivity of the outer cylinder [W / (m K)].
    epsi1 : float
        emissivity of the inner cylinder [-].
    epsi2 : float
        emissivity of the outer cylinder [-].
    """
    A1i = math.pi * D1 * z
    A2o = math.pi * (D2 + 2*t2) * z
    
    # heat balance at the interface between inner-cylinder-air and inner-cylinder-wall
    hb1 = (q_rad(Tg, Tw1i, D1/2, z) + q_conv_in(Tg, Tair, Tw1i, D1, z)) * A1i \
        - Q_cond(Tw1i, Tw1o, D1, D1+2*t1, k1, z)
        
    # heat balance at the interface between outer-cylinder-wall and annulus-air
    hb2 = Q_cond(Tw1i, Tw1o, D1, D1+2*t1, k1, z) - \
        (Q_conv_annulus(Tw1o, Tw2i, D1+2*t1, D2, z) \
         + Q_rad_annulus(Tw1o, Tw2i, D1+2*t1, D2, epsi1, epsi2, z))
            
    # heat balance at the interface between annulus-air and outer-cylinder-wall
    hb3 = (Q_conv_annulus(Tw1o, Tw2i, D1+2*t1, D2, z) +\
           Q_rad_annulus(Tw1o, Tw2i, D1+2*t1, D2, epsi1, epsi2, z)) -\
        Q_cond(Tw2i, Tw2o, D2, D2+2*t2, k2, z)
        
    # heat balance between outer cylinder wall and outer air
    hb4 = Q_cond(Tw2i, Tw2o, D2, D2+2*t2, k2, z) -\
        (q_rad_out(Tw2o, Tair, epsi2) + q_conv_out(Tw2o, Tair, D2+2*t2, z)) * A2o
    
    return hb1, hb2, hb3, hb4

def q_rad(Tg, Tw1i, r, z):
    """
    Returns heat flux per unit area [W/m^2] due to radiative heat transfer
    
    Tg : float
        gas temperature [K]
    Tw1 : float
        inner wall temperature [K]
    r : float
        radius of the inner cylinder [m].
    z : float
        height of the inner cylinder [m].
    """
    
    vol_medium = math.pi*r*r*z
    A_medium = 2*math.pi*r**2 + 2*math.pi*r*z
    A = (0.848 + (9.02*10**-4)*Tg)
    B = (0.9589 + (4.8*10**-6)*Tg)
    
    epsig = math.e**(A + B*math.log(0.2*3.6*(vol_medium/A_medium)))
    alphag = epsig*(Tg/Tw1i)**0.5
    
    return sigma*(epsig*Tg**4 - alphag*Tw1i**4)

def q_conv_in(Tg, Tair, Tw1i, D, z):
    """
    Return heat flux per unit area [W/m^2] due to convective heat transfer
    
    Tg : float
        gas temperature [K].
    Tair : TYPE
        outer air temperature [K].
    Tw1 : float
        inner wall temperature [K].
    D : float
        diameter of the inner cylinder [m].
    z : float
        height of the inner cylinder [m].
    """
    Pr = 0.685 # Prandtl number
    T_ave = (Tg + Tw1i) / 2 # average air temperature at which evaluate fluid properties
    u = (2*g*z*(Tg/Tair-1))**(1/2) #flow velocity
    Re = u*D/nu_air(Tg) # Reynolds number
    Nu = 0.023*(Re**0.8)*(Pr**0.3) # Nusselt number
    return (k_air(T_ave)/D) * Nu * (Tg - Tw1i)

def Q_cond(Twi, Two, Di, Do, k, z):
    """
    Return heat flux [W] due to conductive heat transfer thorugh a cylindrical shell
 
    Twi : float
        inner temperature of the wall [K].
    Two : float
        outer temperature of the wall [K].
    k : float
        thermal conductivity of the metal [W/(K m)].
    z : float
        height of the cylinder [m].
    """
    return 2 * math.pi * k * z * (Twi - Two) / math.log(Do / Di)

def Q_conv_annulus(Tw1, Tw2, D1, D2, z):
    """
    This function was not included in the original code, it has been written to 
    model the heat transfer in the annualr cavity between the inner and the outer
    cylinder. The heat transfer relation herein used come from 
    JP Holman, Heat Transfer, 10th edition
    
    Tw1 : float
        temperature of the innermost wall [K].
    Tw2 : float
        temperature of the outermost wall [K].
    D1 : float
        inner diameter of the inner cylinder [m].
    D2 : float
        inner diameter of the outer cylinder [m].
    z : float
        height of the cylinder [m].
    """
    delta = D2 - D1
    T_ave = (Tw1 + Tw2) / 2 # average air temperature at which evaluate fluid properties
    
    beta = 1 / T_ave #  coefficient of volume expansion
    Pr = 0.685 # Prandtl number
    Gr = g * beta * (Tw1 - Tw2) * delta**3 / nu_air(T_ave)**2
    k = k_air(T_ave)
    GrPr = Gr*Pr
    if GrPr < 2000:
        k_e = k
    elif GrPr > 2000 and GrPr < 200000: 
        k_e = k * 0.197 * (GrPr)**0.25 * (z / delta)**-(1/9)
    else:
        k_e = k * 0.073 * GrPr**(1/3) * (z / delta)**-(1/9)
    return 2 * math.pi * k_e * z * (Tw1 - Tw2) / math.log(D2 / D1)

def Q_rad_annulus(Tw1, Tw2, D1, D2, epsi1, epsi2, z):
    """
    This function was not included in the original code, it has been written to 
    model the heat transfer in the annualr cavity between the inner and the outer
    cylinder. The heat transfer relation herein used come from 
    JP Holman, Heat Transfer, 10th edition
    
    Return radiative heat flux [W/m^2] from inner to outer cylinder
    Tw1 : float
        temperature of the innermost wall [K].
    Tw2 : float
        temperature of the outermost wall [K].
    D1 : float
        inner diameter of the inner cylinder [m].
    D2 : float
        inner diameter of the outer cylinder [m].
    epsi1 : float
        emissivity of the inner wall
    epsi2 : float
        emissivity of the outer wall  
    z : float
        height of the cylinder [m].
    """
    return sigma * math.pi * D1 * z * (Tw1**4 - Tw2**4) / (1 / epsi1 + (1 / epsi2 - 1) * (D1 / D2))

def q_rad_out(Two, Ta, epsi0):
    """
    Return heat flux [W/m^2] due to radiative heat transfer from inner (outer) 
    cylinder and the air in the annular cavity (air outside)
    Two : float
        temperature of the outer wall of the cylinder
    Ta : float
        air temperature
    """
    return sigma*epsi0*(Two**4 - Ta**4)
    
def q_conv_out(Tw, Ta, D, z):
    """
    Return heat flux [W/m^2] due to convective heat transfer from inner (outer)
    cylinder and the air in the annular cavity (air outside)
    """
    T_ave = (Tw + Ta) / 2 # average air temperature at which evaluate fluid properties
    
    beta = 1 / T_ave #  coefficient of volume expansion
    Pr = 0.685 # Prandtl number
    Gr = g*beta*(Tw - Ta)*z**3/(nu_air(T_ave)**2) #Grashoff number
    Nu = 0.59*(Gr*Pr)**0.25
    return k_air(T_ave) / D * Nu * (Ta - Tw)
