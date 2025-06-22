#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 11:07:17 2025

@author: matteocamagna
"""
import math
from utilis import sigma, g

def k_air(T):
    """
    Returns air conductivity at temperature T [K]
    """
    return(0.00031847*T**0.7775) # W / (m K)

def dk_dT(T):
    """
    Returns derivative of air conductivity at temperature T [K]
    """
    return (0.000247610425*T**-0.2225)

def nu_air(T):
    """
    Returns kinematic viscosity of air at temperature T [K]
    """
    return((0.0000644*T**2+0.0631*T-9.54)*10**-6) #

def dnu_dT(T):
    """
    Returns derivative of kinematic viscosity of air at temperature T [K]
    """
    return ((0.0001288*T+0.0631)*10**-6)

def Jacobian_hb(Tw1i, Tw1o, Tw2i, Tw2o, Tg, Tair, D1, t1, D2, t2, z, k1, k2, epsi1, epsi2):
    """
    Computes Jacobian of the heat balance to ease and speed up the numerical solution.
    Takes as input the same parameters of the 'heat_balance' function

    Parameters
    ----------
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

    Returns
    -------
    Jac : list of lists
        Jacobian of the 'heat_balance' function.

    """
    A1i = math.pi * D1 * z
    A2o = math.pi * (D2 + 2*t2) * z
    
    dhb1 = dhb2 = dhb3 = dhb4 = [0, 0, 0, 0]
    
    # non-zero derivatives of heat balance 1
    dhb1[0] = (dq_rad_dTw1i(Tg, Tw1i, D1/2, z) + dq_conv_in_dT(Tg, Tair, Tw1i, D1, z)) * A1i \
        - dQ_cond_dTwi(Tw1i, Tw1o, D1, D1+2*t1, k1, z)
    
    dhb1[1] = - dQ_cond_dTwo(Tw1i, Tw1o, D1, D1+2*t1, k1, z)
    
    # non-zero derivatives of heat balance 2
    dhb2[0] = dQ_cond_dTwi(Tw1i, Tw1o, D1, D1+2*t1, k1, z)
    
    dhb2[1] = dQ_cond_dTwo(Tw1i, Tw1o, D1, D1+2*t1, k1, z) - \
        (dQ_conv_annulus_dTw1(Tw1o, Tw2i, D1+2*t1, D2, z) \
         + dQ_rad_annulus_dTw1(Tw1o, Tw2i, D1+2*t1, D2, epsi1, epsi2, z))
    
    dhb2[2] = (dQ_conv_annulus_dTw2(Tw1o, Tw2i, D1+2*t1, D2, z) \
         + dQ_rad_annulus_dTw2(Tw1o, Tw2i, D1+2*t1, D2, epsi1, epsi2, z))
    
    # non zero-derivatives of heat balance 3
    dhb3[1] = dQ_conv_annulus_dTw1(Tw1o, Tw2i, D1+2*t1, D2, z) \
     + dQ_rad_annulus_dTw1(Tw1o, Tw2i, D1+2*t1, D2, epsi1, epsi2, z)
     
    dhb3[2] = (dQ_conv_annulus_dTw2(Tw1o, Tw2i, D1+2*t1, D2, z) \
         + dQ_rad_annulus_dTw2(Tw1o, Tw2i, D1+2*t1, D2, epsi1, epsi2, z)) \
        - dQ_cond_dTwi(Tw2i, Tw2o, D2, D2+2*t2, k2, z)
    
    dhb3[3] = - dQ_cond_dTwo(Tw2i, Tw2o, D2, D2+2*t2, k2, z)
    
    # non zero-derivatives of heat balance 4
    dhb4[2] = dQ_cond_dTwi(Tw2i, Tw2o, D2, D2+2*t2, k2, z)
    
    dhb4[3] = dQ_cond_dTwo(Tw2i, Tw2o, D2, D2+2*t2, k2, z) -\
        (dq_rad_out_dT(Tw2o, Tair, epsi2) + dq_conv_out_dT(Tw2o, Tair, D2+2*t2, z)) * A2o
    
    Jac = [dhb1, dhb2, dhb3, dhb4]
    return Jac

def dq_rad_dTw1i(Tg, Tw1i, r, z):
    """
    Return derivative of the radiative heat flux inside the inner cylinder, 
    with respect to the temperature of the wall of inner cylinder [W/(m^2 K)]

    Parameters
    ----------
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
    return -3.5*sigma*epsig*Tg**0.5*Tw1i**2.5

def dq_conv_in_dT(Tg, Tair, Tw1i, D, z):
    """
    Return derivative of the convective heat flux inside the inner cylinder, 
    with respect to the temperature of the wall of inner cylinder [W/(m^2 K)].

    Parameters
    ----------
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
    constant = 1 / D * 0.023 * (u*D)**0.8 * Pr**0.3
    num = 1/2*dk_dT(T_ave)*(nu_air(T_ave))**0.8 -\
        0.8*(nu_air(T_ave))**(-0.2)*1/2*dnu_dT(T_ave)*k_air(T_ave)
    den = (nu_air(T_ave))**0.64
    d_k_nu_08_dTw1i = num / den
    return constant * (d_k_nu_08_dTw1i - k_air(T_ave) / (nu_air(T_ave))**0.8)

def dQ_cond_dTwi(Twi, Two, Di, Do, k, z):
    """
    Return derivative of heat flux due to conductive heat transfer thorugh
    a cylindrical shell, with respect to the temperature of the inner wall [W/K].
 
    Twi : float
        inner temperature of the wall [K].
    Two : float
        outer temperature of the wall [K].
    k : float
        thermal conductivity of the metal [W/(K m)].
    z : float
        height of the cylinder [m].
    """
    return 2 * math.pi * k * z / math.log(Do / Di)

def dQ_cond_dTwo(Twi, Two, Di, Do, k, z):
    """
    Return derivative of heat flux due to conductive heat transfer thorugh
    a cylindrical shell, with respect to the temperature of the outer wall [W/K].
 
    Twi : float
        inner temperature of the wall [K].
    Two : float
        outer temperature of the wall [K].
    k : float
        thermal conductivity of the metal [W/(K m)].
    z : float
        height of the cylinder [m].
    """
    return - 2 * math.pi * k * z / math.log(Do / Di)

def d_beta_nu2(T_ave):
    """
    Return derivative of beta/nu^2 with respect to the average temperature 'T_ave'
    """
    beta = 1 / T_ave
    dbeta = - T_ave**-2
    return (dbeta * nu_air(T_ave) - 2 * dnu_dT(T_ave) * beta) / nu_air(T_ave)**3
    
def dQ_conv_annulus_dTw1(Tw1, Tw2, D1, D2, z):
    """
    Returns the derivative of the convective heat flux in the annular cavity
    between inner and outer cylinder, with respect to the the temperature of the wall 
    of the inner cylinder [W/K].
    
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
    GrPr = Gr*Pr
    
    if GrPr < 2000:
        dk_dTw1 = 1/2 * dk_dT(T_ave)
    
    elif GrPr >= 2000 and GrPr <= 200000:
        
        dgr_dTw1 = 0.25 * (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**-0.75 *\
            (1/2 * d_beta_nu2(T_ave) * (Tw1 - Tw2) + beta / nu_air(T_ave)**2)
        
        dk_dTw1 = 0.197 * Pr**0.25 * (g * delta**3)**0.25 * (z / delta)**-(1/9) *\
            (dgr_dTw1 * k_air(T_ave) + (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**0.25 * 1/2 * dk_dT(T_ave))
    else:
        
        dgr_dTw1 = 1/3 * (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**(-2/3) *\
            (1/2 * d_beta_nu2(T_ave) * (Tw1 - Tw2) + beta / nu_air(T_ave)**2)
        
        dk_dTw1 = 0.073 * Pr**(1/3) * (g * delta**3)**(1/3) * (z / delta)**-(1/9) *\
            (dgr_dTw1 * k_air(T_ave) + (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**(1/3) * 1/2 * dk_dT(T_ave))
    
    return 2 * math.pi * z / math.log(D2 / D1) * (dk_dTw1 * (Tw1 - Tw2) + k_air(T_ave))

def dQ_conv_annulus_dTw2(Tw1, Tw2, D1, D2, z):
    """
    Returns the derivative of the convective heat flux in the annular cavity
    between inner and outer cylinder, with respect to the the temperature of the wall 
    of the outer cylinder [W/K].
    
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
    GrPr = Gr*Pr
    
    if GrPr < 2000:
        dk_dTw1 = 1/2 * dk_dT(T_ave)
    
    elif GrPr >= 2000 and GrPr <= 200000:
        
        dgr_dTw1 = 0.25 * (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**-0.75 *\
            (1/2 * d_beta_nu2(T_ave) * (Tw1 - Tw2) - beta / nu_air(T_ave)**2)
        
        dk_dTw1 = 0.197 * Pr**0.25 * (g * delta**3)**0.25 * (z / delta)**-(1/9) *\
            (dgr_dTw1 * k_air(T_ave) + (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**0.25 * 1/2 * dk_dT(T_ave))
    else:
        
        dgr_dTw1 = 1/3 * (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**(-2/3) *\
            (1/2 * d_beta_nu2(T_ave) * (Tw1 - Tw2) - beta / nu_air(T_ave)**2)
        
        dk_dTw1 = 0.073 * Pr**(1/3) * (g * delta**3)**(1/3) * (z / delta)**-(1/9) *\
            (dgr_dTw1 * k_air(T_ave) + (beta * (Tw1 - Tw2) / nu_air(T_ave)**2)**(1/3) * 1/2 * dk_dT(T_ave))
    
    return 2 * math.pi * z / math.log(D2 / D1) * (dk_dTw1 * (Tw1 - Tw2) - k_air(T_ave))

def dQ_rad_annulus_dTw1(Tw1, Tw2, D1, D2, epsi1, epsi2, z):
    """
    Returns the derivative of the radiative heat flux in the annular cavity
    between inner and outer cylinder, with respect to the the temperature of the wall 
    of the inner cylinder [W/K].
    
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
    return sigma * math.pi * D1 * z * 4 * Tw1**3 / (1 / epsi1 + (1 / epsi2 - 1) * (D1 / D2))

def dQ_rad_annulus_dTw2(Tw1, Tw2, D1, D2, epsi1, epsi2, z):
    """
    Returns the derivative of the radiative heat flux in the annular cavity
    between inner and outer cylinder, with respect to the the temperature of the wall 
    of the outer cylinder [W/K].
    
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
    return - sigma * math.pi * D1 * z * 4 * Tw2**3 / (1 / epsi1 + (1 / epsi2 - 1) * (D1 / D2))

def dq_rad_out_dT(Two, Ta, epsi0):
    """
    Return derivative of the heat flux due to radiative heat transfer from outer
    cylinder to the air outside, with respect to the temperature of the wall of
    the outer cylinder [W/(m^2 K)].
    
    Two : float
        temperature of the outer wall of the cylinder [K].
    Ta : float
        air temperature [K].
    epsi0 :
        emissivity of the cylinder [K].
    """
    return 4*sigma*epsi0*Two**3

def dq_conv_out_dT(Tw, Ta, D, z):
    """
    Return derivative of the heat flux due to convective heat transfer from outer
    cylinder to the air outside, with respect to the temperature of the wall of
    the outer cylinder [W/(m^2 K)].
    
    Tw : float
        temperature of the wall of the cylinder [K].
    Ta : float
        outside air temperaure [K].
    D : float
        diameter of the outer cylinder [m].
    z : float
        height of the cylinder [m].
    """
    T_ave = (Tw + Ta) / 2 # average air temperature at which evaluate fluid properties
    
    beta = 1 / T_ave #  coefficient of volume expansion
    Pr = 0.685 # Prandtl number
    dgr = beta / nu_air(T_ave**2) + (Tw - Ta) * 1/2 * d_beta_nu2(T_ave)
    d_kgr_dT = 1/2 * dk_dT(T_ave) * (beta * (Tw - Ta) / nu_air(T_ave)**2)**0.25 +\
        k_air(T_ave) * 0.25 * (beta * (Tw - Ta) / nu_air(T_ave)**2)**-0.75 * dgr
    
    return 0.59 * (Pr * g * z**3)**0.25 / D * (d_kgr_dT * (Tw - Ta) +\
                                               k_air(T_ave) * (beta * (Tw - Ta) / nu_air(T_ave)**2)**0.25)
