# Chapter 8

import sympy
import math

# def Reynolds1D(h: float, p: float, x: float, UA: float, UB: float, n: float): # page 116
#     '''
#     ----------------------------------------------------------------------------------------------------
#     The two fundamental equations of fluid dynamics that describe fkuid flow through a lubricated gap are the Navier-Stokes and continuity equations. Combining these equations for lubricant flow and applying a series of simplifying expressions yields a reduced form fo what is called the Reynolds equation, of which this is the 1D form.

#     h: fluid film thickness
#     p: pressure
#     x: position in the direction of motion
#     UA: relative speed of bottom plate
#     UB: relative speed of top plate
#     n: nu is viscosity
#     '''
#     U: float
#     U = abs(UA-UB)

    

#     return

# def LoadCarryingCapacity(): # page 117
    # return

# def ViscousFrictionForce(): # page 117
#     return

def FilmThicknessInclinedPlane(h0: float, sh: float, x: float, l: float): # page 118
    '''
    ----------------------------------------------------------------------------------------------------
    The inclined plane geometry is defined by the shoulder height sh, minimum film thickness h0, and length l and is calculated as such:
    h(x) = h0 + sh(1 - x/l)
    where 

    h0: film thickness
    sh: shoulder height
    x: 0 < x < 1
    l: length
    '''
    return h0 + sh*(1 - (x/l))

def ReynoldsIntegrated(n: float, U: float, l: float, sh: float, x: float, h0: float): # page 118
    '''
    ----------------------------------------------------------------------------------------------------
    For a fluid with viscosity n and relative speed U, the Reynolds equation can be integrated analytically for inclinded plane geometry to obtain the pressure distribution along the bearing.

    p(x) = ((n*U*l)/(sh**2)) * ((6*(x/l)*(1-(x/l)))/((((h0/sh)+1-(x/l))**2) * (1+2*(h0/sh))))

    n: viscosity
    U: relative speed
    l: distance
    sh: shoulder height
    x: 0 < x < 1
    h0: film thickness
    '''
    return ((n*U*l)/(sh**2)) * ((6*(x/l)*(1-(x/l)))/((((h0/sh)+1-(x/l))**2) * (1+2*(h0/sh))))

def pmax(n: float, U: float, l: float, sh: float, h0: float): # page 118
    '''
    ----------------------------------------------------------------------------------------------------
    The maximum pressure occurs where the pressure gradient is zero, i.e., dp/dx = 0, and can calculated as:
    pmax = (3*n*U*l*sh)/(2*h0*(sh+h0)*(sh+2*h0))

    n: viscosity
    U: relative speed
    l: distance
    sh: shoulder height
    h0: film thickness
    '''
    return (3*n*U*l*sh)/(2*h0*(sh+h0)*(sh+2*h0))

def WPrimeInclinedPlane(n: float, U: float, l: float, sh: float, h0: float): # page 119
    '''
    ----------------------------------------------------------------------------------------------------
    Calculate the load-carrying capacity for an inclined plane.

    n: viscosity
    U: relative speed
    l: distance
    sh: shoulder height
    h0: film thickness
    '''
    return ((n*U*(l**2))/(sh**2)) * (6*math.log10((h0+sh)/h0) - ((12*sh)/(sh - 2*h0)))

def FvPrimeInclinedPlane(n: float, U: float, l: float, sh: float, h0: float): # page 119
    '''
    ----------------------------------------------------------------------------------------------------
    Calculate the viscous friction force acting on an inclined plane.

    n: viscosity
    U: relative speed
    l: distance
    sh: shoulder height
    h0: film thickness
    '''
    return -((n*U*l)/sh) * (4*math.log10(h0/(h0+sh)) + ((6*sh)/(sh - 2*h0)))

def fv(FvPrime: float, WPrime: float): # page 119
    '''
    ----------------------------------------------------------------------------------------------------
    The ratio of the friction force per width Fv' to the load per width W' can be used to calculate the viscous friction coefficient.

    FvPrime: Fv' or the viscous friction force acting on an inclined plane
    WPrime: W' or the load-carrying capacity for an inclined plane
    '''
    return FvPrime/WPrime

def SommerfieldNumber(R: float, c: float, n: float, w: float, pa: float): # page 122
    '''
    ----------------------------------------------------------------------------------------------------
    There is no analytical solution to the Reynolds equation for the jounal bearing geometry. However, solutions can be obtained numerically using the Sommerfield equation. This solution enables calculation of various engineering parameters as functions of a characteristic bearing number, also called the Sommerfield number.

    S = (R/c)^2 * (n*w)/(pa)

    R: nominal bearing radius
    c: clearance
    n: viscosity
    w: angular speed
    pa: applied pressure
    '''
    return ((R/c)**2) * (n*w)/(pa)

def pa(W: float, R: float, l: float): # page 122
    '''
    ----------------------------------------------------------------------------------------------------
    The units of various parameters for the Sommerfield solution must be chosen such that S is dimensionless, and the angular speed w has the unit rev/s. Since contact area is poorly defined for conformal contacts, applied pressure is calculated as the load divided by the projected area:
    pa = W/2Rl

    W: load
    R: nominal bearing radius
    l: length
    '''
    return W/(2*R*l)