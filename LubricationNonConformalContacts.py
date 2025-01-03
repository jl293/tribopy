# Chapter 9

import math

def SRR(UA: list[float], UB: list[float]): # page 132
    '''
    ----------------------------------------------------------------------------------------------------
    A key parameter affecting viscous friction in non-conformal contacts is the relatibve amount of sliding vs rolling. This is quantified by the slide-to-roll ratio.

    UB: list of linear speeds of the bottom plate
    UA: list of linear speeds of the top plate
    '''
    return (2*abs(sum(UA)-sum(UB)))/(abs(sum(UA))+abs(sum(UB)))

def RatioOfRadii(Rax: float, Rbx: float, Ray: float, Rby: float): # page 133
    '''
    ----------------------------------------------------------------------------------------------------
    The ratio of radii of the two contacting bodies k is used to generalize the equations for contacts between various shapes.

    The limiting values for this parameters are k = 1 for point contacts with circular contact patch and k --> ∞ for line contacts where the contact patch is a rectangle.

    Rax: 1/RAx, the effective raadius of body Ax
    Rbx: 1/RBx, the effective raadius of body Bx
    Ray: 1/RAy, the effective raadius of body Ay
    Rby: 1/RBy, the effective raadius of body By
    '''
    Rx: float
    Ry: float

    Rx = 1/(Rax + Rbx)
    Ry = 1/(Ray + Rby)

    return Ry/Rx

def DimensionlessViscosityParameter(a: float, W: float, n0: float, U: float, Rprime: float): # page 134
    '''
    ----------------------------------------------------------------------------------------------------
    The dimensionless viscosity parameter gV

    To determine which type of lubrication is present in a given system, the relative effects of elasticity and viscosity are quantified by two dimensionless parameters.
    '''
    return (a*(W**3))/((n0**2)*(U**2)*(Rprime**4))

def DimensionlessElasticityParameter(W: float, n0: float, U: float, Rprime: float, Eprime: float): # page 135
    '''
    ----------------------------------------------------------------------------------------------------
    The dimensionless elasticity parameter gE
    '''
    return (W**(8/3))/((n0**2)*(U**2)*(Rprime**(10/3))*(Eprime**(2/3)))

def MinimumFilmThickness(h0: float, Rprime: float, W: float, U: float, n0: float): # page 136
    '''
    ----------------------------------------------------------------------------------------------------
    One of the most important parameters for lubrication is the minimum film thickness h0, since it determines the lubrication regime. Minimum film thickness is non-dimensional as h0' to simplify the empirical expressions.

    h0: minimum film thickness
    Rprime: R' is the effective radius
    W: load
    U: relative speed
    n0: ambient viscosity
    '''
    return((h0)/(Rprime))*((W/(Rprime*U*n0))**2)

def PE_MinimumFilmThicnkess(gV: float, gE: float, k: float): # page 136
    '''
    ----------------------------------------------------------------------------------------------------
    '''
    return 3.42*(gV**0.49)*(gE**0.17)*(1-math.e**(-0.0387*(k**(2/math.pi))))