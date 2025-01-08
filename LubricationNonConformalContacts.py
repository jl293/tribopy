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

def gV(a: float, W: float, n0: float, U: float, Rprime: float): # page 134
    '''
    ----------------------------------------------------------------------------------------------------
    The dimensionless viscosity parameter gV

    To determine which type of lubrication is present in a given system, the relative effects of elasticity and viscosity are quantified by two dimensionless parameters.
    '''
    return (a*(W**3))/((n0**2)*(U**2)*(Rprime**4))

def gE(W: float, n0: float, U: float, Rprime: float, Eprime: float): # page 135
    '''
    ----------------------------------------------------------------------------------------------------
    The dimensionless elasticity parameter gE
    '''
    return (W**(8/3))/((n0**2)*(U**2)*(Rprime**(10/3))*(Eprime**(2/3)))

def h0Prime(h0: float, Rprime: float, W: float, U: float, n0: float): # page 136
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

def h0PrimePE(gV: float, gE: float, k: float): # page 136
    '''
    ----------------------------------------------------------------------------------------------------
    Piezoviscous-elastic lubrication (or EHL) occurs when both elasticity and pressure-viscosity effects are significant. This type of lubrication is observed in many common mechanical components, including rolling element bearings and gears. In this regime, the non-dimensional minimum film thickness can be approximated from this function, h0'_PE.

    NOTE: k = 1 for point contacts with circular contact patch and k --> ∞ for line contacts where the contact patch is a rectangle. If you need to calculate for k --> ∞, enter k as 0 in the function

    gE: the dimensionless elasticity parameter
    gV: the dimensionless viscosity parameter
    k: a non-dimensional parameter related to the minimum film thickness.
    '''
    if k == 0:
        k = math.inf
    return 3.42*(gV**0.49)*(gE**0.17)*(1-math.e**(-0.0387*(k**(2/math.pi))))

def h0PrimeIR(k: float): # page 136
    '''
    ----------------------------------------------------------------------------------------------------
    In isoviscous-rigid lubrication, neither viscosity increase nor elastic deformation is significant. In this case, minimum film thickness can be approximated by this function, h0'_IR.

    NOTE: k = 1 for point contacts with circular contact patch and k --> ∞ for line contacts where the contact patch is a rectangle. If you need to calculate for k --> ∞, enter k as 0 in the function

    k: a non-dimensional parameter related to the minimum film thickness.
    '''
    if k == 0:
        k = math.inf
    return 128*k * ((0.131*math.atan(k/2) + 1.683)**2) * (1 + (2/3*k))**-2

def h0PrimePR(gV: float, k: float): # page 136
    '''
    ----------------------------------------------------------------------------------------------------
    In piezoviscous-rigid lubrication, the change of viscosity with pressure affects lubrication behavior. Here, the minimum film thickness is given by the function h0'_PR

    NOTE: k = 1 for point contacts with circular contact patch and k --> ∞ for line contacts where the contact patch is a rectangle. If you need to calculate for k --> ∞, enter k as 0 in the function

    gV: the dimensionless viscosity parameter
    k: a non-dimensional parameter related to the minimum film thickness.
    '''
    if k == 0:
        k = math.inf
    return 1.41*(gV**0.375) * (1 - math.e**(-0.0387*k))

def hoPrimeIE(gE: float, k: float): # page 137
    '''
    ----------------------------------------------------------------------------------------------------
    In isoviscous-elastic lubrication (soft-EHL), observed in contacts between soft bodies, elastically deformation is significant, but the pressure increase with viscosity is negligible. Here, the minimum film thickness is given by the function h0'_IE

    NOTE: k = 1 for point contacts with circular contact patch and k --> ∞ for line contacts where the contact patch is a rectangle. If you need to calculate for k --> ∞, enter k as 0 in the function

    gE: the dimensionless elasticity parameter
    k: a non-dimensional parameter related to the minimum film thickness.
    '''
    if k == 0:
        k = math.inf
    return 8.70*(gE**0.67) * (1 - 0.85*math.e**(-0.31*k**(2/math.pi)))