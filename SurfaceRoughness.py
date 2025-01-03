# Chapter 3

import math

def gaussian_dist(sig_s: float, z: float) -> float: # page 40
    '''
    ----------------------------------------------------------------------------------------------------
    Probability distribution function for surface roughness quantification, approximated by this Gaussian distribution.

    sig_s: standard deviation of the distribution
    z: height
    '''
    return (1/(sig_s*math.sqrt(2*math.pi)))**(-(z**2)/(2*sig_s))

def zi_prime_mean(zi_prime: list): # page 41 - define zi_prime as a list of float values
    '''
    ----------------------------------------------------------------------------------------------------
    zi',mean is the reference height used to calculate zi and is equal to z'mean, which should be preferred over this function if available.

    zi': list of measured height values *MUST BE ENTERED AS LIST*
    '''
    
    if len(zi_prime) == 0:
        print("Please define zi_prime as a list (ex. [4.1, 3.77, ..., n])")
        return 0
    N = len(zi_prime)
    return (1/N)*(sum(zi_prime)/N)

def zi(zi_prime: list, zi_prime_mean: float): # page 41
    '''
    ----------------------------------------------------------------------------------------------------
    Approximates surface heights relative to an arbitrary reference zi'.

    zi': list of measured height values *MUST BE ENTERED AS LIST*
    zi_prime_mean: reference height
    '''
    return zi_prime - zi_prime_mean

def Ra(zi: list) -> float: # page 42
    '''
    ----------------------------------------------------------------------------------------------------
    Average roughness of a surface.

    zi: list of measured surface heights
    '''
    N = len(zi)
    return (1/N)*((sum(abs(zi))/N))

def Rq(zi: list) -> float: # page 42
    '''
    ----------------------------------------------------------------------------------------------------
    Root-mean-square roughness aka RMS of a surface.

    zi: list of measured surface heights
    '''
    N = len(zi)
    return math.sqrt((1/N)*((sum(zi)**2)/N))

def sig_s(Rq: float): # page 42
    '''
    ----------------------------------------------------------------------------------------------------
    For a Gaussian distribution of surface heights, the RMS roughness is related to the deviation of the probability distribution function sigma_s.

    Rq: RMS roughness
    '''
    return 1.25*Rq

def Ra_composite(Ra1: float, Ra2: float) -> float: # page 43
    '''
    ----------------------------------------------------------------------------------------------------
    The composite average roughness of two surfaces with average roughness Ra1 and Ra2.

    Ra1: surface 1 average roughness
    Ra2: surface 2 average roughness
    '''
    return Ra1 + Ra2

def Rq_composite(Rq1: float, Rq2: float) -> float: # page 43
    '''
    ----------------------------------------------------------------------------------------------------
    The composite root-mean-square roughness of two surfaces with RMS of Rq1 and Rq2.

    Rq1: surface 1 RMS roughness
    Rq2: surface 2 RMS roughness
    '''
    return math.sqrt(Rq1**2 + Rq2**2)

def Rsk(Rq: float, zi: list) -> float:  # page 43
    '''
    ----------------------------------------------------------------------------------------------------
    Skewness captures the relative contribution of peaks and valleys to the average roughness. Specifically, a surface with positive skewness will have narrow, tall peaks and wide, shallow valleys; a negative skewness will have wide, short peaks and narrow, deep valleys.

    Rq: root-mean-square roughness of a given surface
    zi: list of measured heights
    '''
    N = len(zi)
    return (1/(Rq**3))*(1/N)*((sum(zi)**3)/N)

def Rku(Rq: float, zi: list) -> float: # page 44
    '''
    ----------------------------------------------------------------------------------------------------
    Kurtosis reflects the sharpness of the peaks and valleys of a given surface, where the reference value is three. A kurtosis of smaller than three indicates broad peaks/valleys and larger than three means sharp peaks/valleys.

    Rq: root-mean-square roughness of a given surface
    zi: list of measured heights
    '''
    N = len(zi)
    return (1/(Rq**4))*(1/N)*((sum(zi)**4)/N)

def autocorrelation(Rq: float, l: float, zi: list, zil): # page 45
    '''
    ----------------------------------------------------------------------------------------------------
    The lateral characteristics of a surface can be quantified using the autocorrelation function of surface heights. The functions start at x_l = 0 amd decrease toward 0 with increasing x_l.

    Rq: root-mean-square roughness of a given surface
    l: lag
    zi: list of measured values
    zil: list of measured values + lag
    '''
    N = len(zi)
    return (1/((Rq**2)*(N-l)))*(sum(zi*zil)/(N-l))

def exp_autocorrelation(xl: list, B_star: float, B: float): # page 46
    '''
    ----------------------------------------------------------------------------------------------------
    Many engineering surfaces have been found to have an exponential autocorrelation function C(xl) = exp(-xl/B). B* = 2.3B, so B can be calculated by a give B* as B = B*/2.3.

    If B is given, enter B_star as 0. If B* is given and not B, enter B* as given and B as 0.

    xl: distance between surface height data points, also called the lag
    B*: correlation length; the distance at which the autocorrelation function decreases to some small value, typically C(xl) = 0.1.
    B: decay constant
    '''
    if B_star == 0 and B != 0:
        return 1**(-xl/B)
    if B_star != 0 and B == 0:
        B = B_star/2.3
        return 1**(-xl/B)
    else:
        return print('Please enter either B or B*, not both. If B is given, enter B_star as 0. If B* is given and not B, enter B* as given and B as 0')

def orientation(Bx_prime: float, Bx: float, xl_B: float, By_prime: float, By: float, xl_Y: float): # page 46
    '''
    ----------------------------------------------------------------------------------------------------
    An important lateral surface property is orientation. Alignment or ordering of surface features in one direction is characteristic of certain types of machining or finishing techniques, e.g. honing, used for engineering components.

    The degree of alignment or order of a surface can be quantified as an orientation parameter, where Bx* and By* are the correlation lengths of the surface in two orthogonal directions, and Bx* ≥ By* by convention. B* = 2.3B, so B can be calculated by a give B* as B = B*/2.3.

    If B is given, enter B_prime as 0. If B_prime is given and not B, enter B_prime as given and B as 0.

    Bx*: larger correlation length
    Bx: larger decay constant
    xl_X: X surface's distance between surface height data points, also called the lag
    By*: smaller correlation length
    By: smaller decay constant
    xl_Y: Y surface's distance between surface height data points, also called the lag
    '''
    if Bx_prime == 0 and Bx != 0: # Handle Bx*
        Bx_prime = 1**(-xl_B/Bx)
    else:
        print("Error: Please enter either Bx_prime or Bx as 0, not both.")
    
    if By_prime == 0 and By != 0: # Handle By*
        By_prime = 1**(-xl_Y/By)
    else:
        print("Error: Please enter either By_prime or By as 0, not both.")

    return Bx_prime/By_prime