# Chapter 1

def friction(F: float, W: float) -> float: # page 13
    ''' 
    ----------------------------------------------------------------------------------------------------
    Calculate the friction coefficient, or the force between two bodies resisting their relative motion.

    F: friction force
    W: normal force
    '''
    return F/W

def lambda_ratio(h0: float, Ra: float) -> float: # page 15
    ''' 
    ----------------------------------------------------------------------------------------------------
    Quantify the relative magnitudes of the film thickness to surface roughness.

    h0: minimum film thickness
    Ra: average roughness of the contacting surfaces
    '''
    return h0/Ra

def hersey_number(w: float, n: float, pa: float) -> float: # page 17
    ''' 
    ----------------------------------------------------------------------------------------------------
    A dimensionless parameter used to quantify the Stribeck curve as friction coefficient vs speed times viscosity divided by pressure.

    w: angular speed
    n: viscosity
    pa: applied pressure
    '''
    return (w*n)/pa