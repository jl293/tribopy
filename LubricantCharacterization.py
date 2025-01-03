# Chapter 7

def ViscosityIndex(L: float, U: float, H: float): # page 107
    '''
    ----------------------------------------------------------------------------------------------------
    The rate of decrease of lubricant viscosity with increasing temperature is critical to its performance as a lubricant. This behavior is characterized by the viscosity index or VO. A good lubricant has a high VI which means that its viscosity decreases slowly with temperature.

    VI is calculated based upon viscosities measured at 40°C and 100°C. Reference the ASTM standard for VI to determine your lubricant's L and H values in units of cSt.
    '''
    return (L-U)/(L-H) * 100
