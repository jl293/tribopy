# Chapter 5

import math

def absolute_visc(self, t: float, y: float): # page 68
    '''
    ----------------------------------------------------------------------------------------------------
    Absolute viscosity η, also called dynamic viscosity, is the ratio of shear stress τ to shear strain y.

    t: shear stress
    y: shear strain
    '''
    return t/y

def kinematic_viscosity_v(self, n: float, p: float): # page 68
    '''
    ----------------------------------------------------------------------------------------------------
    Kinematic viscosity v (nu) is another way to quantify a fluid's resistance to deformation. Note that the same variable v (nu) is used to represent two different quantities. For fluids, v is the kinematic viscosity, and for solids, v is the Poisson's Ratio.

    n: absolute viscosity, η
    p: fluid density
    '''
    return n/p

''' 
page 70: log10log10(v + 0.7) = A - Blog10T
'''

def WilliamsLandelFerry(self, ng: float, c1: float, c2: float, T: float, Tg: float): # page 71
    '''
    ----------------------------------------------------------------------------------------------------
    The Williams-Landel-Ferry equation relates viscosity at a given temperature η(T) to the glass-transistion temperature Tg of the fluid. The term "glass transistion" refers to change of the fluid from viscous to hard and relatively brittle or "glassy."

    ng: viscosity, ηg, at the glass transistion temperature Tg
    c1: typically 17.44
    c2: typically 51.6
    T: given temperature
    Tg: glass transistion temperature
    '''
    return ng*(10**((-c1(T-Tg))/(c2+(T-Tg))))

def Vogel(self, Ca: float, Cb: float, Cc: float, T: float): # page 72
    '''
    ----------------------------------------------------------------------------------------------------
    The Vogel equation is a purely empirical model that is highly accurate because it has three independent constants. Note that the temperature units must be consistent.

    Since there are three constants in this equation, it should be fit using more temperature-viscosity data points than the two-parameter models in this class. Due to the accuracy and closed form of this equation, it is commonly used in numerical solutions for lubrication of non-conformal contacts.

    Ca: constant fit to the measured data
    Cb: constant fit to the measured data
    Cc: constant fit to the measured data
    T: given temperature
    '''
    return Ca**(Cb/(T-Cc))

def Barus(self, n0: float, a: float, p: float): # page 72
    '''
    ----------------------------------------------------------------------------------------------------
    At lower pressures, viscosity can be assumed to be constant. However, as pressure increases to the order of tens of MPa and above, viscosity can increase with increasing pressure. This increase is quantified by the pressure-viscosity coefficient a (alpha). There are various definitions of a, but the simplest is based on the Barus equation.

    Absolute or kinematic viscosity can be used in this formulation as long as the viscosity at p and ambient viscosity have the same units. The pressure-viscosity coefficeint a has units of [1/presure] and is typically on the order of 10 [1/GPa] for lubricating oils.

    The Barus equation is known to be inaccurate for many lubricants and conditions. However, its simple functional form conveys the exponential nature of the pressure-viscosity relationship. Also, the large magnitude of a (alpha) correctly indicates that the exponential term will only become significant at high pressures, consistent with physical observations.

    η0: the viscosity at atmospheric pressure
    a: pressure-viscosity coefficient (alpha)
    p: given pressure
    '''
    return n0**(a*p)

def Roelands(self, np: float, n0: float, p: float, pp: float, z: float): # page 72
    '''
    ----------------------------------------------------------------------------------------------------
    Another model that has been shown to be more accurate than the Barus formulation was developed by Roelands. Roelands found reasonably accurate results for reference values of ηp = 6.31 x 10^-5 Pa-s and pp = -0.196 GPa.

    NOTE: This function returns viscosity n (nu) and a (alpha) values as the tuple [n, a].

    np: viscosity at the intersection of the pressure-viscosity isotherms when extrapolated to negative pressures
    n0: viscosity at given/atmospheric pressure
    p: given pressure
    pp: pressure at the intersection of the pressure-viscosity isotherms when extrapolated to negative pressures
    z: pressure-viscosity index
    '''
    n: float
    a: float
    npt: float

    n = np*((n0/np)**(((pp-p)/pp)**z))
    a = (z*math.log10(n0/np))/-pp
    npt = np*(n0/np)**((((pp-p)/pp)**z)*(()))
    
    return n, a

def Roelands_npt(self, np: float, n0: float, p: float, pp: float, z: float, T0: float, Tinf: float, T: float, Sprime: float): # page 73
    '''
    ----------------------------------------------------------------------------------------------------
    NOTE: This function is an offshoot of the Roeland model.
    
    It is often found that the pressure-viscosity index, Z, is independent of temperature up to about 100°C. However, in practice, the effects of pressure and temperature on viscosity are interrelated. This is captured in this model.

    np: viscosity at the intersection of the pressure-viscosity isotherms when extrapolated to negative pressures
    n0: viscosity at given/atmospheric pressure
    p: given pressure
    pp: pressure at the intersection of the pressure-viscosity isotherms when extrapolated to negative pressures
    z: pressure-viscosity index
    T0: reference temperature
    Tinf: the divergence temperature at which viscosity becomes unbounded and Roelands is identified a reference value of Tinf = -135°C
    Sprime: the slope index that quantifies the rate of decrease of viscosity with temperature. The Roelands equation was designed such that S' = 1.0 for all linear paraffins and this parameter often has a value around unity for lubricants
    '''
    return np*(n0/np)**((((pp-p)/pp)**z)*(((T0-Tinf)/(T-Tinf))**Sprime))

def PSSI(self, vf: float, vs: float, vb: float): # page 74
    '''
    ----------------------------------------------------------------------------------------------------
    The shear stability of a lubricant is characterized by the permanent shear stability index (PSSI).

    The vsheared is measured after the lubricant has been subject to rapid, high shear stress using one of several standard tests.

    vf: vfresh is the kinematic viscosity of the lubricant before shearing at 100°C
    vs: vsheared is the kinematic viscosity after shearing at 100°C
    vb: vbase is the viscosity of the base oil of the lubricant at 100°C
    '''
    return (vf-vs)/(vf-vb) * 100

def Carreau(self, ninf: float, n0: float, y: float, ycr: float, m: float, n: float): # page 74
    '''
    ----------------------------------------------------------------------------------------------------
    Carreau is a more robust model for temporary viscosity loss, having a power law relationship between viscosity and shear rate. It is commonly used in engineering.

    The constants ninf, m, and n are found by fitting the Carreau equation to experimental data, although sometimes m = 2 is used such that there are only two fit terms.

    ninf: viscosity at infinite shear rate, constant found by fitting the model to experimental data
    n0: viscosity at zero shear
    y: shear rate
    ycr: critical shear rate, can be approximated as the shear rate at which high-shear viscosity is 10% lower than the Newtonian value.
    m: constant found by fitting the model to experimental data, sometimes m = 2 suffices
    n: constant found by fitting the model to experimental data
    '''
    return ninf + (n0 - ninf) * (1+(y/ycr)**m)**((n-1)/2)

def SpecificGravity(self, t: float, p: float, ph20: float): # page 78
    '''
    ----------------------------------------------------------------------------------------------------
    Fluid density can be reported as specific gravity, which is the density of the fluid of interest divided by the density of water, so it is unitless. In the petroleum industry, specific gravity can be reported in degrees API (American Petroleum Institute).

    NOTE: This function returns the tuple [s, api], which reflects standard specific gravity and American Petroleum Institute specific gravity, respectively.

    t: given temperature, °C
    p: density of given fluid at given temperature
    ph20: density of water at given temperature
    '''
    p: float
    api: float
    s_api: float

    s = p/ph20
    if t != 15.6: # 15.6°C per API standards
        # convert p and ph20 values to those at 15.6°C --> s1/t1 = s2/t
        p = (p/t)*15.6
        ph20 = (ph20/t)*15.6
        s_api = p/ph20
    else:
        s_api = s
    api = (141.5/s_api)-131.5
    return s, api