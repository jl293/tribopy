''' 
Contains material from:
Chapter 1: Fundamental Concepts
Chapter 2: Solid Materials
Chapter 3: Surface Roughness
Chapter 4: Non-Conformal Contact
Chapter 5: Liquid Properties
'''

import math

####################################################################################

# Chapter 1

class FundamentalConcepts:
    def friction(self, F: float, W: float) -> float: # page 13
        ''' 
        ----------------------------------------------------------------------------------------------------
        Calculate the friction coefficient, or the force between two bodies resisting their relative motion.

        F: friction force
        W: normal force
        '''
        return F/W

    def lambda_ratio(self, h0: float, Ra: float) -> float: # page 15
        ''' 
        ----------------------------------------------------------------------------------------------------
        Quantify the relative magnitudes of the film thickness to surface roughness.

        h0: minimum film thickness
        Ra: average roughness of the contacting surfaces
        '''
        return h0/Ra

    def hersey_number(self, w: float, n: float, pa: float) -> float: # page 17
        ''' 
        ----------------------------------------------------------------------------------------------------
        A dimensionless parameter used to quantify the Stribeck curve as friction coefficient vs speed times viscosity divided by pressure.

        w: angular speed
        n: viscosity
        pa: applied pressure
        '''
        return (w*n)/pa

####################################################################################

# Chapter 2 

class SolidMaterials:
    def shear_modulus(self, E: float, v: float) -> float: # page 25
        ''' 
        ----------------------------------------------------------------------------------------------------
        Also called modulus of rigidity. Provides an approximation via readily available material properties.

        E: elastic modulus
        v: Poisson's ratio
        '''
        return E/(2*(1+v))

    def shear_strenth(self, Y: float) -> float: # page 26
        '''
        ----------------------------------------------------------------------------------------------------
        Provides an approximation of shear strength via readily available material properties.

        Y: tensile yield strength
        '''        
        return Y/2

    def hardness(self, Y: float) -> float: # page 27
        '''
        ----------------------------------------------------------------------------------------------------
        Can be used in lieu of hardness tables, as hardness and yield strength are related.

        Y: yield strength
        '''
        return 3*Y

####################################################################################

# Chapter 3

class SurfaceRoughness:
    def gaussian_dist(self, sig_s: float, z: float) -> float: # page 40
        '''
        ----------------------------------------------------------------------------------------------------
        Probability distribution function for surface roughness quantification, approximated by this Gaussian distribution.

        sig_s: standard deviation of the distribution
        z: height
        '''
        return (1/(sig_s*math.sqrt(2*math.pi)))**(-(z**2)/(2*sig_s))

    def zi_prime_mean(self, zi_prime: list): # page 41 - define zi_prime as a list of float values
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

    def zi(self, zi_prime: list, zi_prime_mean: float): # page 41
        '''
        ----------------------------------------------------------------------------------------------------
        Approximates surface heights relative to an arbitrary reference zi'.

        zi': list of measured height values *MUST BE ENTERED AS LIST*
        zi_prime_mean: reference height
        '''
        return zi_prime - zi_prime_mean

    def Ra(self, zi: list) -> float: # page 42
        '''
        ----------------------------------------------------------------------------------------------------
        Average roughness of a surface.

        zi: list of measured surface heights
        '''
        N = len(zi)
        return (1/N)*((sum(abs(zi))/N))

    def Rq(self, zi: list) -> float: # page 42
        '''
        ----------------------------------------------------------------------------------------------------
        Root-mean-square roughness aka RMS of a surface.

        zi: list of measured surface heights
        '''
        N = len(zi)
        return math.sqrt((1/N)*((sum(zi)**2)/N))

    def sig_s(self, Rq: float): # page 42
        '''
        ----------------------------------------------------------------------------------------------------
        For a Gaussian distribution of surface heights, the RMS roughness is related to the deviation of the probability distribution function sigma_s.

        Rq: RMS roughness
        '''
        return 1.25*Rq

    def Ra_composite(self, Ra1: float, Ra2: float) -> float: # page 43
        '''
        ----------------------------------------------------------------------------------------------------
        The composite average roughness of two surfaces with average roughness Ra1 and Ra2.

        Ra1: surface 1 average roughness
        Ra2: surface 2 average roughness
        '''
        return Ra1 + Ra2

    def Rq_composite(self, Rq1: float, Rq2: float) -> float: # page 43
        '''
        ----------------------------------------------------------------------------------------------------
        The composite root-mean-square roughness of two surfaces with RMS of Rq1 and Rq2.

        Rq1: surface 1 RMS roughness
        Rq2: surface 2 RMS roughness
        '''
        return math.sqrt(Rq1**2 + Rq2**2)

    def Rsk(self, Rq: float, zi: list) -> float:  # page 43
        '''
        ----------------------------------------------------------------------------------------------------
        Skewness captures the relative contribution of peaks and valleys to the average roughness. Specifically, a surface with positive skewness will have narrow, tall peaks and wide, shallow valleys; a negative skewness will have wide, short peaks and narrow, deep valleys.

        Rq: root-mean-square roughness of a given surface
        zi: list of measured heights
        '''
        N = len(zi)
        return (1/(Rq**3))*(1/N)*((sum(zi)**3)/N)

    def Rku(self, Rq: float, zi: list) -> float: # page 44
        '''
        ----------------------------------------------------------------------------------------------------
        Kurtosis reflects the sharpness of the peaks and valleys of a given surface, where the reference value is three. A kurtosis of smaller than three indicates broad peaks/valleys and larger than three means sharp peaks/valleys.

        Rq: root-mean-square roughness of a given surface
        zi: list of measured heights
        '''
        N = len(zi)
        return (1/(Rq**4))*(1/N)*((sum(zi)**4)/N)

    def autocorrelation(self, Rq: float, l: float, zi: list, zil): # page 45
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

    def exp_autocorrelation(self, xl: list, B_star: float, B: float): # page 46
        '''
        ----------------------------------------------------------------------------------------------------
        Many engineering surfaces have been found to have an exponential autocorrelation function C(xl) = exp(-xl/B). B* = 2.3B, so B can be calculated by a give B* as B = B*/2.3.

        If B is given, enter B_star as 0. If B* is given and not B, enter B* as given and B as 1000.1.

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
            return print('Please enter either B or B*, not both. If B is given, enter B_star as 0. If B* is given and not B, enter B* as given and B as 1000.1')

    def orientation(self, Bx_prime: float, By_prime: float): # page 46 !! FIX B logic !!
        '''
        ----------------------------------------------------------------------------------------------------
        An important lateral surface property is orientation. Alignment or ordering of surface features in one direction is characteristic of certain types of machining or finishing techniques, e.g. honing, used for engineering components.

        The degree of alignment or order of a surface can be quantified as an orientation parameter, where Bx* and By* are the correlation lengths of the surface in two orthogonal directions, and Bx* ≥ By* by convention.

        Bx*: larger correlation length
        By*: smaller correlation length
        '''
        return Bx_prime/By_prime

####################################################################################

# Chapter 4

class NonConformalContact:
    def area_journal_bearing(self, R, l): # page 53
        return 2*R*l

    def applied_pressure(self, W, A): # page 53
        return W/A

    def effective_radius(self, RAx, RBx, RAy, RBy): # page 54
        return (1/RAx)+(1/RBx)+(1/RAy)+(1/RBy)

    def effective_elastic_modulus(self, va, vb, Ea, Eb): # page 55
        return (1/2)*(((1-va**2)/Ea)+((1-vb**2)/Eb))

    class point: # page 57
        def a(W, Rprime, Eprime):
            return ((3*W*Rprime)/Eprime)**(1/3)

        def pmax(W, a):
            return (3*W)/(2*math.pi*(a**2))

        def pavg(W, a):
            return W/(math.pi*(a**2))

        def tmax(pmax):
            return pmax/3

        def zmax(a):
            return 0.638*a
        
    class line: # page 57
        def b(W, Rprime, l, Eprime):
            return ((4*W*Rprime)/(math.pi*l*Eprime))**(1/2)

        def pmax(W, b, l):
            return W/(math.pi*b*l)

        def pavg(W, b, l):
            return W/(4*b*l)

        def tmax(pmax):
            return 0.304*pmax

        def zmax(b):
            return 0.786*b
    
    def Aplastic(W, H): # page 59
        return W/H
    
    def Wadh(Rprime, yA, yB): # page 60
        deltaY = 2*((yA*yB)**(1/2))
        return 2*math.pi*Rprime*deltaY
    
    def Ar(Nasperities, Aasperities): # page 61
        return Nasperities*Aasperities
    
    def x_smooth_surface_approx(Rq, Rprime, a): # page 62
        return (Rq*Rprime)/(a**2)
    
    def plasticity_index(Eprime, H, Rq, r): # page 63
        return (Eprime/H)*math.sqrt(Rq/r)
    
####################################################################################

# Chapter 5

class LiquidProperties:
    def n_abs_visc(shear_stress_t, shear_strain_y): # page 68
        return shear_stress_t/shear_strain_y
    
    def kinematic_viscosity_v(n_abs_visc, fluid_density_p): # page 68
        return n_abs_visc/fluid_density_p
    
    ''' 
    page 70: log10log10(v + 0.7) = A - Blog10T
    '''

    def WilliamsLandelFerry(ng, c1, c2, T, Tg): # page 71
        return ng*(10**((-c1(T-Tg))/(c2+(T-Tg))))
    
    def Vogel(Ca, Cb, T, Cc): # page 72
        return Ca**(Cb/(T-Cc))
    
    def Barus(n0, alpha, p): # page 72
        return n0**(alpha*p)
    
    def Roelands(np, n0, p, pp, z): # page 72
        return np*((n0/np)**(((pp-p)/pp)**z))
    
    def alpha_Roelands(z, n0, np, pp): # page 73
        return (z*math.log(n0/np))/(-pp)
    
    def fx_press_temp(np, n0, pp, p, z, T0, Tinf, T, Sprime):
        return 