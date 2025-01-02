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

    def orientation(self, Bx_prime: float, Bx: float, xl_B: float, By_prime: float, By: float, xl_Y: float): # page 46
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

####################################################################################

# Chapter 4

class NonConformalContact:
    def area_journal_bearing(self, R: float, l: float): # page 53
        '''
        ----------------------------------------------------------------------------------------------------
        The most common component that exhibits conformal contact is the journal bearing. For journal bearings, area is approximated from the dimensions of the bearing that would be seen from the side, called the projected area. The projected area A for a journal bearing of length l and radius R is:

        A = 2Rl

        R: radius of given journal bearing
        l: length of given journal bearing
        '''
        return 2*R*l

    def applied_pressure(self, W: float, A: float): # page 53
        '''
        ----------------------------------------------------------------------------------------------------
        The applied pressure of a journal bearing is the load divided by the projected area.

        W: load of givenjournal bearing
        A: projected area of given journal bearing
        '''
        return W/A

    def effective_radius(self, RAx: float, RBx: float, RAy: float, RBy: float): # page 54
        '''
        ----------------------------------------------------------------------------------------------------
        The contact between spheroids A and B is approximated by the effective radius 1/R'

        Radii of curvature for representative shapes are as follows:
         ________________________________
        |Shape       Rx          Ry      |
        |--------------------------------|
        |Cuboid      infinity    infinity|
        |                                |
        |Sphere      R           R       |
        |                                |
        |Cylinder    R           infinity|
        |                                |
        |Ring        -R          infinity|
         ________________________________

        NOTE: If you need to enter an infinity value, enter it as math.inf

        RAx: Radius of curvature for spheroid A, x-direction
        RAy: Radius of curvature for spheroid A, y-direction
        RBx: Radius of curvature for spheroid B, x-direction
        RBy: Radius of curvature for spheroid B, y-direction
        '''
        return (1/RAx)+(1/RBx)+(1/RAy)+(1/RBy)

    def effective_elastic_modulus(self, va: float, vb: float, Ea: float, Eb: float): # page 55
        '''
        ----------------------------------------------------------------------------------------------------
        The elastic deformation of contacting bodies is also a function of their material properties, specifically their elasticity. This is quantified by the elastic modulus E and the Poisson's ratio v. The effective modulus of elasticity between two contacting bodies can be found with these parameters.

        NOTE: Even if the two contacting bodies are made of the same material, the effective elastic modulus is NOT the same as the elastic modulus of the material itself.
        '''
        return (1/2)*(((1-va**2)/Ea)+((1-vb**2)/Eb))

    class point: # page 57
        '''
        ----------------------------------------------------------------------------------------------------
        Point contact refers to a scenario such as two ball bearings touching each other on a single point. See the class "line" for functions pertaining to parallel cylinders contacting each other in a rectangular patch.
        '''
        def a(self, W: float, Rprime: float, Eprime: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Hertz contact equation for half width a: point where two spheres contact one another.

            W: load of contact surface
            R': effective radius
            E': effective modulus of elasticity
            '''
            return ((3*W*Rprime)/Eprime)**(1/3)

        def pmax(self, W: float, a: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the maximum contact pressure between the contacting bodies.

            W: load of contact surface
            a: Hertz half width for point contact
            '''
            return (3*W)/(2*math.pi*(a**2))

        def pavg(self, W: float, a: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the average contact pressure between the contacting bodies.

            W: load of contact surface
            a: Hertz half width for point contact
            '''
            return W/(math.pi*(a**2))

        def tmax(self, pmax: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the magnitude of the maximum subsurface shear stress for point contact.

            pmax: Maximum contact pressure between the contacting bodies.
            '''
            return pmax/3

        def zmax(self, a: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the position od the maximum subsurface shear stress for point contact.

            a: Hertz half width for point contact
            '''
            return 0.638*a
        
    class line: # page 57
        '''
        ----------------------------------------------------------------------------------------------------
        Line contact refers to functions pertaining to parallel cylinders contacting each other in a rectangular patch. See the class "point" for functions that refer to scenarios such as two ball bearings touching each other on a single point. 
        '''
        def b(self, W: float, Rprime: float, l: float, Eprime: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Hertz contact equation for half width b: line contact where two spherical bodies lie parallel with each other.
            '''
            return ((4*W*Rprime)/(math.pi*l*Eprime))**(1/2)

        def pmax(self, W: float, b: float, l: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the average contact pressure between the contacting bodies.

            W: load of contact surface
            b: Hertz half width for line contact
            l: 1/2 of the total length of the contact patch
            '''
            return W/(math.pi*b*l)

        def pavg(self, W: float, b: float, l: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the average contact pressure between the contacting bodies.

            W: load of contact surface
            b: Hertz half width for line contact
            l: 1/2 of the total length of the contact patch
            '''
    
            return W/(4*b*l)

        def tmax(self, pmax: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the magnitude of the maximum subsurface shear stress for line contact.

            pmax: Maximum contact pressure between the contacting bodies.
            '''
            return 0.304*pmax

        def zmax(self, b: float):
            '''
            ----------------------------------------------------------------------------------------------------
            Represents the position od the maximum subsurface shear stress for polineint contact.

            b: Hertz half width for line contact
            '''
            return 0.786*b
    
    def Aplastic(self, W: float, H: float): # page 59
        '''
        ----------------------------------------------------------------------------------------------------
        For a fully plastic contact, the contact area is approximated as the ratio of load to hardness. The onset of fully plastic regime can be approximated by pavg ~~ 2.8Y.

        W: load of contact surface
        H: hardness of softer material
        '''
        return W/H
    
    def Wadh(self, Rprime: float, yA: float, yB: float): # page 60
        '''
        ----------------------------------------------------------------------------------------------------
        Adhesive force will cause the contact area to increase above that calculated using the Hertz equations. There are multiple models for calculating the area of adhesive contacts, but the simplest approach is to add the adhesive force to the applied load and then use the Hertz equations.

        R': effective radius
        yA: surface energy of surface A
        yB: surface energy of surface B
        '''
        deltaY = 2*((yA*yB)**(1/2))
        return 2*math.pi*Rprime*deltaY
    
    def Ar(self, Nasperities: float, Aasperity: float): # page 61
        '''
        ----------------------------------------------------------------------------------------------------
        The total area of contact between asperities is the real contact area, Ar, also known as the true or actual contact area. The real contact area is the sum of the areas of contact between individual asperities.

        Nasperities: the number of asperities in contact
        Aasperity: the average contact area of individual asperities
        '''
        return Nasperities*Aasperity
    
    def x_smooth_surface_approx(self, Rq: float, Rprime: float, a: float): # page 62
        '''
        ----------------------------------------------------------------------------------------------------
        For non-conformal contacts, the Hertz contact equations are often used even for rough surfaces. When taking this approach, it is useful to first evaluate the accuracy (or inaccuracy) of the smooth surface approximation that is inherent to the Hertz models. This is done via this function.

        Rq: composite root-mean-square roughness (Ra can also be used)
        R': effective radius of the bodies in contact
        a: Hertz contact radius
        '''
        return (Rq*Rprime)/(a**2)
    
    def plasticity_index(self, Eprime: float, H: float, Rq: float, r: float): # page 63
        '''
        ----------------------------------------------------------------------------------------------------
        This model assumes that two rough surfaces can be approximated by a smooth rigid plane in contact with a rough elastic surface. The roughness is represented by spherical asperities whose heights follow some statistical distribution. If it is assumed that the asperities deform independently of each other in response to load, the Hertz contact model can be applied to each asperity to calculate local deformation of the areas for individual asperities.

        Greenwood and Williamson developed this criteria to determine if rough surface contact was elastic or plastic. If this function returns a value ψ > 1, the surface will primarily deform plastically, and elastically for ψ < 0.6.

        E': effective modulus of elasticity
        H: hardness of the softer surface
        Rq: composite root-mean-square roughness
        r: average radius of curvature of the surface asperities
        '''
        x = (Eprime/H)*math.sqrt(Rq/r)
        if x > 1:
            print(f'{x} is greater than 1, so the surfaces will mostly deform plastically.')
        if x < 0.6:
            print(f'')
        return x
        
    
####################################################################################

# Chapter 5

class LiquidProperties:
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
    
    def PSSI(self, vf: float, vs: float, vb: float):
        '''
        ----------------------------------------------------------------------------------------------------
        The shear stability of a lubricant is characterized by the permanent shear stability index (PSSI).

        The vsheared is measured after the lubricant has been subject to rapid, high shear stress using one of several standard tests.

        vf: vfresh is the kinematic viscosity of the lubricant before shearing at 100°C
        vs: vsheared is the kinematic viscosity after shearing at 100°C
        vb: vbase is the viscosity of the base oil of the lubricant at 100°C
        '''
        return (vf-vs)/(vf-vb) * 100
    
    def Carreau(self, ninf: float, n0: float, y: float, ycr: float, m: float, n: float):
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

    def SpecificGravity(self, t: float, p: float, ph20: float):
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