# Chapter 4

import math

def area_journal_bearing(R: float, l: float): # page 53
    '''
    ----------------------------------------------------------------------------------------------------
    The most common component that exhibits conformal contact is the journal bearing. For journal bearings, area is approximated from the dimensions of the bearing that would be seen from the side, called the projected area. The projected area A for a journal bearing of length l and radius R is:

    A = 2Rl

    R: radius of given journal bearing
    l: length of given journal bearing
    '''
    return 2*R*l

def applied_pressure(W: float, A: float): # page 53
    '''
    ----------------------------------------------------------------------------------------------------
    The applied pressure of a journal bearing is the load divided by the projected area.

    W: load of givenjournal bearing
    A: projected area of given journal bearing
    '''
    return W/A

def effective_radius(RAx: float, RBx: float, RAy: float, RBy: float): # page 54
    '''
    ----------------------------------------------------------------------------------------------------
    The contact between spheroids A and B is approximated by the effective radius 1/R'. If you need to enter an infinity value, enter it as math.inf. Radii of curvature for representative shapes are as follows:

    |Shape   |    Rx      |    Ry      |
    |--------|------------|------------|
    |Cuboid  |    infinity|    infinity|
    |Sphere  |    R       |    R       |
    |Cylinder|    R       |    infinity|
    |Ring    |    -R      |    infinity|
    
    RAx: Radius of curvature for spheroid A, x-direction
    RAy: Radius of curvature for spheroid A, y-direction
    RBx: Radius of curvature for spheroid B, x-direction
    RBy: Radius of curvature for spheroid B, y-direction
    '''
    return (1/RAx)+(1/RBx)+(1/RAy)+(1/RBy)

def effective_elastic_modulus(va: float, vb: float, Ea: float, Eb: float): # page 55
    '''
    ----------------------------------------------------------------------------------------------------
    The elastic deformation of contacting bodies is also a function of their material properties, specifically their elasticity. This is quantified by the elastic modulus E and the Poisson's ratio v. The effective modulus of elasticity between two contacting bodies can be found with these parameters. Even if the two contacting bodies are made of the same material, the effective elastic modulus is NOT the same as the elastic modulus of the material itself.
    '''
    return (1/2)*(((1-va**2)/Ea)+((1-vb**2)/Eb))

class point: # page 57
    '''
    ----------------------------------------------------------------------------------------------------
    Point contact refers to a scenario such as two ball bearings touching each other on a single point. See the class "line" for functions pertaining to parallel cylinders contacting each other in a rectangular patch.
    '''
    def a(W: float, Rprime: float, Eprime: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Hertz contact equation for half width a: point where two spheres contact one another.

        W: load of contact surface
        R': effective radius
        E': effective modulus of elasticity
        '''
        return ((3*W*Rprime)/Eprime)**(1/3)

    def pmax(W: float, a: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the maximum contact pressure between the contacting bodies.

        W: load of contact surface
        a: Hertz half width for point contact
        '''
        return (3*W)/(2*math.pi*(a**2))

    def pavg(W: float, a: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the average contact pressure between the contacting bodies.

        W: load of contact surface
        a: Hertz half width for point contact
        '''
        return W/(math.pi*(a**2))

    def tmax(pmax: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the magnitude of the maximum subsurface shear stress for point contact.

        pmax: Maximum contact pressure between the contacting bodies.
        '''
        return pmax/3

    def zmax(a: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the position of the maximum subsurface shear stress for point contact.

        a: Hertz half width for point contact
        '''
        return 0.638*a
    
class line: # page 57
    '''
    ----------------------------------------------------------------------------------------------------
    Line contact refers to functions pertaining to parallel cylinders contacting each other in a rectangular patch. See the class "point" for functions that refer to scenarios such as two ball bearings touching each other on a single point. 
    '''
    def b(W: float, Rprime: float, l: float, Eprime: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Hertz contact equation for half width b: line contact where two spherical bodies lie parallel with each other.
        '''
        return ((4*W*Rprime)/(math.pi*l*Eprime))**(1/2)

    def pmax(W: float, b: float, l: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the average contact pressure between the contacting bodies.

        W: load of contact surface
        b: Hertz half width for line contact
        l: 1/2 of the total length of the contact patch
        '''
        return W/(math.pi*b*l)

    def pavg(W: float, b: float, l: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the average contact pressure between the contacting bodies.

        W: load of contact surface
        b: Hertz half width for line contact
        l: 1/2 of the total length of the contact patch
        '''

        return W/(4*b*l)

    def tmax(pmax: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the magnitude of the maximum subsurface shear stress for line contact.

        pmax: Maximum contact pressure between the contacting bodies.
        '''
        return 0.304*pmax

    def zmax(b: float):
        '''
        ----------------------------------------------------------------------------------------------------
        Represents the position od the maximum subsurface shear stress for polineint contact.

        b: Hertz half width for line contact
        '''
        return 0.786*b

def Aplastic(W: float, H: float): # page 59
    '''
    ----------------------------------------------------------------------------------------------------
    For a fully plastic contact, the contact area is approximated as the ratio of load to hardness. The onset of fully plastic regime can be approximated by pavg ~~ 2.8Y.

    W: load of contact surface
    H: hardness of softer material
    '''
    return W/H

def Wadh(Rprime: float, yA: float, yB: float): # page 60
    '''
    ----------------------------------------------------------------------------------------------------
    Adhesive force will cause the contact area to increase above that calculated using the Hertz equations. There are multiple models for calculating the area of adhesive contacts, but the simplest approach is to add the adhesive force to the applied load and then use the Hertz equations.

    R': effective radius
    yA: surface energy of surface A
    yB: surface energy of surface B
    '''
    deltaY = 2*((yA*yB)**(1/2))
    return 2*math.pi*Rprime*deltaY

def Ar(Nasperities: float, Aasperity: float): # page 61
    '''
    ----------------------------------------------------------------------------------------------------
    The total area of contact between asperities is the real contact area, Ar, also known as the true or actual contact area. The real contact area is the sum of the areas of contact between individual asperities.

    Nasperities: the number of asperities in contact
    Aasperity: the average contact area of individual asperities
    '''
    return Nasperities*Aasperity

def x_smooth_surface_approx(Rq: float, Rprime: float, a: float): # page 62
    '''
    ----------------------------------------------------------------------------------------------------
    For non-conformal contacts, the Hertz contact equations are often used even for rough surfaces. When taking this approach, it is useful to first evaluate the accuracy (or inaccuracy) of the smooth surface approximation that is inherent to the Hertz models. This is done via this function.

    Rq: composite root-mean-square roughness (Ra can also be used)
    R': effective radius of the bodies in contact
    a: Hertz contact radius
    '''
    return (Rq*Rprime)/(a**2)

def plasticity_index(Eprime: float, H: float, Rq: float, r: float): # page 63
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