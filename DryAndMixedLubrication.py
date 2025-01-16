# Chapter 10

import math

def fc(Fc: float, W: float): # page 142
    '''
    ----------------------------------------------------------------------------------------------------
    Contact friction force Fc (also called Coulomb friction) increases with normal force W, and the ratio of contact friction to the normal force is the contact friction coefficient, fc.

    Fc: Coulomb or contact friction force
    W: normal force
    '''
    return Fc/W

def fc_sum(fa: float, fd: float): # page 143
    '''
    ----------------------------------------------------------------------------------------------------
    Contact friction coefficient fc can also be approximated with the sum of the adhesive and deformation friction coefficients fa and fd, respectively. The relative effects of adhesion and deformation friction are determined by the materials in contact, the surface roughness, and the operating conditions.

    fa: adhesive friction coefficient
    fd: deformation friction coefficient
    '''
    return fa + fd

def fa_normal(Fa: float, W: float): # page 144
    '''
    ----------------------------------------------------------------------------------------------------
    The adhesive coefficient fa is the ratio of the ratio of adhesive friction force Fa to load W, and the adhesive friction force increase with the size and strength of the asperity junctions.

    Fa: adhesive friction force
    W: load
    '''
    return Fa/W

def fa_SS(Ss: float, Ar: float, W: float): # page 144
    '''
    ----------------------------------------------------------------------------------------------------
    The adhesive coefficient fa is the ratio of the ratio of adhesive friction force Fa to load W, and the adhesive friction force increase with the size and strength of the asperity junctions.

    Ss: shear strength of the weaker material
    Ar: real contact area
    W: load
    '''
    return(Ss*Ar)/W

def fa_plastic(Ss: float, H: float): # page 144
    '''
    ----------------------------------------------------------------------------------------------------
    For plastic contacts, the contact area is simply the ratio of the shear strength of the weaker material and the hardness of the softer material.

    NOTE: fa_plastic is typically ~~ 0.17.

    Ss: shear strength of the softer material
    H: hardness of the softer material
    '''
    return Ss/H

def fa_elastic(Ss: float, EPrime: float, r: float, sigma_s: float): # page 145
    '''
    ----------------------------------------------------------------------------------------------------
    For elastic contacts, adhesive friction is determined by the elasticity and shear strength of the materials in contact as well as the surface roughness. Assuming a Gaussian distribution of surface heights, the adhesive friction for elastic contacts can be estimated via this function fa_elastic.

    Ss: shear strength of the weaker material
    EPrime: effective elastic modulus of the materials
    r: average radius of curvature of the asperities
    sigma_s: standard deviation (sigma_s) of the surface height distribution that can be measured or approximated from the root-mean square roughness
    '''
    return (Ss/EPrime) * (((math.pi*r)/sigma_s)**0.5)

def Fd(H: float, A0side: float):
    '''
    ----------------------------------------------------------------------------------------------------
    Deformation force Fd on a given asperity is the product of hardness of the softer material H and the contact area projected from the side, A0, side.

    H: hardness of the softer material
    A0side: contact area projected from the side
    '''
    return H*A0side

def fd(Fd: float, W: float, H: float, A0side: float, A0top: float, theta: float): # page 146
    '''
    ----------------------------------------------------------------------------------------------------
    fd can be calculated as Fd/W, (H * A0,side) / (H * A0,top), or (2/pi)*cot(theta).
    
    One possible assumed asperity shape is a cone with angle theta. For this shape, the deformation friction coefficient shall be calculated as (2/pi)*cot(theta).

    NOTE: For this function, enter the fewest possible variables and enter 0 for the variables you do not have given/measured.

    Fd: deformation force
    W: load
    H: hardness of the softer material
    A0side: projected area from the side
    A0top: projected area from the top
    theta: angle of asperity in degrees
    '''
    if Fd != 0 and W != 0:
        print("Using fd = Fd/w")
        return Fd/W
    elif (Fd == 0 or W) == 0 and H != 0 and A0top != 0 and A0side != 0:
        print("Using fd = (H * A0,side) / (H * A0,top)")
        return (H * A0side) / (H * A0top)
    elif Fd == 0 and W == 0 and H == 0 and A0top == 0 and A0side == 0 and theta != 0:
        try:
            if theta > 0:
                theta = math.radians(theta)
            print("Using fd = (2/pi)*cot(theta)")
            return (2/math.pi)*(1/math.tan(theta))
        except Exception as e:
            print(f"Error {e} ... NOTE: theta = {theta}")
    else:
        print("Properly define variables and try again.")

def fd_ceramic(C: float, KIC: float, E: float, H: float, W: float): # page 147
    '''
    ----------------------------------------------------------------------------------------------------
    For brittle materials, deformation occurs through fracture, so ceramic deformation friction is a function of fracture toughness K_IC. It has been shown that friction between ceramics increases rapidly with increasing fracture toughness according to the relationship given in this function.

    fd,ceramic = C * ((KIC^2) / (E*sqrt(H*W)))
    
    C: constant of proportionality; obtained by fitting measured friction vs load data by using the equation
    KIC: fracture toughness
    E: elastic modulus of the ceramic being fractured
    H: hardness of the ceramic bering fractured
    W: load
    '''
    return C*((KIC**2)/(E*math.sqrt(H*W)))

def fm(xi: float, fb: float, fv: float): # page 152
    '''
    ----------------------------------------------------------------------------------------------------
    The total friction coefficient in mixed librication fm is due to the contributions of the viscous friction coefficient fv and the boundary friction coefficient fb.

    xi: ξ is the load support ratio
    fb: boundary friction coefficient; typically around 0.1
    fv: viscous friction coefficient
    '''
    return xi*fb + (1-xi)*fv

def xi(l: float): # page 153
    '''
    ----------------------------------------------------------------------------------------------------
    The load support ratio ξ (xi) has been found to vary with λ (lambda). When λ is small, at the transition between boundary and mixed lubrication, boundary friction is dominant and ξ = 1. When λ is large, at the onset of full film lubrication, viscous friction is dominant and ξ = 0.

    l: lambda ratio (λ), see Fundamental Concepts
    '''
    m: float
    n: float
    m = 1.517
    n = 1.327
    return 1/((1 + l**m)**n)

def Pe(U: float, Cp: float, p: float, a: float, k: float): # page 153
    '''
    ----------------------------------------------------------------------------------------------------
    The dimensionless Peclet number is an empirical estimate of temperature rise based on the relative rates of heat generated by friction and heat dissipated away from the contact.

    A low Peclet number arises for slow moving surfaces or for materials that rapidly conduct heat away from the contact, with the opposite being true for high Peclet numbers.

    U: speed
    Cp: specific heat
    p: density
    a: characterization of the contact (e.g. Hertz radius)
    k: thermal conductivity
    '''
    return (U*Cp*p*a)/(2*k)

def AverageFlashTemp(Pe: float, f: float, W: float, UA: float, UB: float, k: float, a: float): # pages 153 and 154
    '''
    ----------------------------------------------------------------------------------------------------
    The average flash temperature at point contacts depends on the magnitude of the Peclet number. The friction coefficient f can be for dry contact, boundary, or mixed conditions (fc, fb, or fm).

    Pe: Peclet number
    f: friction coefficient fc, fb, or fm
    W: load
    UA: speed of bottom surface
    UB: speed of top surface
    k: thermal conductivity
    a: characterization of the contact
    '''
    try:
        a = float(a)
    except Exception as e:
        if a == "Pe":
            pass
        elif a != "Pe":
            print(f'Error {e}. Please enter a valid number or the string "Pe" to calculate the contact characterization value.')    
    if Pe < 0.1:
        print('Using formula for Pe < 0.1')
        return 0.25*((f*W*abs(UA-UB))/(k*a))
    elif 0.1 < Pe < 5:
        try:
            af: float
            af01: float
            af01 = 0.85
            af5: float
            af5 = 0.35
            Pe01: float
            Pe01 = 0.1
            Pe5: float
            Pe5 = 5
            af = af01 + ((Pe-Pe01)*(af5-af01))/(Pe5-Pe01)
            print('Using formula for 0.1 < Pe < 5')
            return (0.25*af)*((f*W*abs(UA-UB))/(k*a))
        except Exception as e:
            print(f'Error: {e}')
    elif Pe > 5:
        try:
            Cp = float(input("Please enter the specific heat of the material's flash temperature being calculated: "))
            p = float(input("Please enter the density of the material's flash temperature being calculated: "))
        except Exception as e:
            print(f'Error: {e}. Please make sure you enter a valid number.')
        print(f'Using formula for Pe > 5')
        return 0.308 * ((f*W*abs(UA-UB))/(k*a)) * (k/(UA*Cp*p*a))**0.5
    
def fr(shape: bool, xi: float, ab: float, R: float): # page 155
    '''
    ----------------------------------------------------------------------------------------------------
    The rolling friction coefficient fr is used for either a cylinder or sphere:
    For a rolling cylinder of radius R with a contact half-width b.
    For a rolling sphere of radius R and contact half-width a

    shape: boolean True for a cylinder, False for a sphere
    xi: load support ratio, ξ; typically ~0.02 for steel surfaces
    ab: Hertz contact half-width, b-value for a cylinder and a-value for a sphere
    R: cylinder radius
    '''
    if shape == True:
        print("Using formula for a cylinder")
        return (2*xi*ab)/(3*math.pi*R)
    elif shape == False:
        print("Using formula for a sphere")
        return (3*xi*ab)/(16*R)