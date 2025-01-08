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

    For this function, enter the fewest possible variables and enter 0 for the variables you do not have given/measured.

    Fd
    W:
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