# Chapter 2 

def shear_modulus(E: float, v: float) -> float: # page 25
    ''' 
    ----------------------------------------------------------------------------------------------------
    Also called modulus of rigidity. Provides an approximation via readily available material properties.

    E: elastic modulus
    v: Poisson's ratio
    '''
    return E/(2*(1+v))

def shear_strenth(Y: float) -> float: # page 26
    '''
    ----------------------------------------------------------------------------------------------------
    Provides an approximation of shear strength via readily available material properties.

    Y: tensile yield strength
    '''        
    return Y/2

def hardness(Y: float) -> float: # page 27
    '''
    ----------------------------------------------------------------------------------------------------
    Can be used in lieu of hardness tables, as hardness and yield strength are related.

    Y: yield strength
    '''
    return 3*Y