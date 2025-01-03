# Chapter 8

import sympy

def Reynolds1D(h: float, p: float, x: float, UA: float, UB: float, n: float):
    '''
    ----------------------------------------------------------------------------------------------------
    The two fundamental equations of fluid dynamics that describe fkuid flow through a lubricated gap are the Navier-Stokes and continuity equations. Combining these equations for lubricant flow and applying a series of simplifying expressions yields a reduced form fo what is called the Reynolds equation, of which this is the 1D form.

    h: fluid film thickness
    p: pressure
    x: position in the direction of motion
    UA: relative speed of bottom plate
    UB: relative speed of top plate
    n: nu is viscosity
    '''
    U: float
    U = abs(UA-UB)

    

    return