# Chapter 12

import math

def fstatic(formula: bool, Fstatic: float, W: float, theta: float): # page 178
    '''
    1. Fstatic/W
    2. tan(theta)

    formula: True for formula 1, False for formula 2
    Fstatic: static friction force
    W: load
    theta: angle at which block starts sliding down the plane
    '''
    if theta > 0:
        theta = math.radians(theta)
    if formula == True:
        return Fstatic/W
    elif formula == False:
        return math.tan(theta)
    else:
        print('Enter "True" for formula 1 or "False" for formula 2.')