# Chapter 12

import math

def fstatic(formula: bool, Fstatic: float, W: float, theta: float): # page 178
    '''
    This is used to calculate the static force coefficient.  The static friction force, Fstatic, and the force in the direction normal to the tilted plane W are related geometrically to the angle theta and the weight of the block times gravity. Then, the static friction coefficient, fstatic, can be calculated with one of two formulas:

    1. Fstatic/W
    2. tan(theta)

    and is essentially a function of the measured angle theta. The inclined plane test does not provide the accuracy or versatility needed in most cases. However, it is a simple and inexpensive was to measure static friction between two flate surfaces and good demonstration of the concept of static friction.

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