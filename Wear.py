# Chapter 11

import math

def Archard(K: float, W: float, L: float, H: float):
    '''
    ----------------------------------------------------------------------------------------------------
    By far the most widely-used model for steady-state wear is the Archard equation, V = K * (WL)/H

    K: wear coefficient
    W: load
    L: sliding distance, enter 0 if calculating V/L = K * W/H
    H: hardness of the softer material
    '''
    if L != 0:
        return K*((W*L)/H)
    elif L == 0:
        return K*(W/H)