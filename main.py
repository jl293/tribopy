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
    def friction(self, F, W): # page 13
        return F/W

    def lambda_ratio(self, h0, Ra): # page 15
        return h0/Ra

    def hersey_number(self, w, n, pa): # page 17
        return (w*n)/pa

####################################################################################

# Chapter 2 

class SolidMaterials:
    def shear_modulus(self, E, v): # page 25
        return E/(2*(1+v))

    def shear_strenth(self, Y): # page 26
        return Y/2

    def hardness(self, Y): # page 27
        return 3*Y

####################################################################################

# Chapter 3

class SurfaceRoughness:
    def gaussian_dist(self, sig_s, z): # page 40
        return (1/(sig_s*math.sqrt(2*math.pi)))**(-(z**2)/(2*sig_s))

    def zprimei_mean(self, zprime_mean, N, zprimei: list): # page 41 - define zprimei as a list of float values
        if len(zprimei) == 0:
            print("Please define zprimei as a list (ex. [4.1, 3.77, ..., n])")
            return 0
        if N == 0:
            print("WARNING: N = 0")
            return 0
        return (1/N)*(sum(zprimei)/N)

    def zi(self, zprimei, zprimei_mean): # page 41
        return zprimei - zprimei_mean

    def Ra(self, N, zi: list): # page 42
        return (1/N)*((sum(abs(zi))/N))

    def Rq(self, N, zi: list): # page 42
        return math.sqrt((1/N)*((sum(zi)**2)/N))

    def sig_s(self, Rq): # page 42
        return 1.25*Rq

    def Ra_composite(self, Ra1, Ra2): # page 43
        return Ra1 + Ra2

    def Rq_composite(self, Rq1, Rq2): # page 43
        return math.sqrt(Rq1**2 + Rq2**2)

    def Rsk(self, Rq, N, zi): # page 43
        return (1/(Rq**3))*(1/N)*((sum(zi)**3)/N)

    def Rku(self, Rq, N, zi): # page 44
        return (1/(Rq**4))*(1/N)*((sum(zi)**4)/N)

    def autocorrelation(self, Rq, N, l, zi, zil): # page 45
        return (1/((Rq**2)*(N-l)))*(sum(zi*zil)/(N-l))

    def exp_autocorrelation(self, xl, B): # page 46
        return 1**(-xl/B)

    def orientation(self, Bx_prime, By_prime): # page 46
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