import numpy as np
import scipy as sc
#CONSTANTES GLOBALES
mu = 4*np.pi*10**(-7) 
f = lambda r,n,L,l: (((r**2)(n**2)*mu*np.pi/(l**2))*((((r**2)+(l**2))**0.5)-r))-L
df = lambda r,n,l: (np.pi*mu*n^2*r*(((r^2+l^2)**0.5)-r)*(2*((r^2+l^2)**0.5)-r))/(l^2*((r^2+l^2)**0.5))

def solver(L,l,n):
    return r
