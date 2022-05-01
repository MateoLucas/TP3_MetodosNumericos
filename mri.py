import numpy as np
import scipy as sc

#CONSTANTES GLOBALES
mu = 4*np.pi*10**(-7) 
tol = 1e-10 #tolerancia

def solver(L,l,n):
    func = lambda r: (((r**2)*(n**2)*mu*np.pi/(l**2))*((((r**2)+(l**2))**0.5)-r))-L
    #dato: bobina mas chica: NANO Tess -> l = 1mm = 1x10**-3 m. l no llega  aser tan chica como para que df de cero, 
    # o que se pueda confundir por cero. 
    dfunc = lambda r: (np.pi*mu*n**2*r*(((r**2+l**2)**0.5)-r)*(2*((r**2+l**2)**0.5)-r))/(l**2*((r**2+l**2)**0.5))
    # Lo que buscamos: x_(k+1) = x_k - f(x_k)/df(x_k)
                  #definimos r inicial
    rk = 16*(10*3) * L/(n**2) 
    rkm1 = rk - (func(rk)/dfunc(rk))
    while(abs(rkm1 - rk) > tol):
        rk = rkm1
        rkm1 = rk - (func(rk)/dfunc(rk))
    return rkm1


print(solver(1e-9,0.2,10))  #BORRAR