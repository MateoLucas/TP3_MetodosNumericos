import numpy as np
import scipy as sc
from matplotlib import pyplot as plt

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


print(solver(0,0.2,100))  #BORRAR

def graph():
    #ploteo de valores
    l= 0.2
    delta_L = 1*(10**-6)
    Ls = np.arange(10**(-9),100*10**(-6), delta_L)
    Ln=len(Ls)
    rs = np.empty((Ln,3))
    for N in [10,100,1000]:         
        for i in range(Ln):
            rs[i,int(np.log10(N)-1)]=solver(Ls[i],l,N)


    radio = plt.figure("Radio")
    ax = radio.add_subplot(3,1,1)
    plt.xlabel("L[H]")
    plt.ylabel("r(L)[m] - N=10")
    plt.plot(Ls,rs[:,0],color = "#40eb34")

    radio.add_subplot(3,1,2)
    plt.xlabel("L[H]")
    plt.ylabel("r(L)[m] - N=100")
    plt.plot(Ls,rs[:,1],color = "#f0fc03")

    radio.add_subplot(3,1,3)
    plt.xlabel("L[H]")
    plt.ylabel("r(L)[m] - N=1000")
    plt.plot(Ls,rs[:,2],color = "#f803fc")

    plt.show()
    return

#graph()