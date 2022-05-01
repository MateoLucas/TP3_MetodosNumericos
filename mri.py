import numpy as np
from matplotlib import pyplot as plt

#CONSTANTES GLOBALES
mu = 4*np.pi*10**(-7) 
tol = 1e-10 #tolerancia
delta = 1e-2

def solver(L,l,n):
    # DEFINICIÓN DE FUNCIONES LAMBDA (f(x) Y f'(x))
    func = lambda r: (((r**2)*(n**2)*mu*np.pi/(l**2))*((((r**2)+(l**2))**0.5)-r))-L
    #dato: bobina mas chica: NANO Tess -> l = 1mm = 1x10**-3 m. l no llega a ser tan chica como para que df de cero, 
    # o que se pueda confundir por cero. 
    dfunc = lambda r: (np.pi*mu*n**2*r*(((r**2+l**2)**0.5)-r)*(2*((r**2+l**2)**0.5)-r))/(l**2*((r**2+l**2)**0.5))
    
    # Se busca realizar la iteración: x_(k+1) = x_k - f(x_k)/df(x_k)
    rk = 16*(10*3) * L/(n**2) # Definimos r inicial (r_0)

    #Primera iteración
    if abs(dfunc(rk)) < tol: # Si la derivada es "0", usamos el metodo de la secante
        div = (func(rk+delta)-func(rk))/delta
        j=2
        while abs(div) < tol: #En caso de que div siga siendo "0" agrandamos el intervalo
             div = (func(rk+(j*delta))-func(rk))/(j*delta)
             j+=1
    else:
        div = dfunc(rk)
    rkm1 = rk - (func(rk)/div)
    
    #Iteración 2 en adelante
    while(abs(rkm1 - rk) > tol):
        if abs(dfunc(rk)) < tol: # Si la derivada es "0", usamos el metodo de la secante
        div = (func(rk+delta)-func(rk))/delta
        j=2
        while abs(div) < tol: #En caso de que div siga siendo 0 agrandamos el intervalo
             div = (func(rk+(j*delta))-func(rk))/(j*delta)
             j+=1
    else:
        div = dfunc(rk)
    rkm1 = rk - (func(rk)/div)
        
    return rkm1



def graph(): 
    # Ploteo de r(L) L ∈ [1pH, 100uH] con l=0.2 m y N= 10,100,1000
    l= 0.2
    delta_L = 1*(10**-6)
    Ls = np.arange(10**(-9),100*10**(-6), delta_L) #Valores de L
    Ln=len(Ls)
    rs = np.empty((Ln,3)) #Se inicializa el array donde se guardaran los r(L)
    
    for N in [10,100,1000]:         
        for i in range(Ln):
            rs[i,int(np.log10(N)-1)]=solver(Ls[i],l,N) #log10(N)-1 es para ubicar la columna correspondiente (0,1,2)


    radio = plt.figure("Radio del inductor - l = 0.2 m")
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
