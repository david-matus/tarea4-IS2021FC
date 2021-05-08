import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

def TrapezoideAnonima(a,b,f,e):
    '''
    Calcula el valor de la integral entre a y b de f.
    
    Parámetros:
    a: float
        valor inicial para evaluar la integral
    b: float
        valor final para evaluar la integral
    f: function
        función a integrar
    e: float
        error porcentual permitido
    
    Returns:
    I: float
        valor numérico de la integral
    '''    
    density=e/100
    I=0
    spacing=(b-a)*density
    xvalue=a
    yvalue=0
    X=[]
    Y=[]
    for m in range(int(1/density)):
        xvalue+=spacing
        yvalue=f(xvalue)
        X.append(xvalue)
        Y.append(yvalue)
    for h in range(len(X)-1):
        semiI=(X[h+1]-X[h])*(Y[h+1]+Y[h])/2
        I+=semiI
    return I

def CalcFourierC(n,L,fcond,intError):
    '''
    Calcula el valor del enésimo coeficiente de Fourier.
    
    Parámetros:
    n: int
        índice del coeficiente de Fourier
    L: float
        largo del intervalo
    fcond: function
        función a integrar
    intError: float
        error porcentual permitido
    
    Returns:
    I: float
        valor numérico de la integral
    '''    
    f=lambda x: fcond(x)*np.sin(n*x/L*np.pi)   
    Csubn=(2/L)*TrapezoideAnonima(0,L,f,intError)
    return Csubn

def CalcMallaResult(endX,endT,xspace,tspace,nFourier,intSpacing,fcond):
    '''
    Genera la matriz solución a partir del método de Fourier.
    
    Parámetros:
    xspace: float
        espaciado en x para evaluar
    tspace: float
        espaciado en t para evaluar
    endX: float
        valor final de x para evaluar
    endT: float
        valor final de y para evaluar
    nFourier: int
        cantidad de coeficientes de Fourier
    intSpacing: float
        espaciado para la integral de los coeficientes
    fcond: function
        función que evalúa la condición inicial

    Returns:
    Zmatrix: nxm matrix
        matriz Z(x,y) solución de la ecuación
    '''
    constFourier=np.zeros(nFourier)
    for i in range(nFourier):
        constFourier[i]=CalcFourierC(i+1,endX,fcond,intError)
    xstep=int(endX/xspace)
    tstep=int(endT/tspace)
    Zmatrix=np.zeros((xstep,tstep))
    for s in range(tstep):
        t=s*tspace
        for r in range(xstep):
            x=r*xspace
            for c in range(nFourier):
                Zmatrix[r,s]+=constFourier[c]*np.sin(c*np.pi*x/endX)*np.exp(-t*0.5*((c*np.pi/endX)**2))
                #Zmatrix[r,s]+=constFourier[c]*np.sin(c*np.pi*x/endX)*np.exp(-t*((c*np.pi)**2))
    return Zmatrix

endX=10
endT=0.1
xspace=0.1
tspace=0.001
nFourier=30
intError=0.01
fcond=lambda x: 2*np.exp(-((x-5)**2)/1.5)
Zmatrix=CalcMallaResult(endX,endT,xspace,tspace,nFourier,intError,fcond)
X = np.arange(0, endX, xspace)
T = np.arange(0, endT, tspace)
T, X = np.meshgrid(T, X)


fig,ax=plt.subplots(subplot_kw={"projection": "3d"})
surf=ax.plot_surface(T,X,Zmatrix)
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter('{x:.02f}')
plt.show()
