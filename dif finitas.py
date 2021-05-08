import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

def MatrixCalc(xspace,yspace,xend,yend,fcond,z0,zf,alpha,itertotal):
    '''
    Genera la matriz solución a partir del método de diferencias finitas.
    
    Parámetros:
    xspace: float
        espaciado en x para evaluar
    yspace: float
    xend: float
        valor final de x para evaluar
    yend: float
        valor final de y para evaluar
    fcond: function
        una función que evalúa la condición inicial
    z0: float
        valor de Z en (0,0)
    zf: float
        valor de z en (xend,0)
    alpha:
        parámetro de la ecuación de calor
    itertotal:
        cantidad de veces en las que se iterará Z
    
    Returns:
    Zmatrix: nxm matrix
        matriz Z(x,y) solución de la ecuación
    '''
    xstep=int(xend/xspace)
    ystep=int(yend/yspace)
    Zmatrix=np.zeros((xstep,ystep))
    r=alpha*yspace/(xspace**2)
    for u in range(xstep):
        Zmatrix[u][0]=fcond(u*xspace)
    for v in range(ystep):
        Zmatrix[0][v]=z0
        Zmatrix[xstep-1][v]=zf
    for itcurrent in range(itertotal):
        for y in range(1,ystep):
            for x in range(1,xstep-1):
                Zmatrix[x][y]=r*(Zmatrix[x+1][y-1]+Zmatrix[x-1][y-1])+(1-2*r)*Zmatrix[x][y-1]
    return Zmatrix
#se definen los parámetros
xend=10
tend=15
xspace=0.1
tspace=0.01
z0=0
zf=0
fcond=lambda x: 2*np.exp(-((x-5)**2)/1.5)
alpha=0.5
itertotal=5
#se calcula la solución
Zmatrix=MatrixCalc(xspace,tspace,xend,tend,fcond,z0,zf,alpha,itertotal)
X=np.arange(0,xend,xspace)
T=np.arange(0,tend,tspace)
T,X=np.meshgrid(T,X)
#se grafica la matriz
fig,ax=plt.subplots(subplot_kw={"projection": "3d"})
surf=ax.plot_surface(T,X,Zmatrix)
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter('{x:.02f}')
plt.show()

