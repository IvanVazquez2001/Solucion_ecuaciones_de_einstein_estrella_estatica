import numpy as np
import matplotlib.pyplot as plt

###### LAS ECUACIONES OBTENIDAS
#### coefiecientes de la funcion incrimento de m(r) y p(r)
def RK(f,x,r):
    k1=h*f(x,r)
    k2=h*f(x+0.5*k1,r+0.5*h)
    k3=h*f(x+0.5*k2,r+0.5*h)
    k4=h*f(x+k3,r+h)
    return (k1+2.0*k2+2.0*k3+k4)/6.0
### funcion m´(r)
def dM(m,r):
    return m*(k*rho*m*r**2+1-m)/r
### funcion p´(r)
def dP(P,r):
    return -(rho+P)*W(n,w,r)/(2.0*n)
### funcion w´(r)=n´´(r)
def dW(n,w,r):
    return (4.0*r*k*P*(m**2.0)*(n**2.0)+2.0*(n**2.0)*(dM(m,r))-2.0*m*n*w+r*n*(dM(m,r))*w+r*m*w**2.0)/(2.0*r*m*n)
    
### funcion w(r)=n´(r)
def W(n,w,r):
    return w

##### la solucion para n

def RK2(f1,f2,n,w,r):
    
    k1=h*f1(n,w,r)
    l1=h*f2(n,w,r)
    k2=h*f1(n+0.5*k1,w+0.5*l1,r+0.5*h)
    l2=h*f2(n+0.5*k1,w+0.5*l1,r+0.5*h)
    k3=h*f1(n+0.5*k2,w+0.5*l2,r+0.5*h)
    l3=h*f2(n+0.5*k2,w+0.5*l2,r+0.5*h)
    k4=h*f1(n+k3,w+l3,r+h)
    l4=h*f2(n+k3,w+l3,r+h)
    return (k1+2.0*k2+2.0*k3+k4)/6.0,(l1+2.0*l2+2.0*l3+l4)/6.0

##### Numero de interaciones
a=0.0
b=.5
N=1000
h=(b-a)/N
#### Valores iniciales
k=8*np.pi
rho=10
P=0.3*rho
m=1.0
n=0.2857
w=0.0
radio=np.arange(a,b,h)
puntos_m=[]
puntos_p=[]
puntos_n=[]
for r in radio:
    if r==0:
        puntos_m.append(m)
        puntos_p.append(P)
        puntos_n.append(n)
    else:
        m = m+RK(dM,m,r)
        n = n+RK2(W,dW,n,w,r)[0]
        w = w+RK2(W,dW,n,w,r)[1]
        P = P+RK(dP,P,r)
        if P<=0:
            P=0
            rho=0
        
        puntos_n.append(n)
        puntos_m.append(m)
        puntos_p.append(P)



fig1=plt.figure()
fig2=plt.figure()

ax1=fig1.add_subplot(111)
ax2=fig2.add_subplot(111)

ax1.plot(radio,puntos_m,label='m(r)',color='blue')
ax1.plot(radio,puntos_n,label='n(r)',color='red')
ax1.plot(radio,np.array(puntos_m)*np.array(puntos_n),label='m(r)n(r)',color='green')
ax1.set_xlabel('Radio')
ax1.legend(loc='upper right')
ax1.grid(True,linestyle='--')

ax2.plot(radio,puntos_p,label='P(r)',color='crimson')
ax2.set_ylabel('P(r)')
ax2.set_xlabel('Radio')
ax2.legend(loc='upper right')
ax2.grid(True,linestyle='--')

plt.show()
