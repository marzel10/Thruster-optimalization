import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def current(V,R,L,E,t,n):
    C=2*E/V**2
    w=np.sqrt(1/L/C-(R/2*L)**2)
    f=(V/w/L*np.exp(-(R/2/L)*t)*np.sin(w*t))**(2/n)
    return f/1000


def get_RLC_current(t,V,L,R,E):
    C = 2 * E / V ** 2
    omega=w=np.sqrt(1/L/C-(R/2*L)**2)
    I = V / (omega * L * 1e-9) * np.exp(-R * 1e-3 / (2 * L * 1e-9) * t) * np.sin(omega * t)
    return I
v=2000
r=24E-3
l=150E-9
e=8
n=0.8

print(sp.integrate.quad(lambda t: (get_RLC_current(t,v,l,r,e) ** 2) ** (1 / n), 0, 10e-6))
t=np.linspace(0,0.00005,1000)

#d=[]
# for i in range(len(t)):
#     d.append(current(v,r,l,e,t[i],n))
d=current(v,r,l,e,t,n)
d1=current(v,r,l,e,t,1)
figure,ax=plt.subplots(1,3)
ax[0].plot(t,d)
ax[1].plot(t,d1)
ax[2].plot(t,current(v,r,l,e,t,0.8))
#plt.xlabel("time[s]"); plt.ylabel("Current[kA]")

plt.show()
