import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def current(V,R,L,E,t,n):
    C=2*E/V**2 #calculate capacity
    w=np.sqrt(1/L/C-(R/2*L)**2) #calculate omega
    f=((V/w/L*np.exp(-(R/2/L)*t)*np.sin(w*t))**2)**(1/n)
    return f/1000 #return current in kA



v=2000
r=24E-3
l=150E-9
e=3
n=2

print(f"The integral is {sp.integrate.quad(lambda t: current(v,r,l,e,t,0.8) , 0,  10e-5,limit=1000)}")

t=np.linspace(0,0.00005,1000)

#d=[]

d=current(v,r,l,e,t,n)
d1=current(v,r,l,e,t,1)
figure,ax=plt.subplots(1,3)
ax[0].plot(t,d)
figure.supxlabel("time[s]"); figure.supylabel("Current[kA]")
ax[1].plot(t,d1)
ax[2].plot(t,current(v,r,l,e,t,0.8))
titles=["n=2","n=1", "n=0.8"]
for i in range(3):
    ax[i].set_title(titles[i])
plt.show()
