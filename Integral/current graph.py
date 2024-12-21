import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def current(V,R,L,E,t,n):
    C=2*E/V**2 #calculate capacity
    w=np.sqrt(1/L/C-(R/2*L)**2) #calculate omega
    f=((V/w/L*np.exp(-(R/2/L)*t)*np.sin(w*t))**2)**(1/n)
    return f/1000#return current in kA


r=34.90568486964674E-3
l=99.98661456654209E-9
e=3.5823264230050493
v=2000


n=2

C = 2 * e / v ** 2  # calculate capacity
w = np.sqrt(1 / l / C - (r / 2 * l) ** 2)  # calculate omega



def Taylor(x,R,L,E,V):
    C = 2 * E / V ** 2  # calculate capacity
    w = np.sqrt(1 / L / C - (R / 2 * L) ** 2)  # calculate omega
    f= w * x
    f+= -(x**2 * R * w) / (2 * L)
    f += (x**3 * (R**2 * w / (8 * L**2) - w**3 / 6))
    f+= -(x**4 * (R**3 * w - 4 * L**2 * R * w**3)) / (48 * L**3)
    f+=  (w * x**5 * (5 * R**4 / L**4 - 40 * R**2 * w**2 / L**2 + 16 * w**4))
    f+= -(x**6 * (48 * L**4 * R * w**5 - 40 * L**2 * R**3 * w**3 + 3 * R**5 * w)) / 1920
    f+=  (x**7 * (-64 * L**6 * w**7 + 336 * L**4 * R * w**5 - 140 * L**2 * R**4 * w**3 + 7 * R**6 * w)) / 11520
    f+= -(x**8 * (322560 * L**6 * R * w**7 + 112 * L**4 * R**3 * w**5 - 28 * L**2 * R**5 * w**3 + 7 * R**7 * w)) / 645120
    f+=  (x**9 * (256 * L**8 * w**9 - 2304 * L**6 * R**2 * w**7 + 2016 * L**4 * R**4 * w**5 - 336 * L**2 * R**6 * w**3 + 9 * R**8 * w)) / 928972800
    f+= - (x**10 * (1280 * L**8 * R * w**9 - 3840 * L**6 * R**3 * w**7 + 2016 * L**4 * R**5 * w**5 - 240 * L**2 * R**7 * w**3 + 5 * R**9 * w)) / 928972800
    f=f*(V/w/L)**2.5
    return f/1000



print(f"The integral is {sp.integrate.quad(lambda t: current(v,r,l,e,t,0.8) , 0,  10e-5,limit=1000)}")
print(f"The integral is {sp.integrate.quad(lambda t: (v/w/l*np.exp(-r/(2*l)*t))**2.5 , 0,  10e-5,limit=1000)}")
t=np.linspace(0,0.00005,1000)
ft=Taylor(t,r,l,e,v)
p=(v/w/l*np.exp(-r/(2*l)*t))**2.5
#d=[]

d=current(v,r,l,e,t,n)
d1=current(v,r,l,e,t,1)
figure,ax=plt.subplots(1,3)
ax[0].plot(t,d)
figure.supxlabel("time[s]"); figure.supylabel("Current[kA]")
ax[1].plot(t,d1)
ax[2].plot(t,current(v,r,l,e,t,0.8),label="funtion")
#ax[2].plot(t,p,color="r",label="without sine")
titles=["n=2","n=1", "n=0.8", "fourier"]
for i in range(3):
    ax[i].set_title(titles[i])
ax[2].legend()
plt.show()
