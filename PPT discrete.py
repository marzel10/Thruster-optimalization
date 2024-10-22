from David_work import ThrusterOptimizer1,print_results #imported the thruster class from david work
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

parameters = {
    "h": [0.5, 5, 1],  # channel height [cm]
    "w": [1, 1, 1],  # electrode width [cm]
    "R": [35, 35, 35],  # effective resistance [mOhm]
    "L": [100, 100, 100],  # effetive inductance [nH]
    "E": [3, 3, 3],  # stored energy [J]
}

constants = {
    "C1": 0.0214,
    "C2": 0.134,
    "C4": 0.154,
    "g": 9.80665,
    "nu0": 1.25663706,
    "um": 19.5856e3,  # magnetosonic speed
    "Rm": 2,  # magnetic Reynolds number
    "n": 0.8,  # solid propellant exponent

    "c": 0.5,

    "T_max": 1.2,  # 0.82 with 1.5 margin
    "I_tot": 116.5,  # 77.6  with 1.5 margin

    "E0A_max": 5,
    "N_shots_max": 1000000,
    "I_bit_max": 400,
}

THR=ThrusterOptimizer1(parameters,constants) #seting up the thruster optimizer
ah=np.linspace(parameters["h"][0],parameters["h"][1],200) #creating an array with 200 values for height starting from 0.5 and ending on 5
aw=np.linspace(parameters["w"][0],parameters["w"][1],200) #doing the same with width
Isp=[] #an array where specific impulses (for each pair of h and w from ah and aw) will be stored
ISP=0 #maximal specific impulse
H=0 #height at maximal Isp
W=0 #width at maximal Isp
ah1=[] # an array where every height inputed into get_resulst is stored (in right order)
aw1=[] # equivalent of ah1 for the width

for i in range(len(ah)): #for every h from ah
    for j in range(len(aw)): #for every w from aw
        pa=[ah[i],aw[j],35,100,3] # thruster parameters
        isp=THR.get_results1(pa)[6] #calculate specific impulse for these parameters
        Isp.append(isp)
        ah1.append(ah[i]) #add used height (in pa) to ah1
        aw1.append(aw[j]) #add used width (in pa) to ah1
        if isp>ISP: #searching for the bigest specific impulse
            ISP=isp
            H=ah[i]
            W=aw[j]
#

print(ISP,H,W) #print optimum parameters
param=[H,W,35,100,3]
opt=THR.get_results1(param) # calculate specifications for optimum parameters
print_results(param,opt) #print all specifications for optimum parameters


#making a scatter plot of every result

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Generate some random data
x = np.random.rand(100)
y = np.random.rand(100)
z = np.random.rand(100)

ax.scatter(ah1, aw1, Isp)

ax.set_xlabel('Height')
ax.set_ylabel('Width')
ax.set_zlabel('Specific Impulse')

plt.show()