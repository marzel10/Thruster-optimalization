import numpy as np
from numpy import unravel_index
import scipy as sp
import matplotlib.pyplot as plt
from scipy import optimize

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
    "I_tot":  116.5,  # 77.6  with 1.5 margin

    "E0A_max": 5,
    "N_shots_max": 1000000,
    "I_bit_max": 200,
}

#Function calculating Isp for the optimalization (from David's code)
def Specific_impulse(x, E=3,constants=constants):
    #x is a vector with all the inputs
    #x[0] is the height of the truster
    #x[1] is the width of the truster
    R=x[2] #Resistance is 3rd entry of the vector x
    L=x[3] #Inductance is 4th entry of the vector x

    #rest of the function goes the same as david's get_results
    C1, _, C3, g, nu0, um, Rm, n, c, T_max, I_tot, _, _, _ = constants.values()
    V = 1720
    C = 2 * E / V ** 2

    omega = np.sqrt(1 / (L * 1e-9 * C) - (R * 1e-3) ** 2 / (4 * (L * 1e-9) ** 2))

    def get_RLC_current(t):
        I = V / (omega * L * 1e-9) * np.exp(-R * 1e-3 / (2 * L * 1e-9) * t) * np.sin(omega * t)
        return I

    m_bit = C1 * (x[0] / x[1]) * (E / R) * 1e3
    L_p = nu0 / np.pi * (3 / 2 + np.log((x[0] / x[1]) / (1 + (c / x[1]))))  # inductance gradient
    Int_GD = sp.integrate.quad(lambda t: (get_RLC_current(t) ** 2) ** (1 / n), 0, 10e-5, limit=1000)[0]
    #print(Int_GD)
    I_EM = E / R * (L_p / 2) * 1e3
    I_GD = C3 * 1e-9 * x[0] / x[1] ** (2 / n - 1) * Int_GD * 1e7
    I_bit = I_EM + I_GD
    I_sp = I_bit * 1e-6 / (m_bit * 1e-9 * g)
    fm = I_EM * 1e-6 / (um * m_bit * 1e-9)  # electromagnetic mass fraction
    fI = I_EM / I_bit  # electromagnetic impulse fraction
    eta = (I_bit * 1e-6) ** 2 / (2 * m_bit * 1e-9 * E)

    # Not sure what to do about plasma inductance
    Rp = (L_p * 10 ** -6 * um / Rm) * 1e3  # plasma resistance
    R0 = R - Rp  # neccessary circuit resistance for effective resistance R

    EpA = E / (x[0] * x[1])  # energy per propellant surface area, to avoid charring

    I_esp = I_bit / E

    nu_max = T_max / I_bit
    P_max = nu_max * E
    N_shots = I_tot / I_bit * 1e6
    M_prop = N_shots * m_bit * 1e-6

    #maximal number of shots and maximal impulse bit was eliminated as it gave weird reults
    #if N_shots>1000000 or EpA>5 or I_bit>200 :
        #I_sp=0

    if EpA>5 : #to avoid charring
       I_sp=0


    return -I_sp #optimalization will find a minimum, but we are looking for maximum so there is - in front of Isp

#Function entirely take from David's code
def get_results( h, w, R, L, E,constants=constants):

    C1, _, C3, g, nu0, um, Rm, n, c, T_max, I_tot, _, _, _ = constants.values()
    V = 1720
    C = 2 * E / V ** 2

    omega = np.sqrt(1 / (L * 1e-9 * C) - (R * 1e-3) ** 2 / (4 * (L * 1e-9) ** 2))

    def get_RLC_current(t):
        I = V / (omega * L * 1e-9) * np.exp(-R * 1e-3 / (2 * L * 1e-9) * t) * np.sin(omega * t)
        return I

    m_bit = C1 * (h / w) * (E / R) * 1e3
    L_p = nu0 / np.pi * (3 / 2 + np.log((h / w) / (1 + (c / w))))  # inductance gradient
    Int_GD = sp.integrate.quad(lambda t: (get_RLC_current(t) ** 2) ** (1 / n), 0, 10e-5, limit=1000)[0]
    #print(Int_GD)
    I_EM = E / R * (L_p / 2) * 1e3
    I_GD = C3 * 1e-9 * h / w ** (2 / n - 1) * Int_GD * 1e7
    I_bit = I_EM + I_GD
    I_sp = I_bit * 1e-6 / (m_bit * 1e-9 * g)
    fm = I_EM * 1e-6 / (um * m_bit * 1e-9)  # electromagnetic mass fraction
    fI = I_EM / I_bit  # electromagnetic impulse fraction
    eta = (I_bit * 1e-6) ** 2 / (2 * m_bit * 1e-9 * E)

    # Not sure what to do about plasma inductance
    Rp = (L_p * 10 ** -6 * um / Rm) * 1e3  # plasma resistance
    R0 = R - Rp  # neccessary circuit resistance for effective resistance R

    EpA = E / (h * w)  # energy per propellant surface area, to avoid charring


    I_esp = I_bit / E

    nu_max = T_max / I_bit
    P_max = nu_max * E
    N_shots = I_tot / I_bit * 1e6
    M_prop = N_shots * m_bit * 1e-6



    return I_EM, I_GD, I_bit, m_bit, fI, fm, I_sp, eta, R0, Rp, EpA, I_esp, nu_max, P_max, N_shots, M_prop

def print_results(input, res):
    print(
        f"E0 = {input[4]} (discharge energy [J])\n"
        f"h = {input[0]} (height [cm])\n"
        f"w = {input[1]} (width [cm])\n"
        f"R = {input[2]} (total effective resistance [mOhm])\n"
        f"L = {input[3]} (total effective inductance [nH])\n"
        f"I_EM = {res[0]} (electromagnetic impulse bit [microNs])\n"
        f"I_GD = {res[1]} (gas dynamic impulse bit [microNs])\n"
        f"I_bit = {res[2]} (total impulse bit [microNs])\n"
        f"m_bit = {res[3]} (total mass bit [microg])\n"
        f"fI = {res[4]} (electromagnetic impulse fraction [-])\n"
        f"fm = {res[5]} (electromagnetic mass fraction [-])\n"
        f"I_sp = {res[6]} (mass specific impulse [s])\n"
        f"eta = {res[7]} (efficiency [-])\n"
        f"R0 = {res[8]} (circuit resistance [mOhm])\n"
        f"Rp = {res[9]} (plasma resistance [mOhm])\n"
        f"E0A = {res[10]} (energy to propellant surface area ratio [J/cm^2])\n"
        f"I_esp = {res[11]} (energy specific impulse [microNs/J])\n"
        f"nu_max = {res[12]} (maximum frequency [Hz])\n"
        f"P_max = {res[13]} (maximum power [W])\n"
        f"N_shots = {res[14]} (number of shots [-])\n"
        f"M_prop = {res[15]} (total propellant mass [g])"
    )


v=2000
r=22
l=45
e=3
n=0.8
bounds=[(0.5,2),(0.5,1.5),(22, 55),(45,170)]  # Bounds for h[cm],w[cm],R[mOhm],L[nH] respectively


#finding minimum via continous analysis
opt=optimize.shgo(Specific_impulse, bounds)
print("Continous optimization goes first, then the discrete one\n")
print(f"Results of the continous optimization: ")
res=np.zeros(5)
res[:4]=opt["x"]
res[4]=e

print_results(np.array(res),get_results(res[0],res[1],res[2],res[3],res[4]))#print results for the optimal configuration

print("\n")

#discrete analysis
H=np.arange(0.5,2,0.01)
W=np.arange(0.5,1.5,0.01)
H,W=np.meshgrid(H,W)


ISPP1=get_results(H,W,r,l,e)[6]
ISPP1=np.where(ISPP1>0,ISPP1,0)

h_index,w_index=unravel_index(np.argmax(ISPP1),np.shape(ISPP1))
print("Results of the discrete optimization: ")
print_results(np.array([H[h_index, w_index],W[h_index, w_index],r,l,e]),get_results(H[h_index, w_index],W[h_index, w_index],r,l,e))#print results for the optimal configuration


#ploting 3d function
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(H, W, ISPP1)
ax.set_xlabel('Height')
ax.set_ylabel('Width')
ax.set_zlabel('Specific impulse')
plt.show()



