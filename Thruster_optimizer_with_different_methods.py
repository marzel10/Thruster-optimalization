import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


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


class ThrusterOptimizer:
    def __init__(self, params, consts, closed_form=False):
        self.params = params
        self.consts = consts

        self.guess = np.array([self.params[k][2] for k in self.params])
        self.bounds = np.array([(self.params[k][0], self.params[k][1]) for k in self.params])

    def update_params(self):
        self.guess = np.array([self.params[k][2] for k in self.params])
        self.bounds = np.array([(self.params[k][0], self.params[k][1]) for k in self.params])

    def cost_function_Isp(self, arr: np.ndarray[8]):
        return - self.get_results(arr)[6]

    def cost_function_Esp(self, arr: np.ndarray[8]):
        return - self.get_results(arr)[11]

    def get_results(self, arr: np.ndarray[8]):
        h, w, R, L, E = arr
        C1, _, C3, g, nu0, um, Rm, n, c, T_max, I_tot, _, _, _ = self.consts.values()
        V = 1720
        C = 2 * E / V ** 2

        omega = np.sqrt(1 / (L * 1e-9 * C) - (R * 1e-3) ** 2 / (4 * (L * 1e-9) ** 2))

        def get_RLC_current(t):
            I = V / (omega * L * 1e-9) * np.exp(-R * 1e-3 / (2 * L * 1e-9) * t) * np.sin(omega * t)
            return I

        m_bit = C1 * (h / w) * (E / R) * 1e3
        L_p = nu0 / np.pi * (3 / 2 + np.log((h / w) / (1 + (c / w))))  # inductance gradient
        Int_GD = sp.integrate.quad(lambda t: (get_RLC_current(t) ** 2) ** (1 / n), 0, 100 * 1e-6, limit=1000)[0]

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

    def investigate_correlation(self, var: str):
        """Function that essentially does a sensitivity study.

        Takes the name of an input variable, and varies it between its min and max values from the params dictionary.
        It keeps all other inputs as their respective guess values. Finally, it plots the results."""

        x_arr = np.linspace(self.params[var][0], self.params[var][1], 100)
        results = {key: np.zeros(100) for key in
                   ['I_EM', 'I_GD', 'I_bit', 'm_bit', 'fI', 'fm', 'I_sp', 'eta', 'R0', 'Rp', 'EpA', 'I_esp',
                    'nu_max', 'P_max', 'N_shots', 'M_prop', 'cost1', 'cost2', 'cost3']}
        arr_dict = {key: self.params[key][2] for key in self.params}

        for i in range(np.size(x_arr)):
            arr_dict[var] = x_arr[i]
            arr = np.fromiter(arr_dict.values(), dtype=float)
            results['I_EM'][i], results['I_GD'][i], results['I_bit'][i], results['m_bit'][i], results['fI'][i], \
                results['fm'][i], results['I_sp'][i], results['eta'][i], results['R0'], results['Rp'], \
                results['EpA'][i], results['I_esp'][i], \
                results['nu_max'][i], results['P_max'][i], results['N_shots'][i], results['M_prop'][
                i] = self.get_results(arr)
            # results['cost1'][i] = self.cost_function_both(arr)[0]
            # results['cost2'][i] = self.cost_function_both(arr)[1]
            # results['cost3'][i] = - self.cost_function_both(arr)[2]

        fig, ax = plt.subplots(4, 2)
        ax[0, 0].plot(x_arr, results["I_esp"], color="red")
        ax[0, 0].set_ylabel(r"$I_{sp,P} \ [\mu Ns / J]$")
        ax[0, 1].plot(x_arr, results["I_sp"], color="green")
        ax[0, 1].set_ylabel("$I_{sp,M}$ [s]")
        ax[1, 0].plot(x_arr, results['m_bit'], color='b')
        ax[1, 0].set_ylabel('$m_{bit} \ [\mu g]$')
        ax[1, 1].plot(x_arr, results['I_bit'], color='b')
        ax[1, 1].set_ylabel('$I_{bit} \ [\mu Ns]$')
        ax[2, 0].plot(x_arr, results['N_shots'], color='b')
        ax[2, 0].set_ylabel('$N_{shots} \ [-]$')
        ax[2, 1].plot(x_arr, results['M_prop'], color='b')
        ax[2, 1].set_ylabel('$M_{prop} \ [g]$')
        ax[3, 0].plot(x_arr, results['nu_max'], color='b')
        ax[3, 0].set_ylabel(r'$\nu_{max} \ [Hz]$')
        ax[3, 1].plot(x_arr, results['P_max'], color='b')
        ax[3, 1].set_ylabel('$P_{max} \ [W]$')
        fig.suptitle(f"Influence of varying {var}")

        plt.tight_layout(rect=(0, 0, 1, 1))
        plt.show()

        plt.plot(x_arr, results['cost1'], label="I_sp")
        plt.plot(x_arr, results['cost2'], label="I_esp")
        plt.plot(x_arr, results['cost3'], label="Combined")
        plt.legend()
        plt.show()

    def optimize_thruster(self, mode,met):
        """Function running the optimization using SciPy.

        It takes the lower and upper bounds, and the initial guess values, from the params dictionary.
        After finding the solution, it prints both the thruster parameters and the performance values.
        """

        if mode == "mass":
            func = self.cost_function_Isp
        elif mode == "power":
            func = self.cost_function_Esp
        # elif mode == "both":
        #     func = self.cost_function_both
        else:
            raise ValueError("Invalid mode.")

        def constraint1(arr):
            res_for_check = self.get_results(arr)
            return self.consts["E0A_max"] - res_for_check[10]

        def constraint2(arr):
            res_for_check = self.get_results(arr)
            return self.consts["N_shots_max"] - res_for_check[14]

        def constraint3(arr):
            res_for_check = self.get_results(arr)
            return self.consts["I_bit_max"] - res_for_check[2]

        cons = [{'type': 'ineq', 'fun': constraint3},
                {'type': 'ineq', 'fun': constraint2},
                {'type': 'ineq', 'fun': constraint1}]

        result = sp.optimize.minimize(func, self.guess, bounds=self.bounds, constraints=cons, method=met)
        return result

    def print_optimized_results(self):
        """Function that prints the thruster parameters and performance values for the optimized configuration."""

        opt_arr = self.optimize_thruster().x
        opt_res = self.get_results(opt_arr)

        print_results(opt_arr, opt_res)

    def print_guess_results(self):
        """Function that prints the thruster parameters and performance values for the initial guess configuration."""

        guess_res = self.get_results(self.guess)

        print_results(self.guess, guess_res)


# First and second list elements are the mininum and maxinum values for the optimization,
# and the third ones are the initial guesses (or treated as fixed values, depending on the context).
parameters = {
    "h": [0.01, 10, 1],  # channel height [cm]
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
    "I_tot":  116.5,  # 77.6  with 1.5 margin

    "E0A_max": 5,
    "N_shots_max": 1000000,
    "I_bit_max": 200,
}



with open("sim results.txt",'w') as file:
    m=["SLSQP", "trust-constr"]

    for i in range(2):
        file.write(f"this is method {m[i]}\n")

        for j in range(1,4):
            parameters["h"][2]=j
            file.write(f"Guess value for height is {parameters["h"][2]}\n")
            optimizer = ThrusterOptimizer(parameters, constants)
            param = optimizer.optimize_thruster("mass",m[i]).x
            res = optimizer.get_results(param)
            #optimizer.investigate_correlation("h")

            file.write(f"E0 = {param[4]} (discharge energy [J])\n"
            f"h = {param[0]} (height [cm])\n"
            f"w = {param[1]} (width [cm])\n"
            f"R = {param[2]} (total effective resistance [mOhm])\n"
            f"L = {param[3]} (total effective inductance [nH])\n"
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
            f"M_prop = {res[15]} (total propellant mass [g])\n")


            file.write(f"this is the end of the method {m[i]}\n\n")

    parameters["w"]=[0.01,10,2.5]
    file.write("after adding width as a degree of freedom\n\n")
    for i in range(2):
        file.write(f"this is method {m[i]}\n")

        for j in range(1,4):
            parameters["h"][2]=j
            file.write(f"Guess value for height is {parameters["h"][2]}\n")
            optimizer = ThrusterOptimizer(parameters, constants)
            param = optimizer.optimize_thruster("mass",m[i]).x
            res = optimizer.get_results(param)
            #optimizer.investigate_correlation("h")

            file.write(f"E0 = {param[4]} (discharge energy [J])\n"
            f"h = {param[0]} (height [cm])\n"
            f"w = {param[1]} (width [cm])\n"
            f"R = {param[2]} (total effective resistance [mOhm])\n"
            f"L = {param[3]} (total effective inductance [nH])\n"
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
            f"M_prop = {res[15]} (total propellant mass [g])\n")


            file.write(f"this is the end of the method {m[i]}\n\n")