import numpy as np
from scipy.optimize import fsolve

# Constants
q = 1.602e-19  # Charge of an electron (C)
k = 1.38e-23  # Boltzmann's constant (J/K)
T = 300  # Temperature in Kelvin (example value)

# Parameters (example values)
I_L = 0.01  # Light-generated current (A)
I_0 = 1e-10  # Reverse saturation current (A)
n = 1.5  # Ideality factor
R_sh = 100  # Shunt resistance (Ohms)

# Equation to solve
def equation(Voc):
    return I_L - I_0 * (np.exp(q * Voc / (n * k * T)) - 1) - Voc / R_sh

# Initial guess for Voc
initial_guess = 0.6  # Volts

# Solve for Voc
Voc_solution = fsolve(equation, initial_guess)

print(f'Solved Voc: {Voc_solution[0]} V')
