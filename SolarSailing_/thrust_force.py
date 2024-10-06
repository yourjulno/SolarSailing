import constants as const
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import initial_params as ip

# скорость истечения топлива
u = 12753 / const.VU # км/день

# максимальная сила тяги
f_max = 20e-6 / (const.MU * const.DU / (const.TU ** 2))

gamma = f_max / u


# coefficients - коэфициенты полинома
def traction_acceleration(m: float, e: np.array) -> np.array:
    global f_max
    return f_max * e / m
