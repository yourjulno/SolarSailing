import constants as const
import numpy as np
import math

c_light = 299_792_458  # m/s
E_Earth = 1372.  # Wt/m^2 = kg / s^3
P_Earth = E_Earth / c_light  # kg / (m * s^2)

# solar pressure at one astronomical unit (AU)
P_au = P_Earth * const.TU ** 2 / (const.AU * 1e3)  # N_in_units / m^2

# sail area
A = 594.364755174347  # m^2

# коэффициент отражения
r = 0.91

# non-Lambertian coefficients
B_f = 0.79  # front
B_b = 0.67  # back
e_f = 0.025  # излучение от передней поверхности
e_b = 0.27

s = 0.89

# управление, psi - угол между v и r
def control(psi: float) -> float:
    tan_psi = np.tan(psi)
    tan_alpha = (-3 + np.sqrt(9 + 8 * tan_psi**2)) / (4 * tan_psi)
    return np.arctan(tan_alpha)

def solar_force(r_i: np.array, v_i: np.array, psi: float, sail_mass: float) -> np.array:
    global r, s, B_f, B_b, e_f, e_b, A, P_au
    c_i = np.cross(r_i, v_i)

    c_i_norm = np.linalg.norm(c_i)
    r_norm = np.linalg.norm(r_i)
    # базисные векторы орбитальной ск в инерциальной ск
    e1 = r_i / r_norm
    e3 = c_i / c_i_norm
    e2 = np.cross(e3, e1)
    # матрица перехода из инерциальной в орбитальную СК
    S = np.array([e1, e2, e3]).T

    angle1 = control(psi)
    angle2 = 0
    # нормаль к поверхности паруса в орбитальной СК
    n = np.array([np.cos(angle1), 0, np.sin(angle2)])
    # направление падения солнечного луча в орбитальной СК
    u = np.array([1, 0, 0])
    cos_a = np.dot(n, u)
    if cos_a < 0:
        print("ERROR IN SAIL ORIENTATION")
    # n = - n
    sin_a = np.sqrt(1 - cos_a ** 2)
    # солнечное давление
    P = P_au * (const.AU_in_units / r_norm) ** 2
    # вектор касательной к парусу

    t = np.cross(n, np.cross(u, n))

    # Force components
    PA = P * A
    f_n = PA * ((1 + r * s) * cos_a**2 + B_f * (1 - s) * r * cos_a + (1 - r) *
                ((e_f * B_f - e_b * B_b) / (e_f + e_b)) * cos_a)
    f_t = PA * (1 - r * s) * cos_a * sin_a
    f_e = f_n * n + f_t * np.array([0, 1, 0])  # Add tangential component
    f_i = np.dot(S, f_e)
    return f_i / sail_mass
