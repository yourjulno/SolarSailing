import constants as const
import numpy as np
import math

c_light = 299_792_458  # m/s
E_Earth = 1372.  # Wt/m^2 = kg / s^3
P_Earth = E_Earth / c_light  # kg / (m * s^2)

# solar pressure at one astronomical unit (AU)
P_au = P_Earth * const.TU ** 2 / (const.AU * 1e3)  # N_in_units / m^2

# sail area
A = 1000.364755174347  # m^2

# коэффициент отражения
r = 0.91

# non-Lambertian coefficients
B_f = 0.79  # front
B_b = 0.67  # back
e_f = 0.025  # излучение от передней поверхности
e_b = 0.27
s = 0.89

a_c = 0.06 / (const.CU * 1000000)


# управление, psi - угол между v и r
def control(psi: float, delta: float) -> float:
    cos_psi = np.cos(psi)
    sqrt_cos_psi = cos_psi * np.sqrt(8 + cos_psi ** 2)
    cos2_a = 8 / 3 * (3 + sqrt_cos_psi) / (12 + cos_psi ** 2 + sqrt_cos_psi)
    cos_a = np.sqrt(cos2_a)
    return np.arccos(cos_a)


def solar_force(r_i: np.array, v_i: np.array, psi: float, sail_mass: float, delta: float, a: float) -> np.array:
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

    angle1 = control(psi, delta)
    # TODO: считать в control черед угол дельта
    angle2 = 0
    # нормаль к поверхности паруса в орбитальной СК
    n = np.array([np.cos(angle1) * np.cos(angle2), np.cos(angle1) * np.sin(angle2), np.sin(angle2)])
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
    m_del_A = P_au / ((1 + r * s) + B_f * (1 - s) * r + (1 - r) *
                      ((e_f * B_f - e_b * B_b) / (e_f + e_b))) * a
    A_del_m = 1 / m_del_A
    f_n = P * ((1 + r * s) * cos_a ** 2 + B_f * (1 - s) * r * cos_a + (1 - r) *
                ((e_f * B_f - e_b * B_b) / (e_f + e_b)) * cos_a)
    f_t = P * (1 - r * s) * cos_a * sin_a
    f_e = f_n * n + f_t * np.array([0, 1, 0])  # Add tangential component
    f_i = np.dot(S, f_e)
    return f_i * A_del_m
