import constants as const
import numpy as np

# solar pressure at one astronomical unit (AU)
P_au = 1

# sail area
A = 100 * 10 ** (-6) / const.AU

# коэффициент отражения
r = 0.99

# non-Lambertian coefficients
B_f = 0.79  # front
B_b = 0.55  # back
e_f = 0.05  # излучение от передней поверхности
e_b = 0.55

def solar_force(r_i: np.array, v_i: np.array) -> np.array:
    c_i = np.cross(r_i, v_i)
    c_i_norm = np.linalg.norm(c_i)
    r_norm = np.linalg.norm(r_i)
    # базисные векторы в связанной СК
    e1 = r_i / r_norm
    e3 = c_i / c_i_norm
    e2 = np.cross(e3, e1)
    # матрица перехода из инерциальной в собственную СК
    S = np.array([e1, e2, e3])
    fi_1 = 0
    fi_2 = 1
    # нормаль к поверхности паруса в собственной СК
    n = np.array([np.cos(fi_1) * np.cos(fi_2), np.cos(fi_1) * np.sin(fi_2), np.sin(fi_1)])
    # направление падения солнечного луча в собственной СК
    u = np.array([1, 0, 0])
    un = np.dot(n, u)
    if un < 0:
        n = - n
    alpha = np.dot(n, u)
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    # солнечное давление
    P = P_au / r_norm ** 2
    # вектор касательной к парусу
    t = np.cross(n, np.cross(u, n))

    # направление отражённого луча
    s = - cos_a * n + sin_a * t
    s_norm = np.linalg.norm(s)
    f_n = P * A * ((1 + r_i * s_norm) * cos_a ** 2 + B_f * (1 - s_norm) * r * cos_a + (1 - r) *
                   (e_f * B_f + e_b * B_f / (e_f + e_b)) * cos_a)

    f_t = P * A * (1 - r * s_norm) * cos_a * sin_a
    f_e = f_t * t + f_n * n
    f_i = np.dot(S, f_e)
    return f_i

