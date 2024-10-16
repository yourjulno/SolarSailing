import constants as const
import numpy as np

# solar pressure at one astronomical unit (AU)
P = 1

# sail area
A = 100 * 10 ** (-6) / const.AU

# коэффициент отражения
r = 0.1

# нормаль
n = np.array([1, 0, 0])

# касательная
t = np.array([0, 1, 0])

# угол падения солнца
alpha = 0
cos_a = np.cos(alpha)
sin_a = np.sin(alpha)
# направление падения солнечного луча
u = cos_a * n + sin_a * t

# направление отражённого луча
s = - cos_a * n + sin_a * t
s_norm = np.linalg.norm(s)

# non-Lambertian coefficients
B_f = 0.79  # front
B_b = 0.55  # back
e_f = 0.05  # излучение от передней поверхности
e_b = 0.55

f_n = P * A * ((1 + r * s_norm) * cos_a**2 + B_f * (1 - s_norm) * r * cos_a + (1 - r) *
               (e_f * B_f + e_b * B_f / (e_f + e_b)) * cos_a)

f_t = P * A * (1 - r * s_norm) * cos_a * sin_a

def solar_force(n: np.array) -> np.array:
    global f_t, f_n
    e = [0, 0, 1]
    t = np.cross(n, e)
    return f_t * t + f_n * n
