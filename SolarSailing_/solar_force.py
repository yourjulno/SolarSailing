import constants as const
import numpy as np

S = 100 * 10 ** (-6) / const.AU

P = 1
A = 1

# коэффициент отражения
r = 0.1

# нормаль
n = np.array([1, 0, 0])

# касательная
t = np.array([0, 1, 0])

# eугол падения солнца
alpha = 10
cos_a = np.cos(alpha)
sin_a = np.sin(alpha)
# направление падения солнечного луча
u = cos_a * n + sin_a * t

# направление отражённого луча
s = - cos_a * n + sin_a * t
s_norm = np.linalg.norm(s)

B_f = 1  # front
B_b = 1  # back
e_f = 1  # излучение от передней поверхности
e_b = 1

f_n = P * A * ((1 + r * s_norm) * cos_a**2 + B_f * (1 - s_norm) * r * cos_a + (1 - r) *
               (e_f * B_f + e_b * B_f / (e_f + e_b)) * cos_a) * n

f_t = P * A * (1 - r * s_norm) * cos_a * sin_a * t

def solar_force() -> np.array:
    global f_t, f_n
    return f_t + f_n
