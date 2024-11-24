import numpy as np
import constants as const
import orbit_conversion as oc
# Начальные условия
r0 = np.array([5 * const.AU_in_units, 0, 0], dtype=np.float64)  # Начальная позиция (1 AU от Солнца)
v0 = np.array(
    [0, -np.sqrt(const.MU_in_units / const.AU_in_units), 0])  # Начальная скорость (орбитальная скорость Земли)


r0_norm = np.linalg.norm(r0)
v_inf = 55 / const.VU
# орбитальные элементы
a = const.MU_in_units / (v_inf ** 2)
r_p = 14 * const.R_sun_in_units
e = 1 + r_p / a
omega = np.pi
Omega = 0
i = 0
nu = np.arccos((a * (1 - e**2) - r0_norm) / (r0_norm * e))


r0_initial, v0_initial = oc.orbital_elements_to_state(a, e, i, Omega, omega, nu, const.MU_in_units)
print(r0_initial)
print(v0_initial)
m0 = [30]

# Время интегрирования
t_span = (365 * 24 * 60 * 60/ const.TU, 0)  # Один год в секундах
t_eval = np.linspace(t_span[0], t_span[1], num=1000)  # Временные точки для оценки
