import numpy as np
import constants as const
# Начальные условия
r0 = np.array([7 * const.AU_in_units, 0, 0], dtype = np.float64)  # Начальная позиция (1 AU от Солнца)
v0 = np.array([0, -np.sqrt(const.MU_in_units / const.AU_in_units), 0])  # Начальная скорость (орбитальная скорость Земли)

m0 = [100]

# Время интегрирования
t_span = (0, 365 * 24 * 60 * 60 * 10/ const.TU)  # Один год в секундах
t_eval = np.linspace(t_span[0], t_span[1], num=1000)  # Временные точки для оценки
