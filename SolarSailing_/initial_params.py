import numpy as np
import constants as const
# Начальные условия
r0 = np.array([const.AU, 0, 0], dtype = np.float64)  # Начальная позиция (1 AU от Солнца)
v0 = np.array([0, np.sqrt(const.MU_sun / const.AU) / 2, 0])  # Начальная скорость (орбитальная скорость Земли)

m0 = [200]

# Время интегрирования
t_span = (0, 365.25 * 24 * 60 * 60 / const.TU)  # Один год в секундах
t_eval = np.linspace(t_span[0], t_span[1], num=10_000)  # Временные точки для оценки
