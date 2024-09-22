# Импорты
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import constants as const
import initial_params as ip
from gravity_force import gravity_force
from traction_force import traction_force
# Функция для вычисления производной
def derivatives(t, y):
    global MU_sun
    r = y[:3]  # Позиция
    v = y[3:6]  # Скорость

    # Ускорение
    a = gravity_force(r) + traction_force()

    return np.concatenate((v, a))


# Объединение начальных условий в один массив
initial_conditions = np.concatenate((ip.r0, ip.v0))

# Решение системы дифференциальных уравнений
solution = solve_ivp(derivatives, ip.t_span, initial_conditions,
                     method='RK45', t_eval=ip.t_eval, atol=1e-12, rtol=1e-12)

# Извлечение координат
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]
vx = solution.y[3]
vy = solution.y[4]
vz = solution.y[5]
times = solution.t  # Времена

# Создаем массивы векторов положений и скоростей
r = np.vstack((x, y, z)).T  # Транспонируем для получения правильной формы
v = np.vstack((vx, vy, vz)).T

# Расчёт первых интегралов
distances = np.sqrt(x ** 2 + y ** 2 + z ** 2)
speeds = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)

# Энергия
total_energy = speeds ** 2 - 2 * const.MU_sun / distances

# Угловой момент
c = np.cross(r, v)

c_norm = np.linalg.norm(c, axis=1)

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

# Отрисовка орбиты
ax.plot(x, y, label='Орбита спутника', color='blue')

# Отрисовка Солнца в начале координат
ax.scatter(0, 0, color='yellow', s=100, label='Солнце')
ax.grid(True)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 15))

# График интеграла энергии
ax1.plot(times, total_energy, label='Энергия', color='purple')
ax1.set_xlabel('Время, дни')
ax1.set_ylabel('Энергия, VU^2')
ax1.set_title('Зависимость интеграла энергии от времени')
ax1.legend()
ax1.grid(True)

# График орбитального момента
ax2.plot(times, c_norm, label='Момент', color='red')
ax2.set_xlabel('Время, дни')
ax2.set_ylabel('Момент')
ax2.set_title('Зависимость орбитального момента от времени')
ax2.legend()
ax2.grid(True)
plt.show()
plt.savefig('momentum_energy')
