# Импорты
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import constants as const
import initial_params as ip
import gravity_force as gf
import thrust_force as tf
import solar_force as sf


# Функция для вычисления производной
def derivatives(t, y):
    r_i = y[:3]  # Позиция
    v_i = y[3:6]  # Скорость

    # Ускорение
    # a_trac = tf.traction_acceleration(m, e)
    a_grav = gf.gravity_acceleration(r_i)

    # Угол между r и v
    psi = np.arctan2(np.linalg.norm(np.cross(r_i, v_i)), np.dot(r_i, v_i))
    a_solar = sf.solar_force(r_i, v_i, psi, ip.m0)

    a = a_grav + a_solar

    # print(np.sign(np.dot(a_solar, v_i)))

    return np.hstack((v_i, a))


# Объединение начальных условий в один массив
initial_conditions = np.hstack((ip.r0, ip.v0))

# Решение системы дифференциальных уравнений
solution = solve_ivp(derivatives, ip.t_span, initial_conditions,
                     method='RK45', t_eval=ip.t_eval, rtol=1e-8, atol=1e-8)

# Извлечение координат
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]
vx = solution.y[3]
vy = solution.y[4]
vz = solution.y[5]
times = solution.t  # Времена

trajectory = solution.y.T

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
# Plot the Sun as a sphere
# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)
# x_sun = 0.1 * np.outer(np.cos(u), np.sin(v))  # Adjust Sun's radius for visualization
# y_sun = 0.1 * np.outer(np.sin(u), np.sin(v))
# z_sun = 0.1 * np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_surface(x_sun, y_sun, z_sun, color='orange', alpha=0.7, label='Sun')
ax.scatter(0, 0, 0, color='yellow', s=100, label='Солнце')
ax.scatter(trajectory[0][0], trajectory[0][1], trajectory[0][2], color='pink', s=100, label='Начальное положение спутника')
# Plot the trajectory
ax.plot(trajectory[:, 0] , trajectory[:, 1] , trajectory[:, 2] , label='Solar Sail')

# Labels and title
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_title('Trajectory')
ax.grid()
plt.show()
# # # Создаем массивы векторов положений и скоростей
# r = np.vstack((x, y, z)).T  # Транспонируем для получения правильной формы
# v = np.vstack((vx, vy, vz)).T
#
# # Расчёт первых интегралов
# distances = np.sqrt(x ** 2 + y ** 2 + z ** 2)
# speeds = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
#
# # Энергия
# total_energy = speeds ** 2 - 2 * const.MU_in_units / distances
#
# # Угловой момент
# c = np.cross(r, v)
# c_norm = np.linalg.norm(c, axis=1)

# # Создание фигуры и 3D-оси
# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')
#
# # Отрисовка орбиты
# ax.plot(x, y, z, label='Орбита спутника', color='blue')
#
# # Отрисовка Солнца в начале координат
# ax.scatter(0, 0, 0, color='yellow', s=100, label='Солнце')
# # # Установка меток на осях
# # ax.set_xticks(np.linspace(x[0], x[-1], 6))
# # ax.set_yticks(np.linspace(y[0], y[-1], 6))
# # ax.set_zticks(np.linspace(z[0], z[-1], 6))
# # Настройка сетки и меток
# ax.grid(True)
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
#
# # Добавление легенды
# ax.legend()
# plt.savefig('images/orbit')
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 15))
#
# # График интеграла энергии
# ax1.plot(times, total_energy, label='Энергия', color='purple')
# ax1.set_xlabel('Время, дни')
# ax1.set_ylabel('Энергия, VU^2')
# ax1.set_title('Зависимость интеграла энергии от времени')
# ax1.legend()
# ax1.grid(True)
#
# # График орбитального момента
# ax2.plot(times, c_norm, label='Момент', color='red')
# ax2.set_xlabel('Время, дни')
# ax2.set_ylabel('Момент')
# ax2.set_title('Зависимость орбитального момента от времени')
# ax2.legend()
# ax2.grid(True)
# plt.savefig('images/momentum_energy')
# plt.show()
