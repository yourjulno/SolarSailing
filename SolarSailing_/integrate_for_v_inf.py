import numpy as np
from derivatives import derivatives
from scipy.integrate import solve_ivp
import constants as const
import orbit_conversion as oc
import matplotlib.pyplot as plt
import events as ev
import initial_params as ip

# Массив различных значений v_inf
v_inf_values = [45, 50, 55, 58, 60]  # Примерные значения для v_inf

# Подготовка пустых массивов для сохранения данных по графикам
a_arr_all = []
omega_arr_all = []
e_arr_all = []
nu_arr_all = []
times_all = []
trajectory_all = []


# Функция для интегрирования для разных v_inf
def integrate_for_v_inf(v_inf):
    # Пересчитываем начальные условия на основе v_inf
    r_p = 15 * const.R_sun_in_units  # Перицентр
    a = const.MU_in_units / (v_inf ** 2)
    e = 1 + r_p / a
    omega = np.pi
    Omega = 0
    i = 0
    nu = np.arccos((a * (1 - e ** 2) - ip.r0_norm) / (ip.r0_norm * e))

    # Преобразуем орбитальные элементы в начальные условия
    r0_initial, v0_initial = oc.orbital_elements_to_state(a, e, i, Omega, omega, nu, const.MU_in_units)

    # Объединение начальных условий в один массив
    initial_conditions = np.hstack((r0_initial, v0_initial))

    # Время интегрирования
    t_span = (365 * 24 * 60 * 60 * 4 / const.TU, 0)  # Один год в секундах
    t_eval = np.linspace(t_span[0], t_span[1], num=1000)  # Временные точки для оценки

    # Решение системы дифференциальных уравнений
    solution = solve_ivp(
        derivatives,
        t_span,
        initial_conditions,
        method='RK45',
        t_eval=t_eval,
        rtol=1e-8,
        atol=1e-8,
        events=[ev.detect_pericenter, ev.detect_apocenter]  # Указываем оба события
    )

    # Извлечение координат
    x = solution.y[0]
    y = solution.y[1]
    z = solution.y[2]
    vx = solution.y[3]
    vy = solution.y[4]
    vz = solution.y[5]
    times = solution.t  # Времена
    trajectory = solution.y.T

    a_arr = []
    omega_arr = []
    e_arr = []
    nu_arr = []

    # Преобразование в орбитальные элементы для всех точек траектории
    for i in range(len(solution.y[0])):
        r = np.array([x[i], y[i], z[i]])
        v = np.array([vx[i], vy[i], vz[i]])
        a_, e_, omega_, nu_ = oc.state_to_orbital_elements(r, v, const.MU_in_units)
        a_arr.append(a_)
        omega_arr.append(omega_ * 180 / np.pi)
        e_arr.append(e_)
        nu_arr.append(nu_)

    # Сохранение результатов для текущего v_inf
    a_arr_all.append(a_arr)
    omega_arr_all.append(omega_arr)
    e_arr_all.append(e_arr)
    nu_arr_all.append(nu_arr)
    times_all.append(times)
    trajectory_all.append(trajectory)


# Интегрируем для всех значений v_inf
for v_inf in v_inf_values:
    integrate_for_v_inf(v_inf / const.VU)

# Построение 3D графиков для всех траекторий
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(0, 0, 0, color='yellow', s=200, label='Солнце')

# Отображаем траектории для разных значений v_inf
colors = ['pink', 'blue', 'green', 'red', 'black']  # Разные цвета для разных траекторий
for i, trajectory in enumerate(trajectory_all):
    ax.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2],
            label=f'Trajectory for v_inf = {v_inf_values[i]} km/s', color=colors[i], linewidth=1)
    # Отображаем начальные положения спутников
    ax.scatter(trajectory[0][0], trajectory[0][1], trajectory[0][2], color=colors[i], s=140,
               label=f'Начальное положение {v_inf_values[i]} км/с')

# Настройка графика
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_xlim(-3, 3)  # Ограничиваем ось X
ax.set_ylim(-3, 3)  # Ограничиваем ось Y
ax.set_zlim(0, 2)  # Ограничиваем ось Z
ax.set_title('Trajectories for different v_inf')
ax.grid(True)
ax.legend()

plt.show()
plt.show()
# Построение графиков для разных значений v_inf
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 15))

# График зависимости большой полуоси от времени для разных v_inf
for i, v_inf in enumerate(v_inf_values):
    ax1.plot(times_all[i], a_arr_all[i], label=f'v_inf = {v_inf} км/с')
ax1.set_xlabel('Время, дни')
ax1.set_ylabel('Большая полуось, км')
ax1.set_title('Зависимость большой полуоси от времени')
ax1.legend()
ax1.grid(True)

# График зависимости аргумента перицентра от времени для разных v_inf
for i, v_inf in enumerate(v_inf_values):
    ax2.plot(times_all[i], omega_arr_all[i], label=f'v_inf = {v_inf} км/с')
ax2.set_xlabel('Время, дни')
ax2.set_ylabel('Аргумент перицентра, градусы')
ax2.set_title('Зависимость аргумента перицентра от времени')
ax2.legend()
ax2.grid(True)

plt.savefig('images/momentum_energy_for_v_inf')
# plt.show()
