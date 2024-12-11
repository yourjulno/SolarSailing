import numpy as np
from derivatives import derivatives
from scipy.integrate import solve_ivp
import constants as const
import orbit_conversion as oc
import matplotlib.pyplot as plt
import events as ev
import plotly.graph_objects as go

# Массив различных значений v_inf
v_inf_values = [20, 30, 45]  # Примерные значения для v_inf

# Подготовка пустых массивов для сохранения данных по графикам
a_arr_all = []
omega_arr_all = []
e_arr_all = []
nu_arr_all = []
times_all = []
trajectory_all = []


# Функция для интегрирования для разных v_inf
def integrate_for_v_inf(v_inf):
    r0_norm = 5 * const.AU_in_units
    # Пересчитываем начальные условия на основе v_inf
    r_p = 10 * const.R_sun_in_units  # Перицентр
    a = -const.MU_in_units / (v_inf ** 2)
    e = 1 - r_p / a
    omega = 180 / 180 * np.pi
    Omega = 0
    i = 0
    nu = -np.arccos((a * (1 - e**2) - r0_norm) / (r0_norm * e))
    print(f"a = {a}")
    print(f"e = {e}")
    print(f"nu = {nu}")
    # Преобразуем орбитальные элементы в начальные условия
    r0_initial, v0_initial = oc.orbital_elements_to_state(a, e, i, Omega, omega, nu, const.MU_in_units)
    print(f"r0 = {r0_initial}")
    print(f"v0 = {v0_initial}")
    print(np.linalg.norm(v0_initial))
    # Объединение начальных условий в один массив
    initial_conditions = np.hstack((r0_initial, v0_initial))

    # Время интегрирования
    t_span = (365 * 24 * 60 * 60 * 10 / const.TU, 0)  # Один год в секундах
    t_eval = np.linspace(t_span[0], t_span[1], num=10000)  # Временные точки для оценки

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


# Отображаем траектории для разных значений v_inf
colors = ['pink', 'blue', 'green', 'red', 'black']  # Разные цвета для разных траекторий

fig = go.Figure()
traces = []


for i, trajectory in enumerate(trajectory_all):
    trace1 = go.Scatter(x = trajectory[:, 0],
                          y = trajectory[:, 1],
                          name = f'Trajectory for v_inf = {v_inf_values[i]} km/s',
                          mode = "lines",
                          line = dict(width = 3,
                                      color = colors[i]))
    traces.append(trace1)

layout = dict(title = dict(text = "TITLE"),
              width = 1800,
              height = 1000,
              margin = dict(l=100, r=100, b=100, t=100),
              xaxis = dict(title = "XLABEL"),
              yaxis = dict(title = "YLABEL"),
              legend = dict(yanchor="top", # Defines whether bot/top of the legend box is the anchor
                            y=0.99,
                            xanchor="left", # Defines whether left/right of the legend box is the anchor
                            x=0.99),
              )

sun = go.Scatter(x = [0],
                  y = [0],
                  name = f'SUN',
                  mode = "markers",
                  marker = dict(size = 50,
                                color = "YELLOW"))
start = go.Scatter(x = [trajectory_all[0][0]],
                  y = [trajectory_all[0][1]],
                  name = f'SATELLITE',
                  mode = "markers",
                  marker = dict(size = 50,
                                color = "BLUE"))
fig.add_trace(sun)
fig.add_trace(start)
fig.add_traces(traces)
fig.update_layout(layout)

fig.show()

# fig = go.Figure()
# traces = []
# for i, trajectory in enumerate(trajectory_all):
#     trace1 = go.Scatter3d(x = trajectory[:, 0],
#                           y = trajectory[:, 1],
#                           z = trajectory[:, 2],
#                           name = f'Trajectory for v_inf = {v_inf_values[i]} km/s',
#                           mode = "lines",
#                           line = dict(width = 3,
#                                       color = colors[i]))
#     traces.append(trace1)
#
# layout = dict(title = dict(text = "TITLE"),
#               width = 2400,
#               height = 800,
#               margin = dict(l=100, r=100, b=100, t=100),
#               xaxis = dict(title = "x, AU"),
#               yaxis = dict(title = "y, AU"),
#               # zaxis = dict(title = "z, AU"),
#               legend = dict(yanchor="top", # Defines whether bot/top of the legend box is the anchor
#                             y=0.99,
#                             xanchor="right", # Defines whether left/right of the legend box is the anchor
#                             x=0.99),
#               )
# sun = go.Scatter3d(x = [0],
#                 y = [0],
#                           z = [0],
#                           name = f'SUN',
#                           mode = "markers",
#                           marker = dict(size = 3,
#                                         color = "YELLOW"))
# fig.add_trace(sun)
# fig.add_traces(traces)
# fig.update_layout(layout)
#
# fig.show()

#
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
ax.set_title('Trajectories for different v_inf')
ax.grid(True)
ax.legend()
plt.show()
# # Построение графиков для разных значений v_inf
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 15))
#
# # График зависимости большой полуоси от времени для разных v_inf
# for i, v_inf in enumerate(v_inf_values):
#     ax1.plot(times_all[i], a_arr_all[i], label=f'v_inf = {v_inf} км/с')
# ax1.set_xlabel('Время, дни')
# ax1.set_ylabel('Большая полуось, км')
# ax1.set_title('Зависимость большой полуоси от времени')
# ax1.legend()
# ax1.grid(True)
#
# # График зависимости аргумента перицентра от времени для разных v_inf
# for i, v_inf in enumerate(v_inf_values):
#     ax2.plot(times_all[i], omega_arr_all[i], label=f'v_inf = {v_inf} км/с')
# ax2.set_xlabel('Время, дни')
# ax2.set_ylabel('Аргумент перицентра, градусы')
# ax2.set_title('Зависимость аргумента перицентра от времени')
# ax2.legend()
# ax2.grid(True)
#
# plt.savefig('images/momentum_energy_for_v_inf')
# # plt.show()
