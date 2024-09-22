
# Импорты
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import constants as const


# Единицы измерения
DU = const.AU # Задано в км
TU = 24 * 60 * 60. # Задано в секундах; сейчас это 1 день
VU = DU / TU

# Константы
MU_sun = 132712440019. * (TU**2/DU**3) # DU^3/TU^2
AU = const.AU / DU


# Функция для вычисления производной
def derivatives(t, y):
    global MU_sun
    r = y[:3]  # Позиция
    v = y[3:6]  # Скорость
    r_magnitude = np.linalg.norm(r)

    # Ускорение
    a = -MU_sun * r / r_magnitude ** 3

    return np.concatenate((v, a))

# Начальные условия
r0 = np.array([AU, 0, 0], dtype = np.float64)  # Начальная позиция (1 AU от Солнца)
v0 = np.array([0, np.sqrt(MU_sun / AU) / 2, 0])  # Начальная скорость (орбитальная скорость Земли)

# Объединение начальных условий в один массив
initial_conditions = np.concatenate((r0, v0))

# Время интегрирования
t_span = (0, 365.25 * 24 * 60 * 60 / TU)  # Один год в секундах
t_eval = np.linspace(t_span[0], t_span[1], num = 10_000)  # Временные точки для оценки


# In[44]:


# Решение системы дифференциальных уравнений
solution = solve_ivp(derivatives, t_span, initial_conditions,
                     method='RK45', t_eval=t_eval, atol = 1e-12, rtol = 1e-12)

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
distances = np.sqrt(x**2 + y**2 + z**2)
speeds = np.sqrt(vx**2 + vy**2 + vz**2)

# Энергия
total_energy = speeds**2 - 2 * MU_sun / distances

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
