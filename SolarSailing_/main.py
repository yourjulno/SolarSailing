import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Константы
G = 6.67430e-11  # Гравитационная постоянная, м^3/(кг*с^2)
M_sun = 1.989e30  # Масса Солнца, кг
AU = 1.496e11  # Один астрономический юнит (расстояние от Земли до Солнца), м

# Начальные условия
r0 = np.array([AU, 0])  # Начальная позиция (1 AU от Солнца)
v0 = np.array([0, 29780])  # Начальная скорость (орбитальная скорость Земли)

# Объединение начальных условий в один массив
initial_conditions = np.concatenate((r0, v0))


# Функция для вычисления производной
def derivatives(t, y):
    r = y[:2]  # Позиция
    v = y[2:]  # Скорость
    r_magnitude = np.linalg.norm(r)

    # Ускорение
    a = -G * M_sun * r / r_magnitude ** 3

    return np.concatenate((v, a))


# Время интегрирования
t_span = (0, 365 * 24 * 60 * 60)  # Один год в секундах
t_eval = np.linspace(t_span[0], t_span[1], num=100000)  # Временные точки для оценки

# Решение системы дифференциальных уравнений
solution = solve_ivp(derivatives, t_span, initial_conditions, method='RK45', t_eval=t_eval)

# Извлечение координат
x = solution.y[0]
y = solution.y[1]
vx = solution.y[2]
vy = solution.y[3]
times = solution.t  # Времена

positions = np.sqrt(x**2 + y**2)
v_perp = np.sqrt(vx**2 + vy**2)
# Вычисление энергии и углового момента
m = 1.0  # Масса спутника (можно считать единичной)

# Кинетическая энергия
kinetic_energy = 0.5 * m * v_perp**2

# Потенциальная энергия
distances = positions
potential_energy = -G * M_sun * m / distances

# Общая энергия
total_energy = kinetic_energy + potential_energy

# Угловой момент
angular_momentum = m * positions * v_perp

# Визуализация графиков зависимости первых интегралов от времени
plt.figure(figsize=(12, 6))
# График энергии
plt.subplot(1, 2, 1)
plt.plot(times, total_energy, label='Общая энергия', color='blue')
plt.plot(times, kinetic_energy, label='Кинетическая энергия', color='green')
plt.plot(times, potential_energy, label='Потенциальная энергия', color='red')
plt.xlabel('Время (с)')
plt.ylabel('Энергия (Дж)')
plt.title('Зависимость энергии от времени')
plt.legend()
plt.grid(True)

# График углового момента
plt.subplot(1, 2, 2)
plt.plot(times, angular_momentum, label='Угловой момент', color='purple')
plt.xlabel('Время (с)')
plt.ylabel('Угловой момент (кг·м²/с)')
plt.title('Зависимость углового момента от времени')
plt.legend()
plt.grid(True)

# Визуализация в 3D
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)

# Отрисовка орбиты
ax.plot(x, y, label='Орбита спутника', color='blue')

# Отрисовка Солнца в начале координат
ax.scatter(0, 0, color='yellow', s=100, label='Солнце')

# Настройка графика
ax.set_xlabel('X (м)')
ax.set_ylabel('Y (м)')
# ax.set_zlabel('Время (с)')
ax.set_title('Траектория спутника вокруг Солнца')
ax.legend()
ax.grid(True)

# Показать график
plt.show()
