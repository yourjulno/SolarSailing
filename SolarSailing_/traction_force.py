import constants as const
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import initial_params as ip

# скорость истечения топлива
u = 12753 * (const.DU / const.TU)  # км/день

# максимальная сила тяги
f_max = 20

gamma = f_max / u
# Определим переменные
t = sp.symbols('t')

# Определим функцию, зависящую от времени
f = t**3
# Найдем производную по времени
f_derivative = sp.diff(f, t)


# Создаем коллаборационный объект
def derivative_function(value):
    return f_derivative.subs(t, value).evalf()  # Подставляем значение и вычисляем


# Обертка для передачи в solve_ivp
def wrapped_derivative(t, y):
    return np.array([derivative_function(t)])


m0 = [1]
# Решите уравнение с использованием метода RK45
mass = solve_ivp(wrapped_derivative, ip.t_span, m0, method='RK45', t_eval=ip.t_eval)


# coefficients - коэфициенты полинома
def traction_force() -> float:
    global gamma, mass
    force = u * gamma
    for m in mass.y[0]:
        return - force / gamma * m
