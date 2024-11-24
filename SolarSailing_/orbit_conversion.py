import numpy as np


def orbital_elements_to_state(a, e, i, Omega, omega, nu, mu):
    # a - большая полуось (м)
    # e - эксцентриситет
    # i - наклонение (рад)
    # Omega - долгота восходящего узла (рад)
    # omega - аргумент перицентра (рад)
    # nu - истинная аномалия (рад)
    # mu - гравитационный параметр (м^3/с^2)

    # 1. Вычисление радиус-вектора r
    r = a * (1 - e**2) / (1 + e * np.cos(nu))  # Радиус-вектор в перигее

    # 2. Вектор в орбитальной системе координат
    r_orbital = np.array([r * np.cos(nu), r * np.sin(nu), 0])

    # 3. Вектор скорости в орбитальной системе координат
    v_orbital = np.array([
        -np.sqrt(mu/a) * np.sin(nu),
        np.sqrt(mu/a) * (e + np.cos(nu)),
        0
    ])

    # 4. Преобразование вектора в инерциальную систему координат
    R = np.array([
        [np.cos(Omega) * np.cos(omega) - np.sin(Omega) * np.sin(omega) * np.cos(i),
         -np.cos(Omega) * np.sin(omega) - np.sin(Omega) * np.cos(omega) * np.cos(i),
         np.sin(Omega) * np.sin(i)],

        [np.sin(Omega) * np.cos(omega) + np.cos(Omega) * np.sin(omega) * np.cos(i),
         -np.sin(Omega) * np.sin(omega) + np.cos(Omega) * np.cos(omega) * np.cos(i),
         -np.cos(Omega) * np.sin(i)],

        [np.sin(omega) * np.sin(i),
         np.cos(omega) * np.sin(i),
         np.cos(i)]
    ])

    # 5. Преобразование радиус-вектора и вектора скорости
    r_inertial = R @ r_orbital
    v_inertial = R @ v_orbital

    return r_inertial, v_inertial
