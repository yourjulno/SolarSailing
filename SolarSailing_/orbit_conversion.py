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

    p = a * (1 - e * e)
    distance = p / (1 + e * np.cos(nu))

    cosTrueAnomaly = np.cos(nu)
    sinTrueAnomaly = np.sin(nu)
    # 2. Вектор в орбитальной системе координат
    r_orbital = np.array([cosTrueAnomaly, sinTrueAnomaly, 0]) * distance
    nu_dot = np.sqrt(p * mu) / distance**2
    r_dot = p * (-e * sinTrueAnomaly * nu_dot) / (1 + e * cosTrueAnomaly)**2
    # 3. Вектор скорости в орбитальной системе координат

    velocity = np.array([r_dot * cosTrueAnomaly - distance * nu_dot * np.sin(nu),
                         r_dot * sinTrueAnomaly + distance * nu_dot * np.cos(nu),
                         0])

    # 4. Преобразование вектора в инерциальную систему координат
    R3_W = np.array([[np.cos(Omega), -np.sin(Omega), 0],
                     [np.sin(Omega), np.cos(Omega), 0],
                     [0, 0, 1]])
    R1_i = np.array([[1, 0, 0], [0, np.cos(i), -np.sin(i)], [0, np.sin(i), np.cos(i)]])

    R3_w = np.array([[np.cos(omega), -np.sin(omega), 0],
                     [np.sin(omega), np.cos(omega), 0],
                     [0, 0, 1]])
    R = R3_W @ R1_i @ R3_w

    # 5. Преобразование радиус-вектора и вектора скорости
    r_inertial = R @ r_orbital
    v_inertial = R @ velocity

    return r_inertial, v_inertial


def state_to_orbital_elements(r, v, mu):
    # r - радиус-вектор (numpy array)
    # v - вектор скорости (numpy array)
    # mu - гравитационный параметр (гравитационная постоянная * масса центрального тела)

    # 1. Вычисление величин
    r_norm = np.linalg.norm(r)  # Длина радиус-вектора
    v_norm = np.linalg.norm(v)  # Длина вектора скорости

    # 2.  (a)
    h = np.cross(r, v)  # Угловой момент
    h_norm = np.linalg.norm(h)

    energy = (v_norm ** 2 / 2) - (mu / r_norm)  # Энергия системы
    a = -mu / (2 * energy)  # Большая полуось

    # 3. Эксцентриситет (e)
    e_vector = (1 / mu) * ((v_norm ** 2 - mu / r_norm) * r - np.dot(r, v) * v)
    e = np.linalg.norm(e_vector)

    # 4. Наклонение (i)
    i = np.arccos(h[2] / h_norm)

    # 5. Долгота восходящего узла (Ω)
    N = h / np.linalg.norm(h)  # Вектор нормали к плоскости орбиты
    N_norm = np.linalg.norm(N)

    if N_norm != 0:
        Omega = np.arctan2(N[1], N[0])
    else:
        Omega = 0  # Если N = 0, то долгота восходящего узла неопределена

    # 6. Аргумент перицентра (ω)
    if e > 0:
        omega = np.arccos(np.dot(N, e_vector) / (N_norm * e))
        if e_vector[2] < 0:
            omega = 2 * np.pi - omega
    else:
        omega = 0  # Если эксцентриситет равен нулю, то ω неопределен

    # 7. Истинная аномалия (ν)
    if e > 0:
        nu = np.arccos(np.dot(e_vector, r) / (e * r_norm))
        if np.dot(r, v) < 0:
            nu = 2 * np.pi - nu
    else:
        nu = 0  # Если эксцентриситет равен нулю, то ν неопределен

    return a, e, omega, nu
