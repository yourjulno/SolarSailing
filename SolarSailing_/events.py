import numpy as np
# Событие: Перицентр
def detect_pericenter(t, y):
    r_i = y[:3]  # Позиция
    v_i = y[3:6]  # Скорость
    vr = np.dot(r_i, v_i) / np.linalg.norm(r_i)  # Радиальная скорость
    return vr  # Условие: радиальная скорость равна нулю


# Настройки для перицентра
detect_pericenter.terminal = False  # Не останавливаем на перицентре
detect_pericenter.direction = -1  # Интересует момент, когда vr меняется с - на +


# Функция события: детектирование апоцентра
def detect_apocenter(t, y):
    r_i = y[:3]  # Позиция
    v_i = y[3:6]  # Скорость

    vr = np.dot(r_i, v_i) / np.linalg.norm(r_i)  # Радиальная скорость
    return vr  # Условие: радиальная скорость равна нулю


# Настройки события
detect_apocenter.terminal = True  # Прекратить интеграцию при апоцентре
detect_apocenter.direction = 1  # Интересует только момент, когда vr меняется с + на -
