import numpy as np
import initial_params as ip
import gravity_force as gf
import solar_force as sf

def derivatives(t, y):
    global a_arr, omega_arr, e_arr, nu_arr
    r_i = y[:3]  # Позиция
    v_i = y[3:6]  # Скорость

    r_i_norm = np.linalg.norm(r_i)
    v_i_norm = np.linalg.norm(v_i)

    # Вычисляем cone angle
    r_ = r_i / r_i_norm
    v_t = v_i - (v_i @ r_) * r_
    c = np.cross(r_i, v_i)
    c_ = c / np.linalg.norm(c)
    x_ = np.cross(r_, c_) # трансверсальное напрваление
    delta = np.pi - np.arctan2(v_t @ c_, v_t @ x_)
    # Ускорение
    # a_trac = tf.traction_acceleration(m, e)
    a_grav = gf.gravity_acceleration(r_i)

    # Угол между r и v
    v_i_minus = - v_i
    psi = np.arccos((r_i @ v_i) / (r_i_norm * v_i_norm))
    a_solar = sf.solar_force(r_i, v_i, psi, ip.m0, delta)

    a = a_grav + a_solar

    # print(np.sign(np.dot(a_solar, v_i)))

    return np.hstack((v_i, a))
