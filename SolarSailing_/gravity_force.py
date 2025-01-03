import constants as const
import numpy as np
def gravity_acceleration(r: np.array) -> np.array:
    r_magnitude = np.linalg.norm(r)
    return - const.MU_in_units * r / (r_magnitude ** 3)
