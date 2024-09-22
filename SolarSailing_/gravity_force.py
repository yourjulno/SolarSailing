import constants as const
import numpy as np
def gravity_force(r: np.array) -> float:
    r_magnitude = np.linalg.norm(r)
    return - const.MU_sun * r / r_magnitude
