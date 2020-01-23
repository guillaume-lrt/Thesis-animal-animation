import numpy as np
import quaternion

def quat_from_alpha_vec(alpha, v):
    a = np.cos(alpha/2)
    b = np.sin(alpha/2)
    return np.quaternion(a, b*v[0], b*v[1], b*v[2])
