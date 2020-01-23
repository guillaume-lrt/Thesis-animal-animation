import numpy as np
import matplotlib as plt

R=0.1
def plot_cylinder(ax, a, b):
    v = b - a
    mag = np.linalg.norm(v)
    v /= mag
    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    
    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    #normalize n1
    n1 /= np.linalg.norm(n1)
    
    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)

    mag = np.linalg.norm(v) 
    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 2)
    theta = np.linspace(0, 2 * np.pi, 100)
    rsample = np.linspace(0, R, 2) 
    #use meshgrid to make 2d arrays
    t, theta2 = np.meshgrid(t, theta)

    rsample, theta = np.meshgrid(rsample, theta)
    X, Y, Z = [a[i] + v[i] * t + R * np.sin(theta2) * n1[i] + R * np.cos(theta2) * n2[i] for i in [0, 1, 2]]
    ax.plot_surface(X, Y, Z, color='blue')
