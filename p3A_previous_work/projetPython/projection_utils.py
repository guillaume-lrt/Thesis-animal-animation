import numpy as np


def projection_on_plane(d, u=None, offset=0):
    """Orthonormal projection on the plane d.x = offset"""
    #Orthonormal projection
    if u is None:
       u = d
    #Let's create a repere
    if (d[2]==0 and d[1]!=0):
       i = np.array([0, -offset/d[1], 1], dtype=np.float32)
    elif (d[2]!=0) :
       i = np.array([0, 1, (offset-d[1])/d[2]], dtype=np.float32)
    else :
       i = np.array([0, 1, 0], dtype=np.float32)
    i = i/np.linalg.norm(i)
    j = np.cross(d,i)

    #Matrice de passage
    M = np.array((i.T, j.T, d.T))
    print('M: ',M)
    def p(x):
       xu = np.cross(x, u)
       dxu = np.cross(d, xu)
       projection = 1/np.dot(u, d)*dxu-offset*u
       return M.dot(projection)
    return p
