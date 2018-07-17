import numpy as np
import matplotlib.pyplot as plt

zTotal = 1  # total depth of the box(stick)
Nz = 100  # space mesh
deltaz = zTotal / Nz  # space step

tTotal = 1000  # total time we will compute
Nt = 1000  # number of time steps
deltat = tTotal / Nt  # time step

kappa = 0.02  # thermal conductivity

'''Function like velocity distribution'''
def u0(PATH='some_value', plot=False):
    ufield = np.zeros((Nt+1, Nz+1), dtype=np.float32)
    if plot:
        plt.plot(np.linspace(0, zTotal, Nz+1), ufield[0, :], color='black')
        plt.xlabel('z')
        plt.ylabel('u')
        plt.title('u0 Velocity field')
        plt.savefig(PATH + '/u_uniform.png')
        plt.cla()
    return ufield
def u_uniform(PATH='some_value', plot=False):
    ufield = np.ones((Nt+1, Nz+1), dtype=np.float32)
    ufield *= -0.05
    if plot:
        plt.plot(np.linspace(0, zTotal, Nz+1), ufield[0, :], color='black')
        plt.xlabel('z')
        plt.ylabel('u')
        plt.title('u_uniform Velocity field')
        plt.savefig(PATH + '/u_uniform.png')
        plt.cla()
    return ufield
def u_step_space(PATH='some_value', plot=False):
    ufield = np.ones((Nt+1, Nz+1), dtype=np.float32)
    for i in range(Nz+1):
        z = i * deltaz
        if z <= zTotal/2:
            ufield[:, i] = 0
        else:
            ufield[:, i] = -0.05
    ufield = ufield * 3
    if plot:
        plt.plot(np.linspace(0, zTotal, Nz+1), ufield[0, :], color='black')
        plt.xlabel('z')
        plt.ylabel('u')
        plt.title('u_step_space Velocity field')
        plt.savefig(PATH + '/u_step_space.png')
        plt.cla()
    return ufield
def u_smooth_step_space(PATH='some_value', plot=False):
    ufield = np.ones((Nt+1, Nz+1), dtype=np.float32)
    for i in range(Nz+1):
        z = i * deltaz
        if z <= zTotal * 0.4:
            ufield[:, i] = 0
        elif z >= zTotal * 0.6:
            ufield[:, i] = -0.05
        else:
            ufield[:, i] = -0.25 * z + 0.1
    ufield = ufield * 3
    if plot:
        plt.plot(np.linspace(0, zTotal, Nz+1), ufield[0, :], color='black')
        plt.xlabel('z')
        plt.ylabel('u')
        plt.title('u_smooth_step_space Velocity field')
        plt.savefig(PATH + '/u_smooth_step_space.png')
        plt.cla()
    return ufield
def u_gaussian_time(PATH='some_value', plot=False):
    ufield = np.ones((Nt+1, Nz+1), dtype=np.float32)
    for n in range(Nt+1):
        t = n * deltat
        ufield[n, :] = - gauss(0.1**2, 0.5, t) * 0.1
    if plot:
        plt.ylim(-0.1, 0.1)
        plt.plot(np.linspace(0, tTotal, Nt+1), ufield[:, 0], color='black')
        plt.xlabel('t')
        plt.ylabel('u')
        plt.title('u_gaussian_time Velocity field')
        plt.savefig(PATH + '/u_gaussian_time.png')
        plt.cla()
    return  ufield

T0 = 0  # Temperature at the top
T1 = 1  # Temperature at the bottom (only useful under Dirichlet B.C.)
p = 0.9188  # Temperature gradient at the bottom (only useful under Neumann B.C.)

# internal heat source H/c
def q_uniform(qs):
    q = np.ones((Nt + 1, Nz + 1), dtype=np.float64) * qs
    return q
def q_continental(qs, hr):

epsilon = 1e-4  # Tolerance scope
MAX = 500  # Max iteration times

outputTimeStep = 0.01
outputSpaceStep = 0.01

'''Function like temperature I.C.'''
def gauss(sigma_2, mu, x):
    """
    A Gaussian like function with max value of 1
    """
    f = np.e**( -(x-mu)**2 / (2*sigma_2) )
    return f

# constant zero
def Tic0():
    T = np.zeros(Nz+1, dtype=np.float32)
    return T

# linear and others
def Tic11():
    z = np.linspace(0, 1, Nz + 1)
    T = z
    return T
def Tic12():
    z = np.linspace(0, 1, Nz + 1)
    T = z**2
    return T

# Guassian 0.5
def Tic21():
    z = np.linspace(0, 1, Nz + 1)
    T = gauss(0.05**2, 0.5, z)
    return T
def Tic22():
    z = np.linspace(0, 1, Nz + 1)
    T = gauss(0.1**2, 0.5, z) / 3
    return T
def Tic23():
    z = np.linspace(0, 1, Nz + 1)
    T = gauss(0.01**2, 0.5, z)
    return T
def Tic24():
    z = np.linspace(0, 1, Nz + 1)
    T = gauss(0.03**2, 0.5, z)
    return T

# Guassian 0.4
def Tic31():
    z = np.linspace(0, 1, Nz + 1)
    T = gauss(0.05**2, 0.4, z)
    return T
def Tic32():
    z = np.linspace(0, 1, Nz + 1)
    T = gauss(0.1**2, 0.4, z) / 3
    return T
def Tic33():
    z = np.linspace(0, 1, Nz + 1)
    T = gauss(0.01**2, 0.4, z) / 3
    return T

# Guassian plus linear and ~
def Tic41():
    z = np.linspace(0, 1, Nz + 1)
    T = -gauss(0.05**2, 0.5, z) / 3 + z
    return T

# Continental geotherm
def Tic_continent():
    z = np.linspace(0, 1, Nz + 1)
    T = 0.9188 * z + 0.0812 * (1-np.e**(-10*z))
    return  T

u = u_uniform
Tic_real = Tic_continent
Tic_guess = Tic11

set_T0_always_positive = False