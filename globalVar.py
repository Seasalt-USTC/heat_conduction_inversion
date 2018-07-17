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
# TODO
def q_continental(qs=56.5, k=3.35, qm=30, hr=10):
    """
    Standard model of exponential distribute with depth radioactive heat source.
    The return value is non-dimensional.
    qs: mW*m^-2
    k: W*m^-1*K^-1
    qm: mW*m^-2
    hr: km
    """
    return 0

epsilon = 1e-4  # Tolerance scope
MAX = 500  # Max iteration times

outputTimeStep = 0.01
outputSpaceStep = 0.01

'''Array like temperature I.C.'''
def gauss(sigma_2, mu, x):
    """
    A Gaussian like function with max value of 1
    """
    f = np.e**( -(x-mu)**2 / (2*sigma_2) )
    return f
z = np.linspace(0, 1, Nz + 1)

# constant
Tic0 = np.zeros(Nz+1, dtype=np.float32)

# linear and others
Tic11 = z
Tic12 = z ** 2

# Guassian 0.5
Tic21 = gauss(0.05**2, 0.5, z)
Tic22 = gauss(0.1**2, 0.5, z) / 3
Tic23 = gauss(0.01**2, 0.5, z)
Tic24 = gauss(0.03 ** 2, 0.5, z)

# Guassian 0.4
Tic31 = gauss(0.05**2, 0.4, z)
Tic32 = gauss(0.1 ** 2, 0.4, z) / 3
Tic33 = gauss(0.01 ** 2, 0.4, z) / 3

# Guassian plus linear and others
Tic41 = -gauss(0.05**2, 0.5, z) / 3 + z

# Continental geotherm
Tic_continent = 0.9188 * z + 0.0812 * (1-np.e**(-10*z))

u = u_uniform
Tic_real = Tic_continent
Tic_guess = Tic11

set_T0_always_positive = False