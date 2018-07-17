import numpy as np
import matplotlib.pyplot as plt

zTotal = 1  # total depth of the box(stick)
Nz = 100  # space mesh
deltaz = zTotal / Nz  # space step

tTotal = 1000  # total time we will compute
Nt = 1000  # number of time steps
deltat = tTotal / Nt  # time step

kappa = 0.02  # thermal conductivity

'''Array like velocity filed'''
u0 = np.zeros((Nt+1, Nz+1), dtype=np.float32)
u_uniform = np.ones((Nt+1, Nz+1), dtype=np.float32)
u_step_space = np.ones((Nt+1, Nz+1), dtype=np.float32)
u_smooth_step_space = np.ones((Nt+1, Nz+1), dtype=np.float32)
u_gaussian_time = np.ones((Nt+1, Nz+1), dtype=np.float32)

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