from utils import *

zTotal = 400  # total depth of the box(stick)  km
Nz =  1 * zTotal  # space mesh
deltaz = zTotal / Nz  # space step

tTotal = 1000  # total time we will compute  Ma
Nt = 1 * tTotal  # number of time steps
deltat = tTotal / Nt  # time step

kappa = 31.6  # thermal diffusivity  km^2 * Ma^-1

Ts = 0  # Temperature at the surface  K
Tb = 1  # Temperature at the bottom (only useful under Dirichlet B.C.)  K
p = 8.96  # Temperature gradient at the bottom (only useful under Neumann B.C.)  K * km^-2

epsilon = 1e-4  # Tolerance scope
MAX = 500  # Max iteration times

outputTimeStep = 0.01
outputSpaceStep = 0.01



def gauss(sigma_2, mu, x):
    """
    A Gaussian like function with max value of 1
    """
    f = np.e**( -(x-mu)**2 / (2*sigma_2) )
    return f
z = np.linspace(0, zTotal, Nz + 1)

'''Array like velocity filed'''
#####################################################################
u0 = np.zeros((Nt+1, Nz+1), dtype=np.float32)
u_uniform = np.ones((Nt+1, Nz+1), dtype=np.float32)

u_step_space = np.ones((Nt+1, Nz+1), dtype=np.float32)
for i in range(Nz + 1):
    zi = i * deltaz
    if zi <= zTotal / 2:
        u_step_space[:, i] = 0
    else:
        u_step_space[:, i] = -0.05

u_smooth_step_space = np.ones((Nt+1, Nz+1), dtype=np.float32)
for i in range(Nz + 1):
        zi = i * deltaz
        if zi <= zTotal * 0.4:
            u_smooth_step_space[:, i] = 0
        elif zi >= zTotal * 0.6:
            u_smooth_step_space[:, i] = -0.05
        else:
            u_smooth_step_space[:, i] = -0.25 * zi + 0.1

u_gaussian_time = np.ones((Nt+1, Nz+1), dtype=np.float32)
for n in range(Nt+1):
    tn = n * deltat
    u_gaussian_time[n, :] = - gauss(0.1**2, 0.5, tn) * 0.1
#####################################################################



'''internal heat source sh=H/c'''
#####################################################################
sh_uniform = np.ones((Nt + 1, Nz + 1), dtype=np.float64) * 0.005
def sh_continental(sh0=23.574, hr=10, u=u0):
    """
    Standard model of exponential distribute with depth radioactive heat source.
    If uplift velocity -u is not equal to zero, then sh0 changes due to erosion.
        u should be uniform in space.
        If u is not uniform in space, only use u at the surface.
    H0: mean heat generation per unit mass at the surface
    c: heat capacity
    hr: At depth z=hr, H is 1/e of its surface value.
    """
    sh = np.ones((Nt + 1, Nz + 1), dtype=np.float64)
    for n in range(Nt + 1):
        shs = sh0 * np.e ** (integrate(u[0:n+1, 0], deltat) / hr)
        for i in range(Nz + 1):
            zi = i * deltaz
            sh[:, i] = shs * np.e ** (-zi/hr)
    return sh
#####################################################################



'''Array like temperature I.C.'''
#####################################################################
# constant
Tic0 = np.zeros(Nz+1, dtype=np.float32)

# linear and others
Tic11 = z / zTotal
Tic12 = (z / zTotal) ** 2

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
Tic41 = -gauss(0.05**2, 0.5, z) / 3 + z / zTotal

# TODO: Continental geotherm
# Continental geotherm
def Tic_continent(p=p, sh0=23.574, hr=10):
    T = Ts + p * z + sh0/kappa * hr**2 * (1 - np.e**(-z/hr))
    return T
#####################################################################



u = u_uniform
Tic_real = Tic11
Tic_guess = Tic12

set_T0_always_positive = False