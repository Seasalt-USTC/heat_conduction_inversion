import numpy as np
import matplotlib.pyplot as plt

'''utils'''
def integrate(f, delta):
    """
    Do integeation on a array with trapezoidal method.
    """
    I = np.sum(f)
    I -= (f[0] + f[-1])/2
    I *= delta
    return I
def norm_abs(a):
    """
    absolute value norm for array a.
    """

    n = np.sum(abs(a))
    return n
def norm_2(a):
    """
    Euclidean norm(L2 norm)
    """

    n = np.sqrt( np.sum(a**2) )
    return n

def solveTri(a, b, c, d, n):
    """
    solve tridiagonal matrix with the central diag b, lower diag a and upper diag c.

    example matrix:
    |b0 c0         | |x0|   |d0|
    |a1 b1 c1      | |x1|   |d1|
    |   a2 b2 c2   | |x2| = |d2|
    |      a3 b3 c3| |x3|   |d3|
    |         a4 b4| |x4|   |d4|
    where n = 5
    """

    x = np.zeros(n, dtype=np.float64)
    cc = np.zeros(n - 1, dtype=np.float64)
    dd = np.zeros(n, dtype=np.float64)

    for i in range(n - 1):
        if i == 0:
            cc[i] = c[i] / b[i]
            dd[i] = d[i] / b[i]
            continue
        cc[i] = c[i] / (b[i] - a[i] * cc[i - 1])
        dd[i] = (d[i] - a[i] * dd[i - 1]) / (b[i] - a[i] * cc[i - 1])
    dd[n-1] = (d[n-1] - a[n-1] * dd[n-1 - 1]) / (b[n-1] - a[n-1] * cc[n-1 - 1])
    for i in range(n - 1, -1, -1):
        if i == n - 1:
            x[i] = dd[i]
            continue
        x[i] = dd[i] - cc[i] * x[i + 1]
    return x
# print (solveTri(1, 2, 3, [1, 1, 1], 3))

def plot_log(logfile, PATH, frame='semilogy'):
    J = []
    with open(logfile) as f:
        for line in f.readlines()[0: -2]:
            J.append( float(line.split()[2]) )
    plt.cla()
    if frame == 'semilogy':
        plt.semilogy(J, label='J')
    elif frame == 'loglog':
        plt.loglog(J, label='J')
    plt.xlabel('k')
    plt.ylabel('J')
    plt.savefig(PATH + '/J.png')

def plot_velocity(u, mode, PATH, zTotal=0, tTotal=0):
    plt.cla()
    if mode == 'space':
        plt.plot(np.linspace(0, zTotal, np.shape(u)[1]), u[0, :], color='black')
        plt.xlabel('z')
        plt.ylabel('u')
        plt.title('u_uniform Velocity field')
        plt.savefig(PATH + '/u.png')
        plt.cla()
    elif mode == 'time':
        plt.plot(np.linspace(0, tTotal, np.shape(u)[0]), u[:, 0], color='black')
        plt.xlabel('t')
        plt.ylabel('u')
        plt.title('u_uniform Velocity field')
        plt.savefig(PATH + '/u.png')
        plt.cla()


