import globalVar
from utils import *

def CN_D(Ts, Tb, kappa, u, Tic, sh=np.zeros((globalVar.Nt + 1, globalVar.Nz + 1), dtype=np.float64)):
    """
    Crank-Nicolson method with Dirichlet B.C.
    Ts as the surface T, Tb as the bottom T.
    """

    Nz = globalVar.Nz
    Nt = globalVar.Nt

    T = np.zeros((Nt + 1, Nz + 1), dtype=np.float64)
    RH = np.zeros(Nz - 1, dtype=np.float64)
    r = globalVar.deltat * kappa / 2 / globalVar.deltaz ** 2

    for n in range(Nt + 1):
        s = globalVar.deltat * u[n, :] / 2 / globalVar.deltaz
        # print(n)
        if n == 0:  # initialize
            for i in range(Nz + 1):
                if i == 0:
                    T[n, i] = Ts
                    continue
                if i == Nz:
                    T[n, i] = Tb
                    continue
                T[n, i] = Tic[i]
            continue
        for i in range(1, Nz):
            RH[i - 1] = r * T[n - 1, i - 1] \
                        + (1 - 2 * r + s[i]) * T[n - 1, i] \
                        + (r - s[i]) * T[n - 1, i + 1] \
                        + (sh[n - 1, i] + sh[n, i])/2 * globalVar.deltat
        d = RH
        d[0] = d[0] + r * Ts
        d[Nz - 1 - 1] = d[Nz - 1 - 1] - (s[Nz - 1] - r) * Tb

        atemp = np.ones(Nz - 1, dtype=np.float64)
        atemp = atemp * (-r)
        btemp = 1 + 2 * r - s[0:Nz-1]
        ctemp = s[0:Nz-1] - r

        Ttemp = solveTri(atemp, btemp, ctemp, d, Nz - 1)
        for i in range(Nz + 1):
            if i == 0:
                T[n, i] = Ts
                continue
            if i == Nz:
                T[n, i] = Tb
                continue
            T[n, i] = Ttemp[i - 1]
    return T
def CN_N(Ts, p, kappa, u, Tic, sh=np.zeros((globalVar.Nt + 1, globalVar.Nz + 1), dtype=np.float64)):
    """
    Crank-Nicolson method with Neumann B.C.
    Ts as the surface T, p as the bottom gradient.
    """

    Nz = globalVar.Nz
    Nt = globalVar.Nt

    T = np.zeros((Nt + 1, Nz + 1), dtype=np.float64)
    RH = np.zeros(Nz, dtype=np.float64)
    r = globalVar.deltat * kappa / 2 / globalVar.deltaz ** 2

    for n in range(Nt + 1):
        s = globalVar.deltat * u[n, :] / 2 / globalVar.deltaz
        # print(n)
        if n == 0:  # initialize
            T[n, :] = Tic
            T[n, 0] = Ts
            continue
        for i in range(1, Nz + 1):
            if i == Nz:
                RH[i - 1] = r * T[n - 1, i - 1] \
                            + (1 - r) * T[n - 1, i] \
                            + (r - s[i]) * globalVar.deltaz * p \
                            + (sh[n - 1, i] + sh[n, i])/2 * globalVar.deltat
                continue
            RH[i - 1] = r * T[n - 1, i - 1] \
                        + (1 - 2 * r + s[i]) * T[n - 1, i] \
                        + (r - s[i]) * T[n - 1, i + 1] \
                        + (sh[n - 1, i] + sh[n, i])/2 * globalVar.deltat
        d = RH
        d[0] = d[0] + r * Ts
        d[Nz - 1] = d[Nz - 1] - (s[Nz] - r) * globalVar.deltaz * p

        atemp = np.ones(Nz, dtype=np.float64)
        atemp = atemp * (-r)
        btemp = 1 + 2 * r - s[0:Nz]
        ctemp = s[0:Nz] - r
        btemp[Nz - 1] = 1 + r

        Ttemp = solveTri(atemp, btemp, ctemp, d, Nz)
        for i in range(Nz + 1):
            if i == 0:
                T[n, i] = Ts
                continue
            T[n, i] = Ttemp[i - 1]
    return T

def CN_D_B_(zTotal, deltaz, tTotal, deltat, Ts, Tb, kappa, u, Tec):
    # Solve the adjoint equation back in time with Crank-Nicolson method with Dirichlet B.C.
    # from Tec back to 0
    # Ts as the surface T, Tb as the bottom T.

    Nz = int(zTotal / deltaz)
    Nt = int(tTotal / deltat)

    T = np.zeros((Nt + 1, Nz + 1), dtype=np.float64)
    RH = np.zeros((Nz - 1), dtype=np.float64)
    r = deltat * kappa / 2 / deltaz ** 2
    s = deltat * u / 2 / deltaz

    for n in range(Nt, -1, -1):
        print(n)
        if n == Nt:  # initialize
            for i in range(Nz + 1):
                if i == 0:
                    T[n, i] = Ts
                    continue
                if i == Nz:
                    T[n, i] = Tb
                    continue
                T[n, i] = Tec(i * deltaz)
            continue
        for i in range(1, Nz):
            RH[i - 1] = r * T[n + 1, i - 1] \
                        + (1 - 2 * r - s) * T[n + 1, i] \
                        + (r + s) * T[n + 1, i + 1]
        d = RH
        d[0] = d[0] + r * Ts
        d[Nz - 2] = d[Nz - 2] + (r + s) * Tb

        atemp = np.ones((Nz - 1), dtype=np.float64)
        btemp = atemp
        ctemp = atemp

        atemp = atemp * (-r)
        btemp = btemp * (1 + 2 * r + s)
        ctemp = ctemp * (-(r + s))

        Ttemp = solveTri(atemp, btemp, ctemp, d, Nz - 1)
        for i in range(Nz + 1):
            if i == 0:
                T[n, i] = Ts
                continue
            if i == Nz:
                T[n, i] = Tb
                continue
            T[n, i] = Ttemp[i - 1]
    return T

def CN_D_B(Ts, Tb, kappa, u, Tec):
    """
    Solve the adjoint equation back in time with Crank-Nicolson method with Dirichlet B.C.
    from Tec back to 0.
    Ts as the surface T, Tb as the bottom T.

    We know from theoretical derivation and the test above that
    the adjoint equation go back in time act exactly like the diffusion equation go forward in time
    but with velocity -u instead of u
    """

    T = CN_D(Ts, Tb, kappa, -u, Tec)
    T = np.flip(T, 0)
    return T
def CN_N_B(Ts, p, kappa, u, Tec):
    """
    Solve the adjoint equation back in time with Crank-Nicolson method with Neumann B.C.
    from Tec back to 0.
    Ts as the surface T, p as the bottom gradient.

    We know from theoretical derivation and the test above that
    the adjoint equation go back in time act exactly like the diffusion equation go forward in time
    but with velocity -u instead of u
    """

    T = CN_N(Ts, p, kappa, -u, Tec)
    T = np.flip(T, 0)
    return T
